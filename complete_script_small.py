import json
import os
import sys
import threading 

ref_directory = "/ebio/ecnv_projects/common_resourses/data/reference_genome/GRCz11/"
cmd_picard = "/ebio/ecnv_projects/common_resourses/data/software/picard-tools-1.130/picard.jar"
cmd_picard2 = "/ebio/ecnv_projects/common_resourses/data/software/picard-tools-1.130/picard2.jar"
cmd_GenomeAnalysisTK = "/ebio/ecnv_projects/common_resourses/data/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar"
cmd_GATK4 = "/ebio/ecnv_projects/common_resourses/data/software/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar"

ref = ref_directory + "Danio_rerio.GRCz11.dna.primary_assembly.fa"
ref_dict = ref_directory + "Danio_rerio.GRCz11.dna.primary_assembly.dict"
snp_db = ref_directory + "danio_rerio.vcf.gz"
chr_big = ref_directory + "chr_big.list"
#chr_big = "chr_big.list"

os.system("cp {}chr_small.list chr_small.list".format(ref_directory))
chr_small = "chr_small.list"

#input files

def check_config(data):
    if "samples" not in data:
        raise Exception("You have not provided sample list")
    if "bam_prefix" not in data:
        raise Exception("You have not provided a prefix for input bam file")


def import_config(filename):
    with open(filename) as f:
        data = json.load(f)
    check_config(data)
    return(data)

def mark_duplicates(sample, input_bam_prefix, markduplicates_postfix):
    input_bam = sample + input_bam_prefix
    output_bam = sample + markduplicates_postfix
    metrics_file = sample + "_dedup_reads_metrics.txt"
    command = "java -jar -Xmx32g {} MarkDuplicates INPUT={} OUTPUT={} METRICS_FILE={}".format(
        cmd_picard, input_bam, output_bam, metrics_file)
    print(command)
    os.system(command)
    command2 = "samtools index {}".format(output_bam)
    print(command2)
    os.system(command2)

def dispatch_samples(samples, target_function, parameter1=None, parameter2=None):
    threads = []
    for sample in samples:
        thread = threading.Thread(target=target_function, args=(sample, parameter1, parameter2))
        threads += [thread]
        thread.start()
    
    for x in threads:
        x.join()

def dispatch_chromosomes(chromosomes, sample, target_function, parameter1=None, parameter2=None):
    threads = []
    for chromosome in chromosomes:
        thread = threading.Thread(target=target_function, args=(sample, chromosome, parameter1, parameter2))
        threads += [thread]
        thread.start()
    
    for x in threads:
        x.join()

def splitncigar(sample, input_bam_postfix, output_bam_postfix):
    input_bam = sample + input_bam_postfix
    output_bam = "{}{}".format(sample, output_bam_postfix)
    command = "java -jar -Xmx32g {} SplitNCigarReads -R {} -I {} -O {}".format(cmd_GATK4, ref, input_bam, output_bam)
    print(command)
    os.system(command)

def gather_bams_rg(sample, bam_postfix, rg_postfix):
 
    input_bam = "{}{}".format(sample, bam_postfix)
    output_bam = "{}{}".format(sample, rg_postfix)
    read_group = "RGID={0} RGLB={0} RGPL=illumina RGPU=unit1 RGSM={0}".format(sample)
    command = "java -jar -Xmx32g {} AddOrReplaceReadGroups I={} OUTPUT={} {} ".format(cmd_picard, input_bam, output_bam, read_group)
    print(command)
    os.system(command)

def bsqr(sample, io_rg_postfix, io_bsqr):
    input_bam = sample + io_rg_postfix
    output_table = sample + ".recal.table"
    output_bam = sample + io_bsqr
    command = "java -jar -Xmx32g {} BaseRecalibrator -I {} -R {} --known-sites {} -O {}".format(cmd_GATK4, input_bam, ref, snp_db, output_table)
    print(command)
    os.system(command)
    command2 = "java -jar -Xmx32g {} ApplyBQSR -R {} -I {} --bqsr-recal-file {} -O {}".format(cmd_GATK4, ref, input_bam, output_table, output_bam)
    print(command2)
    os.system(command2)

def call_variants(sample, chromosome, io_bsqr, io_gvcf):
    input_bam = sample + io_bsqr
    output_gvcf = "gvcf/{}.{}{}".format(sample, chromosome, io_gvcf)
    command = "java -jar -Xmx26g {} HaplotypeCaller -R {} -I  {} -O {} -L {} -ERC GVCF --verbosity ERROR".format(cmd_GATK4, ref, input_bam, output_gvcf, chromosome)
    print(command)
    os.system(command)

def genotype_gvcfs(samples, chromosomes, sample, io_gvcf):
    cohort_gvcfs = []
    output_vcfs = []
    for chromosome in chromosomes:
        input_gvcfs = []
        cohort_gvcf = "cohort.{}.g.vcf.gz".format(chromosome)
        output_vcf = "output.{}.vcf.gz".format(chromosome)
        for sample in samples:
            input_gvcfs.append(" -V gvcf/{}.{}{} ".format(sample, chromosome, io_gvcf))
            command = "java -jar {} CombineGVCFs -R {}  {} -O {}".format(cmd_GATK4, ref, "".join(input_gvcfs), cohort_gvcf)
            print(command)
            os.system(command)
        
        cohort_gvcfs.append(cohort_gvcf)
        command_genotype = "java -jar {} GenotypeGVCFs -R {} -V {} -O {}".format(cmd_GATK4, ref, cohort_gvcf, output_vcf)
        print(command_genotype)
        os.system(command_genotype)
        output_vcfs.append(output_vcf)
    
    return output_vcfs

def process_vcfs(output_vcfs):
    cat_vcfs = " -V ".join(output_vcfs)
    print(cat_vcfs)
    
    command_cat = "java -cp {} org.broadinstitute.gatk.tools.CatVariants -R {} -V {} -out merged_raw_variants.vcf".format(cmd_GenomeAnalysisTK, ref, cat_vcfs)  
    print(command_cat)
    os.system(command_cat)
    
    os.system("rm merged_raw_variants.vcf.idx")
    
    os.system("java -jar {} SortVcf I=merged_raw_variants.vcf O=sorted_merged_raw_variants.vcf SEQUENCE_DICTIONARY={}".format(cmd_picard2, ref_dict))
    
    #remove idx file generated by picard
    os.system("rm sorted_merged_raw_variants.vcf.idx")
    
    # select snps and indels
    select_snps_cmd = "java -jar {} -jdk_deflater -jdk_inflater -T SelectVariants -R {} -V sorted_merged_raw_variants.vcf -selectType SNP -o raw_snps.vcf".format(cmd_GenomeAnalysisTK,ref)
    os.system(select_snps_cmd)
    
    select_indels_cmd = "java -jar {} -jdk_deflater -jdk_inflater -T SelectVariants -R {} -V sorted_merged_raw_variants.vcf -selectType INDEL -o raw_indels.vcf".format(cmd_GenomeAnalysisTK,ref)
    os.system(select_indels_cmd) 
    
    #filter variants
    filter_snps_cmd = "java -jar -Xmx24g {} -T VariantFiltration -R {} -V raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName 'basic_snp_filter' -o raw_merged_filtered_snps.vcf ".format(cmd_GenomeAnalysisTK, ref)
    os.system(filter_snps_cmd)
    
    filter_indels_cmd = "java -jar -Xmx24g {} -T VariantFiltration -R {} -V raw_indels.vcf  --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName 'basic_indel_filter' -o raw_merged_filtered_indels.vcf".format(cmd_GenomeAnalysisTK, ref)
    os.system(filter_indels_cmd)


f_config = sys.argv[1]
print(sys.argv)
config_data = import_config(f_config)

print(config_data)

samples = config_data["samples"]
print(samples)
input_bam_prefix = config_data["bam_prefix"] 

with open(chr_big) as f:
    chromosomes_input = f.readlines()

chromosomes = []
for chrom in chromosomes_input:
    chromosomes.append(chrom.rstrip())

chromosomes.append(chr_small)

os.system("mkdir gvcf")


#Mark Duplicates
io_markduplicates_postfix = "_dedup_reads.bam"
dispatch_samples(samples, mark_duplicates, input_bam_prefix, io_markduplicates_postfix)

# SplitNCigarReads per sample
io_splitncigar_postfix = "_dedup_reads.split.bam"

dispatch_samples(samples, splitncigar, io_markduplicates_postfix, io_splitncigar_postfix)

#Assign read groups to split bams
io_rg_postfix = ".split.RG.bam"

dispatch_samples(samples, gather_bams_rg, io_splitncigar_postfix, io_rg_postfix)

   
# Base recalibration

io_bsqr = ".recal.bam"

dispatch_samples(samples, bsqr, io_rg_postfix, io_bsqr)

# HaplotypeCaller

io_gvcf = ".g.vcf.gz"

for sample in samples:
    dispatch_chromosomes(chromosomes, sample, call_variants, io_bsqr, io_gvcf)

#Variant calls and variant processing

output_vcfs = genotype_gvcfs(samples, chromosomes, sample, io_gvcf)

process_vcfs(output_vcfs)
        


