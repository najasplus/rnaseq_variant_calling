# rnaseq_variant_calling
Variant calling pipeline using GATK 3.8 and 4 and picard from STAR-aligned bam files based on GATK Best-Practices pipeline https://software.broadinstitute.org/gatk/best-practices/workflow?id=11164

The script pre-processes the provided bam files using picard MarkDuplicates, GATK SplitNCigarReads, BaseRecalibrator. The variants are called with GATK4 HaplotypeCaller in GVCF mode.

The output of the script are NCigar-split and recalibrated bam files, vcf files with raw and filtered snps and indels combined for all input samples (sample information preserved).

The script uses zebrafish GRCz11 genome as a reference, this can be modified accoriding to the reference used for producing alignments. 

The script parallelizes the file processing per sample or per chromosomes, therefore chromosome names list should be provided in the path of reference genome. 

# Input

* input_parameters description as .json file. Should include "samples" as a list and "bam_prefix".
Example:

```python
{
	"samples" : ["SampNr1", "SampNr2", "SampNr3"],
	"bam_prefix": ".Aligned.sortedByCoord.out.bam"

}

```

* One or multiple .bam files produced by splice-aware aligner (i.e. STAR). The bam names should start with sample name followed by "bam_prefix" (as in input_paramteres.json file).
Example:

```
SampNr1.Aligned.sortedByCoord.out.bam
SampNr2.Aligned.sortedByCoord.out.bam
SampNr3.Aligned.sortedByCoord.out.bam
```

# Running the script

In the directory where the input bams (or links to them) are located:

```bash
python3 /path/to/script/complete_script_small.py /path/to/input_parameters.json
```



