# SNP calling with GATK-UnifiedGenotyper

## Overview
Here we will use the program [GATK](https://gatk.broadinstitute.org/hc/en-us) to identify short genetic variants (single nucleotide polymorphisms - SNPs, and indels). This step is memory demanding, so I usually split it per chromosome and unplaced scaffolds. However, given that you will run this script in your local computer, we won't split the run but rather will try to do it all in once. As a consequence, this run will probably take several days to complete. If memory problems arise, then we would need to split the analysis per chromosome (44 runs) and unplaced scaffolds (single run).

## Required programs
- GATK version 3.8-0 (download the precompiled jar file to your computer and store it in an appropriate directory)

## Input files
Inside your `data/` directory, create a new directory called `04-vcf-files/` and go there:
```
mkdir 04-vcf-files
cd 04-vcf-files
```
Create a TXT file called `bam_files_UG.list`. This file will list the path to the BAM files generated before. This file should look something like this:
```
/<path>/data/03-bam-files/pool1.sort.MarkDup.bam
/<path>/data/03-bam-files/pool2.sort.MarkDup.bam
/<path>/data/03-bam-files/pool3.sort.MarkDup.bam
/<path>/data/03-bam-files/pool4.sort.MarkDup.bam
/<path>/data/03-bam-files/pool5.sort.MarkDup.bam
/<path>/data/03-bam-files/pool6.sort.MarkDup.bam
```

## Run GATK
Script `02-SNPcalling-UG.sh`:
```bash
#!/bin/bash
# Script to call genetic variants (SNPs and indels) using GATK

# NOTE: Update the <paths> to the files and programs required, as well as
# the RAM memory (`-Xmx60G`) available in your computer

# Set environment variables
REF_FILE='/<path>/data/00-genome/GCF_015228065.1_NSTDA_Pmon_1_genomic.fna'
BAM_LIST='/<path>/data/04-vcf-files/bam_files_UG.list'
GATK_JAR='/<path>/GATK/3.8-0/GenomeAnalysisTK.jar'
OUTFILE_PREFIX='Pmon_6_pools.UG.rawVar'

# Run GATK-UnifiedGenotyper to generate a single VCF file for all pool samples.
java -Xmx60G -jar $GATK_JAR \
-T UnifiedGenotyper \
-R $REF_FILE \
-I $BAM_LIST \
-o $OUTFILE_PREFIX.vcf

```
Good luck!
