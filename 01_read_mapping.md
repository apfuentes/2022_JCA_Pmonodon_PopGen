
# Map reads against the reference genome assembly

## Overview
Here we will align the clean paired-end reads of each pool against the reference genome of *P. monodon* using [BWA](http://bio-bwa.sourceforge.net/). From this, we will get a separate [BAM](https://software.broadinstitute.org/software/igv/BAM) file for each pool. Then we will do some post-processing (mark duplicates, add read groups, generate an index file) to prepare the BAM files for variant calling with [GATK](https://gatk.broadinstitute.org/hc/en-us). Finally, we will generate read mapping quality statistics with [Qualimap](http://qualimap.conesalab.org)

## Required Programs
Make sure these programs are installed in your computer:
- bioinfo-tools
- bwa/0.7.17
- samtools/1.9
- picard/2.20.4
- QualiMap/2.2.1
- java/sun_jdk1.8.0_151 (for Picard)
- java/sun_jdk1.7.0_25 (for Qualimap)

**NOTE:** Please make sure that you have two java versions (JDK 1.8 and 1.7) installed in your machine. The reason is that Picard requires version 1.8 while Qualimap version 1.7. You may need to investigate how to change java version in your machine before running Qualimap and add a command accordingly (see line indicated with >>>>> <<<<<).

## Input file
Let's start by creating a TXT file listing the name of your samples, which should be part of the name of the fastq files (R1 and R2).
Assuming your files are named something like this:
```
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
sample3_R1.fastq.gz
sample3_R2.fastq.gz
sample4_R1.fastq.gz
sample4_R2.fastq.gz
sample5_R1.fastq.gz
sample5_R2.fastq.gz
```
The TXT file called `sample.list`, should look like this (each name in a new line):
```
sample1
sample2
sample3
sample4
sample5
```

## Steps
1. Generate separate bash scripts for each sample
2. Run each of the bash scripts

### 1. Generate separate bash scripts for each sample

**NOTE:** Before running this script, make sure to update the `/path` to the files in your computer as well as the RAM memory (`-Xmx60G`) and threads (`-t 8`) available.

Script `01-create_scripts_read_mapping.sh`:
```bash
#!/bin/bash
# Script to align paired-end reads against a reference genome using BWA-mem,
# post-process the generated BAM file, and generate mapping quality statistics

# Make sure these programs are installed in your computer
#bioinfo-tools
#bwa/0.7.17
#samtools/1.9
#picard/2.20.4
#QualiMap/2.2.1
#java/sun_jdk1.8.0_151 (for Picard)
#java/sun_jdk1.7.0_25 (for Qualimap)

# Make sure to update the paths to the files in your computer as well as the
# RAM memory (`-Xmx60G`) and threads (`-t 8`) available

# Path to the file listing sample names
FASTQ_DIR='/path/data/02-clean-reads'
SAMPLE_LIST='/path/data/03-bam-files/sample.list'
# For testing
#FASTQ_DIR='/proj/snic2020-2-19/private/horse_mackerel/data/02-clean-reads/test-10Kreads'
#SAMPLE_LIST='/proj/snic2020-2-19/private/horse_mackerel/data/02-clean-reads/test-10Kreads/03-bam-files/sample.list'

# Loop to generate the separate scripts
while read -r line; do

	# Path to the fastq files
	R1=$FASTQ_DIR/${line}_R1.fastq.gz
	R2=$FASTQ_DIR/${line}_R2.fastq.gz

	if [ -e $R1 ] && [ -e $R2 ] ; then
		echo "#!/bin/bash" > map_reads_$line.sh
echo "
# Set environment variables
WORK_DIR='/path/data/03-bam-files'
REF_FILE='/path/data/00-genome/GCF_015228065.1_NSTDA_Pmon_1_genomic.fna'
PICARD_JAR='/path/tools/picard.jar'
# For testing
#WORK_DIR='/proj/snic2020-2-19/private/horse_mackerel/data/02-clean-reads/test-10Kreads/03-bam-files'
#REF_FILE='/proj/snic2020-2-19/private/horse_mackerel/data/00-genome/fTraTra1_1.curated_primary.20200310.fa'
#PICARD_JAR='/sw/apps/bioinfo/picard/2.20.4/rackham/picard.jar'

# Directory where the BAM files will be stored.
cd \$WORK_DIR

# Create directory for temporal files.
if [ -d \$WORK_DIR/tmp ]; then echo \"tmp/ exists in WORK_DIR\"; else mkdir \$WORK_DIR/tmp; fi

echo -e \$(date -u) \": Read mapping began...\"

# Map reads to reference genome using several threads [-t 8] and mark split alignment [-M]; then sort reads and generate the bam file index.
bwa mem -M -t 8 -R '@RG\tID:${line}\tSM:${line}' \${REF_FILE} ${R1} ${R2} | samtools view -@ 8 -b -S - > \$WORK_DIR/${line}.bam
samtools sort -@ 8 -T \$WORK_DIR/tmp -o \$WORK_DIR/${line}.sort.bam \$WORK_DIR/${line}.bam && rm \$WORK_DIR/${line}.bam
samtools index -@ 8 \$WORK_DIR/${line}.sort.bam

echo -e \$(date -u) \": Read mapping and bam file indexing ended...\"

# Mark duplicate reads.
if [ -d \$WORK_DIR/MarkDup_metrics ]; then echo \"MarkDup_metrics/ exists\"; else mkdir \$WORK_DIR/MarkDup_metrics; fi
echo -e \$(date -u) \": Mark duplicates began...\"

java -Xmx60G -jar \${PICARD_JAR} MarkDuplicates \
I=\$WORK_DIR/${line}.sort.bam O=\$WORK_DIR/${line}.sort.MarkDup.bam M=\$WORK_DIR/MarkDup_metrics/${line}.MarkDup.txt \
ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT TMP_DIR=\$WORK_DIR/tmp && rm \$WORK_DIR/${line}.sort.bam

echo -e \$(date -u) \": Mark duplicates ended...\"
echo -e \$(date -u) \": Add read groups began...\"

# Add read groups to the bam file.
java -Xmx60G -jar \${PICARD_JAR} AddOrReplaceReadGroups \
I=\$WORK_DIR/${line}.sort.MarkDup.bam O=\$WORK_DIR/${line}.sort.MarkDup.RG.bam \
RGID=${line} RGLB=${line} RGPL=illumina RGPU=${line} RGSM=${line} && rm \$WORK_DIR/${line}.sort.MarkDup.bam

# Create an index for the final bam file.
samtools index -@ 8 \$WORK_DIR/${line}.sort.MarkDup.RG.bam

# Clean up tmp files
rm -R \$WORK_DIR/tmp

echo -e \$(date -u) \": Add read groups and bam index creation ended...\"

# Obtain mapping quality summary statistics.
#module load QualiMap/2.2.1
#unset DISPLAY  # Turn display off to avoid problems with X11 in Uppmax.
#>>>>> Juan, HERE add a command to switch from java/sun_jdk1.8.0_151 => java/sun_jdk1.7.0_25 <<<<<

if [ -d \$WORK_DIR/Qualimap_results ]; then echo \"Qualimap_results/ exists\"; else mkdir \$WORK_DIR/Qualimap_results; fi

echo -e \$(date -u) \": Qualimap began...\"

qualimap --java-mem-size=60G bamqc -bam \$WORK_DIR/${line}.sort.MarkDup.RG.bam -ip -outdir \$WORK_DIR/Qualimap_results/${line} -outformat html

echo -e \$(date -u) \": Qualimap ended...\"

" >> map_reads_$line.sh
			echo $line
			#bash $line.sh

			else
				echo "No reads found for $line"
			fi
		done < "$SAMPLE_LIST"

```
Toy example of one of the generated bash file (there might be a few differences, since it was generated during testing):
```bash
#!/bin/bash

# Set environment variables
WORK_DIR='/proj/snic2020-2-19/private/horse_mackerel/data/02-clean-reads/test-10Kreads/03-bam-files'
REF_FILE='/proj/snic2020-2-19/private/horse_mackerel/data/00-genome/fTraTra1_1.curated_primary.20200310.fa'
PICARD_JAR='/sw/apps/bioinfo/picard/2.20.4/rackham/picard.jar'

#WORK_DIR='/path/data/03-bam-files'
#SAMPLE_LIST='/path/data/03-bam-files/sample.list'
#REF_FILE='/path/data/00-genome/GCF_015228065.1_NSTDA_Pmon_1_genomic.fna'
#PICARD_JAR='/path/tools/picard.jar'

# Directory where the BAM files will be stored.
cd $WORK_DIR

# Create directory for temporal files.
if [ -d $WORK_DIR/tmp ]; then echo "tmp/ exists in WORK_DIR"; else mkdir $WORK_DIR/tmp; fi

echo -e $(date -u) ": Read mapping began..."

# Map reads to reference genome using several threads [-t 8] and mark split alignment [-M]; then sort reads and generate the bam file index.
bwa mem -M -t 8 -R '@RG\tID:1a\tSM:1a' ${REF_FILE} /proj/snic2020-2-19/private/horse_mackerel/data/02-clean-reads/test-10Kreads/1a_R1.fastq.gz /proj/snic2020-2-19/private/horse_mackerel/data/02-clean-reads/test-10Kreads/1a_R2.fastq.gz | samtools view -@ 8 -b -S - > $WORK_DIR/1a.bam
samtools sort -@ 8 -T $WORK_DIR/tmp -o $WORK_DIR/1a.sort.bam $WORK_DIR/1a.bam && rm $WORK_DIR/1a.bam
samtools index -@ 8 $WORK_DIR/1a.sort.bam

echo -e $(date -u) ": Read mapping and bam file indexing ended..."

# Mark duplicate reads.
if [ -d $WORK_DIR/MarkDup_metrics ]; then echo "MarkDup_metrics/ exists"; else mkdir $WORK_DIR/MarkDup_metrics; fi
echo -e $(date -u) ": Mark duplicates began..."

java -Xmx60G -jar ${PICARD_JAR} MarkDuplicates I=$WORK_DIR/1a.sort.bam O=$WORK_DIR/1a.sort.MarkDup.bam M=$WORK_DIR/MarkDup_metrics/1a.MarkDup.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT TMP_DIR=$WORK_DIR/tmp && rm $WORK_DIR/1a.sort.bam

echo -e $(date -u) ": Mark duplicates ended..."
echo -e $(date -u) ": Add read groups began..."

# Add read groups to the bam file.
java -Xmx60G -jar ${PICARD_JAR} AddOrReplaceReadGroups I=$WORK_DIR/1a.sort.MarkDup.bam O=$WORK_DIR/1a.sort.MarkDup.RG.bam RGID=1a RGLB=1a RGPL=illumina RGPU=1a RGSM=1a && rm $WORK_DIR/1a.sort.MarkDup.bam

# Create an index for the final bam file.
samtools index -@ 8 $WORK_DIR/1a.sort.MarkDup.RG.bam

# Clean up tmp files
rm -R $WORK_DIR/tmp

echo -e $(date -u) ": Add read groups and bam index creation ended..."

# Obtain mapping quality summary statistics.
#module load QualiMap/2.2.1
#unset DISPLAY  # Turn display off to avoid problems with X11 in Uppmax.

if [ -d $WORK_DIR/Qualimap_results ]; then echo "Qualimap_results/ exists"; else mkdir $WORK_DIR/Qualimap_results; fi

echo -e $(date -u) ": Qualimap began..."

qualimap --java-mem-size=60G bamqc -bam $WORK_DIR/1a.sort.MarkDup.RG.bam -ip -outdir $WORK_DIR/Qualimap_results/1a -outformat html

echo -e $(date -u) ": Qualimap ended..."

```

### 2. Run each of the bash scripts generated

There should be a single script for each pool sample. While this could be run either in parallel or consecutively in an automatic way, I would recommend to launch each script manually after the other, once you make sure every works.

## Generate a single report of all qualimap reports
Once all the BAM files are generated, combine the Qualimap reports using MultiQC:
```
cd /path/03-bam-files/Qualimap_results
#cd /proj/snic2020-6-128/private/a_obtectus_QTLmap/03-bam-files/Qualimap_results

#module load bioinfo-tools MultiQC/1.12
multiqc .

```
