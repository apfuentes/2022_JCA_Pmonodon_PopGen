# Reference genome of the black tiger shrimp (*Penaeus monodon*)

The assembly covered 2.39 Gb (92.3% of the estimated genome size) and contains 44 pseudomolecules + numerous unplaced scaffolds. More details in the website of the project [here](https://www.biotec.or.th/pmonodon/index.php)

- Files from NCBI:
NSTDA_Pmon_1
https://www.ncbi.nlm.nih.gov/assembly/GCF_015228065.1/
https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Penaeus_monodon/100/
Assembly name	Assembly accession	Submitter	Assembly date	Reference/Alternate	Assembly content
NSTDA_Pmon_1	GCF_015228065.1	SAFE-Aqua	11-05-2020	Reference	45 assembled chromosomes; unplaced scaffolds

- Files from ensembl:
https://rapid.ensembl.org/Penaeus_monodon_gca015228065v1/Info/Index?db=core
Penaeus_monodon-GCA_015228065.1-unmasked.fa.gz     27-Jan-2022 22:58           567515462

Provider: SAFE-Aqua http://www.safeaqua-project.net
Penaeus_monodon-GCA_015228065.1-2021_11-genes.gtf 27-Jan-2022 08:59             7940275
Penaeus_monodon-GCA_015228065.1-2021_11-genes.gff3 27-Jan-2022 20:50             8888429

- Paper fo the genome assembly:
A chromosome-level assembly of the black tiger shrimp (Penaeus monodon) genome facilitates the identification of growth-associated genes https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13357

## Obtain index files for the reference genome
For read mapping and variant calling using BWA and GATK, respectively, it is necessary to generate a group of indexes and dictionaries for the reference genome. For this, I used the bash script called `00-prepare_ref_genome_for_BWA_GATK.sh`:

```bash
#!/bin/bash

# Load required software (only applicable to Uppmax).
#module load bioinfo-tools
#module load bwa/0.7.17
#module load samtools/1.9
#module load picard/2.20.4

# Set environment variables.
REF_FILE='/path/00-genome/GCF_015228065.1_NSTDA_Pmon_1_genomic.fna'
WORK_DIR='/path/00-genome'
PICARD_JAR='/sw/apps/bioinfo/picard/2.20.4/rackham/picard.jar'

cd $WORK_DIR
# Generate BWA file index.
# where -a bwtsw specifies that we want to use the indexing algorithm that is capable of handling the whole human genome.
# Expected Result: This creates a collection of files used by BWA to perform the alignment.
bwa index -a bwtsw ${REF_FILE}

echo "################# bwa index done" ;

# Generate fasta file index.
# This creates a file called reference.fa.fai, with one record per line for each of the contigs in the FASTA reference file. Each record is composed of the contig name, size, location, basesPerLine and bytesPerLine.
samtools faidx ${REF_FILE}

echo "################# samtools faidx done" ;

# Generate sequence dictionary.
# Note that this is the new syntax for use with the latest version of Picard. Older versions used a slightly different syntax because all the tools were in separate jars, so you'd call e.g. java -jar CreateSequenceDictionary.jar directly.
# This creates a file called reference.dict formatted like a SAM header, describing the contents of your reference FASTA file.
# >>>> Update the RAM memory -Xmx48g >>>
java -Xmx48g -jar $PICARD_JAR CreateSequenceDictionary \
REFERENCE=${REF_FILE} \
OUTPUT=${REF_FILE}.dict

echo "################# picard dictionary done" ;

```
**Runtime: 00-00:25:52**

Rename genome dictionary file, because GATK will expect it to be named as `genome.dict` not `genome.fasta.dict`:
```
cd /path/00-genome
mv GCF_015228065.1_NSTDA_Pmon_1_genomic.fna.dict GCF_015228065.1_NSTDA_Pmon_1_genomic.dict
```
