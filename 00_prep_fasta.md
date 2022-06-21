# Reference genome of the black tiger shrimp (Penaeus monodon)

The assembly covered 2.39 Gb (92.3% of the estimated genome size) and contains 44 pseudomolecules, see the website of the project [here](https://www.biotec.or.th/pmonodon/index.php)

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


## Replace contig names in the genome assembly

The contig names in the NCBI assembly are too long, e.g. `NC_051386.1 Penaeus monodon isolate SGIC_2016 chromosome 1, NSTDA_Pmon_1, whole genome shotgun sequence`. Thus, I will replace the fasta headers using the SeqKit package [replace](https://bioinf.shenwei.me/seqkit/usage/#replace) and a key-value file that I created in excel and saved as a TXT file (two columns separated by a tab, 1st column for the original name, and 2nd column for the new name):
```
NC_051386.1 Penaeus monodon isolate SGIC_2016 chromosome 1, NSTDA_Pmon_1, whole genome shotgun sequence	chr_1
NC_051387.1 Penaeus monodon isolate SGIC_2016 chromosome 2, NSTDA_Pmon_1, whole genome shotgun sequence	chr_2
NC_051388.1 Penaeus monodon isolate SGIC_2016 chromosome 3, NSTDA_Pmon_1, whole genome shotgun sequence	chr_3
NC_051389.1 Penaeus monodon isolate SGIC_2016 chromosome 4, NSTDA_Pmon_1, whole genome shotgun sequence	chr_4
NC_051390.1 Penaeus monodon isolate SGIC_2016 chromosome 5, NSTDA_Pmon_1, whole genome shotgun sequence	chr_5
```
To upload both, the assembly and the key file, to the cluster computer, I used:
```
rsync -av --progress /Users/angfu103/Dropbox/Collaborations/Pmonodon_JCA/data/genome/genome_assemblies_genome_fasta/ncbi-genomes-2022-05-14/GCF_015228065.1_NSTDA_Pmon_1_genomic.fna angela@rackham.uppmax.uu.se:/home/angela/Pmon/

rsync -av --progress /Users/angfu103/Dropbox/Collaborations/Pmonodon_JCA/data/genome/genome_assemblies_genome_fasta/ncbi-genomes-2022-05-14/scaffold_names.txt angela@rackham.uppmax.uu.se:/home/angela/Pmon/
```
Info on the SeqKit function:

replace name/sequence by regular expression.

Note that the replacement supports capture variables.
e.g. $1 represents the text of the first submatch.
ATTENTION: use SINGLE quote NOT double quotes in \*nix OS.

Examples: Adding space to all bases.
    seqkit replace -p '(.)' -r '$1 ' -s
Or use the \ escape character.
    seqkit replace -p '(.)' -r '\$1 ' -s
more on: http://bioinf.shenwei.me/seqkit/usage/#replace

Special replacement symbols (only for replacing name not sequence):
    {nr}    Record number, starting from 1
    {kv}    Corresponding value of the key (captured variable $n) by key-value file,
            n can be specified by flag -I (--key-capt-idx) (default: 1)

Usage:
  seqkit replace [flags]

Flags:
  -s, --by-seq                 replace seq
  -h, --help                   help for replace
  -i, --ignore-case            ignore case
  -K, --keep-key               keep the key as value when no value found for the key (only for sequence name)
  -I, --key-capt-idx int       capture variable index of key (1-based) (default 1)
  -m, --key-miss-repl string   replacement for key with no corresponding value
  -k, --kv-file string         tab-delimited key-value file for replacing key with value when using "{kv}" in -r (--replacement) (only for sequence name)
      --nr-width int           minimum width for {nr} in flag -r/--replacement. e.g., formating "1" to "001" by --nr-width 3 (default 1)
  -p, --pattern string         search regular expression
  -r, --replacement string     replacement. supporting capture variables.  e.g. $1 represents the text of the first submatch. ATTENTION: for \*nix OS, use SINGLE quote NOT double quotes or use the \ escape character. Record number is also supported by "{nr}".use ${1} instead of $1 when {kv} given!

Before making the header replacement on the actual fasta file, I will test the commands with toy files:
```
# Make toy files.
nano test.fa
>HiC_scaffold_10
CCCCAAAACCCCATGATCATGGATC
>HiC_scaffold_9
CCCCAAAACCCCATGGCATCATTCA
>HiC_scaffold_4
CCCCAAAACCCCATGTTGCTACTAG
>HiC_scaffold_11
CCCCAAAACCCCATGATCATGGATC
>HiC_scaffold_12
CCCCAAAACCCCATGGCATCATTCA

nano keyValue.txt
HiC_scaffold_10	chr_1
HiC_scaffold_9	chr_2
HiC_scaffold_4	chr_3
HiC_scaffold_11	chr_4
HiC_scaffold_12	chr_5

# Command to replace headers.
seqkit replace -p '(.+)$' -r '{kv}' -k keyValue.txt test.fa > new.test.fa

# Compare before and after files.
cat test.fa
>HiC_scaffold_10
CCCCAAAACCCCATGATCATGGATC
>HiC_scaffold_9
CCCCAAAACCCCATGGCATCATTCA
>HiC_scaffold_4
CCCCAAAACCCCATGTTGCTACTAG
>HiC_scaffold_11
CCCCAAAACCCCATGATCATGGATC
>HiC_scaffold_12
CCCCAAAACCCCATGGCATCATTCA

cat new.test.fa
>chr_1
CCCCAAAACCCCATGATCATGGATC
>chr_2
CCCCAAAACCCCATGGCATCATTCA
>chr_3
CCCCAAAACCCCATGTTGCTACTAG
>chr_4
CCCCAAAACCCCATGATCATGGATC
>chr_5
CCCCAAAACCCCATGGCATCATTCA

```
All good. Now, run the package on the actual fasta file:
```bash
# Load required packages
module load bioinfo-tools
module load SeqKit/0.15.0

# Command to rename the contigs
seqkit replace -p '(.+)$' -r '{kv}' -k scaffold_names.txt GCF_015228065.1_NSTDA_Pmon_1_genomic.fna > tmp.fasta

# Stdout
[INFO] read key-value file: scaffold_names.txt
[INFO] 26875 pairs of key-value loaded
```
Compare the oriiginal and new headers of contigs:
```
cat GCF_015228065.1_NSTDA_Pmon_1_genomic.fna | grep ">" | head
>NC_051386.1 Penaeus monodon isolate SGIC_2016 chromosome 1, NSTDA_Pmon_1, whole genome shotgun sequence
>NC_051387.1 Penaeus monodon isolate SGIC_2016 chromosome 2, NSTDA_Pmon_1, whole genome shotgun sequence
>NC_051388.1 Penaeus monodon isolate SGIC_2016 chromosome 3, NSTDA_Pmon_1, whole genome shotgun sequence
>NC_051389.1 Penaeus monodon isolate SGIC_2016 chromosome 4, NSTDA_Pmon_1, whole genome shotgun sequence
>NC_051390.1 Penaeus monodon isolate SGIC_2016 chromosome 5, NSTDA_Pmon_1, whole genome shotgun sequence
>NC_051391.1 Penaeus monodon isolate SGIC_2016 chromosome 6, NSTDA_Pmon_1, whole genome shotgun sequence
>NC_051392.1 Penaeus monodon isolate SGIC_2016 chromosome 7, NSTDA_Pmon_1, whole genome shotgun sequence
>NC_051393.1 Penaeus monodon isolate SGIC_2016 chromosome 8, NSTDA_Pmon_1, whole genome shotgun sequence
>NC_051394.1 Penaeus monodon isolate SGIC_2016 chromosome 9, NSTDA_Pmon_1, whole genome shotgun sequence
>NC_051395.1 Penaeus monodon isolate SGIC_2016 chromosome 10, NSTDA_Pmon_1, whole genome shotgun sequence

cat tmp.fasta | grep ">" | head
>chr_1
>chr_2
>chr_3
>chr_4
>chr_5
>chr_6
>chr_7
>chr_8
>chr_9
>chr_10

```
All good! Now give a proper name to the tmp file:
```
mv tmp.fasta GCF_015228065.1_NSTDA_Pmon_1_genomic_renamed.fna
```
Download file to laptop:
```
rsync -av --progress angela@rackham.uppmax.uu.se:/home/angela/Pmon/GCF_015228065.1_NSTDA_Pmon_1_genomic_renamed.fna /Users/angfu103/Dropbox/Collaborations/Pmonodon_JCA/data/genome/genome_assemblies_genome_fasta/ncbi-genomes-2022-05-14/
```

## Obtain index files for the reference genome
For read mapping and variant calling using BWA and GATK, respectively, it is necessary to generate a group of indexes and dictionaries for the reference genome. For this, I used the bash script called `00-prepare_ref_genome_for_BWA_GATK.sh`:

```bash
#!/bin/bash

# Load required software.
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9
module load picard/2.20.4

# Set environment variables.
REF_FILE='/path/00-genome/GCF_015228065.1_NSTDA_Pmon_1_genomic_renamed.fna'
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
java -Xmx48g -jar $PICARD_JAR CreateSequenceDictionary \
REFERENCE=${REF_FILE} \
OUTPUT=${REF_FILE}.dict

echo "################# picard dictionary done" ;

```
**Runtime: 00-00:25:52**

Rename genome dictionary file, because GATK will expect it to be named as `genome.dict` not `genome.fasta.dict`:
```
cd /path/00-genome
mv A.obtectus_v2.0.fasta.dict A.obtectus_v2.0.dict
```
