## Reference genome of the black tiger shrimp (Penaeus monodon)

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
