# 2022-04-06

## Compare sequence quality statistics of pools before and after read trimming

Assuming that a quality control of raw reads was already performed (i.e. Illumina adapters and low quality bases were removed with a tool like [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)) and that [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reports were generated for the raw and cleaned reads, let's compare the overall sequence quality statistics of the pools before and after quality control (QC). The aim of this step is to examine whether all pools "behave similarly" and there is no evidence of technical artifacts.

For this we will use [MultiQC](https://multiqc.info). Using the Terminal window, go to the directory where you have the results of FastQC of the **raw reads** and run multiqc from there:
```
cd /path/data/raw_sequences/fastQC_results

multiqc .
```
The program automatically identifies and merges the FastQC reports and should print in screen the progress of the analysis:
```
[WARNING]         multiqc : MultiQC Version v1.8 now available!
[INFO   ]         multiqc : This is MultiQC v1.7
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Searching '.'
[INFO   ]          fastqc : Found 24 reports
[INFO   ]         multiqc : Compressing plot data
[INFO   ]         multiqc : Report      : multiqc_report.html
[INFO   ]         multiqc : Data        : multiqc_data
[INFO   ]         multiqc : MultiQC complete
```
Now, do the same but for the FastQC reports of the **clean reads**:
```
cd /path/data/clean_sequences/fastQC_results

multiqc .
```
Once this is done, please share by email the html files generated in both runs, and we will discuss the results.

Good luck!
