# 2022-04-06

## Explore sequence quality before and after quality control (QC)

Given that you already cleaned the raw reads to remove Illumina adapters and low quality bases, let's compare the overall quality statistics between pools before and after QC. The aim of this is to examine whether all pools "behave similarly" and there is no evidence of technical artifacts.

For this we will use [MultiQC](https://multiqc.info). Using the Terminal window, go to the directory where you have the results of FastQC of the raw reads and run multiqc from there:
```
cd /path/data/raw_sequences/fastQC_results

multiqc .
```
The program automatically identifies and merges the FastQC results and should print in screen the progress of the analysis:
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
Now, do the same but for the FastQC reports of the clean reads:
```
cd /path/data/clean_sequences/fastQC_results

multiqc .
```
Once this is done, please share by email the html files generated in both runs, and we will discuss the results.

Good luck!
