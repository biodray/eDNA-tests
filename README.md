# My first eDNA pipeline

This code purpose is to analysis eDNA data. It was build for the analysis of freshwater fish, at two loci, 12S and cytB. This code run in R, but need external program such as FastQC (), cutadapt (), usearch () and vsearch ().


## Folder

**00_Data** for raw and proceeded data

**01_Results** for all results 

**02_Codes** central codes

**03_Functions** useful functions necessary in more than one code

**04_log** log created by running codes


## Pipelines

### Before starting an analysis

- Put raw read files in 00_Data/01a_RawData
- Put metadata in 00_Data/00_FileInfos

**Change_PARAM.R** to modify repertories and files names

**Rename_RAW.R** to unzip fastq files and rename them

**Create_REF.R** to create a reference dataset for species assignment

### Pipelines

**Process_RAW.R** to transform raw reads into ASV and OTU tables. Must run in UNIX.

**Correct_SEQTAB.R** to correct read numbers based on negative samples

**Assign_SP.R** to assign taxonomy to ASV/OTU tables

**Compare_RESULTS.R** to create graphes and models based on ASV/OTU tables

## Other stuffs

**Blast.R** to use NCBI blast program on a locale reference dataset. Useful to double-check species assignation.

## Stuff to do eventually

- enable doing more analysis under windows (enable local dir for external program)
- Remove yes/no question in Process_RAW.R
- Add more info in Process_RAW.R log file
- Add more log files (within each code)
- Create local metadata, in a more universal version



