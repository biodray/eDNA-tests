# My first eDNA pipeline

This code purpose is to analysis eDNA data. It was build for the analysis of freshwater fish, at two loci, 12S and cytB. This code run in R, but need external programs such as FastQC (), cutadapt (), usearch () and vsearch ().


## Folder

- **00_Data** for raw and proceeded data
- **01_Results** for all results 
- **02_Codes** central codes
- **03_Functions** useful functions necessary in more than one code
- **04_log** log created by running codes


## Pipelines

### Before starting an analysis

- Put raw read files in 00_Data/01a_RawData
- Put metadata in 00_Data/00_FileInfos
- Use **Change_PARAM.R** to modify repertories and files names
- Use **Rename_RAW.R** to unzip fastq files and rename them
- Use **Create_REF.R** to create a reference dataset for species assignment

### Pipelines

Run sequentially these codes:

1. **Process_RAW.R** to transform raw reads into ASV and OTU tables. Must run in UNIX.
2. **Correct_SEQTAB.R** to correct read numbers based on negative samples
3. **Assign_SP.R** to assign taxonomy to ASV/OTU tables
4. **Compare_RESULTS.R** to create graphes and models based on ASV/OTU tables

## Other stuffs

- Use **Blast.R** to use NCBI blast program on a locale reference dataset. Useful to double-check species assignation.

## Stuff to do eventually

- Enable doing more analysis under windows (enable local dir for external program)
- Remove yes/no question in Process_RAW.R
- Add more info in Process_RAW.R log file
- Add more log files (within each code)
- Create local metadata, in a more universal version
- Create a RMD document compiling raw infos + majors results



