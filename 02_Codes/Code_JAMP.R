
# Info --------------------------------------------------------------------






# Library -----------------------------------------------------------------


library(JAMP)

# JAMP generates a new folder for each processing step! Should your files alreay be demultiplexed or you want to start somwhere else in the pipeline with preprocessed reads, you can generate an empty folder and place your files to be processed in "_data"
Empty_folder()

#Remove_last_folder()
files.names <- list.files("./A_Empty_Folder/_data/")


# FastQC

FastQC(exe = "/home/biodray/bin/fastqc")

# Umerge
getwd()

setwd("./eDNA-JAMP")

test <- list.files("./A_Empty_Folder/_data/", pattern = "12s", full.names=T)

U_merge_PE(files = test, fastq_pctid=75)



# trimm primers
Cutadapt(forward="ACTGGGATTAGATACCCC", # mlCOIintF
         reverse="TAGAACAGGCTCCTCTAG", LDist=T) # jgHCO2198


# discard with non target length
Minmax(min=(104-10), max=(104+10))


# SRA = function to download directly from SRA

# discard reads above 1 expected error
U_max_ee(max_ee=1) # usearch


# creating otus
U_cluster_otus(filter=0.01)

