
# Info --------------------------------------------------------------------

# Run the DADA2 and JAMP pipelines in parallele, on eDNA data
# 
# Audrey Bourret
# 2018-10-22
#

start.time <- Sys.time()

# Library -----------------------------------------------------------------

library(gtools)    # for mixedsort
library(stringr)   # for pattern detection

library(tidyverse) # for data manipulation
library(readxl)

library(ggpubr) # on github - for nice graphs
#devtools::install_github("kassambara/ggpubr")

library(dada2); packageVersion("dada2") # Faire mettre cette info dans le log
library(JAMP)

#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/fastqcr")

library(fastqcr)

library(Biostrings)

# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}


# Data --------------------------------------------------------------------

get.value("info.path")

get.value("raw.path")  
get.value("raw_unz.path")   
get.value("raw_unz_rename.path")

get.value("ref.path")

filt_dada2.path <- "./00_Data/02a_Filtered_dada2"
filt_IBIS.path <- "./00_Data/02b_Filtered_IBIS"
filt_JAMP.path  <- "./00_Data/02c_Filtered_JAMP"

get.value("log.path")

get.value("result.path")

DataSample <- read_excel(get.value("Sample.xl"),sheet="DataSample",na="NA",guess_max=100000)
DataSeq    <- read_excel(get.value("Sample.xl"),sheet="DataSeq",na="NA",guess_max=100000)
Amorces    <- read_excel(get.value("Sample.xl"),sheet="Amorces",na="NA",guess_max=100000)
#Inventaire <- read_excel(get.value("Sample.xl"),sheet="DataLac",na="NA",guess_max=100000)

# Format de DataSeq avec les données IBIS (p1-A1)

DataSeq <- DataSeq %>% mutate (IbisID = paste0("p",Plaque,"-",Puit)) %>% 
  left_join(DataSample %>% select(SampleID, NomLac, CatSite), by = "SampleID")


# Log

if(file.exists(file.path(log.path, "Process_RAW.log.txt"))){
  cat("\nWARNING: The log file Process_RAW.log.txt will be erased\n")
}

cat("\n-------------------------\n", 
    "Process raw eDNA data\n",
    date(),
    "\n-------------------------\n", 
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = FALSE, sep = "\n")


# Function to extract file list

# list.raw.files <- function(LOCI, 
#                            PATH, 
#                            FRpattern = c("R1", "R2"), 
#                            STARTpattern = "EP-",
#                            ENDpattern = "_L001_R1_001.fastq"){
#   
#   # Create a list to return all values 
#   res <- list()
#   
#   for(x in LOCI){
#     # List all files based on loci name
#     files.all <- list.files(PATH, pattern = x)
#     
#     # Divide F and R reads
#     files.R1 <- mixedsort(files.all[stringr::str_detect(files.all, pattern = FRpattern[1])])
#     files.R2 <- mixedsort(files.all[stringr::str_detect(files.all, pattern = FRpattern[2])])
# 
#     #
#     names <- files.R1 %>% stringr::str_remove(pattern = STARTpattern) %>% 
#                           stringr::str_remove(pattern = ENDpattern) %>% 
#                           mixedsort()
#     
#     names <- sapply(stringr::strsplit(names, "_"), `[`, 1)
#     
#     # list names
#     files.R1.ln <- paste0(x, ".raw.files.R1")
#     files.R2.ln <- paste0(x, ".raw.files.R2")   
#     names.ln    <- paste0(x, ".names")
#     
#     # Put restults in my list
#     res[[files.R1.ln]] <- files.R1
#     res[[files.R2.ln]] <- files.R2
#     res[[names.ln]]    <- names
#     
#     
#     
#     if(!is.null(log.path)){
#       cat(paste0(x, ": ", length(names), " files found in ", PATH),  
#         file=file.path(log.path, "Process_RAW.log.txt"), 
#         append = T, sep = "\n")
#       
#     }
# 
#     
#   }
#   
#   if(!is.null(log.path)){
#     cat("\n-------------------------\n",  
#         file=file.path(log.path, "Process_RAW.log.txt"), 
#         append = T, sep = "\n")
#     
#   }
#   
#   return(res)
#   
# }
# 
# # add filt path 
# 
# add.filt.files <- function(LOCI, 
#                            PATH,
#                            FILE.LS,
#                            F.PATTERN = "_F_filt.fastq.gz", 
#                            R.PATTERN = "_R_filt.fastq.gz"){
#   
#   # Create a list to return all values 
#   for(x in LOCI){
#     
#     NAMES.LN <-  paste0(x, ".names")
#     NAMES <- FILE.LS[[NAMES.LN]]
#     
#     # Divide F and R reads
#     filt.Fs <- file.path(PATH, paste0(NAMES, F.PATTERN))
#     filt.Rs <- file.path(PATH, paste0(NAMES, R.PATTERN))    
#     
#     # list names
#     files.R1.ln <- paste0(x, ".filt.files.R1")
#     files.R2.ln <- paste0(x, ".filt.files.R2")   
#     
#     # Put restults in my list
#     FILE.LS[[files.R1.ln]] <- filt.Fs
#     FILE.LS[[files.R2.ln]] <- filt.Rs
#   }
#   
#   return(FILE.LS)
#   
# }
# 
# add.filt.OK.files <- function(LOCI, 
#                               PATH,
#                               FILE.LS,
#                               F.PATTERN = "_F_filt.fastq.gz", 
#                               R.PATTERN = "_R_filt.fastq.gz"){
#   
#   # Create a list to return all values 
#   for(x in LOCI){
#     
#     # List all files based on loci name
#     files.all <- list.files(PATH, pattern = x)
#     
#     # Divide F and R reads
#     files.R1 <- mixedsort(files.all[stringr::str_detect(files.all, pattern = F.PATTERN)])
#     files.R2 <- mixedsort(files.all[stringr::str_detect(files.all, pattern = R.PATTERN)])
#     
#     names <- files.R1 %>% stringr::str_remove(pattern =  F.PATTERN)
#     
#     # list names
#     files.R1.ln <- paste0(x, ".filt.files.R1.OK")
#     files.R2.ln <- paste0(x, ".filt.files.R2.OK")
#     
#     names.ln <- paste0(x, ".filt.names")
#     
#     # Put results in my list
#     FILE.LS[[files.R1.ln]] <- files.R1
#     FILE.LS[[files.R2.ln]] <- files.R2
#     FILE.LS[[names.ln]]    <- names    
#   }
#   
#   return(FILE.LS)
#   
# }


# L'arranger pour enlever le merger

get_trackDADA <- function(SUMMARY, MERGER=NULL, SEQTAB, RAW.NAME, FILT.NAMES){
  #Fonction dans la fonction
  getN <- function(x) sum(getUniques(x))
  
  SUMMARY.DF <- as.data.frame(SUMMARY)
  SUMMARY.DF$sample <- RAW.NAME
  
  if(is.null(MERGER)){
    DADA.DF <- as.data.frame(cbind(rowSums(SEQTAB), FILT.NAMES))
    names(DADA.DF) <- c("nonchim", "sample")
    
    res <- SUMMARY.DF %>% full_join(DADA.DF, by = "sample") %>% 
      select(sample, nreads = reads.in, filtered = reads.out, nonchim) 
    
    
  } else {
    DADA.DF <- as.data.frame(cbind(sapply(MERGER, getN), rowSums(SEQTAB), FILT.NAMES))
    names(DADA.DF) <- c("merged", "nonchim", "sample")    
    
    res <- SUMMARY.DF %>% full_join(DADA.DF, by = "sample") %>% 
      select(sample, nreads = reads.in, filtered = reads.out, merged, nonchim) 
    
  }
  
  if(!is.null(MERGER)) {
    res$merged    <- as.numeric(as.character(res$merged))
    res$merged[which(is.na(res$merged))] <- 0
  }
  
  res$nonchim   <- as.numeric(as.character(res$nonchim))
  res$nonchim[which(is.na(res$nonchim))] <- 0  
  
  #res$denoisedF <- as.numeric(as.character(res$denoisedF))
  #res$denoisedR <- as.numeric(as.character(res$denoisedR))
  #res$denoisedF[which(is.na(res$denoisedF))] <- 0
  #res$denoisedR[which(is.na(res$denoisedR))] <- 0
  
  return(res)
  
}

makeSeqTabFromScratch <- function(files, name, path=""){
  files.seq <- list()
  for(x in 1:length(files)){
    #FILE <- path.files[str_detect(path.files, paste0(x,"_"))]
    DNA  <- readDNAStringSet(file.path(path,files[x]))
    #N   <- str_sub(sapply(str_split(names(DNA), "Otu"), `[`, 2),8)
    SEQ <- as.character(DNA)
    SEQ <- stringr::str_replace_all(SEQ, "-", "")
    if(length(SEQ)>0){
      SEQ.df <- as.data.frame(base::table(SEQ))
      names(SEQ.df) <- c("sequence","abundance")
      SEQ.df$abundance <- 1
      files.seq[[name[x]]] <- SEQ.df    
    } else {
      print(paste0("No sequences for sample ",name[x]))
    }
  }
  seqtab <- makeSequenceTable(files.seq, orderBy = "abundance")
  return(seqtab)
}

# Quality assesment - RAW ------------------------------------------------------

#list.raw.files <- function(LOCI, 
#                           PATH, 
#                           FRpattern = c("R1", "R2"), 
#                           STARTpattern = "EP-",
#                           ENDpattern = "_L001_R1_001.fastq")

#all.files <- list.raw.files(LOCI = c("12s", "cytB"), 
#                            PATH = get.value("raw_unz_rename.path"),
#                            FRpattern = c("R1", "R2"), 
#                            STARTpattern = "",
#                            ENDpattern = "R1.fastq")

#str(all.files)



plotQplus <- function(files, locus, pattern) {
  graph.ls <- list()
  
  for(l in locus){
    print(l)
    
    files.red <- files %>% str_subset(l)
    
    for(p in pattern){
      
      files.red2 <- files.red %>% str_subset(p)      
      print(p)
      
      graph <- plotQualityProfile(files.red2, aggregate = T)
      
      graph <- graph + labs(x = "Base pair")
      
      graph.ls[[paste(l,p,sep="-")]] <- graph 
    }
  }
  return(graph.ls)  
}   


if(file.exists(get.value("Qplot.RAW.data"))){
  load(get.value("Qplot.RAW.data"))
  
} else {
  graph.tot.ls   <- plotQplus(list.files(get.value("raw_unz_rename.path"), full.names = T),
                              locus = c("12s", "cytB"), pattern = c("R1", "R2"))
  
  graph.sample.ls <- plotQplus(list.files(get.value("raw_unz_rename.path"), pattern = "Sample", full.names = T),
                               locus = c("12s", "cytB"), pattern = c("R1", "R2"))
  
  save(file = get.value("Qplot.RAW.data"), 
       list = c("graph.tot.ls" , "graph.sample.ls"))
}

# devrait peut-être aller dans un autre script
names(graph.tot.ls)
names(graph.sample.ls)

graphQ.total <- ggarrange(graph.tot.ls[[1]] + labs(title = "12S  - Forward"),
                          graph.tot.ls[[2]] + labs(title = "12S  - Reverse"),
                          graph.tot.ls[[3]] + labs(title = "CYTB - Forward"),
                          graph.tot.ls[[4]] + labs(title = "CYTB - Reverse"),
                          labels = LETTERS[1:4],
                          ncol = 2, nrow = 2)

graphQ.total

graphQ.sample<- ggarrange(graph.sample.ls[[1]] + labs(title = "12S  - Forward"),
                          graph.sample.ls[[2]] + labs(title = "12S  - Reverse"),
                          graph.sample.ls[[3]] + labs(title = "CYTB - Forward"),
                          graph.sample.ls[[4]] + labs(title = "CYTB - Reverse"),
                          labels = LETTERS[1:4],
                          ncol = 2, nrow = 2)

graphQ.sample

# Save graphs

ggsave(filename = "QualityPlotTotal.pdf", 
       path = get.value("result.Q,path"),       
       plot = graphQ.total,
       device = "pdf",
       width = 8, height = 8, units = "in")

ggsave(filename = "QualityPlotTotal.png", 
       path = get.value("result.Q.path"),
       plot = graphQ.total,
       device = "png",
       width = 8, height = 8, units = "in")

ggsave(filename = "QualityPlotSample.pdf", 
       path = get.value("result.Q.path"),       
       plot = graphQ.sample,
       device = "pdf",
       width = 8, height = 8, units = "in")

ggsave(filename = "QualityPlotSample.png", 
       path = get.value("result.Q.path"),
       plot = graphQ.sample,
       device = "png",
       width = 8, height = 8, units = "in")


cat("Graphics done!", "\n-------------------------\n",  
    file=get.value("RAW.log"), 
    append = T, sep = "\n")

# FastQC

if(get_os() %in% c("os","linux")){ # to run only on os and not windows
  
  cmd <- paste("--outdir", get.value("result.FQraw.path"), list.files(get.value("raw_unz_rename.path"), full.names = T))
  system2("fastqc", cmd)
  
} else {cat("Cannot perform FastQC on windows -- sorry!!")}

##from JAMP
#FastQC(files = list.files(get.value("raw_unz_rename.path"), full.names = T), exe = "fastqc")

# from fastqc
#fastqc_install()
#fastqc(fq.dir = get.value("result.FQraw.path"), 
#       qc.dir = get.value("raw_unz_rename.path"),
#       threads = 4)

qc.12s.res  <- qc_JAMP_sum(list.files(get.value("result.FQraw.path") %>% mixedsort(), pattern = "fastqc.zip", full.names = T) %>% str_subset(pattern = "12s"))
qc.cytB.res <- qc_JAMP_sum(list.files(get.value("result.FQraw.path") %>% mixedsort(), pattern = "fastqc.zip", full.names = T) %>% str_subset(pattern = "cytB"))

names(qc.12s.res)
names(qc.cytB.res)

qc_JAMP_plot(qc.12s.res$summary, file = file.path(get.value("result.Q.path"),"FastQC_12s_summary.pdf" ))
qc_JAMP_plot(qc.cytB.res$summary, file = file.path(get.value("result.Q.path"),"FastQC_cytB_summary.pdf" ))


save(file = get.value("FastQC.data"), 
     list = c("qc.12s.res" , "qc.cytB.res"))

# Write the result somewhere
# write.csv(exp, paste(folder, "/FastQC/stats.csv", sep=""))


# RAW to FILT (cutadapt + dada2) ------------------------------------------------------------

# Retry to add cutadapt

file1 <- list.files(get.value("raw_unz_rename.path"), pattern = "12s", full.names=T) %>% str_subset("R1")
file2 <- file1 %>% str_replace("R1","R2")

new.file1 <- file1 %>% str_replace("01c_RawData_unzipped_rename", "02a_Cutadapt") %>% str_replace (".fastq", "_cut.fastq")
new.file2 <- file2 %>% str_replace("01c_RawData_unzipped_rename", "02a_Cutadapt") %>% str_replace (".fastq", "_cut.fastq")

for(x in 1:length(file1)){
  
  cat("\n",file1[x], sep="\n")
  
  cutadapt.12S.pairend <- paste("-g ^ACTGGGATTAGATACCCC",
                                "-G ^TAGAACAGGCTCCTCTAG",
                                "-o", new.file1[x], 
                                "--paired-output", new.file2[x], 
                                file1[x],
                                file2[x],
                                "-f", "fastq",      
                                "--discard-untrimmed", 
                                "--report=minimal",
                                sep = " ") # forward adapter
  
  system2("cutadapt", cutadapt.12S.pairend, stdout="", stderr="") # to R console
  
}
# CYTB

file1 <- list.files(get.value("raw_unz_rename.path"), pattern = "cytB", full.names=T) %>% str_subset("R1")
file2 <- file1 %>% str_replace("R1","R2")

new.file1 <- file1 %>% str_replace("01c_RawData_unzipped_rename", "02a_Cutadapt") %>% str_replace (".fastq", "_cut.fastq")
new.file2 <- file2 %>% str_replace("01c_RawData_unzipped_rename", "02a_Cutadapt") %>% str_replace (".fastq", "_cut.fastq")

for(x in 1:length(file1)){
  
  cat("\n",file1[x], sep="\n")
  
  #F - CYTB
  cutadapt.cytB.F.cmd <- paste("-g ^AAAAAGCTTCCATCCAACATCTCAGCATGATGAAA",
                               #"--minimum-length", "75", #  After truncation
                               "-o", new.file1[x], 
                               file1[x], 
                               "-f", "fastq",
                               "--discard-untrimmed", 
                               "--report=minimal", 
                               sep = " ") # forward adapter
  
  system2("cutadapt", cutadapt.cytB.F.cmd, stdout="", stderr="") # to R consol
  
  #r - CYTB
  cutadapt.cytB.R.cmd <- paste("-g ^AAACTGCAGCCCCTCAGAATGATATTTGTCCTCA", # pas en reverse complement car encore du bon sens
                               #"--minimum-length", "75", #  After truncation
                               "-o", new.file2[x], 
                               file2[x], 
                               "-f", "fastq",
                               "--discard-untrimmed", 
                               "--report=minimal", 
                               sep = " ") # forward adapter
  
  system2("cutadapt", cutadapt.cytB.R.cmd, stdout="", stderr="") # to R consol
  
}


#all.files <-  add.filt.files(LOCI = c("12s", "cytB"), PATH = filt_dada2.path, FILE.LS = all.files) 

#str(all.files)

# Filtrage en soi
# Real filtering

# # function to add filt names
# filt.names <- function(files){
#   new.files <- files %>% str_replace(".fastq", "_filt.fastq.gz") %>%
#                          str_replace(get.value("raw_unz_rename.path"), get.value("filt_dada2.path"))
#   
#   new.files
# }

file1 <- list.files(get.value("filt_cutadapt.path"), pattern = "12s", full.names = T) %>% str_subset("R1")
file2 <- file1 %>% str_replace("R1", "R2")

new.file1 <- file1 %>% str_replace(get.value("filt_cutadapt.path"), get.value("filt_dada2.path")) %>%  str_replace(".fastq", "_filt.fastq")
new.file2 <- file2 %>% str_replace(get.value("filt_cutadapt.path"), get.value("filt_dada2.path")) %>%  str_replace(".fastq", "_filt.fastq")

# 12S

filter.12s.summary <- filterAndTrim(fwd = file1,
                                    filt = new.file1,
                                    rev = file2,
                                    filt.rev = new.file2,
                                    truncQ=10, # minimum Q score, 10 = 90% base call accuracy
                                    truncLen = c(0,0), # Taille min/max des reads
                                    trimLeft= c(0, 0), # Deja enlevé avec cutadapt, sinon c(18,18) 
                                    maxLen = c(Inf,Inf),
                                    minLen = c(50,50), # after trimming and truncation
                                    maxEE=c(1,1),
                                    #orient.fwd = c("ACTGG"), # debut de l'amorce F
                                    compress = FALSE, # Voir si c'est OK de compresser spour JAMP
                                    multithread=FALSE, # TRUE on linux
                                    verbose = TRUE) 


# CytB
file1 <- list.files(get.value("filt_cutadapt.path"), pattern = "cytB", full.names = T) %>% str_subset("R1")
file2 <- file1 %>% str_replace("R1", "R2")

new.file1 <- file1 %>% str_replace(get.value("filt_cutadapt.path"), get.value("filt_dada2.path")) %>%  str_replace(".fastq", "_filt.fastq")
new.file2 <- file2 %>% str_replace(get.value("filt_cutadapt.path"), get.value("filt_dada2.path")) %>%  str_replace(".fastq", "_filt.fastq")



filter.cytB.R1.summary <- filterAndTrim(fwd = file1,
                                        filt = new.file1,
                                        #rev = list.files(get.value("raw_unz_rename.path"), pattern = "cytB", full.names = T) %>% str_subset("R2"),
                                        #filt.rev = list.files(get.value("raw_unz_rename.path"), pattern = "cytB", full.names = T) %>% str_subset("R2") %>% filt.names(),
                                        truncQ = 6,
                                        truncLen = c(75), # Taille min/max des reads
                                        trimLeft= c(0), # Je vais les enlever avec cutadapt - sinon 35
                                        maxLen = c(Inf),
                                        minLen = c(0), # after trimming and truncation
                                        maxEE=c(2),   
                                        compress = FALSE,
                                        multithread=FALSE, # TRUE on linux
                                        verbose = TRUE) 

filter.cytB.R2.summary <- filterAndTrim(fwd = file2,
                                        filt = new.file2,
                                        #rev = list.files(get.value("raw_unz_rename.path"), pattern = "cytB", full.names = T) %>% str_subset("R2"),
                                        #filt.rev = list.files(get.value("raw_unz_rename.path"), pattern = "cytB", full.names = T) %>% str_subset("R2") %>% filt.names(),
                                        truncQ = 6,
                                        truncLen = c(75), # Taille min/max des reads
                                        trimLeft= c(0), # Je vais les enlever avec cutadapt - sinon 34
                                        maxLen = c(Inf),
                                        minLen = c(0), # after trimming and truncation
                                        maxEE=c(2),   
                                        compress = FALSE,
                                        multithread=FALSE, # TRUE on linux
                                        verbose = TRUE) 


#cat(paste0("\n Correlation between 12S and cytB reads after DADA2 filtration : ",
#    round(cor(filter.12s.summary[,"reads.out"], filter.cytB.summary[,"reads.out"], method = "spearman"),2)
#     ),
#    "\n-------------------------\n",
#    file=file.path(log.path, "Process_RAW.log.txt"), 
#    append = T, sep = "\n")


# Add files (beacause some were discard)

#all.files <- add.filt.OK.files(LOCI = c("12s","cytB"), 
#                              PATH = filt_dada2.path,
#                              FILE.LS = all.files,
#                              F.PATTERN = "_F_filt.fastq.gz", 
#                              R.PATTERN = "_R_filt.fastq.gz")

# qUALITY ASSEMENT

if(file.exists(get.value("Qplot.FILT.data"))){
  load(get.value("Qplot.FILT.data"))
  
} else {
  graph.tot.filt.ls   <- plotQplus(list.files(get.value("filt_dada2.path"), full.names = T),
                                   locus = c("12s", "cytB"), pattern = c("R1", "R2"))
  
  graph.sample.filt.ls <- plotQplus(list.files(get.value("filt_dada2.path"), pattern = "Sample", full.names = T),
                                    locus = c("12s", "cytB"), pattern = c("R1", "R2"))
  
  save(file = get.value("Qplot.FILT.data"), 
       list = c("graph.tot.filt.ls" , "graph.sample.filt.ls"))
}

# devrait peut-être aller dans un autre script
names(graph.tot.filt.ls)
names(graph.sample.filt.ls)

graphQ.total.filt <- ggarrange(graph.tot.filt.ls[[1]] + labs(title = "12S  - Forward"),
                               graph.tot.filt.ls[[2]] + labs(title = "12S  - Reverse"),
                               graph.tot.filt.ls[[3]] + labs(title = "CYTB - Forward"),
                               graph.tot.filt.ls[[4]] + labs(title = "CYTB - Reverse"),
                               labels = LETTERS[1:4],
                               ncol = 2, nrow = 2)

graphQ.total.filt

graphQ.sample.filt<- ggarrange(graph.sample.filt.ls[[1]] + labs(title = "12S  - Forward"),
                               graph.sample.filt.ls[[2]] + labs(title = "12S  - Reverse"),
                               graph.sample.filt.ls[[3]] + labs(title = "CYTB - Forward"),
                               graph.sample.filt.ls[[4]] + labs(title = "CYTB - Reverse"),
                               labels = LETTERS[1:4],
                               ncol = 2, nrow = 2)

graphQ.sample.filt

# Save graphs

ggsave(filename = "QualityPlotTotal.filt.pdf", 
       path = get.value("result.Q.path"),       
       plot = graphQ.total.filt,
       device = "pdf",
       width = 8, height = 8, units = "in")

ggsave(filename = "QualityPlotTotal.filt.png", 
       path = get.value("result.Q.path"),
       plot = graphQ.total.filt,
       device = "png",
       width = 8, height = 8, units = "in")

ggsave(filename = "QualityPlotSample.filt.pdf", 
       path = get.value("result.Q.path"),       
       plot = graphQ.sample.filt,
       device = "pdf",
       width = 8, height = 8, units = "in")

ggsave(filename = "QualityPlotSample.filt.png", 
       path = get.value("result.Q.path"),
       plot = graphQ.sample.filt,
       device = "png",
       width = 8, height = 8, units = "in")


# FILT to ASV (denoise with DADA2) --------------------------------------------------

# Calcul du taux d'erreur (DADA2)

err.12s.F <- learnErrors(list.files(get.value("filt_dada2.path"), pattern = "12s", full.names = T) %>% str_subset("R1")) 
err.12s.R <- learnErrors(list.files(get.value("filt_dada2.path"), pattern = "12s", full.names = T) %>% str_subset("R2")) 


#pdf("01_Results/ErrorsRate.dada2.12S.pdf") 
#  plotErrors(err.12s.F, nominalQ=TRUE)
#  plotErrors(err.12s.R, nominalQ=TRUE)
#dev.off()

#msg1 <- "12S: DADA2 error rate plot saved: 01_Results/ErrorsRate.dada2.12S.pdf"  

err.cytB.F <- learnErrors(list.files(get.value("filt_dada2.path"), pattern = "cytB", full.names = T) %>% str_subset("R1")) 
err.cytB.R <- learnErrors(list.files(get.value("filt_dada2.path"), pattern = "cytB", full.names = T) %>% str_subset("R2")) 

#pdf("01_Results/ErrorsRate.dada2.CYTB.pdf") 
#  plotErrors(err.cytB.F, nominalQ=TRUE)
#  plotErrors(err.cytB.R, nominalQ=TRUE)
#dev.off()

#msg1 <- "cytB: DADA2 error rate plot saved: 01_Results/ErrorsRate.dada2.CYTB.pdf"  

#cat(msg1, msg2, "\n-------------------------\n",  
#    file=file.path(log.path, "Process_RAW.log.txt"), 
#    append = T, sep = "\n")

# Déréplication

derep.12s.F <- derepFastq(list.files(get.value("filt_dada2.path"), pattern = "12s", full.names = T) %>% str_subset("R1")) 
#names(derep.12s.F) <- all.files[["12s.filt.names"]] 

derep.12s.R <- derepFastq(list.files(get.value("filt_dada2.path"), pattern = "12s", full.names = T) %>% str_subset("R2")) 
#names(derep.12s.R) <- all.files[["12s.filt.names"]]

derep.cytB.F <- derepFastq(list.files(get.value("filt_dada2.path"), pattern = "cytB", full.names = T) %>% str_subset("R1")) 
#names(derep.cytB.F) <- all.files[["cytB.filt.names"]]

derep.cytB.R <- derepFastq(list.files(get.value("filt_dada2.path"), pattern = "cytB", full.names = T) %>% str_subset("R2")) 
#names(derep.cytB.R) <- all.files[["cytB.filt.names"]]


# Inférence des échantillons

dada.12s.F <- dada(derep.12s.F, 
                   err = err.12s.F, 
                   multithread=FALSE,
                   pool=TRUE)
#HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)

dada.12s.R <- dada(derep.12s.R, 
                   err = err.12s.R, 
                   multithread=FALSE,
                   pool=TRUE)
#HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)


dada.cytB.F <- dada(derep.cytB.F, 
                    err = err.cytB.F, 
                    multithread=FALSE,
                    pool=TRUE)
#HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)

dada.cytB.R <- dada(derep.cytB.R, 
                    err = err.cytB.R, 
                    multithread=FALSE,
                    pool=TRUE)
#HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)



## DADA2 Options: Non-Illumnina sequencing technologies

mergers.12S <- mergePairs(dadaF = dada.12s.F, 
                          derepF = derep.12s.F, 
                          dadaR = dada.12s.R, 
                          derepR = derep.12s.R, 
                          minOverlap = 30, 
                          maxMismatch = 0,
                          returnRejects = TRUE,
                          verbose=TRUE)

#mergers.cytB <- mergePairs(dadaF = dada.cytB.Fs, 
#                          derepF = derep.cytB.Fs, 
#                          dadaR = dada.cytB.Rs, 
#                          derepR = derep.cytB.Rs, 
#                          maxMismatch = 0,
#                          returnRejects = TRUE,
#                          verbose=TRUE,
#                          justConcatenate=TRUE)

# Make sequence table

seqtab.12s.int    <- makeSequenceTable(mergers.12S)
#seqtab.12S.F.int  <- makeSequenceTable(dada.12s.Fs)
#seqtab.12S.R.int  <- makeSequenceTable(dada.12s.Rs)

#seqtab.cytB.int   <- makeSequenceTable(mergers.cytB)
seqtab.cytB.F.int <- makeSequenceTable(dada.cytB.F)
seqtab.cytB.R.int <- makeSequenceTable(dada.cytB.R)

# Remove chimera

seqtab.12s <- removeBimeraDenovo(seqtab.12s.int, method = "consensus", 
                                 multithread = FALSE, verbose = TRUE)

#seqtab.12s.F <- removeBimeraDenovo(seqtab.12S.F.int, method = "consensus", 
#                                 multithread = FALSE, verbose = TRUE)

#seqtab.12s.R <- removeBimeraDenovo(seqtab.12S.R.int, method = "consensus", 
#                                   multithread = FALSE, verbose = TRUE)

#seqtab.cytB <- removeBimeraDenovo(seqtab.cytB.int, method = "consensus", 
#                                  multithread = FALSE, verbose = TRUE)

seqtab.cytB.F <- removeBimeraDenovo(seqtab.cytB.F.int, method = "consensus", 
                                    multithread = FALSE, verbose = TRUE)

seqtab.cytB.R <- removeBimeraDenovo(seqtab.cytB.R.int, method = "consensus", 
                                    multithread = FALSE, verbose = TRUE)   


# Save them somewhere!!!

# Stats on chimera
sum(seqtab.12s)/sum(seqtab.12s.int)
#sum(seqtab.12s.F)/sum(seqtab.12S.F.int)
#sum(seqtab.12s.R)/sum(seqtab.12S.R.int)

#sum(seqtab.cytB)/sum(seqtab.cytB.int)
sum(seqtab.cytB.F)/sum(seqtab.cytB.F.int)
sum(seqtab.cytB.R)/sum(seqtab.cytB.R.int)

cat(paste0("12S: Percentage of ASV remaining after chimera removal: ", 
           round(sum(seqtab.12s)/sum(seqtab.12s.int)*100,1), "%"), 
    paste0("12S.F: Percentage of ASV remaining after chimera removal: ", 
           round(sum(seqtab.12s.F)/sum(seqtab.12S.F.int)*100,1), "%"), 
    paste0("12S.R: Percentage of ASV remaining after chimera removal: ", 
           round(sum(seqtab.12s.R)/sum(seqtab.12S.R.int)*100,1), "%"), 
    paste0("cytB: Percentage of ASV remaining after chimera removal: ", 
           round(sum(seqtab.cytB)/sum(seqtab.cytB.int)*100,1), "%"), 
    paste0("12S.F: Percentage of ASV remaining after chimera removal: ", 
           round(sum(seqtab.cytB.F)/sum(seqtab.cytB.F.int)*100,1), "%"), 
    paste0("12S.R: Percentage of ASV remaining after chimera removal: ", 
           round(sum(seqtab.cytB.R)/sum(seqtab.cytB.R.int)*100,1), "%"),
    "\n-------------------------\n",  
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = T, sep = "\n")

# Save Seqtab

save(file = file.path(result.path, "Seqtab.data"), 
     list = ls(pattern = "seqtab."))

# Get track

sum.12S <- get_trackDADA(SUMMARY = filter.12s.summary,
                         RAW.NAME = all.files[["12s.names"]] %>% str_remove(pattern = "12s-"),
                         #DADAf = dada.12s.Fs.seta, 
                         #DADAr = dada.12s.Rs.seta,
                         MERGER = mergers.12S,
                         SEQTAB = seqtab.12s,
                         FILT.NAMES = all.files[["12s.filt.names"]] %>% str_remove(pattern = "12s-")
)

sum.12S.F <- get_trackDADA(SUMMARY = filter.12s.summary,
                           RAW.NAME = all.files[["12s.names"]] %>% str_remove(pattern = "12s-"),
                           SEQTAB = seqtab.12s.F,
                           FILT.NAMES = all.files[["12s.filt.names"]] %>% str_remove(pattern = "12s-")
)

sum.12S.R <- get_trackDADA(SUMMARY = filter.12s.summary,
                           RAW.NAME = all.files[["12s.names"]] %>% str_remove(pattern = "12s-"),
                           SEQTAB = seqtab.12s.R,
                           FILT.NAMES = all.files[["12s.filt.names"]] %>% str_remove(pattern = "12s-")
)


sum.cytB <- get_trackDADA(SUMMARY = filter.cytB.summary,
                          RAW.NAME = all.files[["cytB.names"]] %>% str_remove(pattern = "cytB-"),
                          #DADAf = dada.12s.Fs.seta, 
                          #DADAr = dada.12s.Rs.seta,
                          MERGER = mergers.cytB,
                          SEQTAB = seqtab.cytB,
                          FILT.NAMES = all.files[["cytB.filt.names"]] %>% str_remove(pattern = "cytB-")
)

sum.cytB.F <- get_trackDADA(SUMMARY = filter.cytB.summary,
                            RAW.NAME = all.files[["cytB.names"]] %>% str_remove(pattern = "cytB-"),
                            SEQTAB = seqtab.cytB.F,
                            FILT.NAMES = all.files[["cytB.filt.names"]] %>% str_remove(pattern = "cytB-")
)
sum.cytB.R <- get_trackDADA(SUMMARY = filter.cytB.summary,
                            RAW.NAME = all.files[["cytB.names"]] %>% str_remove(pattern = "cytB-"),
                            SEQTAB = seqtab.cytB.R,
                            FILT.NAMES = all.files[["cytB.filt.names"]] %>% str_remove(pattern = "cytB-")
)


cat("\n 12S:  Process raw summary : \n",  
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = T, sep = "\n")


write.table(sum.12S, 
            file=file.path(log.path, "Process_RAW.log.txt"),  
            row.names=T, col.names=T, append = T)


cat("\n cytB:  Process raw summary : \n",  
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = T, sep = "\n")


write.table(sum.cytB, 
            file=file.path(log.path, "Process_RAW.log.txt"),  
            row.names=T, col.names=T, append = T)


plot(sum.cytB.F$nonchim, sum.cytB.R$nonchim)
plot(sum.12S.F$nonchim, sum.12S.R$nonchim)

plot(sum.12S.R$nonchim, sum.12S$nonchim)
plot(sum.12S.F$nonchim, sum.12S$nonchim)

plot(sum.cytB.F$nonchim, sum.cytB$nonchim)
plot(sum.cytB.R$nonchim, sum.cytB$nonchim)

sum.12S.wInfo  <- sum.12S  %>% left_join(DataSeq %>% select(SampleID, SeqType, IbisID, NomLac, CatSite), by = c("sample" = "IbisID"))
sum.cytB.wInfo <- sum.cytB %>% left_join(DataSeq %>% select(SampleID, SeqType, IbisID, NomLac, CatSite), by = c("sample" = "IbisID"))


sum.12S.wInfo.res <- sum.12S.wInfo %>%  filter(SeqType %in% c("sample", "mix", "blank")) %>% 
  group_by(SeqType) %>% 
  summarise(perc.all.filt = round(sum(filtered)/sum(nreads),3),
            perc.filt.merged = round(sum(merged)/sum(filtered),3),
            perc.merged.nochim = round(sum(nonchim)/sum(merged),3),
            perc.all.nochim  = round(sum(nonchim)/sum(nreads),3)
  ) %>% as.data.frame()


sum.cytB.wInfo.res <- sum.cytB.wInfo %>%  
  filter(SeqType %in% c("sample", "mix", "blank")) %>% 
  group_by(SeqType) %>% 
  summarise(perc.all.filt = round(sum(filtered)/sum(nreads),3),
            perc.filt.merged = round(sum(merged)/sum(filtered),3),
            perc.merged.nochim = round(sum(nonchim)/sum(merged),3),
            perc.all.nochim  = round(sum(nonchim)/sum(nreads),3)
  ) %>% 
  as.data.frame()

cat("\n 12S:  Process raw summary by SeqType: \n",  
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = T, sep = "\n")


write.table(sum.12S.wInfo.res, 
            file=file.path(log.path, "Process_RAW.log.txt"),  
            row.names=F, col.names=T, append = T)

cat("\n cytB:  Process raw summary by SeqType: \n",  
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = T, sep = "\n")


write.table(sum.cytB.wInfo.res, 
            file=file.path(log.path, "Process_RAW.log.txt"),  
            row.names=F, col.names=T, append = T)

cat("\n-------------------------\n",  
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = T, sep = "\n")



# IBIS: filtered reads (MOTHUR) to ASV ------------------------------------

# 
# 
# all.files[["IBIS.files"]] <- list.files(filt_IBIS.path, pattern = "stability.trim.contigs.good.unique.good.filter.unique.precluster.fn.0.06.rep.pick")
# 
# 
# all.files[["IBIS.files.names"]] <- all.files[["IBIS.files"]] %>% str_remove(pattern = "EP_") %>% 
#                                    str_remove(pattern = "_stability.trim.contigs.good.unique.good.filter.unique.precluster.fn.0.06.rep.pick.fasta") %>% 
#                                    stringr::str_replace_all("_", "-")
# 
# 
# seqtab.12s.IBIS <-  makeSeqTabFromScratch(files = all.files[["IBIS.files"]],
#                                           name = all.files[["IBIS.files.names"]],
#                                           path = filt_IBIS.path) 
# 
# save(file = file.path(result.path, "Seqtab.data"), 
#      list = ls(pattern = "seqtab."))


# FILT to OTU (usearch with JAMP)-------------------------------------------------------------

# Umerge

file1 <- list.files(get.value("filt_dada2.path"), pattern = "12s", full.names=T) %>% str_subset("R1")
file2 <- file1 %>% str_replace("R1", "R2")
new.file <- file1 %>% str_replace(get.value("filt_dada2.path"), get.value("filt_merge.usearch.path")) %>% 
                      str_remove("_R1") %>% 
                      str_replace(".fastq", "_merged.fastq")  

for(x in 1:length(file1)){
  
  print(file1[x])
  
  cmd <- paste("-fastq_mergepairs", file1[x], 
               "-reverse", file2[x],  
               "-fastqout", new.file[x], 
               "-report", new.file[x] %>% str_replace(get.value("filt_merge.usearch.path"),paste0(get.value("filt_merge.usearch.path"), "/log")) %>% str_replace("_merged.fastq", "_log.txt"),
               "-fastq_maxdiffs", "99", 
               "-fastq_pctid", "75", 
               "-fastq_trunctail",  "0",
               "-fastq_minovlen", "30", 
               sep=" ")
  
  system2("usearch", cmd, stdout=F, stderr=F)  
  
}

# Add a min - max because some with 35 pb!!


# Test
files <- c(list.files(get.value("filt_merge.path"), pattern = "12s", full.names = T),
           list.files(get.value("filt_dada2.path"), pattern = "cytB", full.names = T)
           )
new.files <-  files %>% str_replace(get.value("filt_merge.path"), get.value("filt_derep.path")) %>% 
                        str_replace(get.value("filt_dada2.path"), get.value("filt_derep.path")) %>% 
                        str_replace(".fastq","_derep.fasta")

#U_cluster_otus(files= files, filter=0.01)

# Doesn't work - I will do it by hand once again.

# Dereplicated files with vsearch

for(x in 1:length(files)){
  print(files[x])
  
  cmd <- paste("-derep_fulllength", files[x],
             "-output", new.files[x],
             "-sizeout",
             "-sizein",
             sep = " ")
  
  system2("vsearch", cmd, stdout = T, stderr = T)

  gc(verbose = F) # unload memory
    
}

# Create one big file

?paste

cmd1 <- paste(paste(list.files(get.value("filt_derep.path"), full.names = T, pattern = "12s"), collapse = " "), 
             ">",
              paste0(get.value("filt_derep.1file.path"), "/all.12s.fasta"),
             sep = " "
              )

cmd2 <- paste(paste(list.files(get.value("filt_derep.path"), full.names = T, pattern = "cytB") %>% str_subset(pattern = "R1"), collapse = " "), 
              ">",
              paste0(get.value("filt_derep.1file.path"), "/all.cytB.R1.fasta"),
              sep = " "
              )

cmd3 <- paste(paste(list.files(get.value("filt_derep.path"), full.names = T, pattern = "cytB") %>% str_subset(pattern = "R2"), collapse = " "), 
              ">",
              paste0(get.value("filt_derep.1file.path"), "/all.cytB.R2.fasta"),
              sep = " "
              )

system2("cat", cmd1, stdout = T, stderr = T)
system2("cat", cmd2, stdout = T, stderr = T)
system2("cat", cmd3, stdout = T, stderr = T)

# derep these big files


for(x in list.files(get.value("filt_derep.1file.path"), full.names = T)){
  print(x)
  
  cmd <- paste("-derep_fulllength", x,
               "-output", x %>% str_replace(".fasta", "_derep.fasta"),
               "-sizeout",
               "-sizein",
               "-minuniquesize", "2",
               sep = " ")
  
  system2("vsearch", cmd, stdout = T, stderr = T)
  
  gc(verbose = F) # unload memory
}

# Make OTU

files <- list.files(get.value("filt_derep.1file.path"), full.names = T, pattern = "_derep.fasta")

new.files <- files %>% str_replace(get.value("filt_derep.1file.path"), get.value("OTU.usearch")) %>% 
                      str_replace("_derep.fasta", "_OTU.fasta")


for(x in 1:length(files)){

  cmd <- paste("-cluster_otus", files[x], 
               "-otus", new.files[x], 
               "-uparseout ", new.files[x] %>% str_replace(".fasta", "_OTUtab.txt"), 
               "-relabel", "OTU_",
               "-strand ", "plus", 
               sep=" ")
  
  system2("usearch", cmd, stdout=T, stderr=T)  

  #gc(verbose = F)
  
}

# Compare OTU to derep

# 12S

files <- list.files(get.value("filt_derep.path"), pattern = "12s", full.names = T)
new.files <- files %>% str_replace(get.value("filt_derep.path"),get.value("Compare.OTU.usearch")) %>%  
                       str_replace(".fasta",".txt")


for(x in 1:length(files)){
  
  print(files[x])
  
  cmd <- paste("-usearch_global", files[x], 
               "-db", list.files(get.value("OTU.usearch"), pattern = "12s", full.names = T) %>% str_subset("OTU.fasta"), 
               "-strand", "plus",  
               "-id", "0.97",
               "-blast6out", new.files[x], 
               "-maxhits", "1", 
               "-maxaccepts", "1", 
               "-maxrejects", "32",
               sep=" ")

  system2("usearch", cmd, stdout=T, stderr=T)

}

# cytb.R1

files <- list.files(get.value("filt_derep.path"), pattern = "cytB", full.names = T) %>% str_subset("R1")
new.files <- files %>% str_replace(get.value("filt_derep.path"),get.value("Compare.OTU.usearch")) %>%  
  str_replace(".fasta",".txt")

for(x in 1:length(files)){
  
  print(files[x])
  
  cmd <- paste("-usearch_global", files[x], 
               "-db", list.files(get.value("OTU.usearch"), pattern = "cytB.R1", full.names = T) %>% str_subset("OTU.fasta"), 
               "-strand", "plus",  
               "-id", "0.97",
               "-blast6out", new.files[x], 
               "-maxhits", "1", 
               "-maxaccepts", "1", 
               "-maxrejects", "32",
               sep=" ")
  
  system2("usearch", cmd, stdout=T, stderr=T)
  
}

# cytb.R2

files <- list.files(get.value("filt_derep.path"), pattern = "cytB", full.names = T) %>% str_subset("R2")
new.files <- files %>% str_replace(get.value("filt_derep.path"),get.value("Compare.OTU.usearch")) %>%  
  str_replace(".fasta",".txt")

for(x in 1:length(files)){
  
  print(files[x])
  
  cmd <- paste("-usearch_global", files[x], 
               "-db", list.files(get.value("OTU.usearch"), pattern = "cytB.R2", full.names = T) %>% str_subset("OTU.fasta"), 
               "-strand", "plus",  
               "-id", "0.97",
               "-blast6out", new.files[x], 
               "-maxhits", "1", 
               "-maxaccepts", "1", 
               "-maxrejects", "32",
               sep=" ")
  
  system2("usearch", cmd, stdout=T, stderr=T)
  
}

empty <- file.info(new.files)$size==0

new.files[!empty]

list.files(get.value("filt_derep.1file.path"), full.length = T)

getwd()

#Sys.setenv(PATH = paste(Sys.getenv("PATH"), path_to_usearch, sep= .Platform$path.sep))


# Add a cut adapt stuff


Amorces %>% filter(Locus == "12s") %>% pull(Amorce)
Amorces %>% filter(Locus == "cytB") %>% pull(Amorce) 

get.value("filt_cutadapt.path")

# l'intégrer plus tard dans le code ....
# car j'ai besoin de garder le même N reads avant/après ... voir comment trim le fait.







# END of the script -------------------------------------------------------


cat("\n-------------------------\n",
    paste0("Seqtab saved in: ",file.path(result.path, "Seqtab.data")),
    "\n-------------------------\n",  
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = T, sep = "\n")


rm(list = c("msg1","msg2"))

save.image(file.path(log.path,"Process_RAW.Rdata"))

end.time <- round(Sys.time() - start.time,2) 

cat(paste0("Rdata saved: ", file.path(log.path,"Process_RAW.Rdata")),
    paste0("\nTime to run this script: ", end.time, " ",units(end.time)),
    "\n-------------------------\n",  
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = T, sep = "\n")


