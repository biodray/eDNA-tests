
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

library(dada2); packageVersion("dada2") # Faire mettre cette info dans le log

library(Biostrings)

# Data --------------------------------------------------------------------

info.path <- "./00_Data/00_FileInfos"

raw.path     <- "./00_Data/01a_RawData" 
raw_unz.path <- "./00_Data/01b_RawData_unzipped" 

ref.path <- "./00_Data/03_RefSeq"

filt_dada2.path <- "./00_Data/02a_Filtered_dada2"
filt_IBIS.path <- "./00_Data/02b_Filtered_IBIS"
filt_JAMP.path  <- "./00_Data/02c_Filtered_JAMP"

log.path <- "./03_Log" 

seqtab.path <- "./01_Results/01_Seqtab"

Sample.xl <- "DB_Echantillons.xlsx"

DataSample <- read_excel(file.path(info.path,Sample.xl),sheet="DataSample",na="NA",guess_max=100000)
DataSeq    <- read_excel(file.path(info.path,Sample.xl),sheet="DataSeq",na="NA",guess_max=100000)
Amorces    <- read_excel(file.path(info.path,Sample.xl),sheet="Amorces",na="NA",guess_max=100000)
Inventaire <- read_excel(file.path(info.path,Sample.xl),sheet="DataLac",na="NA",guess_max=100000)

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
  
list.raw.files <- function(LOCI, 
                           PATH, 
                           FRpattern = c("R1", "R2"), 
                           STARTpattern = "EP-",
                           ENDpattern = "_L001_R1_001.fastq"){
  
  # Create a list to return all values 
  res <- list()
  
  for(x in LOCI){
    # List all files based on loci name
    files.all <- list.files(PATH, pattern = x)
    
    # Divide F and R reads
    files.R1 <- mixedsort(files.all[stringr::str_detect(files.all, pattern = FRpattern[1])])
    files.R2 <- mixedsort(files.all[stringr::str_detect(files.all, pattern = FRpattern[2])])

    #
    names <- files.R1 %>% stringr::str_remove(pattern = STARTpattern) %>% 
                          stringr::str_remove(pattern = ENDpattern) %>% 
                          mixedsort()
    
    names <- sapply(stringr::strsplit(names, "_"), `[`, 1)
    
    # list names
    files.R1.ln <- paste0(x, ".raw.files.R1")
    files.R2.ln <- paste0(x, ".raw.files.R2")   
    names.ln    <- paste0(x, ".names")
    
    # Put restults in my list
    res[[files.R1.ln]] <- files.R1
    res[[files.R2.ln]] <- files.R2
    res[[names.ln]]    <- names
    
    
    
    if(!is.null(log.path)){
      cat(paste0(x, ": ", length(names), " files found in ", PATH),  
        file=file.path(log.path, "Process_RAW.log.txt"), 
        append = T, sep = "\n")
      
    }

    
  }
  
  if(!is.null(log.path)){
    cat("\n-------------------------\n",  
        file=file.path(log.path, "Process_RAW.log.txt"), 
        append = T, sep = "\n")
    
  }
  
  return(res)
  
}

# add filt path 

add.filt.files <- function(LOCI, 
                           PATH,
                           FILE.LS,
                           F.PATTERN = "_F_filt.fastq.gz", 
                           R.PATTERN = "_R_filt.fastq.gz"){
  
  # Create a list to return all values 
  for(x in LOCI){
    
    NAMES.LN <-  paste0(x, ".names")
    NAMES <- FILE.LS[[NAMES.LN]]
    
    # Divide F and R reads
    filt.Fs <- file.path(PATH, paste0(NAMES, F.PATTERN))
    filt.Rs <- file.path(PATH, paste0(NAMES, R.PATTERN))    
    
    # list names
    files.R1.ln <- paste0(x, ".filt.files.R1")
    files.R2.ln <- paste0(x, ".filt.files.R2")   
    
    # Put restults in my list
    FILE.LS[[files.R1.ln]] <- filt.Fs
    FILE.LS[[files.R2.ln]] <- filt.Rs
  }
  
  return(FILE.LS)
  
}

add.filt.OK.files <- function(LOCI, 
                              PATH,
                              FILE.LS,
                              F.PATTERN = "_F_filt.fastq.gz", 
                              R.PATTERN = "_R_filt.fastq.gz"){
  
  # Create a list to return all values 
  for(x in LOCI){
    
    # List all files based on loci name
    files.all <- list.files(PATH, pattern = x)
    
    # Divide F and R reads
    files.R1 <- mixedsort(files.all[stringr::str_detect(files.all, pattern = F.PATTERN)])
    files.R2 <- mixedsort(files.all[stringr::str_detect(files.all, pattern = R.PATTERN)])
    
    names <- files.R1 %>% stringr::str_remove(pattern =  F.PATTERN)
    
    # list names
    files.R1.ln <- paste0(x, ".filt.files.R1.OK")
    files.R2.ln <- paste0(x, ".filt.files.R2.OK")
    
    names.ln <- paste0(x, ".filt.names")
    
    # Put results in my list
    FILE.LS[[files.R1.ln]] <- files.R1
    FILE.LS[[files.R2.ln]] <- files.R2
    FILE.LS[[names.ln]]    <- names    
  }
  
  return(FILE.LS)
  
}


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



all.files <- list.raw.files(LOCI = c("12s", "cytB"), PATH = raw_unz.path)


str(all.files)


if(file.exists("01_Results/QualityProfile.12S.raw.pdf")){
  
  msg1 <- "12S: Quality assement not done on raw files (previously done)"

} else {
  
  pdf("01_Results/QualityProfile.12S.raw.pdf", width=8, height=8,paper='letter') 

    plotQualityProfile(file.path(raw_unz.path, all.files[["12s.raw.files.R1"]]), aggregate = T)
    plotQualityProfile(file.path(raw_unz.path, all.files[["12s.raw.files.R2"]]), aggregate = T)
  
  dev.off()

  msg1 <- "12S: Quality assement done on raw files: 01_Results/QualityProfile.12S.raw.pdf"  

}


if(file.exists("01_Results/QualityProfile.CYTB.raw.pdf")){
  
  msg2 <- "cytB: Quality assement not done on raw files (previously done)"

} else {

   pdf("01_Results/QualityProfile.CYTB.raw.pdf", width=8, height=8,paper='letter') 

    plotQualityProfile(file.path(raw_unz.path, all.files[["cytB.raw.files.R1"]]), aggregate = T)
    plotQualityProfile(file.path(raw_unz.path, all.files[["cytB.raw.files.R2"]]), aggregate = T)

  dev.off()  
  
  msg2 <- "cytB: Quality assement done on raw files: 01_Results/QualityProfile.cytB.raw.pdf"  

}

cat(msg1, msg2, "\n-------------------------\n",  
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = T, sep = "\n")


# DADA2: raw to ASV ------------------------------------------------------------

all.files <-  add.filt.files(LOCI = c("12s", "cytB"), PATH = filt_dada2.path, FILE.LS = all.files) 

str(all.files)

# Filtrage en soi

Amorces %>% filter(Locus == "12s") %>% pull(Amorce) %>% nchar()
Amorces %>% filter(Locus == "cytB") %>% pull(Amorce) %>% nchar()

# Real filtering


filter.12s.summary <- filterAndTrim(fwd = file.path(raw_unz.path, all.files[["12s.raw.files.R1"]]),
                                         filt = all.files[["12s.filt.files.R1"]],
                                         rev = file.path(raw_unz.path, all.files[["12s.raw.files.R2"]]),
                                         filt.rev = all.files[["12s.filt.files.R2"]],
                                         truncQ=10, # minimum Q score, 10 = 90% base call accuracy
                                         truncLen = c(100,100), # Taille min/max des reads
                                         trimLeft= c(18, 18), # Enlever les seq des amorces, 
                                         maxEE=c(1,1),
                                         #orient.fwd = c("ACTGG"), # debut de l'amorce F
                                         compress = TRUE,
                                         multithread=FALSE, # TRUE on linux
                                         verbose = TRUE) 


filter.cytB.summary <- filterAndTrim(fwd = file.path(raw_unz.path, all.files[["cytB.raw.files.R1"]]),
                                    filt = all.files[["cytB.filt.files.R1"]],
                                    rev = file.path(raw_unz.path, all.files[["cytB.raw.files.R2"]]),
                                    filt.rev = all.files[["cytB.filt.files.R2"]],
                                    truncQ=10, # minimum Q score, 10 = 90% base call accuracy
                                    truncLen = c(100,100), # Taille min/max des reads
                                    trimLeft= c(35, 34), # Enlever les seq des amorces
                                    maxEE=c(1,1),
                                    compress = TRUE,
                                    multithread=FALSE, # TRUE on linux
                                    verbose = TRUE) 


cat(paste0("\n Correlation between 12S and cytB reads after DADA2 filtration : ",
    round(cor(filter.12s.summary[,"reads.out"], filter.cytB.summary[,"reads.out"], method = "spearman"),2)
     ),
    "\n-------------------------\n",
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = T, sep = "\n")


# Add files (beacause some were discard)

all.files <- add.filt.OK.files(LOCI = c("12s","cytB"), 
                              PATH = filt_dada2.path,
                              FILE.LS = all.files,
                              F.PATTERN = "_F_filt.fastq.gz", 
                              R.PATTERN = "_R_filt.fastq.gz")

# qUALITY ASSEMENT

pdf("01_Results/QualityProfile.12S.filter.pdf") 
  
  plotQualityProfile(file.path(filt_dada2.path,all.files[["12s.filt.files.R1.OK"]]), aggregate = T)
  plotQualityProfile(file.path(filt_dada2.path,all.files[["12s.filt.files.R2.OK"]]), aggregate = T)
  
dev.off()  
  
msg1 <- "12S: Quality assement done on DADA2 filtered files: 01_Results/QualityProfile.12S.filter.pdf"  
  


pdf("01_Results/QualityProfile.CYTB.filter.pdf") 
  
  plotQualityProfile(file.path(filt_dada2.path,all.files[["cytB.filt.files.R1.OK"]]), aggregate = T)
  plotQualityProfile(file.path(filt_dada2.path,all.files[["cytB.filt.files.R2.OK"]]), aggregate = T)
  
dev.off()  
  
  msg2 <- "cytB: Quality assement done on DADA2 filtered files: 01_Results/QualityProfile.CYTB.filter.pdf"  


cat(msg1, msg2, "\n-------------------------\n",  
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = T, sep = "\n")


# Calcul du taux d'erreur (DADA2)

err.12s.F <- learnErrors(file.path(filt_dada2.path,all.files[["12s.filt.files.R1.OK"]])) 
err.12s.R <- learnErrors(file.path(filt_dada2.path,all.files[["12s.filt.files.R2.OK"]])) 


pdf("01_Results/ErrorsRate.dada2.12S.pdf") 
  plotErrors(err.12s.F, nominalQ=TRUE)
  plotErrors(err.12s.R, nominalQ=TRUE)
dev.off()

msg1 <- "12S: DADA2 error rate plot saved: 01_Results/ErrorsRate.dada2.12S.pdf"  


err.cytB.F <- learnErrors(file.path(filt_dada2.path,all.files[["cytB.filt.files.R1.OK"]])) 
err.cytB.R <- learnErrors(file.path(filt_dada2.path,all.files[["cytB.filt.files.R2.OK"]])) 

pdf("01_Results/ErrorsRate.dada2.CYTB.pdf") 
  plotErrors(err.cytB.F, nominalQ=TRUE)
  plotErrors(err.cytB.R, nominalQ=TRUE)
dev.off()

msg1 <- "cytB: DADA2 error rate plot saved: 01_Results/ErrorsRate.dada2.CYTB.pdf"  

cat(msg1, msg2, "\n-------------------------\n",  
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = T, sep = "\n")

# Déréplication

derep.12s.Fs <- derepFastq(file.path(filt_dada2.path,all.files[["12s.filt.files.R1.OK"]]))
names(derep.12s.Fs) <- all.files[["12s.filt.names"]] 

derep.12s.Rs <- derepFastq(file.path(filt_dada2.path,all.files[["12s.filt.files.R2.OK"]]))
names(derep.12s.Rs) <- all.files[["12s.filt.names"]]

derep.cytB.Fs <- derepFastq(file.path(filt_dada2.path,all.files[["cytB.filt.files.R1.OK"]]))
names(derep.cytB.Fs) <- all.files[["cytB.filt.names"]]

derep.cytB.Rs <- derepFastq(file.path(filt_dada2.path,all.files[["cytB.filt.files.R2.OK"]]))
names(derep.cytB.Rs) <- all.files[["cytB.filt.names"]]


# Inférence des échantillons

dada.12s.Fs <- dada(derep.12s.Fs, 
                    err = err.12s.F, 
                    multithread=FALSE,
                    pool=TRUE,
                    HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)

dada.12s.Rs <- dada(derep.12s.Rs, 
                    err = err.12s.R, 
                    multithread=FALSE,
                    pool=TRUE,
                    HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)



dada.cytB.Fs <- dada(derep.cytB.Fs, 
                    err = err.cytB.F, 
                    multithread=FALSE,
                    pool=TRUE,
                    HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)

dada.cytB.Rs <- dada(derep.cytB.Rs, 
                    err = err.cytB.R, 
                    multithread=FALSE,
                    pool=TRUE,
                    HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)



## DADA2 Options: Non-Illumnina sequencing technologies

mergers.12S <- mergePairs(dadaF = dada.12s.Fs, 
                               derepF = derep.12s.Fs, 
                               dadaR = dada.12s.Rs, 
                               derepR = derep.12s.Rs, 
                               minOverlap = 50, 
                               maxMismatch = 0,
                               returnRejects = TRUE,
                               verbose=TRUE)

mergers.cytB <- mergePairs(dadaF = dada.cytB.Fs, 
                          derepF = derep.cytB.Fs, 
                          dadaR = dada.cytB.Rs, 
                          derepR = derep.cytB.Rs, 
                          maxMismatch = 0,
                          returnRejects = TRUE,
                          verbose=TRUE,
                          justConcatenate=TRUE)


# Make sequence table

seqtab.12s.int    <- makeSequenceTable(mergers.12S)
seqtab.12S.F.int  <- makeSequenceTable(dada.12s.Fs)
seqtab.12S.R.int  <- makeSequenceTable(dada.12s.Rs)

seqtab.cytB.int   <- makeSequenceTable(mergers.cytB)
seqtab.cytB.F.int <- makeSequenceTable(dada.cytB.Fs)
seqtab.cytB.R.int <- makeSequenceTable(dada.cytB.Rs)


# Remove chimera

seqtab.12s <- removeBimeraDenovo(seqtab.12s.int, method = "consensus", 
                                 multithread = FALSE, verbose = TRUE)

seqtab.12s.F <- removeBimeraDenovo(seqtab.12S.F.int, method = "consensus", 
                                 multithread = FALSE, verbose = TRUE)

seqtab.12s.R <- removeBimeraDenovo(seqtab.12S.R.int, method = "consensus", 
                                   multithread = FALSE, verbose = TRUE)

seqtab.cytB <- removeBimeraDenovo(seqtab.cytB.int, method = "consensus", 
                                  multithread = FALSE, verbose = TRUE)

seqtab.cytB.F <- removeBimeraDenovo(seqtab.cytB.F.int, method = "consensus", 
                                    multithread = FALSE, verbose = TRUE)

seqtab.cytB.R <- removeBimeraDenovo(seqtab.cytB.R.int, method = "consensus", 
                                    multithread = FALSE, verbose = TRUE)   

# Stats on chimera
sum(seqtab.12s)/sum(seqtab.12s.int)
sum(seqtab.12s.F)/sum(seqtab.12S.F.int)
sum(seqtab.12s.R)/sum(seqtab.12S.R.int)

sum(seqtab.cytB)/sum(seqtab.cytB.int)
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
    paste0("\nSeqtab saved in: ",file.path(seqtab.path, "Dada2.Seqtab.data")),
    "\n-------------------------\n",  
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = T, sep = "\n")

# Save Seqtab

save(file = file.path(seqtab.path, "Dada2.Seqtab.data"), 
     list = c("seqtab.12s",
              "seqtab.12s.F",
              "seqtab.12s.R",
              "seqtab.cytB",
              "seqtab.cytB.F",
              "seqtab.cytB.R"))

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



all.files[["IBIS.files"]] <- list.files(filt_IBIS.path, pattern = "stability.trim.contigs.good.unique.good.filter.unique.precluster.fn.0.06.rep.pick")


all.files[["IBIS.files.names"]] <- all.files[["IBIS.files"]] %>% str_remove(pattern = "EP_") %>% 
                                   str_remove(pattern = "_stability.trim.contigs.good.unique.good.filter.unique.precluster.fn.0.06.rep.pick.fasta") %>% 
                                   stringr::str_replace_all("_", "-")


seqtab.12s.IBIS <-  makeSeqTabFromScratch(files = all.files[["IBIS.files"]],
                                          name = all.files[["IBIS.files.names"]],
                                          path = filt_IBIS.path) 


# Save Seqtab

save(file = file.path(seqtab.path, "Dada2.Seqtab.IBIS.data"), 
     list = c("seqtab.12s.IBIS")
     )


cat("\n-------------------------\n",
    paste0("IBIS: Seqtab saved in: ",file.path(seqtab.path, "Dada2.Seqtab.data")),
    "\n-------------------------\n",  
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = T, sep = "\n")




# JAMP: raw to OTU -------------------------------------------------------------

# library("JAMP")


#U_merge_PE(file1 = file.path(raw_unz.path,all.files[["12s.raw.files.R1"]]),  
#           file2 = file.path(raw_unz.path,all.files[["12s.raw.files.R2"]]),
#           exe = "C:\\Users\\BourretA\\Documents\\Programs\\USearch\\usearch11.0.667_win32.exe")


#U_cluster_otus(files= file.path(raw_unz.path,all.files[["cytB.filt.files.R2.OK"]]), filter=0.01)




#Sys.setenv(PATH = paste(Sys.getenv("PATH"), path_to_usearch, sep= .Platform$path.sep))



# END of the script -------------------------------------------------------

rm(list = c("msg1","msg2"))

save.image(file.path(log.path,"Process_RAW.Rdata"))

end.time <- round(Sys.time() - start.time,2) 

cat(paste0("Rdata saved: ", file.path(log.path,"Process_RAW.Rdata")),
    paste0("\nTime to run this script: ", end.time, " ",units(end.time)),
    "\n-------------------------\n",  
    file=file.path(log.path, "Process_RAW.log.txt"), 
    append = T, sep = "\n")


