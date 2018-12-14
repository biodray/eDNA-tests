
# Info --------------------------------------------------------------------

# Run the DADA2 and JAMP pipelines in parallele, on eDNA data
# 
# Audrey Bourret
# 2018-10-22
#

start.time <- Sys.time()

# Library -----------------------------------------------------------------

# Data manipulation

library(tidyverse) # includes ggplot2 dplyr and stringr

library(gtools)    # for mixedsort
library(readxl)

#library(devtools)
#devtools::install_github("kassambara/ggpubr")
library(ggpubr)    # on github - for nice graphs

# Fastq and fasta manipulation

library(Biostrings)

library(dada2); packageVersion("dada2") # Faire mettre cette info dans le log
#library(JAMP) - I don't think this package is still necessary

#devtools::install_github("kassambara/fastqcr")
library(fastqcr)  # to read fastqc ressults

# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}


# Data --------------------------------------------------------------------

#DataSample <- read_excel(get.value("Sample.xl"),sheet="DataSample",na="NA",guess_max=100000)
#DataSeq    <- read_excel(get.value("Sample.xl"),sheet="DataSeq",na="NA",guess_max=100000)
#Amorces    <- read_excel(get.value("Sample.xl"),sheet="Amorces",na="NA",guess_max=100000)
#Inventaire <- read_excel(get.value("Sample.xl"),sheet="DataLac",na="NA",guess_max=100000)

# Format de DataSeq avec les données IBIS (p1-A1)

#DataSeq <- DataSeq %>% mutate (IbisID = paste0("p",Plaque,"-",Puit)) %>% 
#  left_join(DataSample %>% select(SampleID, NomLac, CatSite), by = "SampleID")

# Setting up the log file

if(file.exists(get.value("Raw.log"))){

  switch(menu(title = "Do you want to erase the previous RAW log files?", graphics = T, 
            choice = c("yes", "no")
            )+1,
       # Answer 0
       cat("Nothing done\n"),
       # Answer 1 (yes)
      {cat("\nA new log file Process_RAW.log.txt was created (the previous RAW log was erased)\n\n")
       cat("\n-------------------------\n", 
           "Process raw eDNA data\n",
           date(),
           "\n-------------------------\n", 
           file=get.value("Raw.log"), 
           append = FALSE, sep = "\n")
           }, 
       # Answer 2 (no) 
      {cat ("\nInformation will be append to the log file Process_RAW.log.txt\n\n")
       cat("\n\n-------------------------\n",
            "NEW ANALYSIS PERFORMED", 
            date(),
            "\n-------------------------\n", 
            file=get.value("Raw.log"), 
            append = T, sep = "\n")
      }
        )

} else {
  cat ("\nThe log file Process_RAW.log.txt was created\n\n")  
  cat("\n-------------------------\n", 
      "Process raw eDNA data\n",
      date(),
      "\n-------------------------\n", 
      file=get.value("Raw.log"), 
      append = FALSE, sep = "\n")
}

# Functions ---------------------------------------------------------------

# Function to ugrade plotQualityProfiles

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

# To work with IBIS data - unnecessary now

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

# 1. Quality assesment - RAW (fastqc + dada2) ------------------------------------------------------

# Small code to no redo this part if has been previously done

if(file.exists(get.value("Qplot.RAW.data"))){
  
  switch(menu(title = "Do you want to re-plot RAW data quality?", graphics = T, 
              choice = c("yes", "no")
  )+1,
  # Answer 0
  cat("Nothing done\n"),
  # Answer 1
 {cat("\nRaw data are re-plot (can take some time)\n\n")
   graph.tot.ls   <- plotQplus(list.files(get.value("raw_unz_rename.path"), full.names = T),
                               locus = c("12s", "cytB"), pattern = c("R1", "R2"))
   
   graph.sample.ls <- plotQplus(list.files(get.value("raw_unz_rename.path"), pattern = "Sample", full.names = T),
                                locus = c("12s", "cytB"), pattern = c("R1", "R2"))
   
   save(file = get.value("Qplot.RAW.data"), 
        list = c("graph.tot.ls" , "graph.sample.ls"))  
  
   cat("Raw data quality were replot and saved:",
       get.value("Qplot.RAW.data"),
       "\n-------------------------\n", 
       file=get.value("Raw.log"), 
       append = T, sep = "\n") 
  } , 
  # Answer 2 (no)
 {cat("\nRaw data are not re-plot (old data are reloaded)\n\n")
  load(get.value("Qplot.RAW.data"))
}
  )

} else {
  cat("\nRaw data plot \n\n")
  graph.tot.ls   <- plotQplus(list.files(get.value("raw_unz_rename.path"), full.names = T),
                              locus = c("12s", "cytB"), pattern = c("R1", "R2"))
  
  graph.sample.ls <- plotQplus(list.files(get.value("raw_unz_rename.path"), pattern = "Sample", full.names = T),
                               locus = c("12s", "cytB"), pattern = c("R1", "R2"))
  
  save(file = get.value("Qplot.RAW.data"), 
       list = c("graph.tot.ls" , "graph.sample.ls"))
  
  cat("Raw data quality were plot and saved:",
      get.value("Qplot.RAW.data"),
      "\n-------------------------\n", 
      file=get.value("Raw.log"), 
      append = T, sep = "\n") 
}

# devrait peut-être aller dans un autre script - J'aimerais garder ce script au min

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
switch(menu(title = "Do you want to save RAW data quality plots in PDF and PNG?", graphics = T, 
            choice = c("yes", "no")
            )+1,
       # Answer 0
       cat("Nothing done\n"),
       # Answer 1
       { cat("\nRaw data plot were saved\n\n")
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
       }, 
       # Answer 2 (no)
       cat("\nRaw data quality plot were not saved\n\n")
)


# FastQC

if(get_os() %in% c("os","linux")){ # to run only on os and not windows
  
  file.remove(list.files(get.value("result.FQraw.path"), full.name = T, pattern ="fastqc"))
  
  cmd <- paste("--outdir", get.value("result.FQraw.path"), list.files(get.value("raw_unz_rename.path"), full.names = T))
  
  system2("fastqc", cmd)
  
  # Save info on the log file
  cat("FastQC analysis was performed:",
      get.value("result.FQraw.path"),
      "\n-------------------------\n", 
      file=get.value("Raw.log"), 
      append = T, sep = "\n") 
  
} else {cat("Cannot perform FastQC on windows yet -- sorry!!")}

# Save graphs
switch(menu(title = "Do you want to aggregate FastQC results?", graphics = T, 
            choice = c("yes", "no")
            )+1,
       # Answer 0
       cat("Nothing done\n"),
       # Answer 1
       {
       qc.12s.res  <- qc_JAMP_sum(list.files(get.value("result.FQraw.path") %>% mixedsort(), pattern = "fastqc.zip", full.names = T) %>% str_subset(pattern = "12s"))
       qc.cytB.res <- qc_JAMP_sum(list.files(get.value("result.FQraw.path") %>% mixedsort(), pattern = "fastqc.zip", full.names = T) %>% str_subset(pattern = "cytB"))
       
       qc_JAMP_plot(qc.12s.res$summary, file = file.path(get.value("result.Q.path"),"FastQC_12s_summary.pdf" ))
       qc_JAMP_plot(qc.cytB.res$summary, file = file.path(get.value("result.Q.path"),"FastQC_cytB_summary.pdf" ))
     
       save(file = get.value("FastQC.data"), 
            list = c("qc.12s.res" , "qc.cytB.res"))
       
       cat("FastQC results were aggregated:",
           get.value("result.Q.path"),
           get.value("FastQC.data"),
           "\n-------------------------\n", 
           file=get.value("Raw.log"), 
           append = T, sep = "\n") 
       
       },
       if(file.exists(get.value("FastQC.data"))){
         cat("Old data were reloaded\n")
         load(get.value("FastQC.data"))
         } else {
         cat("Nothing done\n") 
         }
       )

# Write the result somewhere
# write.csv(exp, paste(folder, "/FastQC/stats.csv", sep=""))


# 2. RAW to FILT (cutadapt + dada2) ------------------------------------------------------------

# 2.1 CUTADAPT

if(get_os() %in% c("os","linux")){ 

  # Files
  
  files1 <- list.files(get.value("raw_unz_rename.path"), pattern = "12s", full.names=T) %>% str_subset("R1")
  files2 <- files1 %>% str_replace("R1","R2")
  files3 <- files1 %>% str_replace("12s","cytB")
  files4 <- files2 %>% str_replace("12s","cytB")
  
  cutadapt.files <- function(files){
                             new.files <- files %>% str_replace(get.value("raw_unz_rename.path"), get.value("filt_cutadapt.path")) %>%
                                                    str_replace (".fastq", "_cut.fastq")
                             }                         
    
  new.files1 <- files1 %>% cutadapt.files() 
  new.files2 <- files2 %>% cutadapt.files() 
  new.files3 <- files3 %>% cutadapt.files() 
  new.files4 <- files4 %>% cutadapt.files() 

  # Remove old files
  file.remove(list.files(get.value("filt_cutadapt.path"), full.name = T, pattern =".fastq"))
  file.remove(list.files(get.value("filt_cutadapt.log"), full.name = T, pattern ="_log.txt"))
  
  # 12S
  
  for(x in 1:length(files1)){
    
    cat("\n",files1[x], sep="\n")
    
    cmd <- paste("-g ^ACTGGGATTAGATACCCC",
                                  "-G ^TAGAACAGGCTCCTCTAG",
                                  "-o", new.files1[x], 
                                  "--paired-output", new.files2[x], 
                                  files1[x],
                                  files2[x],
                                  "-f", "fastq",      
                                  "--discard-untrimmed", 
                                  #"--report=minimal",
                                  sep = " ") # forward adapter
    
    A <- system2("cutadapt", cmd, stdout=T, stderr=T) # to R console
    
    # TO UPDaTE
    
    # save a file log
    cat(file = new.files1[x] %>% str_replace(get.value("filt_cutadapt.path"), get.value("filt_cutadapt.log")) %>% 
          str_replace(".fastq","_log.txt") %>% 
          str_remove("_R1"),
          A, # what to put in my file
          append= F, sep = "\n")
    

  }

  cat("12S adapters were removed with cutadapt:",
      get.value("filt_cutadapt.path"),
      get.value("filt_cutadapt.log"),
      "\n-------------------------\n", 
      file=get.value("Raw.log"), 
      append = T, sep = "\n") 
  
  # CYTB

  for(x in 1:length(files3)){
    
    cat("\n",files3[x], sep="\n")
    
    #F - CYTB
    cmd1 <- paste("-g ^AAAAAGCTTCCATCCAACATCTCAGCATGATGAAA",
                                 #"--minimum-length", "75", #  After truncation
                                 "-o", new.files3[x], 
                                 files3[x], 
                                 "-f", "fastq",
                                 "--discard-untrimmed", 
                                 "--report=minimal", 
                                 sep = " ") # forward adapter
    
    A <- system2("cutadapt", cmd1, stdout=T, stderr=T) # to R consol
    
    #r - CYTB
    cmd2 <- paste("-g ^AAACTGCAGCCCCTCAGAATGATATTTGTCCTCA", # pas en reverse complement car encore du bon sens
                                 #"--minimum-length", "75", #  After truncation
                                 "-o", new.files4[x], 
                                 files4[x], 
                                 "-f", "fastq",
                                 "--discard-untrimmed", 
                                 "--report=minimal", 
                                 sep = " ") # forward adapter
    
    B <- system2("cutadapt", cmd2, stdout=T, stderr=T) # to R consol
    
    # Write a log file
    
    cat(file = new.files3[x] %>% str_replace(get.value("filt_cutadapt.path"), get.value("filt_cutadapt.log")) %>% 
                                 str_replace(".fastq","_log.txt"),
        A, # what to put in my file
       append= F, sep = "\n")

    cat(file = new.files4[x] %>% str_replace(get.value("filt_cutadapt.path"), get.value("filt_cutadapt.log")) %>% 
                                 str_replace(".fastq","_log.txt"),
         B, # what to put in my file
         append= F, sep = "\n")
    
  }

  cat("cytB adapters were removed with cutadapt:",
      get.value("filt_cutadapt.path"),
      "\n-------------------------\n", 
      file=get.value("Raw.log"), 
      append = T, sep = "\n") 
  
} else {cat("Cannot perform cutadapt on windows yet -- sorry!!")}


# 2.2 DADA2

files1 <- list.files(get.value("filt_cutadapt.path"), pattern = "12s", full.names = T) %>% str_subset("R1")
files2 <- files1 %>% str_replace("R1", "R2")
files3 <- files1 %>% str_replace("12s","cytB")
files4 <- files2 %>% str_replace("12s","cytB")

dada2.files <- function(files){
                        new.files <- files %>% str_replace(get.value("filt_cutadapt.path"),get.value("filt_dada2.path")) %>%
                                               str_replace (".fastq", "_filt.fastq")
                        }                         

new.files1 <- files1 %>% dada2.files() 
new.files2 <- files2 %>% dada2.files() 
new.files3 <- files3 %>% dada2.files() 
new.files4 <- files4 %>% dada2.files() 

# Remove old files
file.remove(list.files(get.value("filt_dada2.path"), full.name = T, pattern =".fastq"))

# 12S

filter.12s.summary <- filterAndTrim(fwd = files1,
                                    filt = new.files1,
                                    rev = files2,
                                    filt.rev = new.files2,
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

filter.cytB.R1.summary <- filterAndTrim(fwd = files3,
                                        filt = new.files3,
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

filter.cytB.R2.summary <- filterAndTrim(fwd = files4,
                                        filt = new.files4,
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

save(file = get.value("dada2.filt.data"), 
     list = c("filter.12s.summary" , "filter.cytB.R1.summary", "filter.cytB.R2.summary"))

cat("Quality filtering was performed with dada2:",
    get.value("filt_dada2.path"),
    get.value("dada2.filt.data"),
    "\n-------------------------\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 

# qUALITY ASSEMENT

# Small code to no redo this part if has been previously done

if(file.exists(get.value("Qplot.FILT.data"))){
  
  switch(menu(title = "Do you want to re-plot FILT data quality?", graphics = F, 
              choice = c("yes", "no")
  )+1,
  # Answer 0
  cat("Nothing done\n"),
  # Answer 1
  {cat("\nFilt data are re-plot (can take some time)\n\n")
    graph.tot.filt.ls   <- plotQplus(list.files(get.value("filt_dada2.path"), full.names = T),
                                     locus = c("12s", "cytB"), pattern = c("R1", "R2"))
    
    graph.sample.filt.ls <- plotQplus(list.files(get.value("filt_dada2.path"), pattern = "Sample", full.names = T),
                                      locus = c("12s", "cytB"), pattern = c("R1", "R2"))
    
    save(file = get.value("Qplot.FILT.data"), 
         list = c("graph.tot.filt.ls" , "graph.sample.filt.ls"))
    
    cat("Filt data quality were replot and saved:",
        get.value("Qplot.FILT.data"),
        "\n-------------------------\n", 
        file=get.value("Raw.log"), 
        append = T, sep = "\n") 
  } , 
  # Answer 2 (no)
  {cat("\nFilt data are not re-plot (old data are reloaded)\n\n")
    load(get.value("Qplot.FILT.data"))
  }
  )
  
} else {
  cat("\nFilt data quality plot \n\n")
  graph.tot.filt.ls   <- plotQplus(list.files(get.value("filt_dada2.path"), full.names = T),
                                   locus = c("12s", "cytB"), pattern = c("R1", "R2"))
  
  graph.sample.filt.ls <- plotQplus(list.files(get.value("filt_dada2.path"), pattern = "Sample", full.names = T),
                                    locus = c("12s", "cytB"), pattern = c("R1", "R2"))
  
  save(file = get.value("Qplot.FILT.data"), 
       list = c("graph.tot.filt.ls" , "graph.sample.filt.ls"))
  
  cat("Filt data quality were plot and saved:",
      get.value("Qplot.FILT.data"),
      "\n-------------------------\n", 
      file=get.value("Raw.log"), 
      append = T, sep = "\n") 
}


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


# Save graphs
switch(menu(title = "Do you want to save FILT data quality plots in PDF and PNG?", graphics = F, 
            choice = c("yes", "no")
)+1,
# Answer 0
cat("Nothing done\n"),
# Answer 1
{ cat("\nFILT data plot were saved\n\n")
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
  
}, 
# Answer 2 (no)
cat("\nFILT data quality plot were not saved\n\n")
)


# 3. FILT to ASV (denoise with DADA2) --------------------------------------------------

# 3.1 Calcul du taux d'erreur

files1 <- list.files(get.value("filt_dada2.path"), pattern = "12s", full.names=T) %>% str_subset("R1")
files2 <- files1 %>% str_replace("R1","R2")
files3 <- files1 %>% str_replace("12s","cytB")
files4 <- files2 %>% str_replace("12s","cytB")

err.12s.F <- learnErrors(files1) 
err.12s.R <- learnErrors(files2) 

err.cytB.F <- learnErrors(files3) 
err.cytB.R <- learnErrors(files4) 

pdf("01_Results/ErrorsRate.dada2.12S.pdf") 
  plotErrors(err.12s.F, nominalQ=TRUE)
  plotErrors(err.12s.R, nominalQ=TRUE)
dev.off()

pdf("01_Results/ErrorsRate.dada2.CYTB.pdf") 
  plotErrors(err.cytB.F, nominalQ=TRUE)
  plotErrors(err.cytB.R, nominalQ=TRUE)
dev.off()

cat("Dada2 error rate assement was performed and graph were saved:",
    "01_Results/ErrorsRate.dada2.12S.pdf",
    "01_Results/ErrorsRate.dada2.CYTB.pdf",
    "\n-------------------------\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 


# 3.2 Déréplication

derep.12s.F <- derepFastq(files1) 
derep.12s.R <- derepFastq(files2)

derep.cytB.F <- derepFastq(files3) 
derep.cytB.R <- derepFastq(files4) 

cat("Dada2 dereplication was performed.",
    "\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 


# 3.3 Inférence des échantillons

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

cat("Dada2 sample inference was performed.",
    "\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 

# 3.4 Dada2 Merger

mergers.12S <- mergePairs(dadaF = dada.12s.F, 
                          derepF = derep.12s.F, 
                          dadaR = dada.12s.R, 
                          derepR = derep.12s.R, 
                          minOverlap = 30, 
                          maxMismatch = 0,
                          returnRejects = FALSE,
                          verbose=TRUE)

cat("Dada2 samples were merged for 12s.",
    "\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 

# 3.5 Make sequence table

seqtab.12s.int    <- makeSequenceTable(mergers.12S)
#seqtab.12S.F.int  <- makeSequenceTable(dada.12s.Fs)
#seqtab.12S.R.int  <- makeSequenceTable(dada.12s.Rs)

#seqtab.cytB.int   <- makeSequenceTable(mergers.cytB)
seqtab.cytB.F.int <- makeSequenceTable(dada.cytB.F)
seqtab.cytB.R.int <- makeSequenceTable(dada.cytB.R)

# Remove chimera

ASVtab.12s <- removeBimeraDenovo(seqtab.12s.int, method = "consensus", 
                                 multithread = FALSE, verbose = TRUE)

#seqtab.12s.F <- removeBimeraDenovo(seqtab.12S.F.int, method = "consensus", 
#                                 multithread = FALSE, verbose = TRUE)

#seqtab.12s.R <- removeBimeraDenovo(seqtab.12S.R.int, method = "consensus", 
#                                   multithread = FALSE, verbose = TRUE)

#seqtab.cytB <- removeBimeraDenovo(seqtab.cytB.int, method = "consensus", 
#                                  multithread = FALSE, verbose = TRUE)

ASVtab.cytB.F <- removeBimeraDenovo(seqtab.cytB.F.int, method = "consensus", 
                                    multithread = FALSE, verbose = TRUE)

ASVtab.cytB.R <- removeBimeraDenovo(seqtab.cytB.R.int, method = "consensus", 
                                    multithread = FALSE, verbose = TRUE)   

# Write results

write.dada2.res <- function(tab, name, folder){

  file1 <- file.path(folder, paste0("all.", name, "_ASV.fasta"))
  file2 <- file1 %>% str_replace(".fasta", "_ASVtable.txt")
  
  DNA <- DNAStringSet(getSequences(tab))
  names(DNA) <- paste0("ASV_", 1:length(DNA))
  
  writeXStringSet(DNA, file1)
  write.table(tab, file2)
  
}

write.dada2.res(ASVtab.12s, "12s", get.value("ASV.dada2.path"))
write.dada2.res(ASVtab.cytB.F, "cytB.R1", get.value("ASV.dada2.path"))
write.dada2.res(ASVtab.cytB.R, "cytB.R2", get.value("ASV.dada2.path"))

cat("Sequence tables were created and chimeric sequences were removed.",
    "\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 


save(file = get.value("ASVtable.data"), 
     list = ls(pattern = "ASVtab."))


save(file = get.value("dada2.data"), 
     list = c(ls(pattern = "err.12s"),
              ls(pattern = "err.cytB"),
              ls(pattern = "derep.12s"),
              ls(pattern = "derep.cytB"),
              ls(pattern = "dada.12s"),
              ls(pattern = "dada.cytB"),
              ls(pattern = "mergers.12S")
              ))

cat("Data were saved:",
    get.value("ASVtable.data"),
    get.value("dada2.data"),
    "\n-------------------------\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 


load(get.value("ASVtable.data"))



# Save Seqtab somewhere!!


# 4. FILT to OTU (usearch inspired by JAMP)-------------------------------------------------------------

# 4.1 Umerge

files1 <- list.files(get.value("filt_dada2.path"), pattern = "12s", full.names=T) %>% str_subset("R1")
files2 <- files1 %>% str_replace("R1", "R2")
new.files <- files1 %>% str_replace(get.value("filt_dada2.path"), get.value("filt_merge.path")) %>% 
                      str_remove("_R1") %>% 
                      str_replace(".fastq", "_merged.fastq")  

# Remove old files
file.remove(list.files(get.value("filt_merge.path"), full.name = T, pattern =".fastq"))
file.remove(list.files(get.value("filt_merge.path"), full.name = T, pattern =".fastq"))

for(x in 1:length(files1)){
  
  print(files1[x])
  
  cmd <- paste("-fastq_mergepairs", files1[x], 
               "-reverse", files2[x],  
               "-fastqout", new.files[x], 
               "-report", new.files[x] %>% str_replace(get.value("filt_merge.path"),get.value("filt_merge.log")) %>% str_replace("_merged.fastq", "_log.txt"),
               "-fastq_maxdiffs", "99", 
               "-fastq_pctid", "75", 
               "-fastq_trunctail",  "0",
               "-fastq_minovlen", "30", 
               sep=" ")
  
  system2("usearch", cmd, stdout=T, stderr=T)  
  
}

cat("12S sequence were merged with usearch.",
    "\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 

# 4.2 Merged sequences length filtering

files <- list.files(get.value("filt_merge.path"), pattern = "12s", full.names=T)
new.files <- files %>% str_replace(get.value("filt_merge.path"), get.value("filt_min.path")) %>% 
                       str_replace(".fastq", "_min.fastq")

file.remove(list.files(get.value("filt_min.path"), full.name = T, pattern =".fastq"))

for(x in 1:length(files)){
  
  print(files[x])
  
  cmd <- paste(files[x],
               "-o", new.files[x], 
               #"-f", "fastq",  
               "-m", 100, # for minimum
               #"-M", x, # for maximu,
               sep=" ")
  
  A <- system2("cutadapt", cmd, stdout=T, stderr=T)  
  
  # save a log
  cat(file = new.files[x] %>% str_replace(get.value("filt_min.path"), get.value("filt_min.log")) %>% 
        str_replace(".fastq","_log.txt"),
      A, # what to put in my file
      append= F, sep = "\n")

}

cat("Merged sequences were filtered by length.",
    "\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 


# 4.3 Dereplicated files with vsearch

files <- c(list.files(get.value("filt_min.path"), pattern = "12s", full.names = T),
           list.files(get.value("filt_dada2.path"), pattern = "cytB", full.names = T)
           )

new.files <-  files %>% str_replace(get.value("filt_min.path"), get.value("filt_derep.path")) %>% 
                        str_replace(get.value("filt_dada2.path"), get.value("filt_derep.path")) %>% 
                        str_replace(".fastq","_derep.fasta")

file.remove(list.files(get.value("filt_derep.path"), full.name = T, pattern =".fasta"))

for(x in 1:length(files)){
  print(files[x])
  
  cmd <- paste("-derep_fulllength", files[x],
             "-output", new.files[x],
             "-sizeout",
             "-sizein",
             sep = " ")
  
  system2("vsearch", cmd, stdout = T, stderr = T)

}

# 4.4 Create one big file

file.remove(list.files(get.value("filt_derep.1file.path"), full.name = T, pattern =".fasta"))

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


# 4.5 Derep these big files

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

cat("Sequences were dereplicated:",
    get.value("filt_derep.path"),
    get.value("filt_derep.1file.path"),
    "\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 


# 4.6 Make OTU

files <- list.files(get.value("filt_derep.1file.path"), full.names = T, pattern = "_derep.fasta")

new.files <- files %>% str_replace(get.value("filt_derep.1file.path"), get.value("OTU.usearch")) %>% 
                       str_replace("_derep.fasta", "_OTU.fasta")

file.remove(list.files(get.value("OTU.usearch"), full.name = T, pattern ="_OTU"))


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

# 4.7 Compare OTU to derep

usearch_global <- function(files, new.files, DB) {
 
  for(x in 1:length(files)){
  
  print(files[x])
  
  cmd <- paste("-usearch_global", files[x], 
               "-db", DB, 
               "-strand", "plus",  
               "-id", "0.97",
               "-blast6out", new.files[x], 
               "-maxhits", "1", 
               "-maxaccepts", "1", 
               "-maxrejects", "32",
               sep=" ")
  
  A <- system2("usearch", cmd, stdout=T, stderr=T)
  
  # Write a log file for each clusters
  cat(file = new.files[x] %>% str_replace(get.value("Compare.OTU.usearch"), get.value("Compare.OTU.usearch.log")) %>% 
                   str_replace(".txt","_log.txt"),
      paste("usearch", cmd), "\n", A,
      append= F, sep = "\n")
  } 
  
} # fin de la fonction

file.remove(list.files(get.value("Compare.OTU.usearch"), full.name = T, pattern =".txt"))
file.remove(list.files(get.value("Compare.OTU.usearch.log"), full.name = T, pattern =".txt"))

# 12S

usearch_global(files = list.files(get.value("filt_derep.path"), pattern = "12s", full.names = T), 
               new.files = list.files(get.value("filt_derep.path"), pattern = "12s", full.names = T) %>% str_replace(get.value("filt_derep.path"),get.value("Compare.OTU.usearch")) %>%  
                 str_replace(".fasta",".txt"), 
               DB = list.files(get.value("OTU.usearch"), pattern = "12s", full.names = T) %>% str_subset("OTU.fasta"))


# CytB.R1

usearch_global(files = list.files(get.value("filt_derep.path"), pattern = "cytB", full.names = T) %>% str_subset("R1"), 
               new.files = list.files(get.value("filt_derep.path"), pattern = "cytB", full.names = T) %>% str_subset("R1") %>% 
                           str_replace(get.value("filt_derep.path"),get.value("Compare.OTU.usearch")) %>%  
                           str_replace(".fasta",".txt"), 
               DB = list.files(get.value("OTU.usearch"), pattern = "cytB.R1", full.names = T) %>% str_subset("OTU.fasta"))

# CytB.R2

usearch_global(files = list.files(get.value("filt_derep.path"), pattern = "cytB", full.names = T) %>% str_subset("R2"), 
               new.files = list.files(get.value("filt_derep.path"), pattern = "cytB", full.names = T) %>% 
                           str_subset("R2") %>% str_replace(get.value("filt_derep.path"),get.value("Compare.OTU.usearch")) %>%  
                           str_replace(".fasta",".txt"), 
               DB = list.files(get.value("OTU.usearch"), pattern = "cytB.R2", full.names = T) %>% str_subset("OTU.fasta"))


# Create a OTU table

make.OTU.table <- function(files, fasta){

   # ne pas prendre les fichiers vides
   empty <- file.info(files)$size==0
   
   files <- files[!empty]
   
   tab <- data.frame(ID = character())
   
   for(x in 1:length(files)){
     
     data <- read.csv(files[x], sep="\t", header = F, stringsAsFactors = F,
                      col.names = c("query", "otu", "ident","length","mism","gap","qstart","qend","target_s","target_e","e.value","bitscore"))
     
     data <- data %>% mutate(abund = sapply(str_split(query, pattern = "="), `[`,2),
                             abund = as.numeric(as.character(abund)))
     
     temp <- data %>% group_by(otu) %>% summarise(Sum = sum(abund))
     names(temp) <- c("ID", files[x] %>% str_remove(get.value("Compare.OTU.usearch")) %>% 
                                         str_remove("/") %>% str_remove(".txt"))
     
     tab <- tab %>% full_join(temp, by = "ID")

   }

   tab[is.na(tab)] <- 0
   
   tab <- tab[mixedorder(tab$ID),]
   tab[,1:2]
   
   #DNA <- readDNAStringSet(fasta)
   
   #row.names(tab) <- as.character(DNA)
   
   return(tab)

}

files1 <- list.files(get.value("Compare.OTU.usearch"), pattern = "12s", full.names = T)
files2 <- list.files(get.value("Compare.OTU.usearch"), pattern = "cytB", full.names = T) %>% str_subset("R1")
files3 <- list.files(get.value("Compare.OTU.usearch"), pattern = "cytB", full.names = T) %>% str_subset("R2")

fasta1 <- list.files(get.value("OTU.usearch"), pattern = "_OTU.fasta", full.names = T) %>% str_subset("12s")
fasta2 <- list.files(get.value("OTU.usearch"), pattern = "_OTU.fasta", full.names = T) %>% str_subset("cytB.R1")
fasta3 <- list.files(get.value("OTU.usearch"), pattern = "_OTU.fasta", full.names = T) %>% str_subset("cytB.R2")


OTUtab.12s     <- make.OTU.table(files1, fasta1)
OTUtab.cytB.R1 <- make.OTU.table(files2, fasta2)
OTUtab.cytB.R2 <- make.OTU.table(files3, fasta3)

save(file = get.value("OTUtable.data"), 
     list = ls(pattern = "OTUtab."))


cat("OTU tables were created:",
    get.value("OTU.usearch"),
    get.value("Compare.OTU.usearch"),
    get.value("OTUtable.data"),
    "\n-------------------------\n", 
    file=get.value("Raw.log"), 
    append = T, sep = "\n") 


# faire la suite pour le cytB

# then save OTU table

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


# Stats - maybe ailleurs --------------------------------------------------



# sum(seqtab.cytB.R)/sum(seqtab.cytB.R.int)
# 
# cat(paste0("12S: Percentage of ASV remaining after chimera removal: ", 
#            round(sum(seqtab.12s)/sum(seqtab.12s.int)*100,1), "%"), 
#     paste0("12S.F: Percentage of ASV remaining after chimera removal: ", 
#            round(sum(seqtab.cytB.F)/sum(seqtab.cytB.F.int)*100,1), "%"), 
#     paste0("12S.R: Percentage of ASV remaining after chimera removal: ", 
#            round(sum(seqtab.cytB.R)/sum(seqtab.cytB.R.int)*100,1), "%"),
#     "\n-------------------------\n",  
#     file=get.value("Raw.log"), 
#     append = T, sep = "\n")
# 
# sum.12S <- get_trackDADA(SUMMARY = filter.12s.summary,
#                          RAW.NAME = all.files[["12s.names"]] %>% str_remove(pattern = "12s-"),
#                          #DADAf = dada.12s.Fs.seta, 
#                          #DADAr = dada.12s.Rs.seta,
#                          MERGER = mergers.12S,
#                          SEQTAB = seqtab.12s,
#                          FILT.NAMES = all.files[["12s.filt.names"]] %>% str_remove(pattern = "12s-")
# )
# 
# sum.12S.F <- get_trackDADA(SUMMARY = filter.12s.summary,
#                            RAW.NAME = all.files[["12s.names"]] %>% str_remove(pattern = "12s-"),
#                            SEQTAB = seqtab.12s.F,
#                            FILT.NAMES = all.files[["12s.filt.names"]] %>% str_remove(pattern = "12s-")
# )
# 
# sum.12S.R <- get_trackDADA(SUMMARY = filter.12s.summary,
#                            RAW.NAME = all.files[["12s.names"]] %>% str_remove(pattern = "12s-"),
#                            SEQTAB = seqtab.12s.R,
#                            FILT.NAMES = all.files[["12s.filt.names"]] %>% str_remove(pattern = "12s-")
# )
# 
# 
# sum.cytB <- get_trackDADA(SUMMARY = filter.cytB.summary,
#                           RAW.NAME = all.files[["cytB.names"]] %>% str_remove(pattern = "cytB-"),
#                           #DADAf = dada.12s.Fs.seta, 
#                           #DADAr = dada.12s.Rs.seta,
#                           MERGER = mergers.cytB,
#                           SEQTAB = seqtab.cytB,
#                           FILT.NAMES = all.files[["cytB.filt.names"]] %>% str_remove(pattern = "cytB-")
# )
# 
# sum.cytB.F <- get_trackDADA(SUMMARY = filter.cytB.summary,
#                             RAW.NAME = all.files[["cytB.names"]] %>% str_remove(pattern = "cytB-"),
#                             SEQTAB = seqtab.cytB.F,
#                             FILT.NAMES = all.files[["cytB.filt.names"]] %>% str_remove(pattern = "cytB-")
# )
# sum.cytB.R <- get_trackDADA(SUMMARY = filter.cytB.summary,
#                             RAW.NAME = all.files[["cytB.names"]] %>% str_remove(pattern = "cytB-"),
#                             SEQTAB = seqtab.cytB.R,
#                             FILT.NAMES = all.files[["cytB.filt.names"]] %>% str_remove(pattern = "cytB-")
# )
# 
# 
# cat("\n 12S:  Process raw summary : \n",  
#     file=file.path(log.path, "Process_RAW.log.txt"), 
#     append = T, sep = "\n")
# 
# 
# write.table(sum.12S, 
#             file=file.path(log.path, "Process_RAW.log.txt"),  
#             row.names=T, col.names=T, append = T)
# 
# 
# cat("\n cytB:  Process raw summary : \n",  
#     file=file.path(log.path, "Process_RAW.log.txt"), 
#     append = T, sep = "\n")
# 
# 
# write.table(sum.cytB, 
#             file=file.path(log.path, "Process_RAW.log.txt"),  
#             row.names=T, col.names=T, append = T)
# 
# 
# plot(sum.cytB.F$nonchim, sum.cytB.R$nonchim)
# plot(sum.12S.F$nonchim, sum.12S.R$nonchim)
# 
# plot(sum.12S.R$nonchim, sum.12S$nonchim)
# plot(sum.12S.F$nonchim, sum.12S$nonchim)
# 
# plot(sum.cytB.F$nonchim, sum.cytB$nonchim)
# plot(sum.cytB.R$nonchim, sum.cytB$nonchim)
# 
# sum.12S.wInfo  <- sum.12S  %>% left_join(DataSeq %>% select(SampleID, SeqType, IbisID, NomLac, CatSite), by = c("sample" = "IbisID"))
# sum.cytB.wInfo <- sum.cytB %>% left_join(DataSeq %>% select(SampleID, SeqType, IbisID, NomLac, CatSite), by = c("sample" = "IbisID"))
# 
# 
# sum.12S.wInfo.res <- sum.12S.wInfo %>%  filter(SeqType %in% c("sample", "mix", "blank")) %>% 
#   group_by(SeqType) %>% 
#   summarise(perc.all.filt = round(sum(filtered)/sum(nreads),3),
#             perc.filt.merged = round(sum(merged)/sum(filtered),3),
#             perc.merged.nochim = round(sum(nonchim)/sum(merged),3),
#             perc.all.nochim  = round(sum(nonchim)/sum(nreads),3)
#   ) %>% as.data.frame()
# 
# 
# sum.cytB.wInfo.res <- sum.cytB.wInfo %>%  
#   filter(SeqType %in% c("sample", "mix", "blank")) %>% 
#   group_by(SeqType) %>% 
#   summarise(perc.all.filt = round(sum(filtered)/sum(nreads),3),
#             perc.filt.merged = round(sum(merged)/sum(filtered),3),
#             perc.merged.nochim = round(sum(nonchim)/sum(merged),3),
#             perc.all.nochim  = round(sum(nonchim)/sum(nreads),3)
#   ) %>% 
#   as.data.frame()
# 
# cat("\n 12S:  Process raw summary by SeqType: \n",  
#     file=file.path(log.path, "Process_RAW.log.txt"), 
#     append = T, sep = "\n")
# 
# 
# write.table(sum.12S.wInfo.res, 
#             file=file.path(log.path, "Process_RAW.log.txt"),  
#             row.names=F, col.names=T, append = T)
# 
# cat("\n cytB:  Process raw summary by SeqType: \n",  
#     file=file.path(log.path, "Process_RAW.log.txt"), 
#     append = T, sep = "\n")
# 
# 
# write.table(sum.cytB.wInfo.res, 
#             file=file.path(log.path, "Process_RAW.log.txt"),  
#             row.names=F, col.names=T, append = T)
# 
# cat("\n-------------------------\n",  
#     file=file.path(log.path, "Process_RAW.log.txt"), 
#     append = T, sep = "\n")
# 
# 
# 
# 
# 

# WHAT I WILL DO

#ShortRead::countLines(get.value("filt_dada2.path"), pattern = "12s_Sample_Che")/4



# END of the script -------------------------------------------------------

# save.image(file.path(log.path,"Process_RAW.Rdata"))

end.time <- round(Sys.time() - start.time,2) 

cat("\nEND of the raw data processing!",
    #paste0("Rdata saved: ", file.path(log.path,"Process_RAW.Rdata")),
    paste0("\nTime to run this script: ", end.time, " ",units(end.time)),
    "\n-------------------------\n", 
   
    paste("R version", devtools::session_info()$platform$version, sep = ": "),
    paste("OS", devtools::session_info()$platform$os, sep = ": "),
    paste("System", devtools::session_info()$platform$system, sep = ": "),    
    
    "\n~ Important R packages ~",
    paste("dada2", packageVersion("dada2"), sep = ": "),     
    paste("fastqcr", packageVersion("fastqcr"), sep = ": "),     
    paste("Biostrings", packageVersion("Biostrings"), sep = ": "),   

    "\n~ External programs ~",
    
    
    # Add it to the log file
    file=get.value("Raw.log"), 
    append = T, sep = "\n")


if(get_os() %in% c("os","linux")){ # to run only on os and not windows
  
  cat(paste("fastqc", system2("fastqc", "-v", stdout=T, stderr=T), sep = ": "),     
      paste("cutadapt", system2("cutadapt", "--version", stdout=T, stderr=T), sep = ": "),     
      paste("vsearch", system2("vsearch", "-v", stdout=T, stderr=T)[1] , sep = ": "),   
      paste("usearch", system2("usearch", "--version", stdout=T, stderr=T) , sep = ": "),      
      "\n-------------------------\n", 
            # Add it to the log file
      file=get.value("Raw.log"), 
      append = T, sep = "\n")


}

# END OF THE SCRIPT

