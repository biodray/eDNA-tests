# Info --------------------------------------------------------------------

# Correct OTU and ASV table based on negative samples
# 
# Audrey Bourret
# 2018-12-14
#

#start.time <- Sys.time()

# Library -----------------------------------------------------------------

# Data manipulation

library(tidyverse) # includes ggplot2 dplyr and stringr

library(gtools)    # for mixedsort
library(readxl)


#library(devtools)
#devtools::install_github("kassambara/ggpubr")
library(ggpubr)    # on github - for nice graphs

# Fastq and fasta manipulation
#library(Biostrings)


# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}


# Functions ---------------------------------------------------------------

# Function to transform dada2 table on the other side 
SEQtable.df <- function(tab){
  
  new.tab <- data.frame(t(tab))
  new.tab$ID <- paste0("ASV_", 1:nrow(new.tab))
  
  return(new.tab)
}

# Function to remove superflues on a vector of names
simplify.col <- function(x){
  
  x %>% str_remove(".fastq") %>% 
    str_remove("_cut") %>% 
    str_remove("_filt") %>% 
    str_remove("_merged") %>% 
    str_remove("_min") %>% 
    str_remove("_derep") %>% 
    str_replace("X12s", "12s") %>% 
    str_replace("-", ".")
  
}

# Function to add a column to missing values
add.missing <- function(tab, locus, sens = NULL){
  
  if(locus == "12s")  all.files <- list.files(get.value("raw_unz_rename.path"), pattern = "R1")  %>% str_remove("_R1.fastq") %>% str_subset(locus) %>% str_replace("-", ".")
  if(locus == "cytB") all.files <- list.files(get.value("raw_unz_rename.path"), pattern = sens)  %>% str_remove(".fastq") %>% str_subset(locus) %>% str_replace("-", ".")
  
  missing.files <- setdiff(all.files, names(tab))
  
  for(x in missing.files){
    tab[,x] <- 0
    
  }
  
  return(tab)
  
}


# Function to create a heatmap for each sample from a SEQtab
rapid.graph <- function(tab, Sample = T, Tneg = T, Mix = T, maintitle = "Heatmap of sequence frequency"){
  
  Sample1 <- vector()
  Tneg1   <- vector()
  Mix1    <- vector()  
  
  Sample1 <- if(isTRUE(Sample)) names(tab) %>% str_subset("Sample")
  
  Tneg1   <- if(isTRUE(Tneg)) {
    c(names(tab) %>% str_subset("Tneg"),
      names(tab) %>% str_subset("T0"),
      names(tab) %>% str_subset("T1"))
  }
  
  
  Mix1    <- if(isTRUE(Mix)) names(tab) %>% str_subset("Mix")
  
  graph <- tab %>% select("ID",Sample1, Tneg1, Mix1) %>%
    gather(Sample1, Tneg1, Mix1, key = "sample", value = "N") %>% 
    filter(N>=1)  %>% 
    ggplot(aes(x = sample, y = ID, fill = N)) + 
    geom_bin2d() + 
    scale_fill_distiller(palette = "Spectral", trans = "log10") +
    #scale_fill_gradient(low = "darkgray", high = "red", trans = "log") +
    scale_y_discrete(limits=mixedsort(tab$ID), labels = NULL) +
    labs(title= maintitle, x ="Sample", y = "Sequence ID") +
    guides(fill = guide_colourbar(title = "N reads", title.hjust = 0)) + 
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.ticks.y = element_blank())
  
  print(graph)
  
}


# Function to create a heatmap from a SEQtab
plaque.graph <- function(tab, Sample = T, Tneg = T, Mix = T, maintitle = "PCR plate"){
  
  Sample1 <- vector()
  Tneg1   <- vector()
  Mix1    <- vector()  
  
  Sample1 <- if(isTRUE(Sample)) names(tab) %>% str_subset("Sample")
  
  Tneg1   <- if(isTRUE(Tneg)) {
    c(names(tab) %>% str_subset("Tneg"),
      names(tab) %>% str_subset("T0"),
      names(tab) %>% str_subset("T1"))
  }
  
  
  Mix1    <- if(isTRUE(Mix)) names(tab) %>% str_subset("Mix")
  
  #tab <- ASVtab.12s
  
  graph <- 
    
    tab %>% select("ID",Sample1, Tneg1, Mix1) %>%
    gather(Sample1, Tneg1, Mix1, key = "sample", value = "N") %>% 
    group_by(sample) %>% 
    summarise(Sum = sum(N)) %>% 
    mutate(puit = sapply(str_split(sample, "_p"),`[`,2) %>% str_remove("_R1") %>% str_remove("_R2"),
           plaque = str_sub(puit, 1, 1),
           rn = str_sub(puit, 3, 3),
           cn = str_sub(puit, 4,5),
           neg = ifelse(str_detect(sample, pattern = c("Tneg")), "neg.lab", 
                        ifelse(str_detect(sample, pattern = c("T0")), "neg.field",
                               ifelse(str_detect(sample, pattern = c("T1")), "neg.field", NA)))) %>% #View()
    
    
    #filter(plaque == "2")  %>% 
    ggplot(aes(x = cn, y = rn, fill = Sum)) + 
    geom_bin2d() + 
    scale_fill_distiller(palette = "Spectral", trans = "log10") +#+
    scale_shape_discrete(limits= c("neg.field", "neg.lab"), labels = c("field", "lab"))+
    scale_y_discrete("",limits=LETTERS[8:1]) +
    scale_x_discrete("",limits=c(1:12), position = "top") + 
    guides(fill = guide_colourbar(title = "N reads", title.hjust = 0),
           shape = guide_legend("Negative controls")) + 
    labs(title= maintitle) +
    geom_point(aes(shape = neg)) +
    facet_grid(plaque~., labeller = labeller(plaque = label_both)) +
    theme_bw() +
    theme(strip.background = element_rect(colour = "black", fill = "white"))
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  print(graph)
  
}

# Data --------------------------------------------------------------------

load(get.value("ASVtable.data"))
load(get.value("OTUtable.data"))

DataSample <- read_excel(get.value("Sample.xl"),sheet="DataSample",na="NA",guess_max=100000)
DataSeq    <- read_excel(get.value("Sample.xl"),sheet="DataSeq",na="NA",guess_max=100000)
#Amorces    <- read_excel(get.value("Sample.xl"),sheet="Amorces",na="NA",guess_max=100000)
#Inventaire <- read_excel(get.value("Sample.xl"),sheet="DataLac",na="NA",guess_max=100000)

# Format de DataSeq avec les données IBIS (1.1)

DataSeq <- DataSeq %>% mutate (IbisID = paste0(Plaque,".",Puit)) %>% 
  left_join(DataSample %>% select(SampleID, NomLac, CatSite, Temoins), by = "SampleID")


# Data transformation -----------------------------------------------------

# Transform dada2 results as OTU result
ASVtab.12s     <- SEQtable.df(ASVtab.12s)
ASVtab.cytB.R1 <- SEQtable.df(ASVtab.cytB.F)
ASVtab.cytB.R2 <- SEQtable.df(ASVtab.cytB.R)


# Simplifyin column names
names(ASVtab.12s)     <- names(ASVtab.12s) %>% simplify.col() %>% str_remove("_R1")
names(ASVtab.cytB.R1) <- names(ASVtab.cytB.R1) %>% simplify.col()
names(ASVtab.cytB.R2) <- names(ASVtab.cytB.R2) %>% simplify.col()

names(OTUtab.12s)     <- names(OTUtab.12s) %>% simplify.col()
names(OTUtab.cytB.R1) <- names(OTUtab.cytB.R1) %>% simplify.col()
names(OTUtab.cytB.R2) <- names(OTUtab.cytB.R2) %>% simplify.col()


# Add a column for each missing samples

ASVtab.12s     <- add.missing(ASVtab.12s, "12s")
ASVtab.cytB.R1 <- add.missing(ASVtab.cytB.R1, "cytB", "R1")
ASVtab.cytB.R2 <- add.missing(ASVtab.cytB.R2, "cytB", "R2")

OTUtab.12s     <- add.missing(OTUtab.12s, "12s")
OTUtab.cytB.R1 <- add.missing(OTUtab.cytB.R1, "cytB", "R1")
OTUtab.cytB.R2 <- add.missing(OTUtab.cytB.R2, "cytB", "R2")

# Make cool graphs

pdf(file.path(get.value("result.OTUtables"),"SEQheatmap.raw.pdf"), width = 11, height = 8.5)
  
  rapid.graph(ASVtab.12s, maintitle = "Heatmap of ASV frequency - 12S on raw data")
  rapid.graph(ASVtab.cytB.R1, maintitle = "Heatmap of ASV frequency - cytB.R1 on raw data")
  rapid.graph(ASVtab.cytB.R2, maintitle = "Heatmap of ASV frequency - cytB.R2 on raw data")
  
  rapid.graph(OTUtab.12s, maintitle = "Heatmap of OTU frequency - 12S on raw data")
  rapid.graph(OTUtab.cytB.R1, maintitle = "Heatmap of OTU frequency - cytB.R1 on raw data")
  rapid.graph(OTUtab.cytB.R2, maintitle = "Heatmap of OTU frequency - cytB.R2 on raw data")

dev.off()

# Plaque PCR

pdf(file.path(get.value("result.OTUtables"),"PCRheatmap.raw.pdf"), width = 11, height = 8.5)

  plaque.graph(ASVtab.12s, maintitle = "PCR plates heatmap - 12S on ASV raw data")
  plaque.graph(ASVtab.cytB.R1, maintitle = "PCR plates heatmap - cytB.R1 on ASV raw data")
  plaque.graph(ASVtab.cytB.R2, maintitle = "PCR plates heatmap - cytB.R2 on ASV raw data")
  
  plaque.graph(OTUtab.12s, maintitle = "PCR plates heatmap - 12S on OTU raw data")
  plaque.graph(OTUtab.cytB.R1, maintitle = "PCR plates heatmap - cytB.R1 on OTU raw data")
  plaque.graph(OTUtab.cytB.R2, maintitle = "PCR plates heatmap - cytB.R2 on OTU raw data")

dev.off()


# Study on negative values ------------------------------------------------

extract.neg <- function(tab){

  Tnega    <- names(tab) %>% str_subset("Tneg")
  Tfield  <- c(names(tab) %>% str_subset("T0"),
               names(tab) %>% str_subset("T1"))
  
  DATA <- tab %>% select("ID", Tnega, Tfield) %>%
                  gather(Tnega,Tfield, key = "sample", value = "N") %>% 
                  group_by(sample) %>% 
                  summarise(Nreads = sum(N),
                            Nseq = length(N[which(N>0)]))
  
  return(DATA)

}

extract.sample <- function(tab){
  
  Sample <- names(tab) %>% str_subset("Sample")
  
  DATA <- tab %>% select("ID", Sample) %>%
    gather(Sample, key = "sample", value = "N") %>% 
    group_by(sample) %>% 
    summarise(Nreads = sum(N),
              Nseq = length(N[which(N>0)]))
  
  return(DATA)
  
}

compare.neg <- function(tab1, tab2, N, maintitle = "Comparison of ASV and OTU"){
  DATA <- rbind(extract.sample(tab1) %>% mutate(data = "ASV") %>%  mutate(type = "sample"),
                     extract.neg(tab1) %>% mutate(data = "ASV") %>%  mutate(type = "negative"),
                     extract.sample(tab2) %>% mutate(data = "OTU") %>%  mutate(type = "sample"),
                     extract.neg(tab2) %>% mutate(data = "OTU")%>%  mutate(type = "negative")
                     )

  if(N == "Nreads"){
    
    GRAPH <- DATA %>% ggplot(aes(y = Nreads, x = data, col = type)) 
                   
  } else if(N == "Nseq"){
      
    GRAPH <- DATA %>% ggplot(aes(y = Nseq, x = data, col = type))  
        
  } else(stop()) 
    
  GRAPH <- GRAPH + geom_boxplot()+
                   labs(title= maintitle) +
                   theme_bw()
    
  print(GRAPH)

}

new.graph <-ggarrange(compare.neg(ASVtab.12s, OTUtab.12s, "Nreads", maintitle = "12S -  N reads"),
                      compare.neg(ASVtab.12s, OTUtab.12s, "Nseq", maintitle = "12S - N sequences"),
                      compare.neg(ASVtab.cytB.R1, OTUtab.cytB.R1, "Nreads", maintitle = "cytB.R1 - N reads"),
                      compare.neg(ASVtab.cytB.R1, OTUtab.cytB.R1, "Nseq", maintitle = "cytB.R1 - N sequences"),
                      compare.neg(ASVtab.cytB.R2, OTUtab.cytB.R2, "Nreads", maintitle = "cytB.R2 - N reads"),
                      compare.neg(ASVtab.cytB.R2, OTUtab.cytB.R2, "Nseq", maintitle = "cytB.R2 - N sequences"),
                      labels = LETTERS[1:6],
                      ncol = 2, nrow = 3,
                      common.legend = TRUE, legend = "right")

pdf(file.path(get.value("result.OTUtables"),"ComparisonASV-OTU.pdf"), width = 8.5, height = 11)

  new.graph 

dev.off()


# Corrections N reads for samples -----------------------------------------


correct.reads <- function(tab){ 

  Sample <- names(tab) %>% str_subset("Sample")
  Tnega    <- names(tab) %>% str_subset("Tneg")
  Tfield  <- c(names(tab) %>% str_subset("T0"),
               names(tab) %>% str_subset("T1"))
  Mix    <- names(tab) %>% str_subset("Mix")
  
  # Jeux de données liés
  
  DATA <- tab %>% select("ID", Sample, Tnega, Tfield, Mix) %>%
          gather(Sample, Tnega, Tfield, Mix, key = "sample", value = "N") %>% 
          mutate(IbisID = sapply(str_split(sample, "_p"),`[`,2) %>% str_remove("_R1") %>% str_remove("_R2") %>% str_replace("-", ".")) %>%  
          left_join(DataSeq %>% select(IbisID, SampleID, SeqType, Tneg, Temoins), by = "IbisID") 
  
  DATA.Tneg   <- DATA %>% filter(sample %in% c(Tnega)) %>% select(ID, N, SampleID)
  DATA.Tfield <- DATA %>% filter(sample %in% c(Tfield)) %>% select(ID, N, SampleID)
    
  DATA <- DATA %>% left_join(DATA.Tneg, by = c("ID" = "ID", "Tneg" = "SampleID" ), suffix = c("", ".Tneg")) %>% 
           left_join(DATA.Tfield, by = c("ID" = "ID", "Temoins" = "SampleID" ), suffix = c("", ".Tfield")) %>% 
           mutate(Ncor = ifelse(N.Tfield - N.Tneg >0, N - (N.Tneg + (N.Tfield - N.Tneg)), N - N.Tneg),
                  Ncor = ifelse(Ncor < 0, 0 , Ncor)) %>% 
           filter(SeqType %in% c("sample", "dup.sample")) %>%  
           select(ID, sample, Ncor) %>%
           spread(key = sample, value = Ncor)
            
  return(DATA)  
  
}  
  

ASVtab.12s.cor     <- correct.reads(ASVtab.12s)  
ASVtab.cytB.R1.cor <- correct.reads(ASVtab.cytB.R1)  
ASVtab.cytB.R2.cor <- correct.reads(ASVtab.cytB.R2) 

OTUtab.12s.cor     <- correct.reads(OTUtab.12s)  
OTUtab.cytB.R1.cor <- correct.reads(OTUtab.cytB.R1)  
OTUtab.cytB.R2.cor <- correct.reads(OTUtab.cytB.R2) 

# Representation graphique

pdf(file.path(get.value("result.OTUtables"),"SEQheatmap.corrected.pdf"), width = 11, height = 8.5)

  rapid.graph(ASVtab.12s.cor, maintitle = "Heatmap of ASV frequency - 12S on corrected data")
  rapid.graph(ASVtab.cytB.R1.cor, maintitle = "Heatmap of ASV frequency - cytB.R1 on corrected data")
  rapid.graph(ASVtab.cytB.R2.cor, maintitle = "Heatmap of ASV frequency - cytB.R2 on corrected data")
  
  rapid.graph(OTUtab.12s.cor, maintitle = "Heatmap of OTU frequency - 12S on corrected data")
  rapid.graph(OTUtab.cytB.R1.cor, maintitle = "Heatmap of OTU frequency - cytB.R1 on corrected data")
  rapid.graph(OTUtab.cytB.R2.cor, maintitle = "Heatmap of OTU frequency - cytB.R2 on corrected data")

dev.off()

# Plaque PCR

pdf(file.path(get.value("result.OTUtables"),"PCRheatmap.corrected.pdf"), width = 11, height = 8.5)

  plaque.graph(ASVtab.12s.cor, maintitle = "PCR plates heatmap - 12S on ASV corrected data")
  plaque.graph(ASVtab.cytB.R1.cor, maintitle = "PCR plates heatmap - cytB.R1 on ASV corrected data")
  plaque.graph(ASVtab.cytB.R2.cor, maintitle = "PCR plates heatmap - cytB.R2 on ASV corrected data")
  
  plaque.graph(OTUtab.12s.cor, maintitle = "PCR plates heatmap - 12S on OTU corrected data")
  plaque.graph(OTUtab.cytB.R1.cor, maintitle = "PCR plates heatmap - cytB.R1 on OTU corrected data")
  plaque.graph(OTUtab.cytB.R2.cor, maintitle = "PCR plates heatmap - cytB.R2 on OTU corrected data")

dev.off()

# Enregistrer les SEQTAB corrigées
#save()

# The code for mock community should be on another space? Maybe only after assigning species




