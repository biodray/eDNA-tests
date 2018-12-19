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
mixed
# Fastq and fasta manipulation
library(Biostrings)


# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}


# Data --------------------------------------------------------------------


load(get.value("ASVtable.data"))
load(get.value("OTUtable.data"))

DataSample <- read_excel(get.value("Sample.xl"),sheet="DataSample",na="NA",guess_max=100000)
DataSeq    <- read_excel(get.value("Sample.xl"),sheet="DataSeq",na="NA",guess_max=100000)
#Amorces    <- read_excel(get.value("Sample.xl"),sheet="Amorces",na="NA",guess_max=100000)
#Inventaire <- read_excel(get.value("Sample.xl"),sheet="DataLac",na="NA",guess_max=100000)

# Format de DataSeq avec les données IBIS (p1-A1)

DataSeq <- DataSeq %>% mutate (IbisID = paste0(Plaque,".",Puit)) %>% 
  left_join(DataSample %>% select(SampleID, NomLac, CatSite, Temoins), by = "SampleID")


# Transform dada2 to a better format

SEQtable.df <- function(tab){
  
  new.tab <- data.frame(t(tab))
  new.tab$ID <- paste0("ASV_", 1:nrow(new.tab))
  
  return(new.tab)
}

ASVtab.12s     <- SEQtable.df(ASVtab.12s)
ASVtab.cytB.R1 <- SEQtable.df(ASVtab.cytB.F)
ASVtab.cytB.R2 <- SEQtable.df(ASVtab.cytB.R)


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

names(ASVtab.12s)     <- names(ASVtab.12s) %>% simplify.col() %>% str_remove("_R1")
names(ASVtab.cytB.R1) <- names(ASVtab.cytB.R1) %>% simplify.col()
names(ASVtab.cytB.R2) <- names(ASVtab.cytB.R2) %>% simplify.col()

names(OTUtab.12s)     <- names(OTUtab.12s) %>% simplify.col()
names(OTUtab.cytB.R1) <- names(OTUtab.cytB.R1) %>% simplify.col()
names(OTUtab.cytB.R2) <- names(OTUtab.cytB.R2) %>% simplify.col()




# add a column for each sample


add.missing <- function(tab, locus, sens = NULL){

if(locus == "12s")  all.files <- list.files(get.value("raw_unz_rename.path"), pattern = "R1")  %>% str_remove("_R1.fastq") %>% str_subset(locus) %>% str_replace("-", ".")
if(locus == "cytB") all.files <- list.files(get.value("raw_unz_rename.path"), pattern = sens)  %>% str_remove(".fastq") %>% str_subset(locus) %>% str_replace("-", ".")

missing.files <- setdiff(all.files, names(tab))

for(x in missing.files){
  tab[,x] <- 0
  
}

return(tab)

}


ASVtab.12s     <- add.missing(ASVtab.12s, "12s")
ASVtab.cytB.R1 <- add.missing(ASVtab.cytB.R1, "cytB", "R1")
ASVtab.cytB.R2 <- add.missing(ASVtab.cytB.R2, "cytB", "R2")

OTUtab.12s     <- add.missing(OTUtab.12s, "12s")
OTUtab.cytB.R1 <- add.missing(OTUtab.cytB.R1, "cytB", "R1")
OTUtab.cytB.R2 <- add.missing(OTUtab.cytB.R2, "cytB", "R2")



# Mix <- names(ASVtab.12s) %>% str_subset("Sample")
# Sample <- names(ASVtab.12s) %>% str_subset("Sample")
# Tneg   <- c(names(ASVtab.12s) %>% str_subset("Tneg"),
#             names(ASVtab.12s) %>% str_subset("T0"),
#             names(ASVtab.12s) %>% str_subset("T1"))
# 
# ASVtab.12s %>% select("ID",Tneg) %>%
#                gather(Tneg, key = "sample", value = "N") %>% 
#                       filter(N>=1)  %>% 
#                ggplot(aes(x = sample, y = ID, fill = N)) + 
#                       geom_bin2d() + 
#                       scale_fill_distiller(palette = "Spectral") +
#                       scale_y_discrete("ASV", limits=mixedsort(ASVtab.12s$ID), labels = NULL) +
#                       theme(axis.text.x = element_text(angle = 90, hjust = 1))


rapid.graph <- function(tab, Sample = T, Tneg = T, Mix = T){

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
                      scale_y_discrete("Sequences", limits=mixedsort(tab$ID), labels = NULL) +
                      theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(graph)

}


rapid.graph(ASVtab.12s)

rapid.graph(ASVtab.cytB.R1)
rapid.graph(ASVtab.cytB.R2)

rapid.graph(OTUtab.12s)
rapid.graph(OTUtab.cytB.R1)
rapid.graph(OTUtab.cytB.R2)

  
# ?scale_y_discrete
# 
# Sample <- names(OTUtab.12s) %>% str_subset("Sample")
# 
# OTUtab.12s %>% select("ID",Sample) %>%
#                gather(Sample, key = "sample", value = "N") %>% 
#                ggplot(aes(x = sample, y = ID, fill = N)) + 
#                geom_tile() + 
#                scale_fill_distiller(palette = "Spectral")
# 
# 

# PLaque

plaque.graph <- function(tab, Sample = T, Tneg = T, Mix = T){
  
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
    scale_y_discrete("",limits=LETTERS[8:1]) +
    scale_x_discrete("",limits=c(1:12)) + 
    geom_point(aes(shape = neg)) +
      facet_grid(plaque~.) +
      theme_bw()
    #theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  print(graph)
  
}

plaque.graph(ASVtab.12s)

plaque.graph(ASVtab.cytB.R1)

plaque.graph(ASVtab.cytB.R2)


head(ASVtab.cytB.R2[,1:3])


# Commencer différentes tactiques pour enlever des reads - soit le N max, ou encore le N associé à chaque journée

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


plaque.graph(ASVtab.12s.cor)
plaque.graph(ASVtab.cytB.R1.cor)
plaque.graph(ASVtab.cytB.R2.cor)

plaque.graph(OTUtab.12s.cor)
plaque.graph(OTUtab.cytB.R1.cor)
plaque.graph(OTUtab.cytB.R2.cor) 



rapid.graph(ASVtab.12s.cor)
rapid.graph(ASVtab.cytB.R1.cor)
rapid.graph(ASVtab.cytB.R2.cor)

rapid.graph(OTUtab.12s.cor)
rapid.graph(OTUtab.cytB.R1.cor)
rapid.graph(OTUtab.cytB.R2.cor) 
