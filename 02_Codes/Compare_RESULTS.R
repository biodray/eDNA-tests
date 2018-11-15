
# Info --------------------------------------------------------------------

# Compare the differents results
# 
# Audrey Bourret
# 2018-11-14
#


# Library -----------------------------------------------------------------

library(readxl)
library(tidyverse)

# Internal functions
source(file.path("./03_Functions",  list.files("./03_Functions")))


# Data --------------------------------------------------------------------

load(file = file.path(get.value("result.path"), "Compil.data"))
load(file = file.path(get.value("result.path"), "Taxo.data"))

Sample.xl <-  file.path(get.value("info.path"), get.value("Sample.xl"))

DataSample <- read_excel(Sample.xl,sheet="DataSample",na="NA",guess_max=100000)
DataSeq    <- read_excel(Sample.xl,sheet="DataSeq",na="NA",guess_max=100000)
Amorces    <- read_excel(Sample.xl,sheet="Amorces",na="NA",guess_max=100000)
Inventaire <- read_excel(Sample.xl,sheet="DataLac",na="NA",guess_max=100000)

# Format de DataSeq avec les données IBIS (p1-A1)

DataSeq <- DataSeq %>% mutate (IbisID = paste0("p",Plaque,"-",Puit)) %>% 
  left_join(DataSample %>% select(SampleID, NomLac, CatSite), by = "SampleID")

DataSeq

# Functions ---------------------------------------------------------------

count.taxon <- function(TAXO, CLASS = c("Class", "Order", "Family", "Genus", "Species")){
  res1<- vector()
  res2<- vector()
  for(x in 1:length(CLASS)){
    res1[x] <- nrow(TAXO[which(!is.na(TAXO[,CLASS[x]])),])
    res2[x] <- length(unique(na.omit(TAXO[,CLASS[x]])))
  }
    return(cbind(res1, res2))
}

diff.taxon <- function(TAXO1, TAXO2, CLASS = c("Class", "Order", "Family", "Genus", "Species")){
  res1<- vector()
  res2 <- vector()
  for(x in 1:length(CLASS)){
    df <- cbind(TAXO1[,CLASS[x]], TAXO2[,CLASS[x]])
    df <- df[which(!is.na(df[,1]) & !is.na(df[,2])),]
    
    res1[x] <- nrow(df)
    
    df2 <- df[which(df[,1] != df[,2]),]
    
    #row.names(df2) <- NULL
    #print(df2)
    
    res2[x] <- nrow(df2)    
    
  }
  
  return(cbind(res1,res2))
}

RDPvsIDT <- function(RDP, IDT, CLASS = c("Class", "Order", "Family", "Genus", "Species")) {
  res <- data.frame(taxon = CLASS, RDP = NA, RDP.uniq = NA, IDT = NA, IDT.uniq = NA, Both = NA, Diff = NA)
  
  res[,c("RDP", "RDP.uniq")] <- count.taxon(RDP)
  res[,c("IDT", "IDT.uniq")] <- count.taxon(IDT)
  
  res[,c("Both", "Diff")] <- diff.taxon(RDP, IDT)
  
  res$RDP.p <- round(res$RDP / nrow(RDP),3) 
  res$IDT.p <- round(res$IDT / nrow(IDT),3) 
  
  return(res) 

}




#

add.info <- function(OBJ = ls(pattern = "compil.R|Compil.D"), INFO = DataSeq){

  res <- list()
  
  for(x in OBJ){
    DATA <- get(x)
    if(str_detect(x, "RDP")){
      DATA <- DATA[[2]]
    }
    
    DATA$Sample <- paste(sapply(str_split(DATA$Sample, "-"), `[`, 2),
                         sapply(str_split(DATA$Sample, "-"), `[`, 3),
                         sep = "-")
    
    DATA <-  DATA %>% left_join(INFO %>% select(SeqType,IbisID,NomLac,CatSite), by = c("Sample" = "IbisID"))

    res[[x]] <- DATA
  }

  return(res)
}




# Comparison --------------------------------------------------------------


compil.wInfo <- add.info(ls(pattern = "compil.R|Compil.D"), DataSeq)

names(compil.wInfo)


# Assign species

# Check within F vs R vs merged

View(taxo.RDP.12S)
View(taxo.IDT.12S[[2]])




count.taxon(taxo.RDP.12S)

RDPvsIDT(taxo.RDP.12S, taxo.IDT.12S[[2]])
RDPvsIDT(taxo.RDP.12S.IBIS, taxo.IDT.12S.IBIS[[2]])

RDPvsIDT(taxo.RDP.12S.F, taxo.IDT.12S.F[[2]])
RDPvsIDT(taxo.RDP.12S.R, taxo.IDT.12S.R[[2]])

RDPvsIDT(taxo.RDP.cytB.F, taxo.IDT.cytB.F[[2]])
RDPvsIDT(taxo.RDP.cytB.R, taxo.IDT.cytB.R[[2]])

# Check across datasets with the same methods (Biodiv)


