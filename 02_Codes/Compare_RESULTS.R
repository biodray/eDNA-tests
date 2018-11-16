
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
load(file = file.path(get.value("result.path"), "Seqtab.data"))
load(file = file.path(get.value("result.path"), "Compil.data"))
load(file = file.path(get.value("result.path"), "Taxo.data"))

Sample.xl <-  file.path(get.value("info.path"), get.value("Sample.xl"))

DataSample <- read_excel(Sample.xl,sheet="DataSample",na="NA",guess_max=100000)
DataSeq    <- read_excel(Sample.xl,sheet="DataSeq",na="NA",guess_max=100000)

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

add.info <- function(OBJ = ls(pattern = "compil.R|compil.I"), INFO = DataSeq, GATHER = FALSE){

  res <- list()
  
  for(x in OBJ){
    LIST <- get(x)
    
    CLASS <- names(LIST)
    
    for(c in CLASS){
      
      DATA <- LIST[[c]]
    
      DATA$Sample <- paste(sapply(str_split(DATA$Sample, "-"), `[`, 2),
                           sapply(str_split(DATA$Sample, "-"), `[`, 3),
                           sep = "-")
    
      DATA <-  DATA %>% left_join(INFO %>% select(SeqType,SampleID,IbisID,NomLac,CatSite), by = c("Sample" = "IbisID"))

      if(GATHER == TRUE){
    
        GROUP <- names(DATA)[2:(length(names(DATA))-4)]
    
        DATA <- DATA %>% gather(GROUP, key = "groups", value = "N")
    
      }

    res[[paste(x,c, sep=".")]] <- DATA
    }
  }

  return(res)
}




# Comparison --------------------------------------------------------------

# Add sample info

compil.wInfo <- add.info(ls(pattern = "compil.R|compil.I"), DataSeq, GATHER = TRUE)

names(compil.wInfo)


compil.wInfo[[1]]


# fait dans RMD
for(x in names(compil.wInfo)){
  
  cat("\n\n", x, "\n\n")
  
  compil.wInfo[[x]] %>% filter(SeqType %in% c("sample", "mix", "blank")) %>% 
    group_by(SeqType, groups) %>% 
    summarise(Sum = mean(N)) %>% 
    spread(key = SeqType, value = Sum) %>% print()
  
}



# Check within F vs R vs merged


RDPvsIDT(taxo.RDP.12S, taxo.IDT.12S[[2]])
RDPvsIDT(taxo.RDP.12S.IBIS, taxo.IDT.12S.IBIS[[2]])

RDPvsIDT(taxo.RDP.12S.F, taxo.IDT.12S.F[[2]])
RDPvsIDT(taxo.RDP.12S.R, taxo.IDT.12S.R[[2]])

RDPvsIDT(taxo.RDP.cytB.F, taxo.IDT.cytB.F[[2]])
RDPvsIDT(taxo.RDP.cytB.R, taxo.IDT.cytB.R[[2]])



class.diff <- function(TAXO1, TAXO2, CLASS = c("Order", "Family", "Genus", "Species")){
    
  for(x in CLASS){
    cat("\n",x,":","\n\n", sep = "")
    
    
    if(length( setdiff(TAXO1[,x],  TAXO2[,x])>=1)){
      cat("Only in set1: ", setdiff(TAXO1[,x],  TAXO2[,x]), "\n", sep = " ")      
    } else (cat("None only in set1!\n"))
    
    if(length( setdiff(TAXO2[,x],  TAXO1[,x])>=1)){
      cat("Only in set2: ", setdiff(TAXO2[,x],  TAXO1[,x]), "\n", sep = " ")      
    } else (cat("None only in set2!\n"))
    

  }



}


class.diff(taxo.RDP.12S, taxo.IDT.12S[[2]])

class.diff(taxo.RDP.12S, taxo.RDP.12S.F)
class.diff(taxo.RDP.12S, taxo.RDP.12S.R)
class.diff(taxo.RDP.12S.F, taxo.RDP.12S.R)

class.diff( taxo.IDT.12S[[2]], taxo.IDT.12S.F[[2]])
class.diff( taxo.IDT.12S[[2]], taxo.IDT.12S.R[[2]])
class.diff( taxo.IDT.12S.F[[2]], taxo.IDT.12S.R[[2]])

# Ibis
class.diff(taxo.RDP.12S, taxo.RDP.12S.IBIS)
class.diff(taxo.IDT.12S[[2]], taxo.IDT.12S.IBIS[[2]])

# 12S vs cytB

class.diff(taxo.IDT.12S[[2]], taxo.IDT.cytB[[2]])
class.diff(taxo.RDP.12S, taxo.RDP.cytB.F)

# Check across datasets with the same methods (Biodiv)


# Info on blank

names(compil.wInfo)

names(compil.wInfo)[1]

compil.wInfo[[1]] %>% filter(SeqType %in% c("sample", "mix", "blank")) %>% 
  group_by(SeqType, groups) %>% 
  summarise(Sum = mean(N)) %>% 
  spread(key = SeqType, value = Sum)
  




# Penser à enlever les singletons



# Vérifications particulières ---------------------------------------------


# Chrosomus et Luxilus 12S RDP vs IDT

colnames(taxo.IDT.12S[[2]])


res <- as.data.frame(cbind(taxo.IDT.12S[[2]][,c("Family","Genus","Species")], taxo.RDP.12S[,c("Family","Genus","Species")], colSums(seqtab.12s)))

res[,7] <- as.numeric(as.character(res[,7]))
str(res)
View(res)

# Etude des cyprinideae
res <- data.frame(taxo.IDT.12S[[2]], colSums(seqtab.12s))

View(res)                  



# Temoins negatifs --------------------------------------------------------

NAMES <- names(compil.wInfo)[str_detect(names(compil.wInfo), pattern = "12S.Class|12S.Family|12S.Genus|12S.Species")]

for(x in NAMES){
  cat("\n",x,"\n", sep="")
  compil.wInfo[[x]] %>% filter(SeqType %in% c("blank")) %>% 
    group_by(SampleID, groups) %>% 
    summarise(Sum = round(mean(N),1)) %>% 
    spread(key = SampleID, value = Sum) %>% 
    #kable(caption = x) %>% kable_styling(bootstrap_options = c("condensed"), full_width = T) %>% 
    print()
  cat("\n")
  
}


NAMES <- names(compil.wInfo)[str_detect(names(compil.wInfo), pattern = "cytB.Class|cytB.Family|cytB.Genus|cytB.Species")]

for(x in NAMES){
  cat("\n",x,"\n", sep="")
  compil.wInfo[[x]] %>% filter(SeqType %in% c("blank")) %>% 
    group_by(SampleID, groups) %>% 
    summarise(Sum = round(mean(N),1)) %>% 
    spread(key = SampleID, value = Sum) %>% 
    #kable(caption = x) %>% kable_styling(bootstrap_options = c("condensed"), full_width = T) %>% 
    print()
  cat("\n")
  
}

# N samples with Salvelinus


names(compil.wInfo)


NAMES <- names(compil.wInfo)[str_detect(names(compil.wInfo), pattern = "Genus")]

for(x in NAMES){
  cat("\n\n",x,"\n", sep="")
  
  DATA <- compil.wInfo[[x]] %>% filter(groups %in% c("Salvelinus"), SeqType == "sample", N>=1) 

  cat("N samples with Salvelinus:", length(DATA$NomLac), sep ="\n")
  cat("N lakes with Salvelinus:", length(unique(DATA$NomLac)), sep ="\n")

  cat(sort(unique(DATA$NomLac)) , sep=", ")

}



NAMES <- names(compil.wInfo)[str_detect(names(compil.wInfo), pattern = "Genus")]

for(x in NAMES){
  cat("\n\n",x,"\n", sep="")
  
  DATA <- compil.wInfo[[x]] %>% filter(groups %in% c("Perca"), SeqType == "sample", N>=1) 
  
  cat("N samples with Perca:", length(DATA$NomLac), sep ="\n")
  cat("N lakes with Perca:", length(unique(DATA$NomLac)), sep ="\n")
  
  cat(sort(unique(DATA$NomLac)) , sep=", ")
  
}


NAMES <- names(compil.wInfo)[str_detect(names(compil.wInfo), pattern = "Family")]

for(x in NAMES){
  cat("\n\n",x,"\n", sep="")
  
  DATA <- compil.wInfo[[x]] %>% filter(groups %in% c("Cyprinidae"), SeqType == "sample", N>=1) 
  
  cat("N samples with Cyprinidae:", length(DATA$NomLac), sep ="\n")
  cat("N lakes with Cyprinidae:", length(unique(DATA$NomLac)), sep ="\n")
  
  cat(sort(unique(DATA$NomLac)) , sep=", ")
  
}


for(x in NAMES){
  cat("\n\n",x,"\n", sep="")
  
  DATA <- compil.wInfo[[x]] %>% filter(groups %in% c("Gasterosteidae"), SeqType == "sample", N>=1) 
  
  cat("N samples with Gasterosteidae:", length(DATA$NomLac), sep ="\n")
  cat("N lakes with Gasterosteidae:", length(unique(DATA$NomLac)), sep ="\n")
  
  cat(sort(unique(DATA$NomLac)) , sep=", ")
  
}

for(x in NAMES){
  cat("\n\n",x,"\n", sep="")
  
  DATA <- compil.wInfo[[x]] %>% filter(groups %in% c("Catostomidae"), SeqType == "sample", N>=1) 
  
  cat("N samples with Catostomidae:", length(DATA$NomLac), sep ="\n")
  cat("N lakes with Catostomidae:", length(unique(DATA$NomLac)), sep ="\n")
  
  cat(sort(unique(DATA$NomLac)) , sep=", ")
  
}

unique(compil.wInfo[[3]]$groups)       
