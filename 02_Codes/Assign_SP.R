
# Info --------------------------------------------------------------------

# Assign species to different eDNA dataset
# 
# Audrey Bourret
# 2018-11-09 - Update January 11 2019
#

# Library -----------------------------------------------------------------

library(dada2); packageVersion("dada2")
library(DECIPHER); packageVersion("DECIPHER")
library(tidyverse)

# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}


# Data --------------------------------------------------------------------

# Path

#result.path <- get.value("result.path")
#ref.path    <- "./00_Data/03_RefSeq"
#log.path    <- "./03_Log" 

# Seqtab

list.files(get.value("result.data.path"))

#load(get.value("ASVtable.data"))
#load(get.value("OTUtable.data"))
load(get.value("CORRECTEDtable.data"))
load(get.value("CORRECTEDtable.data") %>% str_replace("table.data", "table.2x.data"))

  ls()
# Ref

SEQ.REF <- list.files(get.value("ref.path"), pattern = ".fasta",full.names = T) %>% str_subset("unique")
SEQ.REF


PARAM <-  expand.grid(LOCUS = c("12s", "cytB.R1", "cytB.R2"), 
                      OTU = c("ASV", "OTU"),
                      ASSIGN = c("RDP","IDT"),
                      SP = c("All"))

PARAM$TAB <- paste0(PARAM$OTU, "tab.", PARAM$LOCUS)
PARAM

# REF.12S.SP   <- "All_12S-eco_SP_dup.fasta"  
# REF.12S.taxo <- "All_12S-eco_taxo_dup.fasta"
# 
# REF.CYTB.SP   <- "All_CYTB-Kot_SP_dup.fasta"  
# REF.CYTB.taxo <- "All_CYTB-Kot_taxo_dup.fasta"
# 
# REF.ALL.taxo <- "All_ALL_taxo.fasta"

make.root <- function(FILE){
                      REF <- readDNAStringSet(FILE)
                      names(REF) <- stringr::str_replace(names(REF), "Animalia", "Root")
                      return(REF)
}


# Function to transform df to dada2 table (contrary of SEQtable.df) 
SEQtable.tr <- function(tab){
  
  tab <- tab %>% select(-ID)
  new.tab <- t(tab)
  
  
  return(new.tab)
}



#REF.12S.wROOT  <- make.root(file.path(get.value("ref.path"), REF.12S.taxo))
#REF.CYTB.wROOT <- make.root(file.path(get.value("ref.path"), REF.CYTB.taxo))

#REF.ALL.wROOT  <- make.root(file.path(ref.path, REF.ALL.taxo))


# RDP ---------------------------------------------------------------------


# Function to run assigntaxonomy then assignspecies  
RDP <- function(seqtab, REF.TAXO, REF.SP = NULL){
                
                #seqtab2 <- SEQtable.tr(seqtab)
                
                taxo <- assignTaxonomy(seqtab, REF.TAXO, 
                        taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species_80"), 
                        minBoot=80, tryRC = TRUE)
  
                #taxowSP <-  addSpecies(taxo, REF.SP, allowMultiple = T, tryRC = TRUE)
  
                #for(x in 1:nrow(taxowSP)){
                  
                #  taxowSP[x,"Species"] <- ifelse(!is.na(taxowSP[x, "Genus"])  & !is.na(taxowSP[x, "Species"]), paste(taxowSP[x, "Genus"], taxowSP[x, "Species"]), NA)

                #} 
                
                return(taxo)
}


taxoTAB <- list()

# DO the assignation

for(x in 1:nrow(PARAM)){
  
  if(PARAM[x,"ASSIGN"] == "RDP"){

      print(PARAM[x,"TAB"])
    
      REF.taxo <- SEQ.REF %>% str_subset("unique_wTaxo_dup.fasta") %>% str_subset(toupper(PARAM[x,"LOCUS"])) %>% str_subset(as.character(PARAM[x,"SP"]))
      REF.sp   <- SEQ.REF %>% str_subset("unique_dup.fasta") %>% str_subset(toupper(PARAM[x,"LOCUS"])) %>% str_subset(as.character(PARAM[x,"SP"]))
      DB.sp    <- REF.sp %>% str_replace(get.value("ref.path"), file.path(get.value("ref.path"),"Blast")) %>% str_replace(".fasta", "_DB.fasta")
      
      SEQTAB <- SEQtable.tr(get(PARAM[x,"TAB"]))
  
      # the function get() is use to call an object by its name
      taxoRPD <- RDP(SEQTAB, REF.taxo, REF.sp) 

      # To create a ref DB

      make.blast.db(FILE = REF.sp, 
                    DB = DB.sp)
      
      # To BLAST at 100% identity
       cat(file = file.path(get.value("ref.path"),"Blast", "input.fasta"), append = F)
       for(y in dimnames(SEQTAB)[[2]]){
               cat(paste0(">",y), y, sep = "\n", 
                   file = file.path(get.value("ref.path"),"Blast", "input.fasta"), 
                   append = T)
               }
      
      BLAST.RES <- BLAST.IDENTIC(FILE.TO.BLAST = file.path(get.value("ref.path"),"Blast", "input.fasta"), 
                                 DB =  DB.sp)
      
      taxoRPD2  <- cbind(taxoRPD, BLAST.RES$SP)
      dimnames(taxoRPD2)[[2]][8] <- "Species"
      
      taxoTAB[[paste(PARAM[x,"TAB"], "RDP", as.character(PARAM[x,"SP"]), sep=".")]]<- taxoRPD2
    
  }
}

names(taxoTAB)

View(taxoTAB[[1]])

# IDTAXA ------------------------------------------------------------------


TS.ls <-  list()

# Create training set

for(x in 1:nrow(PARAM)){
  
  if(PARAM[x,"ASSIGN"] == "IDT" & PARAM[x,"OTU"] == "ASV"){ # same training set ASV and OTU
 
   REF.taxo <- SEQ.REF %>% str_subset("unique_wTaxo_dup.fasta") %>% str_subset(toupper(PARAM[x,"LOCUS"]))  %>% str_subset(as.character(PARAM[x,"SP"]))
   print(REF.taxo)  
  
   REF.wROOT <- make.root(REF.taxo)
    
   TS <- LearnTaxa(REF.wROOT, names(REF.wROOT))
   
   TS.ls[[paste(as.character(PARAM[x,"SP"]), as.character(PARAM[x,"LOCUS"]), sep=".")]] <- TS
   
   print(TS$problemSequences)
   print(TS$problemGroups)   
   
   #TS.ls[[paste(as.character(PARAM[x,"SP"]), as.character(PARAM[x,"LOCUS"]), "wRANK", sep=".")]] <-  add.rank.TS(TS)

  }
  
}


names(TS.ls)

# Save training set
TS.ls

#save(file = get.value("IDT.TS.data"), 
#     list = c("TS.ls")
#)

load(get.value("IDT.TS.data"))


plot(TS.ls[["All.12s"]])
plot(TS.ls[["All.cytB.R1"]])
plot(TS.ls[["All.cytB.R2"]])

#plot(TS.ls[["QC.12s"]])
#plot(TS.ls[["QC.cytB.R1"]])
#plot(TS.ls[["QC.cytB.R2"]])


# Idtaxa on all data set

IDT <- function(SEQTAB, TS, ranks = c("Root", "Phylum", "Class", "Order", "Family", "Genus", "Species"),  R = FALSE){
  
# R for reverseComplement - I don't think its usefull anymore
  
  #ranks <- c("Root", "Phylum", "Class", "Order", "Family", "Genus", "Species") # ranks of interest
  
  if(R == FALSE){
    ids <- IdTaxa(DNAStringSet(getSequences(SEQTAB)), TS, strand="top", processors=1, verbose=TRUE) # use all processors
  } else {
    ids <- IdTaxa(reverseComplement(DNAStringSet(getSequences(SEQTAB))), TS, strand="top", processors=1, verbose=TRUE) # use all processors
  }
  
  for(x in 1:length(ids)){
    N <- length(ids[[x]]$taxon)
    ids[[x]]$rank <- ranks[1:N]
  }


  print(plot(ids, TS))
  
  return(ids)
  
}


IDT.dada2 <- function(ids, ranks = c("Root", "Phylum", "Class", "Order", "Family", "Genus", "Species")){
  taxid <- t(sapply(ids, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
  }))
  
  colnames(taxid) <- ranks; rownames(taxid) <- getSequences(SEQTAB)
  
  return(taxid)
  
}


# Do it now!

taxoTAB.IDT <- list()

for(x in 1:nrow(PARAM)){
  
  if(PARAM[x,"ASSIGN"] == "IDT"){
    
    print(paste(PARAM[x,"TAB"], (PARAM[x,"SP"])))
    
    #Get data
    SEQTAB <- SEQtable.tr(get(PARAM[x,"TAB"]))
    
    TS <- TS.ls[[paste(PARAM[x,"SP"], PARAM[x,"LOCUS"], sep=".")]]

    # IDT
    taxoIDT       <- IDT(SEQTAB, TS)
    taxoIDT.dada2 <- IDT.dada2(taxoIDT)             
    
    # Save results
    new.name <- paste(PARAM[x,"TAB"], "IDT", as.character(PARAM[x,"SP"]), sep=".")
    
    taxoTAB[[new.name]]     <- taxoIDT.dada2
    taxoTAB.IDT[[new.name]] <- taxoIDT
    
  }
}


names(taxoTAB)
names(taxoTAB.IDT)

View(taxoTAB[["ASVtab.12s.IDT.QC"]])


# Compare SEQ assignation -------------------------------------------------

names(taxoTAB)



SEQ.PARAM <- expand.grid(ASSIGN = c("RDP", "IDT"),
                         SP = c("All"))

COMP.SEQ <- list()

for(x in c("12s", "cytB.R1", "cytB.R2")){
  
  for(y in c("ASV", "OTU")){
    
     for(z in c("Family", "Genus", "Species")){
       
       TAB <- paste(paste0(y,"tab"),  x, SEQ.PARAM[,"ASSIGN"], SEQ.PARAM[,"SP"], sep=".") 
       
       DATA <- as.data.frame(cbind(taxoTAB[[TAB[1]]][,z], 
                     taxoTAB[[TAB[2]]][,z],
                     taxoTAB[[TAB[3]]][,z],
                     taxoTAB[[TAB[4]]][,z]))
         
       names(DATA) <- TAB
       
       COMP.SEQ[[paste(x,y,z,sep=".")]] <- DATA
       
     }
  }
}

names(COMP.SEQ)

View(COMP.SEQ[[1]])



# Compil datasets ---------------------------------------------------------

TAXO.final <- list()

for(x in c("12s", "cytB.R1", "cytB.R2")){
  
  for(y in c("ASV", "OTU")){
    
    CLASS <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
    
    name <- paste(paste0(y,"tab"), x, "IDT", "All", sep = ".")
    
    print(name)
    
    TAXO.idt <- as.data.frame(taxoTAB[[name]])
    TAXO.rdp <- as.data.frame(taxoTAB[[name %>% str_replace("IDT", "RDP")]])
    
    TAXO <- TAXO.idt %>% cbind(TAXO.rdp %>% select(CLASS))
    names(TAXO) <- c("Root", CLASS, paste0(CLASS, "100"))
    
    # To change 100% identity assignation whenever possible
    TAXO <- TAXO %>%  mutate(Phylum  = ifelse(!is.na(Species), as.character(Phylum),  as.character(Phylum100)),
                             Class   = ifelse(!is.na(Species), as.character(Class),   as.character(Class100)),
                             Order   = ifelse(!is.na(Species), as.character(Order),   as.character(Order100)),
                             Family  = ifelse(!is.na(Species), as.character(Family),  as.character(Family100)),
                             Genus   = ifelse(!is.na(Species), as.character(Genus),   as.character(Genus100)),
                             Species = ifelse(!is.na(Species), as.character(Species), as.character(Species100))
                             )
    
    
    TAXO <- TAXO %>% mutate(TaxoFinal = paste(Root, Phylum, Class, Order, Family, Genus, Species, sep = ";"),
                            TaxoFinal = str_replace_all(TaxoFinal, ";NA", "")) %>% 
                     transmute(SEQ = row.names(TAXO.idt), Assign = TaxoFinal)

    TAXO.final[[name %>% str_remove(".IDT.All")]] <- TAXO

  }
  
}


RES <- data.frame(Level = NULL, N = NULL, DATA = NULL)

for(x in 1:length(TAXO.final)){
  NEW <- str_count(TAXO.final[[x]]$Assign, ";") %>% 
    table() %>% 
    as.data.frame()
  names(NEW) <- c("Level", "N") 
    
  NEW$LOCUS <- names(TAXO.final[x]) %>% str_remove("tab") %>% str_remove("ASV.") %>% str_remove("OTU.")
  NEW$OTU   <- names(TAXO.final[x]) %>% str_remove("tab") %>% str_remove(".12s") %>% str_remove(".cytB.R[:digit:]")
  
  
  RES <- rbind(RES, NEW)  
  
}

RES %>% ggplot(aes(x = as.numeric(as.character(Level)), y = N, fill = LOCUS)) + 
  geom_bar(stat="identity",  position = "dodge") + 
  facet_grid(.~ OTU) + 
  #ylim(0, 200) +
  theme_bw()


names(TAXO.final)


# Add taxo to all seqTAB

for(x in unique(PARAM$TAB)){
   
  print(x)

  TAB  <- get(x)
  TAB2 <- get(paste(x, "cor", sep= "."))
  TAB3 <- get(paste(x, "cor.2x", sep= "."))
  TAXO <- TAXO.final[[x]]
  
  nrow(TAB) == nrow(TAXO)
  
  TAB <- TAB %>% mutate(SEQ = row.names(.))
  
  NEW1 <- TAB %>% full_join(TAXO) %>% 
          select(-SEQ)
  
  NEW2 <- TAB2 %>% full_join(TAB %>% select(ID, SEQ)) %>% 
           full_join(TAXO) %>% 
           select(-SEQ) 
  
  NEW3 <- TAB3 %>% full_join(TAB %>% select(ID, SEQ)) %>% 
           full_join(TAXO) %>% 
           select(-SEQ) 
  
  assign(paste(x, "wTAXO", sep = "."), NEW1)
  assign(paste(x, "cor", "wTAXO", sep = "."), NEW2)
  assign(paste(x, "cor.2x", "wTAXO", sep = "."), NEW3)
}


# Enregistrer les SEQTAB corrigées
save(file = get.value("ALLtable.data"), 
     list = ls() %>% str_subset("wTAXO")) # Pour enlever "SEQtable.df" 


TABLEwTAXO <- ls() %>% str_subset("wTAXO")

rm("TABLEwTAXO")

for(x in ls() %>% str_subset("wTAXO")){
  DATA <- get(x)
  
  SAMPLE <- names(DATA %>% select(-c(ID, Assign)))
  
  NEW <- DATA %>% gather(SAMPLE, key= "Sample", value = Read) %>% 
           group_by(Assign, Sample) %>% 
           summarise(N = sum(Read)) %>% 
           mutate(Level = str_count(Assign, ";"))

  assign(x %>%str_replace("wTAXO", "bySP"), NEW)
  
}



# Save --------------------------------------------------------------------


# Enregistrer les SEQTAB corrigées
save(file = get.value("ALLtable.data"), 
     list = c(ls() %>% str_subset("wTAXO"),
              ls() %>% str_subset("bySP"))
   
     )

