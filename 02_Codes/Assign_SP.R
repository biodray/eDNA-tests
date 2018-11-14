
# Info --------------------------------------------------------------------

# Assign species to different eDNA dataset
# 
# Audrey Bourret
# 2018-11-09
#


# Library -----------------------------------------------------------------

library(dada2); packageVersion("dada2")
library(DECIPHER); packageVersion("DECIPHER")
library(tidyverse)

# Internal functions
source(file.path("./03_Functions",  list.files("./03_Functions")))

# Data --------------------------------------------------------------------

# Path

result.path <- get.value("result.path")
ref.path    <- "./00_Data/03_RefSeq"
#log.path    <- "./03_Log" 

# Seqtab

list.files(get.value("result.path"))

load(file.path(get.value("result.path"), "Seqtab.data"))

# Ref

list.files(get.value("ref.path"))

REF.12S.SP   <- "All_12S-eco_SP_dup.fasta"  
REF.12S.taxo <- "All_12S-eco_taxo_dup.fasta"

REF.CYTB.SP   <- "All_CYTB-Kot_SP_dup.fasta"  
REF.CYTB.taxo <- "All_CYTB-Kot_taxo_dup.fasta"

REF.ALL.taxo <- "All_ALL_taxo.fasta"

make.root <- function(FILE){
                      REF <- readDNAStringSet(FILE)
                      names(REF) <- stringr::str_replace(names(REF), "Animalia;Chordata", "Root")
                      return(REF)
}


REF.12S.wROOT  <- make.root(file.path(get.value("ref.path"), REF.12S.taxo))
REF.CYTB.wROOT <- make.root(file.path(get.value("ref.path"), REF.CYTB.taxo))

#REF.ALL.wROOT  <- make.root(file.path(ref.path, REF.ALL.taxo))


# RDP ---------------------------------------------------------------------


# Function to run assigntaxonomy then assignspecies  
RDP <- function(seqtab, REF.TAXO, REF.SP){
                taxo <- assignTaxonomy(seqtab, REF.TAXO, 
                        taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species_80"), 
                        minBoot=80, tryRC = TRUE)
  
                taxowSP <-  addSpecies(taxo, REF.SP, tryRC = TRUE)
  
                for(x in 1:nrow(taxowSP)){
                  
                  taxowSP[x,"Species"] <- ifelse(!is.na(taxowSP[x, "Genus"])  & !is.na(taxowSP[x, "Species"]), paste(taxowSP[x, "Genus"], taxowSP[x, "Species"]), NA)

                } 
                
                return(taxowSP)
}

taxo.RDP.12S.F <- RDP(seqtab.12s.F, file.path(get.value("ref.path"), REF.12S.taxo), file.path(get.value("ref.path"), REF.12S.SP))
taxo.RDP.12S.R <- RDP(seqtab.12s.R, file.path(get.value("ref.path"), REF.12S.taxo), file.path(get.value("ref.path"), REF.12S.SP))
taxo.RDP.12S   <- RDP(seqtab.12s, file.path(get.value("ref.path"), REF.12S.taxo), file.path(get.value("ref.path"), REF.12S.SP))

taxo.RDP.12S.IBIS   <- RDP(seqtab.12s.IBIS, file.path(get.value("ref.path"), REF.12S.taxo), file.path(get.value("ref.path"), REF.12S.SP))

taxo.RDP.cytB.F <- RDP(seqtab.cytB.F, file.path(get.value("ref.path"), REF.CYTB.taxo), file.path(get.value("ref.path"), REF.CYTB.SP))
taxo.RDP.cytB.R <- RDP(seqtab.cytB.R, file.path(get.value("ref.path"), REF.CYTB.taxo), file.path(get.value("ref.path"), REF.CYTB.SP))
#taxo.RDP.cytB   <- RDP(seqtab.cytB, file.path(ref.path, REF.CYTB.taxo), file.path(ref.path, REF.CYTB.SP))


nrow(taxo.RDP.12S.F)


View(taxo.RDP.cytB.F)


# IDTAXA ------------------------------------------------------------------


# Create training set


TS.12S <- LearnTaxa(REF.12S.wROOT, names(REF.12S.wROOT))
TS.12S

TS.12S$problemSequences
TS.12S$problemGroups

TS.CYTB <- LearnTaxa(REF.CYTB.wROOT, names(REF.CYTB.wROOT))

TS.CYTB$problemGroups
TS.CYTB$problemSequences


TS.CYTB <- LearnTaxa(REF.CYTB.wROOT[-87], names(REF.CYTB.wROOT[-87]))
TS.CYTB

plot(TS.12S)
plot(TS.CYTB)

# Manually add ranks
add.rank.TS <- function(TS){

  ranks <- c("Root", "Class", "Order", "Family", "Genus", "Species") # ranks of interest
  
  add.rank <- vector()

  for(x in 1:length(TS[["taxonomy"]])){
    
    N <- str_count( TS[["taxonomy"]][x], pattern = ";")
  
    add.rank <- c(add.rank, ranks[N])

  }

  TS[["rank"]] <- add.rank

  return(TS)
}


TS.12SwR  <- add.rank.TS(TS.12S)
TS.CYTBwR <- add.rank.TS(TS.CYTB)

# Idtaxa on all data set

IDT <- function(SEQTAB, TS, R = FALSE){
  
  ranks <- c("Root", "Class", "Order", "Family", "Genus", "Species") # ranks of interest
  
  if(R == FALSE){
    ids <- IdTaxa(DNAStringSet(getSequences(SEQTAB)), TS, strand="top", processors=1, verbose=TRUE) # use all processors
  } else {
    ids <- IdTaxa(reverseComplement(DNAStringSet(getSequences(SEQTAB))), TS, strand="top", processors=1, verbose=TRUE) # use all processors
  }
  
  for(x in 1:length(ids)){
    N <- length(ids[[x]]$taxon)
    ids[[x]]$rank <- ranks[1:N]
  }

  taxid <- t(sapply(ids, function(x) {
                         m <- match(ranks, x$rank)
                         taxa <- x$taxon[m]
                         taxa[startsWith(taxa, "unclassified_")] <- NA
                         taxa
                         }))

  colnames(taxid) <- ranks; rownames(taxid) <- getSequences(SEQTAB)

  plot(ids, TS)
  
  return(list(IDT = ids, IDT.dada2 = taxid))
  
  }


taxo.IDT.12S        <- IDT(seqtab.12s, TS.12S, R = FALSE)
taxo.IDT.12S.IBIS   <- IDT(seqtab.12s.IBIS, TS.12S, R = FALSE)

taxo.IDT.12S.F <- IDT(seqtab.12s.F, TS.12S, R = FALSE)
taxo.IDT.12S.R <- IDT(seqtab.12s.R, TS.12S, R = TRUE)


taxo.IDT.cytB   <- IDT(seqtab.cytB, TS.CYTB, R = FALSE)
taxo.IDT.cytB.F <- IDT(seqtab.cytB.F, TS.CYTB, R = FALSE)
taxo.IDT.cytB.R <- IDT(seqtab.cytB.R, TS.CYTB, R = TRUE)



# Compil datasets ---------------------------------------------------------



compil.taxo <- function(TAXO, SEQTAB, CLASS = c("Class", "Order", "Family", "Genus", "Species")){
  compil <- list()
  
  row.names(TAXO) <- NULL
  
  for(x in CLASS){
    TAXO[which(is.na(TAXO[,x])),x] <- "Unknown" 
    colnames(SEQTAB) <- TAXO[,x]     
    
    df <- data.frame(Sample = rownames(SEQTAB))
    
    #Compiler les résultats
    for(n in sort(unique(TAXO[,x]))){
      DATA <- SEQTAB[,which(TAXO[,x] == n)]
      if(is.vector(DATA)){
        df[,n] <-  DATA     
      } else {
        df[,n] <-  rowSums(DATA)     
      }
    }
  
    compil[[x]] <- df
      
  }
  
  return(compil)

}


compil.RDP.12S <- compil.taxo(taxo.RDP.12S, seqtab.12s)
compil.IDT.12S <- compil.taxo(taxo.IDT.12S[[2]], seqtab.12s)

compil.RDP.12S.IBIS <- compil.taxo(taxo.RDP.12S.IBIS, seqtab.12s.IBIS)
compil.IDT.12S.IBIS <- compil.taxo(taxo.IDT.12S.IBIS[[2]], seqtab.12s.IBIS)

compil.RDP.12S.F <- compil.taxo(taxo.RDP.12S.F, seqtab.12s.F)
compil.IDT.12S.F <- compil.taxo(taxo.IDT.12S.F[[2]], seqtab.12s.F)

compil.RDP.12S.R <- compil.taxo(taxo.RDP.12S.R, seqtab.12s.R)
compil.IDT.12S.R <- compil.taxo(taxo.IDT.12S.R[[2]], seqtab.12s.R)


compil.RDP.cytB.F <- compil.taxo(taxo.RDP.cytB.F, seqtab.cytB.F)
compil.IDT.cytB.F <- compil.taxo(taxo.IDT.cytB.F[[2]], seqtab.cytB.F)

compil.RDP.cytB.R <- compil.taxo(taxo.RDP.cytB.R, seqtab.cytB.R)
compil.IDT.cytB.R <- compil.taxo(taxo.IDT.cytB.R[[2]], seqtab.cytB.R)

compil.IDT.cytB <- compil.taxo(taxo.IDT.cytB[[2]], seqtab.cytB)


# Save results ------------------------------------------------------------

# Save compile

save(file = file.path(get.value("result.path"), "Compil.data"), 
     list = ls(pattern = "compil."))

save(file = file.path(get.value("result.path"), "Taxo.data"), 
     list = c(ls(pattern = "taxo.IDT"),ls(pattern = "taxo.RDP")))

save.image(file.path(get.value("log.path"),"Assign_SP.Rdata"))
