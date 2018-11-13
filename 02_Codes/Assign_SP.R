
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

# Data --------------------------------------------------------------------

# Path
seqtab.path <- "./01_Results/01_Seqtab"
ref.path <- "./00_Data/03_RefSeq"

# Seqtab

list.files(seqtab.path)

load(file.path(seqtab.path, "Dada2.Seqtab.data"))
load(file.path(seqtab.path, "Dada2.Seqtab.IBIS.data"))

# Ref

list.files(ref.path)

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


REF.12S.wROOT  <- make.root(file.path(ref.path, REF.12S.taxo))
REF.CYTB.wROOT <- make.root(file.path(ref.path, REF.CYTB.taxo))

REF.ALL.wROOT <- make.root(file.path(ref.path, REF.ALL.taxo))



# RDP ---------------------------------------------------------------------

# Assign species

# Check within F vs R vs merged
# Check across datasets with the same methods (Biodiv)
# Check across methods for the same data set


seqtab.12s[which(seqtab.12s[,"ACTATGCCTAGCCGTAAACTTTGATGAAAACATACAACTGACATCCGCCAGGGGACTATAAGCGCCAGCTTAAAACCCAAAGGACTTGGCGGTGCCTCAGACCCAC"]>0), "ACTATGCCTAGCCGTAAACTTTGATGAAAACATACAACTGACATCCGCCAGGGGACTATAAGCGCCAGCTTAAAACCCAAAGGACTTGGCGGTGCCTCAGACCCAC"] 






B[c("ACTATGCTCAGCCATAAACCTAGATGTCCAACTACAATTAGACGTCCACCCGGGTACTACGAGCATTAGCTTGAAACCCAAAGGACCTGACGGTGCCTTAGACCCCC",
    "ACTATGCTTAGCCTTAAACCCAGATGTATTCTTACACACACATCCGCCCGGGTACTACGAGCATAGCTTAAGACCCAAAGGACTTGGCGGTGTCTCAGACCCAC",
    "ACTATGCCTAGCCATAAACATTGATAGAACAATACAACCTCTATCCGCCAGGGGACTACAAGCATCAGCTTAAAGCCCAAAGGACTTGGCGGTGCTTTAGATCCAC",
    "ACTATGCCTAGCCATAAACATTGGTAGCACATTACACCCACTACCCGCCTGGGAACTACGAGCATCAGCTTGAAACCCAAAGGACTTGGCGGTGCTTTAGATCCAC"),]

seqtab.12s["12s-p2-D10",]

RDP <- function(seqtab, REF.TAXO, REF.SP){
                taxo <- assignTaxonomy(seqtab, REF.TAXO, 
                        taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species_80"), 
                        minBoot=80, tryRC = TRUE)
  
                taxowSP <-  addSpecies(taxo, REF.SP, tryRC = TRUE)
  
                return(taxowSP)
}

taxo.RDP.12S.F <- RDP(seqtab.12s.F, file.path(ref.path, REF.12S.taxo), file.path(ref.path, REF.12S.SP))
taxo.RDP.12S.R <- RDP(seqtab.12s.R, file.path(ref.path, REF.12S.taxo), file.path(ref.path, REF.12S.SP))
taxo.RDP.12S   <- RDP(seqtab.12s, file.path(ref.path, REF.12S.taxo), file.path(ref.path, REF.12S.SP))


taxo.RDP.cytB.F <- RDP(seqtab.cytB.F, file.path(ref.path, REF.CYTB.taxo), file.path(ref.path, REF.CYTB.SP))
taxo.RDP.cytB.R <- RDP(seqtab.cytB.R, file.path(ref.path, REF.CYTB.taxo), file.path(ref.path, REF.CYTB.SP))
#taxo.RDP.cytB   <- RDP(seqtab.cytB, file.path(ref.path, REF.CYTB.taxo), file.path(ref.path, REF.CYTB.SP))



sort(seqtab.cytB.F[,"TTTGGTTCACTCCTAGGCCTATGTTTAGCCACCCAAATTCTTACCGGACTCTTCCTAGCCATACA"])


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
  ranks <- c("root", "class", "order", "family", "genus", "species") # ranks of interest
  add.rank <- vector()

for(x in 1:length(TS[["taxonomy"]])){
  N <- str_count( TS[["taxonomy"]][x], pattern = ";")
  
  add.rank <- c(add.rank, ranks[N])


}

    
  TS[["rank"]] <- add.rank

  return(TS)
}


TS.12SwR <- add.rank.TS(TS.12S)
TS.CYTBwR <- add.rank.TS(TS.CYTB)


names(REF.ALL.wROOT)
names(REF.12S.wROOT)

test <- DNAStringSet(getSequences(taxo.RDP.12S[which(!is.na(taxo.RDP.12S[,"Genus"])),]))
names(test) <- paste(taxo.RDP.12S[which(!is.na(taxo.RDP.12S[,"Genus"])),"Genus"], taxo.RDP.12S[which(!is.na(taxo.RDP.12S[,"Genus"])),"Species"])
test

table(taxo.RDP.12S[,"Genus"], useNA = "always")
    


str(ids)  
      
rapid_tree(test)

ids<- IdTaxa(DNAStringSet(getSequences(seqtab.12s)), TS.12S, strand="top", processors=1, verbose=TRUE) # use all processors
ids<- IdTaxa(DNAStringSet(getSequences(seqtab.12s.IBIS)), TS.12S, strand="top", processors=1, verbose=TRUE) # use all processors
 
plot(ids, TS.12S)


ids <- IdTaxa(reverseComplement(DNAStringSet(getSequences(seqtab.cytB.R))), TS.CYTB, strand="top", processors=1, verbose=TRUE) # use all processors
ids <- IdTaxa(DNAStringSet(getSequences(seqtab.cytB)), TS.CYTB, strand="top", processors=1, verbose=TRUE) # use all processors



ids <- IdTaxa(DNAStringSet(getSequences(seqtab.cytB.F)), TS.ALL, strand="top", processors=1, verbose=TRUE) # use all processors

plot(ids, TS.ALL)



str(plot(ids, TS.CYTB))


assigment <- sapply(ids, function(x) {paste(x$taxon, collapse = ";")})

assigment[str_detect(assigment, "Lampe")]






# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy

# add $rank

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

colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.12s)


View(taxid)



str(TS.12S)



plot(ids, TS.12S)[1]

plot(TS.12S)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)



