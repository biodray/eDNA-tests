
# Info --------------------------------------------------------------------

# Assign species to different eDNA dataset
# 
# Audrey Bourret
# 2018-11-09
#


seqtab.path <- "./01_Results/01_Seqtab"

list.files(seqtab.path)

load()


library(dada2); packageVersion("dada2")
library(DECIPHER); packageVersion("DECIPHER")

# Assign species

A <- assignTaxonomy(seqtab.12s.seta, file.path(ref.path, "All_12S-eco_taxo.fasta"), minBoot=80)
B <- assignTaxonomy(seqtab.12s.setb, file.path(ref.path, "All_12S-eco_taxo.fasta"), minBoot=80)


C <- assignSpecies(seqtab.12s.seta, file.path(ref.path, "All_12S-eco.fasta"))


View(A)
View(B)
View(C)

nrow(A[which(!is.na(A[,"Family"])),])
nrow(B[which(!is.na(A[,"Family"])),])

nrow(C[which(!is.na(C[,"Species"])),])

View(A[which(!is.na(C[,"Species"])),])

rownames(C[which(C[,"Species"]=="salar"),])

#save.image(file="2018-11-06.Rdata")

load("2018-11-06.Rdata")

# Salmo salar
seqtab.12s.seta[,"ACTATGCCTAGCCGTAAACTTTGATGAAAACATACAACTGACATCCGCCAGGGGACTATAAGCGCCAGCTTAAAACCCAAAGGACTTGGCGGTGCCTCAGACCCAC"]

# Luxilus cornutus
seqtab.12s.seta[,"ACTATGCTCAGTCATAAACCCAGACGCTCAACTACAATCAGGCGTCCGCCCGGGTACTACGAGCATTAGCTTGAAACCCAAAGGACCTGACGGTGCCTCAGACCCCC"]

seqtab.12s.seta[,"ACTATGCCTAGCCGTAAACTTTGATAGAAAAATACAACTGATATCCGCCAGGGAACTACAAGCGCCAGCTTAAAACCCAAAGGACTTGGCGGTGCCTCAGACCCAC"]
# Salvelinus

ACTATGCCTAGCCGTAAACTTTGATAGAAAAATACAACTGATATCCGCCAGGGAACTACAAGCGCCAGCTTAAAACCCAAAGGACTTGGCGGTGCCTCAGACCCAC

BwDUP <- assignTaxonomy(seqtab.12s.setb, file.path(ref.path, "All_12S-eco_taxo_dup.fasta"), minBoot=80)
CwDUP  <- assignSpecies(seqtab.12s.seta, file.path(ref.path, "All_12S-eco_SP_dup.fasta"))

D <- addSpecies(BwDUP, file.path(ref.path, "All_12S-eco_SP_dup.fasta"))

View(D)

unname(C)
unname(CwDUP)

cbind(unname(BwDUP)[,7], assigment)
cbind(unname(D[,7:8]), assigment)

View(BwDUP)

nrow(BwDUP[which(!is.na(BwDUP[,"Species"])),])
nrow(B[which(!is.na(B[,"Species"])),])

# Decipher

getSequences(BwDUP[,])

library(DECIPHER); packageVersion("DECIPHER")

## [1] '2.6.0'

dna <- DNAStringSet(getSequences(seqtab.12s.seta)) # Create a DNAStringSet from the ASVs
dna

# Create a training set

ref <- readDNAStringSet(file.path(ref.path, "All_12S-eco_taxo.fasta"))
names(ref) <- str_replace(names(ref), "Animalia;Chordata", "Root")

trainingSet <- LearnTaxa(ref, names(ref))
trainingSet

trainingSet$problemSequences

as.character(ref[[128]])


plot(trainingSet)

ids <- IdTaxa(dna, trainingSet, strand="top", processors=1, verbose=TRUE) # use all processors

plot(ids, trainingSet)

assigment <- sapply(ids, function(x) {paste(x$taxon, collapse = ";")})

assigment[str_detect(assigment, "Salvelinus")]





ranks <- c("phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.12s.seta)


taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)



