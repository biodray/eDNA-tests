# Info --------------------------------------------------------------------

# Code to BLAST various stuffs whenever necessary
# 
# Audrey Bourret
# 2019-02-21
#

# Library -----------------------------------------------------------------

library(Biostrings)
library(tidyverse)

# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}


# Data --------------------------------------------------------------------

# Create a ASV/OTU ref DB

get.value("ref.path")

load(get.value("CORRECTEDtable.data"))
load(get.value("ALLtable.data"))


DNA <- DNAStringSet(c(row.names(ASVtab.12s), row.names(OTUtab.12s)))

names(DNA) <- c(ASVtab.12s$ID, OTUtab.12s$ID)

DNA

writeXStringSet(DNA, file.path(get.value("ref.path"),"Blast","ASVOTU.db.fasta"))


# Blast -------------------------------------------------------------------

# Create DB

get.value("makeblastdb")


cmd <- paste("-in", "ASVOTU.db.fasta",
             "-dbtype", "nucl",
             "-parse_seqids", 
             sep = " ") # forward adapter

system2(get.value("makeblastdb"), cmd, stdout=T, stderr=T) 





# Blast it!


SEQ <- readDNAStringSet(file.path(get.value("ref.path"),"QC_12S-eco_unique_wTAXO.fasta"))
names(SEQ)

SEQ <- readDNAStringSet(file.path(get.value("ref.path"),"LabSeq", "LabSeq.fasta"))
SEQ 
                        
get.value("blastn")


FILES <- file.path(get.value("ref.path"),"LabSeq", "LabSeq.fasta")

cmd <- paste("-db", "ASVOTU.db.fasta",
             "-query", FILES,
             "-outfmt", "7",
             "-out", "results.out", 
             sep = " ") # forward adapter

cmd <- paste("-db", "SP.db.fasta",
             "-query", FILES,
             "-outfmt", "7",
             "-out", "results.out", 
             sep = " ") 
system2(get.value("blastn"), cmd, stdout=T, stderr=T) 



RES <- read.table("results.out")

names(RES) <- c("query acc.ver", "subject acc.ver", "Identity", "AlignmentLength", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")

View(RES)

View(RES %>% filter(Identity >= 99, 
                    AlignmentLength >= 100))

str(RES)



# Function ----------------------------------------------------------------

list.files(file.path(get.value("ref.path"),"blast"), pattern = "ASV", full.names = T)

FILE <- "./00_Data/05_RefSeq/All_12S-eco_unique_dup.fasta"
DB <- FILE %>% str_replace(get.value("ref.path"), file.path(get.value("ref.path"),"Blast")) %>% str_replace(".fasta", "_DB.fasta")  

make.blast.db(FILE, DB) # it works

FILE.TO.BLAST <- file.path(get.value("ref.path"),"LabSeq", "LabSeq.fasta")


BLAST.IDENTIC(FILES.TO.BLAST, DB)


Blast99 <- BLAST.seuil(FILE.TO.BLAST = "./00_Data/05_RefSeq/blast/ASVOTU.db.fasta" ,
            DB = "./00_Data/05_RefSeq/blast/All_12S-eco_unique_dup_DB.fasta",
            seuil = 99)

View(RES)


save(Blast99, file = get.value("Blast99.data"))

