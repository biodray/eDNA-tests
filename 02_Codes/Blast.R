


library(Biostrings)


# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}



# Data --------------------------------------------------------------------


load(get.value("CORRECTEDtable.data"))
load(get.value("ALLtable.data"))

ASVtab.12s.wTAXO %>% filter(ID == "ASV_17")

ASVtab.12s[ ,] %>% View()

DNA <- DNAStringSet(c(row.names(ASVtab.12s), row.names(OTUtab.12s)))

names(DNA) <- c(ASVtab.12s$ID, OTUtab.12s$ID)

DNA

writeXStringSet(DNA, "ASVOTU.db.fasta")


writeXStringSet(readDNAStringSet(file.path(get.value("ref.path"),"QC_12S-eco_unique_dup.fasta")), "SP.db.fasta")


# Blast -------------------------------------------------------------------

# Create DB

get.value("makeblastdb")


cmd <- paste("-in", "ASVOTU.db.fasta",
             "-dbtype", "nucl",
             "-parse_seqids", 
             sep = " ") # forward adapter

system2(get.value("makeblastdb"), cmd, stdout=T, stderr=T) 



cmd <- paste("-in", "SP.db.fasta",
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
