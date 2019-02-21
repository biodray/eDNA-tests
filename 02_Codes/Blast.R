


library(Biostrings)


# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}



# Data --------------------------------------------------------------------


get.value("ref.path")


load(get.value("CORRECTEDtable.data"))
load(get.value("ALLtable.data"))

ASVtab.12s.wTAXO %>% filter(ID == "ASV_17")

ASVtab.12s[ ,] %>% View()

DNA <- DNAStringSet(c(row.names(ASVtab.12s), row.names(OTUtab.12s)))

names(DNA) <- c(ASVtab.12s$ID, OTUtab.12s$ID)

DNA

writeXStringSet(DNA, file.path(get.value("ref.path"),"Blast","ASVOTU.db.fasta"))


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



# Function ----------------------------------------------------------------

list.files(get.value("ref.path"))


FILE <- "./00_Data/05_RefSeq/QC_12S-eco_unique_dup.fasta"
NEW.FILE <- FILE %>% str_replace(get.value("ref.path"), file.path(get.value("ref.path"),"Blast")) %>% str_replace(".fasta", "_BD.fasta")  

# Function to create a BLAST DB

make.blast.db <- function(FILE, DB, EXE = get.value("makeblastdb")){
     
    # Move the fasta file
    DNA <- readDNAStringSet(FILE) 
    writeXStringSet(DNA,DB)

    # Run the commande
    
    cmd <- paste("-in", DB,
                 "-dbtype", "nucl",
                 "-parse_seqids", 
                 sep = " ") 
    
    system2(EXE, cmd, stdout=T, stderr=T) 
  
}

FILE <- "./00_Data/05_RefSeq/All_12S-eco_unique_dup.fasta"
DB <- FILE %>% str_replace(get.value("ref.path"), file.path(get.value("ref.path"),"Blast")) %>% str_replace(".fasta", "_DB.fasta")  

make.blast.db(FILE, DB) # it works


BLAST.IDENTIC <- function(FILES.TO.BLAST, DB, EXE= get.value("blastn")){

   FILE.TO.BLAST <- file.path(get.value("ref.path"),"LabSeq", "LabSeq.fasta")
   
   OUT <- DB %>% str_replace(".fasta", ".result.out")
   
   # Obtain the length of each reference, and SP name
   DNA <- readDNAStringSet(DB)
   ALIGNLENGTH.MIN <- data.frame(ID = keep_ID(DNA@ranges@NAMES),
                                 SP =  keep_sp(DNA@ranges@NAMES),
                                 RealLength = DNA@ranges@width )
   
   # BLAST from command lline
   cmd <- paste("-db", DB,
                "-query", FILE.TO.BLAST,
                "-outfmt", "7",
                "-out", OUT, 
                sep = " ") # forward adapter
   
   system2(EXE, cmd, stdout=T, stderr=T) 
   
   # Load results
   BLAST.RES <- read.table(OUT)
   
   names(BLAST.RES) <- c("query", "db", "Identity", "AlignmentLength", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")
   
   BLAST.RES <- BLAST.RES %>% left_join(ALIGNLENGTH.MIN, by = c("db" = "ID"))
   
   BLAST.RES.100 <- BLAST.RES %>% filter(Identity == 100,
                                         AlignmentLength == RealLength) %>% 
                                  select(query, SP) %>% 
                                  group_by(query) %>% 
                                  summarise(SP = str_flatten(unique(SP), collapse = "/"))
   
   # Prepare results
   
   RES <- data.frame(ID = BLAST.RES %>% pull(query) %>% unique()) %>% 
          left_join(BLAST.RES.100, by = c("ID" = "query"))
   
   return(RES)

}



BLAST.IDENTIC(FILES.TO.BLAST, DB)

