# Blast functions

# More to come

# To create a BLAST DB
make.blast.db <- function(FILE, DB, EXE = get.value("makeblastdb")){
  
#@FILE : fasta files path
#@DB :  name of the new DB
#@EXE : executive path to makeblastdb  
  
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


# To BLAST to 100% identity

BLAST.IDENTIC <- function(FILE.TO.BLAST, DB, EXE= get.value("blastn")){
  
#@FILES.TO.BLAST : fasta files path
#@DB : reference BLAST DB
#@EXE : executive path to blastn
  
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
  
  RES <- data.frame(ID = readDNAStringSet(FILE.TO.BLAST)@ranges@NAMES) %>% 
    left_join(BLAST.RES.100, by = c("ID" = "query"))
  
  return(RES)
  
}



# To BLAST to 100% identity

BLAST.seuil <- function(FILE.TO.BLAST, DB, seuil = 99, EXE= get.value("blastn")){
  
  #@FILES.TO.BLAST : fasta files path
  #@DB : reference BLAST DB
  #@EXE : executive path to blastn
  
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
  
  BLAST.RES.100 <- BLAST.RES %>% filter(Identity >= seuil,
                                        AlignmentLength == RealLength) %>% 
    select(query, SP) %>% 
    group_by(query) %>% 
    summarise(SP = str_flatten(unique(SP), collapse = "/"))
  
  # Prepare results
  
  RES <- data.frame(ID = readDNAStringSet(FILE.TO.BLAST)@ranges@NAMES) %>% 
    left_join(BLAST.RES.100, by = c("ID" = "query"))
  
  return(RES)
  
}


