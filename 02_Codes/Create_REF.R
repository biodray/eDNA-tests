# Description -------------------------------------------------------------

# Code pour creer la liste de reference
# Audrey Bourret
# 31 octobre 2018

# Packages ----------------------------------------------------------------

library(readxl)
library(tidyverse)

library(Biostrings)

library(msa)
library(ape)
library(seqinr)
#library(gtools)

#library(ShortRead)
#library(dada2); packageVersion("dada2")

# Internal functions
source(file.path("./03_Functions",  list.files("./03_Functions")))

# Path --------------------------------------------------------------------

data.path <-"./00_Data/00_FileInfos"
path.group   <- "./00_Data/03_RefSeq"

biodiv.path  <- "S:/Genpop/01-Projets de recherche/Banque reference biodiversite ESTL GSTL 2017-"
path.EXTERNE <- "Sequences/Sequences_EXTERNE"
path.LABO    <- "Sequences/Sequences_LABO_finales"

Sample.xl <- "DB_Echantillons.xlsx"

# Functions ---------------------------------------------------------------

source(file.path(biodiv.path, "R/fct_merge_sequences.R"))

# Ajouter la taxonomie a un fichier fasta
add_taxo <- function(fn, TAXO.REF = REF, fn.new = NULL, verbose = FALSE){

  #@fn       : nom du fichier fasta d'entré
  #@TAXO.REF : Fichier de reference taxonomique
  #@fn.new   : nom du fichier de sortie si différent 
  
   # Read files
  DNA <- readDNAStringSet(fn)
  
  if(str_detect(names(DNA)[1], ";")){
    stop("La taxonomie semble deja ajoutee, action non effectuee", call. = F)
  }
  
  # Extract names
  NAMES <- paste(sapply(str_split(names(DNA), pattern=" "), `[`, 2),
                 sapply(str_split(names(DNA), pattern=" "), `[`, 3))
  
  TAXO <- data.frame(SP = NAMES)
  TAXO <- TAXO %>% left_join(TAXO.REF, by = c("SP" = "Espece_initial"))
  
  FINAL <- paste(TAXO$Kingdom,  
                 TAXO$Phylum, 
                 TAXO$Class,
                 TAXO$Order,
                 TAXO$Family,
                 TAXO$Genus,
                 TAXO$Espece,
                 sep = ";")
  
  if(verbose == TRUE){
    print(sort(unique(FINAL)))    
  }

  
  names(DNA) <- FINAL
  
  if(is.null(fn.new)) {
    fn.new <- fn
  }
  
  cat(paste0("\nRef sequences with taxonomy write in ",fn.new,"\n"))
  
  writeXStringSet(DNA, fn.new, append = FALSE, format = "fasta")
  
  return(list(refseq = DNA))
  
}

# Ajouter id et l'espece a un fichier fasta
add_species <- function(fn, fn.new = NULL, verbose = FALSE, REF = NULL){
  
  #@fn       : nom du fichier fasta d'entré
  #@fn.new   : nom du fichier de sortie si différent 
  
  # Read files
  DNA <- readDNAStringSet(fn)
  
  if(str_detect(names(DNA)[1], ";")){
    stop("La taxonomie semble deja avoir ete ajoutee, action non effectuee", call. = F)
  }
  
  # Extract names
  NAMES1 <- sapply(str_split(names(DNA), pattern=" "), `[`, 1)
  NAMES2 <- paste(sapply(str_split(names(DNA), pattern=" "), `[`, 2),
                  sapply(str_split(names(DNA), pattern=" "), `[`, 3),
                  sep = " ")
  
  if(is.null(REF)){
    NAMES <- paste(NAMES1, NAMES2, sep = " ")
  } else {
    NAMES3 <- vector()
    for(x in NAMES2){
      NAMES3.int <- REF %>% filter(Espece_initial == x) %>% 
                      select(Espece) %>% 
                      pull() %>% 
                      unique()
      if(length(NAMES3.int)!=1){
        NAMES3.int <- x
        cat("no ref for ", x, "\n", sep="")
      }
      
      
      NAMES3 <- c(NAMES3, NAMES3.int)
    }
    
   
    NAMES <- paste(NAMES1, NAMES3, sep = " ")   
  }
  
  if(verbose == TRUE){
    print(unique(NAMES))    
  }

  
  names(DNA) <- NAMES
  
  if(is.null(fn.new)) {
    fn.new <- fn
  }
  
  cat(paste0("\nRef sequences with species write in ",fn.new,"\n"))
  
  writeXStringSet(DNA, fn.new, append = FALSE, format = "fasta")
  
  return(list(refseq = DNA))
  
}

get_uniqueSeq <- function(fn, fn.new = NULL, verbose = F){
  
  #@fn     : nom du fichier fasta d'entré
  #@fn.new : nom du fichier de sortie si différent
  
  # Read ref
  REF <- readDNAStringSet(fn)
  
  # Create a liste of unique occurence
  REF.UNIQUE <- unique(REF)
  
  # Create a dataframe to keep info on duplicated sequences
  REF.DUP <- REF[duplicated(REF)]
  
  if(length(REF.DUP)>0){
    REF.info <- data.frame("Ref" = NA, "Duplicate" = REF.DUP@ranges@NAMES)
    
    for(x in 1:length(REF.DUP)){
      REF.int <- c(REF.DUP[x], REF.UNIQUE)
      REF.info$Ref[x] <- REF.int[duplicated(REF.int)]@ranges@NAMES
    }
    
    if(verbose == TRUE){
      print(REF.info)      
    }

  } else {
    REF.info <- data.frame("Ref", "Duplicate")
    cat("\nNo duplicated sequences\n")
  }
  
  if(is.null(fn.new)) {
    fn.new <- fn
  }
  
  writeXStringSet(REF.UNIQUE, fn.new, 
                  append=FALSE, format="fasta")
  
  cat(paste0("\nUnique ref sequences write in ",fn.new,"\n"))
  
  return(list(refseq = REF.UNIQUE, refdup = REF.info))
} # END function uniqueSeq



# when taxonomy, keep element 7

keep_sp <- function(COL){
  if(str_detect(COL[1], pattern = ";")){
      res <- sapply(str_split(COL, pattern=";"), `[`, 7)
  } else {
      res <- paste(sapply(str_split(COL, pattern=" "), `[`, 2),
                   sapply(str_split(COL, pattern=" "), `[`, 3),
                   sep = " ")
  }
  
  if(is.na(res[1])){
    stop("Il y a deja seulement le nom de l'espece, action non effectuee", call. = F)
  
  
  }
    return(res)
  
}

# apply to each column
capply <- function(df, FUN){
  for(c in 1:ncol(df)){
    COL <- df[,c]
    df[,c] <- FUN(COL)
    
  }
  return(df)
}

# df with 2 columns
keep_dissimilar <- function(df){
  
  new.df <- data.frame(col1 = character(),
                       col2 = character()
  )
  names(new.df) <- names(df)
  
  for(x in 1:nrow(df)){
    if(df[x,1] != df[x,2]){
      new.df <- rbind(new.df, df[x,])
    }
  }
  
  return(unique(new.df))
}

rapid_tree <- function(SEQ, TITLE = "This is a rapide tree"){
  ALIGN <- msa(SEQ, method = "ClustalW")
  TREE.int1 <- msaConvert(ALIGN, type="seqinr::alignment")
  TREE.int2 <- dist.alignment(TREE.int1, "identity")
  TREE.int3 <- njs(TREE.int2)
  
  plot(TREE.int3, main=TITLE)
  
}


# Data --------------------------------------------------------------------

# Obtenir les sequences des amorces
Labo.xls       <- list.files(path = file.path(biodiv.path), pattern = "_LABO_")
Labo.amorces   <- read_excel(file.path(biodiv.path,Labo.xls),sheet="Amorces",na="NA",guess_max=100000)

Labo.amorces 


# Obtenir la liste des fichiers fasta
files.seq      <- c(list.files(file.path(biodiv.path, path.EXTERNE)),
                    list.files(file.path(biodiv.path, path.LABO))
                    )

files.seq [1:6]

# La même liste, mais avec le chemin d'accès
files.seq.wPATH <- c(file.path(biodiv.path,path.EXTERNE, list.files(file.path(biodiv.path, path.EXTERNE))),
                     file.path(biodiv.path,path.LABO, list.files(file.path(biodiv.path, path.LABO)))
                     )
files.seq.wPATH[1:6] 


SP.QC <- read_excel(file.path(data.path,Sample.xl),sheet="Especes",na="NA",guess_max=100000)
SP.QC$QC <- "QC"
SP.QC <- SP.QC %>%  mutate(PC = ifelse(PC == 1, "PC", NA))

SP.QC

# Juste les especes a PC
SP <- SP.QC %>% filter(PC == 1) %>% select(Espece) %>% pull()

# Fichier reference taxo - la boucle permet de verifier si le fichier a bien ete lu
REF <- read_csv(file.path(biodiv.path,"Reference_taxonomie.csv"))
if(ncol(REF) == 1 ) {
  REF <- read_csv2(file.path(biodiv.path,"Reference_taxonomie.csv"))
}
REF

# Ajouter des colonnes a REF pour faire varier la profondeur du jeu de donnees qui m'interresse

REF$All <- "All"
REF <- REF %>% left_join(SP.QC, by = "Espece")


# References --------------------------------------------------------------

# Paramètre de ce que je veux faire
#PARAM <- expand.grid(LOCUS = c("12S"), GROUP = c("All"), PRIMER = c("ecoPrimer"))

PARAM <- expand.grid(LOCUS = c("12S"), GROUP = c("PC", "QC", "All"), PRIMER = c("ecoPrimer", "MiFish"))
PARAM$PRIMER.LAB <- paste(PARAM$LOCUS, str_sub(PARAM$PRIMER, 1,3), sep="-")
PARAM$NMIS <- ifelse(PARAM$PRIMER == "MiFish",5,3)
PARAM

#PARAM <- rbind(PARAM,expand.grid(LOCUS = c("CYTB"), GROUP = c("PC", "QC", "All"), PRIMER = c("Kotcher"), PRIMER.LAB = c("CYTB-Kot"), NMIS = 14))
#PARAM
#PARAM <- expand.grid(LOCUS = c("CYTB"), GROUP = c("PC", "QC", "All"), PRIMER = c("Kotcher"), PRIMER.LAB = c("CYTB-Kot"), NMIS = 14)

# Boucle pour faire rouler le tout

SEQ.INFO <- list()

for(x in 1:nrow(PARAM[,])){
  ID <- paste(as.character(PARAM$GROUP[x]), PARAM$PRIMER.LAB[x], sep="_")
  
  cat(paste0("\n[[Processing ",ID,"]]...\n\n"))

  #Merge sequences
  merge_sequences(LOCUS = c("ADNmt", as.character(PARAM$LOCUS[x])), 
                  GROUP = as.character(PARAM$GROUP[x]), 
                  LIST = files.seq, 
                  LISTwPATH = files.seq.wPATH, 
                  PATH = file.path(path.group), 
                  FN = PARAM$PRIMER.LAB[x],
                  BUFFER = 0, 
                  REF = REF %>% filter(as.character(PARAM$GROUP[x]) == as.character(PARAM$GROUP[x])), 
                  LAB = Labo.amorces %>% filter(NomCommun == PARAM$PRIMER[x]), 
                  KEEP =  FALSE,
                  NMIS = PARAM$NMIS[x],
                  verbose = TRUE)
  
  # Liste des fichiers
  
  FN    <- paste0(as.character(PARAM$GROUP[x]),"_",PARAM$PRIMER.LAB[x],".fasta")
  FN.SP <- paste0(as.character(PARAM$GROUP[x]),"_",PARAM$PRIMER.LAB[x],"_SP.fasta")
  FN.TX <- paste0(as.character(PARAM$GROUP[x]),"_",PARAM$PRIMER.LAB[x],"_taxo.fasta")  
  
  # ADD SP
  add_species(fn = file.path(path.group, FN), fn.new = file.path(path.group, FN), REF = REF)
  
  # Get unique
  SEQ.INFO[[ID]] <- get_uniqueSeq(fn = file.path(path.group, FN), fn.new = file.path(path.group, FN.SP))
  
  # ADD taxo
  add_taxo(fn = file.path(path.group, FN.SP), fn.new = file.path(path.group, FN.TX))
  
  cat(paste0("\nDone!: ", length(SEQ.INFO[[ID]]$refseq), " sequences gardees.\n"))

}


# All_ALL

All.ALL <- readDNAStringSet(files.seq.wPATH)

writeXStringSet(All.ALL, file.path(path.group, "All_ALL.fasta"))

add_species(fn = file.path(path.group, "All_ALL.fasta"), fn.new = file.path(path.group, "All_ALL.fasta"), REF = REF)
add_taxo(fn = file.path(path.group, "All_ALL.fasta"), fn.new = file.path(path.group, "All_ALL_taxo.fasta"))


# Verifier les doublons inter-sp et les ajouter à la main
# Eventuellement ca pourrait être codé aussi



keep_dissimilar(capply(SEQ.INFO[[1]]$refdup, FUN=keep_sp))$Duplicate
SEQ.INFO[[1]]$refdup
SEQ.INFO[[1]]$addSeq <- c("MF621737.1 Salvelinus fontinalis",
                          "MF621744.1 Salvelinus namaycush") 

keep_dissimilar(capply(SEQ.INFO[[2]]$refdup, FUN=keep_sp))$Duplicate
SEQ.INFO[[2]]$refdup
SEQ.INFO[[2]]$addSeq <- c("MF621737.1 Salvelinus fontinalis",
                          "MF621744.1 Salvelinus namaycush",
                          "KP013103.1 Acipenser oxyrinchus",
                          "AP009132.1 Alosa pseudoharengus",
                          "HQ331537.1 Alosa sapidissima",
                          "JQ390060.1 Coregonus clupeaformis",
                          "KM267716.1 Ichthyomyzon fossor",
                          "NC_037014.1 Notropis hudsonius",
                          "KY798500.1 Oncorhynchus mykiss",
                          "AP012101.1 Pimephales notatus",
                          "MF621767.1 Prosopium cylindraceum",
                          "AY216538.1 Hybognathus regius"
                          ) 

keep_dissimilar(capply(SEQ.INFO[[3]]$refdup, FUN=keep_sp))$Duplicate
SEQ.INFO[[3]]$refdup
SEQ.INFO[[3]]$addSeq <- c("MF621737.1 Salvelinus fontinalis",
                          "MF621744.1 Salvelinus namaycush",
                          "KP013103.1 Acipenser oxyrinchus",
                          "AP009132.1 Alosa pseudoharengus",
                          "HQ331537.1 Alosa sapidissima",
                          "KT723026.1 Ammodytes dubius",
                          "JQ390060.1 Coregonus clupeaformis",
                          "KM267716.1 Ichthyomyzon fossor",
                          "NC_037014.1 Notropis hudsonius",
                          "KY798500.1 Oncorhynchus mykiss",
                          "AP012101.1 Pimephales notatus",
                          "MF621767.1 Prosopium cylindraceum",
                          "AY216538.1 Hybognathus regius",
                          "DQ356940.1 Gadus ogac"
                          )

keep_dissimilar(capply(SEQ.INFO[[4]]$refdup, FUN=keep_sp))$Duplicate
SEQ.INFO[[4]]$refdup
SEQ.INFO[[4]]$addSeq <- NA # Ne pas le faire rouler

keep_dissimilar(capply(SEQ.INFO[[5]]$refdup, FUN=keep_sp))$Duplicate
SEQ.INFO[[5]]$refdup
SEQ.INFO[[5]]$addSeq <- c("AB188190.1 Cottus cognatus",
                          "AY372798.1 Percina copelandi",
                          "KM267716.1 Ichthyomyzon fossor") 

keep_dissimilar(capply(SEQ.INFO[[6]]$refdup, FUN=keep_sp))$Duplicate
SEQ.INFO[[6]]$refdup
SEQ.INFO[[6]]$addSeq <- c("AB188190.1 Cottus cognatus",
                          "AY372798.1 Percina copelandi",
                          "KM267716.1 Ichthyomyzon fossor",
                          "DQ536423.1 Lepisosteus osseus") 

keep_dissimilar(capply(SEQ.INFO[[7]]$refdup, FUN=keep_sp))$Duplicate
SEQ.INFO[[7]]$refdup
SEQ.INFO[[7]]$addSeq <- NA # Ne pas le faire rouler

keep_dissimilar(capply(SEQ.INFO[[8]]$refdup, FUN=keep_sp))$Duplicate
SEQ.INFO[[8]]$refdup
SEQ.INFO[[8]]$addSeq <- NA # Ne pas le faire rouler

keep_dissimilar(capply(SEQ.INFO[[9]]$refdup, FUN=keep_sp))$Duplicate
SEQ.INFO[[9]]$refdup
SEQ.INFO[[9]]$addSeq <- NA # Ne pas le faire rouler

# Ajouter les especes larguees (cree un fichier même s'il y en a aucunes)

refSeq.raw.files      <- paste0(names(SEQ.INFO),".fasta")
refSeq.taxo.files     <- paste0(names(SEQ.INFO),"_taxo.fasta")
refSeq.sp.files       <- paste0(names(SEQ.INFO),"_SP.fasta")
refSeq.taxo.new.files <- paste0(names(SEQ.INFO),"_taxo_dup.fasta") 
refSeq.sp.new.files   <- paste0(names(SEQ.INFO),"_SP_dup.fasta") 

for (x in 1:length(SEQ.INFO)){ # pas le 4 car aucun prob
  refSeq.raw   <- readDNAStringSet(file.path(path.group, refSeq.raw.files[x]))
  refSeq.taxo <- readDNAStringSet(file.path(path.group, refSeq.taxo.files[x])) 
  refSeq.sp   <- readDNAStringSet(file.path(path.group, refSeq.sp.files[x]))   
  
  if(is.na(SEQ.INFO[[x]]$addSeq[1])){
    cat(paste0("\nNo sequences added for ", names(SEQ.INFO[x]), "\n"))
    
    refSeq.taxo.new <- refSeq.taxo  
    refSeq.sp.new   <- refSeq.sp      
    
  } else {
    cat(paste0("\nSequences added for ", names(SEQ.INFO[x]), "\n"))
    
    refSeq.toADD <- refSeq.raw[SEQ.INFO[[x]]$addSeq]

    writeXStringSet(refSeq.toADD, "SEQ.fasta.temp", append=FALSE, format = "fasta")
    add_taxo(fn = "SEQ.fasta.temp")
    refSeq.toADDwTAXO <- readDNAStringSet("SEQ.fasta.temp")  

    refSeq.taxo.new <- c(refSeq.taxo,refSeq.toADDwTAXO)
    refSeq.sp.new   <- c(refSeq.sp,refSeq.toADD)  
    
  }

  writeXStringSet(refSeq.taxo.new, file.path(path.group, refSeq.taxo.new.files[x]), append=FALSE, format = "fasta")
  writeXStringSet(refSeq.sp.new, file.path(path.group, refSeq.sp.new.files[x]), append=FALSE, format = "fasta")

  # Enlever le fichier temporaire

}
 
file.remove("SEQ.fasta.temp")  


# Verifier si j'ai tout ce que je veux

for(x in 1:length(SEQ.INFO)){
    print(names(SEQ.INFO)[x])
  CAT <- str_sub(names(SEQ.INFO)[x], 1,2)
  
  if(CAT != "Al"){
      if(file.exists(file.path(path.group, refSeq.sp.new.files[x]))){
        fn <- file.path(path.group, refSeq.sp.new.files[x])        
      } else {
        fn <- file.path(path.group, refSeq.sp.files[x])
      }

      SEQ <- sort(unique(keep_sp(names(readDNAStringSet(fn)))))
      EXP <- SP.QC[!is.na(SP.QC[,CAT]), "Espece"] %>% pull() %>%  sort()
      

      print(setdiff(EXP, SEQ))
  } else {
    print("None!")
  }

}



# Check a la main du comment du pourquoi


DNA <- readDNAStringSet(files.seq.wPATH[str_detect(files.seq.wPATH, "Micropterus dolomieui")])

vmatchPattern(pattern = AMORCE[1], subject = DNA,  max.mismatch = 3)
vmatchPattern(pattern = as.character(reverseComplement(DNAStringSet(AMORCE[2]))), subject = DNA,  max.mismatch = 3)

AMORCE <- Labo.amorces %>% filter(NomCommun == "ecoPrimer") %>% select(Sequence) %>% pull()


length(SEQ)
length(EXP)


names(SEQ.INFO)
str(SEQ.INFO[["PC_12S-eco"]]$refseq)

# Faire un arbre pour voir rapidement ça ressemble a quoi 

pdf(file.path(path.group, "Tree_2018-11-16.pdf"), width=8, height=11,paper='letter') 


DNA1 <-readDNAStringSet(file.path(path.group, "All_12S-eco_SP_dup.fasta"))
DNA2 <-readDNAStringSet(file.path(path.group, "All_12S-eco_taxo_dup.fasta"))

NAMES <- names(DNA1)[str_detect(names(DNA2), pattern = "Cypriniformes")]

rapid_tree(DNA1[NAMES], TITLE = "12S ecoPrimer - Ordre: Cypriniformes")


NAMES <- names(DNA1)[str_detect(names(DNA2), pattern = "Cyprinidae")]

rapid_tree(DNA1[NAMES], TITLE = "12S ecoPrimer - Famille: Cyprinidae")

NAMES <- names(DNA1)[str_detect(names(DNA2), pattern = "Salmonidae")]

rapid_tree(DNA1[NAMES], TITLE = "12S ecoPrimer - Famille: Salmonidae")

DNA1 <-readDNAStringSet(file.path(path.group, "All_12S-MiF_SP_dup.fasta"))
DNA2 <-readDNAStringSet(file.path(path.group, "All_12S-MiF_taxo_dup.fasta"))

NAMES <- names(DNA1)[str_detect(names(DNA2), pattern = "Cypriniformes")]

rapid_tree(DNA1[NAMES], TITLE = "12S MiFish - Ordre: Cypriniformes")


NAMES <- names(DNA1)[str_detect(names(DNA2), pattern = "Cyprinidae")]

rapid_tree(DNA1[NAMES], TITLE = "12S MiFish - Famille: Cyprinidae")

NAMES <- names(DNA1)[str_detect(names(DNA2), pattern = "Salmonidae")]

rapid_tree(DNA1[NAMES], TITLE = "12S MiFish - Famille: Salmonidae")

DNA1 <-readDNAStringSet(file.path(path.group, "All_CYTB-Kot_SP_dup.fasta"))
DNA2 <-readDNAStringSet(file.path(path.group, "All_CYTB-Kot_taxo_dup.fasta"))

## Ne pas faire car trop de séeuqneces trop longue
#NAMES <- names(DNA1)[str_detect(names(DNA2), pattern = "Cypriniformes")]

#rapid_tree(DNA1[NAMES], TITLE = "CytB Kotcher - Ordre: Cypriniformes")


#NAMES <- names(DNA1)[str_detect(names(DNA2), pattern = "Cyprinidae")]

#rapid_tree(DNA1[NAMES], TITLE = "CytB Kotcher  - Famille: Cyprinidae")

NAMES <- names(DNA1)[str_detect(names(DNA2), pattern = "Salmonidae")]

rapid_tree(DNA1[NAMES], TITLE = "CytB Kotcher  - Famille: Salmonidae")

dev.off()



list.files(path.group)

# Compter

length(PC.12Seco$refseq)
length(PC.12SMiF$refseq)

length(All.12Seco$refseq)
length(All.12MiF$refseq)

length(unique(names(All.12Seco$refseq)))
length(unique(names(All.12MiF$refseq)))

# Noms manquant dans MiF
setdiff(names(All.12Seco$refseq), names(All.12SMiF$refseq))

# Noms manquants dans eco
setdiff(names(All.12SMiF$refseq), names(All.12Seco$refseq))


# Old ---------------------------------------------------------------------
# Remettre le code en place car interessant de savoir


sort(unique(REF$Order))

SP.long <- c(REF %>% filter(Espece_initial %in% c(SP %>% pull())) %>% select (Espece_initial) %>% pull(),
             REF %>% filter(Espece_initial %in% c(SP %>% pull())) %>% select (Espece) %>% pull(),
             REF %>% filter(Espece %in% c(SP %>% pull())) %>% select (Espece_initial) %>% pull()
             ) %>% unique()

SP.long 

SP.cypri <- SP.long[which(SP.long %in% c(REF %>% filter(Order %in% c("Cypriniformes", "Cyprinodontiformes")) %>% 
                                                 select (Espece_initial) %>% 
                                                 pull()))]


seq.files <- vector()
seq.res   <- data.frame(SP = SP, N.ADNmt = NA, N.12S = NA,  N.CYTB = NA, N.COI = NA)

for(x in SP){
  seq.x <- list.files(file.path(biodiv.path,"/Sequences/Sequences_EXTERNE"), pattern = x)  
 
  seq.files <- c(seq.files, seq.x)
  
  seq.res[which(seq.res$SP == x), "N.ADNmt"] <- length(seq.x[str_detect(seq.x, pattern = "ADNmt")])
  seq.res[which(seq.res$SP == x), "N.12S"]   <- length(seq.x[str_detect(seq.x, pattern = "12S")])  
  seq.res[which(seq.res$SP == x), "N.CYTB"]   <- length(seq.x[str_detect(seq.x, pattern = "CYTB")])    
  seq.res[which(seq.res$SP == x), "N.COI"]   <- length(seq.x[str_detect(seq.x, pattern = "COI")])
  }

seq.files

seq.res 

seq.res %>% filter(SP %in% SP.cypri)


# Summary -----------------------------------------------------------------

# Combien de sequences par famille

SP.QC.wInfo <- SP.QC %>% dplyr::left_join(REF %>% select(Espece_initial, Class, Order, Family), 
                    by =c("Espece" = "Espece_initial")) #%>% 
                    
SP.QC.wInfo %>% dplyr::group_by(Class, Order, Family) %>% 
                dplyr::summarise(N = length(Espece)) %>% View()


unique(sort(SP.QC.wInfo$Order))


REF %>% filter(Order %in% unique(SP.QC.wInfo$Order)) %>% dplyr::group_by(Order) %>% 
  dplyr::summarise(N = length(Espece)) %>% View()


View(SP.QC)

# Faire un fichier fasta

seq.files.cypri <- vector()


for(x in SP.cypri){
  seq.x <- list.files(file.path(biodiv.path,"/Sequences/Sequences_EXTERNE"), pattern = x)  
  
  seq.files.cypri <- c(seq.files.cypri, seq.x)

}

seq.files.cypri.12S <- c(seq.files.cypri[str_detect(seq.files.cypri, pattern = "12S")],
                         seq.files.cypri[str_detect(seq.files.cypri, pattern = "ADNmt")])

DNA.cypri <- readDNAStringSet(file.path(file.path(biodiv.path,"/Sequences/Sequences_EXTERNE",seq.files.cypri.12S)))

names(DNA.cypri) <- seq.files.cypri.12S

writeXStringSet(DNA.cypri, file.path("./00_Data/03_RefSeq","DNA.cypri.12S.fasta"))
