# Description -------------------------------------------------------------

# Code pour creer la liste de reference
# Audrey Bourret
# 31 octobre 2018 - update january 2019

# Packages ----------------------------------------------------------------

library(readxl)
library(tidyverse)

library(Biostrings)

library(msa)
library(ape)
library(seqinr)
#library(gtools)

#library(devtools)
#devtools::install_github("kassambara/ggpubr")
library(ggpubr)   


#library(ShortRead)
#library(dada2); packageVersion("dada2")

# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}

# Functions ---------------------------------------------------------------

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
Labo.xls       <- list.files(path = get.value("biodiv.path"), pattern = "_LABO_", full.names = T)
Labo.amorces   <- read_excel(Labo.xls,sheet="Amorces",na="NA",guess_max=100000)

Labo.amorces 


# Obtenir la liste des fichiers fasta
files.seq      <- c(list.files(get.value("path.EXTERNE")),
                    list.files(get.value("path.LABO"))
                    )

files.seq [1:6]

# La même liste, mais avec le chemin d'accès
files.seq.wPATH <- c(list.files(get.value("path.EXTERNE"), full.names = T),
                     list.files(get.value("path.LABO"), full.names = T)
                     )

files.seq.wPATH[1:6] 

SP.QC <- read_excel(get.value("Sample.xl"),sheet="Especes",na="NA",guess_max=100000)
SP.QC$QC <- "QC"
SP.QC <- SP.QC %>%  mutate(PC = ifelse(PC == 1, "PC", NA))

SP.QC

# Juste les especes a PC
SP <- SP.QC %>% filter(PC == 1) %>% select(Espece) %>% pull()

# Fichier reference taxo - la boucle permet de verifier si le fichier a bien ete lu
REF <- read_csv(get.value("RefTAXO"))
if(ncol(REF) == 1 ) {
  REF <- read_csv2(get.value("RefTAXO"))
}
REF

# Ajouter des colonnes a REF pour faire varier la profondeur du jeu de donnees qui m'interresse

REF$All <- "All"
REF <- REF %>% left_join(SP.QC, by = "Espece")


# References --------------------------------------------------------------

# Paramètre de ce que je veux faire

PARAM <- expand.grid(LOCUS = c("12S"), GROUP = c("QC", "All"), PRIMER = c("ecoPrimer"))
#PARAM <- expand.grid(LOCUS = c("12S"), GROUP = c("All"), PRIMER = c("ecoPrimer"))
#PARAM <- expand.grid(LOCUS = c("12S"), GROUP = c("PC", "QC", "All"), PRIMER = c("ecoPrimer", "MiFish"))

PARAM$PRIMER.LAB <- paste(PARAM$LOCUS, str_sub(PARAM$PRIMER, 1,3), sep="-")
PARAM$NMIS <- 3
PARAM

PARAM <- rbind(PARAM,expand.grid(LOCUS = c("CYTB"), GROUP = c("QC", "All"), PRIMER = c("Kotcher"), PRIMER.LAB = c("CYTB-Kot"), NMIS = 14))
#PARAM <- rbind(PARAM,expand.grid(LOCUS = c("CYTB"), GROUP = c("PC", "QC", "All"), PRIMER = c("Kotcher"), PRIMER.LAB = c("CYTB-Kot"), NMIS = 14))
#PARAM <- expand.grid(LOCUS = c("CYTB"), GROUP = c("PC", "QC", "All"), PRIMER = c("Kotcher"), PRIMER.LAB = c("CYTB-Kot"), NMIS = 14)
PARAM

# Boucle pour faire rouler le tout

for(x in 1:nrow(PARAM[,])){
  ID <- paste(as.character(PARAM$GROUP[x]), PARAM$PRIMER.LAB[x], sep="_")
  
  cat(paste0("\n[[Processing ",ID,"]]...\n\n"))

  #Merge sequences
  merge_sequences(LOCUS = c("ADNmt", as.character(PARAM$LOCUS[x])), 
                  GROUP = as.character(PARAM$GROUP[x]), 
                  LIST = files.seq, 
                  LISTwPATH = files.seq.wPATH, 
                  PATH = get.value("ref.path"), 
                  FN = PARAM$PRIMER.LAB[x],
                  BUFFER = 0, 
                  REF = REF %>% filter(as.character(PARAM$GROUP[x]) == as.character(PARAM$GROUP[x])), 
                  LAB = Labo.amorces %>% filter(NomCommun == PARAM$PRIMER[x]), 
                  KEEP =  TRUE, # Pour garder ceux qui fit pas
                  NMIS = PARAM$NMIS[x],
                  verbose = TRUE)
  
  # Liste des fichiers
  
  #FN    <- paste0(as.character(PARAM$GROUP[x]),"_",PARAM$PRIMER.LAB[x],".fasta")
  #FN.SP <- paste0(as.character(PARAM$GROUP[x]),"_",PARAM$PRIMER.LAB[x],"_SP.fasta")
  #FN.TX <- paste0(as.character(PARAM$GROUP[x]),"_",PARAM$PRIMER.LAB[x],"_taxo.fasta")  
  
  # ADD SP
  #add_species(fn = file.path(get.value("ref.path"), FN), fn.new = file.path(get.value("ref.path"), FN), REF = REF)
  
  # Get unique
  #SEQ.INFO[[ID]] <- get_uniqueSeq(fn = file.path(get.value("ref.path"), FN), fn.new = file.path(get.value("ref.path"), FN.SP))
  
  # ADD taxo
  #add_taxo(fn = file.path(get.value("ref.path"), FN.SP), fn.new = file.path(get.value("ref.path"), FN.TX))
  
 # cat(paste0("\nDone!: ", length(SEQ.INFO[[ID]]$refseq), " sequences gardees.\n"))

}


# Create alignment to hand edit my stuff

all.files <- list.files(get.value("ref.path"), pattern = "All_", full.names=T) 

ALIGN <- list()

for(x in all.files){
  print(x)  
  name <- str_replace(x, ".fasta", "_align.fasta")
  
  DNA   <- readDNAStringSet(x)
  
  #DNA <- c(DNA, COI.F, COI.R)
  
  if(length(DNA) == 1){
    res <- DNA
    writeXStringSet(res, name)
    
  } else {
    res <- msa(DNA, method = "ClustalW")
    writeXStringSet(res@unmasked, name)
  } 
  
  ALIGN[[x]] <- res

}

# Then open them by hand to have a look... and create an alignement already without my primers

complete.files <- list.files(get.value("ref.path"), pattern = "All_", full.names=T) %>% str_subset("complete.fasta")
complete.files
subset.files <- list.files(get.value("ref.path"), pattern = "QC_", full.names=T) 
subset.files

# Then, select in my subset the ones that I want to keep

for(x in 1:length(complete.files)){

  DNAsub <- readDNAStringSet(subset.files[x])
  DNAcom <- readDNAStringSet(complete.files[x])
  
  # a Few name change because MEGA change my names ...
  names(DNAcom) <- names(DNAcom) %>% str_replace_all("_", " ") %>% str_remove_all(",") %>% str_replace("NC ", "NC_")
  names(DNAsub) <- names(DNAsub) %>% str_replace_all("_", " ") %>% str_remove_all(",") %>% str_replace("NC ", "NC_")
  
  print(length(names(DNAsub)))
  print(length(names(DNAcom)))
  
  new.sub <- intersect(sort(names(DNAsub)), sort(names(DNAcom)))
  
  print(length(new.sub))
  
  writeXStringSet(DNAcom, complete.files[x])
  writeXStringSet(DNAsub[new.sub], subset.files[x] %>% str_replace(".fasta", "_complete.fasta"))

}

# Now its time to changes my files as in my previous loop

complete.files <- list.files(get.value("ref.path"), pattern = "complete.fasta", full.names=T) #%>% str_subset("complete.fasta")
complete.files

SEQ.INFO <- list() # to keep info on duplicated

for(x in complete.files){

  FN    <- x 
  FNnew <- x %>% str_replace("complete", "completewNewNames")  
  FN.SP <- x %>% str_replace("complete", "unique")
  FN.TX <- FN.SP %>% str_replace(".fasta", "_wTAXO.fasta")  

  # ADD SP
  add_species(fn = FN, fn.new = FNnew, REF = REF)

  # Get unique
  SEQ.INFO[[FN.SP]] <- get_uniqueSeq(fn = FNnew, fn.new = FN.SP)
  
  # ADD taxo
  add_taxo(fn = FN.SP, fn.new = FN.TX)

}



# Verifier les doublons inter-sp et les ajouter à la main
# Eventuellement ca pourrait être codé aussi

keep_dissimilar(capply(SEQ.INFO[[1]]$refdup, FUN=keep_sp))$Duplicate
SEQ.INFO[[1]]$refdup

SEQ.INFO[[1]]$addSeq <- c("KM267716.1 Ichthyomyzon fossor",
                          "KP013103.1 Acipenser oxyrinchus",
                          "AF038493.1 Chrosomus neogaeus",
                          
                          "AP012101.1 Pimephales notatus",
                          "NC_037014.1 Notropis hudsonius",
                          "AY216538.1 Hybognathus regius",
                          
                          "KP281293.1 Hypomesus olidus",
                          "JQ390060.1 Coregonus clupeaformis",
                          "MF621767.1 Prosopium cylindraceum",
                          
                          "MF621737.1 Salvelinus fontinalis",
                          "MF621744.1 Salvelinus namaycush",
                          "KP013107.1 Oncorhynchus clarkii",
                          
                          "AP009132.1 Alosa pseudoharengus",
                          "HQ331537.1 Alosa sapidissima",
                          "KX686083.1 Brevoortia tyrannus",
                          
                          "KX929880.1 Anarhichas denticulatus",
                          "MH377821.1 Myoxocephalus scorpius",
                          "KT004432.1 Icelus spatula",
                          
                          "DQ356939.1 Gadus morhua",
                          "DQ356940.1 Gadus ogac",
                          "AM489717.1 Melanogrammus aeglefinus",
                          
                          "KX929917.1 Sebastes mentella",
                          "KX929918.1 Sebastes norvegicus",
                          "KT723026.1 Ammodytes dubius",
                          
                          "EF100182.1 Rajella fyllae",
                          "KF597303.1 Cetorhinus maximus"
                          ) 


keep_dissimilar(capply(SEQ.INFO[[2]]$refdup, FUN=keep_sp))$Duplicate
SEQ.INFO[[2]]$refdup
SEQ.INFO[[2]]$addSeq <- c("DQ678452.1 Sebastes norvegicus"
                          ) 

keep_dissimilar(capply(SEQ.INFO[[3]]$refdup, FUN=keep_sp))$Duplicate
SEQ.INFO[[3]]$refdup
SEQ.INFO[[3]]$addSeq <- c("AF038493.1 Chrosomus neogaeus",
                          "HQ331537.1 Alosa sapidissima",
                          "MF621737.1 Salvelinus fontinalis",
                          
                          "KM267717.1 Ichthyomyzon unicuspis",
                          "KP013107.1 Oncorhynchus clarkii",
                          "KU985081.1 Acipenser fulvescens",
                          
                          "MF621744.1 Salvelinus namaycush",
                          "MF621766.1 Coregonus artedi",
                          "MF621767.1 Prosopium cylindraceum",

                          "MG570409.1 Notropis bifrenatus",
                          "MG570414.1 Notropis heterolepis",
                          "MG570463.1 Alosa aestivalis",
                          
                          "NC_037014.1 Notropis hudsonius"
                          ) 


keep_dissimilar(capply(SEQ.INFO[[4]]$refdup, FUN=keep_sp))$Duplicate
SEQ.INFO[[4]]$refdup
SEQ.INFO[[4]]$addSeq <- NA # Ne pas le faire rouler

# Ajouter les especes larguees (cree un fichier même s'il y en a aucunes)



for (x in names(SEQ.INFO)){ # pas le 4 car aucun prob
  
  # files
  
  refSeq.raw.files      <- x %>% str_replace("unique", "completewNewNames") 
  refSeq.taxo.files     <- x %>% str_replace("unique", "unique_wTaxo") 
  refSeq.sp.files       <- x
  refSeq.taxo.new.files <- refSeq.taxo.files %>% str_replace("unique_wTaxo", "unique_wTaxo_dup")  
  refSeq.sp.new.files   <- refSeq.sp.files %>% str_replace("unique", "unique_dup")  
  
  # DNA
  
  refSeq.raw  <- readDNAStringSet(refSeq.raw.files )
  refSeq.taxo <- readDNAStringSet(refSeq.taxo.files) 
  refSeq.sp   <- readDNAStringSet(refSeq.sp.files  )   
  
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

  writeXStringSet(refSeq.taxo.new, refSeq.taxo.new.files, append=FALSE, format = "fasta")
  writeXStringSet(refSeq.sp.new,   refSeq.sp.new.files, append=FALSE, format = "fasta")

  # Enlever le fichier temporaire

}
 

file.remove("SEQ.fasta.temp")  


# Stats -------------------------------------------------------------------

# Compiler les infos qui pourront etre utiles

final.files <- list.files(get.value("ref.path"), pattern = "unique", full.names=T) 
final.files

# liste des especes manquantes des ref

cat("\n-------------------------\n", 
    "List of missing species\n",
    date(),
    "\n-------------------------\n", 
    file=file.path(get.value("result.ref"), "MissingSP.txt"), 
    append = FALSE, sep = "\n")

for(y in c("PC", "QC")){

  for(x in final.files %>% str_subset("QC") %>% str_subset("unique_dup")){
    print(x)
    
    SEQ <- sort(unique(keep_sp(names(readDNAStringSet(x)))))
    EXP <- SP.QC[!is.na(SP.QC[,y]), "Espece"] %>% pull() %>%  sort()
    
    
    cat("\n", 
        paste(y, x, sep = ": "),
        setdiff(EXP, SEQ), 
        file=file.path(get.value("result.ref"), "MissingSP.txt"), 
        append = TRUE, sep = "\n")
    
  }
}


# N de séquences par fichier

cat("\n-----------------------------\n", 
    "List and size of reference DB\n",
    date(),
    "\n-----------------------------\n", 
    file=file.path(get.value("result.ref"), "REFinfo.txt"), 
    append = FALSE, sep = "\n")

for(x in final.files){
    print(x)
    
    Nseq <- length(readDNAStringSet(x))
    Nsp <- length(unique(keep_sp(names(readDNAStringSet(x)))))
    
    cat("\n", 
        paste(x, sep = ": "),
        paste(Nseq, "sequences", sep= " "), 
        paste(Nsp, "species", sep= " "), 
        file=file.path(get.value("result.ref"), "REFinfo.txt"), 
        append = TRUE, sep = "\n")
    
}


# Liste des espèces qu'on ne peut pas distinguer


cat("\n-----------------------------\n", 
    "List of undistinguable species\n",
    date(),
    "\n-----------------------------\n", 
    file=file.path(get.value("result.ref"), "DUPinfo.txt"), 
    append = FALSE, sep = "\n")

for(x in 1:length(names(SEQ.INFO))){
  
  DATA <- keep_dissimilar(capply(SEQ.INFO[[x]]$refdup, FUN=keep_sp))
  
  cat("\n", 
      paste(names(SEQ.INFO[x]), sep = ": "),
      paste(nrow(DATA),"undistinguable species", sep = " "),
      "",
      file=file.path(get.value("result.ref"), "DUPinfo.txt"), 
      append = TRUE, sep = "\n")

  
  write.table(DATA, 
             file.path(get.value("result.ref"), "DUPinfo.txt"),
             append=T,  row.names = FALSE, col.names = FALSE)
}


# Stats sur les especes par groupe


library(grid) 
library(gridExtra)
library(ggplotify)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


#pdf(file.path(get.value("result.ref"), "testtttt.pdf"), 
#    width = 11, height = 8.5) 


for(x in final.files %>% str_subset("wTaxo_dup")){

  print(x)
  
SEQ <- names(readDNAStringSet(x))
#SEQ

CLASS <- c("Kindom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

DATA <- str_split(SEQ, ";") %>% unlist() %>% matrix(nrow = length(SEQ), ncol=length(CLASS), byrow = T) %>% data.frame()
names(DATA) <- CLASS

theme.legend <- theme(legend.position = "right", 
      legend.text = element_text(size=8),
      legend.justification = "left",
      legend.key.size = unit(10, "points"),
      plot.margin = unit(c(1,1,1,1), "lines"))


GRAPH1 <- DATA %>% group_by(Class, Order, Family) %>% 
                   summarise(N = length(Species)) %>% 
                   ggplot(aes(x="", y = N, fill = Class)) +
                     geom_bar(width = 1, stat = "identity") +
                     coord_polar("y", start=0) +
  guides(fill = guide_legend(ncol = 1, title = NULL)) +
                     labs(title= "N sequences by class")+
  ylab(NULL) + xlab(NULL)+  
  theme_minimal() +
  theme.legend

GRAPH2 <- DATA %>% filter(Class == "Teleostei") %>% 
                   group_by(Class, Order, Family) %>% 
                   summarise(N = length(Species)) %>% 
                   ggplot(aes(x="", y = N, fill = Order)) +
                   geom_bar(width = 1, stat = "identity") +
                   coord_polar("y", start=0) +
  guides(fill = guide_legend(ncol = 1, title = NULL)) +
                   labs(title= "N sequences by order, within Teleostei (class)") +
 ylab(NULL) + xlab(NULL)+  
  theme_minimal() +
  theme.legend




GRAPH3 <- DATA %>% filter(Family == "Cyprinidae") %>% 
                   group_by(Class, Order, Family, Genus, Species) %>% 
                   summarise(N = length(Species)) %>% 
                   ggplot(aes(x="", y = N, fill = Genus)) +
                   geom_bar(width = 1, stat = "identity") +
                   coord_polar("y", start=0) +
                   guides(fill = guide_legend(ncol = 1, title = NULL)) +
                   labs(title= "N sequences by genus, within Cyprinidae (Family)")+
  ylab(NULL) + xlab(NULL)+  
                   theme_minimal() +
  theme.legend 

GRAPH4 <- DATA %>% filter(Family == "Salmonidae") %>% 
                   group_by(Class, Order, Family, Genus, Species) %>% 
                   summarise(N = length(Species)) %>% 
                   ggplot(aes(x="", y = N, fill = Species)) +
                   geom_bar(width = 1, stat = "identity") +
                   coord_polar("y", start=0) +
                  guides(fill = guide_legend(ncol = 1, title = NULL)) +  
                   labs(title= "N sequences by species, within Salmonidae (Family)")+ 
  ylab(NULL) + xlab(NULL)+  
  theme_minimal() +
  theme.legend

GRAPH <- ggarrange(as.ggplot(arrangeGrob(GRAPH1+ theme(legend.position = 'none'), g_legend(GRAPH1), 
                                ncol=2, nrow=1, widths=c(3/4,1/4))),
          as.ggplot(arrangeGrob(GRAPH2+ theme(legend.position = 'none'), g_legend(GRAPH2), 
                                ncol=2, nrow=1, widths=c(3/4,1/4))),
          as.ggplot(arrangeGrob(GRAPH3+ theme(legend.position = 'none'), g_legend(GRAPH3), 
                                ncol=2, nrow=1, widths=c(3/4,1/4))),
          as.ggplot(arrangeGrob(GRAPH4+ theme(legend.position = 'none'), g_legend(GRAPH4), 
                                ncol=2, nrow=1, widths=c(3/4,1/4))),
          labels = LETTERS[1:4],
          legend = "right",
          ncol = 2, nrow = 2
)


ggsave(filename = x %>% str_replace(get.value("ref.path"), get.value("result.ref")) %>%  str_replace("fasta", "pdf"), 
    width = 11, height = 8.5, units = "in",
    plot =  

annotate_figure(GRAPH,
                top = text_grob(x, color = "black", face = "bold", size = 14))
)
#dev.off()

}




