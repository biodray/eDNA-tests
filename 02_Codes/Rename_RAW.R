
# Info --------------------------------------------------------------------

# Simple code to rename Raw files into something little bit easier to work
# with (work better with JAMP)

# Audrey Bourret
# 2018-11-27

# Library -----------------------------------------------------------------

library(tidyverse)
library(readxl)

# Internal functions
# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}

# Data --------------------------------------------------------------------

info.path <- get.value("info.path")
Sample.xl <- get.value("Sample.xl")

raw_unz.path <- get.value("raw_unz.path")

DataSample <- read_excel(get.value("Sample.xl"),sheet="DataSample",na="NA",guess_max=100000)
DataSeq    <- read_excel(get.value("Sample.xl"),sheet="DataSeq",na="NA",guess_max=100000)

# Format de DataSeq avec les données IBIS (p1-A1)

DataSeq <- DataSeq %>% mutate (IbisID = paste0("p",Plaque,"-",Puit)) %>% 
  left_join(DataSample %>% select(SampleID, NomLac, CatSite, Temoins), by = "SampleID")


DataSeq <- DataSeq %>% mutate(Name = ifelse(SeqType %in% c("sample", "dup.sample"), 
                                            paste("Sample",substr(NomLac,1,3), substr(CatSite,1,1), IbisID ,sep="_"),
                                            paste(SampleID,IbisID ,sep="_")))

unique(DataSeq$Name)


# Unzipped files if necessary ---------------------------------------------

# Enlever les fichiers du dossier (unz)

file.remove(list.files(get.value("raw_unz.path"), full.name = T, pattern =".fastq"))

# Copier raw -> unz

file.copy(list.files(get.value("raw.path"), full.name = T, pattern = "fastq"),
          file.path(get.value("raw_unz.path"), list.files(get.value("raw.path"), full.name = F, pattern = "fastq"))
          )

# Unzip (fonctionne sur linux, pê pas sur windows)
system2("gunzip", paste(file.path(list.files(get.value("raw_unz.path"), full.name = T, pattern =".fastq"))))


# Create new files names --------------------------------------------------

# Get the current names of unzipped files

files.names <- list.files(get.value("raw_unz.path"))

new.names <- files.names  %>% stringr::str_remove(pattern = "EP-")

new.names <- paste(sapply(strsplit(new.names, "_"), `[`, 1),
                   sapply(strsplit(new.names, "_"), `[`, 4),
                   sep = "_")

new.names

new.names.final <- vector()

for(x in 1:length(new.names)){
  START  <- sapply(strsplit(new.names[x], "_"), `[`, 1)
  END    <- sapply(strsplit(new.names[x], "_"), `[`, 2) 
  
  START1 <- sapply(strsplit(START, "-"), `[`, 1)
  START2 <- paste(sapply(strsplit(START, "-"), `[`, 2),
                  sapply(strsplit(START, "-"), `[`, 3),
                  sep = "-")
  
  NEW <- DataSeq %>% filter(IbisID == START2) %>% select(Name) %>% pull()
  
  NAME <- paste(START1, NEW, END, sep = "_")
  #print(NAME)
  
  NAME <- paste0(NAME, ".fastq")
  
  new.names.final <- c(new.names.final, NAME) 
  
}

new.names.final



# Change files names ------------------------------------------------------

# Enlever les fichiers précédents (unz)
file.remove(list.files(get.value("raw_unz_rename.path"), full.name = T, pattern =".fastq"))

# Copier unz -> unz_rename (can take some times)

file.copy(list.files(get.value("raw_unz.path"), full.name = T, pattern = "fastq"),
          file.path(get.value("raw_unz_rename.path"), list.files(get.value("raw_unz.path"), full.name = F, pattern = "fastq"))
          )

# Renommer unz_rename

file.rename(file.path(get.value("raw_unz_rename.path"), files.names),
            file.path(get.value("raw_unz_rename.path"), new.names.final)
            )
