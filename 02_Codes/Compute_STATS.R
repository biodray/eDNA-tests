
# Info --------------------------------------------------------------------

# Compute stats for this dataset
# 
# Audrey Bourret
# 2019-01-14 
#

# Library -----------------------------------------------------------------

library(tidyverse)
library(readxl)

# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}

library(gtools)    # for mixedsort

#library(devtools)
#devtools::install_github("kassambara/ggpubr")
library(ggpubr)    # on github - for nice graphs

#Function "Not IN"
`%nin%` = Negate(`%in%`)

# For modelisation
library(lme4)
library(effects)

# Data --------------------------------------------------------------------

# Sample INFO

DataSample <- read_excel(get.value("Sample.xl"),sheet="DataSample",na="NA",guess_max=100000)
DataSample

DataSeq    <- read_excel(get.value("Sample.xl"),sheet="DataSeq",na="NA",guess_max=100000)
DataSeq

DataSeq <- DataSeq %>% mutate (IbisID = paste0(Plaque,".",Puit)) %>% 
  left_join(DataSample %>% select(SampleID, NomLac, CatSite, Nsite, Volume, CorrFiltre), by = "SampleID")
DataSeq

# Lac info

LacSample <- read_excel(get.value("Lac.xl"),sheet="NomLac",na="NA",guess_max=100000)
LacSample

LacInv1 <- read_excel(get.value("Lac.xl"),sheet="Inventaire1996-2003",na="NA",guess_max=100000)
LacInv1 <- na.replace(LacInv1, 0)


LacInv1 <- LacInv1 %>% inner_join(LacSample %>% select(-c(AbrLac)), by = "InvLac") %>% 
           gather(names(.) %>% str_subset(" "), key= "Espece", value = "Presence")

#LacInv1 %>% group_by(NomLac) %>% summarise(Nsp = sum(Presence)) %>% View()

LacInv2 <- read_excel(get.value("Lac.xl"),sheet="Inventaire2018",na="NA",guess_max=100000)
LacInv2 <- LacInv2 %>% gather(names(.) %>% str_subset(" "), key = Espece, value = Value)


LacInv2 %>% spread(key = Mesure, value = Value) %>%  ggplot(aes(x = CPUE, y = BPUE, col = Espece, shape = Peche))+ geom_point(size = 2) +theme_classic() + facet_wrap(~Espece)


LacInv2 <- LacInv2 %>% left_join(LacSample %>% select(NomLac, Surface, Volume)) %>% 
                       mutate(Density.s = Value / Surface,
                              Density.v = Value / Volume)

LacPeche <- read_excel(get.value("Lac.xl"),sheet="StatPeche",na="NA",guess_max=100000)
LacPeche <- LacPeche %>% gather(names(.) %>% str_subset("PUE"), key = "Mesure", value = "Value") 


LacInv1 %>% filter(Presence == 0) %>% pull(Espece) %>% unique()


# Mock community info

Mock.dat <- read_excel(get.value("Sample.xl"),sheet="Mock",na="NA",guess_max=100000)

Mock.dat <- Mock.dat %>% gather(paste0("Mix",1:6), key="Mix", value = "Vol") %>% 
                         mutate(DNA = Vol * 10) # DNA concentration = 10 ng/ul

Mock.dat


# Ref sequences

# Fichier reference taxo - la boucle permet de verifier si le fichier a bien ete lu

# Get the server version if connected, otherwise get a local version
if(file.exists(get.value("RefTAXO"))){
  REF <- read_csv(get.value("RefTAXO"), locale = locale(encoding = "ISO-8859-1"))
  if(ncol(REF) == 1 ) {
    REF <- read_csv2(get.value("RefTAXO"), locale = locale(encoding = "ISO-8859-1"))
  }
} else {
  REF <- read_csv("00_Data/00_FileInfos/Reference_taxonomie.csv", locale = locale(encoding = "ISO-8859-1"))
  if(ncol(REF) == 1 ) {
    REF <- read_csv2("00_Data/00_FileInfos/Reference_taxonomie.csv", locale = locale(encoding = "ISO-8859-1"))
  }
}

REF



# Summary by SP

load(get.value("CORRECTEDtable.data"))
load(get.value("ALLtable.data"))
load(get.value("BasicStats.data"))

ls() %>% str_subset("tab")

# To compute haplotypes - byID.SP

for(x in ls() %>% str_subset("wTAXO")){
  DATA <- get(x)
  
  SAMPLE <- names(DATA %>% select(-c(ID, Assign)))
  
  NEW <- DATA %>% gather(SAMPLE, key= "Sample", value = Read) %>% 
    group_by(ID, Assign, Sample) %>% 
    summarise(N = sum(Read)) %>% 
    mutate(Level = str_count(Assign, ";"))
  
  assign(x %>%str_replace("wTAXO", "byID.SP"), NEW)
  
}

ls() %>% str_subset("byID.SP")

# blast99

load(get.value("Blast99.data"))


#save(file = "Compute_STATS.data", list = ls())
#load("Compute_STATS.data")

# Basic stats -------------------------------------------------------------

library(R.utils)

reads.tab <- expand.grid(Sample = list.files(get.value("raw_unz_rename.path"), pattern = "12s", full.names = F) %>% 
                                  str_subset("R1") %>% str_remove("12s_") %>% str_remove("_R1.fastq"),
                         Locus = c("12s", "cytB"), 
                         Sens = c("R1","R2"))

for(x in 1:nrow(reads.tab)){
  FILES <- file.path(get.value("raw_unz_rename.path"), paste(reads.tab[x,"Locus"], 
                                                              reads.tab[x,"Sample"], 
                                                              paste0(reads.tab[x,"Sens"], ".fastq"), sep = "_"))
  
  print(FILES)
  
  reads.tab[x, "Raw"] <- countLines(FILES)/4
  
}

# Nreads total
sum(reads.tab$Raw)

# Nreads by locus

reads.tab %>% group_by(Locus) %>% 
              summarise(N =  sum(Raw)/1000000)


# After the filtration step # TO DO - check the real name

list.files(get.value("filt_dada2.path"))[1]


for(x in 1:nrow(reads.tab)){
  FILES <- file.path(get.value("filt_dada2.path"), paste(reads.tab[x,"Locus"], 
                                                             reads.tab[x,"Sample"], 
                                                             paste0(reads.tab[x,"Sens"], "_cut_filt.fastq"), 
                                                             sep = "_"))
   if(file.exists(FILES)){
   
     print(FILES)
         
   reads.tab[x, "Filt"] <- countLines(FILES)/4
   } else {
     reads.tab[x, "Filt"] <-      0
   }

  
}

View(reads.tab)

# N ASV/OTU

for(x in 1:nrow(reads.tab)){
  
  SAMPLE <- paste(reads.tab[x, "Locus"], 
                   reads.tab[x, "Sample"] %>% str_replace("-", "."),
                   sep = "_")
  
  if(reads.tab[x, "Locus"] == "12s"){
    
    N.ASV <- sum(ASVtab.12s[,SAMPLE])
    N.OTU <- sum(OTUtab.12s[,SAMPLE])
    
  } else if (reads.tab[x, "Locus"] == "cytB"){
    
    if(reads.tab[x, "Sens"] == "R1"){
      
     N.ASV <- sum(ASVtab.cytB.R1[,paste(SAMPLE, "R1", sep = "_")])
     N.OTU <- sum(OTUtab.cytB.R1[,paste(SAMPLE, "R1", sep = "_")])
     
    } else if(reads.tab[x, "Sens"] == "R2"){
      
     N.ASV <- sum(ASVtab.cytB.R2[,paste(SAMPLE, "R2", sep = "_")])
     N.OTU <- sum(OTUtab.cytB.R2[,paste(SAMPLE, "R2", sep = "_")])
     
     } 
  }
  
 reads.tab[x, "ASV"] <- N.ASV  
 reads.tab[x, "OTU"] <- N.OTU    
  
}




# N ASV/OTU corrected

for(x in 1:nrow(reads.tab)){
  
    SAMPLE <- paste(reads.tab[x, "Locus"], 
                  reads.tab[x, "Sample"] %>% str_replace("-", "."),
                  sep = "_")
  
    if(SAMPLE %>% str_detect("Sample")){
    
  if(reads.tab[x, "Locus"] == "12s"){
    
    N.ASV <- sum(ASVtab.12s.cor[,SAMPLE])
    N.OTU <- sum(OTUtab.12s.cor[,SAMPLE])
    
  } else if (reads.tab[x, "Locus"] == "cytB"){
    
    if(reads.tab[x, "Sens"] == "R1"){
      
      N.ASV <- sum(ASVtab.cytB.R1.cor[,paste(SAMPLE, "R1", sep = "_")])
      N.OTU <- sum(OTUtab.cytB.R1.cor[,paste(SAMPLE, "R1", sep = "_")])
      
    } else if(reads.tab[x, "Sens"] == "R2"){
      
      N.ASV <- sum(ASVtab.cytB.R2.cor[,paste(SAMPLE, "R2", sep = "_")])
      N.OTU <- sum(OTUtab.cytB.R2.cor[,paste(SAMPLE, "R2", sep = "_")])
      
    } 
  }
  
    } else{
      
      N.ASV <- NA
      N.OTU <- NA
    }
    
    
  reads.tab[x, "ASV.cor"] <- N.ASV  
  reads.tab[x, "OTU.cor"] <- N.OTU    
  
}


DataSeq


#save(file = get.value("BasicStats.data"), 
#     list = "reads.tab")

load(get.value("BasicStats.data"))

View(reads.tab) 

reads.tab %>% filter(!is.na(ASV.cor),
                     Locus == "12s") %>% 
              mutate(IbisID = sapply(str_split(Sample, "_"), `[`, 4) %>% str_remove("p") %>% str_replace("-", ".")) %>% 
  left_join(DataSeq %>% select(IbisID, SeqType)) %>%
  filter(SeqType == "sample") %>%
  group_by(Locus, Sens) %>% 
  summarise(Nraw = sum(Raw),
            Nfilt = sum(Filt),
            NASV = sum(ASV),
            NASV.cor = sum(ASV.cor))
  



# Nreads by locus

reads.tab %>% group_by(Locus, Sens) %>% 
  summarise(Nraw =  sum(Raw),
            Nfilt = sum(Filt),
            N.ASV = sum(ASV),
            N.OTU = sum(OTU),
            PRawToFilt = 1 - (Nfilt/Nraw))

reads.tab %>% mutate(ASV.filt = Filt - ASV,
                     ASV.raw = Raw - Filt) %>%
              gather(ASV, ASV.filt, ASV.raw, key = "ASV", value = N.ASV) %>% 
              mutate(P.ASV =N.ASV / Raw) %>% 
              filter(Sens == "R1") %>% 
  
  ggplot(aes(x = Sample, y = P.ASV, fill = ASV)) +
              geom_bar(stat = "identity") +
              scale_y_continuous(name = "N reads") +
              facet_grid(Locus ~.) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  

# Negative samples --------------------------------------------------------


tneg.stat <- function(tab) {
  
  DATA <- tab %>% select(ID,
                         names(.) %>% str_subset("_Tneg"),
                         names(.) %>% str_subset("_T[:digit:]")) %>% 
    gather(c(names(.) %>% str_subset("_Tneg"),
             names(.) %>% str_subset("_T[:digit:]")),
           key = Sample, value = N)
  
  cat("\nN unique haplotypes:",
      DATA %>% filter(N >= 1) %>% select(ID) %>% pull() %>%  unique() %>% length(),
      
      "\nN total reads:",
      DATA %>% select(N) %>% pull() %>%  sum(),
      
      "\nMedian reads by sample:",
      DATA %>% group_by(Sample) %>% summarise(N = sum(N)) %>%  select( N) %>% pull() %>% median(),
      
      "\nMean reads by sample:",
      DATA %>% group_by(Sample) %>% summarise(N = sum(N)) %>%  select( N) %>% pull() %>% mean(),
      
      
      "\nMin reads by sample:",
      DATA %>% group_by(Sample) %>% summarise(N = sum(N)) %>%  select( N) %>% pull() %>% min(),
      
      "\nMax reads by sample:",
      DATA %>% group_by(Sample) %>% summarise(N = sum(N)) %>%  select( N) %>% pull() %>% max(),
      
      "\nMedian reads by haplotype:",
      DATA %>% group_by(ID) %>% summarise(N = sum(N)) %>%  filter(N >=1) %>% select( N) %>% pull() %>% median(),
      
      "\nMean reads by haplotype:",
      DATA %>% group_by(ID) %>% summarise(N = sum(N)) %>%  filter(N >=1) %>% select( N) %>% pull() %>% mean(),
      
      
      "\nMin reads by haplotype:",
      DATA %>% group_by(ID) %>% summarise(N = sum(N)) %>%  filter(N >=1) %>% select( N) %>% pull() %>% min(),
      
      "\nMax reads by haplotype:",
      DATA %>% group_by(ID) %>% summarise(N = sum(N)) %>%  filter(N >=1) %>% select( N) %>% pull() %>% max(),
      
      
      sep = "\n")
  
}

tneg.stat(ASVtab.12s)
tneg.stat(OTUtab.12s)

tneg.stat(ASVtab.cytB.R1)
tneg.stat(OTUtab.cytB.R1)


tneg.cor.stat <- function(tab, tab.cor) {
  
  N1 <- tab %>% select(ID,
                       names(.) %>% str_subset("Sample")) %>% 
    gather(names(.) %>% str_subset("Sample"),
           key = Sample, value = N) %>% 
    pull(N) %>% sum()
  
  N2 <- tab.cor %>% select(ID,
                           names(.) %>% str_subset("Sample")) %>% 
    gather(names(.) %>% str_subset("Sample"),
           key = Sample, value = N) %>% 
    pull(N) %>% sum()
  
  cat("\nP sequence removed:",
      1 - (N2/N1),
      "\nN sequences final",
      N2,
      sep = "\n")
  
  
}


tneg.cor.stat(ASVtab.12s, ASVtab.12s.cor)
tneg.cor.stat(OTUtab.12s, OTUtab.12s.cor)

tneg.cor.stat(ASVtab.cytB.R1, ASVtab.cytB.R1.cor)
tneg.cor.stat(OTUtab.cytB.R1, OTUtab.cytB.R1.cor)



# Stats on assignation ----------------------------------------------------



Summary.Assign.cor <- data.frame(Assign = character(), 
                                 Level = integer(), 
                                 Nread = numeric(),
                                 Nhaplo = numeric(),  
                                 Nsample = integer(),  
                                 Data = character())


for(x in ls() %>% str_subset("cor.byID.SP")){

    DATA <- get(x) %>% filter(N >=1) %>% 
      group_by(Assign, Level) %>% 
    summarise(Nread = sum(N),
              Nhaplo = length(unique(ID)),
              Nsample = length(unique(Sample))) %>% 
    mutate(Data = x %>% str_remove(".cor.byID.SP")) %>% 
    as.data.frame()
  
  
  Summary.Assign.cor <- rbind(Summary.Assign.cor, DATA)
  
}


head( Summary.Assign.cor)

Summary.Assign.cor %>% filter(Level > 2, 
                            str_detect(.$Assign, "Root;Chordata;Teleostei;") == T,
                            Nread > 0) %>% 
  mutate(Assign = str_remove(Assign, "Root;Chordata;Teleostei;"),
         Method = str_sub(Data, 1, 3),
         Locus = Data %>% str_remove(paste0(Method,"tab.")),
         Nread1000 = Nread / 1000) %>% 
  arrange(Assign, Data) %>% 
  ggplot(aes(y= Nread, x = Assign, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Locus ~ ., scale = "free") +
  #scale_y_continuous(trans = "log") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ coord_flip()



Summary.Assign.cor %>% #filter(Level > 2, 
  #        str_detect(.$Assign, "Root;Chordata;Teleostei;") == T) %>% 
  mutate(#Assign = str_remove(Assign, "Root;Chordata;Teleostei;"),
    Method = str_sub(Data, 1, 3),
    Locus = Data %>% str_remove(paste0(Method,"tab."))) %>% 
  arrange(Assign, Data) %>% 
  ggplot(aes(y= Nsample, x = Assign, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Locus ~ .) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Graphe N read par niveau par dataset

Assign.level <- Summary.Assign.cor %>%  mutate(Cat = ifelse(str_detect(.$Assign, "Teleostei") == T,
                                                  "Teleostei",
                                                  "Other"),
                                     Method = str_sub(Data, 1, 3),
                                     Locus = Data %>% str_remove(paste0(Method,"tab."))) %>% 
  group_by(Level, Data, Locus, Method) %>% 
  summarise(Sequences = sum(Nread),
            Haplotypes = sum(Nhaplo))  %>% 
  gather(Sequences, Haplotypes, key = Stat, value = N)



Assign.level %>% filter(Locus %in% c("12s", "cytB.R1")) %>% 
                 group_by(Stat, Locus, Level, Method) %>% 
                 summarise(N = sum(N)) %>% 
                 left_join(Assign.level %>% group_by(Stat, Locus, Method) %>% 
                                            summarise(Ntot = sum(N))) %>% 
                 mutate(Perc = round(N / Ntot,3)) %>% View()

#library(scales)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

taxo.graph <- function(DATA) {
  DATA %>% 
  ggplot(aes(x = "", y = N, fill = as.factor(Level))) + 
  geom_bar(position = "fill", stat = "identity", colour = "black") +
  coord_polar("y", start=0, direction = -1) +
  scale_fill_brewer(palette = "GnBu", 
                    limits = c(0:6),
                    labels = c("Aucun", "Phylum", "Classe", "Ordre", "Famille", "Genre", "Espèce") ) +
  guides(fill = guide_legend(nrow = 2, title = NULL)) +
  labs(title= NULL,
       subtitle = NULL)+
  ylab(NULL) + xlab(NULL)+
  facet_grid(Stat ~ Method) + 
  blank_theme +  theme(axis.text.x=element_blank(),
                       #strip.background = element_rect(colour = "black", fill = "gray"),
                       strip.text.x = element_text(colour = "black", face = "bold"),
                       strip.text.y = element_text(colour = "black", face = "bold"),
                       panel.spacing = unit(1, "lines")) 
  #geom_text(aes(y = N/6 + c(0, cumsum(N)[-length(N)]), 
  #              label = percent(N/100)), size=5)

}  

ggarrange(taxo.graph(Assign.level %>% filter(Locus %in% c("12s"))) + ggtitle("12s"),
          taxo.graph(Assign.level %>% filter(Locus %in% c("cytB.R1")))+ ggtitle("cytB (R1)"),
          labels = LETTERS[1:2],
          nrow=1, ncol=2, common.legend = T,  legend = "bottom")


# Nread.scale <- TEST %>% group_by(Data) %>% 
#   summarise(Nread = sum(Nread)) %>%
#   mutate(Scale = Nread / max(.$Nread)) %>% 
#   select(Data, Scale)
# 
# TEST %>% left_join(Nread.scale ) %>% 
#   ggplot(aes(x=Scale/2, y = Nread, fill = Level, width = Scale)) +
#   geom_bar(position = "fill", stat = "identity") +
#   coord_polar("y", start=0) +
#   guides(fill = guide_legend(ncol = 1, title = "Depth")) +
#   labs(title= "Proportion of reads assigned to each taxonomic depth",
#        subtitle = "From no assignment (0) to species (6)")+
#   ylab(NULL) + xlab(NULL)+
#   facet_grid(Method ~ Locus) + 
#   theme_bw() 






# Comparison cytB R1 and R2 -----------------------------------------------

# On the assignation ... not on the OTU/ASV because not the same!!!
# On Mock community, and on overall data


R1R2.graph <- function(tab.R1, tab.R2) {
  
  # Need the bySP version of tab
  
  Comp.res <- rbind(tab.R1, tab.R2) %>% 
    select(-Level) %>% 
    # Limit to MIX samples
    filter(str_detect(Sample, "Mix")) %>% 
    mutate(Mix = sapply(str_split(Sample, "_"), `[`, 2),
           Puit = sapply(str_split(Sample, "p[:digit:]."), `[`, 2) %>% str_remove("_R[:digit:]"), 
           R = sapply(str_split(Sample, "_"), `[`, 4)) %>% 
    left_join(DataSeq %>% filter(SeqType %in% c("mix", "dup.mix")) %>% select(Puit, SeqType)) %>% 
    # Remove technical duplicates
    filter(SeqType == "mix") %>% 
    select(Assign, N, Mix, R) %>% 
    spread(key = R, value = N, fill = 0) %>% 
    mutate(KEEP = R1 + R2) %>% 
    filter(KEEP > 0) %>% select(-KEEP)
  
  Comp.res
  
  COR.RES <- cor.test(Comp.res$R1, Comp.res$R2, method = "spearman")
  
  Comp.res %>% mutate(R1 = R1 + 1,
                      R2 = R2 +1) %>%  
    ggplot(aes(x = R1, y = R2, col = Mix)) + 
    # Add a 1:1 line
    geom_abline(intercept = 0, colour = "gray", slope = 1, linetype = 1, show.legend = FALSE ) +
    geom_jitter(width = 0.05, height = 0.05) +
    scale_x_continuous(name = "N reads - R1", trans = "log10") +
    scale_y_continuous(name = "N reads - R2",trans = "log10") +
    guides(col = guide_legend(title = NULL, nrow = 3)) +
    annotate("text" , x = 10, y = 10000, 
             label = paste("rho =", round(COR.RES$estimate, 3), ";", 
                           "P", ifelse(COR.RES$p.value < 0.001, "< 0.001", paste("=", round(COR.RES$p.value, 3))))) +
    theme_classic() + 
    theme(legend.justification=c(1 ,0), 
          legend.position=c(1,0),
          legend.background = element_rect(fill = "gray90"))
  
  
}

# Make a cool graph

GRAPH <- ggarrange(R1R2.graph(ASVtab.cytB.R1.bySP, ASVtab.cytB.R2.bySP) + guides(col=FALSE) + ggtitle("cytB - ASV") ,
                   R1R2.graph(OTUtab.cytB.R1.bySP, OTUtab.cytB.R2.bySP) + ggtitle("cytB - OTU"), 
                   
                   labels = LETTERS[1:2],
                   ncol = 2, nrow=1
                   #common.legend = T, legend = "bottom"
                   )

GRAPH

# Save graphes

ggsave(filename = file.path(get.value("result.FINAL"), "CytB.R1R2.pdf"),
       width = 8, height = 4,
       plot = GRAPH
       )

ggsave(filename = file.path(get.value("result.FINAL"), "CytB.R1R2.png"),
       width = 8, height = 4,
       plot = GRAPH
)



# Technical replicates ------------------------------------------------------

# Function to create graphes fro mix samples only, on the entire dataset

Dup.graph <- function(tab){
  
  Dup.res <- tab %>% select(ID, names(.) %>% str_subset("Mix") ) %>%
    gather(names(.) %>% str_subset("Mix"), key = "sample", value = "Nread") %>% 
    mutate(Mix = sapply(str_split(sample, "_"), `[`, 2),
           Puit = sapply(str_split(sample, "p[:digit:]."), `[`, 2) %>% str_remove("_R[:digit:]")) %>% 
    #left_join(Mock.final) %>% 
    #filter(!is.na(DNAfinal)) %>% 
    #mutate(Nread.log = log10(Nread + 1)) %>% 
    left_join(DataSeq %>% filter(SeqType %in% c("mix", "dup.mix")) %>% select(Puit, SeqType)) %>% 
    select(-sample, -Puit) %>% 
    spread(key = SeqType, value = Nread) %>% 
    mutate(KEEP = dup.mix + mix) %>% 
    filter(KEEP > 0) %>% select(-KEEP)
  
  Dup.res
  
  COR.RES <- cor.test(Dup.res$dup.mix, Dup.res$mix, method = "spearman")
  
  Dup.res %>% mutate(dup.mix = dup.mix + 1,
                     mix = mix + 1) %>%  
    ggplot(aes(x = mix, y = dup.mix, col = Mix)) + 
    # Add a 1:1 line
    geom_abline(intercept = 0, colour = "gray", slope = 1, linetype = 1, show.legend = FALSE ) +
    geom_jitter(width = 0.05, height = 0.05) +
    scale_x_continuous(name = "N lectures - échantillon 1", trans = "log10") +
    scale_y_continuous(name = "N lectures - échantillon 2", trans = "log10") +
    guides(col = guide_legend(title = NULL, nrow = 3)) +
    annotate("text" , x = 20, y = 10000, 
             label = paste("rho =", round(COR.RES$estimate, 2), ";", 
                           "P", ifelse(COR.RES$p.value < 0.001, "< 0.001", paste("=", round(COR.RES$p.value, 3))))) +
    theme_classic() + 
    theme(legend.justification=c(1 ,0), 
          legend.position=c(1,0),
          legend.background = element_rect(fill = "gray90"))
  
}


GRAPH <- ggarrange(Dup.graph(ASVtab.12s.wTAXO) + guides(col=FALSE) + ggtitle("12S - ASV") ,
                   Dup.graph(OTUtab.12s.wTAXO) + guides(col=FALSE)+ ggtitle("12S - OTU"), 
                   Dup.graph(ASVtab.cytB.R1.wTAXO) + guides(col=FALSE)+ ggtitle("cytB (R1) - ASV"),
                   Dup.graph(OTUtab.cytB.R2.wTAXO)+ ggtitle("cytB (R1) - OTU"),
                   #Dup.graph(ASVtab.cytB.R1.wTAXO),
                   #Dup.graph(OTUtab.cytB.R2.wTAXO),
                   labels = LETTERS[1:4],
                   ncol = 2, nrow=2
                   #common.legend = T, legend = "bottom"
                   )

print(GRAPH)

ggsave(filename = file.path(get.value("result.FINAL"), "DuplicatedSamples.png"),
       width = 8, height = 8,
       plot = GRAPH
       )

ggsave(filename = file.path(get.value("result.FINAL"), "DuplicatedSamples.pdf"),
       width = 8, height = 8,
       plot = GRAPH
)

# Comparison ASV and OTU -----------------------------------------------

# On the assignation ... as for R1 vs R2 comparison
# On Mock community, and on overall data


ASVOTU.graph <- function(tab.ASV, tab.OTU) {
  
  # Need the bySP version of tab
  
  
  #tab.ASV <- ASVtab.12s.bySP
  #tab.OTU <- OTUtab.12s.bySP
  
  Comp.res <- rbind(tab.ASV %>% mutate(Method = "ASV"), tab.OTU %>% mutate(Method = "OTU")) %>% 
    select(-Level) %>% 
    # Keep only fishes
    filter(str_detect(Assign, "Teleostei")) %>% 
    mutate(Order = sapply(str_split(Assign, ";"), `[`, 5)) %>% 
    filter(!is.na(Order), Order != "Gadidae" ) %>% 
    # Limit to MIX samples
    filter(str_detect(Sample, "Mix")) %>% 
    mutate(Mix = sapply(str_split(Sample, "_"), `[`, 2),
           Puit = sapply(str_split(Sample, "p[:digit:]."), `[`, 2) %>% str_remove("_R[:digit:]")) %>% 
    left_join(DataSeq %>% filter(SeqType %in% c("mix", "dup.mix")) %>% select(Puit, SeqType)) %>% 
    # Remove technical duplicates
    filter(SeqType == "mix") %>% 
    select(Assign, Order,N, Method, Mix) %>% 
    spread(key = Method, value = N, fill = 0) %>% 
    mutate(KEEP = ASV + OTU) %>% 
    filter(KEEP > 0) %>% select(-KEEP)
  
  Comp.res
  
  COR.RES <- cor.test(Comp.res$ASV, Comp.res$OTU, method = "spearman")
  
  Comp.res %>% mutate(ASV = ASV + 1,
                      OTU = OTU +1) %>%  
    ggplot(aes(x = ASV, y = OTU, col = Order)) + 
    # Add a 1:1 line
    geom_abline(intercept = 0, colour = "gray", slope = 1, linetype = 1, show.legend = FALSE ) +
        geom_jitter(width = 0.05, height = 0.05) +
    scale_x_continuous(name = "N lectures - Assignation ASV", trans = "log10", limits= c(0.8,NA)) +
    scale_y_continuous(name = "N lectures - Assignation OTU",trans = "log10", limits= c(0.8,NA)) +
    guides(col = guide_legend(title = "Famille", ncol = 1, title.hjust = 0.5)) +
    annotate("text" , x = 20, y = 10000, 
             label = paste("rho =", round(COR.RES$estimate, 2), ";", 
                           "P", ifelse(COR.RES$p.value < 0.001, "< 0.001", paste("=", round(COR.RES$p.value, 3))))) +
    theme_classic() + 
    theme(#legend.justification=c(1 ,0), 
          #legend.position=c(1,0),
          legend.background = element_rect(fill = "gray90"))
  
  
}

ASVOTU2.graph <- function(tab.ASV, tab.OTU) {
  
  # Need the bySP version of tab
  
  
  #tab.ASV <- ASVtab.12s.cor.bySP
  #tab.OTU <- OTUtab.12s.cor.bySP
  
  Comp.res <- rbind(tab.ASV %>% mutate(Method = "ASV"), tab.OTU %>% mutate(Method = "OTU")) %>% 
    select(-Level) %>% 
    # Keep only fishes
    filter(str_detect(Assign, "Teleostei")) %>% 
    mutate(Order = sapply(str_split(Assign, ";"), `[`, 5)) %>% 
    filter(!is.na(Order), Order != "Gadidae" ) %>% 
    select(Assign, Sample, Order, N, Method) %>%  
    spread(key = Method, value = N, fill = 0) %>% 
    mutate(KEEP = ASV + OTU) %>% 
    filter(KEEP > 0) %>% select(-KEEP)
  
  Comp.res
  
  COR.RES <- cor.test(Comp.res$ASV, Comp.res$OTU, method = "spearman")
  
  Comp.res %>% mutate(ASV = ASV + 1,
                      OTU = OTU +1) %>%  
    ggplot(aes(x = ASV, y = OTU, col = Order)) + 
    # Add a 1:1 line
    geom_abline(intercept = 0, colour = "gray", slope = 1, linetype = 1, show.legend = FALSE ) +
    geom_jitter(width = 0.05, height = 0.05) +
    scale_x_continuous(name = "N lectures - Assignation ASV", trans = "log10", limits= c(0.8,NA)) +
    scale_y_continuous(name = "N lectures - Assignation OTU",trans = "log10", limits= c(0.8,NA)) +
    guides(col = guide_legend(title = "Famille", ncol = 1, title.hjust = 0.5)) +
    annotate("text" , x = 20, y = 10000, 
             label = paste("rho =", round(COR.RES$estimate, 2), ";", 
                           "P", ifelse(COR.RES$p.value < 0.001, "< 0.001", paste("=", round(COR.RES$p.value, 3))))) +
    theme_classic() + 
    theme(#legend.justification=c(1 ,0), 
          #legend.position=c(1,0),
          legend.background = element_rect(fill = "gray90"))
  
  
}

# Make a cool graph

GRAPH <- ggarrange(ASVOTU.graph(ASVtab.12s.bySP, OTUtab.12s.bySP) + guides(col=FALSE) + ggtitle("12S") ,
                   ASVOTU.graph(ASVtab.cytB.R1.bySP, OTUtab.cytB.R1.bySP) + ggtitle("cytB (R1)"), 
                   
                   labels = LETTERS[1:2],
                   ncol = 2, nrow=1
                   #common.legend = T, legend = "bottom"
)


GRAPH <- ggarrange(ASVOTU.graph(ASVtab.12s.bySP, OTUtab.12s.bySP) + ggtitle("12S - Communauté simulée") ,
                                      
                   ASVOTU2.graph(ASVtab.12s.cor.bySP, OTUtab.12s.cor.bySP) + ggtitle("12S - Échantillons") ,
                   
                   ASVOTU.graph(ASVtab.cytB.R1.bySP, OTUtab.cytB.R1.bySP) + ggtitle("cytB (R1) - Communauté simulée"), 
                   
                   ASVOTU2.graph(ASVtab.cytB.R1.cor.bySP, OTUtab.cytB.R1.cor.bySP) + ggtitle("cytB (R1) - Échantillons"), 
                   
                   labels = LETTERS[1:4],
                   ncol = 2, nrow=2,
                   common.legend = T, legend = "right"
)

GRAPH

# Save graphes

ggsave(filename = file.path(get.value("result.FINAL"), "ASVvsOTU.pdf"),
       width = 8, height = 4,
       plot = GRAPH
)

ggsave(filename = file.path(get.value("result.FINAL"), "ASVvsOTU.png"),
       width = 8, height = 4,
       plot = GRAPH
)




# Mock community ----------------------------------------------------------

Mock.final <- Mock.dat %>% left_join(Mock.dat %>% group_by(Mix) %>% summarise(DNAtot = sum(DNA))) %>% 
                           mutate(DNAfinal = (DNA / DNAtot)*2) %>% 
                           select(Species, Mix, DNAfinal)

# Liste des espèces


pull.mock.sp <- function(tab){
tab %>% mutate(Level = str_count(Assign, ";")) %>% 
  filter(Level %in% c(5:6),
         str_detect(.$Assign, "Teleostei") == T) %>% 
  mutate(Assign = str_remove(Assign, "Root;Chordata;Teleostei;")) %>% 
  select(Assign, Level, names(.) %>% str_subset("Mix") ) %>%
  gather(names(.) %>% str_subset("Mix"), key = "sample", value = "Nread") %>% 
  mutate(Mix = sapply(str_split(sample, "_"), `[`, 2),
         Puit = sapply(str_split(sample, "p[:digit:]."), `[`, 2)) %>% 
  left_join(DataSeq %>% filter(SeqType %in% c("mix", "dup.mix")) %>% select(Puit, SeqType)) %>% 
  select(-Puit) %>% 
  filter(Nread >=1) %>% 
  pull(Assign) %>% unique()

}

pull.mock.sp(ASVtab.12s.wTAXO)
pull.mock.sp(OTUtab.12s.wTAXO)

pull.mock.sp(ASVtab.cytB.R1.wTAXO)
pull.mock.sp(OTUtab.cytB.R1.wTAXO)

pull.mock.sp(ASVtab.cytB.R2.wTAXO)
pull.mock.sp(OTUtab.cytB.R2.wTAXO)

# Graph

Mock.graph.data <- data.frame(Assign = character(),
                              Level = integer(), 
                              N = integer(), 
                              Mix = character(), 
                              SeqType = character(), 
                              Name.level = character(),
                              Data = character())

LS <- ls() %>% str_subset(".wTAXO") %>% str_remove(".cor") %>% unique()
LS[LS %>% str_detect("2x")==F]

for(x in LS[LS %>% str_detect("2x")==F]){

DATA <- get(x) %>% mutate(Level = str_count(Assign, ";")) %>% 
                     #filter(Level %in% c(5:6),
                    #        str_detect(.$Assign, "Teleostei") == T) %>% 
                     #mutate(Assign = str_remove(Assign, "Root;Chordata;Teleostei;")) %>% 
                     select(Assign, Level, names(.) %>% str_subset("Mix") ) %>%
                     gather(names(.) %>% str_subset("Mix"), key = "sample", value = "Nread") %>% 
                     group_by(Assign, Level, sample) %>% 
                     summarise(N = sum(Nread)) %>% 
                     mutate(Mix = sapply(str_split(sample, "_"), `[`, 2),
                            Puit = sapply(str_split(sample, "p[:digit:]."), `[`, 2) %>% str_remove("_R[:digit:]")) %>%
                     left_join(DataSeq %>% filter(SeqType %in% c("mix", "dup.mix")) %>% select(Puit, SeqType)) %>% 
                     select(-Puit) %>% mutate(Name.level = NA)

for(y in 1:nrow(DATA)){
  DATA$Name.level[y] <- sapply(str_split(DATA[y,"Assign"], ";"),`[`, pull(DATA[y,"Level"] + 1))
  
  }
  
DATA <- DATA %>% group_by() %>% select(Assign, Level, N, Mix, SeqType, Name.level) %>% 
  mutate(Data = x %>% str_remove(".wTAXO")) %>% as.data.frame()

Mock.graph.data <- rbind(Mock.graph.data, DATA)

}

# Version longue

Mock.graph.data %>% mutate(Method = str_sub(Data,1,3),
                           Locus = Data %>% str_remove(paste0(Method, "tab.")),
                           Taxo = ifelse(Level == 2, "Classe", 
                                               ifelse(Level == 3, "Ordre", 
                                               ifelse(Level == 4, "Famille",
                                               ifelse(Level == 5, "Genre", 
                                               ifelse(Level == 6, "Espèce", NA))))),
                           Taxo = factor(Taxo, levels = c("Classe", "Ordre", "Famille", "Genre", "Espèce" ))
                           ) %>%
  filter(str_detect(Assign, "Teleostei") == T,
         str_detect(Assign, "Gadidae") == F,
         SeqType == "mix",
         Locus %in% c("12s", "cytB.R1")) %>%
  mutate(Locus = ifelse(Locus == "12s", "12s", "cytB (R1)")) %>% 
  filter(N>=1) %>% 
  
  ggplot(aes(x = Mix, y = Name.level, fill = N)) + 
  geom_bin2d() + 
  scale_fill_distiller(palette = "Spectral", trans = "log10") +
  #scale_fill_gradient(low = "darkgray", high = "red", trans = "log") +
  #scale_y_discrete(limits=mixedsort(tab2$Assign)) + #, labels = NULL) +
  labs(title= NULL, x ="Communauté simulée", y = "Assignation") +
  guides(fill = guide_colourbar(title = "N lectures", title.hjust = 0)) + 
  theme_bw()+
  facet_grid(Taxo~Locus + Method, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0)) #+ coord_flip()

# Version courte

graph1 <- Mock.graph.data %>% mutate(Method = str_sub(Data,1,3),
                           Locus = Data %>% str_remove(paste0(Method, "tab.")),
                           Espece = ifelse(str_detect(Name.level, "Salvelinus"), "Salvelinus fontinalis", Name.level),
                           Presence = ifelse(Espece %in% Mock.dat$Species, "Espèce présente", "Espèce absente"),
                           Presence = factor(Presence, levels = c("Espèce présente", "Espèce absente")),
                           Espece = ifelse(str_detect(Name.level, "Salvelinus"), "Salvelinus sp.", Name.level),
                           NwNA = ifelse(N == 0, NA, N)) %>% 
                    filter(Method == "ASV",
                           Locus == "12s",
                           SeqType == "mix",
                           Mix %in% c("Mix1", "Mix2"),
                           str_detect(Assign, "Teleostei"),
                           str_detect(Espece, " "),
                           N >= 1) %>% 
  #bind_rows (expand.grid(Mix = c("Mix1", "Mix2"),
  #                       Espece = "Margariscus margarita")) %>% 
  left_join(REF %>% select(Espece_initial, Class, Order, Family, NomFR),
            by = c("Espece" = "Espece_initial"))  %>% #View()
  #mutate(Espece = ifelse(Espece == "Salvelinus fontinalis", "Salvelinus sp.", Espece),
  #       NomFR =ifelse(Espece == "Salvelinus sp.",  "Salvelinus sp.", NomFR)) %>% 
  filter(Presence == "Espèce présente") %>% 
  ggplot(aes(y = Mix, x = NomFR, fill = NwNA)) + 
  geom_bin2d(col = "darkgray") + 
  scale_fill_distiller(palette = "Spectral", trans = "log10",  na.value = "white", limits= c(1,NA)) +
  #scale_fill_gradient(low = "darkgray", high = "red", trans = "log") +
  #scale_y_discrete(limits=mixedsort(tab2$Assign)) + #, labels = NULL) +
  labs(title= NULL, x =NULL, y = "Communauté simulée") +
  guides(fill = guide_colourbar(title = "N séquences", title.hjust = 0)) + 
  theme_bw()+
  facet_grid(.~ Presence, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        #strip.text.x = element_text(angle = 90),
        strip.background = element_rect(fill="white")) #+ coord_flip()

graph1

ggsave(file.path(get.value("result.FINAL"),"Mock.short.png"), plot = graph1,
       width = 5, height = 3, units = "in")

# Abundance


Mock.abund.data <- Mock.graph.data %>% mutate(Method = str_sub(Data,1,3),
                                              Locus = Data %>% str_remove(paste0(Method, "tab."))) %>%  
                                       filter(Locus %in%  c("12s"),
                                              Method == "ASV") %>% # View()
                                       mutate(Espece = ifelse(str_detect(Name.level, "Salvelinus namaycush"), "Salvelinus fontinalis", Name.level)) %>% #View()
  
                                      filter(Espece %in% c("Salvelinus fontinalis", 
                                                              "Micropterus dolomieu")) %>% #View()
                                       left_join(Mock.final,
                                                 by = c("Espece" = "Species", "Mix" = "Mix")) %>% 

                                          mutate(Espece = ifelse(str_detect(Name.level, "Salvelinus namaycush"), "Salvelinus sp.", Name.level))

Mock.abund.data %>% group_by(Locus, Name.level, Method) %>% 
                    summarise(rho= cor.test(N, DNAfinal, method = "spearman")$estimate,
                              p.value= cor.test(N, DNAfinal, method = "spearman")$p.value,
                              N = length(N))

# Version courte

graph2 <- Mock.abund.data %>% 
  #mutate(Ncor = ifelse(N ==0, 0.8, N)) %>% 
  ggplot(aes(x = DNAfinal, y = N, col = Espece, shape = Espece))+
  geom_smooth(method = "lm", se = F , lty = "dashed", size = 1, fill = "gray", show.legend=FALSE)  +
  #geom_point(size = 2)+
  geom_jitter(width = 0.05, height = 0, size = 2)+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(limits = c(1, 1000), trans = "log10")+
   scale_color_discrete(name = "Légende",
                        breaks = c("Salvelinus sp.", "Micropterus dolomieu"),
                        labels = c("Salvelinus sp.", "Achigan à petite bouche"))+
   scale_shape_discrete(name = "Légende",
                        breaks = c("Salvelinus sp.", "Micropterus dolomieu"),
                        labels = c("Salvelinus sp.", "Achigan à petite bouche"))+
  labs(title= NULL, x ="Concentration d'ADN (ng/ul)", y = "N séquences") +                         
  #facet_grid(. ~ Species) +
  
 # guide_legend(title = "Espèce") +
  theme_classic() + theme(legend.position = c(0.75, 0.15),
                          legend.background = element_rect(colour="darkgray"))

graph2

ggsave(file.path(get.value("result.FINAL"),"Mock.abundancet.png"), plot = graph2,
       width = 5, height = 4, units = "in")


# N haplotypes ------------------------------------------------------------

haplo.mock <- function(tab){
# ASVtab.12s.byID.SP 
  
  tab %>% 
# Limit to MIX samples
filter( str_detect(Assign, "Teleostei") == T, 
        str_detect(Sample, "Mix"),
       Level %in% c(5:6),
      N > 0
       ) %>% 
  mutate(Mix = sapply(str_split(Sample, "_"), `[`, 2),
         Puit = sapply(str_split(Sample, "p[:digit:]."), `[`, 2) %>% str_remove("_R[:digit:]")) %>% 
  left_join(DataSeq %>% filter(SeqType %in% c("mix", "dup.mix")) %>% select(Puit, SeqType)) %>% 
  # Remove technical duplicates
  filter(SeqType == "mix") %>% 
  mutate(Assign2 = str_replace(Assign, ";Salvelinus", ";Salvelinus;Salvelinus fontinalis"),  
         Level = str_count(Assign2, ";")) %>%     
  filter(Level == 6) %>% 
  mutate(Species = sapply(str_split(Assign2, ";"), `[`, 7) ) %>% 
  
  group_by(Species, Mix) %>% 
  summarise(Nhaplo = length(unique(ID))) %>% 
  ggplot(aes(x = Species, y = Nhaplo, fill = Mix)) +
         geom_bar(stat = "identity", position = "dodge") +
         geom_hline(yintercept = 1, col = "black", linetype = 2) + 
         scale_y_continuous(limits = c(0,8))+
         theme_classic() +
         theme(axis.text.x = element_text(angle = 90, hjust = 1))

}



ggarrange(haplo.mock(ASVtab.12s.byID.SP) + labs(title = "12s - ASV"),
          haplo.mock(OTUtab.12s.byID.SP) + labs(title = "12s - OTU"),
          #haplo.mock(ASVtab.cytB.R1.byID.SP),
          #haplo.mock(OTUtab.cytB.R2.byID.SP),
          ncol = 2, nrow = 1, common.legend = T, legend = "right")

ls() %>% str_subset("byID.SP")

Summary.haplo <- data.frame(Assign = character(), Nhaplo = integer(), Data = character())


for(x in ls() %>% str_subset("byID.SP") %>% str_remove(".cor") %>% unique()){
  DATA <- get(x) %>% group_by(Assign, Level) %>% 
                     summarise(Nhaplo = length(unique(ID))) %>%  
                     mutate(Data = x %>% str_remove(".byID.SP")) %>% 
                     as.data.frame()
  
  
  Summary.haplo <- rbind(Summary.haplo, DATA)
  
}

Summary.haplo %>% filter(Assign != "Root", Level > 2) %>% 
  mutate(#Assign = str_remove(Assign, "Root;Chordata;Teleostei;"),
    Method = str_sub(Data, 1, 3),
    Locus = Data %>% str_remove(paste0(Method,"tab."))) %>% 
  arrange(Assign, Data) %>% 
  ggplot(aes(y= Nhaplo, x = Assign, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Locus ~ ., scale = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Trying haplo by lac

ASV.12s.SEQ <- data.frame(ID = ASVtab.12s$ID, SEQ = row.names(ASVtab.12s))

ASVtab.12s.wTAXO %>% left_join(ASV.12s.SEQ) %>%   filter(str_detect(Assign, "Anas"))

ASVtab.12s.wTAXO %>% left_join(ASV.12s.SEQ) %>% 
                         mutate(Assign = Assign %>% str_replace("Salvelinus", "Salvelinus;Salvelinus sp.")) %>% 
                         filter(str_detect(Assign, "Teleostei"),
                                str_detect(Assign, " "), 
                                str_detect(Assign, "Gadus morhua") == FALSE,
                                str_detect(Assign, "Salmo salar") == FALSE
                                ) %>%  
                         mutate(Espece = sapply(str_split(Assign, ";"), `[`, 7)) %>% 
                         gather(names(.) %>% str_subset("Sample"), key = "Sample", value = N) %>% 
                         mutate(Puit = sapply(str_split(Sample, "_"), `[`, 5) %>% str_remove("p")) %>% 
                         left_join(DataSeq %>% select(IbisID, SeqType, NomLac), by = c("Puit" = "IbisID")) %>% 
                         filter(SeqType == "sample") %>% 
                         group_by(NomLac, SEQ, Espece) %>% 
                         summarise(N = sum(N)) %>%
                         mutate(N = ifelse(N == 0, NA, N)) %>% 
                         #filter(N >=1) %>% 
                         #group_by(SEQ) %>% 
                         #summarise(Nlac = length(unique(NomLac)),
                         #         Nreads = sum(N))
                         #spread(SEQ, N) %>%
                         left_join(LacSample %>% select(NomLac, Rotenode, Affluent, Effluent, Bassin, SousBassin, Ordre, Volume),
                                   by = "NomLac") %>% 
                         mutate(NewBassin = ifelse(Bassin == "Isae", paste(Bassin,SousBassin,sep=":"), Bassin),
                                NewBassin = factor(NewBassin, levels = c("Isae:Ecarte", "Isae:Soumire", "Isae:Peche", "Isae:Francais", "Isae:Hamel", "Isae:Isae", "Bouchard", "Wapizagonke", "Aticagamac", "Cinq", "Cauche", "Theode", "St-Maurice", "Mattawin", "Isolé")),
                                NewNomLac = ifelse(is.na(Affluent), paste(NomLac, "*"), NomLac)) %>% 
                         arrange(NewBassin, Ordre) %>% #View()
                         ggplot(aes(x = NewNomLac, y = SEQ, fill = N))+
                         geom_bin2d(col = "darkgray") + 
                         scale_fill_distiller(palette = "Spectral", trans = "log10",  na.value = "white", limits= c(1,NA)) +
                         scale_y_discrete(labels = c(paste0("ASV",1:10))) +
                         labs(title = "Tadam" ) +
                         facet_grid(Espece ~ NewBassin, scale = "free", space = "free") +
                         theme_bw ()+
                         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, aes(size = 1/Volume/10000000)),
                               axis.ticks = element_blank(),
                               strip.text.y = element_text(angle = 0),
                               strip.text.x = element_text(angle = 90),
                               panel.spacing.x = unit(0, "lines"),
                               strip.placement = "outside",
                               strip.background = element_rect(colour = "black", fill = "white"),
                               axis.text.y = element_text(colour="black", hjust = 1),
                               legend.box.just = "left"
                               )
  



DATA <- ASVtab.12s.cor.wTAXO %>% left_join(ASV.12s.SEQ) %>% filter(Assign == "Root;Chordata;Teleostei;Salmoniformes;Salmonidae;Salvelinus") %>% pull(SEQ)
  mutate(Assign = Assign %>% str_replace("Salvelinus", "Salvelinus;Salvelinus sp.")) %>% 
  filter(str_detect(Assign, "Teleostei"),
         str_detect(Assign, " "), 
         str_detect(Assign, "Gadus morhua") == FALSE,
         str_detect(Assign, "Salmo salar") == FALSE
  ) %>%  
  mutate(Espece = sapply(str_split(Assign, ";"), `[`, 7)) %>% 
  gather(names(.) %>% str_subset("Sample"), key = "Sample", value = N) %>% 
  mutate(Puit = sapply(str_split(Sample, "_"), `[`, 5) %>% str_remove("p")) %>% 
  left_join(DataSeq %>% select(IbisID, SeqType, NomLac), by = c("Puit" = "IbisID")) %>% 
  filter(SeqType == "sample") %>% 
  group_by(NomLac, SEQ, Espece) %>% 
  summarise(N = sum(N)) %>%
 # mutate(N = ifelse(N == 0, NA, N)) %>% 
  filter(Espece == "Micropterus dolomieu",
         N >= 1) 
  
library(Biostrings)

DNA <- DNAStringSet(DATA$SEQ)
names(DNA) <- DATA$NomLac

DNA2 <- as.DNAbin(DNA)

library(ape)

DNA.dist <- dist.dna(DNA2)
DNA.dist

library(pegas)

HAPLO <- pegas::haplotype(DNA2)

HAPLO <- sort(HAPLO, what = "label")
str(HAPLO)

(net <- pegas::haploNet(HAPLO))

ind.hap<-with(
  stack(setNames(attr(HAPLO, "index"), rownames(HAPLO))),
  table(hap=ind, individuals=names(DNA)[values])
)



par(mar=c(3,3,2,10))
plot(net, size=c(1,4), scale.ratio=1, pie=ind.hap, show.mutation = 2, labels = FALSE,  main = "Salvelinus network")
legend(0, -3, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=4, xjust = 0.5, yjust = 1)

attr(net, "freq")

par(mar=c(10,3,5,3))
plot(net, size=1 , scale.ratio=1, pie=ind.hap, show.mutation = 1, labels = TRUE,  main = "Ameiurus")
legend(0, -2.5, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=4, xjust = 0.5, yjust = 1)


par(mar=c(10,3,5,3))
plot(net, size=1 , scale.ratio=1, pie=ind.hap, show.mutation = 1, labels = TRUE,  main = "Semotilus")
legend(0, -2.5, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=4, xjust = 0.5, yjust = 1)

par(mar=c(10,3,5,3))
plot(net, size=1 , scale.ratio=1, pie=ind.hap, show.mutation = 1, labels = TRUE,  main = "Perca")
legend(0, -2.5, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=4, xjust = 0.5, yjust = 1)

par(mar=c(10,3,5,3))
plot(net, size=1 , scale.ratio=1, pie=ind.hap, show.mutation = 1, labels = TRUE,  main = "Micropterus")
legend(0, -2.5, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=4, xjust = 0.5, yjust = 1)




# Samples -----------------------------------------------------------------

Sample.data <- data.frame(Assign = character(), 
                          Sample  = character(),
                          N = integer(),
                          Level = integer(),
                          NameAssign = character(),
                          ID = character(),
                          Cat = character(),
                          Nsite = integer(),
                          NomLac = character(),
                          Data = character())

#for(x in ls() %>% str_subset(".cor.2x.wTAXO")){ # For stringent correction
for(x in ls() %>% str_subset(".cor.wTAXO")){
  print(x)
  
DATA <- get(x) %>% gather(names(.) %>% str_subset("Sample"), key = Sample, value = N) %>% 
        group_by(Assign, ID, Sample) %>% 
        summarise(N = sum(N)) %>% 
        
        mutate(Level = str_count(Assign, ";"),
                                 NameAssign = NA,
          Lac = sapply(str_split(Sample, "_"), `[`, 3),
                                                 Cat = sapply(str_split(Sample, "_"), `[`, 4),
                                                Puit = sapply(str_split(Sample, "_"), `[`, 5) %>% str_remove("p")) %>% 
          left_join(Blast99) %>% 
          left_join(DataSeq %>% select(IbisID, Nsite, CatSite, SeqType, NomLac, CorrFiltre, Volume), by = c("Puit" = "IbisID")) %>% 
          filter(SeqType == "sample") %>% 
        mutate(Data = x %>% str_remove(".cor.wTAXO") %>% str_remove(".cor.2x.wTAXO"),
               Cat = CatSite)   
  
  
for(y in 1:nrow(DATA)){
  DATA$NameAssign[y] <- sapply(str_split(DATA[y,"Assign"], ";"),`[`, pull(DATA[y,"Level"] + 1))
  
  }

DATA <- DATA %>% group_by() %>% select(-c(Lac, Puit, SeqType)) %>% 
  as.data.frame()  
 
Sample.data <- rbind(Sample.data, DATA)
 
}  

# When it doesn't work, try this
#for(y in 1:nrow(Sample.data)){
#  Sample.data$NameAssign[y] <- sapply(str_split(Sample.data[y,"Assign"], ";"),`[`, Sample.data[y,"Level"] + 1)
#  
#}

unique(Sample.data$NameAssign)

# Add info from blast 99 results (hereafter NameAssign.99)

Sample.data <- Sample.data %>% mutate(Method = str_sub(Data,1,3),
                                      Locus = Data %>% str_remove(paste0(Method, "tab.")),
                                      #Integrate ALL the 99 blast results
                                      SP = ifelse(Locus == "12s", SP, NA),
                                      #SP = ifelse( str_detect(SP,"/")==F, SP, NA ),
                                      #Identic.SP = ifelse(NameAssign == SP, "Oui", "Non"),
                                      # Integrate only luxilus cornutus  
                                      NameAssign.99 = ifelse(SP == "Luxilus cornutus", SP, NameAssign), 
                                      NameAssign.99 = ifelse(is.na(NameAssign.99), NameAssign, NameAssign.99),
                                      # 
                                      Locus = Data %>% str_remove(paste0(Method, "tab.")),
                                      Taxo = ifelse(Level == 0, "Root",
                                             ifelse(Level == 1, "Phylum",
                                             ifelse(Level == 2, "Classe", 
                                             ifelse(Level == 3, "Ordre", 
                                             ifelse(Level == 4, "Famille",
                                             ifelse(Level == 5, "Genre", 
                                             ifelse(Level == 6, "Espece", NA))))))),
                                      Taxo = factor(Taxo, levels = c("Root","Phylum","Classe", "Ordre", "Famille", "Genre", "Espece" ))
                                      ) %>% 
                              left_join(LacInv1 %>% filter(Espece == "Salvelinus fontinalis") %>% select(NomLac, Location)) %>% 
                              mutate(Location = factor(Location, levels = c("Avant-pays", "Arriere-pays")))

# Basis stats

Sample.data %>% filter(Data == "ASVtab.12s") %>% 
  group_by(Location, NomLac, Cat) %>% 
  summarise(Ntot = length(unique(Sample))) %>% 
  spread(Cat, Ntot, fill = 0) %>% 
  mutate(N = PEL + RIV) %>% View()

Sample.data %>% filter(Data == "ASVtab.12s") %>% 
  group_by(Taxo) %>% 
  summarise(N = sum(N),
            Perc = N / sum(.$N))


# Mean sample by lake
Sample.data %>% filter(Data == "ASVtab.12s") %>% 
  group_by(Location, NomLac, Cat) %>% 
  summarise(Ntot = length(unique(Sample))) %>% 
  spread(Cat, Ntot, fill = 0) %>% 
  mutate(N = PEL + RIV) %>%
  group_by() %>% 
  summarise(Mean = mean(N))


# Potential false positive 
Sample.data %>% filter(Data == "ASVtab.12s") %>% 
  group_by(NameAssign.99) %>%
  summarise(Nech = length(unique(Sample[N>=1])),
            Nlac = length(unique(NomLac[N>=1]))) %>% View()

# Graph of what was done

Sample.data %>% filter(#str_detect(Assign, "Teleostei") == F,
                   #str_detect(Assign, "Gadidae") == F,
                   #str_detect(Assign, "Sebaste") == F,
                   Locus %in% c("12s"),
                   Method %in% c("ASV")) %>%
  group_by(NomLac, Locus, Method, NameAssign, NameAssign.99, Taxo, Location) %>% 
  summarise(Nmed = median(N),
            Nmax = max(N)) %>% #View()
  
  #mutate(Locus = ifelse(Locus == "12s", "12s", "cytB (R1)")) %>% 
  filter(Nmax>=0) %>% # View()
  
  ggplot(aes(x = NomLac, y = NameAssign, fill = Nmax)) + 
  geom_bin2d(col = "gray") + 
  scale_fill_distiller(palette = "Spectral", trans = "log10",  na.value = "White") +
  #scale_fill_gradient(low = "darkgray", high = "red", trans = "log") +
  #scale_y_discrete(limits=mixedsort(tab2$Assign)) + #, labels = NULL) +
  labs(title= NULL, x ="Lac", y = "Assignation") +
  guides(fill = guide_colourbar(title = "N lectures", title.hjust = 0)) + 
  theme_bw()+
  facet_grid(Taxo~ Location, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0)) #+ coord_flip()


# Graph of what was done check first samples

Sample.data %>% filter(#str_detect(Assign, "Teleostei") == F,
  #str_detect(Assign, "Gadidae") == F,
  #str_detect(Assign, "Sebaste") == F,
  Locus %in% c("12s"),
  Method %in% c("ASV"),
  Location == "Avant-pays") %>%
  group_by(NomLac, Sample, Locus, Method, NameAssign, NameAssign.99, Taxo, Location) %>% 
  summarise(N = sum(N, na.rm=T),
            Nmax = max(N)) %>% #View()
  
  #mutate(Locus = ifelse(Locus == "12s", "12s", "cytB (R1)")) %>% 
  filter(Nmax>=0) %>% # View()
  
  ggplot(aes(x = Sample, y = NameAssign, fill = N)) + 
  geom_bin2d(col = "gray") + 
  scale_fill_distiller(palette = "Spectral", trans = "log10",  na.value = "White") +
  #scale_fill_gradient(low = "darkgray", high = "red", trans = "log") +
  #scale_y_discrete(limits=mixedsort(tab2$Assign)) + #, labels = NULL) +
  labs(title= NULL, x ="Lac", y = "Assignation") +
  guides(fill = guide_colourbar(title = "N lectures", title.hjust = 0)) + 
  theme_bw()+
  facet_grid(Taxo~ NomLac, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0)) #+ coord_flip()


# Vérifier ce qu'il se passe dans les lacs Giron, Parker et Pêche

Sample.data %>% filter(#str_detect(Assign, "Teleostei") == F,
  #str_detect(Assign, "Gadidae") == F,
  #str_detect(Assign, "Sebaste") == F,
  Locus %in% c("12s"),
  Method %in% c("ASV"),
  NomLac %in% c("Giron", "Peche (a la)", "Parker")) %>% 
  group_by(NomLac, Locus, Method, NameAssign, NameAssign.99, Taxo, Location, Sample) %>% 
  summarise(N = sum(N)) %>% #View()

  ggplot(aes(x = Sample, y = NameAssign, fill = N)) + 
  geom_bin2d(col = "gray") + 
  scale_fill_distiller(palette = "Spectral", trans = "log10",  na.value = "White") +
  #scale_fill_gradient(low = "darkgray", high = "red", trans = "log") +
  #scale_y_discrete(limits=mixedsort(tab2$Assign)) + #, labels = NULL) +
  labs(title= NULL, x ="Lac", y = "Assignation") +
  guides(fill = guide_colourbar(title = "N lectures", title.hjust = 0)) + 
  theme_bw()+
  facet_grid(Taxo~ NomLac, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0)) #+ coord_flip()


# Graphe pour vois les N min a instaurer
Sample.data %>% mutate(NameAssign.99 = ifelse(NameAssign.99 == "Salvelinus", "Salvelinus sp.", 
                                              ifelse(NameAssign.99 == "Salvelinus namaycush/Salvelinus fontinalis/Salvelinus alpinus", "Salvelinus sp.",NameAssign.99))) %>% 
  filter(str_detect(Assign, "Teleostei"),
         str_detect(NameAssign.99, " "), 
         str_detect(NameAssign.99, "Gadus morhua") == FALSE,
         #str_detect(NameAssign.99, "Salmo salar") == FALSE,
         str_detect(NameAssign.99, "Sebaste") == FALSE,#,
         # Enlever quelques échantillons particulièrement contaminés
         str_detect(Sample, "Sample_Edo_R_p1.G1") == FALSE, #7
         str_detect(Sample, "p1.D3") == FALSE, #19
         str_detect(Sample, "p1.G3") == FALSE, #22
         str_detect(Sample, "p1.B4") == FALSE, # 25
         str_detect(Sample, "p1.B2") == FALSE, #10
         str_detect(Sample, "p1.E2") == FALSE #13
         #Locus %in% c("12s"),
         #Method %in% c("ASV")
  )  %>% #View()
  #Some species with problems
  group_by(NomLac, Locus, Method, Sample, NameAssign.99, Volume, CorrFiltre) %>% 
  summarise(N = sum(N, na.rm = T)) %>%
  filter(N>0, 
         Method == "ASV",
         Locus == "12s",
         NameAssign.99 != "Salmo salar") %>% 
  #filter(NameAssign.99 == "Ameiurus nebulosus", NomLac == "Solitaire") %>% View()
  ggplot(aes(x=N)) + 
  geom_histogram() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~NameAssign.99, scales = "free")





# Proportion des échantilllons avec l'une ou l'autre des espèces

Sample.graph <- Sample.data %>% mutate(NameAssign.99 = ifelse(NameAssign.99 == "Salvelinus", "Salvelinus sp.", 
                                                              ifelse(NameAssign.99 == "Salvelinus namaycush/Salvelinus fontinalis/Salvelinus alpinus", "Salvelinus sp.",NameAssign.99))) %>% 
                                filter(str_detect(Assign, "Teleostei"),
                                       str_detect(NameAssign.99, " "), 
                                       str_detect(NameAssign.99, "Gadus morhua") == FALSE,
                                       #str_detect(NameAssign.99, "Salmo salar") == FALSE, # Potentiellement la ouananiche, on garde
                                       str_detect(NameAssign.99, "Sebaste") == FALSE,#,
                                       # Enlever quelques échantillons particulièrement contaminés
                                       str_detect(Sample, "Sample_Edo_R_p1.G1") == FALSE, #7
                                       str_detect(Sample, "p1.D3") == FALSE, #19
                                       str_detect(Sample, "p1.G3") == FALSE, #22
                                       str_detect(Sample, "p1.B4") == FALSE, # 25
                                       str_detect(Sample, "p1.B2") == FALSE, #10
                                       str_detect(Sample, "p1.E2") == FALSE #13
                                       #Locus %in% c("12s"),
                                       #Method %in% c("ASV")
                                       )  %>% #View()
                                #Some species with problems
                                group_by(NomLac, Locus, Method, Sample, NameAssign.99, Volume, CorrFiltre) %>% 
                                summarise(N = sum(N, na.rm = T)) %>%  
                                #Complete for all species
                                complete(NameAssign.99, NomLac, Method, Locus) %>% 
                                mutate(N = ifelse(is.na(N), 0, N),
                                       # Enlever les N == 1
                                       N = ifelse(N == 1 & Locus == "12s", 0, N),
                                       # Etre plus sévère pour la barbotte
                                       N = ifelse(N < 250 & NameAssign.99 == "Ameiurus nebulosus", 0, N),
                                       Nlog = ifelse(N == 0 , NA, 
                                                    #log2(N)),
                                                     ifelse(N == 1, 0.1, log2(N))),
                                       Nlog.cor = ifelse(is.na(Nlog), NA, Nlog / CorrFiltre / Volume * 1000)
                                       ) %>% #View()
                                # Get one line by sp / lake
                                group_by(NameAssign.99, NomLac, Method, Locus) %>% 
                                summarise(
                                          Ncor = mean(Nlog.cor, na.rm = T),
                                          #Nsample = length(unique(Sample)), # existe plus bas
                                          Nmax = max(N),
                                          NsampleTot = length(N),
                                          NsamplePre = length(N[N>=1]),
                                          PropSample = NsamplePre/NsampleTot,
                                          N = mean(N)) %>% #View() 
                                mutate(PropSample2 = ifelse(PropSample == 0, NA, PropSample),
                                       PresenceADNe = ifelse(NsamplePre>=1,1,NA),
                                       NsamplePre2 = ifelse(NsamplePre >=2,"2 et plus",NA)) %>%
                                #filter(N >=1 & NameAssign.99 == "Ameiurus nebulosus")
                                # Remove grouping factor
                                group_by() %>% 
                                # Comparison trad vs ADNe
                                #mutate(PresenceADNe = ifelse(N > 0, 1, 0)) %>% #View()
                                full_join(LacInv1 %>% select(NomLac, Location, Espece, Presence),
                                          by = c("NomLac" = "NomLac", "NameAssign.99" = "Espece")) %>% 
                                mutate(PresenceADNe = ifelse(is.na (PresenceADNe), 0, PresenceADNe),
                                       DiffInv = ifelse(PresenceADNe == Presence, 
                                                        ifelse(PresenceADNe == 0 , "Absent trad et ADNe", "Present trad et ADNe"),
                                                        ifelse(PresenceADNe == 0, "Present trad seul", "Present ADNe seul"))) %>% #View()
                                # Add species french name
                                left_join(REF %>% select(Espece_initial, Class, Order, Family, NomFR),
                                          by = c("NameAssign.99" = "Espece_initial")) %>% 
                                # Add lac info
                                left_join(LacSample %>% select(NomLac, Rotenode, Affluent, Effluent, Bassin, SousBassin, Ordre, Volume),
                                          by = "NomLac") %>% 
                                mutate(NewBassin = ifelse(Bassin == "Isae", paste(Bassin,SousBassin,sep=":"), Bassin),
                                       NewBassin = factor(NewBassin, levels = c("Isae:Ecarte", "Isae:Soumire", "Isae:Peche", "Isae:Francais", "Isae:Hamel", "Isae:Isae", "Bouchard", "Isae-Ouest", "Wapizagonke", "Aticagamac", "Cinq", "Brier", "Theode", "St-Maurice", "Mattawin", "Isolé")),
                                       NewNomLac = ifelse(is.na(Affluent), paste(NomLac, "*"), NomLac),
                                       Ordre = factor(Ordre),
                                       NewNomLac = factor(NewNomLac, levels=unique(NewNomLac[order(NewBassin, Ordre)]))
                                ) %>% 
                                arrange(NewBassin, Ordre)  

# Keep only
Sample.graph.red <- bind_rows(Sample.graph %>% filter(Method %in% c("ASV", NA),
                                                  Locus %in% c("12s", NA),
                                                  NomFR != "Omble de fontaine"),
                          Sample.graph %>% filter(Method == "ASV",
                                                  Locus == "cytB.R1",
                                                  NomFR == "Omble de fontaine")# %>% View()
                          ) %>%  
                mutate(NomFR = ifelse(NomFR %in% c("Omble de fontaine", "Méné à nageoires rouges"), paste(NomFR, "*"), NomFR),
                       NomFR = ifelse(NomFR == "Saumon atlantique", "Ouananiche", NomFR),
                       NomFR = factor(NomFR, levels = rev(c("Omble de fontaine *", "Omble chevalier", "Touladi", "Salvelinus sp.", "Ouananiche", 
                                                          
                                                            "Barbotte brune", 
                                                            "Meunier noir",
                                                            "Perchaude", "Doré jaune", "Fouille-roche zébré",
                                                            "Achigan à petite bouche", "Crapet de roche", "Crapet-soleil",
                                                            "Épinoche à cinq épines", "Épinoche à neuf épines",
                                                            "Mulet à cornes", "Mulet perlé", "Tête-de-boule", "Ventre rouge du nord", "Méné à nageoires rouges *", "Méné jaune", "Museau noir", "Naseux des rapides", "Ouitouche",
                                                            
                                                            "Chabot à tête plate",
                                                            "Fondule barré",
                                                            "Éperlan arc-en-ciel",
                                                            "Grand brochet")))) %>% 
            left_join(Sample.graph %>% filter(!is.na(NsampleTot)) %>% group_by(NewNomLac) %>% summarise(NsampleTot2 = unique(NsampleTot))) %>% 
            mutate(NewNomLac2 = paste0(NewNomLac, " (",NsampleTot2, ")" ))


# Other species

Sample.other.graph <- Sample.data %>%  filter(Method == "ASV",
                                        Locus == "12s",
                                        NameAssign %in% c("Lithobates pipiens",
                                                          "Lithobates catesbeianus",
                                                          "Eurycea bislineata",
                                                          #"Meleagris gallopavo",
                                                          #"Anas",
                                                           "Anas platyrhynchos",
                                                          "Ondatra zibethicus",
                                                          "Artiodactyla",
                                                          "Homo sapiens"),
                                       # Enlever quelques échantillons particulièrement contaminés
                                       str_detect(Sample, "Sample_Edo_R_p1.G1") == FALSE, #7
                                       str_detect(Sample, "p1.D3") == FALSE, #19
                                       str_detect(Sample, "p1.G3") == FALSE, #22
                                       str_detect(Sample, "p1.B4") == FALSE, # 25
                                       str_detect(Sample, "p1.B2") == FALSE, #10
                                       str_detect(Sample, "p1.E2") == FALSE #13
                                       #Locus %in% c("12s"),
                                       #Method %in% c("ASV")
                                )  %>% #View()
                                #Some species with problems
                                group_by(NomLac, Locus, Method, Sample, NameAssign, Volume, CorrFiltre) %>% 
                                summarise(N = sum(N, na.rm = T)) %>%  
                                #Complete for all species
                                complete(NameAssign, NomLac, Method, Locus) %>% 
                                mutate(N = ifelse(is.na(N), 0, N),
                                       # Enlever les N == 1
                                       N = ifelse(N == 1 & Locus == "12s", 0, N),
                                       # Etre plus sévère pour la barbotte
                                       #N = ifelse(N < 250 & NameAssign.99 == "Ameiurus nebulosus", 0, N),
                                       Nlog = ifelse(N == 0 , NA, 
                                                     #log2(N)),
                                                     ifelse(N == 1, 0.1, log2(N))),
                                       Nlog.cor = ifelse(is.na(Nlog), NA, Nlog / CorrFiltre / Volume * 1000)
                                ) %>% #View()
                                # Get one line by sp / lake
                                group_by(NameAssign, NomLac, Method, Locus) %>% 
                                summarise(
                                  Ncor = mean(Nlog.cor, na.rm = T),
                                  #Nsample = length(unique(Sample)), # existe plus bas
                                  Nmax = max(N),
                                  NsampleTot = length(N),
                                  NsamplePre = length(N[N>=1]),
                                  PropSample = NsamplePre/NsampleTot,
                                  N = mean(N)) %>% #View() 
                                mutate(PropSample2 = ifelse(PropSample == 0, NA, PropSample),
                                       PresenceADNe = ifelse(NsamplePre>=1,1,NA),
                                       NsamplePre2 = ifelse(NsamplePre >=2,"2 et plus",NA)) %>%
                                # Remove grouping factor
                                group_by() %>% 
                                # Comparison trad vs ADNe
                                #mutate(PresenceADNe = ifelse(N > 0, 1, 0)) %>% #View()
                                #full_join(LacInv1 %>% select(NomLac, Location, Espece, Presence),
                                #          by = c("NomLac" = "NomLac", "NameAssign.99" = "Espece")) %>% 
                                #mutate(PresenceADNe = ifelse(is.na (PresenceADNe), 0, PresenceADNe),
                                 #      DiffInv = ifelse(PresenceADNe == Presence, 
                                 #                       ifelse(PresenceADNe == 0 , "Absent trad et ADNe", "Present trad et ADNe"),
                                  #                      ifelse(PresenceADNe == 0, "Present trad seul", "Present ADNe seul"))) %>% #View()
                                # Add species french name
                                #left_join(REF %>% select(Espece_initial, Class, Order, Family, NomFR),
                                #          by = c("NameAssign.99" = "Espece_initial")) %>% 
                                # Add lac info
                                left_join(LacSample %>% select(NomLac, Rotenode, Location, Affluent, Effluent, Bassin, SousBassin, Ordre, Volume),
                                          by = "NomLac") %>% 
                                mutate(NewBassin = ifelse(Bassin == "Isae", paste(Bassin,SousBassin,sep=":"), Bassin),
                                       NewBassin = factor(NewBassin, levels = c("Isae:Ecarte", "Isae:Soumire", "Isae:Peche", "Isae:Francais", "Isae:Hamel", "Isae:Isae", "Bouchard", "Isae-Ouest", "Wapizagonke", "Aticagamac", "Cinq", "Brier", "Theode", "St-Maurice", "Mattawin", "Isolé")),
                                       NewNomLac = ifelse(is.na(Affluent), paste(NomLac, "*"), NomLac),
                                       Ordre = factor(Ordre),
                                       NewNomLac = factor(NewNomLac, levels=unique(NewNomLac[order(NewBassin, Ordre)]))
                                ) %>% 
                                arrange(NewBassin, Ordre)  




# N moyen d'espèce par lac (compte les 0) 

Sample.graph.red %>% filter(N > 0) %>% group_by(NomLac) %>% 
                     summarise(Nsp = length(unique(NameAssign.99))) %>%
                     pull(Nsp) %>% sum()/44

Sample.graph.red %>% filter(N > 0) %>% group_by(NomLac, Location) %>% 
  summarise(Nsp = length(unique(NameAssign.99))) %>%
  ggplot(aes(x = Nsp, fill=Location)) + geom_histogram(position="dodge")

Sample.graph.red %>% filter(N > 0) %>% group_by(NomLac, Location) %>% 
  summarise(Nsp = length(unique(NameAssign.99))) %>%
  group_by(Location, Nsp) %>% 
  summarise(N = n()) %>% 


Sample.graph %>% group_by(NameAssign.99, NsamplePre) %>% 
  summarise(N = n()) %>% 
  group_by() %>% 
  complete(NameAssign.99, NsamplePre, fill = list(0)) %>% #View()
  filter(NsamplePre >=1) %>% 
  ggplot(aes(x=NsamplePre, y=N, fill=factor(NsamplePre))) +
  geom_bar(width = 1, stat = "identity", position = "dodge") +
  #coord_polar("y", start=0) +
  facet_wrap(~NameAssign.99)
 

Sample.graph %>% group_by(NsamplePre) %>% 
  summarise(N = n())


# N group

Sample.graph.red %>% filter(Location == "Avant-pays") %>% 
  group_by(DiffInv) %>% 
  summarise(N = n(),
            Perc = N / 702)

 
SP.presentes <- Sample.graph.red %>% filter(N >= 1) %>% 
                                   pull(NomFR) %>% unique() %>% as.character()




graph3 <- Sample.graph.red  %>% filter(Location == "Avant-pays",
                                       NomFR %in% SP.presentes) %>% 
  ggplot(aes(x = NewNomLac, y = NomFR, fill = factor(PresenceADNe))) + 
  geom_bin2d(col = "gray", na.rm = FALSE) + 
  scale_fill_manual(values = "limegreen", limits = "1") +
  geom_point(aes(shape = factor(Presence)), col = "gray20") +
  
  scale_shape_manual(values = 19, limits = "1", guide = "none") +
  
  labs(title= NULL, x =NULL, y = NULL) +
  guides(fill = FALSE) + 
  theme_bw()+ 
  facet_grid(. ~ NewBassin, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.background = element_rect(fill="white")
        ) 

graph3




graph3.1 <- Sample.graph.red  %>% filter(Location == "Avant-pays",
                                       NomFR %in% SP.presentes) %>% 
  ggplot(aes(x = NewNomLac, y = NomFR, fill = factor(PresenceADNe))) + 
  geom_bin2d(col = "gray", na.rm = FALSE) + 
  scale_fill_manual(values = "white", limits = "1") +
#  geom_point(aes(shape = factor(Presence)), col = "gray20") +
  
#  scale_shape_manual(values = 19, limits = "1", guide = "none") +
  
  labs(title= NULL, x =NULL, y = NULL) +
  guides(fill = FALSE) + 
  theme_bw()+
  facet_grid(. ~ NewBassin, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.background = element_rect(fill="white")
  ) 

graph3.1


graph3.2a <- Sample.graph.red  %>% filter(Location == "Avant-pays",
                                       NomFR %in% SP.presentes) %>% 
                                      mutate(PresenceADNe = ifelse(N >= 10,1,0)) %>% 
  ggplot(aes(x = NewNomLac, y = NomFR, fill = factor(PresenceADNe))) + 
  geom_bin2d(col = "gray", na.rm = FALSE) + 
  scale_fill_manual(values = "limegreen", limits = "1") +
  geom_point(aes(shape = factor(Presence)), col = "gray20") +
  
  scale_shape_manual(values = 19, limits = "1", guide = "none") +
  labs(title= NULL, x =NULL, y = NULL) +
  guides(fill = FALSE) + 
  theme_bw()+
  facet_grid(. ~ NewBassin, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.background = element_rect(fill="white")
  ) 

graph3.2a

graph3.2b <- Sample.graph.red  %>% filter(Location == "Avant-pays",
                                          NomFR %in% SP.presentes) %>% 
  mutate(PresenceADNe = ifelse(N >= 100,1,0)) %>% 
  ggplot(aes(x = NewNomLac, y = NomFR, fill = factor(PresenceADNe))) + 
  geom_bin2d(col = "gray", na.rm = FALSE) + 
  scale_fill_manual(values = "limegreen", limits = "1") +
  geom_point(aes(shape = factor(Presence)), col = "gray20") +
  
  scale_shape_manual(values = 19, limits = "1", guide = "none") +
  labs(title= NULL, x =NULL, y = NULL) +
  guides(fill = FALSE) + 
  theme_bw()+
  facet_grid(. ~ NewBassin, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.background = element_rect(fill="white")
  ) 

graph3.2b


graph3.3 <- Sample.graph.red  %>% filter(Location == "Avant-pays",
                                         NomFR %nin% c("Méné jaune", "Crapet-soleil", "Doré jaune", "Grand brochet")) %>% #View()
  ggplot(aes(x = NewNomLac, y = NomFR, fill = factor(PresenceADNe))) + 
  geom_bin2d(col = "gray", na.rm = FALSE) + 
  scale_fill_manual(values = "limegreen", limits = "1") +
  geom_point(aes(shape = factor(Presence)), col = "gray20") +

  scale_shape_manual(values = 19, limits = "1", guide = "none") +
  
  labs(title= "Avant-Pays", x =NULL, y = NULL) +
  guides(fill = FALSE) + 
  theme_bw()+ 
  facet_grid(. ~ NewBassin, scale = "free", space = "free") +

theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.ticks.y = element_blank(),
      strip.text.x = element_text(angle = 90),
      strip.background = element_rect(fill="white")
) 

graph3.3

  


ggsave(filename = file.path(get.value("result.FINAL"), "PresenteAbsence.AvP.ASV.12S.png"),
       width = 6.5, height = 5,
       plot = graph3
)

ggsave(filename = file.path(get.value("result.FINAL"), "PresenteAbsence.AvP.ASV.12S.vide.png"),
       width = 6.5, height = 5,
       plot = graph3.1
)

ggsave(filename = file.path(get.value("result.FINAL"), "PresenteAbsence.AvP.ASV.12S.min10.png"),
       width = 6.5, height = 5,
       plot = graph3.2a
)

ggsave(filename = file.path(get.value("result.FINAL"), "PresenteAbsence.AvP.ASV.12S.min100.png"),
       width = 6.5, height = 5,
       plot = graph3.2b
)

ggsave(filename = file.path(get.value("result.FINAL"), "PresenteAbsence.COmpTrad.png"),
       width = 10, height = 8,
       plot = graph3.3
)


graph4 <- Sample.graph.red  %>% filter(Location != "Avant-pays",
                                       NomFR %in% SP.presentes) %>% 
  ggplot(aes(x = NewNomLac, y = NomFR, fill = factor(PresenceADNe))) + 
  geom_bin2d(col = "gray", na.rm = FALSE) + 
  scale_fill_manual(values = "limegreen", limits = "1") +
  labs(title= NULL, x =NULL, y = NULL) +
  guides(fill = FALSE) + 
  theme_bw()+
  facet_grid(. ~ NewBassin, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust= 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.background = element_rect(fill="white")
  ) 

graph4


graph4.3 <- Sample.graph.red  %>% filter(Location != "Avant-pays",
                                         NomFR %nin% c("Méné jaune", "Crapet-soleil", "Doré jaune", "Grand brochet")) %>% #View()
  ggplot(aes(x = NewNomLac, y = NomFR, fill = factor(PresenceADNe))) + 
  geom_bin2d(col = "gray", na.rm = FALSE) + 
  scale_fill_manual(values = "limegreen", limits = "1") +
 # geom_point(aes(shape = factor(Presence)), col = "gray20") +
  
  scale_shape_manual(values = 19, limits = "1", guide = "none") +
  
  labs(title= "Arrière-Pays", x =NULL, y = NULL) +
  guides(fill = FALSE) + 
  theme_bw()+ 
  facet_grid(. ~ NewBassin, scale = "free", space = "free") +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.background = element_rect(fill="white")
  ) 

graph4.3

ggsave(filename = file.path(get.value("result.FINAL"), "PresenteAbsence.CompTrad2.png"),
       width = 10, height = 8,
       plot = graph4.3
)


ggsave(filename = file.path(get.value("result.FINAL"), "PresenteAbsence.ArrP.ASV.12S.png"),
       width = 6.5, height = 5,
       plot = graph4
)

ggarrange(graph3 + theme(plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm")) +
                  ggtitle("Avant-pays"),
          graph4 + theme(plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm")) +
            theme(axis.text.y=element_blank())+
            ggtitle("Arrière-pays"),
          #labels = c("A. Avant-pays", "B. Arrière-pays"),
          #vjust = 1, hjust = 1,
          widths = c(1.7, 1),
          ncol=2, align = "hv")

# Then add info on quality

graph5 <- Sample.graph.red  %>% filter(Location == "Avant-pays",
                                       NomFR  %nin% c("Méné jaune", "Crapet-soleil", "Doré jaune", "Grand brochet")) %>% 
                                mutate(Ncor = ifelse(Ncor == 0, NA, Ncor)) %>% #View() #pull(Ncor) %>% max()
  ggplot(aes(x = NewNomLac2, y = NomFR, fill = Ncor, shape = factor(NsamplePre))) + 
  geom_bin2d(col = "gray", na.rm = FALSE) + 
  #scale_fill_distiller(palette = "Spectral", direction = -1, na.value="white", limits = c(2,110)) +
  scale_fill_gradientn(colours =c("skyblue", "navyblue"), na.value="white", limits = c(0,110), breaks= c(1,25,50, 75,100)) + 
  
    #scale_fill_discrete(na.value = "white", guide = "none") +
  
  geom_point(aes(shape = factor(Presence)), col = "black") +
  
  scale_shape_manual(values = 19, limits = "1", guide = "none") +
  
  
  #geom_point(size = 2,  stroke = 1, col = "gray20") +
  #scale_shape_manual(values = c(1), name = NULL, limits = "2 et plus", labels = c("Présent dans plus\nd'un échantillon")) +
  #scale_shape_manual(values = c(49,50,51,52,53,54,55,56,57,58), name = NULL, limits = c("1","2","3", "4", "5", "6", "7", "8","9","10"), guide = "none") +
  
  #scale_fill_gradient(low = "darkgray", high = "red", trans = "log") +
  #scale_y_discrete(limits=mixedsort(tab2$Assign)) + #, labels = NULL) +
  labs(title= NULL, x =NULL, y = NULL, fill = "Indice\nd'abondance") +
  
 # guides(fill = FALSE) + 
  theme_bw()+
  facet_grid(. ~ NewBassin, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.background = element_rect(fill="white"),
        legend.position = "top",
        legend.text = element_text(size=8),
        legend.title = element_text(size=9)#,
        #legend.key.size = unit(0.5, "inch")
        
  ) 


graph5 

graph5.1 <- Sample.graph.red  %>% filter(Location == "Avant-pays",
                                         NomFR  %nin% c("Méné jaune", "Crapet-soleil", "Doré jaune", "Grand brochet")
                                         #NomFR %in% SP.presentes,
                                       #NsampleTot > 1 
                                       
                                       ) %>%  
 # mutate(Ncor = ifelse(Ncor == 0, NA, Ncor),
 #        Ncor = ifelse(PropSample < 0.5, NA, Ncor ),
 #        Ncor = ifelse(NsamplePre ==1, NA, Ncor)) %>% #View() #pull(Ncor) %>% max()
  
  mutate(Ncor = ifelse(Ncor == 0, NA, Ncor),
         #   Ncor = ifelse(PropSample < 0.5, NA, Ncor ),
         Ncor = ifelse(NsamplePre < 2, NA, Ncor)) %>% #View() #pull(Ncor) %>% max()
  
  
  #filter(!is.na(Ncor)) %>% 
  ggplot(aes(x = NewNomLac2, y = NomFR, fill = Ncor, shape = factor(Presence))) + 
  geom_bin2d(col = "gray", na.rm =FALSE) + 
  #scale_fill_distiller(palette = "Reds", direction = 1, na.value="white", limits = c(0,110)) +
  scale_fill_gradient(low = "skyblue", high = "navyblue", na.value="white", limits = c(0,110), breaks= c(1,25,50, 75,100)) + 
  
  #scale_fill_discrete(na.value = "white", guide = "none") +
  geom_point(aes(shape = factor(Presence)), col = "gray20") +
  
  scale_shape_manual(values = 19, limits = "1", guide = "none") +
  #scale_y_discrete(limits=mixedsort(tab2$Assign)) + #, labels = NULL) +
  labs(title= NULL, x =NULL, y = NULL, fill = "Indice\nd'abondance") +
  
  #guides(fill = FALSE) + 
  theme_bw()+
  facet_grid(. ~ NewBassin, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.background = element_rect(fill="white"),
        legend.position = "top",
        legend.text = element_text(size=8),
        legend.title = element_text(size=9)
  ) 

graph5.1 


graph5.2 <- Sample.other.graph  %>% filter(Location == "Avant-pays") %>% 
  mutate(Ncor = ifelse(Ncor == 0, NA, Ncor)) %>% #View() #pull(Ncor) %>% max()
  ggplot(aes(x = NewNomLac2, y = NameAssign, fill = Ncor, shape = factor(NsamplePre))) + 
  geom_bin2d(col = "gray", na.rm = FALSE) + 
  #scale_fill_distiller(palette = "Spectral", direction = -1, na.value="white", limits = c(2,110)) +
  scale_fill_gradientn(colours =c("skyblue", "navyblue"), na.value="white", limits = c(0,110), breaks= c(1,25,50, 75,100)) + 
  
  scale_y_discrete(limits = c("Lithobates pipiens",
                                "Lithobates catesbeianus",
                                "Eurycea bislineata",
                                #"Meleagris gallopavo",
                                #"Anas",
                                "Anas platyrhynchos",
                                "Ondatra zibethicus",
                                "Artiodactyla",
                                "Homo sapiens"), labels= c("Grenouille léopard", "Ouaouaron", "Salamandre à deux lignes", "Anatidés", "Rat musqué", "Artiodactilé", "Humain")) +
  
  labs(title= NULL, x =NULL, y = NULL, fill = "Indice\nd'abondance") +
  
  # guides(fill = FALSE) + 
  theme_bw()+
  facet_grid(. ~ NewBassin, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.background = element_rect(fill="white"),
        legend.position = "top",
        legend.text = element_text(size=8),
        legend.title = element_text(size=9)#,
        #legend.key.size = unit(0.5, "inch")
        
  ) 


graph5.2 

graph6 <- Sample.graph.red  %>% filter(Location == "Arriere-pays",
                                       NomFR  %nin% c("Méné jaune", "Crapet-soleil", "Doré jaune", "Grand brochet", "Ouananiche", "Fouille-roche zébré", "Épinoche à neuf épines", "Méné à nageoires rouges *", "Museau noir", "Ouitouche", "Chabot à tête plate", "Fondule barré", "Éperlan arc-en-ciel")
                                       #NomFR %in% SP.presentes,
                                       ) %>% 
  mutate(Ncor = ifelse(Ncor == 0, NA, Ncor)) %>%# pull(Ncor) %>% max()
  ggplot(aes(x = NewNomLac2, y = NomFR, fill = Ncor, shape = factor(NsamplePre))) + 
  geom_bin2d(col = "gray", na.rm = FALSE) + 
  #scale_fill_distiller(palette = "Reds", direction = 1, na.value="white", limits = c(0,110)) +
  scale_fill_gradient(low = "skyblue", high = "navyblue", na.value="white", limits = c(0,110), breaks= c(1,25,50, 75,100)) + 
  
  #scale_fill_discrete(na.value = "white", guide = "none") +
  
  
  geom_point(aes(shape = factor(Presence)), col = "gray20") +
  
  scale_shape_manual(values = 19, limits = "1", guide = "none") +
  
  
  #geom_point(size = 2,  stroke = 1, col = "gray20") +
  
  #scale_shape_manual(values = c(1), name = NULL, limits = "2 et plus", labels = c("Présent dans plus\nd'un échantillon")) +
  
  #scale_shape_manual(values = c(49,50,51,52,53,54,55,56,57,58), name = NULL, limits = c("1","2","3", "4", "5", "6", "7", "8","9","10"), guide = "none") +
  
  #scale_fill_gradient(low = "darkgray", high = "red", trans = "log") +
  #scale_y_discrete(limits=mixedsort(tab2$Assign)) + #, labels = NULL) +
  labs(title= NULL, x =NULL, y = NULL, fill = "Indice\nd'abondance") +
#  guides(fill = FALSE) + 
  theme_bw()+
  facet_grid(. ~ NewBassin, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.background = element_rect(fill="white"),
        legend.position = "top",
        legend.text = element_text(size=8),
        legend.title = element_text(size=9)
  ) 


graph6

graph6.0 <- Sample.graph.red  %>% filter(Location == "Arriere-pays",
                                         NomFR  %nin% c("Méné jaune", "Crapet-soleil", "Doré jaune", "Grand brochet") #, "Ouananiche", "Fouille-roche zébré", "Épinoche à neuf épines", "Méné à nageoires rouges *", "Museau noir", "Ouitouche", "Chabot à tête plate", "Fondule barré", "Éperlan arc-en-ciel")
                                       #NomFR %in% SP.presentes,
) %>% 
  mutate(Ncor = ifelse(Ncor == 0, NA, Ncor)) %>%# pull(Ncor) %>% max()
  ggplot(aes(x = NewNomLac2, y = NomFR, fill = Ncor, shape = factor(NsamplePre))) + 
  geom_bin2d(col = "gray", na.rm = FALSE) + 
  #scale_fill_distiller(palette = "Reds", direction = 1, na.value="white", limits = c(0,110)) +
  scale_fill_gradient(low = "skyblue", high = "navyblue", na.value="white", limits = c(0,110), breaks= c(1,25,50, 75,100)) + 
  
  #scale_fill_discrete(na.value = "white", guide = "none") +
  
  
  geom_point(aes(shape = factor(Presence)), col = "gray20") +
  
  scale_shape_manual(values = 19, limits = "1", guide = "none") +
  
  
  #geom_point(size = 2,  stroke = 1, col = "gray20") +
  
  #scale_shape_manual(values = c(1), name = NULL, limits = "2 et plus", labels = c("Présent dans plus\nd'un échantillon")) +
  
  #scale_shape_manual(values = c(49,50,51,52,53,54,55,56,57,58), name = NULL, limits = c("1","2","3", "4", "5", "6", "7", "8","9","10"), guide = "none") +
  
  #scale_fill_gradient(low = "darkgray", high = "red", trans = "log") +
  #scale_y_discrete(limits=mixedsort(tab2$Assign)) + #, labels = NULL) +
  labs(title= NULL, x =NULL, y = NULL, fill = "Indice\nd'abondance") +
  #  guides(fill = FALSE) + 
  theme_bw()+
  facet_grid(. ~ NewBassin, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.background = element_rect(fill="white"),
        legend.position = "top",
        legend.text = element_text(size=8),
        legend.title = element_text(size=9)
  ) 


graph6.0


graph6.1 <- Sample.graph.red  %>% filter(Location != "Avant-pays",
                                         NomFR  %nin% c("Méné jaune", "Crapet-soleil", "Doré jaune", "Grand brochet", "Ouananiche", "Fouille-roche zébré", "Épinoche à neuf épines", "Méné à nageoires rouges *", "Museau noir", "Ouitouche", "Chabot à tête plate", "Fondule barré", "Éperlan arc-en-ciel")
                                         #NomFR %in% SP.presentes,
                                         #NsampleTot > 1
                                         ) %>% 
#  mutate(Ncor = ifelse(Ncor == 0, NA, Ncor),
#         Ncor = ifelse(PropSample < 0.5, NA, Ncor ),
#         Ncor = ifelse(NsamplePre ==1, NA, Ncor)) %>% #View() #pull(Ncor) %>% max()
  
  mutate(Ncor = ifelse(Ncor == 0, NA, Ncor),
        #   Ncor = ifelse(PropSample < 0.5, NA, Ncor ),
         Ncor = ifelse(NsamplePre < 2, NA, Ncor)) %>% #View() #pull(Ncor) %>% max()
  
  
  ggplot(aes(x = NewNomLac2, y = NomFR, fill = Ncor, shape = factor(Presence))) + 
  geom_bin2d(col = "gray", na.rm = FALSE) + 
  #scale_fill_distiller(palette = "Reds", direction = 1, na.value="white", limits = c(0,110)) +
  scale_fill_gradient(low = "skyblue", high = "navyblue", na.value="white", limits = c(0,110), breaks= c(1,25,50, 75,100)) + 
  #scale_fill_discrete(na.value = "white", guide = "none") +
  geom_point(aes(shape = factor(Presence)), col = "gray20") +
  
  scale_shape_manual(values = 19, limits = "1", guide = "none") +
  labs(title= NULL, x =NULL, y = NULL, fill = "Indice\nd'abondance") +
  
  #guides(fill = FALSE) + 
  theme_bw()+
  facet_grid(. ~ NewBassin, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.background = element_rect(fill="white"),
        legend.position = "top",
        legend.text = element_text(size=8),
        legend.title = element_text(size=9)
  ) 

graph6.1 

graph6.2 <- Sample.graph.red  %>% filter(Location != "Avant-pays",
                                         NomFR  %nin% c("Méné jaune", "Crapet-soleil", "Doré jaune", "Grand brochet") #, "Ouananiche", "Fouille-roche zébré", "Épinoche à neuf épines", "Méné à nageoires rouges *", "Museau noir", "Ouitouche", "Chabot à tête plate", "Fondule barré", "Éperlan arc-en-ciel")
                                         #NomFR %in% SP.presentes,
                                         #NsampleTot > 1
) %>% 
  #  mutate(Ncor = ifelse(Ncor == 0, NA, Ncor),
  #         Ncor = ifelse(PropSample < 0.5, NA, Ncor ),
  #         Ncor = ifelse(NsamplePre ==1, NA, Ncor)) %>% #View() #pull(Ncor) %>% max()
  
  mutate(Ncor = ifelse(Ncor == 0, NA, Ncor),
         #   Ncor = ifelse(PropSample < 0.5, NA, Ncor ),
         Ncor = ifelse(NsamplePre < 2, NA, Ncor)) %>% #View() #pull(Ncor) %>% max()
  ggplot(aes(x = NewNomLac2, y = NomFR, fill = Ncor, shape = factor(Presence))) + 
  geom_bin2d(col = "gray", na.rm = FALSE) + 
  #scale_fill_distiller(palette = "Reds", direction = 1, na.value="white", limits = c(0,110)) +
  scale_fill_gradient(low = "skyblue", high = "navyblue", na.value="white", limits = c(0,110), breaks= c(1,25,50, 75,100)) + 
  #scale_fill_discrete(na.value = "white", guide = "none") +
  geom_point(aes(shape = factor(Presence)), col = "gray20") +
  
  scale_shape_manual(values = 19, limits = "1", guide = "none") +
  labs(title= NULL, x =NULL, y = NULL, fill = "Indice\nd'abondance") +
  
  #guides(fill = FALSE) + 
  theme_bw()+
  facet_grid(. ~ NewBassin, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.background = element_rect(fill="white"),
        legend.position = "top",
        legend.text = element_text(size=8),
        legend.title = element_text(size=9)
  ) 

graph6.2 

graph6.3 <- Sample.other.graph  %>% filter(Location != "Avant-pays") %>% 
  mutate(Ncor = ifelse(Ncor == 0, NA, Ncor)) %>% #View() #pull(Ncor) %>% max()
  ggplot(aes(x = NewNomLac2, y = NameAssign, fill = Ncor, shape = factor(NsamplePre))) + 
  geom_bin2d(col = "gray", na.rm = FALSE) + 
  #scale_fill_distiller(palette = "Spectral", direction = -1, na.value="white", limits = c(2,110)) +
  scale_fill_gradientn(colours =c("skyblue", "navyblue"), na.value="white", limits = c(0,110), breaks= c(1,25,50, 75,100)) + 
  
  scale_y_discrete(limits = c("Lithobates pipiens",
                              "Lithobates catesbeianus",
                              "Eurycea bislineata",
                              #"Meleagris gallopavo",
                              #"Anas",
                              "Anas platyrhynchos",
                              "Ondatra zibethicus",
                              "Artiodactyla",
                              "Homo sapiens"), labels= c("Grenouille léopard", "Ouaouaron", "Salamandre à deux lignes", "Anatidés", "Rat musqué", "Artiodactilé", "Humain")) +
  
  labs(title= NULL, x =NULL, y = NULL, fill = "Indice\nd'abondance") +
  
  # guides(fill = FALSE) + 
  theme_bw()+
  facet_grid(. ~ NewBassin, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        strip.background = element_rect(fill="white"),
        legend.position = "top",
        legend.text = element_text(size=8),
        legend.title = element_text(size=9)#,
        #legend.key.size = unit(0.5, "inch")
        
  ) 


graph6.3 

ggsave(filename = file.path(get.value("result.FINAL"), "Abondance.AvP.ASV.12S.png"),
       width = 6.5, height = 5.5,
       plot = graph5
)



ggsave(filename = file.path(get.value("result.FINAL"), "Abondance.AvP.Seuil50.ASV.12S.png"),
       width = 6.5, height = 5.5,
       plot = graph5.1
)

ggsave(filename = file.path(get.value("result.FINAL"), "Abondance.ArrP.ASV.12S.png"),
       width = 6, height = 5,
       plot = graph6
)

ggsave(filename = file.path(get.value("result.FINAL"), "Abondance.ArrP.Seuil50.ASV.12S.png"),
       width = 6, height = 5,
       plot = graph6.1
)

graph56 <- ggarrange(graph5.1 + theme(plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm")) +
            ggtitle("Avant-pays"),
          graph6.2 + theme(plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm")) +
            theme(axis.text.y=element_blank())+
            ggtitle("Arrière-pays"),
          #labels = c("A. Avant-pays", "B. Arrière-pays"),
          #vjust = 1, hjust = 1,
          widths = c(1.7, 1),
          ncol=2, align = "hv", common.legend = T, legend = "top") 

graph56

graph56.0 <- ggarrange(graph5 + theme(plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm")) +
                         ggtitle("Avant-pays - Sans minimum"),
                       graph6.0 + theme(plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm")) +
                         theme(axis.text.y=element_blank())+
                         ggtitle("Arrière-pays - Sans minimum"),
                       graph5.1 + theme(plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm")) +
                       ggtitle("Avant-pays - Minimum 2 détections"),
                     graph6.2 + theme(plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm")) +
                       theme(axis.text.y=element_blank())+
                       ggtitle("Arrière-pays - Minimum 2 détections"),
                     #labels = c("A. Avant-pays", "B. Arrière-pays"),
                     #vjust = 1, hjust = 1,
                     widths = c(1.7, 1),
                     ncol=2, nrow=2, align = "hv", common.legend = T, legend = "top") 

graph56.0

graph56.1 <- ggarrange(graph5.2 + theme(plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm")) +
                       ggtitle("Avant-pays"),
                     graph6.3 + theme(plot.margin=unit(c(0.5,0.5, 0.5,0.5),"cm")) +
                       theme(axis.text.y=element_blank())+
                       ggtitle("Arrière-pays"),
                     #labels = c("A. Avant-pays", "B. Arrière-pays"),
                     #vjust = 1, hjust = 1,
                     widths = c(1.7, 1),
                     ncol=2, align = "hv", common.legend = T, legend = "top") 

graph56.1

ggsave(filename = file.path(get.value("result.FINAL"), "Abondance.Total.Seuil50.ASV.12S.png"),
       width = 10, height = 6,
       plot = graph56
)

ggsave(filename = file.path(get.value("result.FINAL"), "Abondance.Total.Seuil50.4panels.ASV.12S.png"),
       width = 10, height = 10,
       plot = graph56.0
)


ggsave(filename = file.path(get.value("result.FINAL"), "Abondance.Total.Others.ASV.12S.png"),
       width = 10, height = 5,
       plot = graph56.1
)


# Liste des espèces échantillonnés dans les lacs du Parc
SP.PARC <- LacInv1 %>% filter(NomLac %in% Sample.graph.red$NomLac %>% unique(),
                              Presence == 1) %>% 
                        pull(Espece) %>% unique()

SP.PARC.Salvelinus <- SP.PARC[str_detect(SP.PARC, "Salvelinus")]

# Calculer le N cases
Sample.graph.red  %>% filter(Location == "Avant-pays",
                             NameAssign.99 %in% SP.PARC) %>%
                      group_by(DiffInv) %>% 
                      summarise(N = n(),
                                Perc = N/546)  

# Calculer le N d'ajout par espèce

Sample.graph.red  %>% filter(Location == "Avant-pays",
                             NameAssign.99 %in% SP.PARC,
                             DiffInv == "Present ADNe seul") %>%
  group_by(NameAssign.99) %>% 
  summarise(N = n(),
            Perc = N/26) %>% 
  arrange(desc(N))


Sample.graph.red  %>% filter(Location == "Avant-pays",
                             NameAssign.99 %in% SP.PARC) %>% #View()
  ggplot(aes(x = NewNomLac, y = NomFR, shape = DiffInv)) + 
  geom_bin2d(aes(fill = Ncor), col = "gray", size = 0.5) + 
  geom_point(size = 3, fill = "yellow",  stroke = 1, col = "gray30") + 
  
  scale_fill_distiller(palette = "Reds",
                       direction = 1,
                       #trans = "log10",  
                       na.value = "White", limits = c(1,110)) +
  
  scale_shape_manual(values = c(3,6), name = "Comparaison avec\nl'inventaire traditionnel", limits = c("Present ADNe seul", "Present trad seul"), labels = c("Ajout", "Manquant")) +
  #scale_colour_manual(values = c("gray", "red"), name = "Traitement à la roténode", limits = c("0", "1"), labels = c("Non", "Oui")) +
  
  labs(title= NULL, x =NULL, y = NULL) +
  #guides(fill = guide_colourbar(title = "N moyen\nlectures", title.hjust = 0)) +
  guides(fill = guide_colourbar(title = "Indice d'abondance\nd'ADNe", title.hjust = 0)) +
  theme_bw ()+
  facet_grid(. ~ NewBassin, 
             scale = "free", space = "free", 
             labeller = label_value,
             switch = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, aes(size = 1/Volume/10000000)),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90),
        panel.spacing = unit(0, "lines"),
        strip.placement = "outside",
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.text.y = element_text(colour="black", hjust = 1),
        legend.box.just = "left") + 
  ggtitle("Avant-Pays")

graph7.1 <- Sample.graph.red  %>% filter(Location == "Avant-pays",
                                         #NameAssign.99 %in% SP.PARC,             
                                         NameAssign.99 %in% c(SP.PARC.Salvelinus)
                                         ) %>% #View()
  ggplot(aes(x = NewNomLac, y = NomFR, shape = DiffInv)) + 
  geom_bin2d(aes(fill = Ncor), col = "gray", size = 0.5) + 
#  geom_point(size = 3, fill = "yellow",  stroke = 1, col = "gray30") + 
  
  scale_fill_distiller(palette = "Reds",
                       direction = 1,
                       #trans = "log10",  
                       na.value = "White", limits = c(0,110)) +
  
  geom_point(aes(shape = factor(Presence)), col = "gray20") +
  
  scale_shape_manual(values = 19, limits = "1", guide = "none") +
  
  
 # scale_shape_manual(values = c(3,6), name = "Comparaison avec\nl'inventaire traditionnel", limits = c("Present ADNe seul", "Present trad seul"), labels = c("Ajout", "Manquant")) +
  #scale_colour_manual(values = c("gray", "red"), name = "Traitement à la roténode", limits = c("0", "1"), labels = c("Non", "Oui")) +
  
  labs(title= NULL, x =NULL, y = NULL) +
  #guides(fill = guide_colourbar(title = "N moyen\nlectures", title.hjust = 0)) +
  guides(fill = guide_colourbar(title = "Indice d'abondance\nd'ADNe", title.hjust = 0)) +
  theme_bw ()+
  facet_grid(. ~ NewBassin, 
             scale = "free", space = "free", 
             labeller = label_value,
             switch = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, aes(size = 1/Volume/10000000)),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90),
        panel.spacing = unit(0, "lines"),
        strip.placement = "outside",
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.text.y = element_text(colour="black", hjust = 1),
        legend.position = "none") 

graph7.1


graph7.1.1 <- Sample.graph.red  %>% filter(Location == "Avant-pays",
                                         #NameAssign.99 %in% SP.PARC,             
                                         NameAssign.99 %in% c(SP.PARC.Salvelinus)) %>% 
                                           mutate(Ncor = 1) %>% 

  ggplot(aes(x = NewNomLac, y = NomFR, shape = DiffInv)) + 
  geom_bin2d(aes(fill = Ncor), col = "gray", size = 0.5) + 
  # geom_point(size = 3, fill = "yellow",  stroke = 1, col = "gray30") + 
  
  scale_fill_distiller(palette = "Reds",
                       direction = 1,
                       #trans = "log10",  
                       na.value = "White", limits = c(2,110)) +
  
  scale_shape_manual(values = c(3,6), name = "Comparaison avec\nl'inventaire traditionnel", limits = c("Present ADNe seul", "Present trad seul"), labels = c("Ajout", "Manquant")) +
  #scale_colour_manual(values = c("gray", "red"), name = "Traitement à la roténode", limits = c("0", "1"), labels = c("Non", "Oui")) +
  
  labs(title= NULL, x =NULL, y = NULL) +
  #guides(fill = guide_colourbar(title = "N moyen\nlectures", title.hjust = 0)) +
  guides(fill = guide_colourbar(title = "Indice d'abondance\nd'ADNe", title.hjust = 0)) +
  theme_bw ()+
  facet_grid(. ~ NewBassin, 
             scale = "free", space = "free", 
             labeller = label_value,
             switch = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, aes(size = 1/Volume/10000000)),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90),
        panel.spacing = unit(0, "lines"),
        strip.placement = "outside",
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.text.y = element_text(colour="black", hjust = 1),
        legend.position = "none") 

graph7.1.1

graph7.1.2 <- Sample.graph.red  %>% filter(Location == "Avant-pays",
                                         #NameAssign.99 %in% SP.PARC,             
                                         NameAssign.99 %in% c(SP.PARC.Salvelinus)
) %>% #View()
  ggplot(aes(x = NewNomLac, y = NomFR, shape = DiffInv)) + 
  geom_bin2d(aes(fill = Ncor), col = "gray", size = 0.5) + 
#  geom_point(size = 3, fill = "yellow",  stroke = 1, col = "gray30") + 
  
  scale_fill_distiller(palette = "Reds",
                       direction = 1,
                       #trans = "log10",  
                       na.value = "White", limits = c(0,110)) +
  
  scale_shape_manual(values = c(3,6), name = "Comparaison avec\nl'inventaire traditionnel", limits = c("Present ADNe seul", "Present trad seul"), labels = c("Ajout", "Manquant")) +
  #scale_colour_manual(values = c("gray", "red"), name = "Traitement à la roténode", limits = c("0", "1"), labels = c("Non", "Oui")) +
  
  labs(title= NULL, x =NULL, y = NULL) +
  #guides(fill = guide_colourbar(title = "N moyen\nlectures", title.hjust = 0)) +
  guides(fill = guide_colourbar(title = "Indice d'abondance\nd'ADNe", title.hjust = 0)) +
  theme_bw ()+
  facet_grid(. ~ NewBassin, 
             scale = "free", space = "free", 
             labeller = label_value,
             switch = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, aes(size = 1/Volume/10000000)),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90),
        panel.spacing = unit(0, "lines"),
        strip.placement = "outside",
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.text.y = element_text(colour="black", hjust = 1),
        legend.position = "none") 

graph7.1.2

graph7.2 <- Sample.graph.red  %>% filter(Location == "Avant-pays",
                             NameAssign.99 %in% SP.PARC,             
                             NameAssign.99 %nin% c(SP.PARC.Salvelinus),
                             NomFR %nin% c("Fouille-roche zébré",
                                           "Épinoche à neuf épines",
                                           "Museau noir",
                                           "Ouitouche",
                                           "Chabot à tête plate",
                                           "Fondule barré",
                                           "Éperlan arc-en-ciel")) %>% #View()
  ggplot(aes(x = NewNomLac, y = NomFR, shape = DiffInv)) + 
  geom_bin2d(aes(fill = Ncor), col = "gray", size = 0.5) + 
  #geom_point(size = 3, fill = "yellow",  stroke = 1, col = "gray30") + 
  
  scale_fill_distiller(palette = "Reds",
                       direction = 1,
                       #trans = "log10",  
                       na.value = "White", limits = c(0,110)) +
  geom_point(aes(shape = factor(Presence)), col = "gray20") +
  
  scale_shape_manual(values = 19, limits = "1", guide = "none") +
  
  #scale_shape_manual(values = c(3,6), name = "Comparaison avec\nl'inventaire traditionnel", limits = c("Present ADNe seul", "Present trad seul"), labels = c("Ajout", "Manquant")) +
  #scale_colour_manual(values = c("gray", "red"), name = "Traitement à la roténode", limits = c("0", "1"), labels = c("Non", "Oui")) +
  
  labs(title= NULL, x =NULL, y = NULL) +
  #guides(fill = guide_colourbar(title = "N moyen\nlectures", title.hjust = 0)) +
  guides(fill = guide_colourbar(title = "Indice d'abondance\nd'ADNe", title.hjust = 0)) +
  theme_bw ()+
  facet_grid(. ~ NewBassin, 
             scale = "free", space = "free", 
             labeller = label_value,
             switch = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, aes(size = 1/Volume/10000000)),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90),
        panel.spacing = unit(0, "lines"),
        strip.placement = "outside",
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.text.y = element_text(colour="black", hjust = 1),
        legend.position = "none") 

graph7.2


ggsave(filename = file.path(get.value("result.FINAL"), "CompSalvelinus.AvP.ASV.12S.png"),
       width = 6.5, height = 3.5,
       plot = graph7.1
)

ggsave(filename = file.path(get.value("result.FINAL"), "CompSalvelinus.AvP.ASV.12S.vide.png"),
       width = 6.5, height = 3.5,
       plot = graph7.1.1
)

ggsave(filename = file.path(get.value("result.FINAL"), "CompSalvelinus.AvP.ASV.12S.min.png"),
       width = 6.5, height = 3.5,
       plot = graph7.1.2
)

ggsave(filename = file.path(get.value("result.FINAL"), "CompAutres.AvP.ASV.12S.png"),
       width = 6.5, height = 5,
       plot = graph7.2
)


Sample.graph.red  %>% 
  group_by(NewNomLac, Location, NewBassin) %>% 
  summarise(Volume = mean(Volume/1e-9),
            Nsample = max(NsampleTot, na.rm=T)) %>% 
  filter(Location == "Avant-pays") %>% #View()
  ggplot(aes(x = NewNomLac, y = "Volume", size = Volume, fill = factor(Nsample))) + 
  geom_bin2d(fill = "White", col = "gray", size = 0.5) + 
  geom_point(colour = "black", shape = 21) +
  scale_y_discrete(labels = "Caractéristiques\ndes lacs")+
  scale_fill_brewer(name = "N échantillons", palette = "Blues", na.value = "White", limits = c(2:7), breaks = c(2:7)) +
  scale_size(name = "Volume") +
  #scale_shape_discrete(limits = c("Present ADNe seul", "Present trad seul")) +
  labs(title= NULL, x =NULL, y = NULL) +
  theme_bw ()+
  facet_grid(. ~ NewBassin, 
             scale = "free", space = "free", 
             labeller = label_value,
             switch = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, aes(size = 1/Volume/10000000)),
        axis.text.y = element_text(colour="black", hjust = 1, face = "bold"),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90),
        panel.spacing = unit(0, "lines"),
        strip.placement = "outside",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.background = element_rect(colour = "black", fill = "white"),
        panel.grid = element_blank(),
        legend.box.just = "top",
        legend.box = "horizontal"
        #legend.position =  "none"
        #panel.grid.major = element_blanck()),  panel.grid.minor = element_line(colour = NULL)
  ) #+ coord_flip()


Sample.graph.red  %>% filter(Location == "Avant-pays",
                             NameAssign.99 %in% SP.PARC) %>% #View()
  ggplot(aes(x = NewNomLac, y = NomFR, shape = DiffInv)) + 
  geom_bin2d(aes(fill = DiffInv), col = "gray", size = 0.5) + 
  #geom_point(size = 3, fill = "yellow",  stroke = 1, col = "gray30") + 
  scale_fill_manual(values = c("green1", "green1", "yellow1", "yellow1"), 
                    limits = c("Absent trad et ADNe", "Present trad et ADNe", "Present trad seul", "Present ADNe seul"), 
                    labels = c("Similaire", "Similaire", "Different", "Different"),
                    guide = "none")+

  labs(title= NULL, x =NULL, y = NULL) +
  #guides(fill = guide_colourbar(title = "N moyen\nlectures", title.hjust = 0)) +
  #guides(fill = guide_colourbar(title = "Indice d'abondance\nd'ADNe", title.hjust = 0)) +
  theme_bw ()+
  facet_grid(. ~ NewBassin, 
             scale = "free", space = "free", 
             labeller = label_value,
             switch = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, aes(size = 1/Volume/10000000)),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90),
        panel.spacing = unit(0, "lines"),
        strip.placement = "outside",
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.text.y = element_text(colour="black", hjust = 1),
        legend.box.just = "left",
        legend.position =  "none") + 
  ggtitle("Avant-Pays")





# Modelisation ------------------------------------------------------------


Sample.mod.data.int <- Sample.data %>% select(Assign, Sample, N, NameAssign.99, Method, Locus) %>% 
  
                                   mutate(NameAssign.99 = ifelse(NameAssign.99 == "Salvelinus", "Salvelinus sp.", 
                                                                                  ifelse(NameAssign.99 == "Salvelinus namaycush/Salvelinus fontinalis/Salvelinus alpinus", "Salvelinus sp.",NameAssign.99))) %>% 
                                    filter(str_detect(Assign, "Teleostei"),
                                           str_detect(NameAssign.99, " "), 
                                           str_detect(NameAssign.99, "Gadus morhua") == FALSE,
                                           #str_detect(NameAssign.99, "Salmo salar") == FALSE,
                                           str_detect(NameAssign.99, "Sebaste") == FALSE         
                                           #Locus %in% c("12s"),
                                           #Method %in% c("ASV")
                                    )  %>% 
                                   # Correct the number of observed reads
                                   
                                   select(NameAssign.99, Sample, Method, Locus, N) %>% #View()
                                   group_by(NameAssign.99, Sample, Method, Locus) %>% 
                                  summarise(N = sum(N)) %>%    
                                   complete(NameAssign.99, Sample) %>% 
                                   group_by() %>% 
                                    mutate(N = ifelse(is.na(N), 0, N),
                                          PresenceADNe = ifelse(N > 0, 1, 0),
                                          #Locus = sapply(str_split(Sample, "_"), `[`, 1),
                                          Puit = sapply(str_split(Sample, "_"), `[`, 5) %>% str_remove("p")) %>% #View()
                                   left_join(DataSeq %>% select(IbisID, Nsite, CatSite, SeqType, NomLac, CorrFiltre, Volume), by = c("Puit" = "IbisID")) %>% 
                                   left_join(LacInv1 %>% select(NomLac, Location, Espece, Presence),
                                             by = c("NomLac" = "NomLac", "NameAssign.99" = "Espece")) %>% 
                                   mutate(PresenceADNe = ifelse(is.na (PresenceADNe), 0, PresenceADNe),
                                          DiffInv = ifelse(PresenceADNe == Presence, 
                                                           ifelse(PresenceADNe == 0 , "Absent trad et ADNe", "Present trad et ADNe"),
                                                           ifelse(PresenceADNe == 0, "Present trad seul", "Present ADNe seul"))) %>% 
                                   left_join(REF %>% select(Espece_initial, Class, Order, Family, NomFR),
                                             by = c("NameAssign.99" = "Espece_initial")) %>% 
                                   left_join(LacSample %>% select(NomLac, Rotenode, Affluent, Effluent, Bassin, SousBassin, Ordre, Volume, Surface),
                                             by = "NomLac") %>% 
                                   mutate(NewBassin = ifelse(Bassin == "Isae", paste(Bassin,SousBassin,sep=":"), Bassin),
                                          NewBassin = factor(NewBassin, levels = c("Isae:Ecarte", "Isae:Soumire", "Isae:Peche", "Isae:Francais", "Isae:Hamel", "Isae:Isae", "Bouchard", "Wapizagonke", "Aticagamac", "Cinq", "Cauche", "Theode", "St-Maurice", "Mattawin", "Isolé")),
                                          NewNomLac = ifelse(is.na(Affluent), paste(NomLac, "*"), NomLac)) %>% 
                                   arrange(NewBassin, Ordre) %>% 
                                   mutate(NewNomLac = factor(NewNomLac, levels = unique(NewNomLac)), 
                                          Family = factor(Family, levels = c("Salmonidae", "Cyprinidae", "Catostomidae", "Gasterosteidae", "Cottidae", "Fundulidae", "Osmeridae", "Centrarchidae", "Ictaluridae", "Percidae", "Esocidae"))
                                   )
 
Sample.mod.data.int$Locus %>% unique()
   
Sample.mod.data  <- bind_rows(Sample.mod.data.int %>% filter(Method %in% c("ASV", NA),
                                                      Locus %in% c("12s", NA),
                                                      NomFR != "Omble de fontaine"),
                              Sample.mod.data.int %>% filter(Method == "ASV",
                                                      Locus == "cytB.R1",
                                                      NomFR == "Omble de fontaine")) %>% 
                                filter(#Method %in% c("ASV"),
                            Location == "Avant-pays",
                            #NomFR != "Omble de fontaine",
                            SeqType == "sample") %>% 
                     mutate(DiffInv01 = ifelse(DiffInv %in% c("Absent trad et ADNe", "Present trad et ADNe"), 1,
                                        ifelse(DiffInv %in% c("Present ADNe seul", "Present trad seul"),0,NA)),
                            #DiffInv01.pres = ifelse(DiffInv %in% c("Present trad et ADNe"), 1,
                            #                   ifelse(DiffInv %in% c("Present ADNe seul", "Present trad seul"),0,NA)),
                            #CompTrad =  ifelse(DiffInv %in% c("Present trad et ADNe"), 1,
                            #                   ifelse(DiffInv %in% c("Present trad seul"),0,NA)), 
                            PresenceADNe = factor(PresenceADNe),
                            NameAssign.99 = factor(NameAssign.99),
                            Bassin = factor(Bassin),
                            CatSite = factor(CatSite),
                            VolSample.std = scale(Volume.x),
                            CorrFiltre.std = scale(CorrFiltre),
                            CorrFiltre.fct = factor(CorrFiltre),
                            VolLac.std = scale(Volume.y)) 

Sample.mod.data$NomFR %>% unique()
Sample.mod.data$NomLac %>% unique()
names(Sample.mod.data)



sum(Sample.mod.data$DiffInv01) / nrow(Sample.mod.data)
sum(DATA2$DiffInv01) / nrow(DATA)

mSample.m0 <- glmer(DiffInv01 ~ NameAssign.99 + 
                        VolSample.std +  
                        CorrFiltre.fct + #I(CorrFiltre.std^2)+
                        VolLac.std + I(VolLac.std^2)+ 
                        CatSite + 
                        Bassin +
                        #CatSite:VolSample.std +  
                        NameAssign.99:CatSite +
                        NameAssign.99:Bassin+
                  
                        (1|NomLac) + (1|Sample), 
            family = "binomial", 
            control=glmerControl(optimizer="bobyqa"),
            data = Sample.mod.data)

summary(mSample.m0)


plot(allEffects(m1))



# Maintenant en regroupant les lacs

Sample.mod.data2 <- Sample.mod.data %>% group_by(NomLac, NameAssign.99, VolLac.std, Bassin, Presence) %>% 
                                        summarise(Nsample = length(unique(Sample)),
                                                  Npel = length(unique(Sample[CatSite == "PEL"])),
                                                  Nriv = length(unique(Sample[CatSite == "RIV"])),
                                                  Nread = sum(N),
                                                  Nsample1 = length(unique(Sample[N>=1])),
                                                  VolSample = sum(VolSample.std),
                                                  CorrFiltre = sum(CorrFiltre),
                                                  PresenceADNe = sum(as.numeric(as.character(PresenceADNe)))
                                                  )    %>% 
                                        group_by() %>% 
                                        mutate(PresenceADNe = ifelse(PresenceADNe >=1, 1 ,0),
                                               DiffInv = ifelse(PresenceADNe == Presence, 
                                                                ifelse(PresenceADNe == 0 , "Absent trad et ADNe", "Present trad et ADNe"),
                                                                ifelse(PresenceADNe == 0, "Present trad seul", "Present ADNe seul")),
                                               DiffInv01 = ifelse(DiffInv %in% c("Absent trad et ADNe", "Present trad et ADNe"), 1,
                                                                  ifelse(DiffInv %in% c("Present ADNe seul", "Present trad seul"),0,NA)),
                                               VolLac.std = as.numeric(VolLac.std),
                                               Nsample.std = scale(Nsample),
                                               VolSample.std = scale(VolSample),

                                               CorrFiltre.std = scale(CorrFiltre))  %>% 
                                        as.data.frame()

str(Sample.mod.data2)  
        
mLake.m0 <- glmer(DiffInv01 ~ #NameAssign.99 + 
                              #VolLac.std + #I(VolLac.std^2) +
                              #Nsample + #I(Nsample^2) +
                              #Bassin +

                              #Nsample:VolLac.std +
                              #CatSite:VolSample.std +  
                              #NameAssign:CatSite +
                              #NameAssign:Bassin+
                              #Nsample +
                              #NameAssign:Nsample +
                              (1|NomLac), 
                 family = "binomial", 
                 control=glmerControl(optimizer="bobyqa"),
                 data = Sample.mod.data2)


mLake.m0 <- glmer(DiffInv01 ~ NameAssign.99 + 
                    VolLac.std + #I(VolLac.std^2) +
                    Nsample.std + #I(Nsample^2) +
                    #Bassin +
                    factor(Presence) +
                    NameAssign.99:VolLac.std +
                    #Nsample:VolLac.std +
                    #CatSite:VolSample.std +  
                    #NameAssign:CatSite +
                    #NameAssign:Bassin+
                    #Nsample +
                    #NameAssign:Nsample +
                    (1|NomLac), 
                  family = "binomial", 
                  control=glmerControl(optimizer="bobyqa"),
                  data = Sample.mod.data2)



mLake.m0 <- glmer(PresenceADNe ~ #NameAssign.99 + 
                    VolLac.std + #I(VolLac.std^2) +
                    #Nsample.std + #I(Nsample^2) +
                    VolSample.std +
                    #CorrFiltre.std +
                    #VolLac.std:Nsample +
                    #Bassin +
                    #factor(Presence) +
                    #NameAssign.99:VolLac.std +
                    #Nsample:VolLac.std +
                    #CatSite:VolSample.std +  
                    #NameAssign:CatSite +
                    #NameAssign:Bassin+
                    #Nsample +
                    #NameAssign:Nsample +
                    (1|NomLac), 
                  family = "binomial", 
                  control=glmerControl(optimizer="bobyqa"),
                  data = Sample.mod.data2 %>% filter(NameAssign.99 == "Salvelinus sp.",
                                                     Presence == 1))

mLake.m0 <- glm(PresenceADNe ~ #NameAssign.99 + 
                    VolLac.std #+ #I(VolLac.std^2) +
                    #Nsample.std + #I(Nsample^2) +
                    #VolSample.std# +
                    #CorrFiltre.std
                    #VolLac.std:Nsample +
                    #Bassin +
                    #factor(Presence) +
                    #NameAssign.99:VolLac.std +
                    #Nsample:VolLac.std +
                    #CatSite:VolSample.std +  
                    #NameAssign:CatSite +
                    #NameAssign:Bassin+
                    #Nsample +
                    #NameAssign:Nsample +
                  #(1|NomLac), 
                  #family = "binomial", 
                  #control=glmerControl(optimizer="bobyqa"),
                  ,data = Sample.mod.data2 %>% filter(NameAssign.99 == "Salvelinus sp.",
                                                     Presence == 1))

summary(mLake.m0)

plot(allEffects(mLake.m0))

# Modèle pour les CPUE/BPUE

Sample.mod.data.Etal <- Sample.mod.data %>% filter(NomLac %in% LacInv2$NomLac,
                                                   NameAssign.99 %in% LacInv2$Espece#,
                                                   #str_detect(Sample, "Sample_Edo_R_p1.G1") == FALSE, #7
                                                   #str_detect(Sample, "p1.D3") == FALSE, #19
                                                   #str_detect(Sample, "p1.G3") == FALSE, #22
                                                   #str_detect(Sample, "p1.B4") == FALSE, # 25
                                                   #str_detect(Sample, "p1.B2") == FALSE, #10
                                                   #str_detect(Sample, "p1.E2") == FALSE #13) %>%          
                                                   ) %>% 
                    left_join(LacInv2 %>% filter(Peche == "Verveux",
                                                 Mesure == "BPUE") %>% 
                                          select(NomLac, Espece, Value, Density.s), 
                              by = c("NomLac" = "NomLac", "NameAssign.99" = "Espece")) %>% 
                         mutate(Density.s.std = scale(Density.s),
                                Value.std = scale(Value),
                                Surface.std = scale(Surface),
                                NameAssign.99 = factor(NameAssign.99),
                                VolSample.fct = factor(VolSample.std),
                                VolFiltre = scale(VolSample.std / CorrFiltre))


Sample.mod.data.Etal$Bassin %>% unique()
Sample.mod.data.Etal$NomLac %>% unique()

Sample.mod.data.Etal$Surface.std %>% unique()
Sample.mod.data.Etal$Value.std %>% unique()

Sample.mod.data.Etal$NameAssign.99 %>% unique()
Sample.mod.data.Etal %>% group_by(VolFiltre, CatSite) %>% summarise(N = n())
Sample.mod.data.Etal %>% group_by(CorrFiltre) %>% summarise(N = n())

plot(Sample.mod.data.Etal %>% filter(Value > 0) %>% select(Value, Surface))

mEtal01.m0 <- glmer(PresenceADNe ~ 
                      VolFiltre +
                      Surface.std + 
                      
                      CatSite + 
                      
                      Value.std + 
                      
                      Value.std:Surface.std +

                      (1|NomLac) + 
                      (1|NameAssign.99) + 
                      (1|Sample), 
                    family = "binomial", 
                    control=glmerControl(optimizer="bobyqa"),
                    data = Sample.mod.data.Etal %>% filter(Value > 0))

summary(mEtal01.m0)

mEtal01.m1 <- glmer(PresenceADNe ~ 
                      VolFiltre +
                      Surface.std + 
                      
                      CatSite + 
                      
                      Value.std + 
                      
                     
                      (1|NomLac) + 
                      (1|NameAssign.99) + 
                      (1|Sample), 
                    family = "binomial", 
                    control=glmerControl(optimizer="bobyqa"),
                    data = Sample.mod.data.Etal %>% filter(Value > 0))

summary(mEtal01.m1)

mEtal01.m2 <- glmer(PresenceADNe ~ 
                      VolFiltre +
                      
                      CatSite + 
                      
                      Value.std + 
                      
                      (1|NomLac) + 
                      (1|NameAssign.99) + 
                      (1|Sample), 
                    family = "binomial", 
                    control=glmerControl(optimizer="bobyqa"),
                    data = Sample.mod.data.Etal %>% filter(Value > 0))

summary(mEtal01.m2)

mEtal01.m3 <- glmer(PresenceADNe ~ 
                      VolFiltre +
                      
                      CatSite + 

                      (1|NomLac) + 
                      (1|NameAssign.99) + 
                      (1|Sample), 
                    family = "binomial", 
                    control=glmerControl(optimizer="bobyqa"),
                    data = Sample.mod.data.Etal %>% filter(Value > 0))

summary(mEtal01.m3)

plot(All)

pred.m1 <- effect("VolFiltre", mEtal01.m3, xlevels=20, confidence.level=0.95)
pred.m1

resp.m1 <- summary(pred.m1, type="response")

names(resp.m1$effect)

plot(x=100, y=100, xlim = c(-1.6, 1.6), ylim = c(0,1), xlab = "Volume d'eau utilisé (std)", ylab = "Probabilité de détection par l'ADNe" )
lines(x = names(resp.m1$effect), y = resp.m1$effect, col = "black", lwd=1.5)
lines(x = names(resp.m1$effect), y = resp.m1$lower, col = "darkgray", lwd=1.5, lty = "dashed")
lines(x = names(resp.m1$effect), y = resp.m1$upper, col = "darkgray", lwd=1.5, lty = "dashed")

lines(x = colnames(resp.m1$effect), y = resp.m1$effect[2,], col = "blue")
lines(x = colnames(resp.m1$effect), y = resp.m1$lower[2,], col = "blue", lty = "dashed")
lines(x = colnames(resp.m1$effect), y = resp.m1$upper[2,], col = "blue", lty = "dashed")

plot(x = colnames(resp.m1$effect), y = resp.m1$effect[2,])


graph7.3 <- Sample.graph.red  %>% filter(Location == "Avant-pays",
                                         NomLac %in% Sample.mod.data.Etal$NomLac,
                                         NameAssign.99 %in% Sample.mod.data.Etal$NameAssign.99,
                                         NomFR != "Omble de fontaine *") %>% #View()
  ggplot(aes(x = NewNomLac, y = NomFR, shape = DiffInv)) + 
  geom_bin2d(aes(fill = Ncor), col = "gray", size = 0.5) + 
  #geom_point(size = 3, fill = "yellow",  stroke = 1, col = "gray30") + 
  
  scale_fill_distiller(palette = "Reds",
                       direction = 1,
                       #trans = "log10",  
                       na.value = "White", limits = c(0,110)) +
  
  geom_point(aes(shape = factor(Presence)), col = "gray20") +
  
  scale_shape_manual(values = 19, limits = "1", guide = "none") +
  
  
  #scale_shape_manual(values = c(3,6), name = "Comparaison avec\nl'inventaire traditionnel", limits = c("Present ADNe seul", "Present trad seul"), labels = c("Ajout", "Manquant")) +
  #scale_colour_manual(values = c("gray", "red"), name = "Traitement à la roténode", limits = c("0", "1"), labels = c("Non", "Oui")) +
  
  labs(title= NULL, x =NULL, y = NULL) +
  #guides(fill = guide_colourbar(title = "N moyen\nlectures", title.hjust = 0)) +
  guides(fill = guide_colourbar(title = "Indice d'abondance\nd'ADNe", title.hjust = 0)) +
  theme_bw ()+
  facet_grid(. ~ NewBassin, 
             scale = "free", space = "free", 
             labeller = label_value,
             switch = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, aes(size = 1/Volume/10000000)),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90),
        panel.spacing = unit(0, "lines"),
        strip.placement = "outside",
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.text.y = element_text(colour="black", hjust = 1),
        legend.position = "none") 

graph7.3

graph7.4 <- Sample.graph.red  %>% filter(Location == "Avant-pays",
                                         NomLac %in% Sample.mod.data.Etal$NomLac,
                                         NameAssign.99 %in% Sample.mod.data.Etal$NameAssign.99,
                                         NomFR != "Omble de fontaine *") %>% 
  left_join(LacInv2 %>% filter(Peche == "Verveux",
                               Mesure == "BPUE") %>% 
              select(NomLac, Espece, Value, Density.s), 
            by = c("NomLac" = "NomLac", "NameAssign.99" = "Espece")) %>% 
  filter(DiffInv %in% c(c("Present trad seul","Present trad et ADNe"))) %>% 
  #View()
  ggplot(aes(x = NewNomLac, y = NomFR, shape = DiffInv)) + 
  geom_bin2d(aes(fill = Ncor), col = "gray", size = 0.5) + 
  #geom_point(size = 3, fill = "yellow",  stroke = 1, col = "gray30") + 
  
  scale_fill_distiller(palette = "Reds",
                       direction = 1,
                       #trans = "log10",  
                       na.value = "White", limits = c(0,110)) +
  
  #scale_shape_manual(values = c(3,6), name = "Comparaison avec\nl'inventaire traditionnel", limits = c("Present ADNe seul", "Present trad seul"), labels = c("Ajout", "Manquant")) +
  #scale_colour_manual(values = c("gray", "red"), name = "Traitement à la roténode", limits = c("0", "1"), labels = c("Non", "Oui")) +
  
  labs(title= NULL, x =NULL, y = NULL) +
  #guides(fill = guide_colourbar(title = "N moyen\nlectures", title.hjust = 0)) +
  guides(fill = guide_colourbar(title = "Indice d'abondance\nd'ADNe", title.hjust = 0)) +
  theme_bw ()+
  facet_grid(. ~ NewBassin, 
             scale = "free", space = "free", 
             labeller = label_value,
             switch = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, aes(size = 1/Volume/10000000)),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90),
        panel.spacing = unit(0, "lines"),
        strip.placement = "outside",
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.text.y = element_text(colour="black", hjust = 1),
        legend.position = "none",
        panel.background = element_rect(fill = 'darkgray'),
        panel.grid.major = element_blank()
        ) 

graph7.4

graph7.4.1 <- Sample.graph.red  %>% filter(Location == "Avant-pays",
                                         NomLac %in% Sample.mod.data.Etal$NomLac,
                                         NameAssign.99 %in% Sample.mod.data.Etal$NameAssign.99,
                                         NomFR != "Omble de fontaine *") %>% 
  #mutate(Ncor = 1) %>% 
  left_join(LacInv2 %>% filter(Peche == "Verveux",
                               Mesure == "BPUE") %>% 
              select(NomLac, Espece, Value, Density.s), 
            by = c("NomLac" = "NomLac", "NameAssign.99" = "Espece")) %>% 
  #filter(Value > 0) %>% 
  #View()
  ggplot(aes(x = NewNomLac, y = NomFR, shape = DiffInv)) + 
  geom_bin2d(aes(fill = DiffInv), col = "gray", size = 0.5) + 
  #geom_point(size = 3, fill = "yellow",  stroke = 1, col = "gray30") + 
    
  scale_fill_manual(values = c("white", "white"),
                    limits = c("Present trad seul","Present trad et ADNe")) +
  
  #scale_shape_manual(values = c(3,6), name = "Comparaison avec\nl'inventaire traditionnel", limits = c("Present ADNe seul", "Present trad seul"), labels = c("Ajout", "Manquant")) +
  #scale_colour_manual(values = c("gray", "red"), name = "Traitement à la roténode", limits = c("0", "1"), labels = c("Non", "Oui")) +
  
  labs(title= NULL, x =NULL, y = NULL) +
  #guides(fill = guide_colourbar(title = "N moyen\nlectures", title.hjust = 0)) +
 # guides(fill = guide_colourbar(title = "Indice d'abondance\nd'ADNe", title.hjust = 0)) +
  theme_bw ()+
  facet_grid(. ~ NewBassin, 
             scale = "free", space = "free", 
             labeller = label_value,
             switch = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, aes(size = 1/Volume/10000000)),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90),
        panel.spacing = unit(0, "lines"),
        strip.placement = "outside",
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.text.y = element_text(colour="black", hjust = 1),
        legend.position = "none",
        panel.background = element_rect(fill = 'darkgray'),
        panel.grid.major = element_blank()
  ) 

graph7.4.1

ggsave(filename = file.path(get.value("result.FINAL"), "CompEtal.FauxPositif.AV.ASV.12S.png"),
       width = 5, height = 5,
       plot = graph7.3
)

ggsave(filename = file.path(get.value("result.FINAL"), "CompEtal.AV.ASV.12S.png"),
       width = 5, height = 5,
       plot = graph7.4
)

ggsave(filename = file.path(get.value("result.FINAL"), "CompEtal.AV.ASV.12S.vide.png"),
       width = 5, height = 5,
       plot = graph7.4.1
)



res <- vector()

for(x in 1:50000){
 
  
val <- DATA2 %>% #filter(DiffInv = "Absent trad et ADNe") %>% 
                 mutate(PresenceADNe = ifelse(Nread >= x, 1, 0),
                 DiffInv = ifelse(PresenceADNe == Presence, 
                                  ifelse(PresenceADNe == 0 , "Absent trad et ADNe", "Present trad et ADNe"),
                                  ifelse(PresenceADNe == 0, "Present trad seul", "Present ADNe seul")),
                 DiffInv01 = ifelse(DiffInv %in% c("Absent trad et ADNe", "Present trad et ADNe"), 1,
                                    ifelse(DiffInv %in% c("Present ADNe seul", "Present trad seul"),0,NA))) %>% 
                summarise(voila = sum(DiffInv01)) %>% pull()
                          
  
#val <- DATA2 %>% filter(Nsample1 >= x) %>% 
#                 pull(DiffInv01) %>% sum() 

res[x] <- val / nrow(DATA2)

}

plot(y= res[1:15000], x=1:15000, type = "l", xlab= "N min reads", ylab = "% concordance trad / ADNe")



sum(DATA2$DiffInv01) / nrow(DATA2)

# Trying to model with density

DATA3 <- DATA %>% filter(NomLac %in% LacInv2$NomLac,
                          NameAssign %in% LacInv2$Espece,
                          NameAssign %nin% c("Luxilus cornutus", "Margariscus margarita")
                          ) %>% 
                   left_join(LacInv2 %>% filter(Peche == "Verveux", Mesure == "CPUE") %>% 
                                         select(NomLac, Value, Density.s), by = "NomLac") %>% 
  mutate( NameAssign = factor(NameAssign),
          Bassin = factor(Bassin),
          CatSite = factor(CatSite),
          VolSample.std = scale(Volume.x),
          CorrFiltre.std = scale(CorrFiltre),
          VolLac.std = scale(Volume.y),
          Density.s.std = scale(Density.s),
          CPUE.std = scale(Value)) #%>% pull(NameAssign) %>% unique()



m1 <- glmer(PresenceADNe ~ NameAssign +
            #VolLac.std +# CPUE.std +
              Density.s.std +
              Density.s.std:NameAssign +
              #NameAssign:CPUE.std +
              #CatSite:VolSample.std +  
              #NameAssign:CatSite +
              #NameAssign:Bassin+
              
              (1|NomLac) + (1|Sample), 
            family = "binomial", 
            control=glmerControl(optimizer="bobyqa"),
            data = DATA3)


summary(m1)

plot(allEffects(m1))


# TO do the same but with PEL vs RIV

CompPAtrad.graph.fct <- function(data) {   
    
  bind_rows(data %>% filter(Method %in% c("ASV", NA),
                                  Locus %in% c("12s", NA),
                                  NomFR != "Omble de fontaine"),
            data %>% filter(Method == "ASV",
                                  Locus == "cytB.R1",
                                  NomFR == "Omble de fontaine")# %>% View()
  ) %>%  
    mutate(NomFR = ifelse(NomFR == "Omble de fontaine", paste(NomFR, "*"), NomFR),
           NomFR = factor(NomFR, levels = rev(c("Omble de fontaine *", "Omble chevalier", "Touladi", "Salvelinus sp.", 
                                                "Mulet à cornes", "Mulet perlé", "Tête-de-boule", "Ventre rouge du nord", "Méné à nageoires rouges", "Méné jaune", "Museau noir", "Naseux des rapides", "Ouitouche",
                                                "Meunier noir",
                                                "Épinoche à cinq épines", "Épinoche à neuf épines",
                                                "Chabot à tête plate",
                                                "Fondule barré",
                                                "Éperlan arc-en-ciel",
                                                "Achigan à petite bouche", "Crapet de roche", "Crapet-soleil",
                                                "Barbotte brune", 
                                                "Perchaude", "Doré jaune", "Fouille-roche zébré",
                                                "Grand brochet")
           ))) %>% 
    filter(Location == "Avant-pays") %>% 
    #Method %in% c("ASV", NA),
    #Locus %in% c("12s", NA)) %>%# View()
    ggplot(aes(x=NewNomLac, y = NomFR, fill = DiffInv)) +
    geom_bin2d(col = "black") +
    labs(title= NULL, x =NULL, y = NULL) +
    #scale_y_discrete(position = "right") +
    scale_fill_manual(limits = c("Present trad et ADNe",  "Absent trad et ADNe", "Present trad seul", "Present ADNe seul"),  
                      values = c("green3", 
                                 "darkseagreen1", 
                                 "dodgerblue3", #"darkorange1", 
                                 "goldenrod1")) + 
    guides(fill = guide_legend(title = NULL)) + 
    theme_bw ()+
    facet_grid(Family ~ NewBassin, 
               scale = "free", space = "free", 
               labeller = label_value,
               switch = NULL) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.ticks = element_blank(),
          strip.text.y = element_text(angle = 0),
          strip.text.x = element_text(angle = 90),
          panel.spacing = unit(0, "lines"),
          strip.placement = "outside",
          strip.background = element_rect(colour = "black", fill = "white"),
          axis.text.y = element_text(colour="black", hjust = 1)
    ) #+ coord_flip()
  

}



# The three graphes                     
ggarrange(CompPAtrad.graph.fct(Sample.data %>% CompPAtrad.data.fct()),
          ggarrange(CompPAtrad.graph.fct(Sample.data %>% filter(CatSite == "RIV") %>% CompPAtrad.data.fct()) + theme(legend.position="none"),
                    CompPAtrad.graph.fct(Sample.data %>% filter(CatSite == "PEL") %>% CompPAtrad.data.fct() %>% filter(NomLac %nin% c("Berube", "Etienne", "Hamel Ouest"))) + theme(legend.position="none"),
                    ncol =2, labels = c("B. Échantillons riverains", "C. Échantillons pélagiques"),  hjust = 0),
          nrow=2, labels = "A. Tous les échantillons", hjust = 0 
   )


CompPAtrad.fct(CompPAtrad.cat %>% filter(CatSite != "PEL"))
CompPAtrad.fct(CompPAtrad.cat %>% filter(CatSite %in% c("RIV", NA)))

CompPAtrad.cat %>% View()

CompPAtrad %>% filter(NameAssign == "Salvelinus fontinalis",
                      Location == "Avant-pays") %>% 
  ggplot(aes(y=NomLac, x = Method, fill = DiffInv)) +
  geom_bin2d(col = "black") +
  labs(title= NULL, x ="Lac", y = "Assignation") +
  scale_fill_manual(limits = c("Present trad et ADNe", "Present trad seul", "Present ADNe seul", "Absent trad et ADNe"),  
                    values = c("green", "red", "yellow", "white")) + 
  theme_bw()+
  facet_grid(. ~ Locus, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0)) #+ coord_flip()

# # Idem mais pour les salvelidées
# 
# Sample.data %>% filter(str_detect(NameAssign, "Salvelinus"),
#                        Location == "Avant-pays") %>% 
#   group_by(NameAssign, Location, NomLac, Method, Locus) %>% 
#   summarise(N = sum(N)) %>% 
#   mutate(PresenceADNe = ifelse(N >=1, 1, 0)) %>% 
#   select(-c(N)) %>% 
#   spread(NameAssign, value = PresenceADNe, fill = 0) %>%
#   rename(Salvelinus = "Salvelinus.ADNe", "Salvelinus fontinalis" = "Salvelinus fontinalis.ADNe") %>% 
#   left_join(LacInv1 %>% filter(str_detect(Espece, "Salvelinus")) %>% 
#                         spread(Espece, Presence) %>% 
#                         select(-c(InvLac, Rotenode, Location)),
#                         by = c("NomLac" = "NomLac")) %>% 
#   gather(names(.) %>% str_subset("Salvelinus"), key = Salvelinus, value = Presence) %>% 
#   mutate(Salvelinus = factor(Salvelinus, levels = c("Salvelinus fontinalis", "Salvelinus namaycush", "Salvelinus alpinus", "Salvelinus fontinalis.ADNe", "Salvelinus.ADNe"))) %>% 
#   filter(Locus %in% c("12s", "cytB.R1"),
#        Method %in% c("ASV", "OTU"),
#        Location == "Avant-pays") %>% 
#   ggplot(aes(y=NomLac, x = Salvelinus, fill = factor(Presence))) +
#   geom_bin2d(col = "black") +
#   labs(title= NULL, x ="", y = "Lac") +
#   scale_fill_manual(limits = c("0", "1"),  
#                     values = c("white","green")) + 
#   theme_bw()+
#   facet_grid(Locus ~ Method, scale = "free", space = "free") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#         axis.ticks.y = element_blank(),
#         strip.text.y = element_text(angle = 0))
#   
# LacInv1 %>% head()  
# 
# 
# # ## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# ## biocLite("BiocUpgrade") ## you may need this
# biocLite("ggtree")
# 
# 
# 
#   pull(NameAssign) %>% unique()tion

# Comparison riverain, pélagique

Sample.data %>% filter(Data == "ASVtab.12s") %>% select(Sample, Cat) %>% group_by(Cat) %>% summarise(Ntot = length(unique(Sample)))


Comp.RivPel <-Sample.data %>% filter(Locus %in% c("12s"),
                                     Method %in% c("ASV"),#Data == "ASVtab.12s",
                       N>=1,
                       Cat %in% c("PEL", "RIV"),
                       str_detect(Assign, "Teleostei") == T,
                       Level >= 5) 

Comp.RivPel$Famille <- sapply(str_split(Comp.RivPel$Assign, ";"),`[`,5)


Comp.RivPel <-Comp.RivPel %>% 
  group_by(Locus, Method, Cat, Level, Famille, NameAssign) %>% 
  summarise(Nesp = length(unique(Sample))) %>% 
  spread(Cat, Nesp) 

Comp.RivPel <- na.replace(Comp.RivPel,0)

Comp.RivPel$P.rel <- Comp.RivPel$PEL / 67 *100
Comp.RivPel$R.rel <- Comp.RivPel$RIV / 95 * 100

#Comp.RivPel %>% View()

Comp.RivPel %>% ggplot(aes(P.rel, R.rel, col = Famille)) + 
                geom_point(size = 2) + 
                geom_abline(slope =1, intercept = 0) + 
                geom_text(aes(label=NameAssign), check_overlap = TRUE, vjust=0, nudge_x = 0.01, size=3, hjust = 0)+
                theme_bw() +
                scale_x_continuous(limits=c(0,60))+
                scale_y_continuous(limits=c(0,60))+
                xlab("Échantillons pélagiques (%)") + ylab("Échantillons riverains (%)")
               # facet_grid(Locus ~ Method)




# Correlation with inventaire -----------------------------------------------------

# Etalonnage

graph8 <- Sample.mod.data.Etal %>%  mutate(N = ifelse(is.na(N), 0, N),
                                 Nlog = ifelse(N == 0 , NA, 
                                               ifelse(N == 1, 0.5, log2(N))),
                                 Nlog.cor = ifelse(is.na(Nlog), 0, Nlog / CorrFiltre / Volume.x * 1000)) %>% 
                             filter(#Value > 0,
                                    NomLac %nin% c("Giron", "Parker"),
                                   NomFR %nin% c("Perchaude", "Omble de fontaine")) %>% #View()
  ggplot(aes(x = Value, y = Nlog.cor, col = CatSite)) +
  geom_smooth(method = "lm", se = FALSE) +                              
  geom_count(show.legend = F) + 
  #stat_cor(method = "spearman", cex= 3.5) +
                                facet_wrap(~ NomFR, scale = "free", nrow=2) +
  scale_colour_manual(values = c("blue","red"), name = "Type d'échantillon", limits = c("RIV", "PEL"), labels = c("Riverain", "Pélagique")) +
  #scale_colour_manual(values = c("gray", "red"), name = "Traitement à la roténode", limits = c("0", "1"), labels = c("Non", "Oui")) +
  
  labs(title= NULL, x = "BPUE", y = "Indice d'abondance d'ADNe") +
  #guides(fill = guide_colourbar(title = "N moyen\nlectures", title.hjust = 0)) +
  guides(fill = guide_colourbar(title = "Indice d'abondance\nd'ADNe", title.hjust = 0)) +
  
  theme_bw() + theme(strip.background = element_rect(fill="white"),
                     legend.position = "top")


graph8


ggsave(filename = file.path(get.value("result.FINAL"), "CompAbondance.ASV.12S.png"),
       width = 7, height = 5,
       plot = graph8
)


# Pêche
library(GGally)
ggpairs( LacPeche %>% spread(Mesure, Value) %>% select(BPUEhm, BPUEjp, CPUEhm, CPUEjp))

LacPeche$NomLac %>% unique()

Sample.mod.data.Peche <- Sample.mod.data %>% filter(NomLac %in% LacPeche$NomLac,
                                                    NameAssign.99 %in% c(LacPeche$Espece %>% unique, "Salvelinus sp.")) %>% #pull(NomLac) %>% unique()
   mutate(NameAssign.99.join = NameAssign.99,
          NameAssign.99.join = str_replace(NameAssign.99.join, "Salvelinus sp.", "Salvelinus fontinalis")) %>% 
  left_join(LacPeche %>% filter(Mesure == "BPUEjp") %>% 
              select(NomLac, Espece, Value), 
            by = c("NomLac" = "NomLac", "NameAssign.99.join" = "Espece")) %>% 
  mutate(Value = ifelse(is.na(Value), 0, Value),
         Marker = ifelse(Locus == "12s", "12S - Salvelinus sp.", "cytB - Omble de fontaine"))


graph9 <- Sample.mod.data.Peche %>%  mutate(N = ifelse(is.na(N), 0, N),
                                           Nlog = ifelse(N == 0 , NA, 
                                                         ifelse(N == 1, 0.5, log2(N))),
                                           Nlog.cor = ifelse(is.na(Nlog), 0, Nlog / CorrFiltre / Volume.x * 1000)) %>% 
  #group_by(NomFR, CatSite, NomLac)%>%
  #summarise(Nlog.cor = mean(Nlog.cor),
  #          Value = mean(Value)) %>% 
  filter(NomFR != "Achigan à petite bouche",
         NomLac %nin% c("Giron", "Parker", "Peche (a la)"),
         Nlog.cor > 0
    ) %>% #View()
  ggplot(aes(x = Value, y = Nlog.cor, col = CatSite)) +
  geom_smooth(method = "lm", se = FALSE) +                              
  geom_count(show.legend = F) + 
  #stat_cor(method = "spearman", cex= 3.5) +
  facet_wrap(~ Marker, nrow=1) +
  scale_colour_manual(values = c("blue","red"), name = "Type d'échantillon", limits = c("RIV", "PEL"), labels = c("Riverain", "Pélagique")) +
  #scale_colour_manual(values = c("gray", "red"), name = "Traitement à la roténode", limits = c("0", "1"), labels = c("Non", "Oui")) +
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2,2.5), labels = c("0,0", "0,5", "1,0", "1,5", "2,0", "2,5"))+  
  labs(title= NULL, x = "BPUE (heures moyennes)", y = "Indice d'abondance d'ADNe") +
  #guides(fill = guide_colourbar(title = "N moyen\nlectures", title.hjust = 0)) +
  guides(fill = guide_colourbar(title = "Indice d'abondance\nd'ADNe", title.hjust = 0)) +
  
  theme_bw() + theme(strip.background = element_rect(fill="white"),
                     legend.position = "top")

graph9

ggsave(filename = file.path(get.value("result.FINAL"), "CompPeche.ASV.12S.png"),
       width = 5, height = 3,
       plot = graph9
)


RES.COR.INV2 <- expand.grid(Espece = c( "Ameiurus nebulosus",
                                        "Ambloplites rupestris",
                                        "Chrosomus eos",
                                        "Culaea inconstans",
                                        "Perca flavescens",
                                        "Salvelinus fontinalis",
                                        "Semotilus atromaculatus"#,
                                        #"Margariscus margarita",
                                        #"Luxilus cornutus"
                                        ),
                            Peche = c("Verveux", "Alaska"),
                            Mesure = c("CPUE", "BPUE"),
                            Cat = c("RIV", "PEL", "T"),
                            CorMax = NA,
                            CorMax.pvalue = NA,
                            CorMed = NA,
                            CorMed.pvalue = NA,
                            CorMean = NA,
                            CorMean.pvalue = NA)



for(x in 1:nrow(RES.COR.INV2)){

  DATA <- Sample.data %>% filter(NomLac %in% LacInv2$NomLac,
                       #Cat == "R",
                       Data %in% "ASVtab.12s",
                       NameAssign %in% c(LacInv2$Espece, "Salvelinus")) %>% 
                mutate(Espece = ifelse(NameAssign == "Salvelinus", "Salvelinus fontinalis", NameAssign)) %>% 
                select(Sample, NomLac, Cat, Nsite, Espece, N) %>%
                left_join(LacInv2) %>% 
                filter(Peche == as.character(RES.COR.INV2[x,"Peche"]),
                       Mesure == as.character(RES.COR.INV2[x,"Mesure"]),
                       Espece == as.character(RES.COR.INV2[x,"Espece"])) %>% 
                group_by(NomLac)
              

  if(RES.COR.INV2[x,"Cat"] == "T") {
    
  DATA <- DATA %>%  summarise(Nmax = max(N),
                              Nmed = median(N),
                              Nmean = mean(N),
                              Value = unique(Value)) 
    
  } else {
    DATA <- DATA %>% filter(Cat == RES.COR.INV2[x,"Cat"]) %>%  
      summarise(Nmax = max(N),
                Nmed = median(N),
                Nmean = mean(N),
                Value = unique(Value)) 
    
  }
  
  RES1 <- cor.test(DATA$Nmax, DATA$Value, method = "spearman")
  RES2 <- cor.test(DATA$Nmed, DATA$Value, method = "spearman")
  RES3 <- cor.test(DATA$Nmean, DATA$Value, method = "spearman")
  
  RES.COR.INV2[x,"CorMax"]        <- RES1$estimate
  RES.COR.INV2[x,"CorMax.pvalue"] <- RES1$p.value
  RES.COR.INV2[x,"CorMed"]        <- RES2$estimate
  RES.COR.INV2[x,"CorMed.pvalue"] <- RES2$p.value
  RES.COR.INV2[x,"CorMean"]        <- RES3$estimate
  RES.COR.INV2[x,"CorMean.pvalue"] <- RES3$p.value
}


plot(RES.COR.INV2$CorMax, RES.COR.INV2$CorMean)

RES.COR.INV2.graph <- RES.COR.INV2 %>% gather(names(.) %>% str_subset("CorM"), key= Stat, value = Value)

RES.COR.INV2.graph$Methode <- paste(RES.COR.INV2$Peche, RES.COR.INV2$Mesure, RES.COR.INV2$Cat)

RES.COR.INV2.graph %>% filter(Stat %in% c("CorMax", "CorMed", "CorMean")) %>% 
  
  ggplot(aes(x = Espece, y = Methode, fill = Value)) + 
  geom_bin2d(color = "gray")+
  scale_fill_distiller(palette = "Spectral") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(Cat ~
             Stat, scale= "free")
 
# Avec les corrections

Sample.data.cor <- Sample.data %>% left_join(LacSample  %>% select(NomLac, Surface, Volume) %>% rename (Volume = "LacVol")) %>% 
                                   mutate(Nlog = ifelse(N == 0 , 0, log2(N)),
                                          Nlog.cor = Nlog / CorrFiltre / Volume * 1000,
                                          Nlog.cor.vol = Nlog.cor / LacVol / 1000000000) %>% 
                                   gather(Nlog, Nlog.cor, Nlog.cor.vol, key = Correction, value = Ncor)


RES.COR.LOG.INV2 <- expand.grid(Espece = c( "Ameiurus nebulosus",
                                        "Ambloplites rupestris",
                                        "Chrosomus eos",
                                        "Culaea inconstans",
                                        "Perca flavescens",
                                        "Salvelinus fontinalis",
                                        "Semotilus atromaculatus"#,
                                        #"Margariscus margarita",
                                        #"Luxilus cornutus"
                                        ),
                                        Peche = c("Verveux", "Alaska"),
                                        Mesure = c("CPUE", "BPUE"),
                                        Cat = c("RIV", "PEL", "T"),
                                        Correction = c("Nlog.cor",  "Nlog.cor.vol"),
                                        CorMax = NA,
                                        CorMax.pvalue = NA,
                                        CorMed = NA,
                                        CorMed.pvalue = NA,
                                        CorMean = NA,
                                        CorMean.pvalue = NA)


for(x in 1:nrow(RES.COR.LOG.INV2)){
  
  DATA <- Sample.data.cor %>%   filter(NomLac %in% LacInv2$NomLac,
         #Cat == "R",
         Data %in% "ASVtab.12s",
         NameAssign %in% c(LacInv2$Espece, "Salvelinus"),
         Correction == RES.COR.LOG.INV2[x,"Correction"]) %>% 
    mutate(Espece = ifelse(NameAssign == "Salvelinus", "Salvelinus fontinalis", NameAssign)) %>% 
    select(Sample, NomLac, Cat, Nsite, Espece, Ncor) %>%
    left_join(LacInv2) %>% 
    filter(Peche == as.character(RES.COR.LOG.INV2[x,"Peche"]),
           Mesure == as.character(RES.COR.LOG.INV2[x,"Mesure"]),
           Espece == as.character(RES.COR.LOG.INV2[x,"Espece"])) %>% 
    group_by(NomLac)
  
  
  if(RES.COR.LOG.INV2[x,"Cat"] == "T") {
    
    DATA <- DATA %>%  summarise(Nmax = max(Ncor),
                                Nmed = median(Ncor),
                                Nmean = mean(Ncor),
                                Value = unique(Value)) 
    
  } else {
    DATA <- DATA %>% filter(Cat == RES.COR.LOG.INV2[x,"Cat"]) %>%  
      summarise(Nmax = max(Ncor),
                Nmed = median(Ncor),
                Nmean = mean(Ncor),
                Value = unique(Value)) 
    
  }
  
  RES1 <- cor.test(DATA$Nmax, DATA$Value, method = "spearman")
  RES2 <- cor.test(DATA$Nmed, DATA$Value, method = "spearman")
  RES3 <- cor.test(DATA$Nmean, DATA$Value, method = "spearman")
  
  RES.COR.LOG.INV2[x,"CorMax"]        <- RES1$estimate
  RES.COR.LOG.INV2[x,"CorMax.pvalue"] <- RES1$p.value
  RES.COR.LOG.INV2[x,"CorMed"]        <- RES2$estimate
  RES.COR.LOG.INV2[x,"CorMed.pvalue"] <- RES2$p.value
  RES.COR.LOG.INV2[x,"CorMean"]        <- RES3$estimate
  RES.COR.LOG.INV2[x ,"CorMean.pvalue"] <- RES3$p.value
}



plot(RES.COR.LOG.INV2$CorMax, RES.COR.LOG.INV2$CorMean)

RES.COR.LOG.INV2.graph <- RES.COR.LOG.INV2 %>% gather(names(.) %>% str_subset("CorM"), key= Stat, value = Value.cor)

RES.COR.LOG.INV2.graph$Methode <- paste(RES.COR.LOG.INV2$Peche, RES.COR.LOG.INV2$Mesure, RES.COR.LOG.INV2$Cat)


RES.COR.LOG.INV2 %>% filter(#Stat %in% c("CorMax", "CorMed", "CorMean"),
                           Correction == "Nlog.cor") %>% View()

RES.COR.LOG.INV2.graph %>% filter(Stat %in% c("CorMax", "CorMed", "CorMean"),
                                  Correction == "Nlog.cor") %>% 
  
  ggplot(aes(x = Espece, y = Methode, fill = Value.cor)) + 
  geom_bin2d(color = "gray")+
  scale_fill_distiller(palette = "Spectral") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(Cat ~
            Stat + Correction, scale= "free")


RES.COR.LOG.INV2.graph %>% spread(Correction, Value.cor) %>% ggplot(aes(x = Nlog.cor, y = Nlog.cor.vol)) + geom_point()

RES.COR.LOG.INV2.graph %>% spread(Correction, Value.cor) %>%  
  mutate(diff = Nlog.cor - Nlog.cor.vol) %>% #head()
  #group_by(Peche, Mesure, Cat) %>% summarise(Mean = mean(diff, na.rm=T)) %>% View()
  ggplot(aes(x = diff)) + geom_histogram()

  
Sample.data.cor %>% filter(Correction == "Nlog.cor.vol",
                           Method == "ASV",
                           Locus == "12s",
                           NameAssign %in% c(LacInv2$Espece),
                           NomLac %in% LacInv2$NomLac) %>%
  mutate(Espece = NameAssign) %>% 
  left_join(LacInv2 %>% 
  filter(#Peche == "Verveux",
         Mesure == "BPUE")) %>% #View()
  mutate(Method = paste(Peche, Mesure, sep = "-")) %>% 
  group_by(NomLac, Cat, Espece, Method, Value) %>% 
  summarise(Ncor = mean(Ncor)) %>% 
  ggplot(aes(x=Value, y = Ncor, col = Cat)) + 
  geom_jitter() + 
  geom_smooth(method = "lm", se = F, na.rm=T)+
 # facet_wrap( ~ Espece, scale = "free") + theme_bw()
  facet_wrap( ~ Method + Espece, ncol = 6,  scale = "free")  + 
  theme_bw() + 
  theme(strip.text = element_text(size = 8))
  


head(Sample.data.cor)


Sample.data.cor %>% filter(Ncor>=1) %>%  ggplot(aes(Ncor)) + geom_histogram() + facet_grid(Locus ~ Method)


Sample.data.cor %>% spread(Correction, Ncor) %>% ggplot(aes(Nlog, Nlog.cor)) + geom_count()





LacPeche$Mesure %>% unique() 


RES.COR.PECHE <- expand.grid(Espece = c("Salvelinus fontinalis"#,
                                        #"Ambloplites rupestris",
                                        #"Chrosomus eos",
                                        #"Culaea inconstans",
                                        #"Perca flavescens",
                                        #"Salvelinus fontinalis",
                                        #"Semotilus atromaculatus"#,
                                        #"Margariscus margarita",
                                        #"Luxilus cornutus"
),
#Peche = c("Verveux", "Alaska"),
Mesure = c("CPUEjp", "CPUEhm", "BPUEjp", "BPUEhm"),
Cat = c("R", "P", "T"),
Cor = NA,
Cor.pvalue = NA)

for(x in 1:nrow(RES.COR.PECHE)){
  
  DATA <- Sample.data %>% filter(NomLac %in% LacPeche$NomLac,
                                 #Cat == "R",
                                 Data %in% "ASVtab.12s",
                                 Name.level %in% c(LacInv2$Espece, "Salvelinus")) %>% 
    mutate(Espece = ifelse(Name.level == "Salvelinus", "Salvelinus fontinalis", Name.level)) %>% 
    select(Sample, NomLac, Cat, Nsite, Espece, N) %>%
    left_join(LacPeche %>% select(-Lac)) %>% 
    filter(#Peche == as.character(RES.COR.INV2[x,"Peche"]),
           Mesure == as.character(RES.COR.PECHE[x,"Mesure"]),
           Espece == as.character(RES.COR.PECHE[x,"Espece"])) %>% 
    group_by(NomLac)
  
  
  if(RES.COR.PECHE[x,"Cat"] == "T") {
    
    DATA <- DATA %>%  summarise(Nmed = median(N),
                                Value = unique(Value)) 
    
  } else {
    DATA <- DATA %>% filter(Cat == RES.COR.PECHE[x,"Cat"]) %>%  
      summarise(Nmed = median(N),
                Value = unique(Value)) 
    
  }
  
  RES <- cor.test(DATA$Nmed, DATA$Value, method = "spearman")
  
  RES.COR.PECHE[x,"Cor"] <- RES$estimate
  RES.COR.PECHE[x,"Cor.pvalue"] <- RES$p.value
  
}


RES.COR.PECHE %>% View()

RES.COR.PECHE$Methode <- paste(RES.COR.PECHE$Mesure, RES.COR.PECHE$Cat)

RES.COR.PECHE %>% ggplot(aes(x = Mesure, y = Cat, fill = Cor)) + 
  geom_bin2d()+
  scale_fill_distiller(palette = "Spectral") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#

DATA <- Sample.data %>% filter(NomLac %in% LacPeche$NomLac,
                               #Cat == "R",
                               Data %in% "ASVtab.12s",
                               Name.level %in% c("Salvelinus")) %>% 
  mutate(Espece = ifelse(Name.level == "Salvelinus", "Salvelinus fontinalis", Name.level),
         N01 = ifelse(N>=1,1,0)) %>% 
  select(Sample, NomLac, Cat, Nsite, Espece, N, N01) %>%
  left_join(LacPeche %>% select(-Lac)) %>% 
  spread(Mesure, Value) 


DATA %>% head()

DATA$Cat <- factor(DATA$Cat)

mod1 <- glmer(N01 ~ Cat + BPUEhm + (1|NomLac), family = "binomial", data = DATA)

summary(mod1)


plot(allEffects(mod1))



# Samples - RIV vs PEL ----------------------------------------------------

ASVtab.12s.cor.bySP

DataSeq


Comp.res <- ASVtab.12s.cor.bySP %>% mutate(Lac = sapply(str_split(Sample, "_"), `[`, 3),
                                           Cat = sapply(str_split(Sample, "_"), `[`, 4),
                                           Puit = sapply(str_split(Sample, "_"), `[`, 5) %>% str_remove("p")) %>% 
                                    left_join(DataSeq %>% select(IbisID, Nsite), by = c("Puit" = "IbisID"))  %>% 
                                    select(Assign, N, Lac, Cat, Nsite) %>% 
                                    spread(key = Nsite, value = N, fill = 0) %>% 
  mutate(KEEP = ASV + OTU) %>% 
  filter(KEEP > 0) %>% select(-KEEP)spread()
  
  


# Phase I -----------------------------------------------------------------

LacInv3 %>% ggplot(aes(Nsp.x, Nsp.y)) + geom_count()

cor.test(LacInv3$Nsp.x, LacInv3$Nsp.y, method = "spearman")

# N species

# Pour les lacs d'avant-pays, compter le N minimal d'espèce, puis modèle en fonction inventaire trad


DATA <- Sample.data %>% mutate(Method = str_sub(Data,1,3),
                       Locus = Data %>% str_remove(paste0(Method, "tab."))) %>% 
   filter(str_detect(Assign, "Teleostei") == T,
                       Level >= 6,
                       Locus %in% c("12s"),
                       Method %in% c("OTU")) %>%
  mutate(Genus = sapply(str_split(Name.level, " "), `[`, 1)) %>% 
  filter(Genus != "Sebaste",
         N>=1) %>% 
  group_by(NomLac) %>% summarise(Nsp.ADNe = length(unique(Genus))) %>% 
  left_join(LacInv2, by = "NomLac")

DATA

DATA %>% 
  ggplot(aes(x = Nsp.y, y = Nsp.ADNe)) + geom_count()

DATA %>% 
  ggplot(aes(x = Nsp.ADNe)) + geom_histogram()


mod <- lm(Nsp.ADNe ~ Nsp, data = DATA)
summary(mod)

cor.test(DATA$Nsp.ADNe, DATA$Nsp.x, method = "spearman")
cor.test(DATA$Nsp.ADNe, DATA$Nsp.y, method = "spearman")


# PCA


DATA <- ASVtab.12s.cor.bySP %>% mutate(Lac = sapply(str_split(Sample, "_"), `[`, 3),
                                       Cat = sapply(str_split(Sample, "_"), `[`, 4),
                                       Puit = sapply(str_split(Sample, "_"), `[`, 5) %>% str_remove("p")) %>% 
                                left_join(DataSeq %>% select(IbisID, NomLac, CatSite, Nsite), by = c("Puit" = "IbisID"))  %>% 
                                filter(NomLac %in% LacInv2$NomLac,
                                      Level == 6,
                                       str_detect(Assign, "Teleostei")==T) %>%
                                mutate(Species = sapply(str_split(Assign, ";"), `[`, 7)) %>%
                                group_by() %>% 
                                select(Sample, Species, N) %>% 
                                mutate(N = ifelse(N>=1,1,0)) %>% 
                                spread(key = Species, value = N) 
DATA$KEEP <- rowSums(DATA[,-1])
  
DATA <- DATA %>% filter(KEEP >= 1) %>% select(-c(KEEP))

library(vegan)

str(DATA)

  
DATA.jac <- vegdist(DATA[,-1], method="jac") 

summary(DATA.jac)


library(phangorn)

?pcoa
DATA.jac.pcoa<-pcoa(DATA.jac)


RES <- cbind(DATA[,1], DATA.jac.pcoa$vectors)

RES <-  RES %>% mutate(Lac = sapply(str_split(Sample, "_"), `[`, 3),
                       Cat = sapply(str_split(Sample, "_"), `[`, 4))

RES %>% ggplot(aes(x = Axis.1, y =Axis.2, col = Lac, shape= Cat)) + 
  geom_point() + 
  #facet_wrap(~Lac) + 
  theme_bw()

# Extraction des résultats
DATA.jac.pcoa

# Représentation graphique
biplot.pcoa(DATA.jac.pcoa)



DATA.jac.upgma <- upgma(DATA.jac)

DATA.jac.upgma[["tip.label"]] <- paste(RES$Lac, RES$Cat, sep = "-")


plot(DATA.jac.upgma)



summary(DATA.cca)

View(DATA)


summary(DATA.jac.pcoa)

  mutate(Lac = sapply(str_split(Sample, "_"), `[`, 3),
                                           Cat = sapply(str_split(Sample, "_"), `[`, 4),
                                           Puit = sapply(str_split(Sample, "_"), `[`, 5) %>% str_remove("p")) %>% 
  left_join(DataSeq %>% select(IbisID, Nsite), by = c("Puit" = "IbisID"))  %>% 
  select(Assign, N, Lac, Cat, Nsite) %>% 
  spread(key = Nsite, value = N, fill = 0)



# Salvelinus detection ----------------------------------------------------


DATA <-  Sample.data %>% mutate(NameAssign = ifelse(NameAssign == "Salvelinus", "Salvelinus fontinalis", NameAssign)) %>% 
                  filter(str_detect(Assign, "Teleostei"),
  str_detect(NameAssign, " "),
  str_detect(NameAssign, "Gadus morhua") == FALSE,
  str_detect(NameAssign, "Salmo salar") == FALSE )%>%   
  select(NomLac, NameAssign, CatSite, Volume, CorrFiltre, N) %>% 
  mutate(N = ifelse(N>=1,1,0)) %>% 
  left_join(LacInv1 %>% select(NomLac, Presence))

head(DATA)
  
  
DATA$CatSite <- factor(DATA$CatSite)
DATA$CorrFiltre <- factor(DATA$CorrFiltre)
DATA$Volume <- factor(DATA$Volume)
DATA$Presence <- factor(DATA$Presence)

library(lme4)    

mod1 <- glmer(N ~ CatSite + Volume + CorrFiltre + Presence + NameAssign + 
                  NameAssign:CatSite + NameAssign:Presence +
                   (1|NomLac), 
              family = "binomial", data = DATA,
              control=glmerControl(optimizer="bobyqa"))

mod1 <- glmer(N ~ CatSite + (1|NomLac), family = "binomial", data = DATA)  

summary(mod1)    
    
library(effects)

plot(allEffects(mod1))


# Check blast 99 ----------------------------------------------------------


ASVtab.12s.byID.SP %>% group_by(ID, Assign) %>% 
                       summarise(N = sum(N)) %>% 
                       left_join(Blast99) %>% View()


ASVtab.12s.cor.wTAXO %>% View()


head(Blast99)

