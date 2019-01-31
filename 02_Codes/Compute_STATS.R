
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


LacPeche <- read_excel(get.value("Lac.xl"),sheet="StatPeche",na="NA",guess_max=100000)
LacPeche <- LacPeche %>% gather(names(.) %>% str_subset("PUE"), key = "Mesure", value = "Value") 


LacInv1 %>% filter(Presence == 0) %>% pull(Espece) %>% unique()

ggarrange(LacPeche %>% ggplot(aes(x = CPUEjp, y = CPUEhm, col = Espece))+ geom_point(size = 2) +theme_classic(),
          LacPeche %>% ggplot(aes(x = BPUEjp, y = BPUEhm, col = Espece))+ geom_point(size = 2) +theme_classic(),
          LacPeche %>% ggplot(aes(x = CPUEjp, y = BPUEjp, col = Espece))+ geom_point(size = 2) +theme_classic(),
          LacPeche %>% ggplot(aes(x = CPUEhm, y = BPUEhm, col = Espece))+ geom_point(size = 2) +theme_classic(),
          ncol =2 , nrow =2, common.legend = T, legend = "right")

# Mock community info

Mock.dat <- read_excel(get.value("Sample.xl"),sheet="Mock",na="NA",guess_max=100000)

Mock.dat <- Mock.dat %>% gather(paste0("Mix",1:6), key="Mix", value = "Vol") %>% 
                         mutate(DNA = Vol * 10) # DNA concentration = 10 ng/ul

Mock.dat


# Ref sequences

# Fichier reference taxo - la boucle permet de verifier si le fichier a bien ete lu
REF <- read_csv(get.value("RefTAXO"))
if(ncol(REF) == 1 ) {
  REF <- read_csv2(get.value("RefTAXO"))
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


View(reads.tab)

#save(file = get.value("BasicStats.data"), 
#     list = "reads.tab")

names(ASVtab.12s)



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


for(x in ls() %>% str_subset(".wTAXO") %>% str_remove(".cor") %>% unique()){

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


# Abundance


Mock.abund.data <- Mock.graph.data %>% mutate(Method = str_sub(Data,1,3),
                                              Locus = Data %>% str_remove(paste0(Method, "tab."))) %>%  
                                       filter(Locus %in%  c("12s", "cytB.R1"),
                                              Name.level %in% c("Salvelinus fontinalis",
                                                                "Salvelinus", 
                                                                "Micropterus dolomieu")) %>%
                                       mutate(Species = ifelse(Name.level == "Salvelinus", "Salvelinus fontinalis", Name.level)) %>% 
                                       left_join(Mock.final)
                                                 #by = c("Name.level" = "Species", "Mix" = "Mix"))

Mock.abund.data %>% filter(!(Locus == "cytB.R1" & Name.level =="Salvelinus" )) %>% 
                    group_by(Locus, Name.level, Method) %>% 
                    summarise(rho= cor.test(N, DNAfinal, method = "spearman")$estimate,
                              p.value= cor.test(N, DNAfinal, method = "spearman")$p.value,
                              N = length(N))

Mock.abund.data %>% filter(!(Locus == "cytB.R1" & Name.level =="Salvelinus" )) %>%
                    mutate(Ncor = N + 1) %>% 
                     ggplot(aes(x = DNAfinal, y = Ncor, col = Method, shape = Method))+
                           geom_smooth(method = "lm", se = T , lty = "dashed", size = 0.5, fill = "gray85")  +
                           geom_jitter(width = 0.05, height = 0.05, size = 2)+
                           scale_x_continuous(trans = "log10")+
                           scale_y_continuous(trans = "log10")+
  labs(title= NULL, x ="Concentration d'ADN (ng/ul)", y = "N lectures") +                         
  facet_grid(Locus ~ Species) +
                           theme_bw()










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
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

}

ggarrange(
haplo.mock(ASVtab.12s.byID.SP),
haplo.mock(OTUtab.12s.byID.SP),
haplo.mock(ASVtab.cytB.R1.byID.SP),
haplo.mock(OTUtab.cytB.R2.byID.SP),
ncol = 2, nrow = 2)

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


# Samples -----------------------------------------------------------------

Sample.data <- data.frame(Assign = character(), 
                          Sample  = character(),
                          N = integer(),
                          Level = integer(),
                          NameAssign = character(),
                          Cat = character(),
                          Nsite = integer(),
                          NomLac = character(),
                          Data = character())

for(x in ls() %>% str_subset(".cor.wTAXO")){

  print(x)
  
DATA <-
  
  get(x) %>% gather(names(.) %>% str_subset("Sample"), key = Sample, value = N) %>% 
        group_by(Assign, Sample) %>% 
        summarise(N = sum(N)) %>% 
        
        mutate(Level = str_count(Assign, ";"),
                                 NameAssign = NA,
          Lac = sapply(str_split(Sample, "_"), `[`, 3),
                                                 Cat = sapply(str_split(Sample, "_"), `[`, 4),
                                                Puit = sapply(str_split(Sample, "_"), `[`, 5) %>% str_remove("p")) %>% 
          left_join(DataSeq %>% select(IbisID, Nsite, CatSite, SeqType, NomLac, CorrFiltre, Volume), by = c("Puit" = "IbisID")) %>% 
          filter(SeqType == "sample") %>% 
        mutate(Data = x %>% str_remove(".cor.wTAXO"),
               Cat = CatSite)   
  
  
for(y in 1:nrow(DATA)){
  DATA$NameAssign[y] <- sapply(str_split(DATA[y,"Assign"], ";"),`[`, pull(DATA[y,"Level"] + 1))
  
  }

DATA <- DATA %>% group_by() %>% select(-c(Lac, Puit, SeqType)) %>% 
  as.data.frame()  
 
Sample.data <- rbind(Sample.data, DATA)
 
}  


Sample.data <- Sample.data %>% mutate(Method = str_sub(Data,1,3),
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

head(Sample.data)

# Graph of what was done

Sample.data %>% filter(str_detect(Assign, "Teleostei") == T,
                   str_detect(Assign, "Gadidae") == F,
                   str_detect(Assign, "Sebaste") == F,
                   Locus %in% c("12s"),
                   Method %in% c("ASV")) %>%
  group_by(NomLac, Locus, Method, NameAssign, Taxo, Location) %>% 
  summarise(Nmed = median(N),
            Nmax = max(N)) %>%
  
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

# Comparison trad vs eDNA

head(Sample.data)

CompPAtrad <- Sample.data %>% mutate(NameAssign = ifelse(NameAssign == "Salvelinus", "Salvelinus fontinalis", NameAssign)) %>% 
  filter(str_detect(Assign, "Teleostei"),
         str_detect(NameAssign, " "),
         str_detect(NameAssign, "Gadus morhua") == FALSE,
         str_detect(NameAssign, "Salmo salar") == FALSE) %>% 
  group_by(NameAssign, Location, NomLac, Method, Locus) %>% 
  summarise(N = sum(N)) %>% 
  mutate(PresenceADNe = ifelse(N >=1, 1, 0)) %>% 
  full_join(LacInv1 %>% select(NomLac, Location, Espece, Presence),
            by = c("NomLac" = "NomLac", "Location" = "Location", "NameAssign" = "Espece")) %>% 
  mutate(PresenceADNe = ifelse(is.na (PresenceADNe), 0, PresenceADNe),
         DiffInv = ifelse(PresenceADNe == Presence, 
                          ifelse(PresenceADNe == 0 , "Absent trad et ADNe", "Present trad et ADNe"),
                   ifelse(PresenceADNe == 0, "Present trad seul", "Present ADNe seul"))) %>% 
 left_join(REF %>% select(Espece_initial, Class, Order, Family, NomFR),
                         by = c("NameAssign" = "Espece_initial")) 


  #group_by(Method, Locus, DiffInv) %>% summarise(N = length(NomLac)) %>% View()
  #View()
CompPAtrad %>%   filter(Location == "Avant-pays",
         Method %in% c("ASV", NA),
         Locus %in% c("12s", NA)) %>% #View()
  ggplot(aes(x=NomLac, y = NomFR, fill = DiffInv)) +
  geom_bin2d(col = "black") +
  labs(title= NULL, x ="Lac", y = "Assignation") +
  #scale_y_discrete(position = "right") +
  scale_fill_manual(limits = c("Present trad et ADNe", "Present trad seul", "Present ADNe seul", "Absent trad et ADNe"),  
                    values = c("green", "red", "yellow", "white")) + 
  guides(fill = guide_legend(title = NULL)) + 
  theme_bw ()+
  facet_grid(Order + Family ~ ., 
             scale = "free", space = "free", 
             labeller = label_value,
             switch = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0),
        panel.spacing = unit(0, "lines"),
        strip.placement = "outside",
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.text.y = element_text(colour="black", hjust = 1)) #+ coord_flip()


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

# Idem mais pour les salvelidées

Sample.data %>% filter(str_detect(NameAssign, "Salvelinus"),
                       Location == "Avant-pays") %>% 
  group_by(NameAssign, Location, NomLac, Method, Locus) %>% 
  summarise(N = sum(N)) %>% 
  mutate(PresenceADNe = ifelse(N >=1, 1, 0)) %>% 
  select(-c(N)) %>% 
  spread(NameAssign, value = PresenceADNe, fill = 0) %>%
  rename(Salvelinus = "Salvelinus.ADNe", "Salvelinus fontinalis" = "Salvelinus fontinalis.ADNe") %>% 
  left_join(LacInv1 %>% filter(str_detect(Espece, "Salvelinus")) %>% 
                        spread(Espece, Presence) %>% 
                        select(-c(InvLac, Rotenode, Location)),
                        by = c("NomLac" = "NomLac")) %>% 
  gather(names(.) %>% str_subset("Salvelinus"), key = Salvelinus, value = Presence) %>% 
  mutate(Salvelinus = factor(Salvelinus, levels = c("Salvelinus fontinalis", "Salvelinus namaycush", "Salvelinus alpinus", "Salvelinus fontinalis.ADNe", "Salvelinus.ADNe"))) %>% 
  filter(Locus %in% c("12s", "cytB.R1"),
       Method %in% c("ASV", "OTU"),
       Location == "Avant-pays") %>% 
  ggplot(aes(y=NomLac, x = Salvelinus, fill = factor(Presence))) +
  geom_bin2d(col = "black") +
  labs(title= NULL, x ="", y = "Lac") +
  scale_fill_manual(limits = c("0", "1"),  
                    values = c("white","green")) + 
  theme_bw()+
  facet_grid(Locus ~ Method, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0))
  
LacInv1 %>% head()  


# ## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade") ## you may need this
biocLite("ggtree")



  pull(NameAssign) %>% unique()tion

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


# Stats data

Sample.data %>% filter(Data == "ASVtab.12s") %>% 
                group_by(Location, NomLac, Cat) %>% 
                summarise(Ntot = length(unique(Sample))) %>% 
                spread(Cat, Ntot, fill = 0) %>% 
                mutate(N = PEL + RIV) %>% View()





# Correlation with inventaire -----------------------------------------------------

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


Sample.data.cor <- Sample.data %>% mutate(Nlog = ifelse(N == 0 , 0, log2(N)),
                       Nlog.cor = Nlog / CorrFiltre / Volume * 1000) %>% 
                  gather(Nlog, Nlog.cor, key = Correction, value = Ncor)

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
                                        Correction = c("Nlog",  "Nlog.cor"),
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


RES.COR.LOG.INV2.graph %>% spread(Correction, Value.cor) %>% ggplot(aes(x = Nlog, y = Nlog.cor)) + geom_point()

RES.COR.LOG.INV2.graph %>% spread(Correction, Value.cor) %>%  
  mutate(diff = Nlog - Nlog.cor) %>% #head()
  group_by(Peche, Mesure, Cat) %>% summarise(Mean = mean(diff, na.rm=T))
  ggplot(aes(x = diff)) + geom_histogram()

  
Sample.data.cor %>% filter(Correction == "Nlog.cor",
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

