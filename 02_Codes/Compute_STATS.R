
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
  left_join(DataSample %>% select(SampleID, NomLac, CatSite, Nsite), by = "SampleID")

# Lac info

LacSample <- read_excel(get.value("Lac.xl"),sheet="NomLac",na="NA",guess_max=100000)
LacSample

LacInv1 <- read_excel(get.value("Lac.xl"),sheet="Inventaire1996-2003",na="NA",guess_max=100000)
LacInv1 <- na.replace(LacInv1, 0)
LacInv1

LacInv2 <- read_excel(get.value("Lac.xl"),sheet="Inventaire2018",na="NA",guess_max=100000)
LacInv2

LacPeche <- read_excel(get.value("Lac.xl"),sheet="StatPeche",na="NA",guess_max=100000)
LacPeche

LacInv1 <- LacInv1 %>% inner_join(LacSample %>% select(-c(AbrLac)), by = "InvLac") %>% 
                       gather(names(.) %>% str_subset(" "), key= "Espece", value = "Presence")

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
                     mix = mix +1) %>%  
    ggplot(aes(x = mix, y = dup.mix, col = Mix)) + 
    # Add a 1:1 line
    geom_abline(intercept = 0, colour = "gray", slope = 1, linetype = 1, show.legend = FALSE ) +
    geom_jitter(width = 0.05, height = 0.05) +
    scale_x_continuous(name = "N reads - sample 1", trans = "log10") +
    scale_y_continuous(name = "N reads - sample 2",trans = "log10") +
    guides(col = guide_legend(title = NULL, nrow = 3)) +
    annotate("text" , x = 10, y = 10000, 
             label = paste("rho =", round(COR.RES$estimate, 3), ";", 
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

GRAPH

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
    # Limit to MIX samples
    filter(str_detect(Sample, "Mix")) %>% 
    mutate(Mix = sapply(str_split(Sample, "_"), `[`, 2),
           Puit = sapply(str_split(Sample, "p[:digit:]."), `[`, 2) %>% str_remove("_R[:digit:]")) %>% 
    left_join(DataSeq %>% filter(SeqType %in% c("mix", "dup.mix")) %>% select(Puit, SeqType)) %>% 
    # Remove technical duplicates
    filter(SeqType == "mix") %>% 
    select(Assign, N, Method, Mix) %>% 
    spread(key = Method, value = N, fill = 0) %>% 
    mutate(KEEP = ASV + OTU) %>% 
    filter(KEEP > 0) %>% select(-KEEP)
  
  Comp.res
  
  COR.RES <- cor.test(Comp.res$ASV, Comp.res$OTU, method = "spearman")
  
  Comp.res %>% mutate(ASV = ASV + 1,
                      OTU = OTU +1) %>%  
    ggplot(aes(x = ASV, y = OTU, col = Mix)) + 
    # Add a 1:1 line
    geom_abline(intercept = 0, colour = "gray", slope = 1, linetype = 1, show.legend = FALSE ) +
        geom_jitter(width = 0.05, height = 0.05) +
    scale_x_continuous(name = "N reads - ASV", trans = "log10") +
    scale_y_continuous(name = "N reads - OTU",trans = "log10") +
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

GRAPH <- ggarrange(ASVOTU.graph(ASVtab.12s.bySP, OTUtab.12s.bySP) + guides(col=FALSE) + ggtitle("12S") ,
                   ASVOTU.graph(ASVtab.cytB.R1.bySP, OTUtab.cytB.R1.bySP) + ggtitle("cytB (R1)"), 
                   
                   labels = LETTERS[1:2],
                   ncol = 2, nrow=1
                   #common.legend = T, legend = "bottom"
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
                           Taxo = ifelse(Level == 2, "Ordre", 
                                               ifelse(Level == 3, "Classe", 
                                               ifelse(Level == 4, "Famille",
                                               ifelse(Level == 5, "Genre", 
                                               ifelse(Level == 6, "Espece", NA))))),
                           Taxo = factor(Taxo, levels = c("Ordre", "Classe", "Famille", "Genre", "Espece" ))
                           ) %>%
  filter(str_detect(Assign, "Teleostei") == T,
         SeqType == "mix",
         Locus %in% c("12s", "cytB.R1")) %>%
  filter(N>=1) %>% 
  
  ggplot(aes(x = Mix, y = Name.level, fill = N)) + 
  geom_bin2d() + 
  scale_fill_distiller(palette = "Spectral", trans = "log10") +
  #scale_fill_gradient(low = "darkgray", high = "red", trans = "log") +
  #scale_y_discrete(limits=mixedsort(tab2$Assign)) + #, labels = NULL) +
  labs(title= NULL, x ="Sample", y = "Assigment") +
  guides(fill = guide_colourbar(title = "N reads", title.hjust = 0)) + 
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
                           geom_smooth(method = "loess", se = T , lty = "dashed", size = 0.5, fill = "gray85")  +
                           geom_jitter(width = 0.05, height = 0.05, size = 2)+
                           scale_x_continuous(trans = "log10")+
                           scale_y_continuous(trans = "log10")+
                           facet_grid(Species~Locus) +
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
                          Name.level = character(),
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
                           Name.level = NA,
    Lac = sapply(str_split(Sample, "_"), `[`, 3),
                                           Cat = sapply(str_split(Sample, "_"), `[`, 4),
                                           Puit = sapply(str_split(Sample, "_"), `[`, 5) %>%
      

      str_remove("p")) %>% 
                                    left_join(DataSeq %>% select(IbisID, Nsite, SeqType, NomLac), by = c("Puit" = "IbisID")) %>% 
    filter(SeqType == "sample") %>% mutate(Data = x %>% str_remove(".cor.wTAXO"))   
  
  
for(y in 1:nrow(DATA)){
  DATA$Name.level[y] <- sapply(str_split(DATA[y,"Assign"], ";"),`[`, pull(DATA[y,"Level"] + 1))
  
  }

DATA <- DATA %>% group_by() %>% select(-c(Lac, Puit, SeqType)) %>% 
  as.data.frame()  
 
Sample.data <- rbind(Sample.data, DATA)
 
}  


Sample.data %>% mutate(Method = str_sub(Data,1,3),
                       Locus = Data %>% str_remove(paste0(Method, "tab.")),
                       Taxo = ifelse(Level == 2, "Ordre", 
                                               ifelse(Level == 3, "Classe", 
                                               ifelse(Level == 4, "Famille",
                                               ifelse(Level == 5, "Genre", 
                                               ifelse(Level == 6, "Espece", NA))))),
                           Taxo = factor(Taxo, levels = c("Ordre", "Classe", "Famille", "Genre", "Espece" ))
                           ) %>%
  filter(str_detect(Assign, "Teleostei") == T,
         #SeqType == "mix",
         Locus %in% c("12s"),
         Method %in% c("ASV")) %>%#head()
  group_by(Name.level, NomLac, Method, Locus, Taxo) %>% 
    summarise(N = sum(N)) %>% 
  #complete(Name.level, Cat, Nsite, NomLac, Method, Locus, Taxo, fill = list(0)) %>%  #View()
#  filter(NomLac == "Alphonse", Name.level == "Culaea") %>% 
  
  filter(N>=1) %>% 
  
  ggplot(aes(x = NomLac, y = Name.level, fill = N)) + 
  geom_bin2d() + 
  scale_fill_distiller(palette = "Spectral", trans = "log10") +
  #scale_fill_gradient(low = "darkgray", high = "red", trans = "log") +
  #scale_y_discrete(limits=mixedsort(tab2$Assign)) + #, labels = NULL) +
  labs(title= NULL, x ="Sample", y = "Assigment") +
  guides(fill = guide_colourbar(title = "N reads", title.hjust = 0)) + 
  theme_bw()+
  facet_grid(Taxo~Locus + Method, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0)) #+ coord_flip()
                                    


assign.graph(tab = ASVtab.12s.wTAXO, Sample = T, Tneg = F, Mix = F, 
             maintitle = "12s ASVs - Samples", 
             Nlevel = c(2:6), 
             subAssign = "Root;Chordata;Teleostei;")




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
  
  




