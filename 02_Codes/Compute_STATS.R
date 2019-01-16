
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

DataSeq    <- read_excel(get.value("Sample.xl"),sheet="DataSeq",na="NA",guess_max=100000)
DataSeq

# Mock community info

Mock.dat <- read_excel(get.value("Sample.xl"),sheet="Mock",na="NA",guess_max=100000)

Mock.dat <- Mock.dat %>% gather(paste0("Mix",1:6), key="Mix", value = "Vol") %>% 
                         mutate(DNA = Vol * 10) # DNA concentration = 10 ng/ul

Mock.data

# Summary by SP

load(get.value("ALLtable.data"))

ls() %>% str_sub("tab")

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

ls() %>% str_sub("byID.SP")


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


Mock.res <- ASVtab.12s.wTAXO %>% mutate(Level = str_count(Assign, ";")) %>% 
                     filter(Level %in% c(5:6),
                            str_detect(.$Assign, "Teleostei") == T) %>% 
                     mutate(Assign = str_remove(Assign, "Root;Chordata;Teleostei;")) %>% 
                     select(Assign, Level, names(.) %>% str_subset("Mix") ) %>%
                     gather(names(.) %>% str_subset("Mix"), key = "sample", value = "Nread") %>% 
                     mutate(Mix = sapply(str_split(sample, "_"), `[`, 2),
                            Puit = sapply(str_split(sample, "p[:digit:]."), `[`, 2)) %>% 
                     mutate(Assign = str_replace(Assign, ";Salvelinus", ";Salvelinus;Salvelinus fontinalis"),  
                            Level = str_count(Assign, ";")) %>%     
                     filter(Level == 3) %>% 
                     mutate(Species = sapply(str_split(Assign, ";"), `[`, 4) ) %>% 
                     group_by(Species, Mix, Puit) %>%
                     summarise(Nread = sum(Nread)) %>% 
                     left_join(Mock.final) %>% 
                     filter(!is.na(DNAfinal)) %>% 
                     mutate(Nread.log = log10(Nread + 1)) %>% 
                     left_join(DataSeq %>% filter(SeqType %in% c("mix", "dup.mix")) %>% select(Puit, SeqType)) %>% 
                     select(-Puit)
                     
                   
cor.test(Mock.res %>% filter(Species == "Salvelinus fontinalis") %>% select(Nread) %>% pull(),
         Mock.res %>% filter(Species == "Salvelinus fontinalis") %>% select(DNAfinal) %>% pull()
         )

cor.test(Mock.res %>% filter(Species == "Micropterus dolomieu") %>% select(Nread.log) %>% pull(),
         Mock.res %>% filter(Species == "Micropterus dolomieu") %>% select(DNAfinal) %>% pull(),
         method = "spearman")

cor.test(Mock.res %>% select(Nread) %>% pull(),
         Mock.res %>% select(DNAfinal) %>% pull(),
         method = "spearman")

par(mfrow = c(1,1))

plot(Mock.res %>% select(Nread.log) %>% pull(),
         Mock.res %>% select(DNAfinal) %>% pull())

Mock.res %>% filter(Species == "Salvelinus fontinalis")

Mock.res %>% ggplot(aes(x= DNAfinal, y = Nread.log, col = Species)) + 
  geom_jitter() + geom_smooth(method = 'loess') + facet_wrap(~Species, nrow = 2) + theme_bw()


Mock.res %>% ggplot(aes(x= DNAfinal, y = Nread.log)) + 
  geom_jitter() + geom_smooth(method = 'loess') #+ facet_wrap(~Species, nrow = 2) + theme_bw()


Mock.res %>% ggplot(aes(x= Nread)) + 
  geom_histogram()

Mock.res$Species <- as.factor(Mock.res$Species)

library(lme4)
library(lmerTest)
library(effects)

mod1 <- lmer(Nread.log ~ DNAfinal + Species + (1|Mix), data = Mock.res)

plot(allEffects(mod1))

summary(mod1)


E1 <- resid(mod1)
F1 <- fitted(mod1)
plot(x = F1, 
     y = E1, 
     xlab = "Fitted Values",
     ylab = "Normalized residuals")
abline(h = 0, lty = 2)



# Function to create a heatmap for each sample from a SEQtab
assign.graph <- function(tab, 
                         Sample = T, Tneg = T, Mix = T, 
                         maintitle = "Heatmap of species assigment", 
                         Nlevel = c(0:6), 
                         subAssign = "Root;"){
  
  Sample1 <- vector()
  Tneg1   <- vector()
  Mix1    <- vector()  
  
  Sample1 <- if(isTRUE(Sample)) names(tab) %>% str_subset("Sample")
  
  Tneg1   <- if(isTRUE(Tneg)) {
    c(names(tab) %>% str_subset("Tneg"),
      names(tab) %>% str_subset("T0"),
      names(tab) %>% str_subset("T1"))
  }
  
   Mix1    <- if(isTRUE(Mix)) names(tab) %>% str_subset("Mix")
 
  tab2 <- tab %>% mutate(Level = str_count(Assign, ";")) %>% 
                 filter(Level %in% Nlevel,  
                 str_detect(.$Assign, subAssign) == T) %>% 
                 mutate(Assign = str_remove(Assign, subAssign)) %>% 
    select(Assign,Sample1, Tneg1, Mix1) %>%
    gather(Sample1, Tneg1, Mix1, key = "sample", value = "N") %>% 
    filter(N>=1)  %>%
    group_by(Assign, sample) %>% 
    summarise(N = sum(N)) %>% 
    arrange(desc(Assign))
  
  graph <- tab2 %>% 
    ggplot(aes(x = sample, y = Assign, fill = N)) + 
    geom_bin2d() + 
    scale_fill_distiller(palette = "Spectral", trans = "log10") +
    #scale_fill_gradient(low = "darkgray", high = "red", trans = "log") +
    #scale_y_discrete(limits=mixedsort(tab2$Assign)) + #, labels = NULL) +
    labs(title= maintitle, x ="Sample", y = "Assigment") +
    guides(fill = guide_colourbar(title = "N reads", title.hjust = 0)) + 
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.ticks.y = element_blank())
  
  print(graph)
  
}


GRAPH1 <- assign.graph(tab = ASVtab.12s.wTAXO, Sample = F, Tneg = F, Mix = T, 
             maintitle = "12s ASVs - Mock community", 
             Nlevel = c(1:6), 
             subAssign = "Root;Chordata;")

GRAPH2 <- assign.graph(tab = OTUtab.12s.wTAXO, Sample = F, Tneg = F, Mix = T, 
             maintitle = "12s OTUs - Mock community", 
             Nlevel = c(1:6), 
             subAssign = "Root;Chordata;")



GRAPH3 <- assign.graph(tab = ASVtab.cytB.R1.wTAXO, Sample = F, Tneg = F, Mix = T, 
                       maintitle = "cytB (R1) ASVs - Mock community", 
                       Nlevel = c(1:6), 
                       subAssign = "Root;Chordata;")

GRAPH4 <- assign.graph(tab = OTUtab.cytB.R1.wTAXO, Sample = F, Tneg = F, Mix = T, 
                       maintitle = "cytB (R1) OTUs - Mock community", 
                       Nlevel = c(1:6), 
                       subAssign = "Root;Chordata;")

GRAPH5 <- assign.graph(tab = ASVtab.cytB.R2.wTAXO, Sample = F, Tneg = F, Mix = T, 
                       maintitle = "cytB (R2) ASVs - Mock community", 
                       Nlevel = c(1:6), 
                       subAssign = "Root;Chordata;")

GRAPH6 <- assign.graph(tab = OTUtab.cytB.R2.wTAXO, Sample = F, Tneg = F, Mix = T, 
                       maintitle = "cytB (R2) OTUs - Mock community", 
                       Nlevel = c(1:6), 
                       subAssign = "Root;Chordata;")


ggarrange(GRAPH1, GRAPH2,
          labels = LETTERS[1:2],
          ncol = 1, nrow = 2,
          common.legend = TRUE, legend = "right")

ggarrange(GRAPH3, GRAPH4,
          labels = LETTERS[3:4],
          ncol = 1, nrow = 2,
          common.legend = TRUE, legend = "right")

ggarrange(GRAPH5, GRAPH6,
          labels = LETTERS[5:6],
          ncol = 1, nrow = 2,
          common.legend = TRUE, legend = "right")

pdf(file.path(get.value("result.FINAL"),"MockCommunity_heatmap.pdf"), width = 8, height = 6)

  GRAPH1
  GRAPH2
  GRAPH3
  GRAPH4
  GRAPH5
  GRAPH6

dev.off()






### Mes échantillons

assign.graph(tab = ASVtab.12s.wTAXO, Sample = T, Tneg = F, Mix = F, 
             maintitle = "12s ASVs - Samples", 
             Nlevel = c(2:6), 
             subAssign = "Root;Chordata;Teleostei;")










# N haplotypes ------------------------------------------------------------


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



# Summary -----------------------------------------------------------------



Summary.bySP.cor <- data.frame(Assign = character(), Level = integer(), Nread = numeric(), Nsample = integer(), Data = character())


for(x in ls() %>% str_subset("cor.bySP")){
  DATA <- get(x) %>% group_by(Assign, Level) %>% 
    summarise(Nread = sum(N),
              Nsample = length(which(N > 0))) %>% 
    mutate(Data = x %>% str_remove(".cor.bySP")) %>% 
    as.data.frame()
  
  
  Summary.bySP.cor <- rbind(Summary.bySP.cor, DATA)
  
}


Summary.bySP.cor %>% filter(Level > 2, 
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



Summary.bySP.cor %>% #filter(Level > 2, 
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

TEST <- Summary.bySP.cor %>%  mutate(Cat = ifelse(str_detect(.$Assign, "Teleostei") == T,
                                                  "Teleostei",
                                                  "Other"),
                                     Method = str_sub(Data, 1, 3),
                                     Locus = Data %>% str_remove(paste0(Method,"tab."))) %>% 
  group_by(Level, Data, Locus, Method) %>% 
  summarise(Nread = sum(Nread),
            Nsample = sum(Nsample)) 

Nread.scale <- TEST %>% group_by(Data) %>% 
  summarise(Nread = sum(Nread)) %>%
  mutate(Scale = Nread / max(.$Nread)) %>% 
  select(Data, Scale)


TEST %>% left_join(Nread.scale ) %>% 
  ggplot(aes(x=Scale/2, y = Nread, fill = Level, width = Scale)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_polar("y", start=0) +
  guides(fill = guide_legend(ncol = 1, title = "Depth")) +
  labs(title= "Proportion of reads assigned to each taxonomic depth",
       subtitle = "From no assignment (0) to species (6)")+
  ylab(NULL) + xlab(NULL)+
  facet_grid(Method ~ Locus) + 
  theme_bw() 



