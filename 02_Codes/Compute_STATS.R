
# Info --------------------------------------------------------------------

# Compute stats for this dataset
# 
# Audrey Bourret
# 2019-01-14 
#

# Library -----------------------------------------------------------------

library(tidyverse)


# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}



# Data --------------------------------------------------------------------

load(get.value("ALLtable.data"))


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
                            str_detect(.$Assign, "Root;Chordata;Teleostei;") == T) %>% 
                     mutate(Assign = str_remove(Assign, "Root;Chordata;Teleostei;"),
                            Method = str_sub(Data, 1, 3),
                            Locus = Data %>% str_remove(paste0(Method,"tab."))) %>% 
                     arrange(Assign, Data) %>% 
                     ggplot(aes(y= Nread, x = Assign, fill = Method)) +
                              geom_bar(stat = "identity", position = "dodge") +
                              facet_grid(Locus ~ .) +
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
  
                group_by() %>%  summarise(max = max(Nread)) %>% 
                pull()
               
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


# Mock community ----------------------------------------------------------

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


library(gtools)    # for mixedsort

#library(devtools)
#devtools::install_github("kassambara/ggpubr")
library(ggpubr)    # on github - for nice graphs




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







# Duplicated samples ------------------------------------------------------





# N haplotypes ------------------------------------------------------------



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




