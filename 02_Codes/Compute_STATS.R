
# Info --------------------------------------------------------------------

# Compute stats for this dataset
# 
# Audrey Bourret
# 2019-01-14 
#

# Library -----------------------------------------------------------------


# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}



# Data --------------------------------------------------------------------

load(get.value("ALLtable.data"))


ASVtab.12s.cor.bySP %>% group_by(Assign, Level) %>% 
                        summarise(Nread = sum(N),
                                  Nsample = length(which(N > 0))) %>% 
                        arrange(Assign) %>% View()



# Graphe N read par niveau par dataset
