# Info --------------------------------------------------------------------

# Change parameters if necessery
# 
# Audrey Bourret
# 2018-11-14
#

# Library -----------------------------------------------------------------

#  
source(file.path("./03_Functions",  list.files("./03_Functions")))


PARAM <- data.frame(x = character(), value = character())

# Specific functions ------------------------------------------------------


"add.param<-" <- function(PARAM, ..., value){
  PARAM <- rbind(PARAM, value)
  names(PARAM) <- c("x", "value")
  PARAM$x     <- as.character(PARAM$x)
  PARAM$value <- as.character(PARAM$value)
  PARAM

}


# Path --------------------------------------------------------------------

add.param(PARAM) <- c("info.path", "./00_Data/00_FileInfos")
add.param(PARAM) <- c("ref.path", "./00_Data/03_RefSeq")
add.param(PARAM) <- c("result.path", "./01_Results/01_data")
add.param(PARAM) <- c("log.path", "./04_Log")


# Files -------------------------------------------------------------------

add.param(PARAM) <- c("Sample.xl", "DB_Echantillons.xlsx")




# Save Parameters ---------------------------------------------------------

PARAM

write.csv2(PARAM, file = file.path("./04_Log", "Options.csv"), row.names=F)

get.value("result.path")



