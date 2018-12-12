# Info --------------------------------------------------------------------

# Change parameters if necessery
# 
# Audrey Bourret
# 2018-11-14
#

# Library -----------------------------------------------------------------

# Internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}


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

# Data
add.param(PARAM) <- c("info.path", "./00_Data/00_FileInfos")

add.param(PARAM) <- c("raw.path", "./00_Data/01a_RawData")
add.param(PARAM) <- c("raw_unz.path", "./00_Data/01b_RawData_unzipped")
add.param(PARAM) <- c("raw_unz_rename.path", "./00_Data/01c_RawData_unzipped_rename")

add.param(PARAM) <- c("filt_cutadapt.path", "./00_Data/02a_Cutadapt")
add.param(PARAM) <- c("filt_dada2.path", "./00_Data/02b_Filtered_dada2")

add.param(PARAM) <- c("ASV.dada2.path", "./00_Data/03_ASV")

add.param(PARAM) <- c("filt_merge.path", "./00_Data/04a_Merged_usearch")
add.param(PARAM) <- c("filt_derep.path", "./00_Data/04b_Derep_vsearch")
add.param(PARAM) <- c("filt_derep.1file.path", "./00_Data/04c_Derep_united_vsearch")
add.param(PARAM) <- c("OTU.usearch", "./00_Data/04d_OTU")
add.param(PARAM) <- c("Compare.OTU.usearch", "./00_Data/04e_Compared_OTU")

#add.param(PARAM) <- c("filt_maxEE.usearch.path", "./00_Data/04b_maxEE_usearch")

                      
#add.param(PARAM) <- c("filt_JAMP.path", "./00_Data/02c_Filtered_JAMP")

add.param(PARAM) <- c("ref.path", "./00_Data/05_RefSeq")



# Results

add.param(PARAM) <- c("result.path", "./01_Results")
add.param(PARAM) <- c("result.data.path", "./01_Results/01_data")
add.param(PARAM) <- c("result.Q.path", "./01_Results/02_ReadsQ")
add.param(PARAM) <- c("result.FQraw.path", "./01_Results/02_ReadsQ/FastQC_Raw")
add.param(PARAM) <- c("log.path", "./04_Log")




# Update file before continue to ensure next section will do OK!
write.csv2(PARAM, file = file.path("./04_Log", "Options.csv"), row.names=F)

# Files -------------------------------------------------------------------

add.param(PARAM) <- c("Sample.xl", paste(get.value("info.path"),"DB_Echantillons.xlsx", sep = "/"))

add.param(PARAM) <- c("Raw.log", paste(get.value("log.path"),"Process_RAW.log.txt", sep = "/"))



# Data --------------------------------------------------------------------

add.param(PARAM) <- c("Qplot.RAW.data", paste(get.value("result.data.path"),"Qplot.RAW.data", sep = "/"))
add.param(PARAM) <- c("FastQC.data", paste(get.value("result.data.path"),"FastQC.data", sep = "/"))
add.param(PARAM) <- c("Qplot.FILT.data", paste(get.value("result.data.path"),"Qplot.FILT.data", sep = "/"))






# Save Parameters ---------------------------------------------------------

PARAM

write.csv2(PARAM, file = file.path("./04_Log", "Options.csv"), row.names=F)

get.value("result.path")



