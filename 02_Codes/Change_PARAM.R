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
add.param(PARAM) <- c("filt_cutadapt.log", "./00_Data/02a_Cutadapt/log")

add.param(PARAM) <- c("filt_dada2.path", "./00_Data/02b_Filtered_dada2")

add.param(PARAM) <- c("ASV.dada2.path", "./00_Data/03_ASV")

add.param(PARAM) <- c("filt_merge.path", "./00_Data/04a_Merged_usearch")
add.param(PARAM) <- c("filt_merge.log", "./00_Data/04a_Merged_usearch/log")
add.param(PARAM) <- c("filt_min.path", "./00_Data/04b_Min_cutadapt")
add.param(PARAM) <- c("filt_min.log", "./00_Data/04b_Min_cutadapt/log")
add.param(PARAM) <- c("filt_derep.path", "./00_Data/04c_Derep_vsearch")
add.param(PARAM) <- c("filt_derep.1file.path", "./00_Data/04d_Derep_united_vsearch")
add.param(PARAM) <- c("OTU.usearch", "./00_Data/04e_OTU")
add.param(PARAM) <- c("Compare.OTU.usearch", "./00_Data/04f_Compared_OTU")
add.param(PARAM) <- c("Compare.OTU.usearch.log", "./00_Data/04f_Compared_OTU/log")

#add.param(PARAM) <- c("filt_maxEE.usearch.path", "./00_Data/04b_maxEE_usearch")

                      
#add.param(PARAM) <- c("filt_JAMP.path", "./00_Data/02c_Filtered_JAMP")

add.param(PARAM) <- c("ref.path", "./00_Data/05_RefSeq")



add.param(PARAM) <- c("biodiv.path", "S:/Genpop/01-Projets de recherche/08-Banque reference biodiversite ESTL GSTL 2017-")
add.param(PARAM) <- c("path.EXTERNE", "S:/Genpop/01-Projets de recherche/08-Banque reference biodiversite ESTL GSTL 2017-/Sequences/Sequences_EXTERNE")
add.param(PARAM) <- c("path.LABO",  "S:/Genpop/01-Projets de recherche/08-Banque reference biodiversite ESTL GSTL 2017-/Sequences/Sequences_LABO_finales")


# Results

add.param(PARAM) <- c("result.path", "./01_Results")
add.param(PARAM) <- c("result.data.path", "./01_Results/01_data")
add.param(PARAM) <- c("result.Q.path", "./01_Results/02_ReadsQ")
add.param(PARAM) <- c("result.FQraw.path", "./01_Results/02_ReadsQ/FastQC_Raw")
add.param(PARAM) <- c("result.OTUtables", "./01_Results/03_OTUtables")
add.param(PARAM) <- c("result.ref", "./01_Results/04_Ref")

add.param(PARAM) <- c("result.FINAL", "./01_Results/05_Final")

add.param(PARAM) <- c("log.path", "./04_Log")




# Update file before continue to ensure next section will do OK!
write.csv2(PARAM, file = file.path("./04_Log", "Options.csv"), row.names=F)

# Files -------------------------------------------------------------------

add.param(PARAM) <- c("Sample.xl", paste(get.value("info.path"),"DB_Echantillons.xlsx", sep = "/"))
add.param(PARAM) <- c("Lac.xl", paste(get.value("info.path"),"DB_Lacs.xlsx", sep = "/"))

add.param(PARAM) <- c("RefTAXO", paste(get.value("biodiv.path"),"Reference_taxonomie.csv", sep = "/"))

add.param(PARAM) <- c("Raw.log", paste(get.value("log.path"),"Process_RAW.log.txt", sep = "/"))



# External programs -------------------------------------------------------


add.param(PARAM) <- c("blastn", "C:/Users/BourretA/Documents/Programs/ncbi-blast-2.8.1+/bin/blastn.exe")
add.param(PARAM) <- c("makeblastdb", "C:/Users/BourretA/Documents/Programs/ncbi-blast-2.8.1+/bin/makeblastdb.exe")



# Data --------------------------------------------------------------------

# RAW
add.param(PARAM) <- c("Qplot.RAW.data", paste(get.value("result.data.path"),"Qplot.RAW.data", sep = "/"))
add.param(PARAM) <- c("Qplot.FILT.data", paste(get.value("result.data.path"),"Qplot.FILT.data", sep = "/"))
add.param(PARAM) <- c("FastQC.data", paste(get.value("result.data.path"),"FastQC.data", sep = "/"))
add.param(PARAM) <- c("dada2.filt.data", paste(get.value("result.data.path"),"dada2.filt.data", sep = "/"))

add.param(PARAM) <- c("dada2.data", paste(get.value("result.data.path"),"dada2.data", sep = "/"))
add.param(PARAM) <- c("ASVtable.data", paste(get.value("result.data.path"),"ASVtable.data", sep = "/"))

add.param(PARAM) <- c("OTUtable.data", paste(get.value("result.data.path"),"OTUtable.data", sep = "/"))

add.param(PARAM) <- c("CORRECTEDtable.data", paste(get.value("result.data.path"),"CORRECTEDtable.data", sep = "/"))

add.param(PARAM) <- c("ALLtable.data", paste(get.value("result.data.path"),"ALLtable.FINAL.data", sep = "/"))

# IDT training set
add.param(PARAM) <- c("IDT.TS.data", paste(get.value("result.data.path"),"IDT.TS.data", sep = "/"))

# Blast99
add.param(PARAM) <- c("Blast99.data", paste(get.value("result.data.path"),"Blast99.data", sep = "/"))

# Other stats
add.param(PARAM) <- c("BasicStats.data", paste(get.value("result.data.path"),"BasicStats.data", sep = "/"))


# Save Parameters ---------------------------------------------------------

PARAM

write.csv2(PARAM, file = file.path("./04_Log", "Options.csv"), row.names=F)

get.value("result.path")



