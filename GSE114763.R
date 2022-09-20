
getwd()
setwd("/home/mda/Documents/data")
getwd()

# Loading libraries
install.packages("HelpersMG")
library(GEOquery)
library(ChAMP)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(magrittr)
library(HelpersMG)
library(limma)

#Making a directory for the files 

dir.create("GSE114763")
WORKING_DIR="GSE114763"
ARRAY_DATA="GSE114763.tar"
DEST=paste(WORKING_DIR,"/",ARRAY_DATA,sep="")

# Getting the GEO file
gse_s <- getGEO("GSE114763",
                GSEMatrix = F,
                getGPL = F,
                destdir = "GSE114763/")

#Finding GEO ids in the file
GSM_IDs <- gse_s@header$sample_id
class(GSM_IDs)


# Getting pheno characters into a df

## In here Initializing pheno variable as null so that it can be passed via the loop.
## getGPL in here will be false, if not data@header$characteristics_ch1 row will be duplicated

pheno <- NULL

for (g in GSM_IDs) {
  data<-getGEO(g,
               GSEMatrix = TRUE,
               AnnotGPL = FALSE,
               getGPL = F)
  pheno<- rbind(pheno,
                data@header$characteristics_ch1)  
}


# Extracting first word from a string in character 1 and naming those extracted words as column names 

test <- pheno
colnames <- as.character(sapply(test[1,],
                                str_extract,
                                pattern = ".*(?=:)"))


test <- pheno
column_names <-str_extract(test[1,],pattern = ".*(?=:)")
colnames(test)<-column_names

# Adding GSM numbers to the pheno characters
test<-test%>%
  as_tibble() %>%
  mutate_all(str_extract,pattern="(?<=:[[:blank:]]).*")%>%
  mutate(GEO_accession = GSM_IDs)


################Load raw data ################
getwd()
setwd("/home/mda/Documents/data/GSE114763/")

#Download raw tar files
system('wget -O idats.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE114763&format=file"')
file.rename("idats.tar.gz",ARRAY_DATA)
untar(exdir = "GSE114763/IDAT", tarfile = ARRAY_DATA)


############################
#Since Series matrix is a file preprocess by the authors , avoid using series matrix for analysis. 
#Use always idat files. 
#Download series matrix file 
#getwd()
#setwd("/home/mda/Documents/data/GSE100825_s/")
#if (! file.exists("GSE100825_series_matrix.txt.gz") ) { 
  #URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100825/matrix/GSE100825_series_matrix.txt.gz"
  #download.file(URL,destfile = "GSE100825_series_matrix.txt.gz")
#}
#  when the network is too slow
#system('wget -O GSE100825_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100825/matrix/GSE100825_series_matrix.txt.gz"')

#Download signal matrix file

download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100825/suppl/GSE100825_signal_intensities.txt.gz'",destfile = "GSE100825_SIgnalIntensityMatrix.txt.gz")
R.utils::gunzip("GSE100825_signal_intensities.txt.gz",overwrite=F)
# Connection issue so i used Wget()
system('wget -O GSE100825_SIgnalIntensityMatrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100825/suppl/GSE100825_signal_intensities.txt.gz"')

## getting the raw data and cpg sites
raw_heading <- read.table("GSE100825_SIgnalIntensityMatrix.txt.gz",
                          nrows = 1)

# Reading signal intensity/ raw data 
## checking for the content of first row in the raw data table

raw_heading <- read.table("GSE100825_SIgnalIntensityMatrix.txt.gz",
                          nrows = 1)
#skipping the first line in the row of the raw data
raw_s <- read.table("GSE100825_SIgnalIntensityMatrix.txt.gz",
                    skip = 1)
CpGs <- pull(raw_s,
             var=1)

#######Shortcut for getting CpGs 
Raw<-read.delim2("GSE100825_SIgnalIntensityMatrix.txt.gz", header = TRUE, sep = "\t", dec = ",")
Cpg<-Raw$ID_REF
#####################################################################

#Read GPL FILE 
Matrix_table<-read.csv("IDAT/GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz", header = T, sep = ",")

######################################################################################
gse <- getGEO(filename="GSE100825_series_matrix.txt.gz")
dim(gse)
sample_metadata <- pData(phenoData(gse))
dim(sample_metadata)
targets <- sample_metadata
setwd("/home/mda/Documents/data/")
files <- list.files(WORKING_DIR,pattern = "GSM",recursive = TRUE)
files
mybase <- unique(gsub("_Red.idat.gz" ,"", gsub("_Grn.idat.gz", "" ,files)))
mybase
mybase <- paste(WORKING_DIR,"/",mybase,sep="")
class(mybase)
targets$Basename <- mybase
head(targets)
dim(targets)