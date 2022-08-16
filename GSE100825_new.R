getwd()
setwd("/home/mda/Documents/data")
getwd()

# Loading libraries

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(magrittr)


#Making a directory for the files 

dir.create("GSE100825_s")

# Getting the GEO file
gse_s <- getGEO("GSE100825",
              GSEMatrix = F,
              getGPL = F,
              destdir = "GSE100825_s/")

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


#Load raw data 

raw_heading <- read.table("GSE100825_SIgnalIntensityMatrix.txt.gz",
                          nrows = 1)

# File was not able to find in the directory. Therefore, signal intensity file was downloaded
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100825/suppl/GSE100825_signal_intensities.txt.gz",destfile = "cpg.txt.gz")
R.utils::gunzip("cpg.txt.gz",overwrite=F)

### side note- signal intensity was not found cz raw data file was not extracted. Anyhow, signal intensity file was downloaded
# Signal intensity file was moved to the directory
file.copy("cpg.txt.gz", "GSE100825_s/")

# Reading signal intensity/ raw data 
## checking for the content of first row in the raw data table

raw_heading <- read.table("cpg.txt",
                          nrows = 1)
#skipping the first line in the row of the raw data
raw_s <- read.table("cpg.txt.gz",
                  skip = 1)
CpGs <- pull(raw_s,
             var=1)

#Shortcut for getting CpGs 
Raw<-read.delim2("cpg.txt", header = TRUE, sep = "\t", dec = ",")
Cpg<-Raw$ID_REF


#Load Sentrix ID and Position and add to pheno table
library(readxl)
batch <-("fam.xml") %>%
  rename(`twin id`=`Sample name`)
test <- left_join(test,
                  batch)
if(!require('methylumi')) {
  install.packages('methylumi')
  library('methylumi')
}

#Download series_matrix
if (! file.exists("GSE100825_series_matrix.txt.gz") ) { 
  URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100825/matrix/GSE100825_series_matrix.txt.gz"
  download.file(URL,destfile = "GSE100825_series_matrix.txt.gz")
}
R.utils::gunzip("GSE100825_series_matrix.txt.gz",overwrite=F)

Matrix_table <- read.delim2("GSE100825_series_matrix.txt.gz", header = TRUE, sep = "\t", dec = ",")
head(Matrix_table)


#RAW IDAT FILES 

#Raw data file was not getting downloaded so used wget
ARRAY_DATA="GSE100825_RAW.tar"
WORKING_DIR="GSE100825_s"
DEST=paste(WORKING_DIR,"/",ARRAY_DATA,sep="")

if ( !dir.exists("GSE100825_s/IDAT") ) {
  dir.create("GSE100825_s/IDAT")  
if ( !file.exists(ARRAY_DATA)  ) {
    system('wget -O idats.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100825&format=file')
    file.rename("idats.tar.gz",ARRAY_DATA)
    untar(exdir = "GSE100825_s/IDAT", tarfile = ARRAY_DATA)
  }
}          


download.file("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL21nnn/GPL21145/suppl/GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz",destfile = "GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz")
R.utils::gunzip("GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz",overwrite=F)
lp<-read.csv("GPL21145_MethylationEPIC_15073387_v-1-0.csv", header = T, sep = ",")
getwd()
R.utils::gunzip("GSE100825.soft.gz",overwrite=F)

library(Biobase)
lets_see<-read_file("GSE100825.soft")

head(lets_see)
str(lets_see)





