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
GSM_IDs
column_names
test

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
  as.tibble() %>%
  mutate_all(str_extract,pattern="(?<=:[[:blank:]]).*")%>%
  mutate(GEO_accession = GSM_IDs)


#Load raw data 

raw_heading <- read.table("GSE100825_SIgnalIntensityMatrix.txt.gz",
                          nrows = 1)

# File was not able to find in the directory. Therefore, signal intensity file was downloaded
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100825/suppl/GSE100825_signal_intensities.txt.gz",destfile = "cpg.txt.gz")
u<-R.utils::gunzip("cpg.txt.gz",overwrite=T)
read_file("cpg.txt")
### side note- signal intensity was not found cz raw data file was not extracted. Anyhow, signal intensity file was downloaded
# Signal intensity file was moved to the directory
file.copy("cpg.txt.gz", "GSE100825_s/")

getwd()
# Reading signal intensity/ raw data 
## checking for the content of first row in the raw data table
read.table("cpg.txt",sep = " ")
raw_heading_s <- read.table("cpg.txt",
                          nrows = 1)
#skipping the first line in the row of the raw data
raw_s <- read.table("cpg.txt.gz",
                  skip = 1)
CpGs <- pull(raw_s,
             var=1)




