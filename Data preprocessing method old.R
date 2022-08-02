
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("minfi")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChAMP")


library(minfi)
library(GEOquery)
library(ChAMP)
library(tidyverse)
library(ggplot2)




WORKING_DIR="GSE100825_NEW"
ARRAY_DATA="GSE100825_RAW.tar"
DEST=paste(WORKING_DIR,"/",ARRAY_DATA,sep="")

if(!dir.exists(WORKING_DIR)){
  dir.create(WORKING_DIR)
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100825/suppl/GSE100825_RAW.tar",
                destfile=DEST) 
  untar(exdir = "IDAT", tarfile = WORKING_DIR)
}
SERIES_MATRIX=paste(WORKING_DIR,"/","GSE100825_series_matrix.txt.gz",sep="")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100825/matrix/GSE100825_series_matrix.txt.gz", 
              destfile=SERIES_MATRIX)

gse <- getGEO(filename=SERIES_MATRIX)

baseDir <- "GSE100825_NEW"
sample_metadata <- pData(phenoData(gse))
targets <- sample_metadata

#use pdata to find the phenotypic data in the set
targets <-pData(gse)

#finding methyl_array sample sheet
methyl_array<- read.metharray.sheet(baseDir)

#therefore alternative approach

baseDir <- "GSE100825_NEW"
R.utils::gunzip("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL21nnn/GPL21145/suppl/GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz",overwrite=TRUE, remove=FALSE)
waaaa <- read.metharray.sheet(baseDir,pattern="csv")


files <- list.files(WORKING_DIR,pattern = "GSM",recursive = TRUE)
mybase <- unique(gsub("_Red.idat.gz" ,"", gsub("_Grn.idat.gz", "" ,files)))
mybase <- paste(WORKING_DIR,"/",mybase,sep="")
targets$Basename <- mybase
head(targets)



files <- list.files(WORKING_DIR,pattern = "GSM",recursive = TRUE)
mybase <- unique(gsub("_Red.idat.gz" ,"", gsub("_Grn.idat.gz", "" ,files)))
mybase <- paste(WORKING_DIR,"/",mybase,sep="")
targets$Basename <- mybase
head(targets)

targets$Basename <-mybase
rgSet <- read.metharray.exp(targets = targets)
?read.metharray.exp
unique(targets$Basename)
