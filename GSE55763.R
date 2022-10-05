library(GEOquery)
library(ChAMP)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(magrittr)
library(HelpersMG)
library(limma)
library(minfi)

dir.create("GSE55763")
WORKING_DIR="GSE55763"
ARRAY_DATA="GSE55763.tar"
DEST=paste(WORKING_DIR,"/",ARRAY_DATA,sep="")

# Getting the GEO file
gse_s <- getGEO("GSE55763",
                GSEMatrix = F,
                getGPL = F,
                destdir = "GSE55763/")


gse<-getGEO("GSE55763")



#Finding GEO ids in the file
GSM_IDs <- gse_s@header$sample_id

#Finding the number of GSM Ids in the GSM_IDs
str(GSM_IDs)

#Lets see the str of the gse file
str(gse_s)

################Load raw data ################
getwd()
setwd("/home/mandhri/Data preprocess/Datapreprocess/GSE55763/")
system('wget -O idats.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55763&format=file"')
  file.rename("idats.tar.gz",ARRAY_DATA)
  untar(exdir = "GSE55763/IDAT", tarfile = ARRAY_DATA)

baseDir <- "GSE55763"  
R.utils::gunzip("/home/mandhri/Data preprocess/Datapreprocess/GSE55763/GSE55763/IDAT/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz",overwrite=TRUE, remove=FALSE)
targets <- read.metharray.sheet(baseDir,pattern="csv$")
  
Matrix_table<-read.csv("GSE55763/IDAT/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz", header = T, sep = ",",skip = 7)

pheno=NULL
pheno<- c()

for (i in 1:length(gse_s@gsms)) {
  print(gse_s@gsms[[i]]@header$characteristics_ch1)
  pheno[[i]] <- gse_s@gsms[[i]]@header$characteristics_ch1
  }
length(pheno)
pheno<-as.data.frame.list(pheno)
pheno<-rbind(GSM_IDs,pheno)
pheno <- as.data.frame(t(pheno))

#
dim(pheno)

#setting name for the 1st column which doesn't have the name 

pheno<-setNames(cbind(rownames(pheno), pheno, row.names = NULL), 
         c("COL1", "V1", "V2","V3","V4","V5"))

#Removing the first column and keeping the rest
dim(pheno)
pheno<-pheno[,2:6]

# Splitting the column with age to extract just the age

age<-sapply(strsplit(pheno$V5, " "), "[[", 2)
class(age)
pheno<-cbind(pheno,age)

# Splitting the column with sex to extract just the sex
sex<-sapply(strsplit(pheno$V4, " "), "[[", 2)
class(sex)
pheno<-cbind(pheno,sex)

# Splitting the column with tissue to extract just the tissue
tissue<-sapply(strsplit(pheno$V2, " "), "[[", 3)
tissue
pheno<-cbind(pheno,tissue)


##
pheno= subset(pheno, select = -c(V2,V3,V4,V5))

summary(pheno)


#checking for NA values 
is.na(pheno)
newdata <- na.omit(pheno)
pheno$sex
sum(is.na(pheno))
## Confirmed that there is no missing values

#Checking the number of females and males in the sample
n_occur <- data.frame(table(pheno$sex))
n_occur[n_occur$Freq > 1,]


unique(pheno$V3)
###
#Need to extract the description category to see the sampling type


sampling_group<-c()
for (x in 1:length(gse_s@gsms)) {
  print(gse_s@gsms[[x]]@header$description) 
  sampling_group[[x]] <- gse_s@gsms[[x]]@header$description 
}

#Change the sampling_group list to a df and combine with previous pheno df
pheno1<-as.matrix(sampling_group)
pheno1<-as.data.frame(sampling_group)
names(pheno1)[names(pheno1) == "V1"]<- "sampling_grp"
dim(pheno1)
pheno2<-as.data.frame(t(pheno1))


##Renaming the first column 
pheno2<-setNames(cbind(rownames(pheno2), pheno2, row.names = NULL), 
                c("COL1", "V1"))

#Removing the first column and keeping the rest
dim(pheno2)
pheno2<-pheno2[,2]

final_pheno<-cbind(pheno,pheno2)


library(arsenal)
install.packages("arsenal")
#checking for missing values in the final_pheno
sum(is.na(final_pheno))
summary(is.na(final_pheno))

##Comparing two df to check whether there is any variation 
comparedf(reconfirm,final_pheno)
summary(comparedf(reconfirm,final_pheno))

### No missing values or replicated samples found in both the files. Therefore, this dataset can not have any duplicate 
#or missing values at this stage##########