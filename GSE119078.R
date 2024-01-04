

# Loading libraries
library(GEOquery)
library(ChAMP)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(magrittr)
library(HelpersMG)
library(limma)
library(readxl)
library(sva)
#Making a directory for the files 
WORKING_DIR="GSE119078"
baseDir<-WORKING_DIR
dir.create("GSE119078")
ARRAY_DATA="GSE119078.tar"
DEST=paste(WORKING_DIR,"/",ARRAY_DATA,sep="")

#Annotation

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
myann <- data.frame(ann450k[,c("UCSC_RefGene_Name","Regulatory_Feature_Group")])
promoters <- grep("Prom",myann$Regulatory_Feature_Group)

# Getting the GEO file
gse_s <- getGEO("GSE119078",
              GSEMatrix = F,
              getGPL = F,
              destdir = "GSE119078/") 

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

#Download raw tar files
system('wget -O idats.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE119078&format=file"')
file.rename("idats.tar.gz",ARRAY_DATA)
untar(exdir = "GSE119078/IDAT", tarfile = ARRAY_DATA)
#open up the gpl file from the untar files
R.utils::gunzip("/home/mandhri/nerw/GSE119078/IDAT/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz",overwrite=T)


############################
#Since Series matrix is a file preprocess by the authors , avoid using series matrix for analysis. 
#Use always idat files. moreover , no signal matrix file was found. Anyways you can use gpl file from idat files to get signal matrix file.
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

#download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100825/suppl/GSE100825_signal_intensities.txt.gz'",destfile = "GSE100825_SIgnalIntensityMatrix.txt.gz")
#R.utils::gunzip("GSE100825_signal_intensities.txt.gz",overwrite=F)
# Connection issue so i used Wget()
#system('wget -O GSE100825_SIgnalIntensityMatrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100825/suppl/GSE100825_signal_intensities.txt.gz"')

#####################################################################
setwd("~/nerw/GSE114763/IDAT")
#Read GPL FILE 
Matrix_table<-read.csv("/home/mandhri/nerw/GSE114763/IDAT/GPL13534_HumanMethylation450_15017482_v.1.1.csv", header = T, sep = ",")

######################################################################################
gse_s <- getGEO(filename="GSE119078_series_matrix.txt.gz")
gse_s <- getGEO("GSE119078",
                GSEMatrix = T,
                getGPL = F,
                destdir = "GSE119078/")

gse_s <- getGEO(filename="GSE119078_series_matrix.txt.gz")

dim(gse_s)
sample_metadata <- pData(phenoData(gse_s))
dim(sample_metadata)
targets <- sample_metadata
setwd("~/nerw/GSE119078/IDAT")

files <- list.files(getwd(),pattern = "GSM",recursive = TRUE)
files
mybase <- unique(gsub("_Red.idat.gz" ,"", gsub("_Grn.idat.gz", "" ,files)))
mybase
#mybase <- paste(WORKING_DIR,"/",mybase,sep="")
class(mybase)
targets$Basename <- mybase
head(targets)
dim(targets)
# Make sure you have the files in the same directory. Due to having files in different wd(), i had to spend 2 days trying to find the error. Sigh!!!!. I forgot the rule 1
getwd()
rgset<-read.metharray.exp(targets = targets)
rgset

#Reading sample sheet (GPL)

baseDir <- "GSE100825_s"
sample_sheet_GPL<-read.metharray.sheet(baseDir)

# Lets get Sentrix ID and Sentrix Position from column names
Basename
sentrixIDs<-targets$description.1
sentrixIds<- str_extract(mybase,"(?<=_)[0-9]+")
class(sentrixIds)
targets$SentrixID <-sentrixIds
sentrixPositions <- sapply(strsplit(mybase, "_"), function(x) x[3])
sentrixPositions
targets$sentrixPositions <-sentrixPositions
head(targets)

#Add sample name to merge with matrices
test <- test %>%
  mutate(Sentrix_ID = sentrixIds,
         Sentrix_Position = sentrixPositions,
         Basename=mybase,
         Tissue=targets$source_name_ch1)


write_csv(test,
          file = "GSE119078 phenotypes.csv")



###############

library(minfi)

#method1
rgset
mset<-preprocessRaw(rgset)
rset<-ratioConvert(mset, what = "both", keepCN = TRUE)
grset<-mapToGenome(rset)
predictedSex<-getSex(grset,cutoff = -2)$predictedSex
targets$predictedSex<-predictedSex


# Matching predicted sex column with sex/ gender column 
test <- test %>%
  mutate(sex_match = ifelse((predictedSex == 'M' & gender == 'Male') | 
                              (predictedSex == 'F' & gender == 'Female'), 
                            'match', 'mismatch'))

#method2- graphical rep

## sex chromosome include

detp<- detectionP(rgset)
keep<-rowSums(detp<0.01) ==ncol(rgset)
include_sexchrom<-mset[keep,]
dim(detp)
dim(rgset)


# exclude SNP probes
mSetSw <- mset[keep,]
mSetSw <- mapToGenome(mSetSw)
mSetSw_nosnp <- dropLociWithSnps(mSetSw)
dim(mSetSw)
dim(mSetSw_nosnp)
mSetSw <- mSetSw_nosnp

# exclude sex chromosomes
keep <- !(featureNames(mSetSw) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
mSetFlt <- mSetSw[keep,]
mSetFlt[1:6,1:5]
dim(mSetFlt)


## Extracting Beta and M-values


# include sex chromosomes
meth <- getMeth(mSetSw)
unmeth <- getUnmeth(mSetSw)
Mval <- log2((meth + 100)/(unmeth + 100))
beta <- getBeta(mSetSw)

# exclude sex chromosomes
meth <- getMeth(mSetFlt)
unmeth <- getUnmeth(mSetFlt)
Mval_flt <- log2((meth + 100)/(unmeth + 100))
beta_flt <- getBeta(mSetFlt)

# MDS plots to depict sex chromosome seperation
colnames(Mval) <- sapply(strsplit(colnames(Mval),"_"),"[[",1)
colnames(Mval_flt) <- sapply(strsplit(colnames(Mval_flt),"_"),"[[",1)

cgx <- rownames(Locations[which(Locations$chr %in% "chrX"),])
cgy <- rownames(Locations[which(Locations$chr %in% "chrY"),])

mvx <- Mval[which(rownames(Mval) %in% cgx),]
mvy <- Mval[which(rownames(Mval) %in% cgy),]

targets_m <- rownames(subset(targets,`gender:ch1`=="Male"))
str(targets_m)
targets_f <- rownames(subset(targets,`gender:ch1`=="Female"))
str(targets_f)

Mvalm <- Mval[,colnames(Mval) %in% targets_m]
Mvalf <- Mval[,colnames(Mval) %in% targets_f]

mvxm <- Mvalm[which(rownames(Mvalm) %in% cgx),]
mvym <- Mvalm[which(rownames(Mvalm) %in% cgy),]

mvxf <- Mvalf[which(rownames(Mvalf) %in% cgx),]
mvyf <- Mvalf[which(rownames(Mvalf) %in% cgy),]


plot(colMeans(mvx),colMeans(mvy),col="gray")
points(colMeans(mvxm),colMeans(mvym),col="lightblue",pch=19,cex=1.5)
points(colMeans(mvxf),colMeans(mvyf),col="pink",pch=19,cex=1.5)
text(colMeans(mvx),colMeans(mvy),labels = colnames(mvx),cex=0.75)
dev.off()

##################Making beta-matrix

myLoad <- champ.load(directory=getwd(),arraytype="450K")

#Give the column names to the beta matrix

colnames(myLoad$beta)<-mybase
# Use the place where you have your IDAT files and phenotype file 

sheets <- excel_sheets("/home/mandhri/Cross reactive and SNP probes from Pidsley et al.xlsx")
beta = myLoad$beta
for (s in sheets)
{
  qcprobes <- read_excel('/home/mandhri/Cross reactive and SNP probes from Pidsley et al.xlsx',
                         sheet=s)%>%
    pull(ProbeID)
  beta <- beta[setdiff(rownames(beta),qcprobes),]
}


#Produce quality control graphs to look at the data
champ.QC(beta = beta,
         pheno =myLoad$pd$gender,
         dendrogram = FALSE,
         resultsDir="./CHAMP_QCimages/")


#Normalization of Type I and Type II probes
Champ_Norm <- champ.norm(beta=beta)

## No beta values missing , so no need to impute  and run QC again
#impute = champ.impute(beta = beta,pd = pheno)
#champ.QC(beta = impute$beta,pheno = impute$pd$sex,dendrogram = FALSE,resultsDir = "./CHAMP_QCimages")

## Incase there is missing values, then first get an estimate of the missing values and then 
# decide whether you want to impute or remove the missing cpgs.


#num_NA <- sum(is.na(beta))
#num_Inf <- sum(is.infinite(beta))

#cat("Number of NA values:", num_NA, "\n")
#cat("Number of Infinite values:", num_Inf, "\n")

#probeNA <- rowMeans(is.na(beta)) * 100
#summary(probeNA)

#sampleNA <- colMeans(is.na(beta)) * 100  # in percentage
#summary(sampleNA)


write.table(beta, "/mnt/vol1/buccal/GSE119078_beta.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")



write.table(test,"/mnt/vol1/buccal/GSE119078_phenotype.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

#testing the files
f<-read.delim("/mnt/vol1/buccal/GSE119078_beta.txt")

f<-read.delim("/mnt/vol1/buccal/GSE119078_phenotype.txt")


champ.SVD(beta=Champ_Norm,
              pd=myLoad$pd%>%dplyr::select(Slide,Basename,gender,age,disease.status,Array,Tissue),
          resultsDir="./CHAMP_SVDimages/")

#THe above didnt work cz of dimensions of the champ_norm needs to be 1 not 2 to run champ.svd.
class(Champ_Norm) %>% length() 


champ.SVD(beta=Champ_Norm%>% as.data.frame(),
          pd=myLoad$pd%>%dplyr::select(Slide,Basename,gender,age,disease.status,Array,Tissue),
          resultsDir="./CHAMP_SVDimages/")

## If you look at the screeplot generated by champ.svd, you need to search for the principal components hat create the elbow in the plot.
# In your plot, the elbow appears to be around Component3 or Component4, suggesting that the first 3 or 4 components are the most significant in the datset. 

#Correcting for batch effect

M <- logit2(Champ_Norm)

myCombat=ComBat(dat=as.matrix(M), #it outputs an M-value matrix adjusted for batch
                batch=test$Sentrix_ID,
                mod=NULL)

myCombat=ilogit2(myCombat)

champ.QC(beta = myCombat,
         pheno = test$gender,
         dendrogram = FALSE,
         resultsDir="./CHAMP_QCimages/combat_QCimages/")
getwd()

write.table(myCombat, "/mnt/vol1/buccal/GSE119078_beta.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

#Run SVD again
champ.SVD(beta=as.data.frame(myCombat),
          pd=pheno,
          resultsDir="./CHAMP_SVDimages/batch_corrected/")
