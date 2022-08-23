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

dir.create("GSE100825_s")
WORKING_DIR="GSE100825_s"
ARRAY_DATA="GSE100825_s.tar"
DEST=paste(WORKING_DIR,"/",ARRAY_DATA,sep="")

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


################Load raw data ################
getwd()
setwd("/home/mda/Documents/data/GSE100825_s/")

#Download raw tar files
system('wget -O idats.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100825&format=file"')
file.rename("idats.tar.gz",ARRAY_DATA)
untar(exdir = "GSE100825_s/IDAT", tarfile = ARRAY_DATA)

#Download series matrix file 
getwd()
setwd("/home/mda/Documents/data/GSE100825_s/")
if (! file.exists("GSE100825_series_matrix.txt.gz") ) { 
  URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100825/matrix/GSE100825_series_matrix.txt.gz"
  download.file(URL,destfile = "GSE100825_series_matrix.txt.gz")
}
#  when the network is too slow
system('wget -O GSE100825_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100825/matrix/GSE100825_series_matrix.txt.gz"')

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

sample_metadata <- pData(phenoData(gse))
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

# Make sure you have the files in the same directory. Due to having files in different wd(), i had to spend 2 days trying to find the error. Sigh!!!!. I forgot the rule 1
getwd()
rgset<-read.metharray.exp(targets = targets)
rgset

#Reading sample sheet (GPL)

baseDir <- "GSE100825_s"
sample_sheet_GPL<-read.metharray.sheet(baseDir)

# Lets get Sentrix ID and Sentrix Position from column names

sentrixIDs_post<-targets$description.1
sentrixIds<- str_extract(sentrixIDs_post,"^[^_:]+")
class(sentrixIds)
targets$SentrixID <-sentrixIds
Sentrix_positions<-str_extract(sentrixIDs_post,"_.+")
Sentrix_positions
targets$Sentrix_positions <-Sentrix_positions
targets

#Add sample name to merge with matrices
test <- test %>%
  mutate(Sentrix_ID = sentrixIds,
         Sentrix_Position = Sentrix_positions)


write_csv(test,
          file = "GSE100825 phenotypes.csv")


#Extract meth, unmeth and detP signals
meth <- raw_s[,seq(3,ncol(raw_s),by=3)] %>% as.matrix()
rownames(meth) <- CpGs
colnames(meth) <- test$GEO_accession
unmeth <- raw_s[,seq(2,ncol(raw_s),by=3)] %>% as.matrix()
rownames(unmeth) <- CpGs
colnames(unmeth) <- test$GEO_accession
detP<- raw_s[,seq(4,ncol(raw_s),by=3)] %>% as.matrix()
rownames(detP) <- CpGs
colnames(detP) <- test$GEO_accession


##### For the raw data did in short cut

#methylated 
methylated <-as.matrix(Raw[,seq(3,ncol(Raw),by=3)])
rownames(methylated )<-Cpg
colnames(methylated )<-GSM_IDs
methylated 
#Unmethylated
unmethylated<-as.matrix(Raw[,seq(2,ncol(Raw),by=3)])
rownames(unmethylated)<-Cpg
colnames(unmethylated)<-GSM_IDs
unmethylated

# detP signals
p_val<-as.matrix(Raw[,seq(4,ncol(Raw),by=3)])
rownames(p_val)<-Cpg
colnames(p_val)<-GSM_IDs
p_val

## Shortcut for detection of p val from rgset
p_val<-detectionP(rgset)
######


#Now put data into an methylset object
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annotation <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19",force = T)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
names(annotation) <-"annotation"

methylset=MethylSet(Meth=meth,
                    Unmeth=unmeth,
                    annotation = annotation)

methylsetR<- MethylSet(Meth = methylated,
                       Unmeth = unmethylated,
                       annotation = annotation)
RSet <- ratioConvert(methylset,
                     what = "both",
                     keepCN = TRUE)

RsetR<- ratioConvert(methylsetR,
                     what="both",
                     keep=T)
  
GRset <- mapToGenome(RSet)

mapping_to_genome<-mapToGenome(RsetR)
?mapToGenome
predictedSex <- getSex(GRset,
                       cutoff = -2)$predictedSex

#Troubleshoot gender in this 
test <- test %>%
  mutate(gender = ifelse(gender=="Male",
                         "M",
                         "F"))
test$gender==predictedSex #All sexes match!
targets$`Sex:ch1`

###############

library(minfi)

#method1
rgset
mset<-preprocessRaw(rgset)
rset<-ratioConvert(mset, what = "both", keepCN = TRUE)
grset<-mapToGenome(rset)
predictedSex<-getSex(grset,cutoff = -2)$predictedSex


#method2
## sex chromosome include

detp<- detectionP(rgset)
keep<-rowSums(detp<0.01) ==ncol(rgset)
include_sexchrom<-mset[keep,]

## SNP exclude
include_sexchrom1<-mapToGenome(include_sexchrom)
include_sexchrom1_nosnp<-dropLociWithSnps(include_sexchrom1)
dim(include_sexchrom1_nosnp)

## sex chromsome exclude
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
xyprobes<- anno$Name[anno$chr %in% c("chrX","chrY")]
xyprobes
nosex_chromo<-include_sexchrom1[which(!rownames(include_sexchrom1) %in% xyprobes)]


## beta and m values with sex chromosome 
meth_with_sexchrom<-getMeth(include_sexchrom1)
unmeth_with_sexchrom<-getUnmeth(include_sexchrom1)
Mval_include_sexchrom1<-log2((meth_with_sexchrom+100)/(unmeth_with_sexchrom+100))
beta_include_sexchrom1<-getBeta(include_sexchrom1)
dim(Mval_include_sexchrom1)


## beta and m values with sex chromosome 
meth_no_sexchromo<-getMeth(nosex_chromo)
unmeth_no_sexchromo<-getUnmeth(nosex_chromo)
Mval_no_sexchromo<-log2((meth_no_sexchromo+100)/(unmeth_no_sexchromo+100))
beta_no_sexchromo<-getBeta(nosex_chromo)
dim(Mval_no_sexchromo)

#let's go graphic with MDS plots to depict sex chromosome seperation

targets$`Sex:ch1`<-as.factor(targets$`Sex:ch1`)
library(limma)
plotting<-plotMDS(Mval_include_sexchrom1, labels=targets$`Sex:ch1`, col=colours,cex=0.7)

#Interestingly, top probe (GSM2694066) is a male participant but is more towards female side.Need to do a sex diagnostic


