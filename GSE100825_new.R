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


#Now put data into an methylset object
library(minfi)
annotation <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19",force = T)

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
 ?ratioConvert 
  
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


