getwd()
setwd("/home/mandhri/Data preprocess/Datapreprocess/GSE55763/GSE55763/")
# Reconfirming once again for missing values#####
phenotypes<- c()
gsm_nos <- c()
sample<- c()
for (i in 1:length(gse_s@gsms)) {
  phenotypes[[i]] <- gse_s@gsms[[i]]@header$characteristics_ch1
  gsm_nos[[i]] <- gse_s@gsms[[i]]@header$geo_accession
  sample[[i]]<-gse_s@gsms[[i]]@header$description
}

## For samples
samples<-as.data.frame.Date(sample)
dim(samples)
names(samples)[names(samples) == "sample"]<- "sampling_grp"

## For GSM NUMBER 
Gsm_num<-as.data.frame.Date(gsm_nos)
dim(Gsm_num)
names(Gsm_num)[names(Gsm_num) == "Gsm_num"]<- "GSM_ID"

## For Phenotypes
phenotypes<-as.data.frame(phenotypes)
dim(phenotypes)
phenotypes<-as.data.frame(t(phenotypes))

###
phenotypes<-setNames(cbind(rownames(phenotypes), phenotypes, row.names = NULL), 
                c("COL1", "V1", "V2","V3","V4"))

#Renamed and removed the first column 
phenotypes<-phenotypes[,2:5]
phenotypes= subset(phenotypes, select = -c(V2) )
getwd

# Splitting the column with age to extract just the age

phenotypes<-sapply(strsplit(phenotypes$V4, " "), "[[", 3)
class(age)
phenotypes<-cbind(phenotypes,age)

# Splitting the column with sex to extract just the sex
sex<-sapply(strsplit(phenotypes$V3, " "), "[[", 2)
class(sex)
phenotypes<-cbind(phenotypes,sex)

# Splitting the column with tissue to extract just the tissue
tissue<-sapply(strsplit(phenotypes$V2, " "), "[[", 1)
tissue
phenotypes<-cbind(phenotypes,tissue)

########

reconfirm<-cbind(Gsm_num,phenotypes,samples)
reconfirm= subset(reconfirm, select = -c(V1,V3,V4) )

summary(pheno12)


