getwd()
setwd("/home/mandhri/Documents/R data/")
getwd()

# Loading libraries

library(GEOquery)
library(ChAMP)
library(tidyverse)
library(ggplot2)

#Making a directory for the files 

dir.create("GSE100825")

ARRAY_DATA="GSE74548/GSE74548_RAW.tar"
# only download it if it is not present on the system
if ( !dir.exists("GSE74548/IDAT") ) {
  dir.create("GSE74548/IDAT")
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE74548&format=file",
                destfile = ARRAY_DATA)
  untar(exdir = "GSE74548/IDAT", tarfile = ARRAY_DATA)
}



# Getting the GEO file
gse <- getGEO("GSE100825",
              GSEMatrix = F,
              getGPL = F)
?getGEO

GSM_IDs <- gse@header$sample_id

pheno <- rbind()
for (g in GSM_IDs)
{
  data <- getGEO(g,
                 GSEMatrix = TRUE,
                 AnnotGPL = FALSE,
                 getGPL = FALSE)
  
  pheno <- rbind(pheno,
                 data@header$characteristics_ch1)
}

test <- pheno
colnames <- as.character(sapply(test[1,],
                                str_extract,
                                pattern = ".*(?=:)"))



colnames(test)
colnames(test) <- colnames
test <- test %>%
  as_tibble%>%
  mutate_all(str_extract,
             pattern = "(?<=:[[:blank:]]).*")%>%
  mutate(`GEO accession` = GSM_IDs,
         Sample_Name = GSM_IDs)

