library(minfi)
library(GEOquery)
library(ChAMP)
library(tidyverse)
library(ggplot2)


dir.create("GSE100825_old")


#file was not getting downloaded so used wget
if ( !dir.exists("GSE100825_old/IDAT") ) {
  dir.create("GSE100825_old/IDAT")
  if ( !file.exists(ARRAY_DATA)  ) {
    system('wget -O idats.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100825&format=file"')
    file.rename("idats.tar.gz",ARRAY_DATA)
    untar(exdir = "GSE100825_old/IDAT", tarfile = ARRAY_DATA)
  }
}

# Finding methyl array sheets
baseDir <- "GSE100825_old"
targets <- read.metharray.sheet(baseDir)

#Since it is not the correct array datasheet, 

if (! file.exists("GSE100825_series_matrix.txt.gz") ) { 
  URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100825/matrix/GSE100825_series_matrix.txt.gz"
  download.file(URL,destfile = "GSE100825_series_matrix.txt.gz")
}
gse<- getGEO(filename = "GSE100825_series_matrix.txt.gz")



