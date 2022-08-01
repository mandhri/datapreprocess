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

baseDir <- "."
sample_metadata <- pData(phenoData(gse))
targets <- sample_metadata

targets <- pData(phenoData(gse))
targets <- targets[order(rownames(targets)),]
mybase <- unique(gsub("_Red.idat.gz" ,"", gsub("_Grn.idat.gz", "" ,list.files("./GSE100825",pattern = "GSM",recursive = TRUE))))
mybase <- paste("GSE100825/", mybase, sep = "")
head(targets)

targets$basename <-mybase
rgSet <- read.metharray.exp(targets = targets)
mSet <- preprocessRaw(rgSet)
library("minfi")
