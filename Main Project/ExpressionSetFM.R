source("http://bioconductor.org/biocLite.R")
biocLite("BioUpgrade")
library(Biobase)
biocLite('affy')
library('affy')
biocLite('gahgu95av2.db')
library(gahgu95av2.db)
library(convert)
C<-getwd()
C=sub("Main Project", "Data", C)
setwd(C)
data=read.table("datasetA_scans.txt", header = TRUE, sep="\t")
data=data[c(191:205,211:225),]
opis=opis[c(191:205,211:225),]
opis=read.AnnotatedDataFrame("datasetA_scans.txt", sep="\t", header = TRUE, row.names=4, stringsAsFactors= FALSE)

sampleNames(opis)=paste(sampleNames(opis), ".CEL", sep="")

data_Affy=ReadAffy(filenames=sampleNames(opis), verbose=T)
data_Affy@cdfName=paste("ga",data_Affy@cdfName,sep="")
s