source("http://bioconductor.org/biocLite.R")
biocLite()
library(Biobase)
biocLite('affy')
library(affy)

biocLite('gahgu95av2.db')

library(gahgu95av2.db)
install.packages(pkgs='gplots')
library(gplots)
biocLite("convert")
library(convert)

dane=read.table("datasetA_scans.txt",header = TRUE, sep="\t")
dane=dane[c(191:194,211:214),]


opis=read.AnnotatedDataFrame("datasetA_scans.txt",sep="\t",header = TRUE,row.names=4,stringsAsFactors=F)
opis=opis[c(191:194,211:214),]
sampleNames(opis)=paste(sampleNames(opis),".CEL",sep="")
data_Affy=ReadAffy(filenames=sampleNames(opis),verbose=T)
data_Affy@cdfName=paste("ga",data_Affy@cdfName,sep="")


#normalizacjia
RMA=rma(data_Affy)
dataRMA=exprs(RMA)

experiment=new("MIAME",name="Dane mikromacierzowe",lab="IO",title="dane testowe")