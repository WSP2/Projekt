# czarne, długie, kręcone zęby
# ścieżka Projekt/ i do każdego folderu blebleble
# ftp://157.158.14.226 wsp wsp2015
source('http://www.bioconductor.org/biocLite.R')
# biocLite('Biobase')
# library(Biobase)
# biocLite('affy')
# library(affy)
# biocLite('annotate')
# library(annotate)
# biocLite('hgu95av2.db')
# library(hgu95av2.db)
# biocLite('gahgu95av2.db')
# library(gahgu95av2.db)
# biocLite('gcrma')
# library('gcrma')

# d1=read.affybatch('Main Project/Cells/CL2001031606AA.CEL')
# # image(d1[,1])
# # hist(d1,main="bleble")
# # MAplot(d1)
# d1r=rma(d1)
# d1m=mas5(d1)
# # ai<-compute.affinities(cdfName(d1))
# # d1g=gcrma(d1,affinity.info=ai,type="affinities")



# szki=read.table('http://www.broadinstitute.org/mpr/publications/projects/LUNG/datasetA_scans.txt',header=T
#                 ,sep ="\t")

dataDirectory <- system.file("extdata", package="Biobase")
exprsFile <- file.path("http://www.broadinstitute.org/mpr/publications/projects/LUNG/datasetA_scans.txt")
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t",as.is=TRUE))

c=exprs[exprs[,2]=='CARCINOID',]; c1=c[sample(nrow(c),15),] # wybranie losowych 15 z klasy CARCINOID
# c=szki[szki$CLASS=='CARCINOID',]; c1=c[sample(nrow(c),15),] # wybranie losowych 15 z klasy CARCINOID
paste(sprintf("%s",c1[,4]),".CEL",sep="")
# d1=read.affybatch(paste("Main Project/Cells/",sprintf("%s",c1[,4]),".CEL",sep=""))
d=read.affybatch(paste("C:/Users/Anna/Documents/WSP/",
                        sprintf("%s",c1[,4]),".CEL",sep=""),verbose=T)
drma=rma(d)
de=exprs(drma)
