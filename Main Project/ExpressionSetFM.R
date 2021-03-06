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
pathname <- file.path(C, "datasetA_scans.txt")
data = read.table(pathname, header = TRUE,sep = "\t" )
data = data[c(191:193, 211:213),]
opis = read.AnnotatedDataFrame(pathname, sep = "\t", header = TRUE, row.names = 4, stringsAsFactors = F)
opis = opis[c(191:193, 211:213),]
sampleNames(opis) = paste(sampleNames(opis), ".CEL", sep = "")
data_Affy = ReadAffy(filenames=sampleNames(opis), verbose = TRUE)
data_Affy@cdfName = paste("ga",data_Affy@cdfName, sep="")
RMA = rma(data_Affy)
dataRMA = exprs(RMA)
experiment = new("MIAME", name="Dane mikromacierzowe", lab ="IO", title = "dane testowe",
                 abstract = "Przyklad", url = "http://www.bioconductor.org",
                 other = list(notes = "inne"))
ExprSet = new("ExpressionSet", expr = dataRMA, phenoData = opis,
              experimentData = experiment, annotation="gahgu95av2.db")
expr_sort = sort(rowMeans(exprs(ExprSet)),index.return=T)
feat_num = dim(ExprSet)[1]
cutoff = round(dim(ExprSet)[1]*0.025)
ind_clear = expr_sort$ix[c(1:cutoff,(feat_num-cutoff):feat_num)]
ExprSet = ExprSet[-ind_clear,]
View(ExprSet)
