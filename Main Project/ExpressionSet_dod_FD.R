# dodane (Filipczyk, Dyduch)
source("http://bioconductor.org/biocLite.R")
biocLite("BioUpgrade")
library(Biobase)
biocLite('affy')
library('affy')
biocLite('gahgu95av2.db')
library('gahgu95av2.db')
install.packages(pkgs='gplots')
library(gplots)
biocLite("convert")
library(convert)

###wczytanie plikoww .CEL -> obiekt AffyBatch
C<-getwd()
C=sub("Main Project", "Data", C)
pathname <- file.path(C, "datasetA_scans.txt")
data = read.table(pathname, header = T,sep = "\t" )
data = data[c(191:205, 211:225),]

###summary(data)

opis = read.AnnotatedDataFrame(pathname, sep = "\t", header = T, row.names = 4, stringsAsFactors = F)
opis = opis[c(191:205, 211:225),]
###look=head(pData(opis),n=nrow(opis)) #opis probek
sampleNames(opis) = paste(sampleNames(opis), ".CEL", sep = "")
data_Affy = ReadAffy(filenames=sampleNames(opis), verbose = T)
##annotacje Ferrariego

data_Affy@cdfName = paste("ga",data_Affy@cdfName, sep="")
data_Affy@annotation = paste("ga",data_Affy@annotation, sep="")

##################################################################
### normalizacja danych
dane_n = expresso(data_Affy, bgcorrect.method='rma', normalize.method='quantiles', pmcorrect.method='pmonly', summary.method='medianpolish') #normalizacja danych

#### dodane - stworzenie obiektu exprSet
assdat = dane_n@assayData$exprs #assayData
pdat = read.table('moje_dane.txt', row.names=4, header=TRUE, sep='\t') #pData
colnames(assdat)=rownames(pdat)
metdat = data.frame(labelDescription=c('Gene simple annotation', 'CARCINOID/NORMAL', 'Sample name'), row.names=names(pdat)) #metaData
phdat = new('AnnotatedDataFrame', data=pdat, varMetadata=metdat) #phenoData
expd = new('MIAME', name='Dane mikromacierzowe', lab='WSP', title='Porownanie carc i normal') #opis eksperymentu experimentDesc
obiekt = new('ExpressionSet', exprs=assdat, phenoData=phdat, experimentData=expd, annotation = 'gahgu95av2') #obiekt ExprSet
View(obiekt)

save(obiekt, file="obiekt_nf.RData") #zapisanie obiektu
rm(list=ls()) #wyczyszczenie workspace
load(file='obiekt_nf.RData') #wczytanie obiektu

##### usuniecie 5% (najw/najn) carc
ldof = round(dim(exprs(obiekt))[1] * 0.05) #5% liczby sond -> do odfiltrowania
posortowane_carc = sort(rowMeans(exprs(obiekt)[,1:15]), index.return=TRUE) #posortowanie rosnaco wg sredniej ekspresji poszczegolnych sond dla macierzy carc
exprs(obiekt) = exprs(obiekt)[posortowane_carc$ix,]
exprs(obiekt) = exprs(obiekt)[c((ldof+1):(dim(exprs(obiekt))[1]-ldof)),]

# usuniecie 5% (najw/najn) normal
ldof = round(dim(exprs(obiekt))[1] * 0.05) #5% liczby sond -> do odfiltrowania
posortowane_normal = sort(rowMeans(exprs(obiekt)[,16:30]), index.return=TRUE) #posortowanie rosnaco wg sredniej ekspresji poszczegolnych sond dla macierzy normal
exprs(obiekt) = exprs(obiekt)[posortowane_normal$ix,]
exprs(obiekt) = exprs(obiekt)[c((ldof+1):(dim(exprs(obiekt))[1]-ldof)),]

save(obiekt, file="obiekt.RData") # zapis obiektu
rm(list=ls()) #wyczyszczenie workspace
load(file='obiekt.RData') #wczytanie odfiltrowanego obiektu

View(obiekt)
#save(obiekt, file="ExprSet.RData")

# wykresy
analiza_pca = prcomp(exprs(obiekt))
klasy = obiekt@phenoData@data$CLASS
kolory = c('blue', 'pink')
plot(analiza_pca$x[,1], analiza_pca$x[,2], col=kolory[klasy], xlab='pc1', ylab='pc2',pch=20)
legend("topright", legend=levels(klasy), col=kolory, pch=20)
slupkowy = (analiza_pca$sdev[1:5]/sum(analiza_pca$sdev))*100
barplot(slupkowy, names.arg=c('1','2','3','4','5'), main='Udzial skladowych w calkowitej zmiennosci', xlab='Skladowe', ylab='Udzial [%]')



