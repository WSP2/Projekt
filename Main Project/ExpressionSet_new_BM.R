source("http://bioconductor.org/biocLite.R")
#biocLite("BioUpgrade")
library(Biobase)
#biocLite('affy')
library('affy')
#biocLite('gahgu95av2.db')
library(gahgu95av2.db)

#install.packages(pkgs='gplots')
library(gplots)
#biocLite("convert")
library(convert)

###wczytanie plików .CEL -> obiekt AffyBatch
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

#normalizacjia RMA
RMA=rma(data_Affy)
dataRMA=exprs(RMA)


#opis eksperymentu
experiment = new("MIAME", name="Arindam Bhattacharjee1,William G. Richards,Jane Staunton,Cheng Li, Stefano Monti, Priya Vasa, Christine Ladd, Javad Beheshti, Raphael Bueno, Michael Gillette, Massimo Loda, Griffin Weber, Eugene J. Mark, Eric S. Lander, Wing Wong, Bruce E. Johnson, Todd R. Golub, David J. Sugarbaker,Matthew Meyerson1", 
                  title = "Classification of Human Lung Carcinomas by mRNA Expression Profiling Reveals Distinct Adenocarcinoma Sub-classes",
                  lab="1.Departments of Adult Oncology and Pediatric Oncology, Dana–Farber Cancer Institute, Harvard Medical School, 2.Department of Biostatistics, Harvard School of Public Health, Boston, 3. Departments of Surgery and Pathology, Brigham and Women's Hospital, Boston, 4.Whitehead Institute/Massachusetts Institute of Technology Center for Genome Research, Cambridge, 5. Department of Pathology, Massachusetts General Hospital",
                  pubMedIds="PMID:11707567[PubMed - indexed for MEDLINE], PMCID: PMC61120",
                  samples=list(data$Sample),
                  contact = "Arindam_Bhattacharjee@dfci.harvard.edu, staunton@genome.wi.mit.edu, Matthew_Meyerson@dfci.harvard.edu,golub@genome.wi.mit.edu",
                  abstract = "We have generated a molecular taxonomy of lung carcinoma, the leading cause of cancer death in the United States and worldwide. Using oligonucleotide microarrays, we analyzed mRNA expression levels corresponding to 12,600 transcript sequences in 186 lung tumor samples, including 139 adenocarcinomas resected from the lung. Hierarchical and probabilistic clustering of expression data defined distinct sub-classes of lung adenocarcinoma. Among these were tumors with high relative expression of neuroendocrine genes and of type II pneumocyte genes, respectively. Retrospective analysis revealed a less favorable outcome for the adenocarcinomas with neuroendocrine gene expression. The diagnostic potential of expression profiling is emphasized by its ability to discriminate primary lung adenocarcinomas from metastases of extra-pulmonary origin. These results suggest that integration of expression profile data with clinical parameters could aid in diagnosis of lung cancer patients",
                  url = "www.broadinstitute.org/mpr/lung/",other = list(notes="Publication date: 11/13/2001","Publication URL: http://www.pnas.org/cgi/content/full/191502998v1"))

#obiekt ExpressionSet
ExprSet = new("ExpressionSet", expr = dataRMA, phenoData = opis,
              experimentData = experiment, annotation="gahgu95av2.db")

#usunięcie 5% sond o najwyższej i najliższej średniej ekspresji
expr_sort = sort(rowMeans(exprs(ExprSet)),index.return=T)
feat_num = dim(ExprSet)[1]
cutoff = round(dim(ExprSet)[1]*0.025)
ind_clear = expr_sort$ix[c(1:cutoff,(feat_num-cutoff):feat_num)]
ExprSet = ExprSet[-ind_clear,]
View(ExprSet)

