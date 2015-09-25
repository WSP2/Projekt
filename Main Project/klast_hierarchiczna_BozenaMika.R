install.packages("cluster", dependencies = TRUE)
library(cluster)

#wybór metody pomiaru odleglości dla klasteryzacji hierarchicznej- podaje użytkownik
fun <- function(metody) {
  m=as.numeric(readline("distance measure: 
      1.euclidean,  
      2.maximum, 
      3.manhattan, 
      4.canberra, 
      5.binary, 
      6.minkowski      
                        give the number:"))
  met=metody[m]
  return(met)
  
}
metody=c("euclidean","maximum","manhattan","canberra","binary","minkowski")

metoda_odl=fun(metody)

# agglomeration method -podaje użutkownik
fun2 <-function(metody) {
  m=as.numeric(readline("distance measure: 
  1.ward.D, 
  2.ward.D2, 
  3.single, 
  4.complete, 
  5.average, 
  6.mcquitty,
  7.median,
  8.centroid,  
                  give the number:"))  #podaje użytkownik
  met=metody[m]
  return(met)
  
}
ag=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")

metoda_hier=fun2(ag)

#użytkownik podaje czy chce sam dokonać podziału na klastry i zaznaczyć je na dendrogramie
c=as.numeric(readline("do you want to give the number of clusters? 1-Yes, 0-NO:"))

# podział na ile klastrów?- podaje użytkownik 
kl=as.numeric(readline("divide for (give number) clusters:"))

#użytkownik podaje czy chce poziomy dendrogram
den=readline('rotates plot by 90 degrees T/F :') 


#klasteryzacja
kl_hierarch <-function(ExprSet,metoda_odl,metoda_hier,c,kl,den) {

#pobranie danych 
dane=ExprSet@assayData$exprs 
dane=aperm(exprs(ExprSet))
#pobranie przynależności do klasy
sample=ExprSet@phenoData@data$Sample
KLASY=ExprSet@phenoData@data$CLASS


d=dist(dane, method = metoda_odl) #pomiar odleglości
fit=hclust(d, method=metoda_hier)
ka=as.dendrogram(fit)

tiff(filename="dendrogram_R.tif",compression='none')

# dendogram - plot
nP=list(col = c(1,NA), cex = c(0.7, 0.8), pch =21:22,
    bg =  c("light blue", "cyan3"),lab.cex = 0.75, lab.col = "black")

plot(ka,nodePar=nP,edgePar=list(col="azure4", lwd = 1.2),horiz=den)

#podział na kl klastrów
if (c==1){
rect.hclust(fit,k=2, border=c(2:(kl+2)))  #wyrysowanie na dendogramie kl klastrów

clast=as.data.frame(cutree(fit, k=kl)) #podział próbek na kl klastrów- podsumowanie
pods=as.data.frame(cbind(sample,KLASY,clast)) #tabela wyników klasteryzacji
colnames(pods)=c('Sample','prawidłowa klasa','klaster')
dev.off()
return(pods)
} else{
  dev.off()
}

}

##gdy użytkownik podaje, że chce sam podać ilość klastrów kl to results zawiera 
##tabele podsumowującą; jeśli c==0 - brak tabeli, tylko sam dendrogram jest zapisywany do pliku

results=kl_hierarch(ExprSet,metoda_odl,metoda_hier,c,kl,den)

