install.packages("cluster", dependencies = TRUE)
library(cluster)

#pobranie danych
dane=ExprSet@assayData$exprs 
dane=aperm(dane) 
#pobranie przynależności do klasy
sample=ExprSet@phenoData@data$Sample
KLASY=ExprSet@phenoData@data$CLASS

#### klasteryzacja hierarchiczna: na podstawie: http://www.statmethods.net/advstats/cluster.html

#wybór metody pomiaru odleglości
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

d=dist(dane, method = metoda_odl) #pomiar odleglości

# agglomeration method 
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

fit=hclust(d, method=metoda_hier)

# dendogram
plot(fit,cex =0.75,)

#różne opcje podziału na klastry:
#samodzielne zaznaczanie na dendogramie poddziału na klastry
x=identify(fit) 

# ile klastrów?- podaje użytkownik (dzięki uprzedmiemu samodzielnemu zaznaczeniu można sprawdzić jak ukladają się próbki i dostosować ilość klastrów)
kl=as.numeric(readline("divide for (give number) clusters:"))


#podział na kl klastrów
rect.hclust(fit, k=kl, border=c(2:(kl+1)))  #wyrysowanie na dendogramie kl klastrów

clast=as.data.frame(cutree(fit, k=kl)) #podział próbek na kl klastrów- podsumowanie
pods=as.data.frame(cbind(sample,KLASY,clast)) #tabela wyników klasteryzacji
colnames(pods)=c('Sample','prawidłowa klasa','klaster')


