
#testy_stat<- function(ExprSet,t){

#wyodr?bnienie dw?ch grup - potrzebne do testu t. 
# sprawdzenie przynależnosci do konkretnego klastra i zapisywanie do 
#osobnej macierzy.

 is.matrix(results)
  quest =as.numeric(readline("Aby wykonać testy powinna zostać przeprowadzona klasteryzacja? Czy chcesz ją zrobić teraz? T/F"))
  if (quest==T)
    
  #carcinoid
  dlugosc = c(length(results[,1]))
  dlugosc_kol = c(length(results[1,]))
  
 dane2_t = matrix(0,15,3)
 dane2_t2 = matrix(0,15,3)
 
   for (i in 1:dlugosc){ 
     for (j in 1:dlugosc_kol){
    if (results[i, 3]==1)
  dane2_t[i,j]= matrix(results[i,j])
  #normal
  else
    dane2_t2[1:15,j] = matrix(results[i,j])
  
     }
}
#przypisanie podzielonych klastrow do danych z RMA
 #dalam tu na sztywno, nie wiem jak to przypisac, zeby 
 # pobieralo mi wynik z klastra do RMAdata, moze ktos cos ogarnie?
 
dataRMA_cancer= data2RMA[1:15,]
dataRMA_normal=data2RMA[16:30,]

#prawilny test t dla prob niezaleznych
  var.test(dataRMA_cancer,dataRMA_normal) #sprawdzmy czy jest roznorodnosc wariancji
  t = t.test(dataRMA_cancer,dataRMA_normal) # to chyba trzeba by by?o zzapisa? do jakiego? pliku
  
#wilcoxon? robic?

#heatmapa
#zapisywanie do jpg-a
jpeg(filename="heatmap.jpg")
hmp= heatmap(dataRMA_cancer,dataRMA_normal)
dev.off()
