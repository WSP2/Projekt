
#testy_stat<- function(ExprSet,t){

#wyodr?bnienie dw?ch grup

  dataRMA = t(dataRMA)
  #carcinoid
 dane2_t = matrix(0,15,3)
   for (i in 1:30){ 
     for (j in 1:3){
    if (results[i, 3]==1)
  dane2_t[i,j]= results[i,j]
  #normal
  else
  dane2_t2[i,j] = results[i,j]
     }
}

#prawilny test t dla prob niezale?nych

  var.test(dane2_t,dane2_t2) #sprawd?my czy jest r??norodno?? wariancji
  t = t.test(dane2_t,dane2_t2)
  #wilcoxon?
  
wyniki = testy_stat(ExprSet, t)
#write.table(file='testy_statystyczne.txt', wyniki)

#heatmapa

jpeg(filename="heatmap_")
hmp= heatmap(dane2_t, dane2_t2, Rowv=NULL)
dev.off()
