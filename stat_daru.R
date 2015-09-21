#to jest tylko proba 

testy_stat<- function(ExprSet,t) {
  
#wyodr?bnienie dw?ch grup

  dataRMA = t(dataRMA)
  #carcinoid
  dane2_t=dataRMA[,(1:15)]
  #normal
  dane2_t2 = dataRMA[,(16:30)]
   
#prawilny test t dla prob niezale?nych
  var.test(dane2_t,dane2_t2) #sprawd?my czy jest r??norodno?? wariancji
  t = t.test(dane2_t,dane2_t2)
}
wyniki = testy_stat(ExprSet, t)
#write.table(file='testy_statystyczne.txt', wyniki)
#coś mi nie commituje