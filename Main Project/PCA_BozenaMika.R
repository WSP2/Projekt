##podaje użytkownik
odp=readline('log transformation T/F :') #czy dane mają być w skali logarytmicznej?

PCA<- function(ExprSet,odp) {
 
dane=exprs(ExprSet)
dane2=t(dane)

pca2=prcomp(dane2,retx=odp,scale=T) #PCA

jpeg(filename="PCA_ScreePlot.jpg")
plot(pca2,main="PCA")
dev.off()

jpeg(filename="PCA_ScreePlot2.jpg")
plot(pca2,type='l',main="PCA")
dev.off()

su=summary(pca2) 

scores=pca2$x

jpeg(filename="PCA_new.jpg")
color=ifelse(pData(ExprSet)$CLASS=='CARCINOID','green','red')
plot(scores[,1:2],col=color, main="PCA")
dev.off()

return(su)
  
}

var=PCA(ExprSet,odp)  ##wykresy są zapisane w folderze, a zmienna var zawiera podsumowanie
