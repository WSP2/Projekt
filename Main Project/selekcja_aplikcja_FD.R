# fun. selekcja i apl. okienkowa (Filipczyk, Dyduch)
require(tcltk) || stop("tcltk support is absent")

library('Biobase')
library('affy')
biocLite('estrogen')
library('estrogen')
biocLite('hgu95av2')
library('hgu95av2')
biocLite('hgu95av2.db')
library('hgu95av2.db')
biocLite('hgu95av2cdf')
library('hgu95av2cdf')
biocLite('gahgu95av2.db')
library('gahgu95av2.db')
biocLite('heatmap.plus')
library('heatmap.plus')
biocLite('vcd')
library('vcd')

selekcja <- function(obiekt, korekcja, wybor_sond, ilosc_sond, kolumna_sortujaca, separator){

  ekspresja_carc= 2^(exprs(obiekt)[,which(obiekt@phenoData@data$CLASS=='CARCINOID')])
  ekspresja_normal = 2^(exprs(obiekt)[,which(obiekt@phenoData@data$CLASS=='NORMAL')])
  
  statystyka_t = 0
  p_wartosc = 0
  r = dim(obiekt@assayData$exprs)[1]
  
  for (i in 1:r){
    test = t.test(ekspresja_carc[i,], ekspresja_normal[i,]) #test t dla grup carc i normal (kazda sonda osobno)
    statystyka_t[i] = test$statistic
    p_wartosc[i] = test$p.value
  }
  p_wartosc_skorygowana = p.adjust(p_wartosc, korekcja) #korekcja p-wartosci wskazanc metode
  
  #annotacja symboli genow
  identyfikator_ferrari = featureNames(obiekt)
  annotacje_symboli <- gahgu95av2SYMBOL 
  mapa_symboli <- mappedkeys(annotacje_symboli)
  lista_symboli <- as.list(annotacje_symboli[mapa_symboli])
  annotowane_symbole = lista_symboli[identyfikator_ferrari]
  
  #annotacja nazwy genow
  annotacje_nazw <- gahgu95av2GENENAME
  mapa_nazw <- mappedkeys(annotacje_nazw)
  lista_nazw <- as.list(annotacje_nazw[mapa_nazw])
  annotowane_nazwy = lista_nazw[identyfikator_ferrari]

  if (wybor_sond==1) { #wybor sond 1. - funkcja zwraca zadana liczbe sond
    posortowane_sondy = sort(p_wartosc_skorygowana, index.return=TRUE)
    posortowane_probki = posortowane_sondy$ix[1:ilosc_sond]
    
  } else { #wyb?r sond 2. - funkcja zwraca sondy spod progu
    posortowane_probki = which(p_wartosc_skorygowana < ilosc_sond)
  }
  
  FerrariID = identyfikator_ferrari[posortowane_probki]
  wartosc_statystyki_t = statystyka_t[posortowane_probki]
  wartosc_nieskorygowana_p = p_wartosc[posortowane_probki]
  wartosc_skorygowana_p = p_wartosc_skorygowana[posortowane_probki]
  srednia_carc = as.vector(rowMeans(ekspresja_carc)[posortowane_probki])
  srednia_normal = as.vector(rowMeans(ekspresja_normal)[posortowane_probki])
  FoldChange = srednia_carc/srednia_normal
  
  GeneSymbol = annotowane_symbole[posortowane_probki]
  GeneSymbol[unlist(lapply(GeneSymbol, is.null))] = NA
  GeneSymbol = as.vector(unlist(GeneSymbol))
  
  GeneName = annotowane_nazwy[posortowane_probki]
  GeneName[unlist(lapply(GeneName, is.null))] = NA
  GeneName = as.vector(unlist(GeneName))
  
  tabelka = data.frame(FerrariID,wartosc_statystyki_t,wartosc_nieskorygowana_p,wartosc_skorygowana_p,srednia_carc,srednia_normal,FoldChange,GeneSymbol,GeneName)
  colnames(tabelka) = c('FerrariID', 'wartosc statystyki t', 'p-wartosc', 'skorygowana p-wartosc', 'srednia CARC', 'srednia NORMAL', 'FoldChange', 'GeneSymbol', 'GeneName')
  sortowanie_tablicy = tabelka[kolumna_sortujaca] #sortowanie tablicy wg zaqdanej kolumny
  posortowane_sondy_tabelka = sort(as.matrix(sortowanie_tablicy), index.return=TRUE)
  tabelka_posortowana = tabelka[posortowane_sondy_tabelka$ix,]
  
  write.table(tabelka_posortowana, file='wyniki.txt', sep=separator, row.names=FALSE, col.names=TRUE, quote=FALSE) #zapis tabelki do pliku tekstowego
  
  png('heatmapa.png') #zapis heatmapy do pliku
  kolory = ifelse(obiekt@phenoData@data$CLASS=='CARC', 'blue', 'red')
  heatmap(obiekt@assayData$exprs[posortowane_probki,], ColSideColors=kolory)
  legend('topleft',legend=c('CARC', 'NORMAL'), fill=c('blue', 'red'))
  dev.off()
  
  
  return(tabelka_posortowana)
  
}

aplikacja <- function(){

  
  korekcja_var <- tclVar("")
  wybor_sond_var <- tclVar("")
  ilosc_sond_var <- tclVar("")
  kolumna_sortujaca_var <- tclVar("")
  separator_var <- tclVar("")
  
  tt <- tktoplevel()
  tkwm.title(tt,"Aplikacja okienkowa")
  korekcja.entry <- tkentry(tt, textvariable=korekcja_var)
  wybor_sond.entry <- tkentry(tt, textvariable=wybor_sond_var)
  ilosc_sond.entry <- tkentry(tt, textvariable=ilosc_sond_var)
  kolumna_sortujaca.entry <- tkentry(tt, textvariable=kolumna_sortujaca_var)
  separator.entry <- tkentry(tt, textvariable=separator_var)
  
  reset <- function() {
    tclvalue(korekcja_var)<-""
    tclvalue(wybor_sond_var)<-""
    tclvalue(ilosc_sond_var)<-""
    tclvalue(kolumna_sortujaca_var)<-""
    tclvalue(separator_var)<-""
  }
  
  reset.but <- tkbutton(tt, text="Wyczysc", command=reset)
  
  submit <- function() {
    korekcja <- (tclvalue(korekcja_var))
    wybor_sond <- as.numeric(tclvalue(wybor_sond_var))
    ilosc_sond <- as.numeric(tclvalue(ilosc_sond_var))
    kolumna_sortujaca <- as.numeric(tclvalue(kolumna_sortujaca_var))
    separator_wybor <- as.numeric(tclvalue(separator_var))
    
    if (separator_wybor==1){
      separator = '\t'
    } else {
      if (separator_wybor==2){
        separator = ';'
      } else {
        separator = ' '
      }
    }
    
    nazwa_obiektu <- tclvalue(tkgetOpenFile())
    load(nazwa_obiektu)
    
    selekcja(obiekt, korekcja, wybor_sond, ilosc_sond, kolumna_sortujaca, separator)
  }
  
  
  submit.but <- tkbutton(tt, text="Wybierz obiekt i licz", command=submit)
  
  quit.but <- tkbutton(tt, text = "Zamknij R", 
                       command = function() {
                         q(save = "no")
                         tkdestroy(tt)
                       }
  )
  
  tkgrid(tklabel(tt,text="\n\n"))
  tkgrid(tklabel(tt,text="Wprowadz dane do obliczen"), columnspan=3)
  tkgrid(tklabel(tt,text="\n\n"))
  tkgrid(tklabel(tt,text="Korekcja p-value: BH, bonferroni, BY, fdr, hochberg, holm, hommel, none"), korekcja.entry)
  tkgrid(tklabel(tt,text="Wybor sond: 1 - podana ilosc, 2 - ponizej progu"), wybor_sond.entry)
  tkgrid(tklabel(tt,text="Ilosc sond lub prog odciecia"), ilosc_sond.entry)
  tkgrid(tklabel(tt,text="Numer kolumny do sortowania"), kolumna_sortujaca.entry)
  tkgrid(tklabel(tt,text="Separator: 1 - tab, 2 - srednik, 3 - spacja"), separator.entry)
  tkgrid(tklabel(tt,text="\n\n"))
  tkgrid(submit.but, columnspan=3)
  tkgrid(tklabel(tt,text="\n"))
  tkgrid(reset.but, columnspan=3)
  tkgrid(quit.but, columnspan=3)
  tkgrid(tklabel(tt,text="\n\n"))
  

  tkwait.window(tt)
  
}

aplikacja()

