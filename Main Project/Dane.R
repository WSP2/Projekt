# TODO: Add comment
# 
# Author: Filip
##############################################################################
#install.packages("R.matlab")
library(R.matlab)
C<-getwd()
C=sub("Main Project", "Data", C)
pathname <- file.path(C, "genes_list.mat")
data <- readMat(pathname) 

