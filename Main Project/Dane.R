# TODO: Add comment
# 
# Author: Filip
##############################################################################
#install.packages("R.matlab")
library(R.matlab)
pathname <- file.path("~/R/Projekt/Data", "genes_list.mat")
data <- readMat(pathname) 

