#first generate data: TRUE and RANDOM DATA
library(QUIC)
library(igraph)
#parameters
t_Zvals<-read.table("~/Dropbox (CRG ADV)/Personal_Estefania/Network/data/labchip_cleaned_idata.tab")
sampleData <- as.matrix(t_Zvals)
setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/data/")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/CRobCor.R")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/functions_fdr_noMC.R")
#FUNCTION:
estimateFDR(
  inputM=sampleData,
  NumRandomM=10,
  start=0.3,
  end=0.6,
  interval=0.05)
######################################  
#optimized in 0.5



