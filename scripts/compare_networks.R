library("impute") #Availabe from bioconductor
library("glasso")
library("QUIC")
library("RColorBrewer")
library("igraph")
library("SDMTools")
library("scales")
library("png") # for reading in PNGs
library(data.table)
library(stringr)
library(matrixStats)
library(corrplot)
setwd("/home/emancini/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/")
source("scripts/CentralityRanking.R")
source("scripts/CreateGraph.R")
source("scripts/CRobCor.R")
source("scripts/Vscale.R")
source("scripts/CreateGraph.R")
source("scripts/original_createGraph.R")
#comparison networks:
orig_events<-
read.table("/home/emancini/Dropbox (CRG ADV)/Personal_Estefania/GosiaAndEstefi/AS_317_OK/Original_Network_events.tsv", header=T, sep="\t")
head(orig_events); dim(orig_events)
ii<-match(orig_events$EVENT, rownames(data))
dfOrigEvents<-data[ii,]; head(dfOrigEvents)
head(data)
eventsIDs<-paste(rownames(data)[ii],
                 1:length(rownames(data)[ii]), sep="_")
rownames(dfOrigEvents)<-paste(data$GENE[ii], eventsIDs, sep=":")
onlyValues<-dfOrigEvents[,7:ncol(dfOrigEvents)]

#labchips
t_Zvals<-read.table(file="/home/emancini/Dropbox (CRG ADV)/Personal_Estefania/Network/comparison_old_new/profiles/newNames_t_Zvals.tab", sep="\t")#OK
sampleData <- as.matrix(t_Zvals)
###################################################3
head(sampleData); dim(sampleData)#35*267
originalEvents<-rownames(sampleData)
write(colnames(sampleData))
write(originalEvents)
getwd()
#######################################################
CLabChips <- CRobCor(sampleData)
gList<-CreateGraph(CLabChips,rho=0.48)   # Rho was set to achieve FDR < 5%
g<-gList[[1]]
memb<-gList[[2]]
values<-gList[[3]]
mean(values[1,])#4.09
deg <- degree(g, mode="all"); mean(deg)
deg.dist <- degree_distribution(g, cumulative=T, mode="all")
V(g);E(g);
df<-get.edgelist(g); 
colnames(df)<-c("source", "target")
write.table(df, file="edgelist_048_labchip.tab", col.names = NA, quote = F, sep="\t")
#196 KNOCK DOWNS
#################################################################
############################
############################
psi<-read.table( "/home/emancini/Dropbox (CRG ADV)/Personal_Estefania/Network/comparison_old_new/profiles/psi_only_labChipKDS.tab",  sep="\t")
psi[1:5,1:5]
dim(psi)
write(colnames(psi))
SparseRows<-which(rowMeans(is.na(psi))*100>30)#Row numbers with >30% missing values
length(SparseRows)#5 columns  
SparseCols<-which(colMeans(is.na(psi))*100>30)#Row numbers with >30% missing values
length(SparseCols)
psiCleaned<-psi[-c(SparseRows), -c(SparseCols)]
dim(psiCleaned)#26,266
#pisa la original
write(colnames(psiCleaned))
tpsi<-t(psiCleaned)
psiI<-impute.knn(as.matrix(tpsi),
                 k=round(sqrt(nrow(tpsi))*0.25),
                 rowmax=0.5)$data
psiI
psifinal=t(psiI)
psifinal[1:5,grep("^AA", colnames(psifinal))]
controls<-rowMeans(psifinal[,grep("^AA", colnames(psifinal))])
delta<-psifinal-controls
#remove controls from matrix before scaling
delta<-delta[,-grep("^AA", colnames(delta)) ]
head(delta); dim(delta);delta[1,1:10]
#Knock Down (Column) Scaling. 
#  Column Scaling -> Only Shape not Scale of KD effect is important. This is essential since for example not all KDs have the same efficiecny 
deltascaled<-scale(delta)    #scaled by KDs, maintain the profiles
eventsscaled<-t(scale(t(deltascaled)))
write.table(eventsscaled,"eventscaled.tab", sep="\t")
Cdouble <- CRobCor(eventsscaled)
gListDouble<-CreateGraph(Cdouble,rho=0.55) 
gList<-gListDouble
# Rho was set to achieve FDR < 5%
##############################################
g<-gList[[1]]
memb<-gList[[2]]
values<-gList[[3]]
mean(values[1,])#1.71
V(g)#200
E(g)#397
deg <- degree(g, mode="all")
mean(deg)#3.97
hist(values[1,])#
summary(values[1,])
df<-get.edgelist(g); 
colnames(df)<-c("source", "target")
write.table(df, file="edgelist_055_rnaseq.tab", col.names = NA, quote = F, sep="\t")
######################################################

