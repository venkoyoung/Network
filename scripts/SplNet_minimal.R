setwd("~/Dropbox (CRG)/Personal_Estefania/Network/")
##################################
source("https://bioconductor.org/biocLite.R")
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
#######################################
########## Functions
source("scripts/CentralityRanking.R")
source("CreateGraph.R")
source("scripts/CRobCor.R")
source("scripts/Vscale.R")
#Calculation of custom robust correlation. This tries to discriminate between technical outliers (artefacts) and true biological outliers.
#See supplemental methods of Papasaikas, Tejedor et al for details.
#Arguments are M -> a n x p matrix (n features -e.g ASEs- in the rows, p variables -e.g KDs- in the columns)
#and ct a scalar defining the threshold for finding close neighbors of variables
############ DATA FILTERING AND PREPROCESSING:
#original table:
#########################################################################
#read.data
setwd("~/Dropbox (CRG)/Personal_Estefania/Network/")
source("scripts/CentralityRanking.R")
source("scripts/CreateGraph.R")
source("scripts/CRobCor.R")
source("scripts/Vscale.R")
#########################################################################library("impute") #Availabe from bioconductor
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
library(impute)
#be careful because original data is trasposed respect to this data
load("compatibledata/psifinalOnlyLabChipKD.RData")
dim(psifinal)
setwd("comparison_old_new/")
psi<-read.table("compatibledata/psi_only_labChipKDS.tab", stringsAsFactors = F, 
                colClasses = c("character",rep("numeric", 271)), 
                row.names = 1)
dim(psi)
SparseRows<-which(rowMeans(is.na(psi))*100>30)#Row numbers with >30% missing values
head(SparseRows)
#5 columns  
SparseCols<-which(colMeans(is.na(psi))*100>30)#Row numbers with >30% missing values
psiCleaned<-psi[-c(SparseRows), -c(SparseCols)]
class(psiCleaned)#26,266
#pisa la original
tpsi<-t(psiCleaned); dim(tpsi)
psiI<-impute.knn(as.matrix(tpsi),
                 k=round(sqrt(nrow(tpsi))*0.25),
                 rowmax=0.5)$data
hist(rowMeans(is.na(psiI))*100)# desaparecen los NAs
#####################compute deltaPSI###############
psifinal=t(psiI); 
head(psifinal); dim(psifinal); psifinal[1:5,260:266]
controls<-rowMeans(psifinal[,260:266])
sdcontrols<-rowSds(psifinal[,260:266])
delta<-psifinal-controls
head(delta); dim(delta)
####################################################
delta[1:5,1:10]
dim(delta)#26 266
#Knock Down (Column) Scaling. 
#  Column Scaling -> Only Shape not Scale of KD effect is important. This is essential since for example not all KDs have the same efficiecny 
deltascaled<-scale(delta)    #scaled by KDs, maintain the profiles
#Calculate Robust Correlation Matrix
C <- CRobCor(deltascaled)
save(C, file="Original_Events_RNASeq_C.RData")
write.table(C, file="Original_Events_RNASeq_C.tab", sep="\t", col.names = NA)
load("Original_Events_RNASeq_C.RData")
#Create graph, calculate centrality measures, identify communities and plot
source("original_createGraph.R")
getwd()
gList<-CreateGraph(C,rho=0.65)   # Rho was set to achieve FDR < 5%
g<-gList[[1]]
memb<-gList[[2]]
values<-gList[[3]]
mean(values[1,])#6.1  
V(g)
E(g)
deg <- degree(g, mode="all")
mean(deg)
hist(values[1,])#6.1  
summary(values[1,])
getwd()
save(gList, file="gList_06_26-266-RNA-Seq.RData")
##############################################################
#exporta parameters
deg.dist <- degree_distribution(g, cumulative=T, mode="all")
png("cumulativeFreq_RNA-Seq_06_26-266.png")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
  xlab="Degree", ylab="Cumulative Frequency")
dev.off()
########################################3
write.table(t(values), file="RNA_Seq_06_26-266topology.tab", sep="\t")
#add community information:
comps <-memb$membership; head(comps)#dpmde esta cada gen, em que grupo
table(comps)
write.table(comps, file="RNA_Seq_06_26-266memberships.txt", col.names=NA)
########################################
#match commnuities/descriptores
topo<-t(values)
topo
names(comps)
row.names(topo)[1]
names(comps)[1]
both<-row.names(topo)%in%names(comps)
ii<-match(row.names(topo),names(comps))
names<-names(comps)[ii]
commu<-comps[ii]
names[is.na(names)]<-"notCommu"
commu[is.na(commu)]<-0
row.names(topo)
write.table(data.frame(topo, names, commu), file="RNA_Seq_allInfo062.tab", sep="\t")
########################################
#######################compare original kds 
