setwd("~/Dropbox (CRG)/Personal_Estefania/Network/")
getwd()
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
psi<-read.table("compatibledata/GeneEventsKds.tab", stringsAsFactors = F)
dim(psi)
#pisa la original
#####################compute deltaPSI###############
####################################################
#Knock Down (Column) Scaling. 
#  Column Scaling -> Only Shape not Scale of KD effect is important. This is essential since for example not all KDs have the same efficiecny 
hclust(psi)
deltascaled<-scale(psi)    #scaled by KDs, maintain the profiles
#Calculate Robust Correlation Matrix
C <- CRobCor(deltascaled)
dim(C)
gList<-CreateGraph(C,rho=0.75)   # Rho was set to achieve FDR < 5%
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
