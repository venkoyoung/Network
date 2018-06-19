setwd("~/Dropbox (CRG)/Personal_Estefania/Network/comparison_old_new/")
source("scripts/CentralityRanking.R")
source("scripts/CreateGraph.R")
source("scripts/CRobCor.R")
source("scripts/Vscale.R")
source("scripts/CreateGraph.R")
source("original_createGraph.R")
setwd("rnaseq/")
##################################
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
#######################################
########## Functions
#Calculation of custom robust correlation. This tries to discriminate between technical outliers (artefacts) and true biological outliers.
#See supplemental methods of Papasaikas, Tejedor et al for details.
#Arguments are M -> a n x p matrix (n features -e.g ASEs- in the rows, p variables -e.g KDs- in the columns)
#and ct a scalar defining the threshold for finding close neighbors of variables
############ DATA FILTERING AND PREPROCESSING:
#original table:
#viene del renames
geneEventsKdsNoDup<-read.table("profiles/GeneEventsKds.tab")
delta<-geneEventsKdsNoDup
#remove controls from matrix before scaling
head(delta); dim(delta)
deltascaled<-scale(delta)    #scaled by KDs, maintain the profiles
eventsscaled<-t(scale(t(deltascaled)))
dim(eventsscaled)
dim(deltascaled)
#scaled by KDs, maintain the profiles)
#Calculate Robust Correlation Matrix
getwd()
C <- CRobCor(deltascaled)
pdf("byKDs_corrPlot_GE_alphabet.pdf")
corrplot(C, method = "color", order = "alphabet", tl.cex = 0.1)
dev.off()
pdf("byKds_corrPlot_ge_hclust.pdf")
corrplot(C, method = "color", order = "hclust",  col = brewer.pal(n = 8, name = "RdBu"), tl.cex = 0.1)
dev.off()

Cdouble <- CRobCor(eventsscaled)
dim(eventsscaled)
write.table(eventsscaled, "ge_evenstscaled.tab")
getwd()
pdf("Double_corrPlot_GE_alphabet.pdf")
corrplot(Cdouble, method = "color", order = "alphabet", tl.cex = 0.1)
dev.off()
setwd("ge/")
pdf("Double_corrPlot_GE_hclust.pdf")
corrplot(Cdouble, method = "color", order = "hclust",  col = brewer.pal(n = 8, name = "RdBu"), tl.cex = 0.1)
dev.off()
write.table(eventsscaled,"profiles/gene_double_scaled.tab", sep="\t")

#columns=events
#Calculate Robust Correlation Matrix
save(C, file="Original_Events_GE_C_singlescaled.RData")
write.table(C, file="Original_Events_GE_C_singlescaled.tab", sep="\t", col.names = NA)
save(Cdouble, file="Original_Events_GE_C_doublescaled.RData")
write.table(Cdouble, file="Original_Events_GE_C_doublescaled.tab", sep="\t", col.names = NA)
#Create graph, calculate centrality measures, identify communities and plot
gListDouble<-CreateGraph(Cdouble,rho=0.65) 
gList<-gListDouble
# Rho was set to achieve FDR < 5%
##############################################
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
save(gList, file="gList_6_GE-Seq.RData")
getwd()
write_graph(g, format="edgelist", file="GE")
##############################################################
#exporta parameters
deg.dist
?degree_distribution
deg.dist <- degree_distribution(g, cumulative=T, mode="all")
png("cumulativeFreq_GE_065.png")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
  xlab="Degree", ylab="Cumulative Frequency")
dev.off()
########################################3
write.table(t(values), file="GE_06topology.tab", sep="\t")
#add community information:
comps <-memb$membership; head(comps)#dpmde esta cada gen, em que grupo
table(comps)
write.table(comps, file="GE_06memberships.txt", col.names=NA)
########################################
g <- make_ring(10)
degree(g)
g2 <- sample_gnp(1000, 10/1000)
degree_distribution(g2)
deg2 <- degree(g2, mode="all")
mean(deg2)
deg.dist2 <- degree_distribution(g2, cumulative=T, mode="all")
plot( x=0:max(deg2), y=1-deg.dist2, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")

