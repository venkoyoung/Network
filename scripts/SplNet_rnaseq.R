setwd("~/Dropbox (CRG)/Personal_Estefania/Network/comparison_old_new/")
source("../scripts/CentralityRanking.R")
source("../scripts/CreateGraph.R")
source("../scripts/CRobCor.R")
source("../scripts/Vscale.R")
source("../scripts/CreateGraph.R")
source("../original_createGraph.R")
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
gc()
all331<-read.table("~/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/AS_331/FcorrectedPsiValues_NA3_newNames_b.tab",
                       header = T, sep="\t", row.names = 1);dim(all331)
#remove bad samples
data <- all331[, ! colnames(all331) %in% 
                  c("AA2","AA1","CCDC12","C1orf55",
                    "C1orf55_b","CDC5L","HFM1","LENG1",
                    "RBM17","PPIL1","SRRT"),
                  drop = F]
#select the events:
orig_events<-read.table("~/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/AS_317_OK/Original_Network_events.tsv .csv", header=T, sep="\t")
  head(orig_events)
  ii<-match(orig_events$EVENT, rownames(data))
  dfOrigEvents<-data[ii,]; head(dfOrigEvents)
  head(data)
  eventsIDs<-paste(rownames(data)[ii],
 1:length(rownames(data)[ii]), sep="_")
  rownames(dfOrigEvents)<-paste(data$GENE[ii], eventsIDs, sep=":")
onlyValues<-dfOrigEvents[,7:ncol(dfOrigEvents)]

write.table(onlyValues,"onlyValuesOriginalEventsRNASEQ.txt", sep="\t")
getwd()
onlyValues<-read.table("../320onlyValuesOriginalEventsRNASEQ.txt")
psi<-onlyValues
#match the names of the KDs:
####################
colnames(onlyValues); dim(onlyValues)
length(which(onlyValues=="NAold"))#575
length(which(onlyValues=="NAnew3"))#833
length(which(onlyValues=="NAnew2"))#0
#############################################################
dim(onlyValues)[1]*dim(onlyValues)[2]#9920
dim(psiO)[1]*dim(psiO)[2]
575/9920*100#5.76
833/9920*100#8.39
(575+833)/9920*100#15%

psi[psi=="NAold"]<-NA
psi[psi=="NAnew3"]<-NA
#desde aca vengo con el script renames
psi<-psifinal
length(which(is.na(psi)))
nas_by_rows<-rowSums(is.na(psi)); 
nas_by_kds<-colSums(is.na(psi))
barplot(sort(nas_by_rows, decreasing = T),ylim=c(1,350), xlim=c(1,30) )
barplot(sort(nas_by_kds, decreasing = T),ylim=c(1,30), xlim=c(1,330) )
#########################################################################
#be careful because original data is trasposed respect to this data
SparseRows<-which(rowMeans(is.na(psi))*100>30)#Row numbers with >30% missing values
head(SparseRows)
#5 columns  
SparseCols<-which(colMeans(is.na(psi))*100>30)#Row numbers with >30% missing values
psiCleaned<-psi[-c(SparseRows), -c(SparseCols)]
dim(psiCleaned)#26,315
#pisa la original
plot(density(colMeans(is.na(psi))*100))
plot(density(rowMeans(is.na(psi))*100))
plot(density(colMeans(is.na(psiCleaned))*100))
plot(density(rowMeans(is.na(psiCleaned))*100))
tpsi<-t(psiCleaned)
psiI<-impute.knn(as.matrix(tpsi),
                 k=round(sqrt(nrow(tpsi))*0.25),
                 rowmax=0.5)$data
hist(rowMeans(is.na(psiI))*100)# desaparecen los NAs
#####################compute deltaPSI###############
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
dim(eventsscaled)
dim(deltascaled)
#scaled by KDs, maintain the profiles)
#Calculate Robust Correlation Matrix
C <- CRobCor(deltascaled)
pdf("byKDs_corrPlot_rnaseq_alphabet.pdf")
corrplot(C, method = "color", order = "alphabet", tl.cex = 0.1)
dev.off()
getwd()
pdf("byKds_corrPlot_rnaseq_hclust.pdf")
corrplot(C, method = "color", order = "hclust",  col = brewer.pal(n = 8, name = "RdBu"), tl.cex = 0.1)
dev.off()

Cdouble <- CRobCor(eventsscaled)
pdf("Double_corrPlot_rnaseq_alphabet.pdf")
corrplot(Cdouble, method = "color", order = "alphabet", tl.cex = 0.1)
dev.off()
pdf("Double_corrPlot_rnaseq_hclust.pdf")
corrplot(Cdouble, method = "color", order = "hclust",  col = brewer.pal(n = 8, name = "RdBu"), tl.cex = 0.1)
dev.off()
write.table(eventsscaled, "profiles/eventsscaled.tab", sep="\t")
#pisa la original
### Transpose matrix to ASEs (rows) x KDs (columns)
#columns=events
#Calculate Robust Correlation Matrix
save(C, file="Original_Events_RNASeq_C_singlescaled.RData")
write.table(C, file="Original_Events_RNASeq_C_singlescaled.tab", sep="\t", col.names = NA)
save(Cdouble, file="Original_Events_RNASeq_C_doublescaled.RData")
write.table(Cdouble, file="Original_Events_RNASeq_C_doublescaled.tab", sep="\t", col.names = NA)
#Create graph, calculate centrality measures, identify communities and plot
load("Original_Events_RNASeq_C_doublescaled.RData")
gListDouble<-CreateGraph(Cdouble,rho=0.6) 
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
save(gList, file="gList_6_RNA-Seq.RData")
getwd()
write_graph(g, format="edgelist", file="rnaseq")
##############################################################
#exporta parameters
deg.dist <- degree_distribution(g, cumulative=T, mode="all")
png("cumulativeFreq_RNA-Seq_06.png")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
  xlab="Degree", ylab="Cumulative Frequency")
dev.off()
########################################3
write.table(t(values), file="RNA_Seq_06topology.tab", sep="\t")
#add community information:
comps <-memb$membership; head(comps)#dpmde esta cada gen, em que grupo
table(comps)
write.table(comps, file="RNA_Seq_06memberships.txt", col.names=NA)
########################################


