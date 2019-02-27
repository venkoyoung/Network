library(data.table)
library(ggplot2)
all<-read.table("~/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/lastVtHsa/VT_Hg19_NA3_NAI.tab",
                header = T, sep="\t", row.names = 1);
dim(all)
#remove bad samples
data <- all[, ! colnames(all) %in% 
                  c("AA2",
                    "AA1",
                    "CCDC12",
                    "C1orf55",
                    "C1orf55_b",
                    "CDC5L",
                    "HFM1",
                    "LENG1",
                    "RBM17",
                    "PPIL1",
                    "SRRT"),
                  drop = F]
#remove 11 samples
dim(data)
#select the events:
orig_events<-read.table("~/Dropbox (CRG)/Personal_Estefania/GosiaAndEstefi/AS_317_OK/Original_Network_events.tsv", header=T, sep="\t")
head(orig_events)
ii<-match(orig_events$EVENT, rownames(data)); is.na(ii)
dfOrigEvents<-data[ii,]; 
head(dfOrigEvents);dim(dfOrigEvents)
eventsIDs<-paste(rownames(data)[ii],1:length(rownames(data)[ii]), sep="_")
rownames(dfOrigEvents)<-paste(data$GENE[ii], eventsIDs, sep=":")
onlyValues<-dfOrigEvents[,7:ncol(dfOrigEvents)]
psi<-onlyValues; head(onlyValues); dim(onlyValues)
#match the names of the KDs:
####################
psi[psi=="NAold"]<-NA
psi[psi=="NAnew3"]<-NA
head(psi)
length(which(is.na(psi)))#1405
dim(psi)
#desde aca vengo con el script renames
####################################################################
nas_by_rows<-rowSums(is.na(psi)); 
nas_by_kds<-colSums(is.na(psi))
barplot(sort(nas_by_rows, decreasing = T),ylim=c(1,350), xlim=c(1,30) )
barplot(sort(nas_by_kds, decreasing = T),ylim=c(1,30), xlim=c(1,330) )
#########################################################################
#be careful because original data is trasposed respect to this data
SparseRows<-which(rowMeans(is.na(psi))*100>30)#Row numbers with >30% missing values
length(SparseRows)#5 columns  
SparseCols<-which(colMeans(is.na(psi))*100>30)#Row numbers with >30% missing values
length(SparseCols)
psiCleaned<-psi[-c(SparseRows), -c(SparseCols)]
dim(psiCleaned)#26,315
#pisa la original
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
########################################
write.table(eventsscaled,"eventscaled.tab", sep="\t")
Cdouble <- CRobCor(eventsscaled)
pdf("Double_corrPlot_rnaseq_hclust.pdf")
corrplot(Cdouble, method = "color", order = "hclust",  col = brewer.pal(n = 8, name = "RdBu"), tl.cex = 0.1)
dev.off()
#pisa la original
### Transpose matrix to ASEs (rows) x KDs (columns)
#columns=events
#Calculate Robust Correlation Matrix
save(Cdouble, file="Original_Events_RNASeq_C_doublescaled.RData")
write.table(Cdouble, file="Original_Events_RNASeq_C_doublescaled.tab", sep="\t", col.names = NA)
#Create graph, calculate centrality measures, identify communities and plot
gListDouble<-CreateGraph(Cdouble,rho=0.6) 
gList<-gListDouble
# Rho was set to achieve FDR < 5%
##############################################
g<-gList[[1]]
memb<-gList[[2]]
values<-gList[[3]]
mean(values[1,])#1.71
V(g)
E(g)
deg <- degree(g, mode="all")
mean(deg)#2.9
hist(values[1,])#
summary(values[1,])
save(gList, file="gList_6_RNA-Seq.RData")
write_graph(g, format="edgelist", file="rnaseq")
df<-get.edgelist(g); 
colnames(df)<-c("source", "target")
write.table(df, file="edgelist_06_rnaseq.tab", col.names = NA, quote = F, sep="\t")

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




