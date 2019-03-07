library(pheatmap)
library(tidyr)
library(data.table)
library(ggplot2)
library(limma)
library(data.table)
library(stringr)
library(ggplot2)
library(reshape) 
##################################################################
setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/selectedEventsHs2/randomA3_30_15/")
edgelists<-read.table("edgelists.tab", header = F, stringsAsFactors = F ,sep="\t")
edgelists[1:5,]
finalDF<-data.frame()
linksList<-list()
length(linksList)
nodesList<-list()
for (i in 1:nrow(edgelists))
{

    file<-edgelists$V1[i]
    print(file)
    name<-gsub("_edgelist.tab", "", edgelists$V1[i])
    print(name)
    edgeListLC<-read.table(file);
    print(head(edgeListLC));
    print(dim(edgeListLC))
    ###########################################################
    #Links
    l11<-paste(edgeListLC$source, edgeListLC$target, sep="-")
    linksList[[i]]<-l11
    names(linksList)[[i]]<-name
    
    ###########################################################
    #nodes
    n1<-unique( c(as.character(edgeListLC$source),
              as.character(edgeListLC$target)))
    nodesList[[i]]<-n1
    names(nodesList)[[i]]<-name
    
  }
##################################################################
allNodes<-unique(unlist(nodesList)); length(allNodes)
allNodesDF<-data.frame(matrix(0, ncol=length(nodesList), nrow=length(allNodes)), row.names = allNodes)
head(allNodesDF)
#select KD:
for (i in 1:length(nodesList))
{
  ii<-match(nodesList[[i]], rownames(allNodesDF))
  colnames(allNodesDF)[i]<-names(nodesList)[i]
  allNodesDF[ii,i]<-1  
}

##################################################################
#finalMatrixCOmparison
dim(allNodesDF)
KDsSum<-rowSums(allNodesDF)
KDsSum
sort(KDsSum)[1:10]
#network de 100 nodos:
KDs<-sort(KDsSum, decreasing = T)[1:100]#MIN 81
summary(KDs)
KDs[1:10]
KDsSum[1:5]
topNodes<-allNodesDF[names(KDs),]
dim(topNodes)
#esto se puede plotear como heatmap
pdf("presence_ausence_ES_150_30_100_min98.pdf",width = 16,  height = 20)
pheatmap(topNodes,
         fontsize = 6, 
         fontsize_row = 10, 
         fontsize_col = 6,
         color=c("red","white"),
         main="presence_ausence_A3_150_30_100_min98")
dev.off()
#################################################################
allLinks<-unique(unlist(linksList)); length(allLinks)  
allLinksDF<-data.frame(matrix(0, 
                              ncol=length(linksList), 
                              nrow=length(allLinks)), 
                       row.names = allLinks)
head(allLinksDF)
dim(allLinksDF)#34528 *100
for (i in 1:length(linksList))
{
  ii<-match(linksList[[i]], rownames(allLinksDF))
  colnames(allLinksDF)[i]<-names(linksList)[i]
  allLinksDF[ii,i]<-1  
}
head(allLinksDF)
dim(allLinksDF)
##################################################################
head(allLinksDF)
#filtered<-allLinksDF[rowSums(allLinksDF)>20,]
filtered<-allLinksDF
dim(filtered)#343*100
resultingEdgeList<-data.frame(matrix(unlist(strsplit(rownames(filtered),"-")), ncol=2, byrow = T))
head(resultingEdgeList)
rownames(resultingEdgeList)<-rownames(filtered)
resultingEdgeList$occur<-rowSums(filtered)
colnames(resultingEdgeList)<-c("from","to","occur")
dim(resultingEdgeList)
write.table(resultingEdgeList,"30_150_100_A3_oc.txt",sep="\t", col.names = NA, quote = F)
##################################################################
g<-simplify(graph_from_data_frame(resultingEdgeList[,1:2], directed = FALSE),
              remove.multiple = TRUE, remove.loops = TRUE)
df2<-as.data.frame(get.edgelist(g))
head(df2)
dim(df2)
g
df2$consensus<-rep("30_600_all", 41251)
head(df2)
colnames(df2)<-c("source", "target","dset")
write.table(df2,"30_600_all_rho05.txt",sep="\t", col.names = NA, quote = F)
