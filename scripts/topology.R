#first generate data: TRUE and RANDOM DATA
library(QUIC)
library(igraph)
library("RColorBrewer")
library(scales)
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/CRobCor.R")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/functions_fdr_noMC.R")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/CreateGraph.R")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/CentralityRanking.R")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/Vscale.R")
library(igraph)
############################################################
getMeanDeg<-function(start, end, gap, Cdouble, name)
{
  
  df<-data.frame(rho=seq(start,end,gap))
  for (i in 1:length(seq(start,end,gap))) {
    print(i)
    gListDouble<-CreateGraph(Cdouble,rho=df$rho[i]) #
    gList<-gListDouble
    g<-gList[[1]]
    deg <- degree(g, mode="all")
    df$md[i]<-mean(deg)
    print(mean(deg))
    file=paste(paste("MeanDegree", name, sep="_"), "png", sep=".")
    png(file)
    plot(df$rho, df$md, xlim=c(0,1), ylim=c(0,max(df$md)), type = "l", main=paste("Mean degree", name))
    abline(h=6)
    dev.off()
  
}
  
return(df)
}

###############################################################################
NetworkDesc<-function(g, name, rho)
  {
df<-data.frame()
dfSmall<-data.frame()
df[1,1]<-"NumOfEdges"
df[1,2]<-length(E(g))#547
print(paste("NumOfEdges",length(E(g)), sep=": " ))
df[2,1]<-"NumOfVertices"
df[2,2]<-length(V(g))
print(paste("NumOfVertices",length(V(g)), sep=": " ))
##################################################################
deg <- degree(g, mode="all")
deg.dist <- degree_distribution(g, cumulative=T, mode="all")
df[3,1]<-"MeanDegree"
df[3,2]<-mean(deg)#3.97
print(paste("MeanDegree",round(mean(deg),2), sep=": " ))
file1<-paste(name,"cumulativeFreq.png" , sep="_")
print(file1)
########################################################
png(file1)
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency", main=name)
dev.off()
##################################################################
file1df<-paste(name,"cumulativeFreq.tab" , sep="_")
print(file1df)
df1<-data.frame( deg=0:max(deg), cum=1-deg.dist)
write.table(df1, file=file1df, sep="\t", col.names = NA)
########################################################
file2<-paste(name, "HistDegree.png", sep="_")
print(file2)
png(file2)
hist(deg, main=name)
dev.off()
#densidad:
df[4,1]<-"EdgeDensity"
df[4,2]<-edge_density(g, loops = FALSE)
print(paste("EdgeDensity",round(edge_density(g, loops = FALSE),2), sep=": " ))

shortest.paths(g)
diameter(g,directed=FALSE,unconnected=FALSE,weights=NULL)
df[4,1]<-"Diameter"
df[4,2]<- diameter(g)
print(paste("Diameter", round(diameter(g),2), sep=": " ))


df[5,1]<-"Diameter 2"
df[5,2]<-length(get.diameter(g))
print(paste("Diameter 2", length(get.diameter(g)), sep=": " ))

#transitivity
df[6,1]<-"Ave Transitivity"
df[6,2]<- round(transitivity(g, type="average"),2)
print(paste("Ave Transitivity", round(transitivity(g, type="average"),2), sep=": " ))

df[7,1]<-"centr_degree"
df[7,2]<- round(centr_degree(g)$centralization,2)
print(paste("centr_degree", round(centr_degree(g)$centralization,2), sep=": " ))

df[8,1]<-"centr_clo"
df[8,2]<- round(centr_clo(g)$centralization,2)
print(paste("centr_clo", round(centr_clo(g)$centralization,2), sep=": " ))

df[8,1]<-"centr_betw"
df[8,2]<- round(centr_betw(g)$centralization,2)
print(paste("centr_betw", round(centr_betw(g)$centralization,2), sep=": " ))

df[9,1]<-"centr_eigen"
df[9,2]<- round(centr_eigen(g, directed = FALSE)$centralization, 2)
print(paste("centr_eigen", round(centr_eigen(g)$centralization,2), sep=": " ))

df[10,1]<-"Motifs3"
df[10,2]<- count_motifs(g, 3)
print(paste("Motifs 3",  count_motifs(g, 3), sep=": " ))

df[11,1]<-"Motifs4"
df[11,2]<- count_motifs(g, 4)
print(paste("Motifs 4",  count_motifs(g, 4), sep=": " ))

file3<-paste(name, "plot.pdf", sep="_")
pdf(file3)
print(file3)
plot.igraph(g,
            layout=layout.fruchterman.reingold,
            vertex.size=1,
            vertex.label.cex=.5,
            edge.arrow.size=.5)
dev.off()
#####################################################
file4<-paste(name, "df.tab", sep="_")
write.table(df, file4, sep="\t")

file4<-paste(name, "df_small.tab", sep="_")
dfSmall<-cbind(rho=rho, t(df[1:5,]))
head(dfSmall)
write.table(dfSmall, file4, sep="\t")

###################################################
df2<-get.edgelist(g)
colnames(df2)<-c("source", "target")
file5<-paste(name, "edgelist.tab", sep="_")
write.table(df2, file5, col.names = NA, quote = F, sep="\t")

DG<-degree(g)
NDG<-DG / df[1,2]
PR<-page.rank(g)$vector  		#Pagerank Score
OPR<-order(PR,decreasing=TRUE)
BC<-betweenness(g)		#Betweeness Centrality
OBC<-order(BC,decreasing=TRUE)	

CC<-closeness(g)		#Closeness Centrality
OCC<-order(CC,decreasing=TRUE)	

  RANK<-rank(PR)+rank(BC)+rank(CC)
ORANK<-order(RANK,decreasing=TRUE)	#Ordering of Nodes based on aggregate score.

file6<-paste(name, "CentralityRanking.tab", sep="_")

df3<-data.frame(DG,NDG,PR,OPR,BC,OBC,CC,OCC,RANK,ORANK)
head(df3)
write.table(df3, file6, col.names = NA, quote = F, sep="\t")
print(getwd())
##############################################################
#normalized degree
}
#####################################################################
getwd()
setwd("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/rangesDSets/MIC/")
scaledfiles<-read.table("eventscaled.tab", header = F, stringsAsFactors = F ,sep="\t")
scaledfiles
colnames(scaledfiles)<-"file"
scaledfiles$rho<-rep(0.5, nrow(scaledfiles))
  for (i in 1:nrow(scaledfiles))
      {
        file<-scaledfiles$file[i]  
        print(file)
        name<-gsub("_eventscaled.tab", "", scaledfiles$file[i])
    #    name<-paste("OR", name, sep="_")
        print(name)
        rho<-scaledfiles$rho[i]
        print(rho)
        M<-read.table(file); print(dim(M))
        
        Cdouble <- CRobCor(M)
        
        png(paste(name, "density.png", sep="_"))
        plot(density(Cdouble), xlim=c(-1,1), col="red")
        abline(v=0)
        dev.off()
        
        gListDouble<-CreateGraph(Cdouble,rho) 
        gList<-gListDouble
        g<-gList[[1]]
        NetworkDesc(g, name, rho)
          
      }
    
