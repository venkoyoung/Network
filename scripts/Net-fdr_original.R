#first generate data: TRUE and RANDOM DATA
library(QUIC)
library(igraph)
library(corrplot)
library("RColorBrewer")
library(scales)

#parameters
t_Zvals<-read.table("~/Dropbox (CRG ADV)/Personal_Estefania/Network/data/labchip_cleaned_idata.tab")
sampleData <- as.matrix(t_Zvals)
setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/data/")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/CRobCor.R")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/functions_fdr_noMC.R")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/CreateGraph.R")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/CentralityRanking.R")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/Vscale.R")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/")

#FUNCTION:
estimateFDR(
  inputM=sampleData,
  NumRandomM=10,
  start=0.3,
  end=0.6,
  interval=0.05)
######################################  
#build the network
CLabChips <- CRobCor(sampleData)
getwd()
pdf("corrPlot_labchips_alphabet.pdf")
  corrplot(CLabChips, method = "color", order = "alphabet", tl.cex = 0.1)
dev.off()


pdf("corrPlot_labchips_hclust.pdf")
  corrplot(CLabChips, method = "color", order = "hclust",  col = brewer.pal(n = 8, name = "RdBu"), tl.cex = 0.1)
dev.off()
save(CLabChips, file="C_labchip.RData")
class(CLabChips)
write.table(CLabChips, file="C_labchip.tab", sep="\t", col.names = NA)
####plot###############
gList<-CreateGraph(CLabChips,rho=0.48)   # Rho was set to achieve FDR < 5%
g<-gList[[1]]
memb<-gList[[2]]
values<-gList[[3]]
mean(values[1,])#6.1  
deg <- degree(g, mode="all")
deg.dist <- degree_distribution(g, cumulative=T, mode="all")
V(g)
E(g)
mean(deg)
hist(values[1,])#6.1  
head(values)
summary(values[1,])
##############################################################
#exporta parameters
setwd("../NumOfM100/")
png("cumulativeFreqLabChip_048.png")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")
  dev.off()
########################################3
write.table(t(values), file="048topology.tab", sep="\t")
comps <-memb$membership; head(comps)#dpmde esta cada gen, em que grupo
write.table(comps, file="048memberships.txt", col.names=NA)
########################################
#save grapho, as nodes:
write_graph(g, format="edgelist", file="edgelist_048.tab")
names(V(g))
df<-get.edgelist(g); 
colnames(df)<-c("source", "target")
write.table(df, file="edgelist_048.tab", col.names = NA, quote = F, sep="\t")
############################################################################

