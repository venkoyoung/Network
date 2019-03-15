
#ppi validation
library(igraph)
setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/consenso/")
##################################################################
#consensus network
t1<-read.table("30_600_100_rho05.txt");head(t1); dim(t1)
g1<-graph_from_data_frame(t1, directed = F)
length(E(g1))#97
length(V(g1))#48

sp<-distance_table(g1, directed = F)
sp<-shortest_paths(g1,3)
sp$vpath
 600/20
 barplot(seq(20,600, 20), cex.names = 1.5, besides=T)
         cex.names=seq(20,600, 20))
 ?barplot
 ?seq
 barplot(1:30)
 df<-data.frame(dose=seq(20,600, 20), len=seq(20,600, 20))
 library(ggplot2)
 dev.off()
 p<-
 svg("subsampling600.svg")  
   ggplot(data=df, aes(x=dose, y=len)) +
   geom_bar(stat="identity")+
   theme_minimal()+
   geom_text(aes(label=len), vjust=2, color="white", size=3)
   dev.off()
   getwd()
 
   
   140/20
 p
 