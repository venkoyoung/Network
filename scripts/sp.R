
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
