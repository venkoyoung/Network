##### define layouts and coord for network 1 (consensus) #####
set.seed(1)
g1
layg1 <- layout.fruchterman.reingold(g1)
rownames(layg1)<-V(g1)$name
#define coords for network 2
g2
set.seed(2)
layg2 <- layout.fruchterman.reingold(g2)
rownames(layg2)<-V(g2)$name
# overwrite coords for shared nodes
#############################################################
ii<-match(rownames(layg2), row.names(layg1))
ii<-ii[!is.na(ii)]
layg2<-layg1[ii,]
xlim <- range(c(layg1[,1], layg2[,1]))
ylim <- range(c(layg1[,2], layg2[,2]))
V(g2)$color <- "red"
E(g2)$color <- "red"
pdf("both_overlap_Cons.pdf", width =11,  height = 8)
plot(g1, 
     layout=layg1, 
     xlim=xlim, 
     ylim=ylim, 
     rescale=FALSE,
     vertex.label.dist=3,
     vertex.label.cex=0.8,
     vertex.label.font=1,
     vertex.size=20)
plot(g2 , 
     vertex.size=20, 
     layout=layg2, 
     xlim=xlim, 
     ylim=ylim, 
     rescale=FALSE,
     vertex.label=NA,
     add = TRUE)
dev.off()
#############################################################
pdf("both_par_Consensus.pdf", width = 11,  height =8)
par(mfrow=c(1,2))
plot(g1, 
     layout=layg1, 
     xlim=xlim, 
     ylim=ylim, 
     rescale=FALSE,
     vertex.label.dist=5,
     vertex.label.cex=0.8,
     vertex.label.font=20,
     vertex.size=20)
dev.off()
plot(g2 , 
     vertex.size=20, 
     layout=layg2, 
     xlim=xlim, 
     ylim=ylim, 
     rescale=FALSE,
     vertex.label.dist=5,
     vertex.label.cex=0.8,
     vertex.label.font=20,
     vertex.size=50)
dev.off()
#############################################################
