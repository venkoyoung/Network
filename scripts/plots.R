install.packages("pheatmap")
library(pheatmap)
setwd("Documents/surfaces/")
data<-read.csv("optimalRho.csv", sep="\t", header =T, row.names = 1)
data[data>100]<-100
colnames(data)<-  gsub("X", "", colnames(data))
tp<-read.csv("optimalTP.csv", sep="\t", header =T, row.names = 1)
colnames(tp)<-  gsub("X", "", colnames(tp))
#######################################
cols <- heat.colors(30)
par(mfrow = c(2, 1))
pdf("rho_curves_mix.pdf", width =12, height = 8)
plot(NA, ylim = c(0, 100), xlim = c(1, 15))
for(i in 1:30) lines(data[[i]], col = cols[i], lwd=2)
legend("bottomright", legend=(colnames(data)), lwd=2, col=cols)
abline(h=c(5,10))
#######################################
plot(NA, ylim = c(0, 6000), xlim = c(1, 15))
for(i in 1:30) lines(tp[[i]], col = cols[i], lwd=2)
legend("bottomright", legend=(colnames(tp)), lwd=2, col=cols)
dev.off()
#######################################
#division
dataCorr<-data
dataCorr[dataCorr>5]<-100
ratio<-tp/dataCorr
n=20
cc<-rev(heat.colors(n, alpha = 1))
pdf("heatmaps_ratios_mix.pdf", width =8, height = 4)
pheatmap(data,  cluster_rows = FALSE,  cluster_cols = FALSE, color=cc )
pheatmap(tp,    cluster_rows = FALSE,  cluster_cols = FALSE, color=cc )
pheatmap(dataCorr,  cluster_rows = FALSE,  cluster_cols = FALSE, color=cc )
pheatmap(ratio,    cluster_rows = FALSE,  cluster_cols = FALSE, color=cc )
dev.off()
