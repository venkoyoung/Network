library(pheatmap)
setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/selectedEventsHs2/convIR//")
scaledfiles<-read.table("IR.eventscaled.csv", header = T, stringsAsFactors = F ,sep="\t")
############################################################
numEvent<- matrix(unlist(strsplit(as.character(scaledfiles$file), "_")), byrow = T , ncol=4)[,3]
Event<- matrix(unlist(strsplit(as.character(scaledfiles$file), "_")), byrow = T , ncol=4)[1,1]
############################################################
for (i in 1:nrow(scaledfiles))
{
file<-scaledfiles$file[i]  
print(file)
name<-gsub("_eventscaled.tab", "", scaledfiles$file[i])
rho<-scaledfiles$rho[i]
print(rho)
M<-read.table(file); print(dim(M))
##################################check correlation density
Cdouble <- CRobCor(M)
#################################################3
png(paste(name, "heatmap.png", sep="_"),
    width = 1000, height = 1000)
pheatmap(Cdouble,  fontsize_col= 5, fontsize_row=5,
         color = colorRampPalette(rev(brewer.pal(n =7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList,
         treeheight_row = 0,
         treeheight_col = 0)
dev.off()
dev.off()
pdf(paste(name, "heatmap.pdf", sep="_"),
    width = 20, height = 20)

pheatmap(Cdouble,  fontsize_col= 5, fontsize_row=5,
         color = colorRampPalette(rev(brewer.pal(n =7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList,
         treeheight_row = 0,
         treeheight_col = 0)
png(paste(name, "density.png", sep="_"))
plot(density(Cdouble), xlim=c(-1,1), col="red")
abline(v=0)
dev.off()

}
####################################################################################
