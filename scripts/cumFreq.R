setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/selectedEventsHs2/random/exons/")
freqFiles<-read.table("ES_cumulative_files.txt",
                      header = F, stringsAsFactors = F ,sep="\t")
colnames(freqFiles)<-"file"
df<-data.frame()
for (i in 1:nrow(freqFiles))
{
  file<-freqFiles$file[i]  
  print(file)
  name<-gsub("_cumulativeFreq.tab", "", scaledfiles$file[i])
  print(name)
  df<-rbind(df,read.table(file, header = T, row.names = 1))
}
###############################################################
dfAg<-aggregate(df$cum, by=list(df$deg), FUN=mean)
colnames(dfAg)<-c("deg","cum")
###############################################################
plotname<-unlist(strsplit(name, "_"))[1]
fileName<-paste(plotname, "CumFreqAll.pdf", sep="_")
pdf(fileName)
plot(df$deg,df$cum, 
     xlim = c(0,max(df$deg)),
     ylim = c(0,1), 
     pch=20, cex=0.2,
     xlab="Degree", ylab="Cumulative Freq",
     main=paste("Cumulative Freq",plotname))
lines(dfAg$deg,dfAg$cum, col="red",cex=2)
dev.off()
###############################################################
meanA3<-dfAg; dim(meanA3)
meanA5<-dfAg; dim(meanA5)
meanIR<-dfAg; dim(meanIR)
meanES<-dfAg; dim(meanES)


meanA3$col<-rep("red", nrow(meanA3))
meanA5$col<-rep("blue", nrow(meanA5))
meanIR$col<-rep("yellow", nrow(meanIR))
meanES$col<-rep("green", nrow(meanES))

finalDF<-rbind(meanA3, meanA5, meanIR, meanES)

pdf("CumulativeFreq.pdf")
plot(finalDF$deg, finalDF$cum, col=finalDF$col, cex=1, pch=20, 
     xlab="Degree", ylab="Cumulative Freq",
     main=paste("Cumulative Freq"))
legend(60, 0.4, 
       legend=c("A3", "A5","IR","ES"),
       col=c("red", "blue","yellow","green"), lty=1, cex=0.8)
dev.off()
