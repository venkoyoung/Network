setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/selectedEventsHs2/convTotalSum/tabs/")
topoFiles<-read.table("ORconvTotalSum_smallDF.txt",
                      header = F, stringsAsFactors = F ,sep="\t")
colnames(topoFiles)<-"file"
topoFiles
numEvent<-  matrix(unlist(strsplit(as.character(topoFiles$file), "_")), byrow = T , ncol=6)[,4]
Event<- matrix(unlist(strsplit(as.character(topoFiles$file), "_")), byrow = T , ncol=5)[1,2]
plotname<-paste(Event, "20-600_OR.png", sep="_")
df<-data.frame()
for (i in 1:nrow(topoFiles))  {
    print(i)
  file<-topoFiles$file[i]  
  print(file)
  name<-gsub("_df_small.tab", "", topoFiles$file[i])
  print(name)
  df<-rbind(df, read.table(file, header = F, skip=2) )
}
df
rownames(df)<-topoFiles$file
df$V1<-NULL
colnames(df)<-c("rho","NumOfEdges","NumOfVertices","MeanDegree","Diameter","Diameter2")
df$ne<-     matrix(unlist(strsplit(as.character(topoFiles$file), "_")), byrow = T , ncol=6)[,4]
df<-df[order(as.numeric(df$ne)),]
head(df)
######################
dfA3<-df; dim(dfA3); #7,7
dfA3$color<-rep("blue", nrow(dfA3))
dfA5<-df; dim(dfA5); #7,7
dfA5$color<-rep("red", nrow(dfA5))
dfES<-df; dim(dfES);# 7,7
dfES$color<-rep("green", nrow(dfES))
dfIR<-df; dim(dfIR)# 7,7
dfIR$color<-rep("orange", nrow(dfIR))
#cargo todos lo DFs y luego los comnbino y ploteo todo junto
finalDF<-rbind(dfA5, dfA3, dfIR, dfES)
head(finalDF)
######################
#Plots
setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/selectedEventsHs2/allConv/")
write.table(finalDF,"finalDF.tab", col.names = NA, sep="\t")
finalDF$color<-as.factor(finalDF$color)
ntrees <- levels(finalDF$color)
# get the range for the x and y axis 
finalDF$ne<-as.numeric(finalDF$ne)
xrange <- range(finalDF$ne) 
png("NumberOfEdges.png")
yrange <- range(as.numeric(finalDF$NumOfEdges)) 
# set up the plot 
  plot(xrange, yrange, type="n", xlab="Num of events", ylab="Num of edges" )  
  # add lines 
  for (i in 1:length(ntrees)) { 
    print(ntrees[i])
    etype <- subset(finalDF, color==ntrees[i]) 
    lines(etype$ne, etype$NumOfEdges, type="l", lwd=1.5,
   col=ntrees[i]) 
  } 

  # add a title and subtitle 
title("Number of edges")
# add a legend 
legend("topright", 
       legend=c("A3", "A5","IR","ES"),
       col=c("red", "blue","orange","green"), lty=1, cex=0.8)
dev.off()
################################################################
#repito para el resto de los parametros:
head(finalDF)
png("NumberOfVertices.png")
yrange <- range(as.numeric(finalDF$NumOfVertices)) 
# set up the plot 
plot(xrange, yrange, type="n", xlab="Num of events", ylab="Num of vertices" )  
# add lines 
for (i in 1:length(ntrees)) { 
  print(ntrees[i])
  etype <- subset(finalDF, color==ntrees[i]) 
  lines(etype$ne, etype$NumOfVertices, type="l", lwd=1.5,
        col=ntrees[i]) 
} 
# add a title and subtitle 
title("Number of vertices")
# add a legend 
legend("topright", 
       legend=c("A3", "A5","IR","ES"),
       col=c("red", "blue","orange","green"), lty=1, cex=0.8)
dev.off()
################################################################
head(finalDF)
png("MeanDegree.png")
yrange <- range(as.numeric(finalDF$MeanDegree)) 
# set up the plot 
plot(xrange, yrange, type="n", xlab="Num of events", ylab="Mean Degree" )  
# add lines 
for (i in 1:length(ntrees)) { 
  print(ntrees[i])
  etype <- subset(finalDF, color==ntrees[i]) 
  lines(etype$ne, etype$MeanDegree, type="l", lwd=1.5,
        col=ntrees[i]) 
} 
# add a title and subtitle 
title("Mean Degree")
# add a legend 
legend("topright", 
       legend=c("A3", "A5","IR","ES"),
       col=c("red", "blue","orange","green"), lty=1, cex=0.8)
dev.off()

###############################################################
head(finalDF)
  png("Diameter.png")
  yrange <- range(as.numeric(finalDF$Diameter2)) 
  # set up the plot 
  plot(xrange, yrange, type="n", xlab="Num of events", ylab="Diameter" )  
  # add lines 
  for (i in 1:length(ntrees)) { 
    print(ntrees[i])
    etype <- subset(finalDF, color==ntrees[i]) 
    lines(etype$ne, etype$Diameter2, type="l", lwd=1.5,
          col=ntrees[i]) 
  } 
  # add a title and subtitle 
  title("Diameter")
  # add a legend 
  legend("topright", 
         legend=c("A3", "A5","IR","ES"),
         col=c("red", "blue","orange","green"), lty=1, cex=0.8)
  dev.off()
#####################################
#### Total
finalDF<-df
head(finalDF)
fixrho<-finalDF
OR<-df
#####################################
finalDF$color<-rep("red", nrow(df))
# get the range for the x and y axis 
head(finalDF)
finalDF$ne<-as.numeric(finalDF$ne)
xrange <- range(finalDF$ne) 
png("NumberOfEdges_OR.png")
yrange <- range(as.numeric(finalDF$NumOfEdges)) 
# set up the plot 
plot(finalDF$ne, finalDF$NumOfEdges,
     xlim=xrange, ylim=yrange, xlab="Num of events", ylab="Num of edges",
     type="l", lwd=1.5) 
lines(OR$ne, OR$NumOfEdges, col="red")
title("Number of edges")
dev.off()
################################################################
#repito para el resto de los parametros:
png("NumberOfVertices_OR.png")
yrange <- range(as.numeric(finalDF$NumOfVertices)) 
  # set up the plot 
plot(finalDF$ne, 
     finalDF$NumOfVertices, 
     type="l", lwd=1.5,
     xlim=xrange, 
     ylim=yrange, 
     xlab="Num of events",
     ylab="Num of vertices" )  
title("Number of vertices")
lines(OR$ne, OR$NumOfVertices, col="red")
# add a legend 
dev.off()
################################################################
png("MeanDegree_OR.png")
yrange <- range(as.numeric(finalDF$MeanDegree)) 
plot(finalDF$ne, finalDF$MeanDegree, type="l", lwd=1.5,
  xlim=  xrange, 
  ylim= yrange, xlab="Num of events", ylab="Mean Degree" )  
  title("Mean Degree")
  lines(OR$ne, OR$MeanDegree, col="red")
  dev.off()
###############################################################
png("Diameter_OR.png")
yrange <- range(as.numeric(OR$Diameter2)) 
# set up the plot 
plot(finalDF$ne, finalDF$Diameter2,
       xlim= xrange, 
     ylim=c(1,30), 
       xlab="Num of events", ylab="Diameter" , type="l", lwd=1.5)
title("Diameter")
lines(OR$ne, OR$Diameter2, col="red")
  # add a legend 
  dev.off()