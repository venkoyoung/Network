library(ggplot2)
setwd("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/rangesDSets/MIC/tabs/")
  topoFiles<-read.table("smallDF_all.txt",
                      header = F, stringsAsFactors = F ,sep="\t")
colnames(topoFiles)<-"file"
topoFiles
numEvent<-  matrix(unlist(strsplit(as.character(topoFiles$file), "_")), byrow = T , ncol=7)[,4]
Event<- matrix(unlist(strsplit(as.character(topoFiles$file), "_")), byrow = T , ncol=7)[,1]
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
head(df)
dim(df)
rownames(df)<-topoFiles$file
df$V1<-NULL
colnames(df)<-c("rho","NumOfEdges","NumOfVertices","MeanDegree","Diameter","Diameter2")
df$nrep<- matrix(unlist(strsplit(as.character(topoFiles$file), "_")), byrow = T , ncol=7)[,5]
df$ne<-   as.numeric(matrix(unlist(strsplit(as.character(topoFiles$file), "_")), byrow = T , ncol=7)[,4])
df$event<- matrix(unlist(strsplit(as.character(topoFiles$file), "_")), byrow = T , ncol=7)[,1]
df$range<- matrix(unlist(strsplit(as.character(topoFiles$file), "_")), byrow = T , ncol=7)[,3]
table(df$range)
#################################################
dfMIC<-df[df$event=="MIC",]
dfMIC0<-df[df$range=="range0",]
dim(dfMIC0)
dfMIC5<-df[df$range=="range5",]
dim(dfMIC5)
dfMIC15<-df[df$range=="range15",]
dim(dfMIC15)
dfMIC0$ne
dfMIC5$ne
dfMIC15$ne

plot(dfMIC5$ne, dfMIC5$NumOfEdges)
plot(dfMIC15$ne, dfMIC15$NumOfEdges)
lines(meanMIC15$Group.1,meanMIC15$x)
meanMIC0<-aggregate(dfMIC0$NumOfEdges, by=list(dfMIC0$ne), FUN=mean)
meanMIC5<-aggregate(dfMIC5$NumOfEdges, by=list(dfMIC5$ne), FUN=mean)
meanMIC15<-aggregate(dfMIC15$NumOfEdges, by=list(dfMIC15$ne), FUN=mean)


dir.create("plots")
pdf("../plots/MIC_all_Events_NumofEdges.pdf", width = 11.69,  height = 8.27)
  ggplot(df) + 
  geom_point(aes(x = ne, y = NumOfEdges, color=range), size=0.1)+
  theme_classic()   +
  geom_line(data=meanMIC0, aes(x=Group.1, y=x))+
    geom_line(data=meanMIC5, aes(x=Group.1, y=x))+    
    geom_line(data=meanMIC15, aes(x=Group.1, y=x))+    
  scale_x_continuous(breaks = seq(0, 500, by = 20), name="Number of events",limits = c(0,600))+
  #scale_y_continuous(breaks = seq(0, 50000, by = 1000), name="Number of Edges",limits = c(0,10000))+
  scale_y_continuous( name="Number of Edges")+
  geom_hline(yintercept=c(0, 1000), color='black', size=0.3)+
  geom_vline(xintercept=0, color='black', size=0.3)+
  theme(
  axis.text.x = element_text(angle = -90, hjust = 0.5, size = 5))  
  dev.off()
####################################################################
mean_V_MIC0<-aggregate(dfMIC0$NumOfVertices, by=list(dfMIC0$ne), FUN=mean)
mean_V_MIC5<-aggregate(dfMIC5$NumOfVertices, by=list(dfMIC5$ne), FUN=mean)
mean_V_MIC15<-aggregate(dfMIC15$NumOfVertices, by=list(dfMIC15$ne), FUN=mean)
  
pdf("plots/MIC_all_NumberOfVertices.pdf", width = 11.69,  height = 8.27)
ggplot(df) + 
  geom_point(aes(x = ne, y = NumOfVertices, color=range), size=0.1)+
  theme_classic()   +
  geom_line(data=mean_V_MIC0, aes(x=Group.1, y=x))+
  geom_line(data=mean_V_MIC5, aes(x=Group.1, y=x))+
  geom_line(data=mean_V_MIC15, aes(x=Group.1, y=x))+
  scale_x_continuous(breaks = seq(0, 500, by = 20), name="Number of events", limits = c(0,600))+
  scale_y_continuous(breaks = seq(0, 300, by = 20), name="Number of Vertices (KDs)",
                     limits = c(0,300))+
  theme(  axis.text.x = element_text(angle = -90, hjust = 0.5, size = 5)) +
  geom_hline(yintercept=0, color='black', size=0.3)+
  geom_vline(xintercept=0, color='black', size=0.3)
dev.off()

#################################################################
pdf("plots/Mix_NumofEdges.pdf")
ggplot(dfMix) + geom_point(aes(x = ne, y = NumOfEdges), size=0.1)+
theme_classic()   +
geom_line(data=meanTotalNE, aes(x=Group.1, y=x))+
scale_x_continuous(breaks = seq(0, 600, by = 20), name="Number of events")+
scale_y_continuous( name="Number of Edges")+
theme(
  axis.text.x = element_text(angle = -90, hjust = 0.5, size = 8)) 
dev.off()
#################################################################
pdf("plots/Mix_NumofVertices.pdf")
ggplot(dfMix) + geom_point(aes(x = ne, y = NumOfVertices), size=0.1)+
  theme_classic()   +
  geom_line(data=meanTotalNV, aes(x=Group.1, y=x))+
  scale_x_continuous(breaks = seq(0, 600, by = 20), name="Number of events")+
  scale_y_continuous(breaks = seq(0, 300, by = 20), name="Number of Vertices (KDs)")+
  theme(
  axis.text.x = element_text(angle = -90, hjust = 0.5, size = 8)) 
dev.off()
################################################################
pdf("plots/all_NumofVertices.pdf")
ggplot(df) + geom_point(aes(x = ne, y = NumOfVertices,  color=event), size=0.3)+
  theme_classic()   +
  geom_line(data=meanTotalNV, aes(x=Group.1, y=x), color="purple")+
  geom_line(data=meanExNV, aes(x=Group.1, y=x), color="green")+
  geom_line(data=meanInNV, aes(x=Group.1, y=x), color="blue")+
  geom_line(data=meanA5NV, aes(x=Group.1, y=x), color="yellow")+
  geom_line(data=meanA3NV, aes(x=Group.1, y=x), color="red")+
  scale_x_continuous(breaks = seq(0, 600, by = 20), name="Number of events")+
  scale_y_continuous(breaks = seq(0, 300, by = 20), name="Number of Vertices (KDs)")+
  theme( axis.text.x = element_text(angle = -90, hjust = 0.5, size = 8)) 
dev.off()
##################################################################
pdf("plots/all_NumofEdges_lines.pdf")
ggplot(df) +
# geom_point(aes(x = ne, y = NumOfEdges,  color=event), size=0.3)+
  theme_classic()   +
  geom_line(data=meanTotalNE, aes(x=Group.1, y=x), color="purple")+
  geom_line(data=meanExNE, aes(x=Group.1, y=x), color="green")+
  geom_line(data=meanInNE, aes(x=Group.1, y=x), color="blue")+
  geom_line(data=meanA5NE, aes(x=Group.1, y=x), color="yellow")+
  geom_line(data=meanA3NE, aes(x=Group.1, y=x), color="red")+
  scale_x_continuous(breaks = seq(0, 600, by = 20), name="Number of events")+
  #scale_y_continuous(name="Number of Edges")+
  scale_y_continuous( name="Number of Edges", limits = c(0,7000))+
  theme( axis.text.x = element_text(angle = -90, hjust = 0.5, size = 8)) 
dev.off()


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