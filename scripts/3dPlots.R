library(plotly)
setwd("Dropbox (CRG)/Personal_Estefania/Network/3dplots/")
datos<-read.table("lrtTable_331_high.txt",  header=T)
plot(datos$logCPM,
     datos[,2], xlim=c(0,15), ylim=c(-5,5) , type="p", pch=20)
head(datos); dim(datos)
colnames(datos)  <-sub("logFC.group","",colnames(datos))
#
sortedAll<-sort(apply(datos[,-c(324:327)], 2, sd))
names(sortedAll)[1:10]
length(sortedAll)#322
orderDatos<-datos[,names(sortedAll)]; dim(orderDatos)
orderDatos$logCPM<-datos$logCPM
#################################################
subset<-data.frame(orderDatos[,1:5], logCPM=datos$logCPM)
head(subset)
colnames(subset); dim(subset)
logfc<-unlist(subset[,1:5]); length(logfc)
cpm<-rep(subset$logCPM, 5); length(cpm)
kdsNames<-colnames(subset); length(kdsNames)
kd<-rep(kdsNames[1:5], each=nrow(subset)) ; length(kd)
kdNumber<-rep(1:5, each=nrow(subset)) ; length(kdNumber)
dataK<-data.frame(  x = cpm,  y = logfc,  cut = kdNumber);head(dataK)
####################################################
#full:
dim(orderDatos)
logfcA<-unlist(orderDatos[,1:100]); length(logfcA)
cpmA<-rep(subset$logCPM, 100); length(cpmA)
kdsNamesA<-colnames(orderDatos[1:100]); length(kdsNamesA)
kdA<-rep(kdsNamesA, each=nrow(orderDatos)) ; length(kdA)
kdNumberA<-rep(1:100, each=nrow(orderDatos)) ; length(kdNumberA)
dataA<-data.frame(  x = cpmA,  
                    y = logfcA,  cut = kdNumberA);head(dataA)

dim(dataA)
