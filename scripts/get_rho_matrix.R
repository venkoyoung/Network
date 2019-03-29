setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/selectedEventsHs2/results")
rhofiles<-read.table("NumberOfEdgesAlt3.tab", header = T, stringsAsFactors = F ,sep="\t")
df1<-read.table(rhofiles$file[1], header = T, row.names = 1)
name<-gsub("\\.tab", "",rhofiles$file[1])
Oname<-gsub("NumberOfEdges_", "",name)
######################################################
dfRho<-data.frame(df1$ratio)
rownames(dfRho)<-df1$rho
dfTE<-data.frame(df1$true)
rownames(dfTE)<-df1$rho
NumOfEvents<- unlist(strsplit(Oname, "_"))[4]
colnames(dfRho)[1]<-as.numeric(NumOfEvents)
colnames(dfTE)[1]<-as.numeric(NumOfEvents)

for (i in 2:nrow(rhofiles))
{
name<-gsub("\\.tab", "",rhofiles$file[i])
Oname<-gsub("NumberOfEdges_", "",name)
print(i)
print(Oname)
NumOfEvents<- unlist(strsplit(Oname, "_"))[4]
print(NumOfEvents)
rhocurve<-read.table(rhofiles$file[i], header = T, row.names = 1)
#######################################################
dfRho<-cbind(dfRho, rhocurve$ratio)
dfTE<-cbind(dfTE,rhocurve$true)
colnames(dfRho)[i]<-as.numeric(NumOfEvents)
colnames(dfTE)[i]<-as.numeric(NumOfEvents)
}
#######################################################
dfRho<-dfRho*100
dfRho<-dfRho[,order(as.numeric(colnames(dfRho)))]
dfTE<-dfTE[,order(as.numeric(colnames(dfTE)))]	
write.table(dfRho, "optimalRho_Matrix_alt3.tab", col.names=NA)
write.table(dfTE, "optimalTE_Matrix_alt3.tab",col.names=NA)
