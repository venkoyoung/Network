setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/selectedEventsHs2/convTotalSum/")
rhofiles<-read.table("NumberOfEdges.tab", header = T, stringsAsFactors = F ,sep="\t")
head(rhofiles)
df<-data.frame() 
dim(df)
i<-0
for (i in 1:nrow(rhofiles))
{
name<-gsub("\\.tab", "",rhofiles$file[i])
Oname<-gsub("NumberOfEdges_", "",name)
print(i)
print(Oname)
rhocurve<-read.table(rhofiles$file[i], header = T, row.names = 1)
plotfile<-paste(Oname, "ratio.png", sep="_")
if ( length(which(rhocurve$ratio*100<5))==0  ) {
  print(Oname)
  MinRho<-rhocurve$rho[which.min(rhocurve$ratio)]   }
else {   
  MinRho<-min(rhocurve$rho[rhocurve$ratio*100<=5]) 
  print("OK")
}
#######################################################
print(paste(Oname, MinRho, sep=" : "))  
df<-rbind(df, data.frame(Oname, MinRho ))  

png(plotfile)
plot(rhocurve$rho, rhocurve$ratio*100, type = "l", col="blue", ylim=c(0,100), xlab="rho", ylab="ratio", main=Oname)
abline(h=5, v=MinRho)
legend("topleft", paste("Min rho =",MinRho )  )
dev.off()
}
#######################################################
#######################################################
df$Oname
df$ne<- matrix(unlist(strsplit(as.character(df$Oname), "_")), byrow = T , ncol=4)[,4]
df$event<-matrix(unlist(strsplit(as.character(df$Oname), "_")), byrow = T , ncol=4)[,2]
table(df$event)
write.table(df, file="rho_vs_numOfEvents.txt", sep="\t", col.names = NA)
#######################################################
rhoA3<-df[df$event=="Alt3",]
rhoA3<-rhoA3[order(as.numeric(rhoA3$ne)),]
png("A3_NumOfEventsVsRho.png")
plot(rhoA3$ne,rhoA3$MinRho, type = "l",  xlab="Num of events", ylab="rho", main=unique(rhoA3$event))
dev.off()        
#######################################################
rhoA5<-df[df$event=="Alt5",]
rhoA5<-rhoA5[order(as.numeric(rhoA5$ne)),]
png("A5_NumOfEventsVsRho.png")
plot(rhoA5$ne,rhoA5$MinRho, type = "l",  xlab="Num of events", ylab="rho", main=unique(rhoA5$event))
dev.off()        
#######################################################
rhoEX<-df[df$event=="exons",]
rhoEX<-rhoEX[order(as.numeric(rhoEX$ne)),]
png("EX_NumOfEventsVsRho.png")
plot(rhoEX$ne,rhoEX$rMinRho, type = "l",  xlab="Num of events", ylab="rho", main=unique(rhoEX$event))
dev.off()        
#######################################################
rhoIR<-df[df$event=="introns",]
rhoIR<-rhoIR[order(as.numeric(rhoIR$ne)),]
png("IR_NumOfEventsVsRho.png")
plot(rhoIR$ne,rhoIR$MinRho, type = "l",  xlab="Num of events", ylab="rho", main=unique(rhoIR$event))
dev.off()        
#all together
png("NumOfEventsVsRho.png")
plot(rhoA3$ne,rhoA3$rho, type = "b",  xlab="Num of events", ylab="rho", main="All types", col="blue", ylim=c(0.3,0.7), xlim=c(20,140))
lines(rhoA5$ne,rhoA5$rho, type = "b", col="red")
lines(rhoEX$ne,rhoEX$rho, type = "b", col="green")
lines(rhoIR$ne,rhoIR$rho, type = "b",  col="yellow")
legend("topright", 
       c("A3","A5","EX","IR"),
       lty = 1, 
       col = c("blue", "red","green","yellow"))
dev.off()
#########################################
  
######################################################
rhoTotal<-df[df$event=="total",]
rhoTotal<-rhoTotal[order(as.numeric(rhoTotal$ne)),]
head(rhoTotal)
tail(rhoTotal)
xrange<-  range(as.numeric(rhoTotal$ne) )
yrange<-range(as.numeric(rhoTotal$MinRho) )
png("total_NumOfEventsVsRho.png")
plot(
  rhoTotal$ne,
  rhoTotal$MinRho, type = "b",  xlab="Num of events", ylab="rho", 
  main=unique(rhoTotal$event), xlim=xrange, ylim=yrange)
dev.off()        
#######################################################

