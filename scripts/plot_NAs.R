 #finalMatrixCOmparison
library(data.table)
library(stringr)
setwd("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/newNAs/")
getwd()
luisaData <- read.csv("INCLUSION_LEVELS_FULL-Hsa7-hg19.tab", 
                header=T, 
                sep="\t", 
                stringsAsFactors = F, row.names =2 )
dim(luisaData);luisaData[1:5,1:5]
onlyValues<-luisaData[,6:ncol(luisaData)]
Qvals <- grep("\\.Q$", colnames(onlyValues))
length(Qvals)
qualValues <- onlyValues[, Qvals]; dim(qualValues)
head(qualValues)
psiValues  <- onlyValues[,-Qvals]; dim(psiValues)
head(psiValues)

finalValues<-matrix(0, ncol=ncol(psiValues), nrow=nrow(psiValues))
colnames(finalValues)<-colnames(onlyValues)[-Qvals]
dim(finalValues);head(finalValues)
dim(qualValues);head(qualValues)
#now compare the matrixes

#first replace Ns values
Ns_comma<-t(apply(qualValues, 1, function(x)
{ str_count(x, "N,")}))
NsByKds<-colSums(Ns_comma); 
barplot(NsByKds, xlim=c(0,8), ylim=c(0,400000))
NsByEvents<-rowSums(Ns_comma); 
plot(density(NsByEvents, breaks = 20, ylim=c(0,250000)))
indexes_N<-which(Ns_comma !=0, arr.in=TRUE)
finalValues[indexes_N]<-Ns_comma[indexes_N];head(finalValues)
originalNA<-is.na(psiValues); length(originalNA)#2549043
correctedPsiVal<-psiValues; dim(psiValues)[1]*dim(psiValues)[2]#2549043
indexes_NewN2<-which(Ns_comma ==2, arr.in=TRUE); length(indexes_NewN2)/2#113576
indexes_NewN3<-which(Ns_comma ==3, arr.in=TRUE);length(indexes_NewN3)/2#726229
#correctedPsiVal[indexes_NewN1]<-"NAnew1"
correctedPsiVal[indexes_NewN2]<-"NAnew2"
correctedPsiVal[indexes_NewN3]<-"NAnew3"
correctedPsiVal[originalNA]<-"NAold"
fullNAS<-cbind(luisaData[, 1:5], correctedPsiVal)
write.table(fullNAS,"Luisa_PsiValues_NA3_NA2_newNames.tab", sep="\t", col.names = NA )
#####################################################

##################################################################
sf3b1h3 <- read.csv("INCLUSION_LEVELS_FULL-Hsa12-hg19.tab", 
                      header=T, 
                      sep="\t", 
                      stringsAsFactors = F, row.names = 2)
dim(sf3b1h3);sf3b1h3[1:5,6:10]
onlyValuesSF3<-sf3b1h3[,6:ncol(sf3b1h3)]
#####################################################
#then replace NAs, in order not to overlap
indexes_NA<-which(is.na(psiValues), arr.in=TRUE); #NAS
(dim(indexes_NA)[1]*dim(indexes_NA)[2]) /(dim(psiValues)[1]*dim(psiValues)[2])*100#67%
finalValues[indexes_NA]<-NA
######################################
#select the events:
orig_events<-read.table("/home/emancini/Dropbox (CRG)/GosiaAndEstefi/AS_317_OK/Original_Network_events.tsv .csv", header=T, sep="\t")
head(orig_events)
rownames(finalValues)<-row.names(onlyValues)
dim(finalValues); dim(fullNAS)

ii<-match(orig_events$EVENT, rownames(finalValues))

dfOrigEvents<-finalValues[ii,]; head(dfOrigEvents)
eventsIDs<-paste(rownames(finalValues)[ii],
                              1:length(rownames(finalValues)[ii]), sep="_")
rownames(dfOrigEvents)<-paste(fullNAS$GENE[ii], eventsIDs, sep=":")
#############################################################
#SD controles
colnames(dfOrigEvents)<-colnames(finalValues)
forPlot<-dfOrigEvents
#row.names(forPlot)<-rownames(dfOrigEvents)
row.names(forPlot)<-1:31
colnames(forPlot)<-colnames(dfOrigEvents)
length(which(is.na(forPlot)))#681
forPlot[is.na(forPlot)] <-4
forPlot[1:5,1:9]
#SD controles
forPlotMelted<-melt(forPlot)
head(forPlotMelted)
table(forPlotMelted$value)
#xlabels
eventName<-row.names(dfOrigEvents)
head(dfOrigEvents)
head(forPlotMelted)
table(forPlotMelted$value)
####################################
#pdf('withNs_hor_large.pdf', width = 80, height = 10)
pdf("withNs_hor_A4.pdf", width = 11.69,  height = 8.27)
  ggplot(forPlotMelted, 
         aes(x = Var2, y = Var1, fill = factor(value)))  +
    theme_minimal() +
    geom_tile(position = "identity",  
              color = "grey") +
    scale_y_continuous(breaks=c(1:31), labels=eventName,
                       limits=c(0,32) , expand = c(0,-0.5)) +
    scale_x_discrete(position = "top") +
    scale_fill_manual(values=c("white","blue","green","pink","red"))+
    theme(
      axis.text.x = element_text(angle = 90, hjust=0, vjust=-5,size=3),
      axis.ticks.length=unit(0,"cm"),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border=element_blank(),
      legend.position="none",
      axis.line=element_blank()
    )
  dev.off()

  
  