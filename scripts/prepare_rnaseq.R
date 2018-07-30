library(data.table)
library(ggplot2)
all<-read.table("~/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/lastVtHsa/VT_Hg19_NA3_NAI.tab",
                header = T, sep="\t", row.names = 1);
dim(all)
#remove bad samples
data <- all[, ! colnames(all) %in% 
                  c("AA2",
                    "AA1",
                    "CCDC12",
                    "C1orf55",
                    "C1orf55_b",
                    "CDC5L",
                    "HFM1",
                    "LENG1",
                    "RBM17",
                    "PPIL1",
                    "SRRT"),
                  drop = F]
#remove 11 samples
dim(data)
#select the events:
orig_events<-read.table("~/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/AS_317_OK/Original_Network_events.tsv .csv", header=T, sep="\t")
head(orig_events)
ii<-match(orig_events$EVENT, rownames(data))
dfOrigEvents<-data[ii,]; 
head(dfOrigEvents);dim(dfOrigEvents)
eventsIDs<-paste(rownames(data)[ii],1:length(rownames(data)[ii]), sep="_")
rownames(dfOrigEvents)<-paste(data$GENE[ii], eventsIDs, sep=":")
onlyValues<-dfOrigEvents[,7:ncol(dfOrigEvents)]
psi<-onlyValues; head(onlyValues); dim(onlyValues)
#match the names of the KDs:
####################
colnames(onlyValues); dim(onlyValues)
length(which(onlyValues=="NAold"))#575
length(which(onlyValues=="NAnew3"))#833
length(which(onlyValues=="NAI"))#833
#############################################################
dim(onlyValues)[1]*dim(onlyValues)[2]#9920
dim(psi)[1]*dim(psi)[2]
575/9920*100#5.76
833/9920*100#8.39
(575+833)/9920*100#15%
psi[psi=="NAold"]<-NA
psi[psi=="NAnew3"]<-NA
head(psi)
length(which(is.na(psi)))#575
dim(psi)
#desde aca vengo con el script renames
#llamo a la funcion:
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/plot_NAs.R")
psi<-as.matrix(psi)
setwd("~")
####################################################################
 