#first generate data: TRUE and RANDOM DATA
library(QUIC)
library(igraph)
library(corrplot)
library("RColorBrewer")
library(scales)
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/CRobCor.R")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/functions_fdr_noMC.R")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/CreateGraph.R")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/CentralityRanking.R")
source("~/Dropbox (CRG ADV)/Personal_Estefania/Network/Network/scripts/Vscale.R")

###############################################################################
getMeanDeg<-function(start, end, gap, Cdouble, name)
{
  
  df<-data.frame(rho=seq(start,end,gap))
  for (i in 1:length(seq(start,end,gap))) {
    print(i)
    gListDouble<-CreateGraph(Cdouble,rho=df$rho[i]) #
    gList<-gListDouble
    g<-gList[[1]]
    deg <- degree(g, mode="all")
    df$md[i]<-mean(deg)
    print(mean(deg))
    file=paste(paste("MeanDegree", name, sep="_"), "png", sep=".")
    png(file)
    plot(df$rho, df$md, xlim=c(0,1), ylim=c(0,max(df$md)), type = "l", main=paste("Mean degree", name))
    abline(h=6)
    dev.off()
    
  }
  
  return(df)
}

###############################################################################
all<-
read.delim("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/lastVtHsa/selected/Alt3_30.txt",
header = T,
sep="\t",  
stringsAsFactors = F, 
dec = ",", 
row.names = 1);
dim(all); head(all)
data<-all
###############################################################
onlyValues<-data[,17:ncol(data)]
dim(onlyValues)
onlyValues$SD<-NULL
onlyValues$RANGE<-NULL
head(onlyValues); 
#match the names of the KDs:
psi<-onlyValues; 
colnames(psi)[colnames(psi)=="LENG1_b"]<-"LENG1"
colnames(psi)[colnames(psi)=="RBM17con"]<-"RBM17"
colnames(psi)[colnames(psi)=="HFM1_b"]<-"HFM1"
colnames(psi)[colnames(psi)=="CCDC12_b"]<-"CCDC12"
colnames(psi)[colnames(psi)=="CDC5L_b"]<-"CDC5L"
colnames(psi); dim(psi)#30*305
write(colnames(psi))
####################
#? psiCompatible -> colnames labchip
#
labchipKDs<-read.table("/home/emancini/Dropbox (CRG ADV)/Personal_Estefania/Network/comparison_old_new/compatibledata/newNames_t_Zvals.tab", 
                       stringsAsFactors = F, 
                       nrows = 1,header=T)
psiCompatible<-colnames(labchipKDs)
psiCompatible[1:10]
length(psiCompatible)#267
class(psiCompatible)
ii<-match(psiCompatible, colnames(psi)) 
psiCompatible[which(is.na(ii))]
ii<-ii[!is.na(ii)]
colnames(psi)
length(ii)#264
#? correcto o no? dejamos todos los KDds o solo los que son compatibles con la primera network?
length(psiCompatible)#267 3 son NA
psifinal<-psi[,ii]
psifinal<-psi
dim(psifinal)#305
####################################################################
delta<-psifinal
#remove controls from matrix before scaling
#  Column Scaling -> Only Shape not Scale of KD effect is important. This is essential since for example not all KDs have the same efficiecny 
deltascaled<-scale(delta)    #scaled by KDs, maintain the profiles
eventsscaled<-t(scale(t(deltascaled)))
########################################
getwd()
setwd("/home/emancini/Dropbox (CRG)/Personal_Estefania/Network/subsets")
write.table(eventsscaled,"best30Alt3_all_eventscaled.tab", sep="\t")
getwd()
# Rho was set to achieve FDR < 5%
#para cada nueva nwtwork tengo que hacer esto:
#IR All
M<-read.table("best30IR_all_eventscaled.tab"); head(M)
Cdouble <- CRobCor(M)
res<-getMeanDeg(start=0.5,1,gap=0.05,Cdouble, name="IR")
gListDouble<-CreateGraph(Cdouble,rho=0.85) 
gList<-gListDouble
g<-gList[[1]]

#IR comp
M<-read.table("best30IR_comp_eventscaled.tab"); head(M)
Cdouble <- CRobCor(M)
gListDouble<-CreateGraph(Cdouble,rho=0.48) 
gList<-gListDouble
g<-gList[[1]]
###Alt5 All
M<-read.table("best30Alt5_all_eventscaled.tab"); head(M)
Cdouble <- CRobCor(M)
res<-getMeanDeg(start=0.5,0.9,gap=0.05,Cdouble, name="Alt5")
gListDouble<-CreateGraph(Cdouble,rho=0.78) 
gList<-gListDouble
g<-gList[[1]]
###Alt5 Comp
M<-read.table("best30Alt5_comp_eventscaled.tab"); head(M)
Cdouble <- CRobCor(M)
gListDouble<-CreateGraph(Cdouble,rho=0.45) 
gList<-gListDouble
g<-gList[[1]]
##############################################
#repeat new exons
M<-read.table("newExons_eventscaled.tab"); head(M)
Cdouble <- CRobCor(M)
gListDouble<-CreateGraph(Cdouble,rho=0.55) 
gList<-gListDouble
g<-gList[[1]]
##############################################################
#repeat 
M<-read.table("newExonsAll_eventscaled.tab"); head(M)
Cdouble <- CRobCor(M)
gListDouble<-CreateGraph(Cdouble,rho=0.55) 
gList<-gListDouble
g<-gList[[1]]
############################################
#merge all (escaled) events in 1 big matrix 120 events:
Mscaled<-rbind(ME,MI,M5,M3)
dim(Mscaled)
write.table(Mscaled,"best120_all_eventscaled.tab", sep="\t")
#randomly select 30 events for this matrix
MscaledSampled<-Mscaled[sample(1:120, 30),]
head(MscaledSampled)
write.table(MscaledSampled,"best30_random_all_eventscaled.tab", sep="\t")
############################################
#armar la matrix total y luego escalarla:
e30<-
  read.delim("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/lastVtHsa/selected/exons_30.txt",
             header = T,
             sep="\t",  
             stringsAsFactors = F, 
             dec = ",", 
             row.names = 1);
dim(e30);head(e30)
i30<-
  read.delim("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/lastVtHsa/selected/introns_30.txt",
             header = T,
             sep="\t",  
             stringsAsFactors = F, 
             dec = ",", 
             row.names = 1);
dim(i30);head(i30)
a530<-
  read.delim("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/lastVtHsa/selected/Alt5_30.txt",
             header = T,
             sep="\t",  
             stringsAsFactors = F, 
             dec = ",", 
             row.names = 1);
dim(a530);head(a530)
a330<-
  read.delim("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/lastVtHsa/selected/Alt3_30.txt",
             header = T,
             sep="\t",  
             stringsAsFactors = F, 
             dec = ",", 
             row.names = 1);
dim(a330);head(a330)
###############################################################
data<-rbind(e30,i30,a530,a330)
onlyValues<-data[,17:ncol(data)]
dim(onlyValues)
onlyValues$SD<-NULL
onlyValues$RANGE<-NULL
#match the names of the KDs:
psi<-onlyValues; 
colnames(psi)[colnames(psi)=="LENG1_b"]<-"LENG1"
colnames(psi)[colnames(psi)=="RBM17con"]<-"RBM17"
colnames(psi)[colnames(psi)=="HFM1_b"]<-"HFM1"
colnames(psi)[colnames(psi)=="CCDC12_b"]<-"CCDC12"
colnames(psi)[colnames(psi)=="CDC5L_b"]<-"CDC5L"
colnames(psi); dim(psi)#120*305
####################
psifinal<-psi
dim(psifinal)#120*305
delta<-psifinal
deltascaled<-scale(delta)    #scaled by KDs, maintain the profiles
eventsscaled<-t(scale(t(deltascaled)))
########################################
setwd("/home/emancini/Dropbox (CRG)/Personal_Estefania/Network/subsets")
write.table(eventsscaled,"best120_all_eventscaled_together.tab", sep="\t")
eventsscaledSampled<-eventsscaled[sample(1:120, 30),]
head(eventsscaledSampled); dim(eventsscaledSampled)#30*305
######################################################
#original
M<-read.table("eventscaled.tab"); head(M)
Cdouble <- CRobCor(M)
res<-getMeanDeg(start=0.4,0.7,gap=0.05,Cdouble, name="RNASEQ_comp")
gListDouble<-CreateGraph(Cdouble,rho=0.52) #
gList<-gListDouble
g<-gList[[1]]

#3alt all
M<-read.table("best30Alt3_all_eventscaled.tab"); head(M)
Cdouble <- CRobCor(M)
res<-getMeanDeg(start=0.5,0.8,gap=0.05,Cdouble, name="Alt3")
gListDouble<-CreateGraph(Cdouble,rho=0.75) #
gList<-gListDouble
g<-gList[[1]]

#3alt comp
M<-read.table("best30Alt3_comp_eventscaled.tab"); head(M)
Cdouble <- CRobCor(M)
gListDouble<-CreateGraph(Cdouble,rho=0.5) #
gList<-gListDouble
g<-gList[[1]]

##120
M<-read.table("best120_all_eventscaled.tab"); head(M)
dim(M)
Cdouble <- CRobCor(M)
gListDouble<-CreateGraph(Cdouble,rho=0.4) 
gList<-gListDouble
g<-gList[[1]]

setwd("/home/emancini/Dropbox (CRG)/Personal_Estefania/Network/subsets")
M<-read.table("best120_all_eventscaled_together.tab"); head(M)
dim(M)
Cdouble <- CRobCor(M)
gListDouble<-CreateGraph(Cdouble,rho=0.45) 
gList<-gListDouble
g<-gList[[1]]

setwd("/home/emancini/Dropbox (CRG)/Personal_Estefania/Network/subsets")
M<-read.table("best30_random_all_eventscaled.tab"); head(M)
dim(M)
Cdouble <- CRobCor(M)
gListDouble<-CreateGraph(Cdouble,rho=0.58) 
gList<-gListDouble
g<-gList[[1]]

setwd("/home/emancini/Dropbox (CRG)/Personal_Estefania/Network/subsets")
M<-read.table("best30_random_all_eventscaled_together.tab"); head(M)
dim(M)
Cdouble <- CRobCor(M)
res<-getMeanDeg(start=0.5,0.8,gap=0.05,Cdouble, name="Random")
gListDouble<-CreateGraph(Cdouble,rho=0.66) 
gList<-gListDouble
g<-gList[[1]]
