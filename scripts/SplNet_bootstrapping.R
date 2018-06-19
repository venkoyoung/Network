setwd("~/Dropbox (CRG)/Network/")
library(vegan)
library("impute") #Availabe from bioconductor
library("glasso")
library("QUIC")
library("RColorBrewer")
library("igraph")
library("SDMTools")
########## Functions
source("scripts/CentralityRanking.R")
source("scripts/CRobCor.R")
source("scripts/Vscale.R")
#Network Data. These are already transformed to Z-scores + Pvalues. This file is generated from the prepare_z_extra_data_2016.pl script.
Ztable<-read.table("LABCHIPS/labchip37_26_06_13_bare_Z_P.tsv",as.is=TRUE,na.strings="NAn",header=TRUE,row.names=1,sep="\t",quote=""); 
Zmatrix<-as.matrix(Ztable[,seq(1,(ncol(Ztable)-1),2)]);  
Pmatrix<-as.matrix(Ztable[,seq(2,ncol(Ztable),2)]);		
#Remove Uninformative genes. These are genes whose inclusion does not change in the HeLa context .
#Columns events
Zmatrix<-(Zmatrix[,-c(20,24,28,36)]);             #Do not include VEGFA (20th), SYK (24th), CND1 (28th), ???GADD45A (31st), LMNA (36)
Pmatrix<-as.matrix(Pmatrix[,-c(20,24,28,36)]);		#Do not include VEGFA (20th), SYK (24th), CCND1 (28th), ???GADD45A (31st), LMNA (36)
# Remove Sparse KD data and impute remaining missing values:
Zvals<-Zmatrix; head(Zvals); dim(Zvals)
hist(rowMeans(is.na(Zvals))*100)
nais<-is.na(Zvals); Zvals[nais]
SparseRows<-which(rowMeans(is.na(Zvals))*100>30)#Row numbers with >30% missing values
if (length(SparseRows)>0){
  Z<-Z[-c(SparseRows),]}
ZvalsI<-impute.knn(as.matrix(Zvals),k=round(sqrt(nrow(Zvals))*0.25),rowmax=0.5)$data
t_Zvals=t(ZvalsI)
#Remove conditions not corresponding to KDs (e.g untransfected, mocks...). This is based on a collection of suffixes/prefixes, since data kept being added to the Labchip this kept expanding...
NonKDs=c("Untransfected",
         "untransfected",
         "Mock","CN ","CN_",
         "CL_","GFP_",
         "Control","_CN",
         "DMSO","UNT_","CTRL_",
         "ctr_4h","ctr_24h",
         "ctr_37h","Contr_",
         "Scrambled","0h","rep")
for (i in 1: length(NonKDs)){
  word=NonKDs[i]
  if (length(grep(word,rownames(Zvals) ) )>0  )  {t_Zvals<-t_Zvals[1:ncol(Zvals),-c(grep(word,colnames(t_Zvals)))]       }   
}
#Knock Down (Column) Scaling. 
t_Zvals<-scale(t_Zvals); dim(t_Zvals) 
C<-CRobCor(t_Zvals)
TPList<-vector("list", 1)
TPList[[1]]<-C

var1<-seq(from=0.3,  to=0.7, by=0.02)
TrueAdjMatrix<-getADJRHO(TPList); length(TrueAdjMatrix)#10
##################################################
totalTP<-vector()
length(totalTP)
class(totalTP)
##################################################
for (i in 1:length(TrueAdjMatrix))  {
  ll<-lapply(TrueAdjMatrix[[i]],
             function(x)
             {
               TruePos<-length(which( unlist(x)>0 ) )
               return(TruePos)
             }
  )   
  totalTP<-c(totalTP, ll)
}
plot(var1,unlist(totalTP), type="l")
##################################################
length(which((unlist(TrueAdjMatrix[[1]])>0)))

###############################################
MList <- vector("list", 10)
for (i in 1:10) {
  MList[[i]]<-t(apply(t_Zvals,1, sample))
}
################################################
RList <- vector("list", length(MList))
RList<-lapply(MList, CRobCor)
#################################################
RandomeList<-getEdgesByRHO(RList); length(RandomeList)#10
RandomAdjMatrix<-getADJRHO(RList); length(RandomAdjMatrix)#10
#############compare matrix by matrix for each RHO


final<-vector()
for (i in 1:length(RandomAdjMatrix))  {
      ll<-lapply(RandomAdjMatrix[[i]],
             function(x)
             {
               
           falsePos<- length( which((unlist(x)>0 )) %in% (unlist(TrueAdjMatrix[[i]]) >0)  )
           return(falsePos)
            }
       )   
      final<-c(final, ll)
}
tpByho<-matrix(ncol=10, nrow=21, data=as.numeric(unlist(final)), byrow = T)
class(tpByho)
head(tpByho)
####################################################################
final2<-vector()
for (i in 1:length(RandomAdjMatrix))  {
  ll<-lapply(RandomAdjMatrix[[i]],
             function(x)
             {
               
               falsePos<- length(which((unlist(x)>0 )))
               
               return(falsePos)
             }
  )   
  final2<-c(final2, ll)
}
tpByho2<-matrix(ncol=10, nrow=21, data=as.numeric(unlist(final2)), byrow = T)
class(tpByho2)
rowMeans(tpByho2)
rowMeans(tpByho)

plot(var1,rowMeans(tpByho2), type="l", col="red")
lines(var1,rowMeans(tpByho), type="l", col="blue")
lines(var1,unlist(totalTP), col="green")
ratio=rowMeans(tpByho)/unlist(totalTP)
plot(var1,ratio, col="green", type="l")
min(ratio)

####################################################################
rr<-unlist(lapply(RandomeList, mean))
tt<-unlist(lapply(TPList))
plot(seq(from=0.3,  to=0.7, by=0.02) , rr )
plot(rep(seq(from=0.3,  to=0.7, by=0.02), each=10) , unlist(RandomeList))
var1<-seq(from=0.3,  to=0.7, by=0.02)
plot(var1, rr, col="blue", type="l")
lines(var1, tt, col="red")
#random cae mas rapido
###########################################
plot(var1,  rr/tt, type="l")
abline(h=c(0,0.05), v=0)
######################################################################
#Ahora testeamos con la matriz de adyacencias
