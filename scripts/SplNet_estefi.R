setwd("~/Dropbox (CRG)/Personal_Estefania/Network/")
##################################
source("https://bioconductor.org/biocLite.R")
biocLite("impute")
install.packages("glasso")
install.packages("QUIC")
install.packages("RColorBrewer")
install.packages("irlba")
install.packages("igraph")
install.packages("SDMTools")
install.packages("png") # for reading in PNGs
install.packages("scales")
library("impute") #Availabe from bioconductor
library("glasso")
library("QUIC")
library("RColorBrewer")
library("igraph")
library("SDMTools")
library("scales")
library("png") # for reading in PNGs
########## Functions
source("scripts/CentralityRanking.R")
source("scripts/CreateGraph.R")
source("scripts/CRobCor.R")
source("scripts/Vscale.R")
#Calculation of custom robust correlation. This tries to discriminate between technical outliers (artefacts) and true biological outliers.
#See supplemental methods of Papasaikas, Tejedor et al for details.
#Arguments are M -> a n x p matrix (n features -e.g ASEs- in the rows, p variables -e.g KDs- in the columns)
#and ct a scalar defining the threshold for finding close neighbors of variables
############ DATA FILTERING AND PREPROCESSING:
#Network Data. These are already transformed to Z-scores + Pvalues. 
#This file is generated from the prepare_z_extra_data_2016.pl script.
Ztable<-read.table("LABCHIPS/labchip37_26_06_13_bare_Z_P.tsv",
                   as.is=TRUE,
                   na.strings="NAn",
                   header=TRUE,
                   row.names=1,
                   sep="\t",
                   quote=""); 
#Splicing Annotation Data:
dim(Ztable)
#Split data to a Z-values and a P-values matrix:
Zmatrix<-as.matrix(Ztable[,seq(1,(ncol(Ztable)-1),2)]);  
head(Zmatrix)
dim(Zmatrix)
Pmatrix<-as.matrix(Ztable[,seq(2,ncol(Ztable),2)]);		
head(Pmatrix)
colnames(Zmatrix)
  #Remove Uninformative genes. These are genes whose inclusion does not change in the HeLa context .
head(Zmatrix); nrow(Zmatrix)
#rows: KD
#Columns events
  colnames(Zmatrix)[c(20,24,28,36)]
#Do not include VEGFA (20th), SYK (24th), CND1 (28th), ???GADD45A (31st), LMNA (36)
#not analyzed because they are not expressed in HEla or they dont change
Zmatrix<-(Zmatrix[,-c(20,24,28,36)]);        
Pmatrix<-as.matrix(Pmatrix[,-c(20,24,28,36)]);		#Do not include VEGFA (20th), SYK (24th), CCND1 (28th), ???GADD45A (31st), LMNA (36)
# Remove Sparse KD data and impute remaining missing values:
Zvals<-Zmatrix; head(Zvals)
hist(rowMeans(is.na(Zvals))*100)
nais<-is.na(Zvals); Zvals[nais]
SparseRows<-which(rowMeans(is.na(Zvals))*100>30)#Row numbers with >30% missing values
if (length(SparseRows)>0){
  Z<-Z[-c(SparseRows),]
}
head(SparseRows)
nrow(Zvals)#269
#pisa la original
ZvalsI<-impute.knn(as.matrix(Zvals),k=round(sqrt(nrow(Zvals))*0.25),rowmax=0.5)$data
hist(rowMeans(is.na(ZvalsI))*100)# desaparecen los NAs
ZvalsI[nais]
### Transpose matrix to ASEs (rows) x KDs (columns)
#we have rows=KD factors
#columns=events
t_Zvals=t(ZvalsI)
head(Zvals)
colMeans(t_Zvals)
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
nrow(Zvals)
for (i in 1: length(NonKDs)){
  word=NonKDs[i]
  if (length(grep(word,rownames(Zvals) ) )>0  )  {t_Zvals<-t_Zvals[1:ncol(Zvals),-c(grep(word,colnames(t_Zvals)))]       }   
}

#Knock Down (Column) Scaling. 
t_Zvals<-scale(t_Zvals)    
#  Column Scaling -> Only Shape not Scale of KD effect is important. This is essential since for example not all KDs have the same efficiecny 
colMeans(t_Zvals)
#Calculate Robust Correlation Matrix
CLabChips <- CRobCor(t_Zvals)
head(C)
save(C, file="C.RData")
class(C)
write.table(C, file="C.tab", sep="\t", col.names = NA)
#Create graph, calculate centrality measures, identify communities and plot
getwd()
setwd("comparison_old_new/")
source("CreateGraph.R")
gList<-CreateGraph(C,rho=0.48)   # Rho was set to achieve FDR < 5%
g<-gList[[1]]
memb<-gList[[2]]
values<-gList[[3]]
mean(values[1,])#6.1  
V(g)
E(g)
mean(deg)
hist(values[1,])#6.1  
summary(values[1,])
##############################################################
#exporta parameters
deg <- degree(g, mode="all")
deg.dist <- degree_distribution(g, cumulative=T, mode="all")
png("cumulativeFreq_048.png")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
  xlab="Degree", ylab="Cumulative Frequency")
dev.off()
########################################3
write.table(t(values), file="048topology.tab", sep="\t")
#add community information:
comps <-memb$membership; head(comps)#dpmde esta cada gen, em que grupo
write.table(comps, file="048memberships.txt", col.names=NA)
########################################3
#match commnuities/descriptores
topo<-t(values)
head(topo)
names(comps)
row.names(topo)[1]
names(comps)[1]
both<-row.names(topo)%in%names(comps)
ii<-match(row.names(topo),names(comps))
names<-names(comps)[ii]
commu<-comps[ii]
names[is.na(names)]<-"notCommu"
commu[is.na(commu)]<-0
row.names(topo)
write.table(data.frame(topo, names, commu), file="allInfo.tab", sep="\t")
########################################3
