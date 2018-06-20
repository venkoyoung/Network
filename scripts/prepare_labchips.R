setwd("~/Dropbox (CRG)/Personal_Estefania/Network/comparison_old_new/")
library(impute)
##################################
Ztable<-read.table("labchip/labchip37_26_06_13_bare_Z_P.tsv",
                   as.is=TRUE,
                   na.strings="NAn",
                   header=TRUE,
                   row.names=1,
                   sep="\t",
                   quote=""); 
#Splicing Annotation Data:
#Split data to a Z-values and a P-values matrix:
Zmatrix<-as.matrix(Ztable[,seq(1,(ncol(Ztable)-1),2)]) ;  #impares # valores
Pmatrix<-as.matrix(Ztable[,seq(2,ncol(Ztable),2) ] );		 #pares
#rows: KD
#Columns events
#Remove Uninformative genes. These are genes whose inclusion does not change in the HeLa context .
#Do not include VEGFA (20th), SYK (24th), CND1 (28th), ???GADD45A (31st), LMNA (36)
#not analyzed because they are not expressed in HEla or they dont change
colnames(Zmatrix)[c(20,24,28,36)]
Zmatrix<-(Zmatrix[,-c(20,24,28,36)]);        
Pmatrix<-as.matrix(Pmatrix[,-c(20,24,28,36)]);		
# Remove Sparse KD data and impute remaining missing values:
Zvals<-Zmatrix; 
dim(Zvals)#269
#impute ROWS-> for KNOCK DOWNS (not for events)
SparseRows<-which(rowMeans(is.na(Zvals))*100>30)#Row numbers with >30% missing values

if (length(SparseRows)>0){
  Zvals<-Zvals[-c(SparseRows),]
}
dim(Zvals)#268
#pisa la original
ZvalsI<-impute.knn(as.matrix(Zvals),
                   k=round(sqrt(nrow(Zvals))*0.25),rowmax=0.5)$data
hist(rowMeans(is.na(Zvals))*100)# desaparecen los NAs
hist(rowMeans(is.na(ZvalsI))*100)# desaparecen los NAs
### Transpose matrix to ASEs (rows) x KDs (columns)
#we have rows=KD factors
#columns=events
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
nrow(Zvals)#268
for (i in 1: length(NonKDs)){
  word=NonKDs[i]
  if (length(grep(word,rownames(Zvals) ) )>0  )  
  {
    t_Zvals<-t_Zvals[1:ncol(Zvals),-c(grep(word,colnames(t_Zvals)))]       
    }   
}
#final cleaned table:
dim(t_Zvals)#39*266
write.table(t_Zvals, file="../data/labchip_cleaned_idata.tab")
