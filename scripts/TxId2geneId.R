library(limma)
#benjamin predicion
dfBenInt<-read.table("Re%3a_list_of_genes/output-ensg-1-sigma-input-estefi-ensgs.txt", header=T);head(ben3)
conversion<-read.table("~/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/ensembl/completeConversionTable_ENSEMBL_VAST.txt",
                       sep="\t", header=T, stringsAsFactors = F)
source2vast<-match(as.character(dfBenInt$rbp_ensg), 
                   as.character(conversion$FGeneID)); length(which(is.na(source2vast)))#0
target2vast<-match(as.character(dfBenInt$transcript_ensg),
                   as.character(conversion$FGeneID)); length(which(is.na(target2vast)))#91521

dfBenInt$sourceVT<-conversion$VT[source2vast]
dfBenInt$targetVT<-conversion$VT[target2vast]

dfNet<-data.frame(dfBenInt$sourceVT, dfBenInt$targetVT); dim(dfNet)
dfNet<-na.omit(dfNet); dim(dfNet)
##########################################################
iclipData<-read.table("Re%3a_list_of_genes/output-ensg-eclip-input-estefi-ensgs.txt", header=T, sep="\t")
source2vast<-match(as.character(iclipData$rbp_ensg), 
                   as.character(conversion$FGeneID)); length(which(is.na(source2vast)))#0
target2vast<-match(as.character(iclipData$transcript_ensg),
                   as.character(conversion$FGeneID)); length(which(is.na(target2vast)))#91521
iclipData$sourceVT<-conversion$VT[source2vast]
iclipData$targetVT<-conversion$VT[target2vast]
dfIclip<-data.frame(iclipData$sourceVT, iclipData$targetVT); dim(iclipData)
iclipData<-na.omit(dfIclip); dim(dfIclip)
head(iclipData)
