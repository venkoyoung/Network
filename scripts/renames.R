library(limma)
setwd("~/Dropbox (CRG)/Personal_Estefania/Network/comparison_old_new")
t_Zvals<-read.table("t_Zvals_LabChip.tab"); head(t_Zvals)#ya estan corregidos los nombres
dim(t_Zvals)
#######################################################################################
all331<-read.table("~/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/AS_331/FcorrectedPsiValues_NA3_newNames_b.tab",
                   header = T, sep="\t", row.names = 1)
#Splicing Annotation Data:
dim(all331)
data <- all331[, ! colnames(all331) %in% 
                 c("AA2","AA1","CCDC12","C1orf55",
                   "C1orf55_b","CDC5L","HFM1","LENG1",
                   "RBM17","PPIL1","SRRT"),
               drop = F]
dim(data)
data[1:5,1:8]
all331[1:5,8:10]
#remove bad samples
colnames(all331)
#select the events:
orig_events<-read.table("~/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/AS_317_OK/Original_Network_events.tsv .csv", header=T, sep="\t")
head(orig_events)
ii<-match(orig_events$EVENT, rownames(data))
dfOrigEvents<-data[ii,]; head(dfOrigEvents)
head(data)
eventsIDs<-paste(rownames(data)[ii],
                 1:length(rownames(data)[ii]), sep="_")
rownames(dfOrigEvents)<-paste(data$GENE[ii], eventsIDs, sep=":")
onlyValues<-dfOrigEvents[,7:ncol(dfOrigEvents)]
getwd()
dim(onlyValues)
write.table(onlyValues,"320onlyValuesOriginalEventsRNASEQ.txt", sep="\t")
psi<-onlyValues
head(psi)
dim(psi)
dim(t_Zvals)
###############################################################
allNames<-union(colnames(psi),colnames(t_Zvals))
auxdf<-data.frame(rnaseq=allNames%in%colnames(psi),
                  labchip=allNames%in%colnames(t_Zvals), row.names = allNames)
vennDiagram(auxdf)
####################################################
conversion<-read.table("~/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/ensembl/conversionTable_samples_ENSEMBL_VAST.txt", 
                       sep="\t", 
                       header=T, stringsAsFactors = F)
ii<-match(colnames(t_Zvals), conversion$SampleoOriginalName)
conversion$SampleoOriginalName[conversion$SampleoOriginalName =="ESRP1"]
write.table(conversion[(ii),], "conversion_labchip_vast.txt", sep="\t")
labNewNames<-t_Zvals# viene del chip

length(as.character(conversion$SampleoOriginalName[ii]))
colnames(labNewNames)<-as.character(conversion$VT[ii])
colnames(labNewNames)[is.na(colnames(labNewNames))]<-colnames(t_Zvals)[ is.na(ii) ]
colnames(labNewNames)[colnames(labNewNames)=="RBM35A"]<-"ESRP1"
colnames(labNewNames)[colnames(labNewNames)=="U5.116KD"]<-"EFTUD2"
colnames(labNewNames)[colnames(labNewNames)=="U5.200KD"]<-"SNRNP200"
colnames(labNewNames)[colnames(labNewNames)=="IMP.3"]<-"IGF2BP3"
head(labNewNames); dim(labNewNames)
write.table(labNewNames, file="profiles/newNames_t_Zvals.tab", sep="\t")#OK
###########################################################
allNames<-union(colnames(psi),colnames(labNewNames))
auxdf<-data.frame(rnaseq=allNames%in%colnames(psi),
                  labchip=allNames%in%colnames(labNewNames), row.names = allNames)
vennDiagram(auxdf)
###########################################################
colnames(psi)[colnames(psi)=="LENG1_b"]<-"LENG1"
colnames(psi)[colnames(psi)=="RBM17con"]<-"RBM17"
colnames(psi)[colnames(psi)=="HFM1_b"]<-"HFM1"
colnames(psi)[colnames(psi)=="CCDC12_b"]<-"CCDC12"
colnames(psi)[colnames(psi)=="CDC5L_b"]<-"CDC5L"
colnames(psi)[colnames(psi)=="CCDC12_b"]<-"CCDC12"
colnames(onlyValues); dim(onlyValues)
length(which(onlyValues=="NAold"))#575
length(which(onlyValues=="NAnew3"))#833
psi[psi=="NAold"]<-NA
dim(psi)[1]*dim(psi)[2]
psi[psi=="NAnew3"]<-NA
ll<-data.frame(matrix(unlist(apply(psi,1,as.numeric)), ncol=320, byrow = T))
colnames(ll)<-colnames(psi)
rownames(ll)<-row.names(psi)
ll$IKm<-ave(ll$IK, ll$IKcon, na.rm = T)
ll$SF3B1m <- ave(ll$SF3B1, ll$SF3B1con, na.rm = T)
ll$SMU1m  <- ave(ll$SMU1, ll$SMU1con, na.rm = T)
ll$SRPK2m  <- ave(ll$SRPK2, ll$SRPK2_b, na.rm = T)
ll$XAB2m  <- ave(ll$XAB2, ll$XAB2_b, na.rm = T)
ll$PRPF8m <- ave(ll$PRPF8, ll$PRPF8con, na.rm = T)
ll$CWC22m <- ave(ll$CWC22, ll$CWC22_b, na.rm = T)
###########################################################
#remover los origininales:
psifinal<-ll[, ! colnames(ll) %in% 
               c(  "SF3B1","SF3B1con",
                    "SMU1","SMU1con",
                    "SRPK2","SRPK2_b",
                    "XAB2","XAB2_b",
                   "PRPF8","PRPF8con",
                   "CWC22","CWC22_b"),
                 drop = F]
colnames(psifinal)<-gsub("m","",colnames(psifinal))
###########################################################
head(psifinal)
allNames<-union(colnames(psifinal),colnames(labNewNames))
auxdf<-data.frame(rnaseq=allNames%in%colnames(psifinal),
                  labchip=allNames%in%colnames(labNewNames), row.names = allNames)
vennDiagram(auxdf)
#select by columns:
auxdf[rowSums(auxdf)==2,]
kdsI<-rownames(auxdf)[rowSums(auxdf)==2]
psifinal[,kdsI]
length(kdsI)
psifinal<-cbind(psifinal[,kdsI], psi[,1:7])
dim(psifinal)
length(which(psifinal=="NAold"))#0
length(which(psifinal=="NAnew3"))#0
dim(t_Zvals)
dim(psifinal)
###########################################################
write.table(psifinal, "profiles/psi_only_labChipKDS.tab", quote = F, sep="\t")
getwd()
save(psifinal,file="psifinalOnlyLabChipKD.RData")
#desde aca sigo con el scrit SplNetRNASEQ
################
genes<-matrix(unlist(strsplit(row.names(psifinal), ":")), ncol=2, byrow = T)[,1]
ge220<-read.table("/home/emancini/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/GE_hg19_noBadSamples/lrt_table_hg19.tab", header=T, row.names = 1)
ge220[1:5,1:5]
#rename colnames:
colnames(ge220)<-gsub("logFC.group","",colnames(ge220))
#rename rownames:
dim(ge220)
ge220[1:5,]
conversion[1:5,1:5]
compconversion<-read.table("~/Dropbox (CRG)/Personal_Gosia/Shared/GosiaAndEstefi/ensembl/completeConversionTable_ENSEMBL_VAST.txt", 
                       sep="\t", 
                       header=T, stringsAsFactors = F)
head(compconversion)
ii<-match(rownames(ge220), as.character(compconversion$FGeneID))
ge220$symbol<-compconversion$VT[ii]
dim(ge220)
ge220[1:5,300:313]  
#select the rows;
gg<-match(genes, compconversion$VT)
ee<-match(compconversion$FGeneID[gg], row.names(ge220))
geneEvents<-ge220[ee[!is.na(ee)],]
###############################################
colnames(geneEvents)
colnames(geneEvents)[colnames(geneEvents)=="LENG1_b"]<-"LENG1"
colnames(geneEvents)[colnames(geneEvents)=="HFM1_b"]<-"HFM1"
colnames(geneEvents)[colnames(geneEvents)=="CCDC12_b"]<-"CCDC12"
colnames(geneEvents)[colnames(geneEvents)=="CDC5L_b"]<-"CDC5L"
colnames(geneEvents)[colnames(geneEvents)=="PPIL1"]
kds<-match(colnames(labNewNames),
           colnames(geneEvents)  )
cbind(colnames(labNewNames), kds)
dim(geneEvents)
geneEventsKds<-geneEvents[ ,kds[!is.na(kds)]]
dim(geneEventsKds)
dim(geneEventsKds)
ss<-match(rownames(geneEventsKds), compconversion$FGeneID)
which(is.na(ss))#e should remove them, because are duplicates
geneId<-rownames(geneEventsKds)
geneId[20]<-"ENSG00000115414"
geneId[24]<-"ENSG00000067225"
geneId[25]<-"ENSG00000113648"
geneId[27]<-"ENSG00000089685"
geneEventsKdsNoDup<-geneEventsKds[-which(is.na(ss)),]
gsy<-match(rownames(geneEventsKdsNoDup), compconversion$FGeneID)
row.names(geneEventsKdsNoDup)<-compconversion$VT[gsy]
setwd("../ge/")
write.table(geneEventsKdsNoDup, file="profiles/GeneEventsKds.tab", sep="\t", quote = F)
getwd()
dim(geneEventsKdsNoDup)
colnames(geneEventsKdsNoDup)
allNamesGE<-union(union(colnames(psifinal),colnames(labNewNames)), colnames(geneEventsKdsNoDup))
auxdf<-data.frame(rnaseq=allNamesGE%in%colnames(psifinal),
                  labchip=allNamesGE%in%colnames(labNewNames),
                  GE=allNamesGE%in%colnames(geneEventsKdsNoDup),row.names = allNamesGE)
vennDiagram(auxdf)
auxdf[rowSums(auxdf)==1,]
#######################################we have to select this 264 common in all the samples:
#GE is OK geneEventsKdsNoDup
#AS is OK psifinal
#lapchip remove: #PPIL1 C1orf55 
#table: labNewNames
allNamesGE<-union(union(colnames(deltascaled),colnames(labNewNames)), colnames(geneEventsKdsNoDup))
auxdf<-data.frame(rnaseq=allNamesGE%in%colnames(deltascaled),
                  labchip=allNamesGE%in%colnames(labNewNames),
                  GE=allNamesGE%in%colnames(geneEventsKdsNoDup),row.names = allNamesGE)
vennDiagram(auxdf)




