#first generate data: TRUE and RANDOM DATA
###############################################################################
#1:read suppa delta PSI table
setwd("/home/emancini/Dropbox (CRG)/Personal_Estefania/Network/selectedEventsHs2/")
dir.create("/home/emancini/Dropbox (CRG)/Personal_Estefania/Network/selectedEventsHs2/randomA3_30_15")
###############################################################################
file<-"Alt3_150.txt"
nramdomCycles=100
nramdomPoints=30
Oname<-gsub("\\.txt", "", file)
print(Oname)
data<-  read.delim(file,
                 header = T,
                 sep="\t",  
                 stringsAsFactors = F, 
                 dec = ",", row.names = 1);
###############################################################
onlyValues<-data[,18:ncol(data)]
onlyValues$SD<-NULL
onlyValues$RANGE<-NULL
  dim(onlyValues)#30*305
  head(onlyValues); 
  #match the names of the KDs:
psi<-onlyValues; 
  colnames(psi)[colnames(psi)=="LENG1_b"]<-"LENG1"
  colnames(psi)[colnames(psi)=="RBM17con"]<-"RBM17"
  colnames(psi)[colnames(psi)=="HFM1_b"]<-"HFM1"
  colnames(psi)[colnames(psi)=="CCDC12_b"]<-"CCDC12"
  colnames(psi)[colnames(psi)=="CDC5L_b"]<-"CDC5L"
  colnames(psi); dim(psi)#30*305

  set.seed(1234)
  for (i in 1:nramdomCycles)
  {
  ii<-sample(1:nrow(psi), nramdomPoints)
  print(ii)
  delta<-psi[ii,]
  deltascaled<-scale(delta)    #scaled by KDs, maintain the profiles
  eventsscaled<-t(scale(t(deltascaled)))
  #####################################################
  dir2print<- "/home/emancini/Dropbox (CRG)/Personal_Estefania/Network/selectedEventsHs2/randomA3_30_15/"
  name<-paste(Oname, i, sep = "_")
  file<-paste(dir2print, paste(name, "eventscaled.tab", sep="_"), sep="/")
  fileQsub<-paste(name, "eventscaled.tab", sep="_")
  fileQ<-paste(dir2print, paste(name, "qsub.sh", sep="."), sep="/") 
     
  qsub<-paste(
  paste(
      paste("qsub -q long-sl7 -V -cwd  -N fdr  -l virtual_free=100G -l h_rt=06:00:00 -b y -pe smp 10 Rscript --vanilla ../Network/scripts/Net-fdr_CL.R -s 0.3 -e 0.9 -i 0.05 -r 50 -b /no_backup/jvalcarcel/emancini/Network/scripts/ -c 10 -f /no_backup/jvalcarcel/emancini/Network-testData/",
            fileQsub, sep=""),
        "-n", sep=" "), 
    name, sep=" ")
    write.table(eventsscaled,file,  sep="\t")
    write( qsub,fileQ)
    }
  
