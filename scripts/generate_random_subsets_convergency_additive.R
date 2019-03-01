  #first generate data: TRUE and RANDOM DATA
###############################################################################
#1:read suppa delta PSI table
setwd("/home/emancini/Dropbox (CRG)/Personal_Estefania/Network/selectedEventsHs2/convTotalSum/")
dir.create("/home/emancini/Dropbox (CRG)/Personal_Estefania/Network/selectedEventsHs2/convTotalSum")
dir2print<- "/home/emancini/Dropbox (CRG)/Personal_Estefania/Network/selectedEventsHs2/convTotalSum/"
###############################################################################
file<-"total_600.txt"

###############################################################################
Oname<-gsub("\\.txt", "", file)
print(Oname)
data<-  read.delim(file,
                 header = T,
                 sep="\t",  
                 stringsAsFactors = F, 
                 dec = ",", row.names = 1);
dim(data)
colnames(data)
nramdomPoints=round(nrow(data)/20)
###############################################################
onlyValues<-data[,18:ncol(data)]
onlyValues$SD<-NULL
onlyValues$RANGE<-NULL
onlyValues[1:5,]
head(onlyValues); 
dim(onlyValues)#30*305
#match the names of the KDs:
psi<-onlyValues; 
colnames(psi)[colnames(psi)=="LENG1_b"]<-"LENG1"
colnames(psi)[colnames(psi)=="RBM17con"]<-"RBM17"
colnames(psi)[colnames(psi)=="HFM1_b"]<-"HFM1"
colnames(psi)[colnames(psi)=="CCDC12_b"]<-"CCDC12"
colnames(psi)[colnames(psi)=="CDC5L_b"]<-"CDC5L"
colnames(psi); dim(psi)#30*305
psiO<-psi
set.seed(1234)
for (i in 1:nramdomPoints)
  {
      total<-1:nrow(psi)
      if (i==1){
      ii<-sample(total, 20)
      delta<-psi[ii,]
      deltascaled<-scale(delta)    #scaled by KDs, maintain the profiles
      eventsscaled<-t(scale(t(deltascaled)))
      #####################################################
      name<-paste(Oname, nrow(delta), sep = "_")
      file<-paste(dir2print, paste(name, "eventscaled.tab", sep="_"), sep="/")
      fileQsub<-paste(name, "eventscaled.tab", sep="_")
      fileQ<-paste(dir2print, paste(name, "qsub.sh", sep="."), sep="/") 
      qsub<-paste(
        paste(
          paste("qsub -q long-sl7 -V -cwd  -N fdr  -l virtual_free=100G -l h_rt=06:00:00 -b y -pe smp 10 Rscript --vanilla ../Network/scripts/Net-fdr_CL.R -s 0.2 -e 0.9 -i 0.05  -r 50 -b /no_backup/jvalcarcel/emancini/Network/scripts/ -c 10 -f /no_backup/jvalcarcel/emancini/Network-testData/",
                fileQsub, sep=""),
          "-n", sep=" "), 
        name, sep=" ")
      write.table(eventsscaled,file,  sep="\t")
      write( qsub,fileQ)  
      total<-total[-ii]
      print(i)
      print(ii)
      i=i+1
    }
  else {

    bb<-sample(total, 20)
    ii<-c(ii,bb)
    delta<-psi[ii,]
    deltascaled<-scale(delta)    #scaled by KDs, maintain the profiles
    eventsscaled<-t(scale(t(deltascaled)))
    #####################################################
    name<-paste(Oname, nrow(delta), sep = "_")
    file<-paste(dir2print, paste(name, "eventscaled.tab", sep="_"), sep="/")
    fileQsub<-paste(name, "eventscaled.tab", sep="_")
    fileQ<-paste(dir2print, paste(name, "qsub.sh", sep="."), sep="/") 
    qsub<-paste(
      paste(
        paste("qsub -q long-sl7 -V -cwd  -N fdr  -l virtual_free=100G -l h_rt=06:00:00 -b y -pe smp 10 Rscript --vanilla ../Network/scripts/Net-fdr_CL.R -s 0.2 -e 0.9 -i 0.05  -r 50 -b /no_backup/jvalcarcel/emancini/Network/scripts/ -c 10 -f /no_backup/jvalcarcel/emancini/Network-testData/",
              fileQsub, sep=""),
        "-n", sep=" "), 
      name, sep=" ")
    write.table(eventsscaled,file,  sep="\t")
    write( qsub,fileQ)  
    total<-total[-ii]
    i=i+1
    print(i)
    print(ii)
  }
}
getwd()
    
