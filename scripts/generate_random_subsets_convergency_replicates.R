#first generate data: TRUE and RANDOM DATA
###############################################################################
#1:read suppa delta PSI table
setwd("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/rangesDSets/MIC/")
dir2print<- "/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/rangesDSets/"
###############################################################################
rangedfiles<-read.table("ranges.txt", header = F, stringsAsFactors = F ,sep="\t")
head(rangedfiles)
dim(rangedfiles)
colnames(rangedfiles)<-"file"
  for (i in 1:nrow(rangedfiles))
  {
    file<-rangedfiles$file[i]  
    print(file)
    nramdomCycles=10
    Oname<-gsub("\\.txt", "", file)
    print(Oname)
    data<-  read.delim(file,
                   header = T,
                   sep="\t",  
                   stringsAsFactors = F, 
                   dec = ",", row.names = 1);
    
    nramdomPoints=seq(20,nrow(data)-10,20 )
  #############################################################
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
  colMeans(psi)  
  set.seed(1234)
  
  for (i in 1:length(nramdomPoints))
  {
    print (nramdomPoints[i])
    for (r in 1:nramdomCycles)
    {
      print(r)
  
      ii<-sample(nrow(psi), nramdomPoints[i])
  #    print(ii)
      delta<-psi[ii,]
      deltascaled<-scale(delta)    #scaled by KDs, maintain the profiles
      eventsscaled<-t(scale(t(deltascaled)))
      #####################################################
      name<-paste(Oname, nramdomPoints[i], sep = "_")
      name<-paste(name, r, sep = "_")
      print(name)
      file<-paste(dir2print, paste(name, "eventscaled.tab", sep="_"), sep="/")
      fileQsub<-paste(name, "eventscaled.tab", sep="_")
      fileQ<-paste(dir2print, paste(name, "qsub.sh", sep="."), sep="/") 
      qsub<-paste(
      paste(
      paste("qsub -q short-sl7 -V -cwd  -N fdr  -l virtual_free=100G -l h_rt=05:00:00 -b y -pe smp 10 Rscript --vanilla ../Network/scripts/Net-fdr_CL.R -s 0.2 -e 0.9 -i 0.05  -r 50 -b /no_backup/jvalcarcel/emancini/Network/scripts/ -c 10 -f /no_backup/jvalcarcel/emancini/Network-testData/",
                fileQsub, sep=""),
            "-n", sep=" "), 
        name, sep=" ")
  
       write.table(eventsscaled,file,  sep="\t")
       write( qsub,fileQ)
  }
  }
  }
getwd()
    
