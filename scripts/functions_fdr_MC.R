randomRows = function(df,n){
  return(df[sample(nrow(df),n),])
}

#functions
createRandomMatrix<-function(TrueM, NumRandomM)
{
  if (is.numeric(TrueM))
  {
    MList <- vector("list", NumRandomM)
    for (i in 1:NumRandomM) {
      MList[[i]]<-t(apply(TrueM,1, sample))
    }
    return(MList)
  }
  else
  {
    print("You are not entering a Numerical matrix")
  }
  
}
####################################################################
getEdges<- function(C, rho)
{
  RHO<-matrix(data=rho,nrow=nrow(C),ncol=ncol(C)) 
  Cquic<-QUIC(C,rho=RHO,tol=1e-02,maxIter=100)
  InvCov<-Cquic$X
  ADJ<-abs(InvCov)
  GRAO<-graph.adjacency(ADJ,mode=c("max"),weighted=TRUE,diag=FALSE)
  g<-induced.subgraph(GRAO,which(degree(GRAO)>0))
  print(length(E(g)))
  return(length(E(g)))
}
######################################################################
#agregar from to como parametros
getEdgesByRHO<-function(RList,start,end, interval, ncores)
{
  
  MList<-vector("list", length(seq(start, end, interval)))
  for(rho in seq(start,
                 end, interval))
  {
    res<-unlist(mclapply(RList,  getEdges,  rho, mc.cores = cores))
    MList[[i]] <- res
    i <- i + 1
  }
  return(MList)
}

######################################################################
getADJ<- function(C, rho)
{
  RHO<-matrix(data=rho,nrow=nrow(C),ncol=ncol(C)) 
  Cquic<-QUIC(C,rho=RHO,tol=1e-02,maxIter=100)
  InvCov<-Cquic$X
  ADJ<-abs(InvCov)
  return(ADJ)
}
#########################################################################
getADJRHO<-function(RList, start,end, interval, ncores)
{
  MList<-vector("list", length(seq(start, end , interval)))
  i=1
  for(rho in seq(start, end , interval))
  {
    res<-mclapply(RList,  getADJ,  rho, mc.cores = ncores)
    MList[[i]] <- res
    i <- i + 1
  }
  return(MList)
}
#######################################################################
getPositiveValues<-function(ADJList, ncores){
  totalPos<-vector()
  for (i in 1:length(ADJList))  {
    ll<-mclapply(ADJList[[i]],
                 function(x)
                 {
                   PosValues<-length(which( unlist(x)>0 ) )
                   return(PosValues)
                 }, mc.cores = ncores
    )   
    totalPos<-c(totalPos, ll)
  }
  return(totalPos)
}
######################################################################
#tengo una lista de listas
getEdgesBySample<-function(rrlist, rho)
{
  RRList<-vector("list")
  i=1
  for (i in 1:length(rrlist))
  {
  res<-mean(unlist(lapply(rrlist[[i]],  getEdges,  rho)))
  RRList[[i]] <- res
  i <- i + 1
  }
  return(RRList)
}
#############################################################
getEdgesByRHO<-function(RList,start,end, interval, ncores)
{
  MList<-vector("list", length(seq(start, end, interval)))
  i=1
  for(rho in seq(start,
                 end, interval))
  {
    res<-unlist(mclapply(RList,  getEdges,  rho, mc.cores = ncores))
    MList[[i]] <- res
    i <- i + 1
  }
  return(MList)
}

#########################Plots###########################################

plots<-function(
  start,
  end,
  interval,
  trueEdges,
  trueMEAN,
  randomEdges,
  randomMEAN,
  name  )
{


ratioNumberOfEdges<-randomEdges/ trueEdges*100
ratioPositiveValues<-randomMEAN/trueMEAN*100
eje1<-seq(start,end,interval)
name<-name
###########################################################################################
file1<-paste(
       paste(
       paste("NumberOfEdges_", paste(start,end,sep="-"),sep="") 
            , name, sep="_"),
              ".png", sep="")

  png(file1)
  plot(eje1,
       trueEdges, 
       main="Number of edges real / random data vs rho" ,
       type="l", 
       xlab="rho", 
       ylab="Number of edges",
       col="blue",
       ylim=c(min(trueEdges, randomEdges), max(trueEdges, randomEdges)),
       xlim=c(start, end))
      lines(eje1, randomEdges, col="red")
      legend("topright", 
         c("real data","random data"), lty = 1, col = c("blue","red"))
  dev.off()
message("Plot 1: Number of edges real / random data vs rho")
###########################################################################################
file2<-paste(
        paste(
        paste("RatioNumberOfEdges_",
        paste(start,end,sep="-"), sep=""), 
        name, sep="_"),
              ".png", sep="")

  png(file2)
  plot(eje1, (randomEdges/trueEdges)*100, 
       main="Ratio edges in random data vs real data",
       type="l")
  abline(h=c(0,5,10))
  dev.off()
message("Plot 2: Ratio edges in random data vs real data")
###########################################################
file3<-paste(paste(paste("RealPositiveValues_",paste(start,end,sep="-"), sep=""), name, sep="_"),".png", sep="")

png(file3)
plot(eje1,
     trueMEAN, 
     type="l", col="blue", 
     main="Adj positive values and edges in real data", 
     xlab="values", ylab="rho",
     ylim=c(min(trueMEAN, trueEdges), max(trueMEAN, trueEdges)),
     xlim=c(start, end))
lines(eje1,trueEdges, col="red")
legend("topright", 
       c("edges","positives values"),
       lty = 1, col = c("red","blue"))
dev.off()
message("Plot 3: Adj positive values and edges in real data")
##############################################
file4<-paste(paste(paste("RandomPositiveValues_",paste(start,end,sep="-"), sep=""), name, sep="_"),".png", sep="")
png(file4)
plot(eje1,randomMEAN, 
     type="l", col="blue", 
     main="Adj positive values and edges in random data",
     ylim=c(min(randomMEAN, randomEdges), max(randomMEAN, randomEdges)),
     xlim=c(start, end)
)
lines(eje1,randomEdges, col="red")
legend("topright", c("edges","positives values"), lty = 1, col = c("red","blue"))
dev.off()
message("Plot 4: Adj positive values and edges in random data")
##############################################
file5<-paste(paste(paste("Edges_vs_Positive_",paste(start,end,sep="-"), sep=""), name, sep="_"),".png", sep="")
png(file5)
  plot(eje1,trueMEAN, 
       type="l", col="blue", 
       main="Adj Matrix positive values and edges in real and random data",
       xlab="rho", 
       ylab="values" ,
       ylim=c(min(randomMEAN,trueMEAN, randomMEAN,randomEdges), 
              max(randomMEAN,trueMEAN, randomMEAN,randomEdges)),
       xlim=c(start, end))
      lines(eje1,trueEdges, col="green")#t
    lines(eje1,randomMEAN, type="l", col="red")#
    lines(eje1,randomEdges, col="black")#
    legend("topright", c("positives values real data",
                       "edges real values",
                       "positives values random data", 
                       "edges random data"), 
         lty = 1, 
         col = c("blue", "green","red", "black"))
dev.off()
message("Plot 5: Adj Matrix positive values and edges in real and random data")
#################################################  
#ratios:

file6<-paste(paste(paste("Ratios_Edges_vs_Positives_",paste(start,end,sep="-"), sep=""), name, sep="_"),".png", sep="")

png(file6)
plot(eje1,  
     ratioNumberOfEdges, 
     type="l", 
     main="Ratios random/real", 
     col="red", 
     xlab="rho", 
     ylab="ratio",
     ylim=c(min(ratioNumberOfEdges, ratioPositiveValues), 
            max(ratioNumberOfEdges, ratioPositiveValues)),
     xlim=c(start, end))
lines(eje1, ratioPositiveValues, col="blue", type="l") 
legend("topright", 
       c("ratio number of edges",
         "ratio positive values"), 
       lty = 1, 
       col = c("red", "blue"))
dev.off()
message("Plot 6: Ratios random/real")
}


##################################################
estimateFDR<-function(
  inputM=sampleData,
  NumRandomM=10,
  start=0.3,
  end=1,
  interval=0.05,
  ncores=1,
  name)
{
  TrueList<-vector("list",1)
  TrueList[[1]]<-sampleData
  TrueList<-rep(TrueList, NumRandomM)
  TPList<-lapply(TrueList, CRobCor)# 
  
  RandomList<-createRandomMatrix(TrueM = sampleData,  NumRandomM)
  RList <- vector("list", length(RandomList))
  RList<-lapply(RandomList, CRobCor)# paraleliza OK
  TPList<-rep(TPList, NumRandomM)
  RandomEdgesList<-getEdgesByRHO(RList, start, end, interval, ncores)
  print(paste("Lenght of random edge list:",length(RandomEdgesList))) 
  TrueEdgesList<-getEdgesByRHO(TPList, start, end, interval, ncores)
  print(paste("Lenght of true edge list:", length(TrueEdgesList)   )  ) 
  trueEdges<-unlist(lapply(TrueEdgesList, mean))
  randomEdges<- unlist(lapply(RandomEdgesList, mean))
  ##################################################
  #Using adjacency matrix
  TrueAdjMatrix<-getADJRHO(TPList, start, end , interval, ncores)
  print(paste("Lenght of true adj matrix list", length(TrueAdjMatrix) )  ) 
  RandomAdjMatrix<-getADJRHO(RList, start, end , interval, ncores)
  message(paste("Lenght of random adj matrix list",  length(RandomAdjMatrix)  )  ) 
  ##################################################
  #from here divide by 2
  tp<-unlist(getPositiveValues(TrueAdjMatrix, ncores))/2
  trueValues<-matrix(nrow=length(seq(start, end, interval)), data=tp, byrow = T)
  trueMEAN<-rowMeans(trueValues)
  rr<-unlist(getPositiveValues(RandomAdjMatrix, ncores))/2
  randomValues<-matrix(nrow=length(seq(start, end, interval)), data=rr, byrow = T)
  randomMEAN<-rowMeans(randomValues)
  ##################################################
  #add an small correcction in order to avoid infinite values
  message("Finish computation. Lets plot.")
  
print( start)
print( end)
print(interval)
print(trueEdges+0.01)
print(trueMEAN+0.01)
print(randomEdges+0.01)
print(randomMEAN+0.01)
print(name)#all together

 plots( start,  end, interval, 
         trueEdges=trueEdges+0.01,
         trueMEAN=trueMEAN+0.01,
         randomEdges=randomEdges+0.01,
         randomMEAN=randomMEAN+0.01,
	 name)#all together
  ##################################################
  
}


