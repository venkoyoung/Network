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
  print(dim(C))
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
getEdgesByRHO<-function(RList,start,end, interval)
{
  
  MList<-vector("list", length(seq(start, end, interval)))
  for(rho in seq(start,
                 end, interval))
  {
    res<-unlist(lapply(RList,  getEdges,  rho))
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
getADJRHO<-function(RList, start,end, interval)
{
  MList<-vector("list", length(seq(start, end , interval)))
  i=1
  for(rho in seq(start, end , interval))
  {
    res<-lapply(RList,  getADJ,  rho)
    MList[[i]] <- res
    i <- i + 1
  }
  return(MList)
}
#######################################################################
getPositiveValues<-function(ADJList){
  totalPos<-vector()
  for (i in 1:length(ADJList))  {
    ll<-lapply(ADJList[[i]],
               function(x)
               {
                 PosValues<-length(which( unlist(x)>0 ) )
                 return(PosValues)
               }
    )   
    totalPos<-c(totalPos, ll)
  }
  return(totalPos)
}
######################################################################
PositivesByCell<-function(TrueValuesMatrixList,
                          RandomValuesMatrixList, cores)
{
  if(length(TrueValuesMatrixList) == length(RandomValuesMatrixList))
  {
    final<-vector()
    for (i in 1:length(RandomValuesMatrixList))  {
      ll<-mclapply(RandomValuesMatrixList[[i]],mc.cores = cores,
                 function(x)
                 {
                   
                   falsePos<- length(which(which(unlist(x)  >0  )  %in% which(unlist(TrueValuesMatrixList[[i]]) >0 )     ) )
                   return(falsePos)
                 }
      )   
      final<-c(final, ll)
    }
    return(final)
  }
  else
  {
    print("List are not of the same length")
  }
}
#######################################################################
randomRowsList<-function(mmList)
{
  rrowsList<-vector("list")
  i=1
  for(n in 5:nrow(mmList[[1]]) )
  {
    print(n)
    res<-lapply(mmList,  randomRows,  n)
    print(dim(res[[1]]))
    cres<-lapply(res, CRobCor)
    rrowsList[[i]] <- cres#lista de matrices
    print(length(cres))
    i <- i + 1
  }
  return(rrowsList)
}
#esto me da una lista de listas de matrices con 
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
getEdgesBySampleByRho<-function(rrlist, rho,  start, end , interval)
{
  SampleList<-vector("list", length(rrlist))
  i=1
  for (i in 1:length(rrlist))
  {
    RandomEdgesList<-getEdgesByRHO(rrlist[[i]], start, end, interval)
    meanByRHO<-
      
    tpByho<-matrix(nrow=length(seq(0.3, 0.6,.02)),
                   data=as.numeric(unlist(RandomEdgesList)),
                   byrow = T)
    
    SampleList[[i]] <-  rowMeans(tpByho)
    i <- i + 1
  }
  return(SampleList)
 }
#############################################################
getEdgesByRHOmc<-function(RList,start,end, interval, cores)
{
  MList<-vector("list", length(seq(start, end, interval)))
  i=1
  for(rho in seq(start,
                 end, interval))
  {
    res<-unlist(mclapply(RList,  getEdges,  rho, mc.cores = cores))
    MList[[i]] <- res
    i <- i + 1
  }
  return(MList)
}
#############################################################
getEdgesByRHO<-function(RList,start,end, interval)
{
  MList<-vector("list", length(seq(start, end, interval)))
  i=1
  for(rho in seq(start,
                 end, interval))
  {
    res<-unlist(lapply(RList,  getEdges,  rho))
    MList[[i]] <- res
    i <- i + 1
  }
  return(MList)
}

