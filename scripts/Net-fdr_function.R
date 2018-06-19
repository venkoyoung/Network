 fdrPlots<- function (matrix, NR, start, end, interval)
 {
 #############################################
 #functions
   ###### Scale vector in 0-1 range
   Vscale <- function(v, a, b) {
     v <- v-min(v) ; v <- v/max(v) ; v <- v * (b-a) ; v+a
   }
  #############################################
   CRobCor <- function (M,ct=0.5) {
     #rows: eventos
     #columns: factores
     Mcov<-cov((M),use="all.obs")
     #evalua con Pearson
     Mcor<-cov2cor(Mcov)
     #Scaling a covariance matrix into a correlation one can be achieved in many ways,
     #mathematically most appealing by multiplication with a diagonal matrix from left and right, 
     #or more efficiently by using sweep(.., FUN = "/") twice. 
     #The cov2cor function is even a bit more efficient, and provided mostly for didactical reasons.  
     DELTA<-array(data=0,
                  dim=c(ncol(M),
                        ncol(M),
                        nrow(M)),
                  dimnames=list(colnames(M),
                                colnames(M),
                                rownames(M)))
     # p x p x n dimensional array (obviously ofr large numbers of n this becomes impossible to compute)
     MRobCor<-Mcor;
     for (ev in 1:nrow(M)) {
       marg_M<-M[-c(ev),];
       marg_M<-scale(marg_M) 
       marg_Mcor<-cov(marg_M,use="all.obs")
       marg_Mcor<-cov2cor(marg_Mcor)
       DELTA[,,ev]<-abs(marg_Mcor-Mcor)		#marginal Delta Covariance
     }
     
     INFL<-apply(DELTA,3,sum)/(ncol(M)^2-ncol(M)); #Calculate "influence" of outliers
     
     for (p_i in 1:(ncol(M)-1)){
       for (p_j in (p_i+1):ncol(M)){
         if (abs(Mcor[p_i,p_j])<0.25) next               
         cta=ct-0.15 /( 0.5*length(which(abs(Mcor[p_i,])>ct))+ 0.5*length(which(abs(Mcor[p_j,])>ct)) )
         WT=apply(DELTA[p_i,  c(p_j,which(abs(Mcor[p_i,])>cta)),],  2,sum)/(length(which(abs(Mcor[p_i,])>cta)) )    #Calculate Weights
         WT=WT+apply(DELTA[p_j,  c(p_i,which(abs(Mcor[p_j,])>cta)),],  2,sum)/(length(which(abs(Mcor[p_j,])>cta)) )    #Calculate Weights 
         CRij=cov.wt(M[,c(p_i,p_j)], ((1/WT^(0.8))),cor=TRUE)$cor  #Weighted correlation
         MRobCor[p_i,p_j]=CRij[1,2]
         MRobCor[p_j,p_i]=CRij[2,1]
       }
     }
     
     return(MRobCor)
   }
   ###############################################################
   
   # CENTRALITY/AUTHORITY MEASURES. Return a matrix of authority scores and corresponding node rankings
   CentralityRanking<-function(g){
     DG<-degree(g)
     PR<-page.rank(g)$vector  		#Pagerank Score
     OPR<-order(PR,decreasing=TRUE)
     
     BC<-betweenness(g)		#Betweeness Centrality
     OBC<-order(BC,decreasing=TRUE)	
     
     CC<-closeness(g)		#Closeness Centrality
     OCC<-order(CC,decreasing=TRUE)	
     
     RANK<-rank(PR)+rank(BC)+rank(CC)
     ORANK<-order(RANK,decreasing=TRUE)	#Ordering of Nodes based on aggregate score.
     return(rbind(DG,PR,OPR,BC,OBC,CC,OCC,RANK,ORANK))
   }
##############################################################
   ##### Network Reconstruction, Plotting and Analysis. Input is a covariance matrix C and a regularization paramether rho
   CreateGraph<-function(C,rho){
     RHO<-matrix(data=rho,nrow=nrow(C),ncol=ncol(C)) #Create reqularization matrix, es la matriz para regularizar 
     #InvCov<-glasso(s=C,rho=RHO,maxit=20000,penalize.diagonal=TRUE)$wi  #Graphical Lasso estimation of sparse inverse covariance matrix
     Cquic<-QUIC(C,rho=RHO,tol=1e-02,maxIter=100)
     rownames(Cquic$X)<-rownames(C)
     colnames(Cquic$X)<-colnames(C)
     #InvCov<-QUIC(C,rho=RHO,tol=1e-02,maxIter=100,msg=0)$X  #0.1 for C4/ 0.3 for C3 / 0.5 for C2
     InvCov<-Cquic$X
     ADJ<-abs(InvCov) #Weighted adjacency matrix (taking absolute values)
     TADJ<-InvCov     #"True" Weighted Adjacency matrix (keeping the sign)
     
     GRAO<-graph.adjacency(ADJ,mode=c("max"),weighted=TRUE,diag=FALSE)
     TGRAO<-graph.adjacency(TADJ,mode=c("max"),weighted=TRUE,diag=FALSE)
     labels=rownames(C)
     GRAO<-set.vertex.attribute(GRAO,"labels",value=labels)
     ### Remove isolated nodes:
     g<-induced.subgraph(GRAO,which(degree(GRAO)>0)) ; head(g)
     Tg<-induced.subgraph(TGRAO,which(degree(GRAO)>0)); head (Tg) 
     ### Calculate Centrality scores and node rankings
     Centr_G<-CentralityRanking(GRAO)
     ######## COMMUNITY DETECTION #########
     fg_C<-fastgreedy.community(g, modularity=TRUE) #infomap now available in igraph is better...
     im_C<-infomap.community(g, modularity=TRUE) #infomap now available in igraph is better...
     memb=fg_C
     names(memb$membership)<-get.vertex.attribute(g,"labels")
     V(g)$memb<-memb$membership
     csize=table(memb$membership)
     ##################################### Plotting Parameters for the Graph Object
     comps <-memb$membership; head(comps)#dpmde esta cada gen, em que grupo
     colbar <- c(brewer.pal(8,"Dark2")[2:3],
                 brewer.pal(8,"Set1"),
                 brewer.pal(8,"Accent")[1],
                 brewer.pal(8,"Pastel2"),
                 brewer.pal(9,"Pastel1"),
                 brewer.pal(12,"Paired"),
                 brewer.pal(8,"Dark2"),
                 brewer.pal(8,"Set2"))
     V(g)$color <- colbar[comps]
     ##### ADD TRANSPARENCIES #####
     E(g)$width<-Vscale(E(g)$weight,1.5,9)
     V(g)$color <- colbar[comps]# color according int community
     V(g)$color<-paste(V(g)$color,"44",sep="")
     E(g)$color<-ifelse(E(Tg)$weight<0,"#55885599","#BB111199")
     V(g)$label_color="#222222"
     l<-layout.fruchterman.reingold(g)
     size<-rescale(betweenness(g),   c(5, 30) )
     
     ##########
     gList<-list(g, memb, Centr_G)
     return(gList)
   }
###############################################
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
                           RandomValuesMatrixList)
 {
   if(length(TrueValuesMatrixList) == length(RandomValuesMatrixList))
   {
     final<-vector()
     for (i in 1:length(RandomValuesMatrixList))  {
       ll<-lapply(RandomValuesMatrixList[[i]],
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
#############################################################
#first assign class ofdata for a small number of rows
t_Zvals<-matrix
NumRandomM=NR
#param: t_Zvals, number of matrix to generarate, number of nodes to use.
TrueList<-vector("list",1)
TrueList[[1]]<-t_Zvals
TrueList<-rep(TrueList, NumRandomM)
length(TrueList)
TPList<-lapply(TrueList, CRobCor)# paraleliza OK
RandomList<-createRandomMatrix(TrueM = t_Zvals,  NumRandomM)
RList <- vector("list", length(RandomList))
RList<-lapply(RandomList, CRobCor)# paraleliza OK
TPList<-rep(TPList, NumRandomM)
##################################################
#ya tengo las 2 listas
#paso a extraer los valores
#1.- Number of nodes
start
end
interval

RandomEdgesList<-getEdgesByRHO(RList, start=start, end=end, interval=interval)
TrueEdgesList<-getEdgesByRHO(TPList, start=start, end=end, interval=interval)
##################################################
#PLOTS
eje1<-seq(start,end,interval)
trueEdges<-unlist(lapply(TrueEdgesList, mean))
randomEdges<- unlist(lapply(RandomEdgesList, mean))
file1<-paste( paste(start, end, sep=""),"NumEdges.png", sep="_" )

png(file1)
plot(eje1,trueEdges, main="Number of edges real / random data vs rho" ,
     type="l", xlab="rho",ylab="number of edges", col="blue")
lines(eje1, randomEdges, col="red")
legend(.6,2000 , c("real data","random data"), 
       lty = 1, col = c("blue","red"))
dev.off()
##################################################################
file2<-paste( paste(start, end, sep=""),"NumEdgesRatio.png", sep="_" )
png(file2)
plot(eje1, (randomEdges/trueEdges)*100, main="Ratio edges in random data vs real data",type="l", xlab="rho", ylab="ratio")
abline(h=c(0,5,10), v=c(0,0.5))
dev.off()
##################################################
#Using adjacency matrix 
TrueAdjMatrix<-getADJRHO(TPList, start=start, end=end, interval=interval)
RandomAdjMatrix<-getADJRHO(RList, start=start, end=end, interval=interval)
##################################################
#from here divide by 2
tp<-unlist(getPositiveValues(TrueAdjMatrix))/2
#trueValues<-matrix(nrow=length(seq(0.3, 0.7 , .02)), data=tp, byrow = T)
trueValues<-matrix(nrow=length(seq(start, end, interval)), data=tp, byrow = T)
trueMEAN<-rowMeans(trueValues)
rr<-unlist(getPositiveValues(RandomAdjMatrix))/2
#randomValues<-matrix(nrow=length(seq(0.3, 0.7 , .02)), data=rr, byrow = T)
randomValues<-matrix(nrow=length(seq(start, end, interval)), data=rr, byrow = T)
randomMEAN<-rowMeans(randomValues)
#number of positives
file3<-paste( paste(start, end, sep=""),"TruePosValues.png", sep="_" )
png(file3)
plot(eje1,trueMEAN, type="l", col="blue",
     main="Adj positive values  and edges in real data", 
     ylab="Pos values", xlab="rho", ylim=c(0,1600))
lines(eje1,trueEdges, col="red")
legend(.6,1000 , c("edges","positives values"), lty = 1, col = c("red","blue"))
dev.off()
##################################################
file4<-paste( paste(start, end, sep=""),"RandomPosValues.png", sep="_" )
png(file4)
plot(eje1,randomMEAN, type="l", 
     col="blue", 
     main="Adj positive values and edges in random data")
lines(eje1,randomEdges, col="red")
legend(.6,300 , c("edges","positives values"), lty = 1, col = c("red","blue"))
dev.off()
##################################################
#alltogether
file5<-paste( paste(start, end, sep=""),"Allvalues.png", sep="_" )
png(file5)
plot(eje1,trueMEAN, type="l", col="blue", 
     main="Adj Matrix positive values and edges in real and random data",
     xlab="rho", ylab="values" )
lines(eje1,trueEdges, col="green")#t
lines(eje1,randomMEAN, type="l", col="red")#
lines(eje1,randomEdges, col="black")#
legend(.65,1000 , c("positives values real data",
                   "edges real values",
                   "positives values random data", 
                   "edges random data"), 
       lty = 1, 
       col = c("blue", "green","red", "black"))
dev.off()
#ratios:
ratioNumberOfEdges<-randomEdges/ trueEdges*100
ratioPositiveValues<-randomMEAN/trueMEAN*100

file6<-paste( paste(start, end, sep=""),"AllRatios.png", sep="_" )
png(file6)
plot(eje1,  ratioNumberOfEdges, type="l", 
     main="Ratios random/real", col="red", xlab="rho", ylab="ratio")
lines(eje1, ratioPositiveValues , col="blue", type="l") 
legend(.6,20 , c("ratio number of edges",
                 "ratio positive values"), 
       lty = 1, 
       col = c("red", "blue"))
dev.off()
#############compare matrix by matrix for each RHO
comparison<-PositivesByCell(TrueAdjMatrix, RandomAdjMatrix)
tpByho<-matrix(nrow=length(seq(start, end , interval)),
               data=as.numeric(unlist(comparison)), byrow = T)
falsePositives<-rowMeans(tpByho)
################################################################
file7<-paste( paste(start, end, sep=""),"cell.png", sep="_" )
png(file7)
plot(eje1, falsePositives, type="l", col="red")
dev.off()
################################################################
}