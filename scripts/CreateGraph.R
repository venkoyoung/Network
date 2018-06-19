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
  #size<-rescale(betweenness(g),  c(min(betweenness(g)), max(betweenness(g)) ))
  size<-rescale(betweenness(g),   c(5, 30) )
  
  ##########
  gList<-list(g, memb, Centr_G)
  return(gList)
}


###############################################

