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
