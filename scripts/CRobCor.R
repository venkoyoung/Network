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



