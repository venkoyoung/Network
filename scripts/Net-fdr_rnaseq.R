setwd("../")
source("../scripts/functions_fdr_noMC.R")
t_Zvals<-eventsscaled
sampleData <- t_Zvals
#first assign class ofdata for a small number of rows
NumRandomM=10
#param: t_Zvals, number of matrix to generarate, number of nodes to use.
TrueList<-vector("list",1)
TrueList[[1]]<-t_Zvals
TrueList<-rep(TrueList, NumRandomM)
length(TrueList)
TPList<-lapply(TrueList, CRobCor)# paraleliza OK
length(TPList)
RandomList<-createRandomMatrix(TrueM = t_Zvals,  NumRandomM)
RList <- vector("list", length(RandomList))
RList<-lapply(RandomList, CRobCor)# paraleliza OK
TPList<-rep(TPList, NumRandomM)
##################################################
#ya tengo las 2 listas
#paso a extraer los valores
#1.- Number of nodes
RandomEdgesList<-getEdgesByRHO(RList, start=0.3, end=0.9, interval=0.02);length(RandomEdgesList)#21
TrueEdgesList<-getEdgesByRHO(TPList, start=0.3, end=0.9, interval=0.02); length(TrueEdgesList)#21#OK
#not using mclapply:
lapply(RandomEdgesList, mean)
lapply(TrueEdgesList, mean)
#PLOTS
eje1<-seq(0.3,0.9,0.02)
trueEdges<-unlist(lapply(TrueEdgesList, mean))
randomEdges<- unlist(lapply(RandomEdgesList, mean))

png("Number_of_edges0309.png")
plot(eje1,trueEdges, main="Number of edges real / random data vs rho" ,
     type="l", xlab="rho",ylab="number of edges", col="blue", ylim=c(0, max(trueEdges)))
lines(eje1, randomEdges, col="red")
legend(.6,2000 , c("real data","random data"), lty = 1, col = c("blue","red"))
dev.off()
##################################################################
png("ratio_number_of_edges0509.png")
plot(eje1, (randomEdges/trueEdges)*100, main="Ratio edges in random data vs real data",type="l")
abline(h=c(0,5,10), v=c(0,0.5))
dev.off()
##################################################
#Using adjacency matrix 
TrueAdjMatrix<-getADJRHO(TPList, start=0.3, end=0.9 , interval=0.02); length(TrueAdjMatrix)#10
RandomAdjMatrix<-getADJRHO(RList, start=0.3, end=0.9 , interval=0.02); length(RandomAdjMatrix)#10
##################################################
#from here divide by 2
tp<-unlist(getPositiveValues(TrueAdjMatrix))/2
#trueValues<-matrix(nrow=length(seq(0.3, 0.7 , .02)), data=tp, byrow = T)
trueValues<-matrix(nrow=length(seq(0.3, 0.9 , 0.02)), data=tp, byrow = T)
trueMEAN<-rowMeans(trueValues)

rr<-unlist(getPositiveValues(RandomAdjMatrix))/2
#randomValues<-matrix(nrow=length(seq(0.3, 0.7 , .02)), data=rr, byrow = T)
randomValues<-matrix(nrow=length(seq(0.3, 0.9 , 0.02)), data=rr, byrow = T)
randomMEAN<-rowMeans(randomValues)
#number of positives
png("real_data_positiveValues0509.png")
plot(eje1,trueMEAN, type="l", col="blue", main="Adj positive values  and edges in real data", 
     ylab="values", xlab="rho", ylim=c(0,4000))
lines(eje1,trueEdges, col="red")
legend(.6,1000 , c("edges","positives values"), lty = 1, col = c("red","blue"))
dev.off()
##################################################
png("random_data_positiveValues0509.png")
plot(eje1,randomMEAN, type="l", col="blue", main="Adj positive values and edges in random data", ylim=c(0,500))
lines(eje1,randomEdges, col="red")
legend(.6,300 , c("edges","positives values"), lty = 1, col = c("red","blue"))
dev.off()
##################################################
#alltogether
png("numberOfEdges_positiveValues0309.png")
plot(eje1,trueMEAN, type="l", col="blue", 
#     main="Adj Matrix positive values and edges in real and random data",
     xlab="rho", ylab="values",
     ylim=c(0,max(trueMEAN)) )
lines(eje1,trueEdges, col="green")#t
lines(eje1,randomMEAN, type="l", col="red")#
lines(eje1,randomEdges, col="black")#
legend(.55,2000 , c("positives values real data",
                   "edges real values",
                   "positives values random data", 
                   "edges random data"), 
       lty = 1, 
       col = c("blue", "green","red", "black"))
dev.off()
#ratios:  
ratioNumberOfEdges<-randomEdges/ trueEdges*100
ratioPositiveValues<-randomMEAN/trueMEAN*100
plot(eje1, randomMEAN, type="l", col="red")
lines(eje1, trueMEAN, col="blue")
plot(eje1, (randomMEAN/ trueMEAN)*100, type="l")
plot(eje1, (randomEdges/ trueEdges)*100)
abline(h=5)
png("ratios_numberOfEdges_positiveValues.png")
plot(eje1,  ratioNumberOfEdges, type="l", main="Ratios random/real", col="red", xlab="rho", ylab="ratio", ylim=c(0,100))
lines(eje1, ratioPositiveValues , col="blue", type="l")
legend(.6,300, legend=c("ratio number of edges","ratio positive values"),
       lty = 1,  col = c("red", "blue"))
       
dev.off()
#############compare matrix by matrix for each RHO
comparison<-PositivesByCell(TrueAdjMatrix, RandomAdjMatrix, cores=6)
comparison<-PositivesByCell(TrueAdjMatrix, RandomAdjMatrix)
tpByho<-matrix(nrow=length(seq(0.3, 0.7 , .02)),
               data=as.numeric(unlist(comparison)), byrow = T)
falsePositives<-rowMeans(tpByho)
length(unlist(TrueAdjMatrix[[1]]))/2
################################################################
png("cell_by_cell.png")
plot(eje1, falsePositives, type="l", col="red")
dev.off()
################################################################
plot(rep(eje1, each=10), as.numeric(unlist(comparison)))
plot(eje1, falsePositives/35644*100, type="l", col="red")
abline(h=1)
################################################################
###subsampling rows from total to total/2
##3rho fijo de 0.5
#genero una lista de 10 matrices iguales reales
library(parallel)
#first generate data: TRUE and RANDOM DATA
sampleData <- read.csv("t_Zvals.tab", header = TRUE, nrows = 5, sep="\t")
#first assign class ofdata for a small number of rows
classes <- sapply(sampleData, class); table(classes)
#now read the full data
largeData <- data.matrix(read.csv("t_Zvals.tab", header = TRUE, colClasses = classes, sep="\t"))
TrueList<-rep(TrueList, NumRandomM); length(TrueList)
truevalList<-TrueList
#genero una lista de 10 matrices iguales a partir de una matrix aleatoria
randomval<-list()
randomval[[1]]<-RandomList[[1]];length(randomval)

#############################################################
#subsampleo cada matrix desde 18 hasta 35 filas, el numero de matrices que tenga la lista. 
#Sobre cada matrix estimo correlacion
source("scripts/functions_fdr.R")
TrueList<-vector("list",1)
TrueList[[1]]<-largeData
NumRandomM=10000
TrueList<-rep(TrueList, NumRandomM)
truevalList<-TrueList
#genero una lista de 10 matrices iguales a partir de una matrix aleatoria
RandomList<-createRandomMatrix(TrueM = largeData,  NumRandomM)
randomval<-list()
randomval[[1]]<-RandomList[[1]];length(randomval)
randomvalList<-rep(randomval, NumRandomM)

trlist<-randomRowsList(truevalList, cores=8); length(truevalList)# ya calculó la correlación
rrlist<-randomRowsList(randomvalList, cores=8); length(randomvalList)
#randomval<-randomRowsList(RandomList)
#obtengo el numero de edges promedio para cada grafo estimado a partir de cada matrix de correlacion
edgesByRRows<-getEdgesBySample(rrlist, rho=.54, cores=8)#OK   
edgesByTRows<-getEdgesBySample(trlist, rho=.54, cores=8)
#############################################################
plot(5:35,unlist(edgesByRRows), type="l")
plot(5:35,unlist(edgesByTRows), type="l", col="blue",    xlim=range(5:35), 
     ylim=range(0:3000) )
lines(5:35,unlist(edgesByRRows), col="red")
legend(20,2000 , c("real data",
                   "random data"), 
       lty = 1, 
       col = c("blue", "red"))
#############################################################
  plot(5:35,unlist(edgesByRRows)/unlist(edgesByTRows),
       type="l", col="blue",    xlim=range(5:35))
#para esto hay que hacer un loop para obtener el rho optimo para cada numero de eventos
length(rrlist)
rrlist[1:5]
edgesByRRowsbyRho<-getEdgesBySampleByRho(rrlist[1:5], 
                                         start=0.3, 
                                         end=0.7, 
                                         interval=0.02,
                                         cores=6)#OK   
edgesByTRowsbyRho<-getEdgesBySampleByRho(trlist[1:5],
                                         start=0.3,
                                         end=0.7,
                                         interval=0.02,
                                         cores=6)
##################################################################################
plot(lines(seq(0.3, 0.6,.02),SampleList[[1]]), type="l",    xlim=range(0.3, 0.6), 
     ylim=range(0:3000) )
colores<-rev(heat.colors(length(SampleList), alpha = 1))
##################################################################################
for (l in 1:length(SampleList))
{
  lines(seq(0.3, 0.6,.02),SampleList[[l]], col=colores[l])
}
legend("topright",
       legend=c(5:35), 
       lty = 1, 
       col = colores)
##################################################################################
