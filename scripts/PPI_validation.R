#ppi validation
library(igraph)
setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/NumOfM100/")
dir()
setwd("../interactions/")
##################################################################
edgeListLC<-read.table("edgelist_048_labchip.tab");head(edgeListLC); dim(edgeListLC)
labchipG<-graph_from_data_frame(edgeListLC, directed = F)
length(E(labchipG))#547
length(V(labchipG))#196
##################################################################
edgeListLR<-read.table("edgelist_055_rnaseq.tab");head(edgeListLR); dim(edgeListLR)
networkR<-graph_from_data_frame(edgeListLR, directed = F)
length(E(networkR))#389
length(V(networkR))#197
###################################################################
# Make the plot
l <- layout_with_fr(g)
l<-layout_randomly(g)
l <- layout_in_circle(g)
l <- layout_on_sphere(g)
l <- layout_with_fr(g)
plot(g, layout=l)
plot(g,ceb, vertex.color="black",
     vertex.frame.color="white",
     vertex.label.cex=0.5,
     vertex.size=1,
     vertex.label.color="black",
     main=" RNA-Seq", 
     edge.arrow.size=.5,
          layout=l)
ceb <- cluster_edge_betweenness(g)
##################################################################
PPi<-read.table("ppi_all_names.txt", sep="\t", header=T)
PPi_DF<-PPi[,7:8]; dim(PPi_DF)
ppiG<-graph_from_data_frame(PPi_DF, directed = F)
length(E(ppiG))#length
length(V(ppiG))#7798
##################################################################
setwd("../interactions/")
enzyme<-read.table("../interactions/enzyme_substrate.tsv", sep="\t", header=T); head(enzyme)
enzymeDF<-data.frame(enzyme$enzyme_genesymbol, enzyme$substrate_genesymbol)
enzymeN<-graph_from_data_frame(enzymeDF, directed = F)
length(E(enzymeN))#710
length(V(enzymeN))#204
##################################################################
int<-read.table("interactions.tsv", sep="\t", header=T); head(int)
intDF<-data.frame(int$source_genesymbol, int$target_genesymbol)
intN<-graph_from_data_frame(intDF, directed = F)
length(E(intN))#1930
length(V(intN))#1127
intNEdges<-get.edgelist(intN); 
write.table(intNEdges, file="edgelist_intN", col.names = NA, quote = F, sep="\t")
#################################################################
resuting<-intersection(labchipG, ppiG,keep.all.vertices = F)
length(E(resuting))
length(V(resuting))
dfRes<-get.edgelist(resuting)
write.table(dfRes, file="edgelist_resulting_labchip.tab", col.names = NA, quote = F, sep="\t")
#################################################################
resuting2<-intersection(networkR, ppiG,keep.all.vertices = F)
length(E(resuting2))#49
length(V(resuting2))#159
dfRes2<-get.edgelist(resuting2)
write.table(dfRes2, file="edgelist_resulting_rnaseq.tab", col.names = NA, quote = F, sep="\t")
#################################################################
labVsR<-intersection(labchipG, networkR,keep.all.vertices = F)
length(E(labVsR))#64
length(V(labVsR))#140
LVR<-get.edgelist(labVsR); 
write.table(LVR, file="edgelist_resulting_labVRNASeq.tab", col.names = NA, quote = F, sep="\t")
#################################################################
validation<-intersection(labVsR, ppiG,keep.all.vertices = F)
length(E(validation))#29
length(V(validation))#117
valida<-get.edgelist(validation); 
write.table(valida, file="edgelist_resulting_validation.tab", col.names = NA, quote = F, sep="\t")
#################################################################
#other sources of validation:
labVsEnzyme<-intersection(labchipG, enzymeN,keep.all.vertices = F)
length(E(labVsEnzyme))#0
length(V(labVsEnzyme))#57
#################################################################
rnaseqVsEnzyme<-intersection(networkR, enzymeN,keep.all.vertices = F)
length(E(rnaseqVsEnzyme))#0
length(V(rnaseqVsEnzyme))#58
############################
labVsInt<-intersection(labchipG, intN,keep.all.vertices = F)
length(E(labVsInt))#9
length(V(labVsInt))#135
#################################################################
rnaseqVsInt<-intersection(networkR, intN,keep.all.vertices = F)
length(E(rnaseqVsInt))#5
length(V(rnaseqVsInt))#140
############################
validationE<-intersection(labVsR, enzymeN,keep.all.vertices = F)
length(E(validationE))#0
length(V(validationE))#39
############################
validationInt<-intersection(labVsR, intN,keep.all.vertices = F)
length(E(validationInt))#4
length(V(validationInt))#102
validaI<-get.edgelist(validationInt); 
write.table(validaI, file="edgelist_resulting_validationInt.tab", col.names = NA, quote = F, sep="\t")
############################
rbpN<-graph_from_data_frame(dfNet, directed = F)
length(unique(dfNet$dfBenInt.sourceVT))#293
length(E(rbpN))#11553360
length(V(rbpN))#12943
rbpsEdges<-get.edgelist(rbpN); 
write.table(rbpsEdges, file="edgelist_rbps.txt", col.names = NA, quote = F, sep="\t")
############################
labVsRBP<-intersection(labchipG, rbpN,keep.all.vertices = F)
length(E(labVsRBP))#105
length(V(labVsRBP))#1 94
validaLabchipRBP<-get.edgelist(labVsRBP); 
write.table(validaLabchipRBP, file="edgelist_resulting_validationLabChip_RBPs.tab", col.names = NA, quote = F, sep="\t")
#################################################################
rnaseqVsRBP<-intersection(networkR, rbpN,keep.all.vertices = F)
length(E(rnaseqVsRBP))#97
length(V(rnaseqVsRBP))#196
validaRBP<-get.edgelist(rnaseqVsRBP); 
write.table(validaRBP, file="edgelist_resulting_validationRNASEQ_RBPs.tab", col.names = NA, quote = F, sep="\t")
############################
validationRBP<-intersection(labVsR, rbpN,keep.all.vertices = F)
length(E(validationRBP))#18
length(V(validationRBP))#140
validaBoth<-get.edgelist(validationRBP); 
write.table(validaBoth, file="edgelist_resulting_validationRNASEQ_Labchip_RBPs.tab", col.names = NA, quote = F, sep="\t")
############################
#ICLIP:
iclip<-graph_from_data_frame(iclipData, directed = F)
length(E(iclip))#724121
length(V(iclip))#11439
clipEdges<-get.edgelist(iclip); 
write.table(validaLabchipClip, file="edgelist_iclip", col.names = NA, quote = F, sep="\t")

labVsCLIP<-intersection(labchipG, iclip,keep.all.vertices = F)
length(E(labVsCLIP))#82
length(V(labVsCLIP))#191
validaLabchipClip<-get.edgelist(labVsCLIP); 
write.table(validaLabchipClip, file="edgelist_resulting_validationLabChip_iclip.tab", col.names = NA, quote = F, sep="\t")

#################################################################
rnaseqVsCLIP<-intersection(networkR,  iclip,keep.all.vertices = F)¡
length(E(rnaseqVsCLIP))#50
length(V(rnaseqVsCLIP))#192
validaRNASEQCLIP<-get.edgelist(rnaseqVsCLIP); 
write.table(validaRNASEQCLIP, file="edgelist_resulting_validationRNASEQ_CLIP.tab", col.names = NA, quote = F, sep="\t")
############################
validationBothClip<-intersection(labVsR, iclip,keep.all.vertices = F)
length(E(validationBothClip))#19
length(V(validationBothClip))#138
validaBothCLIP<-get.edgelist(validationBothClip); 
write.table(validaBothCLIP, file="edgelist_resulting_validation_bothCLIP.tab", col.names = NA, quote = F, sep="\t")
############################
############comparison iclip int
RBPvsICLIP<- intersection(rbpN, iclip,keep.all.vertices = F)
V(rbpN)
length(E(RBPvsICLIP))#44418
length(V(RBPvsICLIP))#9057

RBPvsPPI<- intersection(rbpN, ppiG,keep.all.vertices = F)
length(E(RBPvsPPI))#1084
length(V(RBPvsPPI))#5936

PPIvsICLIP<- intersection(ppiG, iclip,keep.all.vertices = F)
length(E(PPIvsICLIP))#553
length(V(PPIvsICLIP))#5924
############################
#only exons vs eventos del network
exons<-g
exonsvsICLIP<- intersection(exons, iclip,keep.all.vertices = F)
length(E(exonsvsICLIP))#50 #75 #111
length(V(exonsvsICLIP))#142 #198 #231
V(iclip)=="CW22"
###########################
exonsvsppI<- intersection(exons, ppiG,keep.all.vertices = F)
length(E(exonsvsppI))#41 #49 #59
length(V(exonsvsppI))#126 #165 #193
############################
exonsvRnaSeq<- intersection(exons, networkR,keep.all.vertices = F)
length(E(exonsvRnaSeq))#30 #33 #39
length(V(exonsvRnaSeq))#115 #156 #179
############################
exonsVsLabchip<-intersection(exons, labchipG,keep.all.vertices = F)
length(E(exonsVsLabchip))#65
length(V(exonsVsLabchip))#173
############################
#0.6, 0.55, 0.5

