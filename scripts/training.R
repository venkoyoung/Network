net1 <- graph_from_literal(D-A:B:F:G, A-C-F-A, B-E-G-B, A-B, F-G,
                           H-F:G, H-I-J)
net1DF<-get.edgelist(net1); 
write.table(net1DF, file="net1DF.txt", col.names = NA, quote = F, sep="\t")

net2 <- graph_from_literal(D-A:F:Y, B-A-X-F-H-Z, F-Y)
net2DF<-get.edgelist(net2); 
write.table(net2DF, file="net2DF.txt", col.names = NA, quote = F, sep="\t")


net1 <- simplify(net1)%>%
set_edge_attr("weight", value = 1:length(E(net1))) %>%
set_edge_attr("color", value = "red")  %>%
set_edge_attr("label", value = LETTERS[1:16])
plot(net1, edge.width = E(net1)$weight)
plot(net1,layout=layout.fruchterman.reingold,
     edge.width=E(net1)$weight)

net3<-graph.union(net1, net2)
g = graph.union(g1, g2, byname=TRUE, keep.x.vertex.attributes=TRUE)

library(igraph)
par(mfrow=c(1,2))
g1 <- barabasi.game(10)
g2 <- barabasi.game(5)
plot(g1)
plot(g2, edge.color='black', vertex.color='green', add=T)
ecount(net1)
ecol(net1)


g <- make_(ring(10), with_vertex_(name = LETTERS[1:10]))
plot(intersection(E(g)[1:6], E(g)[5:9]))
get.edgelist(net2)

gg<-matrix(c(1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,3, 2,3, 3,5), ncol=2)
g <- graph(  c(1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,3, 2,3, 3,5),
  n=5, directed=FALSE )
V(g) #Vertex sequence;
unique(c(1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,3, 2,3, 3,5))
E(g) #Edge sequence
plot(g, layout=layout.circle, vertex.label=c("a","b","c","d","e"))

net1 <- graph_from_literal(D-A:B:F:G, A-C-F-A, B-E-G-B, A-B, F-G,H-F:G, H-I-J)
net2 <- graph_from_literal(D-A:F:Y, B-A-X-F-H-Z, F-Y)
net3<-net1 + net2
E(net3)$curved <- TRUE
E(net3)$width <- count.multiple(net3)
length(E(net3))
plot(net3, layout=layout.circle)

net4<-interesectnet1+net3
plot(net4, layout=layout.circle, edge.curved=TRUE)

plot(net3, layout=layout.circle)
net3

library(igraph)
g <- graph.ring(10, directed=TRUE, mutual=TRUE)
E(g)$color <- sample(1:5, ecount(g), replace=TRUE)
plot(, 
     layout=layout.circle, 
     edge.curved=TRUE, edge.arrow.mode=1)

E(net)$width <- 1.5
plot(net, edge.color=c("dark red", "slategrey")[(E(net)$type=="hyperlink")+1],
     vertex.color="gray40", layout=layout_in_circle, edge.curved=.3)


multigtr <- graph( edges=c(1,2, 1,2, 1,2), n=2 )
l <- layout_with_kk(multigtr)
# Let's just plot the graph:
plot(multigtr, vertex.color="lightsteelblue", vertex.frame.color="white",
     vertex.size=40, vertex.shape="circle", vertex.label=NA,
     edge.color=c("gold", "tomato", "yellowgreen"), edge.width=10,
     edge.arrow.size=3, edge.curved=0.1, layout=l)

plot(multigtr, vertex.color="lightsteelblue", vertex.frame.color="white",
     vertex.size=40, vertex.shape="circle", vertex.label=NA,
     edge.color=c("gold", "tomato", "yellowgreen"), edge.width=10,
     edge.arrow.size=3, edge.curved=curve_multiple(multigtr), layout=l)

net1; E(net1)
net2<-net1
E(net2)

net1
net1 <- add_edges(net1, E(net2))
?add_edges
###############################################################
install.packages("GGally")
install.packages("network")
install.packages("geomnet")
library(GGally)
library(network)
library(geomnet)
library(ggplot2)

data(madmen, package = 'geomnet')
head(madmen)
rownames(madmen$vertices) <- madmen$vertices$label
mm.net %v% "gender" <- as.character(
  madmen$vertices[ network.vertex.names(mm.net), "Gender"])
mm.net <- network(madmen$edges[, 1:2], directed = FALSE)
mm.net
mm.col <- c("female" = "red", "male" = "blue")
set.seed(10052016)

ggnet2(mm.net, 
#       color = mm.col[ mm.net %v% "gender" ],
       labelon = TRUE, 
       label.color = mm.col[ mm.net %v% "gender" ],
       size = 2, 
       vjust = -0.6, 
       mode = "kamadakawai", 
       label.size = 3)

MMnet <- fortify(as.edgedf(madmen$edges), madmen$vertices)
head(MMnet)
set.seed(10052016)
ggplot(data = MMnet, aes(from_id = from_id, to_id = to_id)) +
  geom_net(aes(colour = Gender), layout.alg = "kamadakawai",
           size = 2, labelon = TRUE, vjust = -0.6, ecolour = "grey60",
           directed =FALSE, fontsize = 3, ealpha = 0.5) +
  scale_colour_manual(values = c("#FF69B4", "#0099ff")) +
  xlim(c(-0.05, 1.05)) +
  theme_net() +
  theme(legend.position = "bottom")

library(ggnetwork)
install.packages("ggnetwork")
set.seed(10052016)
########################################
data(email, package = 'geomnet')
# create node attribute data
em.cet <- as.character(email$nodes$CurrentEmploymentType)
names(em.cet) = email$nodes$label
class(em.cet)
em.cet
# remove the emails sent to all employees
edges <- subset(email$edges, nrecipients < 54)
head(edges)
# create network
em.net <- edges[, c("From", "to") ]; head(em.net)
em.net <- network(em.net, directed = TRUE)
em.net
# create employee type node attribute
head(em.cet)#vector
em.net %v% "curr_empl_type" <- em.cet[ network.vertex.names(em.net) ]

set.seed(10312016)
ggnet2(em.net, 
       color = "curr_empl_type",
       size = 4, 
       palette = "Set1", 
       arrow.gap = 0.02,
       arrow.size = 5, 
       edge.alpha = 0.25,
       mode = "fruchtermanreingold",
       edge.color = c("color", "grey50"),
       color.legend = "Employment Type") +
      theme(legend.position = "bottom")
#############################################################
# use em.net created in ggnet2step
#####################try to reproduce
net1E<-as.data.frame(get.edgelist(net1) ); class(net1E)
colnames(net1E)<-c("From","To")
net1E2<-rbind(net1E,net1E)
em.cet<-c(rep(1, 4), rep(2, 5), 3)
names(em.cet)<-unique(union(net1E$From, net1E$To))
em.net12 <- network(net1E2, directed = TRUE)
em.net12 %v% "curr_empl_type" <- em.cet[ network.vertex.names(em.net12) ]
set.seed(10312016)

ggnet2(em.net12, 
       color = "curr_empl_type",
       size = 4, 
       palette = "Set1", 
       arrow.gap = 0.02,
       arrow.size = 5, 
       edge.alpha = 0.25,
       mode = "fruchtermanreingold",
       edge.color = c("color", "grey50"),
       color.legend = "Employment Type") +
  theme(legend.position = "bottom")



email$edges <- email$edges[, c(1,5,2:4,6:9)]
# use em.net created in ggnet2step
em.net
set.seed(10312016)
ggplot(ggnetwork(em.net,
                 arrow.gap = 0.02,
                 layout = "fruchtermanreingold"),
                aes(x, y, xend = xend, yend = yend)) +
  geom_edges(    aes(color = curr_empl_type),
    alpha = 0.25,    arrow = arrow(length = unit(5, "pt"),
                  type = "closed"), curvature = 0.05) +
  geom_nodes(aes(color = curr_empl_type),size = 4) +
  scale_color_brewer("Employment Type",palette = "Set1") +
  theme_blank() +  theme(legend.position = "bottom")



###################################################
## Load R libraries
library(igraph)

# Set adjacency matrix
g <- matrix(c(0,1,1,1, 1,0,1,0, 1,1,0,1, 0,0,1,0),nrow=4,ncol=4,byrow=TRUE)

# Set adjacency matrix to graph object
g <- graph.adjacency(g,mode="directed")
# Add node attribute label and name values
V(g)$name <- c("n1","n2","n3","n4")
# Set subgraph members
c <- c("n1","n2","n3")
# Add edge attribute id values
E(g)$id <- seq(ecount(g))
# Extract supgraph
ccsg <- induced.subgraph(graph=g,vids=c)
# Extract edge attribute id values of subgraph
ccsgId <- E(ccsg)$id
# Set graph and subgraph edge and node colors and sizes
E(g)$color="grey"
E(g)$width=2
E(g)$arrow.size=1
E(g)$arrow.width=1
E(g)[ccsgId]$color <- "#DC143C" # Crimson
E(g)[ccsgId]$width <- 2
V(g)$size <- 4
V(g)$color="#00FFFF" # Cyan
V(g)$label.color="#00FFFFyan"
V(g)$label.cex <-1.5
V(g)[c]$label.color <- "#DC143C" # Crimson
V(g)[c]$color <- "#DC143C" # Crimson

# Set seed value
set.seed(40041)
# Set layout options
l <- layout.fruchterman.reingold(g)
# Plot graph and subgraph
plot.igraph(x=g,layout=l)

#####################################

