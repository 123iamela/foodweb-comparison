# Food web analysis Arctic & Atlantic

# load FW matrix and PA per polygon
#poly.pasp<- read.table("poly_pasp.txt", h=T)
TS<- read.table("TS.txt", h=T)
pw<- read.table("PW.txt", h=T)
FG<- read.table("FG.txt", h=T)
Regions<- read.table("Regions.txt", h=T)


#convert pairwise list to matrix

fw<- matrix(0, dim(TS)[1], dim(TS)[1])
row.names(fw)<-TS$TS
colnames(fw)<-TS$TS

# Convert from pairwise list to matrix
for (i in 1:dim(fw)[2]){
ID<-0
 for (j in 1:dim(pw)[1]){
	if(pw$PREDATOR[j]==as.vector(TS$TS[i])) {ID<-c(ID,which(row.names(fw)==pw$PREY[j]))}
	}
fw[ID,i]<-1
}

# Change matrix to presence/absence
fw1<-1*(fw>0)		# Multiply by 1 to convert the true and false statements in fw to 1 and 0 (presence/absence)
fw2<- fw1

# Load packages
library(igraph)		# Network analysis
library(NetIndices)	# Trophic level and omnivory 


# Taxa in Boreal and Arctic
Boreid<-which(Regions[1,]==1)
Arctid.b<-which(Regions[2,]==1)
Arctid.a<-which(Regions[3,]==1)

# Atlantic and Arctic food webs
B.FW<- fw1[Boreid, Boreid] #Boreal food web
Ab.FW<- fw1[Arctid.b, Arctid.b] # Arctic food web before
Aa.FW<-fw1[Arctid.a, Arctid.a] # Arctic food web after

# Which taxa belong to which functional group
FG<-cbind(rownames(fw1),FG)
FG.B<- FG[which(FG[,1]%in%rownames(B.FW)),]
FG.Ab<- FG[which(FG[,1]%in%rownames(Ab.FW)),] 
FG.Aa<- FG[which(FG[,1]%in%rownames(Aa.FW)),] 

# Creates an edge list of interactions
B.FW.edge.list<-graph.adjacency(B.FW)
Ab.FW.edge.list<-graph.adjacency(Ab.FW)
Aa.FW.edge.list<-graph.adjacency(Aa.FW)

# Trophic level as a plotting parameter.

TL.B.FW<-TrophInd(B.FW)
TL.Ab.FW<-TrophInd(Ab.FW)
TL.Aa.FW<-TrophInd(Aa.FW)

TL.B<- TL.B.FW$TL-1
TL.B[1]<-0
#TL.B[2]<-0
TL.Ab<-TL.Ab.FW$TL-1
TL.Ab[1]<-0
#TL.Ab[2]<-0
TL.Aa<-TL.Aa.FW$TL-1
TL.Aa[1]<-0
#TL.Aa[2]<-0
	

Cent.metrics.B<- data.frame(
deg.B=degree(B.FW.edge.list),
bet.B=betweenness(B.FW.edge.list) ,
eig.B=evcent(B.FW.edge.list)$vector ,
clo.B=closeness(B.FW.edge.list)
)

norm.deg.B<- Cent.metrics.B[,1]/dim(Cent.metrics.B)[1]
plot(Cent.metrics.B[,3]~norm.deg.B)
points(Cent.metrics.Ab[,3]~norm.deg.Ab, col="red")
points(Cent.metrics.Aa[,3]~norm.deg.Aa, col="grey")

power.law.fit(Cent.metrics.B[,2], xmin=NULL, start=2, force.continuous=FALSE, implementation=c("plfit", "R.mle"))

Cent.metrics.B[-c(1:10),4]
Cent.metrics.B[-c(1:10),4]>0.00005
FG.B

Cent.metrics.Ab<- data.frame(
deg.A=degree(Ab.FW.edge.list),
bet.A=betweenness(Ab.FW.edge.list) ,
eig.A=evcent(Ab.FW.edge.list)$vector ,
clo.A=closeness(Ab.FW.edge.list)
)

norm.deg.Ab<- Cent.metrics.Ab[,1]/dim(Cent.metrics.Ab)[1]
plot(Cent.metrics.Ab[,4]~norm.deg.Ab)
plot(Cent.metrics.Ab[,4]~TL.Ab)

Cent.metrics.Aa<- data.frame(
deg.A=degree(Aa.FW.edge.list),
bet.A=betweenness(Aa.FW.edge.list) ,
eig.A=evcent(Aa.FW.edge.list)$vector ,
clo.A=closeness(Aa.FW.edge.list)
)


norm.deg.Aa<- Cent.metrics.Aa[,1]/dim(Cent.metrics.Aa)[1]
plot(Cent.metrics.Aa[,2]~norm.deg.Aa)

h<-Cent.metrics.Aa[,4]>0.00005

FG.Aa[,1][h]

# creates a two-column matrix identifying the x and y values for each node.
layout.matrix.B<-matrix(nrow=length(V(B.FW.edge.list)), ncol=2)
layout.matrix.B[,1]<-runif(length(V(B.FW.edge.list))) # randomly assign along x-axis
layout.matrix.B[,2]<-TL.B # y-axis value based on trophic level
 
layout.matrix.Ab<-matrix(nrow=length(V(Ab.FW.edge.list)), ncol=2)
layout.matrix.Ab[,1]<-runif(length(V(Ab.FW.edge.list)))
layout.matrix.Ab[,2]<-TL.Ab
 
layout.matrix.Aa<-matrix(nrow=length(V(Aa.FW.edge.list)), ncol=2)
layout.matrix.Aa[,1]<-runif(length(V(Aa.FW.edge.list)))
layout.matrix.Aa[,2]<-TL.Aa

# plot food webs according to Trophic Level
SP.B <- spinglass.community(B.FW.edge.list, gamma=1) # simulated annealing
SP.Ab <- spinglass.community(Ab.FW.edge.list, gamma=1) # simulated annealing  
SP.Aa <- spinglass.community(Aa.FW.edge.list, gamma=1) # simulated annealing  

groups.B<- FG.B[,3] # 1=phytopl, 2=zoopl, 3=benthos, 4=fish, 5=birds, 6=mammals
groups.Ab<- FG.Ab[,3]# 1=phytopl, 2=zoopl, 3=benthos, 4=fish, 5=birds, 6=mammals
groups.Aa<- FG.Aa[,3]# 1=phytopl, 2=zoopl, 3=benthos, 4=fish, 5=birds, 6=mammals
colbar.FG<- c("lightgreen","cyan","orange","blue","yellow","lightpink", "black")

all.shapes <- setdiff(vertex.shapes(), "")
shapes <- c("circle", "square", "rectangle", "triangle")
shapes <- c("circle", "square", "rectangle", "triangle", "circle")

par(mar=c(.1,.1,.1,.1),mfrow=c(1,1))


plot(B.FW.edge.list, layout=layout.matrix.B, vertex.size = 1/5*degree(B.FW.edge.list)+4, edge.arrow.size=.10, main ="", vertex.color=colbar.FG[groups.B],vertex.label = "", edge.color="grey70", edge.width=0.05)
#plot(B.FW.edge.list, layout=layout.circle, vertex.size = 1/5*degree(B.FW.edge.list)+4, edge.arrow.size=.10, main ="", vertex.color=colbar.FG[groups.B],vertex.label = "", edge.color="black")
plot(Ab.FW.edge.list, layout=layout.matrix.Ab, vertex.size = 1/5*degree(Ab.FW.edge.list)+4, edge.arrow.size=.10, main ="", vertex.color=colbar.FG[groups.Ab],vertex.label = "", edge.color="grey70", edge.width=0.05)

plot(Ab.FW.edge.list, layout=layout.matrix.Ab, vertex.size = 1/5*degree(Ab.FW.edge.list)+4, edge.arrow.size=.10, main ="", vertex.color=colbar.FG[groups.Ab],vertex.label = "", edge.color="grey70", edge.width=0.05)
plot(Aa.FW.edge.list, layout=layout.matrix.Aa, vertex.size = 1/5*degree(Aa.FW.edge.list)+4, edge.arrow.size=.10, main ="", vertex.color=colbar.FG[groups.Aa],vertex.label = "", edge.color="grey70", edge.width=0.05)










