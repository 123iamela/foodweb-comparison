library(igraph)

spingB<-spinglass.community(B.FW.edge.list, spins=25, gamma=1)
spingB.mem<- spingB$membership
spingAB<-spinglass.community(Ab.FW.edge.list, spins=25, gamma=1)
spingAB.mem<- spingAB$membership
spingAA<-spinglass.community(Aa.FW.edge.list, spins=25, gamma=1)
spingAA.mem<- spingAA$membership

Cod<-which(FG.B=="GAD_MOR",arr.ind=TRUE) # find col and row no of species
Haddock<-which(FG.B=="MEL_AEG",arr.ind=TRUE) # find col and row no of species

# Degree of all species
deg<- as.numeric(degree(B.FW.edge.list))

# within-module degree
l<-vector()
memMod<-vector()
for (i in 1:dim(FG.B)[1]){
		
	sp.in.deg <- V(B.FW.edge.list)[nei(i, "in")]
	sp.out.deg<- V(B.FW.edge.list)[nei(i, "out")]
	mem.sp<-spingB.mem[i]
	k<- length(which(spingB.mem[c(sp.in.deg, sp.out.deg)]==mem.sp))
	mem<- which(spingB.mem==mem.sp)
	
	for (m in 1:length(mem)){
	mem.in.deg <- V(B.FW.edge.list)[nei(mem[m], "in")]
	mem.out.deg<- V(B.FW.edge.list)[nei(mem[m], "out")]
	memMod.id<- length(which(spingB.mem[c(mem.in.deg, mem.out.deg)]==mem.sp))
	memMod[m]<- memMod.id
	}
	
	k.ave<- mean(memMod)
	k.sd<- sd(memMod)
	l[i]<- (k-k.ave)/k.sd
}

# among module connectivity
r<- vector()
for (i in 1:dim(FG.B)[1]){
	
	d<-degree(B.FW.edge.list)[i]
	sp.in.deg <- V(B.FW.edge.list)[nei(i, "in")]
	sp.out.deg<- V(B.FW.edge.list)[nei(i, "out")]
	mod.sp<-table(spingB.mem[c(sp.in.deg, sp.out.deg)])
	mod.no<-as.numeric(names(mod.sp))
	mod<-rep(0, length(unique(spingB.mem)))
	mod[mod.no]<-c(mod.sp)
	r[i]<- 1-(sum((mod/d)^2))
}

write.table(l, "l.emp.txt")
write.table(r, "r.emp.txt")


# Plot
colbar<- c("red","lightgreen","orange","blue","magenta","black")
colbar.FG<- c("grey","lightgreen","cyan","orange","blue","magenta","lightpink")
groups.B<- FG.B[,2] # 1=phytopl, 2=zoopl, 3=benthos, 4=fish, 5=birds, 6=mammals

par(mfrow=c(1,1))
par(oma=c(2.7,1.3,0.7,1)) # change outer figure margins
par(mai=c(1,1,0.7,0.7)) # size of figure margins

plot(l~r, type="n", axes=T ,tck=F, lwd=2, ann=F, cex.axis=1.2, xlim=c(0,0.8), ylim=c(-1.5,5))
lines(c(0.625,0.625), c(-1.5, 5), col="grey")
lines(c(0,0.8), c(2.5, 2.5), col="grey")
points(r, l, col=colbar.FG[groups.B], pch=20, cex=2)
mtext(2, text="Within module degree", line=4,font=2, cex=1.2)
mtext(1, text="Among module connectivity",line=4, font=2, cex=1.2)
axis(1, tck=0.02, lwd=1,las=2,lty=1, labels=F, xlim=c(0,1))
axis(2,tck=0.02, labels=FALSE)
points(r[128], l[128], cex=4, col="blue", pch=20)
points(r[144], l[144], cex=4, col="blue", pch=20)

#### Arctic before

# within-module degree
ll<-vector()
memMod2<-vector()
for (i in 1:dim(FG.Ab)[1]){
		
	sp.in.deg <- V(Ab.FW.edge.list)[nei(i, "in")]
	sp.out.deg<- V(Ab.FW.edge.list)[nei(i, "out")]
	mem.sp<-spingAB.mem[i]
	k<- length(which(spingAB.mem[c(sp.in.deg, sp.out.deg)]==mem.sp))
	mem<- which(spingAB.mem==mem.sp)
	
	for (m in 1:length(mem)){
	mem.in.deg <- V(Ab.FW.edge.list)[nei(mem[m], "in")]
	mem.out.deg<- V(Ab.FW.edge.list)[nei(mem[m], "out")]
	memMod2[m]<- length(which(spingAB.mem[c(mem.in.deg, mem.out.deg)]==mem.sp))
	}

	k.ave<- mean(memMod2)
	k.sd<- sd(memMod2)
	ll[i]<- (k-k.ave)/k.sd
}

# among module connectivity
rr<- vector()
for (i in 1:dim(FG.Ab)[1]){
	
	d<-degree(Ab.FW.edge.list)[i]
	sp.in.deg <- V(Ab.FW.edge.list)[nei(i, "in")]
	sp.out.deg<- V(Ab.FW.edge.list)[nei(i, "out")]
	mod.sp<-table(spingB.mem[c(sp.in.deg, sp.out.deg)])
	mod.no<-as.numeric(names(mod.sp))
	mod<-rep(0, length(unique(spingAB.mem)))
	mod[mod.no]<-c(mod.sp)
	rr[i]<- 1-(sum((mod/d)^2))
}

# Plot
colbar<- c("red","lightgreen","orange","blue","magenta","black")
colbar.FG<- c("grey","lightgreen","cyan","orange","blue","magenta","lightpink")
groups.Ab<- FG.B[,2] # 1=phytopl, 2=zoopl, 3=benthos, 4=fish, 5=birds, 6=mammals

#par(mfrow=c(1,1))
#par(oma=c(2.7,1.3,0.7,1)) # change outer figure margins
#par(mai=c(1,1,0.7,0.7)) # size of figure margins

plot(ll~rr, type="n", axes=T ,tck=F, lwd=2, ann=F, cex.axis=1.2, xlim=c(0,0.8), ylim=c(-1.5,5))
lines(c(0.625,0.625), c(-1.5, 5), col="grey")
lines(c(0,0.8), c(2.5, 2.5), col="grey")
points(rr, ll, col=colbar.FG[groups.Ab], pch=20, cex=2)
mtext(2, text="Within module degree", line=4,font=2, cex=1.2)
mtext(1, text="Among module connectivity",line=4, font=2, cex=1.2)
axis(1, tck=0.02, lwd=1,las=2,lty=1, labels=F)
axis(2,tck=0.02, labels=FALSE)

#points(r[128], l[128], cex=4, col="blue", pch=20)
#points(r[144], l[144], cex=4, col="blue", pch=20)


# Degree of all species
# within-module degree
lll<-vector()
memMod3<-vector()
for (i in 1:dim(FG.Aa)[1]){
		
	sp.in.deg <- V(Aa.FW.edge.list)[nei(i, "in")]
	sp.out.deg<- V(Aa.FW.edge.list)[nei(i, "out")]
	mem.sp<-spingAA.mem[i]
	k<- length(which(spingAA.mem[c(sp.in.deg, sp.out.deg)]==mem.sp))
	mem<- which(spingAA.mem==mem.sp)
	
	for (m in 1:length(mem)){
	mem.in.deg <- V(Aa.FW.edge.list)[nei(mem[m], "in")]
	mem.out.deg<- V(Aa.FW.edge.list)[nei(mem[m], "out")]
	memMod3[m]<- length(which(spingAA.mem[c(mem.in.deg, mem.out.deg)]==mem.sp))
	}

	k.ave<- mean(memMod2)
	k.sd<- sd(memMod2)
	lll[i]<- (k-k.ave)/k.sd
}

# among module connectivity
rrr<- vector()
for (i in 1:dim(FG.Aa)[1]){
	
	d<-degree(Aa.FW.edge.list)[i]
	sp.in.deg <- V(Aa.FW.edge.list)[nei(i, "in")]
	sp.out.deg<- V(Aa.FW.edge.list)[nei(i, "out")]
	mod.sp<-table(spingAA.mem[c(sp.in.deg, sp.out.deg)])
	mod.no<-as.numeric(names(mod.sp))
	mod<-rep(0, length(unique(spingAA.mem)))
	mod[mod.no]<-c(mod.sp)
	rrr[i]<- 1-(sum((mod/d)^2))
}

# Plot

colbar<- c("red","lightgreen","orange","blue","magenta","black")
colbar.FG<- c("grey","lightgreen","cyan","orange","blue","magenta","lightpink")
groups.AA<- FG.Aa[,2] # 1=phytopl, 2=zoopl, 3=benthos, 4=fish, 5=birds, 6=mammals

plot(lll~rrr, type="n", axes=T ,tck=F, lwd=2, ann=F, cex.axis=1.2, xlim=c(0,0.8), ylim=c(-1.5,5))
lines(c(0.625,0.625), c(-1.5, 5), col="grey")
lines(c(0,0.8), c(2.5, 2.5), col="grey")
points(rrr, lll, col=colbar.FG[groups.AA], pch=20, cex=2)
mtext(2, text="Within module degree", line=4,font=2, cex=1.2)
mtext(1, text="Among module connectivity",line=4, font=2, cex=1.2)
axis(1, tck=0.02, lwd=1,las=2,lty=1, labels=F)
axis(2,tck=0.02, labels=FALSE)
points(rrr[which(FG.Aa=="GAD_MOR",arr.ind=TRUE)[1]], lll[which(FG.Aa=="GAD_MOR",arr.ind=TRUE)[1]], cex=4, col="blue", pch=20)
points(rrr[which(FG.Aa=="MEL_AEG",arr.ind=TRUE)[1]], lll[which(FG.Aa=="MEL_AEG",arr.ind=TRUE)[1]], cex=4, col="blue", pch=20)

### Module colored
par(mfrow=c(1,3))
colmod<-c("red","green","blue","black"

plot(l~r, type="n", axes=T ,tck=F, lwd=2, ann=F, cex.axis=1.2, xlim=c(0,0.8), ylim=c(-1.5,5))
lines(c(0.625,0.625), c(-1.5, 5), col="grey")
lines(c(0,0.8), c(2.5, 2.5), col="grey")
points(r, l, col=colmod[spingB.mem], pch=20, cex=2)
mtext(2, text="Within module degree", line=4,font=2, cex=1.2)
mtext(1, text="Among module connectivity",line=4, font=2, cex=1.2)
axis(1, tck=0.02, lwd=1,las=2,lty=1, labels=F, xlim=c(0,1))
axis(2,tck=0.02, labels=FALSE)
points(r[128], l[128], cex=4, col="blue", pch=20)
points(r[144], l[144], cex=4, col="blue", pch=20)

plot(ll~rr, type="n", axes=T ,tck=F, lwd=2, ann=F, cex.axis=1.2, xlim=c(0,0.8), ylim=c(-1.5,5))
lines(c(0.625,0.625), c(-1.5, 5), col="grey")
lines(c(0,0.8), c(2.5, 2.5), col="grey")
points(rr, ll, col=colmod[spingAB.mem], pch=20, cex=2)
mtext(2, text="Within module degree", line=4,font=2, cex=1.2)
mtext(1, text="Among module connectivity",line=4, font=2, cex=1.2)
axis(1, tck=0.02, lwd=1,las=2,lty=1, labels=F)
axis(2,tck=0.02, labels=FALSE)

plot(lll~rrr, type="n", axes=T ,tck=F, lwd=2, ann=F, cex.axis=1.2, xlim=c(0,0.8), ylim=c(-1.5,5))
lines(c(0.625,0.625), c(-1.5, 5), col="grey")
lines(c(0,0.8), c(2.5, 2.5), col="grey")
points(rrr, lll, col=colmod[spingAA.mem], pch=20, cex=2)
mtext(2, text="Within module degree", line=4,font=2, cex=1.2)
mtext(1, text="Among module connectivity",line=4, font=2, cex=1.2)
axis(1, tck=0.02, lwd=1,las=2,lty=1, labels=F)
axis(2,tck=0.02, labels=FALSE)
points(rrr[which(FG.Aa=="GAD_MOR",arr.ind=TRUE)[1]], lll[which(FG.Aa=="GAD_MOR",arr.ind=TRUE)[1]], cex=4, col="blue", pch=20)
points(rrr[which(FG.Aa=="MEL_AEG",arr.ind=TRUE)[1]], lll[which(FG.Aa=="MEL_AEG",arr.ind=TRUE)[1]], cex=4, col="blue", pch=20)


#### Among module connectivity as a function of degree
# Boreal
plot(r ~ as.vector(deg.B.in),  type="n", xlim=c(0,90),axes=T ,tck=F, lwd=2, ann=F, cex.axis=1.2)
lines(c(0, 90), c(0.625,0.625), col="darkgrey", lwd=2)
text(r~ as.vector(deg.B.in), labels=as.character(names.B), cex=0.5, col=colbar[groups.B]) 
mtext(2, text="Among module connectivity", line=4,font=2, cex=1.2)
mtext(1, text="In-degree",line=4, font=2, cex=1.2)
axis(1, tck=0.02, lwd=1,las=2,lty=1, labels=F, xlim=c(0,1))
axis(2,tck=0.02, labels=FALSE)

plot(r~ deg.B.out, type="n", xlim=c(0,90),axes=T ,tck=F, lwd=2, ann=F, cex.axis=1.2)
lines(c(0, 90), c(0.625,0.625), col="darkgrey", lwd=2)
text(r~ as.vector(deg.B.out), labels=as.character(names.B), cex=0.5, col=colbar[groups.B]) 
mtext(2, text="Among module connectivity", line=4,font=2, cex=1.2)
mtext(1, text="Out-degree",line=4, font=2, cex=1.2)
axis(1, tck=0.02, lwd=1,las=2,lty=2, labels=F, xlim=c(0,1))
axis(2,tck=0.02, labels=FALSE)


# Arctic
plot(rr~ as.vector(deg.A.in), col=colbar.FG[groups.Ab], type="n", xlim=c(0,60),axes=T ,tck=F, lwd=2, ann=F, cex.axis=1.2)
lines(c(0, 60), c(0.625,0.625), col="darkgrey", lwd=2)
text(rr~ as.vector(deg.A.in), labels=as.character(names.Ab), cex=0.5, col=colbar[groups.Ab]) 
mtext(2, text="Among module connectivity", line=4,font=2, cex=1.2)
mtext(1, text="In-degree",line=4, font=2, cex=1.2)
axis(1, tck=0.02, lwd=1,las=2,lty=1, labels=F, xlim=c(0,1))
axis(2,tck=0.02, labels=FALSE)

plot(rr~ as.vector(deg.A.out), col=colbar.FG[groups.Ab], type="n", xlim=c(0,60),axes=T ,tck=F, lwd=2, ann=F, cex.axis=1.2)
lines(c(0, 60), c(0.625,0.625), col="darkgrey", lwd=2)
text(rr~ as.vector(deg.A.out), labels=as.character(names.Ab), cex=0.5, col=colbar[groups.Ab]) 
mtext(2, text="Among module connectivity", line=4,font=2, cex=1.2)
mtext(1, text="Out-degree",line=4, font=2, cex=1.2)
axis(1, tck=0.02, lwd=1,las=2,lty=1, labels=F, xlim=c(0,1))
axis(2,tck=0.02, labels=FALSE)


# Arctic after
plot(rrr~ as.vector(deg.Aa.in), col=colbar.FG[groups.Aa], type="n", xlim=c(0,80),axes=T ,tck=F, lwd=2, ann=F, cex.axis=1.2)
lines(c(0, 80), c(0.625,0.625), col="darkgrey", lwd=2)
text(rrr~ as.vector(deg.Aa.in), labels=as.character(names.Aa), cex=0.5, col=colbar[groups.Aa]) 
mtext(2, text="Among module connectivity", line=4,font=2, cex=1.2)
mtext(1, text="In-degree",line=4, font=2, cex=1.2)
axis(1, tck=0.02, lwd=1,las=2,lty=1, labels=F, xlim=c(0,1))
axis(2,tck=0.02, labels=FALSE)

plot(rrr~ as.vector(deg.Aa.out), col=colbar.FG[groups.Aa], type="n", xlim=c(0,80),axes=T ,tck=F, lwd=2, ann=F, cex.axis=1.2)
lines(c(0, 80), c(0.625,0.625), col="darkgrey", lwd=2)
text(rrr~ as.vector(deg.Aa.out), labels=as.character(names.Aa), cex=0.5, col=colbar[groups.Aa]) 
mtext(2, text="Among module connectivity", line=4,font=2, cex=1.2)
mtext(1, text="Out-degree",line=4, font=2, cex=1.2)
axis(1, tck=0.02, lwd=1,las=2,lty=1, labels=F, xlim=c(0,1))
axis(2,tck=0.02, labels=FALSE)






































fish.list<- c("GAD_MOR","MEL_AEG")
gammaRegional.B<- matrix(0,2,10)
rownames(gammaRegional.B)<-c("GAD_MOR","MEL_AEG")
colnames(gammaRegional.B)<-seq(0.5, 1.4, 0.1)

gammaLocal.B<- matrix(0,2,10)
rownames(gammaLocal.B)<-c("GAD_MOR","MEL_AEG")
colnames(gammaLocal.B)<-seq(0.5, 1.4, 0.1)


for (i in 1:length(fish.list)){
	for (j in 1:10){
	spin.gamma<-0.4+j/10
	spingB <- spinglass.community(B.FW.edge.list, spins=25, gamma=spin.gamma)
	
	pos.fish<-match(fish.list[i], FG.B[,1])		#find position of cod	
	mem.fish<-spingB$membership[pos.fish]	

	sp.in.deg <- V(B.FW.edge.list)[nei(pos.fish, "in")]
	sp.out.deg<- V(B.FW.edge.list)[nei(pos.fish, "out")]
	
	k<- length(which(spingB$membership[c(sp.in.deg, sp.out.deg)]==mem.fish))
	mem<- which(spingB$membership==mem.fish)
	k.ave<- mean(deg[mem])
	k.sd<- sd(deg[mem])

	d<-degree(B.FW.edge.list)[pos.fish]
	mod.fish<-table(spingB$membership[c(sp.in.deg, sp.out.deg)])
	mod.no<-as.numeric(names(mod.fish))
	mod<-rep(0, length(unique(spingB$membership)))
	mod[mod.no]<-c(mod.fish)
	
	gammaLocal.B[i,j]<- (k-k.ave)/k.sd
	gammaRegional.B[i,j]<- 1-(sum((mod/d)^2))
}
}

par(mfrow=c(1,1))
plot(gammaRegional.B[1,], gammaLocal.B[1,], type="n", xlab="among module connectivity", ylab="within module degree", main="Cod", cex=2)
text(gammaRegional.B[1,], gammaLocal.B[1,],label=seq(0.5, 1.4, 0.1))
plot(gammaRegional.B[2,], gammaLocal.B[2,], type="n", xlab="among module connectivity", ylab="within module degree", main="Haddock", cex=2)
text(gammaRegional.B[2,], gammaLocal.B[2,],label=seq(0.5, 1.4, 0.1))











### Fish in Boreal food web

xx<-c(0,1,0,-1,0)
yy<-c(1,0,-1,0,0)
layout<-cbind(xx,yy)

fish.list<- c("GAD_MOR", "MEL_AEG", "SEB_MEN", "MAL_VIL")

par(mfrow=c(2,2))

pdf()  
for(i in 1:length(fish.list)){
pos.fish<-match(fish.list[i], FG.B[,1])		#find position of cod
mem.fish<-spingB.mem[pos.fish]				#which module is cod in
fish.in.deg <- V(B.FW.edge.list)[ nei( pos.fish, "in") ]
fish.out.deg<- V(B.FW.edge.list)[ nei( pos.fish, "out") ]

deg.size.fish<-table(spingB.mem[c(fish.in.deg, fish.out.deg)])
node.size.fish<-table(spingB.mem)
fish.mtx<-matrix(0,length(node.size.fish)+1,length(node.size.fish)+1)
idx<-as.numeric(names(deg.size.fish))
fish.mtx[,5]<- rep(0,5)
fish.mtx[idx,5]<- 1
node.size.fish<-c(as.vector(node.size.fish),30)
fish.edge.list<-graph.adjacency(fish.mtx)

col.mod<-colors()[c(652,468,259,636)]
plot.igraph(fish.edge.list,vertex.label=c(1,2,3,4,mem.fish), vertex.label.color= "black", vertex.label.cex=1.2,vertex.size=node.size.fish,
edge.arrow.size=.1,edge.width=deg.size.fish/2,layout=layout,vertex.color=col.mod[c(1,2,3,4,mem.fish)], main=fish.list[i])
}
dev.off() 

### Fish in Arctic now food web


xx<-c(0,1,0,-1,0)
yy<-c(1,0,-1,0,0)
layout<-cbind(xx,yy)

fish.list<- c("GAD_MOR", "MEL_AEG", "SEB_MEN", "MAL_VIL")

par(mfrow=c(2,2))

pdf()  
for(i in 1:length(fish.list)){
pos.fish<-match(fish.list[i], FG.A.now[,1])		#find position of cod
mem.fish<-spingA.now.mem[pos.fish]				#which module is cod in
fish.in.deg <- V(A.now.edge.list)[ nei( pos.fish, "in") ]
fish.out.deg<- V(A.now.edge.list)[ nei( pos.fish, "out") ]

deg.size.fish<-table(spingA.now.mem[c(fish.in.deg, fish.out.deg)])
node.size.fish<-table(spingA.now.mem)
fish.mtx<-matrix(0,length(node.size.fish)+1,length(node.size.fish)+1)
idx<-as.numeric(names(deg.size.fish))
fish.mtx[,5]<- rep(0,5)
fish.mtx[idx,5]<- 1
node.size.fish<-c(as.vector(node.size.fish),30)
fish.edge.list<-graph.adjacency(fish.mtx)

col.mod<-colors()[c(652,468,259,636)]
plot.igraph(fish.edge.list,vertex.label=c(1,2,3,4,mem.fish), vertex.label.color= "black", vertex.label.cex=1.2,vertex.size=node.size.fish,
edge.arrow.size=.1,edge.width=deg.size.fish/2,layout=layout,vertex.color=col.mod[c(1,2,3,4,mem.fish)], main=fish.list[i])
}
dev.off() 














#Arctic # Boregadus saida
pos.cod<-match("BOR_SAI", FG.A[,1])		#find position of cod
mem.cod<-spingA.mem[pos.cod]				#which module is cod in

cod.in.deg <- V(A.FW.edge.list)[ nei( 118, "in") ]
cod.out.deg<- V(A.FW.edge.list)[ nei( 118, "out") ]

hist(spingA.mem[c(cod.in.deg, cod.out.deg)], breaks = seq(0, 4, by = 1))

# which fish are in module
list.mem.fish<- spingA.mem[FG.A$FG==4]
list.fish<- FG.A[,1][FG.A$FG==4]

fish.aff.gamma0.5<-cbind(as.character(list.fish),list.mem.fish)

write.table(fish.aff.gamma0.5, "fish.affiliation.gamma0.5.txt")

