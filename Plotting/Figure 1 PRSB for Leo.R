### Food web Network Figures ###
#FG.B is a list of spcies with their functional affiliation
#B.FW.edge.list is the edgelists of your food webs matrices created in igraph

# creates a two-column matrix identifying the x and y values for each node.
layout.matrix.B <- matrix(nrow=length(V(b)), ncol=2)

spingB.mem <- mod_by_red[[2]]$membership

grB1 <- runif(length(which(spingB.mem==3)), min = 1.0, max = 1.3 )
grB2 <- runif(length(which(spingB.mem==1)), min = 0, max = 0.3)
grB3 <- runif(length(which(spingB.mem==2)), min = 0.5, max = 0.8)
grB4 <- runif(length(which(spingB.mem==4)), min = 1.5, max = 1.8)

spingB1 <- replace(spingB.mem, spingB.mem==3, grB1)
spingB2 <- replace(spingB1, spingB.mem==1, grB2)
spingB3 <- replace(spingB2, spingB.mem==2, grB3)
spingB4 <- replace(spingB3, spingB.mem==4, grB4)

layout.matrix.B[,1] <- spingB4 # randomly assign along x-axis
layout.matrix.B[,2] <- troph.net2$TL # y-axis value based on trophic level NECESITO CALCULAR EL TL Y QUE COINCIDA CON EL ORDEN DE LAS SPP POR MODULO

kb <- read.delim(file = "Data/bmod_rol_k.txt", stringsAsFactors = FALSE)

groups.B <- kb$Roln # 5=hubconn, 6=modconn, 7=modspe
#groups.Ab<- FG.Ab[,2]# 1=non-organic,2=phytopl, 3=zoopl, 4=benthos, 5=fish, 6=birds, 7=mammals
colbar.FG<- c("cyan","orange","blue")

par(mai=c(0,0,0,0))

par(mar=c(.5,.5,.5,.5),mfrow=c(1,1))

# Plot
plot(b, layout=layout.matrix.B, vertex.size = 1/5*degree(b)+4, vertex.frame.color = "gray70", edge.arrow.size=F,  vertex.color=colbar.FG[groups.B],vertex.label = NA, edge.color="grey70", edge.width=0.15)

#plot(Ab.FW.edge.list, layout=layout.matrix.Ab, vertex.size = 1/5*degree(Ab.FW.edge.list)+4, vertex.frame.color ="gray50",edge.arrow.size=0.001, vertex.color=colbar.FG[groups.Ab],vertex.label = NA, edge.color="grey40", edge.width=0.15)
#plot(Aa.FW.edge.list, layout=layout.matrix.Aa, vertex.size = 1/5*degree(Aa.FW.edge.list)+4, vertex.frame.color = "gray50",edge.arrow.size=0.001, vertex.color=colbar.FG[groups.Aa],vertex.label = NA, edge.color="grey70", edge.width=0.05)

#
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("topleft", legend=c("non-living","plankton","zooplankton", "benthos", "fish", "sea birds", "sea mammals"), xpd = TRUE, horiz = F, bty = "n", pch = c(rep(16, 7)), col = colbar.FG, cex=1)
inset = c(0, 0)
#c("non-living","plankton","zooplankton", "benthos", "fish", "sea birds", "sea mammals")