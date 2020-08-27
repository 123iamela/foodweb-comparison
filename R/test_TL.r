# Porque dan distintos los niveles troficos si se calcula a partir de la matriz GLV
#
b1 <- readNetwork("Data/Beagle_FW.txt", edgeListFormat=1)
b <- simplify(b1)
E(b1)[which_multiple(b1)]
E(b1)[which_loop(b1)]

difference(b1,b)

cbRed <- list(make_empty_graph(n=vcount(b)), make_empty_graph(n=vcount(b)), b) 

# Crea la matriz de adjacencia GLV

cbRed <- fromIgraphToMgraph(cbRed, c("empty","empty","Antagonistic"))
A <- toGLVadjMat(cbRed, c("empty","empty","Antagonistic"))

# Prueba de que sea una matriz de interaccion antagonica unica = Red trofica

calcPropInteractionsGLVadjMat(A, rep(1,times=nrow(A)))==c(0,1,0,0,0)

# Comprueba que los predadores y las presas estan bien codificados

A[57,58] # Presa    -1
A[58,57] # Predador  1 

V(b)[57]
V(b)[58]
V(b)[nei(V(b)[57],"out")]
V(b)[nei(V(b)[57],"in")]

V(b)[nei(V(b)[58],"in")]

A_Beagle <- A


gP1  <-  ifelse(A_Beagle<0, 1, 0)
diag(gP1) <- 0
tlp1 <- TrophInd(gP1)$TL
tlb <- TrophInd(get.adjacency(b,sparse=F))$TL
round(tlb-tlp1,5)
gP <- get.adjacency(b,sparse=F)
gP[57,58]
gP[58,57]
#### TA MAL gP!!!!!!!!!!!!!!!!!!
all(gP==gP1)
which(gP1 != gP, arr.ind=TRUE)



gP[57,58]
gP1[57,58]
A_Beagle[57,58]
A_Beagle[58,57]
