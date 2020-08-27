## Codigo para generar matrices aleatorias y calcular el "characteristic path length" 
# y el "clustering coefficient" 

library("igraph")
library("xlsx")
library("plyr")

## Si queremos cargar una red lo debemos hacer como una matriz de adyacencia 
# (nombre de nodos en 1ra fila y columna) que luego sera convertida en grafo para analizarla

red.e <- read.csv ("C:/Users/usuario/Desktop/Redes_empiricas/Beagle_Channel_FW.csv", row.names = 1, header = TRUE)
red.e <- as.matrix(red.e) 
red.e <- graph_from_adjacency_matrix(red.e) 


## Calculo del "characteristic path length" y el "clustering coefficient"  
# type="local" o "global" en transitivity ("calculo del clustering coefficient") 
# para nodos individuales o para toda la red

cha.path.e <- average.path.length(red.e)
clus.coef.e <- transitivity(red.e, type = "global") 

## Creacion de 1000 redes aletorias con S y L igual a la red empirica
# Calculo del "characteristic path length" y el "clustering coefficient" 

redes.r <- lapply(1:1000, function (x) erdos.renyi.game(33, 183, type="gnm")) 
#cha.path.e = 2.30
#clus.coef.e = 0.30

cha.path.r <- c() 
for (i in 1:1000){
  cha.path.r <- c(cha.path.r,average.path.length(redes.r[[i]]))
}

clus.coef.r <- c() 
for (i in 1:1000){
  clus.coef.r <- c(clus.coef.r,transitivity(redes.r[[i]], type="Global"))
}

## Calculo de cocientes empiric/random para "characteristic path length" (lambda) 
# y "clustering coefficient" (gamma)  

lambda <- cha.path.e/cha.path.r
lambda.mean <- mean(lambda)
lambda.mean
lambda.sd <- sd(lambda)
lambda.sd

gamma <- clus.coef.e/clus.coef.r
gamma.mean <- mean(gamma)
gamma.mean
gamma.sd <- sd(gamma)
gamma.sd

## Exportar datos a una planilla .xlsx

cha.path.res <- write.xlsx(cha.path.r, "c:/Users/usuario/Desktop/cha.data.xlsx")
clus.coef.res <- write.xlsx(clus.coef.r, "c:/Users/usuario/Desktop/clus.data.xlsx")

## Histograma de "characteristic path length" y "clustering coefficient" en redes 
# aleatorias con S y C igual a Potter Cove. Ver en Navia et al. (2016).

hist(cha.path.r, breaks = 6, freq = F, col = "gray", 
     ylim = c(0, 50), xlim = c(1.7, 2.8), 
     xlab = "Characteristic Path Length (CPL)", ylab = "Frequency", main = NA)
bold <- rgb(0, 0, 0)
arrows(1.82, 15, 1.82, 0, length = 0.15, angle = 30, lty = 1, lwd = 2, col = bold)
text(1.97, 17, labels = "CPL obs = 1.82", cex = 0.9, font = 2)

hist(clus.coef.r, breaks = 6, freq = F, col = "gray",
     ylim = c(0, 50), xlim = c(0, 0.12), 
     xlab = "Clustering Coefficient (CC)", ylab = "Frequency", main = NA)
arrows(0.096, 20, 0.096, 8, length = 0.15, angle = 30, lty = 1, lwd = 2, col = bold)
text(0.098, 25, labels = "CC obs = 0.096", cex = 0.9, font = 2)

## Grafico bivariado con x=lambda e y=gamma para cada red trofica

marinas <- read.table ("c:/Users/usuario/Desktop/Marine_fw.csv", sep=",", header = T)
plot(marinas$xmean, marinas$ymean, type = "p", cex = 0.9, col = c(marinas$food.web), pch = 16, 
     xlab="CPL / CPL random", ylab="CC / CC random", xlim=c(0,2), ylim=c(0,10))
legend("topright", title = "FOOD WEBS", legend = as.vector(marinas$food.web), 
       pch=16, cex = 0.6, col = unique(marinas$food.web))
arrows(marinas$xmean-marinas$xsd, marinas$ymean, marinas$xmean+marinas$xsd, marinas$ymean, 
       length=0.025, angle=90, code=3)
arrows(marinas$xmean, marinas$ymean-marinas$ysd, marinas$xmean, marinas$ymean+marinas$ysd, 
       length=0.025, angle=90, code=3)
with(subset(marinas,food.web=="Potter Cove"), text(xmean,ymean,food.web, cex = 0.8, pos = 1))

# Colorear de manera transparente el rectangulo de small-world (L/Lrandom = 1, C/Crandom > 1)
# y agregar texto

col2rgb("blue", alpha=TRUE)
bluetrans <- rgb(0, 0, 255, 127, maxColorValue = 255)
rect(0.75, 3, 1.25, 8, col = bluetrans, lty = 3)
bold <- rgb(0, 0, 0)
text(1, 7.5, "Small-world", col = bold)

## Grafico log-log de S (riqueza) VS C/Crandom.
# Redes troficas marinas
marinas <- read.table ("c:/Users/usuario/Desktop/Marine_fw.csv", sep=",", header = T)
log.size.m <- log10(marinas$S)
log.clus.m <- log10(marinas$ymean)
plot(log.size.m, log.clus.m, type = "p", cex = 0.9, col = c(marinas$food.web), pch = 16, 
     xlab="S richness", ylab="C / C random", xlim=c(0,2.5), ylim=c(-0.5,0.5))
legend("topleft", title = "FOOD WEBS", legend = as.vector(marinas$food.web), 
       pch=16, cex = 0.6, col = unique(marinas$food.web))

fit.m <-lm(log.clus.m ~ log.size.m)
summary(fit.m)
segments(1, -0.4111+0.2497*1, 2.5, -0.4111+0.2497*2.5)
legend(1.25, -0.3, legend = "y = -0.41 + 0.25x", bty = "n", cex = 0.85)
legend(1.25, -0.37, legend = "R squared = 0.16", bty = "n", cex = 0.85)

# Redes troficas marinas y no marinas
todas <- read.table ("c:/Users/usuario/Desktop/All_fw.csv", sep=",", header = T)
log.size <- log10(todas$S)
log.clus <- log10(todas$ymean)
plot(log.size, log.clus, type = "p", cex = 0.9, pch = c(todas$ecosystem),
     xlab = "S richness", ylab = "C / C random", xlim = c(0,2.5), ylim = c(-0.5,0.5))
legend("topleft", title = "FOOD WEBS", c("Marine","Non-marine"),
       pch = c(21,24), cex = 0.6)

# Ajuste a un modelo lineal 
fit <- lm(log.clus ~ log.size)
summary(fit)
segments(1, -0.4692+0.2051*1, 2.5, -0.4692+0.2051*2.5)
