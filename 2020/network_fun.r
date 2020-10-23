# Plot modules vs trophic level
#
#
plot_modules_TL <- function(redl,modulos)
{
  #
  # Calculate Trophic Position
  #
  require(NetIndices)
  troph.net2<-TrophInd(get.adjacency(redl,sparse=F))
  layout.matrix.1<-matrix(
    nrow=length(V(redl)),  # Rows equal to the number of vertices
    ncol=2
  )
  
  # 
  # Add colors to nodes 
  #
  require(RColorBrewer)
  
  colTL <-as.numeric(cut(troph.net2$TL,11))
  colnet <- brewer.pal(11,"RdYlGn")
  V(redl)$color <- colnet[12-colTL]
  
  #
  # Plot modules
  #
  layout.matrix.1[,2]<-jitter(troph.net2$TL,0.1) # y-axis value based on trophic level
  layout.matrix.1[,1]<-jitter(modulos$membership,1) # randomly assign along x-axis
  
#  png("Figures/RegionalFoodWeb.png",width=6,height=6,units="in",res=600)
#  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
  
  plot(redl, vertex.color=vertex_attr(redl)$cor, #vertex.label=NA,
       #vertex.size=log(3*igraph::degree(redl)),
       edge.width=.3,edge.arrow.size=.4, 
       edge.color="grey50",
       edge.curved=0.3, layout=layout.matrix.1)
  
}

#' Calculate motif counts for observed network and CI for erdos-renyi random networks and Z-scores 
#'
#' @param red igraph network object
#' @param redes list with networks used as a null model with the same number nodes and links
#' @param ncores number of cores for parallel computation, if null do not use a parallel 
#'
#' @return data.frame with all the results
#' @export
#'
#' @examples
calc_motif_zscores <- function(red, redes, ncores=NULL)
{
    Size <- vcount(red)
    Links <- ecount(red)

    if(!is.null(ncores)) {
      cn <-parallel::detectCores()
      if(cn>ncores)
        cn <- ncores
      else
        cn <- cn
      # cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug
      cl <- parallel::makeCluster(cn)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
    } else {
      foreach::registerDoSEQ()
    }
    
    ind <- data.frame()

    ind <- foreach(g=redes,.combine='rbind',.inorder=FALSE,.packages='igraph') %dopar% 
    {
      mot <- triad_census(g)
      mot[4] # Exploitative competition
      mot[5] # Apparent competition
      mot[6] # Tri-trophic chain
      mot[9] # Omnivory

     data.frame(explComp=mot[4],apprComp=mot[5],triTroph=mot[6],omnivory=mot[9])
    }

    # 99% confidence interval
    #
    
    qEC <- quantile(ind$explComp,c(0.005,0.995))
    qAC <- quantile(ind$apprComp,c(0.005,0.995))
    qTT <- quantile(ind$triTroph,c(0.005,0.995))
    qOM <- quantile(ind$omnivory,c(0.005,0.995))
    
    # Calculate motif for the original network
    obs <- triad_census(red)
    
    zEC <- (obs[4] - mean(ind$explComp))/sd(ind$explComp)
    zAC <- (obs[5] - mean(ind$apprComp))/sd(ind$apprComp)
    zTT <- (obs[6] - mean(ind$triTroph))/sd(ind$triTroph)
    zOM <- (obs[9] - mean(ind$omnivory))/sd(ind$omnivory)
    
    zscores <- data.frame(explComp=obs[4],apprComp=obs[5],triTroph=obs[6],omnivory=obs[9],zEC=zEC,zAC=zAC,zTT=zTT,zOM=zOM,EClow=qEC[1],EChigh=qEC[2],AClow=qAC[1],AChigh=qAC[2],TTlow=qTT[1],TThigh=qTT[2],OMlow=qOM[1],OMhigh=qOM[2])
    
    return(list(zscores=zscores,simuls=ind))         
}    

# Calculation of the clustering coefficients and average path for random network simulations
#
#
calc_modularity_random<- function(red, nsim=1000){
  
  t <- calc_topological_indices(red)
    
  redes.r <- lapply(1:nsim, function (x) {
    e <- erdos.renyi.game(t$Size, t$Links, type="gnm",directed = TRUE)
    while(components(e)$no>1)
      e <- erdos.renyi.game(t$Size, t$Links, type="gnm",directed = TRUE)
    
    return(e) }
    )

  ind <- data.frame()
  require(doParallel)
  cn <-detectCores()
#  cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug 
  cl <- makeCluster(cn)
  registerDoParallel(cl)
  
  ind <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE,.packages='igraph') %dopar% 
    {
    m<-cluster_spinglass(redes.r[[i]])
    modl <- m$modularity
    ngrp <- length(m$csize)
    clus.coef <- transitivity(redes.r[[i]], type="Global")
    cha.path  <- average.path.length(redes.r[[i]])
    data.frame(modularity=modl,ngroups=ngrp,clus.coef=clus.coef,cha.path=cha.path)
  }
  stopCluster(cl)
  ind <- ind %>% mutate(gamma=t$Clustering/clus.coef,lambda=t$PathLength/cha.path,SWness=gamma/lambda)
  # 99% confidence interval
  #
  qSW <- quantile(ind$SWness,c(0.005,0.995))
  qmo <- quantile(ind$modularity,c(0.005,0.995))
  qgr <- quantile(ind$ngroups,c(0.005,0.995))
  mcc <- mean(ind$clus.coef)
  mcp <- mean(ind$cha.path)
  mmo <- mean(ind$modularity)
  mgr <- mean(ind$ngroups)
  mSW <- mean(t$Clustering/mcc*mcp/t$PathLength)
  mCI <- 1+(qSW[2]-qSW[1])/2  
  return(data_frame(rndCC=mcc,rndCP=mcp,rndMO=mmo,rndGR=mgr,SWness=mSW,SWnessCI=mCI,MOlow=qmo[1],MOhigh=qmo[2],
                    GRlow=qgr[1],GRhigh=qgr[2]))         
}

#' Calc incoherence z-score and confidence interval under a random Erdos-Renyi 
#' networks with the condition of at least one basal node
#'
#' @param g igraph network object
#' @param ti trophic level vector
#'
#' @return
#' @export
#'
#' @examples
calc_incoherence_z <- function(g,ti=NULL,nsim=1000) {

    t <- calc_topological_indices(g)
    
    redes.r <- lapply(1:nsim, function (x) {
        e <- erdos.renyi.game(t$Size, t$Links, type="gnm",directed = TRUE)
        basal <- length(V(e)[degree(e,mode="in")==0])
        while(components(e)$no>1 | basal==0){
          e <- erdos.renyi.game(t$Size, t$Links, type="gnm",directed = TRUE)
          basal <- length(V(e)[degree(e,mode="in")==0])
        }
      return(e) }
    )
    
    ind <- data.frame()
    require(doParallel)
    cn <-detectCores()
    # #  cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug 
    cl <- makeCluster(cn)
    registerDoParallel(cl)
    
    ind <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE,.packages='igraph',.export = 'calc_incoherence') %do% 
    {
      m<-calc_incoherence(redes.r[[i]])
      data.frame(Q=m$Q,mTI=m$mTI)
    }
    stopCluster(cl)
    
    qQ <- quantile(ind$Q,c(0.005,0.995))
    qTI <- quantile(ind$mTI,c(0.005,0.995))
    rndQ <- mean(ind$Q)
    rndTI <- mean(ind$mTI)

    m <- calc_incoherence(g,ti)
    
    zQ <-  (m$Q- rndQ)/sd(ind$Q)
    zTI <- (m$mTI - rndTI)/sd(ind$mTI) # the same as sd(ind$mTI)
    #
    return(data_frame(rndQ=rndQ,rndTI=rndTI,Qlow=qQ[1],Qhigh=qQ[2],
                      TIlow=qTI[1],TIhigh=qTI[2],zQ=zQ,zTI=zTI))         
    
    
}
  

#' Calc topological roles 
#'
#' @param g an Igraph object with the network 
#' @param nsim  number of simulations with different community
#'
#' @return a  data frame with two fields: within_module_degree, among_module_conn
#' @export
#'
#' @examples
calc_topological_roles <- function(g,nsim=1000)
{
  
  toRol <- data.frame()
  require(doParallel)
  cn <-detectCores()
  #  cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug 
  cl <- makeCluster(cn)
  registerDoParallel(cl)
  
  toRol <- foreach(idx=1:nsim,.combine='rbind',.inorder=FALSE,.packages='igraph') %dopar% 
  {
    # within-module degree
    #
    # Standarized Within module degree z-score 
    #
    m<-cluster_spinglass(g)
    spingB.mem<- m$membership
    
    l<-vector()
    memMod<-vector()
    
    for (i in 1:vcount(g)){
      
      sp.in.deg <- V(g)[nei(i, "in")]
      sp.out.deg<- V(g)[nei(i, "out")]
      mem.sp<-spingB.mem[i]
      k<- length(which(spingB.mem[c(sp.in.deg, sp.out.deg)]==mem.sp))
      mem<- which(spingB.mem==mem.sp)
      
      for (m in 1:length(mem)){
        mem.in.deg <- V(g)[nei(mem[m], "in")]
        mem.out.deg<- V(g)[nei(mem[m], "out")]
        memMod.id<- length(which(spingB.mem[c(mem.in.deg, mem.out.deg)]==mem.sp))
        memMod[m]<- memMod.id
      }
      
      k.ave<- mean(memMod)
      k.sd<- sd(memMod)
      l[i]<- (k-k.ave)/k.sd
    }
    
    # among module connectivity
    r<- vector()
    for (i in 1:vcount(g)){
      
      d<-degree(g)[i]
      sp.in.deg <- V(g)[nei(i, "in")]
      sp.out.deg<- V(g)[nei(i, "out")]
      mod.sp<-table(spingB.mem[c(sp.in.deg, sp.out.deg)])
      mod.no<-as.numeric(names(mod.sp))
      mod<-rep(0, length(unique(spingB.mem)))
      mod[mod.no]<-c(mod.sp)
      r[i]<- 1-(sum((mod/d)^2))
    }
    return(data.frame(node=1:vcount(g),within_module_degree=l, among_module_conn=r))
    
  }
  stopCluster(cl)
  
  # toRol %>% group_by(node) %>% summarize(wtmLowCI=quantile(within_module_degree,0.005,na.rm=TRUE),
  #                                        wtmHiCI=quantile(within_module_degree,0.995,na.rm=TRUE),
  #                                        amcLowCI=quantile(among_module_conn,0.005,na.rm=TRUE),
  #                                        amcHiCI=quantile(among_module_conn,0.995,na.rm=TRUE),
  #                                        within_module_degree=mean(within_module_degree,na.rm=TRUE),
  #                                        among_module_conn=mean(among_module_conn,na.rm=TRUE))
  return(toRol)
}


#' Plot topological roles
#'
#' @param tRoles Calculated topological roles with the function calc_topological_roles
#' @param g Igraph network object   
#' @param spingB Igraph community object
#'
#' @return
#' @export
#'
#' @examples
plot_topological_roles <- function(tRoles,g,spingB){
  
  
  spingB.mem<- spingB$membership
  
  l <- tRoles$within_module_degree
  r <- tRoles$among_module_conn
  # Plot
  require(RColorBrewer)
  colbar.FG <- brewer.pal(length(spingB$csize),"Dark2")
  groups.B<- spingB.mem # 1=phytopl, 2=zoopl, 3=benthos, 4=fish, 5=birds, 6=mammals
  
  par(mfrow=c(1,1))
  par(oma=c(2.7,1.3,0.7,1)) # change outer figure margins
  par(mai=c(1,1,0.7,0.7)) # size of figure margins
  yran <- range(l)
  xran <- range(r)
  plot(l~r, type="n", axes=T ,tck=F, lwd=2, ann=F, cex.axis=1.2, xlim=xran, ylim=yran)
  lines(c(0.625,0.625), yran, col="grey")
  lines(xran, c(2.5, 2.5), col="grey")
  points(r, l, col=colbar.FG[groups.B], pch=20, cex=2)
  mtext(2, text="Within module degree", line=4,font=2, cex=1.2)
  mtext(1, text="Among module connectivity",line=4, font=2, cex=1.2)
  axis(1, tck=0.02, lwd=1,las=2,lty=1, labels=F, xlim=c(0,1))
  axis(2,tck=0.02, labels=FALSE)
  
  # Which are the module hubs: many links within its own module.
  #
  modhub <- which(l>2.5)
  modhub <- modhub[which(l>2.5) %in% which(r<=0.625)]
  modlbl <- unlist(vertex_attr(g,index=modhub))
  if(is.null(modlbl))
    modlbl <- modhub
  hub_conn <- data.frame()
  
  if(length(modhub)) {
    text(r[modhub],l[modhub],labels = modlbl,cex=0.7,pos=3)
    hub_conn <- data.frame(type="modhub",node=modhub,name=modlbl)  
  }
  
  #points(r[modhub], l[modhub], cex=4, col="blue", pch=20)
  
  # Which are the hub connectors: high within and between-module connectivity
  #                              and are classified super-generalists
  #
  modhub <- which(l>2.5)
  modhub <- modhub[which(l>2.5) %in% which(r>0.625)]
  modlbl <- unlist(vertex_attr(g,index=modhub))
  if(is.null(modlbl))
    modlbl <- modhub
  if(length(modhub)) {
    text(r[modhub],l[modhub],labels = modlbl,cex=0.7,pos=3)
  }
  
  #points(r[modhub], l[modhub], cex=4, col="blue", pch=20)
  if(length(modhub)){
    hub_conn <- rbind(hub_conn, data.frame(type="hubcon",node=modhub,name=modlbl))  
  }
  
  # Which are the module specialist: Few links and most of them within its own module
  #
  modhub <- which(l<=2.5)
  modhub <- modhub[which(l<=2.5) %in% which(r<=0.625)]
  modlbl <- unlist(vertex_attr(g,index=modhub))
  if(is.null(modlbl))
    modlbl <- modhub
  
  hub_conn <- rbind(hub_conn, data.frame(type="modspe",node=modhub,name=modlbl))  
  
  # Which are the module connectors: Few links and between modules
  #
  modhub <- which(l<=2.5)
  modhub <- modhub[which(l<=2.5) %in% which(r>0.625)]
  modlbl <- unlist(vertex_attr(g,index=modhub))
  if(is.null(modlbl))
    modlbl <- modhub
  
  hub_conn <- rbind(hub_conn, data.frame(type="modcon",node=modhub,name=modlbl))  
  
}


#' Calculate average topological roles doing nsimStep simulations and repeating until there is no
#' differences using an Anderson Darling test.
#'
#' @param g igraph object with the network 
#' @param nsimStep number of repeated simulations until testing for differences, the minimun number of simulations is nsimStep*2
#'
#' @return data.frame with topological roles averaged over n*nsimStep repetitions
#' @export
#'
#' @examples
calc_avg_topological_roles <- function(g, net_name,nsimStep){
  
  tR1 <- calc_topological_roles(g,nsimStep)                       # 30 simulations is enough to obtain stable topological roles
  tsim <- nsimStep
  topoRoles_mWA_temp  <- tR1 %>% group_by(node) %>% summarize(wtmLowCI=quantile(within_module_degree,0.005,na.rm=TRUE),
                                              wtmHiCI=quantile(within_module_degree,0.995,na.rm=TRUE),
                                              amcLowCI=quantile(among_module_conn,0.005,na.rm=TRUE),
                                              amcHiCI=quantile(among_module_conn,0.995,na.rm=TRUE),
                                              within_module_degree=mean(within_module_degree,na.rm=TRUE),
                                              among_module_conn=mean(among_module_conn,na.rm=TRUE))
  
  
  print(tsim)
  # Loop
  while(TRUE){
    tR1 <- bind_rows(tR1, calc_topological_roles(g,nsimStep))
    
    saveRDS(tR1,"TopoRoles_mWA_temp.rds")
    
    #tR1 <- readRDS("TopoRoles_mWA_temp.rds")
    
    tR  <- tR1 %>% group_by(node) %>% summarize(wtmLowCI=quantile(within_module_degree,0.005,na.rm=TRUE),
                                                wtmHiCI=quantile(within_module_degree,0.995,na.rm=TRUE),
                                                amcLowCI=quantile(among_module_conn,0.005,na.rm=TRUE),
                                                amcHiCI=quantile(among_module_conn,0.995,na.rm=TRUE),
                                                within_module_degree=mean(within_module_degree,na.rm=TRUE),
                                                among_module_conn=mean(among_module_conn,na.rm=TRUE))
    
    
    require(kSamples)
    t1 <- ad.test(list(tR$among_module_conn,topoRoles_mWA_temp$among_module_conn),method="simulated",nsim=1000)
    t2 <- ad.test(list(tR$within_module_degree,topoRoles_mWA_temp$within_module_degree),method="simulated",nsim=1000)
    topoRoles_mWA_temp <- tR %>% mutate(Network=net_name)
    tsim <- tsim +nsimStep
    print(tsim)
    if(t1$ad[1,4]>0.1 && t2$ad[1,4]>0.1) break()
    
  }
  
  return(topoRoles_mWA_temp)
}

#' Plot topological roles by trophic level and module
#'
#' @param netFrame dataframe with all the networks 
#' @param netName String with name of the food web to analyse
#' @param deadNodes Vector of strings with name of dead nodes to calculate trophic level
#' @param modulObj Igraph community object with the module organization of the food web
#' @param topoFrame dataframe with topological role and node index
#' @param legendPos position of the legend "topleft", "topright" or if "" no legend.
#' @param redl igraph object, if it is not null the network is taken from it    
#'
#' @return
#' @export
#'
#' @examples

plotTopoRolesByTLByMod <- function(netFrame,netName,deadNodes,modulObj,topoFrame,legendPos="",redl=NULL){
  # 
  # Local 
  #
  if(is.null(redl)){
    
    dtot1 <- as.matrix(netFrame %>% filter(Network==netName) %>% dplyr::select(Prey_name,Predator_name))
    redl <- graph_from_edgelist(dtot1, directed  = T)
    redl <- simplify(redl)
  }
  require(NetIndices)
  if(deadNodes=="")
      troph.net2<-TrophInd(get.adjacency(redl,sparse=F))
  else
      troph.net2<-TrophInd(get.adjacency(redl,sparse=F),Dead=deadNodes)
  layout.matrix.1<-matrix(
    nrow=length(V(redl)),  # Rows equal to the number of vertices
    ncol=2
  )
  
  # 
  # Add colors with topological roles to nodes 
  #
  require(RColorBrewer)
  colnet <- brewer.pal(4,"RdYlGn")
  
  hc <- topoFrame %>% mutate(type = factor(type)) %>% filter(Network==netName) %>% arrange(node) %>% mutate(col= as.numeric(type), TL=troph.net2[,1]) 
  V(redl)$color <- colnet[hc$col]
  
  # Transform y-axis coordinates
  #
  maxnew <- max(hc$TL)
  minnew <- min(hc$TL)
  maxold <- 1
  minold <- -1
  t2 <- function(x) (maxold-minold)/(maxnew -minnew)*(x - maxnew)+maxold 
  
  
  #
  # Plot modules
  #
  layout.matrix.1[,2]<-jitter(troph.net2$TL,0.4) # y-axis value based on trophic level
  layout.matrix.1[,1]<-jitter(modulObj$membership,1) # randomly assign along x-axis
  
  
  plot(redl, vertex.color=vertex_attr(redl)$cor,vertex.label=NA,
       vertex.size=log(3*igraph::degree(redl)),
       edge.width=.3,edge.arrow.size=.2, 
       edge.color=add.alpha("grey50",0.5),
       edge.curved=0.3, layout=layout.matrix.1)
  
  
  axis(side=2,at=t2(1:4),labels=1:4,  las=1, col = NA, col.ticks = 1)
  
  legstr <- levels(hc$type)
  legstr <- c("Hub conn.", "Mod. Conn.", "Mod. Hubs", "Mod. Spec.")
  if(legendPos!="")
      legend(legendPos, pch=19, col=colnet, legend= legstr)
  
}


# Add alpha to base plot colors
# 
#
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}





#' Get species names in degree/preys outdegree/predators of one topological role ("hubcon","modspe","modcon","modhub")
#'
#' @param netFrame dataframe with all the networks 
#' @param netName String with name of the food web to analyse
#' @param deadNodes Vector of strings with name of dead nodes to calculate trophic level
#' @param topoFrame dataframe with topological role and node index
#'
#' @return
#' @export
#'
#' @examples

getTopoRolesTLdegree <- function(netFrame,netName,deadNodes,topoFrame,topoType=NULL){
  # 
  # Igraph object from dataframe
  #
  dtot1 <- as.matrix(netFrame %>% filter(Network==netName) %>% dplyr::select(Prey_name,Predator_name))
  redl <- graph_from_edgelist(dtot1, directed  = T)
  redl <- simplify(redl)
  
  require(NetIndices)
  TL<-TrophInd(get.adjacency(redl,sparse=F),Dead=deadNodes)

  if(!is.null(topoType)){
    topoFrame %>% filter(type==topoType,Network==netName) %>% rowwise() %>% mutate( preys=degree(redl,node,mode=c("in")), predators= degree(redl,node,mode=c("out")), trophLevel=TL[node,1])
  } else {
    topoFrame %>% filter(Network==netName) %>% group_by(type) %>% rowwise() %>% mutate( preys=degree(redl,node,mode=c("in")), predators= degree(redl,node,mode=c("out")), trophLevel=TL[node,1])
  }
}
  
#' Title Plot net assembly model time series of S and L only the last timeW steps are ploted
#'
#' @param AA output of a net assembly model
#' @param timeW time window used
#' @param fname file name to save the plot 
#'
#' @return
#' @export
#'
#' @examples
plot_NetAssemblyModel <- function(AA,timeW,fname=NULL){
  require(viridis)
  colnet <- viridis(3)
  
  tf <- length(AA$L)
  if(tf<timeW) stop("timeW parameter must be less than the time of the simulation")
  
  dfA <- data.frame(S=AA$S[(tf-timeW):tf],L=as.numeric(AA$L[(tf-timeW):tf]),T=c((tf-timeW):tf))
  dfA$C <- dfA$L/(dfA$S*dfA$S)
  if(is.null(fname)){
    gS <- ggplot(dfA, aes(x=T,y=S)) + geom_line(colour=colnet[1]) + theme_bw() + geom_hline(yintercept = mean(dfA$S),linetype = 2,colour="grey50")
    print(gS)
    gL <- ggplot(dfA, aes(x=T,y=L)) + geom_line(colour=colnet[2]) + theme_bw() + ylab("L") + geom_hline(yintercept = mean(dfA$L),linetype = 2,colour="grey50")
    print(gL)
    gC <- ggplot(dfA, aes(x=T,y=C)) + geom_line(colour=colnet[3]) + theme_bw() + ylab("C") + geom_hline(yintercept = mean(dfA$C),linetype = 2,colour="grey50")
    print(gC)
    return(list(gS=gS,gL=gL,gC=gC))
  } else {
    require(cowplot)
    g1 <- ggplot(dfA, aes(x=T,y=S)) + geom_line() + theme_bw() + geom_hline(yintercept = mean(dfA$S),linetype=3)
    g2 <- ggplot(dfA, aes(x=T,y=L)) + geom_line() + theme_bw() + ylab("L") + geom_hline(yintercept = mean(dfA$L),linetype=3)
    g3 <- plot_grid(g1,g2,labels = c("A","B"),align = "h")
    save_plot(fname,g3,base_width=8,base_height=5,dpi=600)
  }
}


#' Title Plot net assembly model S and L average by a moving window to check if equilibrium is reached
#'
#' @param AA output of a net assembly model
#' @param timeW time window used
#' @param fname file name to save the plot 
#'
#' @return
#' @export
#'
#' @examples
plot_NetAssemblyModel_eqw <- function(AA,timeW,fname=NULL){
  
  tf <- length(AA$L)
  df <- data.frame(S=AA$S,L=as.numeric(AA$L),T=c(1:tf))
  grandS <- mean(df$S[timeW:nrow(df)])
  grandL <- mean(df$L[timeW:nrow(df)])
  
  df$gr <- rep(1:(nrow(df)/timeW), each = timeW)
  df <- df %>% group_by(gr) %>% summarise(mS=mean(S),sdS=sd(S), mL=mean(L), sdL=sd(L),time=max(T))
  
  g1 <- ggplot(df,aes(y=mS,x=time,colour=time))+ theme_bw() + geom_point() + geom_errorbar(aes(ymin=mS-sdS,ymax=mS+sdS)) + scale_color_distiller(palette = "RdYlGn",guide=FALSE)+ geom_hline(yintercept =grandS,linetype=3 )
  g2 <- ggplot(df,aes(y=mL,x=time,colour=time))+ theme_bw() + geom_point() + geom_errorbar(aes(ymin=mL-sdL,ymax=mL+sdL))+ scale_color_distiller(palette = "RdYlGn",guide=FALSE)+ geom_hline(yintercept =grandL,linetype=3 )
  
  if(is.null(fname)){
    print(g1)
    print(g2)
  } else {
    require(cowplot)
    g3 <- plot_grid(g1,g2,labels = c("A","B"),align = "h")
    save_plot(fname,g3,base_width=8,base_height=5,dpi=600)
  }
    
  return(list(df=df,g1=g1,g2=g2))
}


#' Estimation of z-scores using Meta-Web assembly model as a null 
#'
#' @param red This is the reference network as an igraph object
#' @param Adj Adyacency matrix for the meta-web
#' @param mig Migration parameter of the meta-Web assembly model
#' @param ext Exctinction parameter of the meta-Web assembly model
#' @param ti  trophic level vector 
#' @param nsim number of simulations
#'
#' @return
#' @export
#'
#' @examples
calc_modularity_metaWebAssembly<- function(red, Adj, mig,ext,nsim=1000,ti=NULL){
  
  t <- calc_topological_indices(red)
  final_time <- 500  # Final time used in simulations of the meta-web assembly
  
  ind <- data.frame()
  require(doParallel)
  cn <-detectCores()
  #  cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug 
  cl <- makeCluster(cn)
  registerDoParallel(cl)
  ind <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE,.packages=c('MetaWebAssemblyModels','igraph'), 
                 .export = c('Adj','ext','mig','final_time','calc_incoherence')) %dopar% 
  {
    AA <- metaWebNetAssembly(Adj,mig,1,ext,final_time)
    g <- graph_from_adjacency_matrix( AA$A, mode  = "directed")
    # Select only a connected subgraph graph 
    dg <- components(g)
    g <- induced_subgraph(g, which(dg$membership == which.max(dg$csize)))
    mmm<-cluster_spinglass(g)
    modl <- mmm$modularity
    ngrp <- length(mmm$csize)
    clus.coef <- transitivity(g, type="Global")
    cha.path  <- average.path.length(g)
    mmm<-calc_incoherence(g)

    data.frame(modularity=modl,ngroups=ngrp,clus.coef=clus.coef,cha.path=cha.path,Q=mmm$Q,mTI=mmm$mTI)
  }
  stopCluster(cl)
  ind <- ind %>% mutate(gamma=t$Clustering/clus.coef,lambda=t$PathLength/cha.path,SWness=gamma/lambda)
  # 99% confidence interval
  #
  qSW <- quantile(ind$SWness,c(0.005,0.995))
  qmo <- quantile(ind$modularity,c(0.005,0.995))
  qgr <- quantile(ind$ngroups,c(0.005,0.995))
  mcc <- mean(ind$clus.coef)
  mcp <- mean(ind$cha.path)
  mmo <- mean(ind$modularity)
  mgr <- mean(ind$ngroups)
  mSW <- mean(t$Clustering/mcc*mcp/t$PathLength)
  mCI <- 1+(qSW[2]-qSW[1])/2  

  qQ <- quantile(ind$Q,c(0.005,0.995))
  qTI <- quantile(ind$mTI,c(0.005,0.995))
  mdlQ <- mean(ind$Q)
  mdlTI <- mean(ind$mTI)
  
  m <- calc_incoherence(red,ti)
  
  zQ <-  (m$Q- mdlQ)/sd(ind$Q)
  zTI <- (m$mTI - mdlTI)/sd(ind$mTI) # the same as sd(ind$mTI)
  
  return(data_frame(mdlCC=mcc,mdlCP=mcp,mdlMO=mmo,mdlGR=mgr,SWness=mSW,SWnessCI=mCI,MOlow=qmo[1],MOhigh=qmo[2],
                    GRlow=qgr[1],GRhigh=qgr[2], mdlQ=mdlQ,mdlTI=mdlTI,Qlow=qQ[1],Qhigh=qQ[2],
                                                 TIlow=qTI[1],TIhigh=qTI[2],zQ=zQ,zTI=zTI))         
}




#' Calculate motif counts for observed network and CI for meta-web assembly model networks and Z-scores 
#'
#' @param red igraph network object
#' @param Adj Adyacency matrix for the meta-web
#' @param mig Migration parameter of the meta-Web assembly model
#' @param ext Exctinction parameter of the meta-Web assembly model
#' @param nsim number of simulation to calculate random networks with the same nodes and links
#'
#' @return data.frame with all the results
#' @export
#'
#' @examples
calc_motif_metaWebAssembly<- function(red, Adj, mig, ext, nsim=1000)
{
  ind <- data.frame()
  require(doParallel)
  cn <-detectCores()
  #  cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug 
  cl <- makeCluster(cn)
  registerDoParallel(cl)
  
  final_time <- 500  # Final time used in simulations of the meta-web assembly
  
  ind <- data.frame()
  require(doParallel)
  cn <-detectCores()
  #  cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug 
  cl <- makeCluster(cn)
  registerDoParallel(cl)
  ind <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE,.packages=c('MetaWebAssemblyModels','igraph'), 
                 .export = c('Adj','ext','mig','final_time')) %dopar% 
  {
    AA <- metaWebNetAssembly(Adj,mig,1,ext,final_time)
    g <- graph_from_adjacency_matrix( AA$A, mode  = "directed")
    # Select only a connected subgraph graph 
    dg <- components(g)
    g <- induced_subgraph(g, which(dg$membership == which.max(dg$csize)))
    
    mot <- triad_census(g)
    mot[4] # Exploitative competition
    mot[5] # Apparent competition
    mot[6] # Tri-trophic chain
    mot[9] # Omnivory
    
    data.frame(explComp=mot[4],apprComp=mot[5],triTroph=mot[6],omnivory=mot[9])
  }
  stopCluster(cl)
  # 99% confidence interval
  #
  
  qEC <- quantile(ind$explComp,c(0.005,0.995))
  qAC <- quantile(ind$apprComp,c(0.005,0.995))
  qTT <- quantile(ind$triTroph,c(0.005,0.995))
  qOM <- quantile(ind$omnivory,c(0.005,0.995))
  
  # Calculate motif for the original network
  obs <- triad_census(red)
  
  zEC <- (obs[4] - mean(ind$explComp))/sd(ind$explComp)
  zAC <- (obs[5] - mean(ind$apprComp))/sd(ind$apprComp)
  zTT <- (obs[6] - mean(ind$triTroph))/sd(ind$triTroph)
  zOM <- (obs[9] - mean(ind$omnivory))/sd(ind$omnivory)
  
  return(data_frame(explComp=obs[4],apprComp=obs[5],triTroph=obs[6],omnivory=obs[9],zEC=zEC,zAC=zAC,zTT=zTT,zOM=zOM,EClow=qEC[1],EChigh=qEC[2],AClow=qAC[1],AChigh=qAC[2],TTlow=qTT[1],TThigh=qTT[2],OMlow=qOM[1],OMhigh=qOM[2]))         
}    



#' Title Plot 5 simulations of net assembly model time series of S and L only the last timeW steps are ploted
#'
#' @param metaW meta-web adjacency matrix 
#' @param m     migration
#' @param q     probability of link
#' @param a     extinction
#' @param timeW time window used
#'
#' @return
#' @export
#'
#' @examples
plot_NetAssemblyModel_sims <- function(metaW,m, q, a, tf,timeW){
  require(viridis)

  if(tf<timeW) stop("timeW parameter must be less than the time of the simulation")
  
  dfA <- data.frame()
  
  for(n in 1:5){
    AA <- metaWebNetAssembly(metaW,m,q,a,tf)
    tdfA <- data.frame(S=AA$S[(tf-timeW):tf],L=as.numeric(AA$L[(tf-timeW):tf]),T=c((tf-timeW):tf))
    tdfA$C <- tdfA$L/(tdfA$S*tdfA$S)
    tdfA$sim <- n
    dfA <- bind_rows(dfA,tdfA)
  }
  gS <- ggplot(dfA, aes(x=T,y=S,colour=sim)) + geom_point() + theme_bw() + geom_hline(yintercept = mean(dfA$S),linetype = 2,colour="grey50") + scale_color_viridis(guide=FALSE)
  print(gS)
  gL <- ggplot(dfA, aes(x=T,y=L,colour=sim)) + geom_point() + theme_bw() + ylab("L") + geom_hline(yintercept = mean(dfA$L),linetype = 2,colour="grey50") + scale_color_viridis(guide=FALSE)
  print(gL)
  gC <- ggplot(dfA, aes(x=T,y=C,colour=sim)) + geom_point() + theme_bw() + ylab("C") + geom_hline(yintercept = mean(dfA$C),linetype = 2,colour="grey50") + scale_color_viridis(guide=FALSE)
  print(gC)
  return(list(gS=gS,gL=gL,gC=gC))
}




generateER <- function(ig,nsim=1000){
  if(!is_igraph(ig))
    stop("Parameter ig must be an igraph object")
  
  size <- vcount(ig)
  links <- ecount(ig)
  
  er <-lapply(1:nsim, function (x) {
    e <- sample_gnm(size, links, directed = TRUE)
    
    # Check that the ER networks has only one connected component
    #
    while(components(e)$no>1)
      e <- sample_gnm(size, links,directed = TRUE)
    
    return(e) })
  
}

calcZscore <- function(mProp,nullDist,name,plevel=0.01)
{
  
  qQ <- quantile(nullDist,c(plevel/2,1-plevel/2))
  rndQ <- mean(nullDist)
  zQ <-  (mProp - rndQ)/sd(nullDist)
  isSig <- (mProp<qQ[1] | mProp>qQ[2])
  return(data_frame(propName=name,propMean=mProp,nullMean=rndQ,lowCI=qQ[1],highCI=qQ[2],Z=zQ,isSig))
  
}

#' Calculate migration as a linear function of trophic level with mean equal to meanMigr
#'
#'
#'
#' @param adj Adjacency matrix
#' @param meanMigr mean migration probability
#'
#' @return
#' @export
#'
#' @examples
calcMigrAsTL <- function(adj,meanMigr){
  a <- TrophInd(adj)$TL
  b <- meanMigr*max(a)/mean(a)
  return(a/max(a)*b)
  
}

#' Plot a simulation of the metaWebNetAssemblyGLV function from time=901-1000
#'
#' @param glv result of one simulation
#' @param timeRange time range to plot the simulation 
#' @param meanTreshold select species with mean greater than this number
#' @param lowMeanTreshold select species with mean lower or equal than this number
#' @param type character, "l" usual time series plot with lines, "r" raster type with abundance as colors
#'             and time running vertically.   
#'
#' @return
#' @export
#'
#' @examples
plotGLVmodel <- function( glv, timeRange=NULL, meanTreshold=NULL,lowMeanTreshold=NULL,type="l" ){
  require(tidyr)
  require(ggplot2)
  require(dplyr)

  df <- as.data.frame(t(glv$STime))
  df$Time <- 1:(nrow(df))
  require(tidyr)
  df <- gather(df,key="Species",value="N", -Time)
  if(!is.null(timeRange))
    df <- filter(df, Time>=min(timeRange) & Time<=max(timeRange))
  if(!is.null(meanTreshold)){
    
    dfm <- df %>% group_by(Species) %>% summarise(N=mean(N)) %>% filter(N>meanTreshold)
    df <- inner_join(df,dfm,by="Species") %>% rename(N=N.x)
  }
  if(!is.null(lowMeanTreshold)){
    
    dfm <- df %>% group_by(Species) %>% summarise(N=mean(N)) %>% filter(N<=lowMeanTreshold)
    df <- inner_join(df,dfm,by="Species") %>% rename(N=N.x)
  }
  if(type=="l")
    ggplot(df, aes(Time,N,colour=Species)) + geom_line() +guides(colour=FALSE) + scale_colour_viridis_d() + theme_bw()
  else if(type=="r")
    ggplot(df, aes(Species,Time,fill=N)) + geom_raster() + scale_fill_viridis_c() + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 0))
}



#' Simulates the Generalized Lotka Volterra model from a metaweb varying intensity of interactions and migration
#'
#' @param adjM        GLV Lotka-Volterra adjacency matrix that acts as a meta-web
#' @param intens      vector of Intensities of interactions to simulate
#' @param migrat      vector of migration probabilities 
#' @param nsim        Number of simulations for each combination of intensity and migration values
#' @param ncores      if >0 use parallel computation with package doParallel 
#' @param fixedM      if TRUE the migration probability is equal for each species
#'                    if FALSE the migration is centered at the value but dependente on the trophic level
#' @param selfLim     numeric vector, maximal self limitation values for mutualistic,basal,predator species
#'                    this is the diagonal of the Lotka-Volterra matrix a lower value implies a greter 
#'                    carring capacity.  
#' @param eliminateLowMean numeric value, eliminate from the calculation of C and S the species with mean 
#'                         lower than this value
#'
#' @return 
#' @export
#'
#' @examples
#' 
simulateGLVIntMigr <- function(adjM, intens, migrat,nsim=10,ncores=0,fixedM=TRUE,selfLim=c(0,0.001,0),
                               intGrowRate=NULL,eliminateLowMean=4,pesca=NULL,preserveInt=FALSE){
  require(doParallel)
  require(igraph)
  require(EcoNetwork)
  require(MetaWebAssemblyModels)
  require(dplyr)
  dfSC <- tibble()
  
  stopifnot(length(selfLim)==3)
  
  if(ncores) {
    cn <-parallel::detectCores()
    if(ncores>cn)
      ncores <- cn
    # cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else {
    foreach::registerDoSEQ()
  }
  
  sdf <- as.data.frame(expand.grid(intens,migrat))
  names(sdf) <- c('intens','migrat')
  sdf <- sdf[rep(seq_len(nrow(sdf)), nsim), ]
  
  tl <- TrophInd( ifelse(adjM<0,1,0) )$TL
  
  dfSC <- foreach(i=1:nrow(sdf),.combine='rbind',.inorder=FALSE,
                  .packages=c('MetaWebAssemblyModels','EcoNetwork')) %dopar%
  {
    # Self limitation term 
    # 0 for Mutualistic (NO MUTUALISTIC SPECIES)
    # 0.001 for basal species
    # 0 for antagonistic (No self limitation for predators)
    adjS <- adjM*sdf$intens[i]
    A1 <- generateGLVparmsFromAdj(adjS,0.1,sdf$intens[i],selfLim,migrMin = 0.0,preserveInt=preserveInt)
    
    # Growth rate as parameter or random
    if(!is.null(intGrowRate))
    {
      rsign <- sign(A1$r)
      A1$r <- intGrowRate*rsign
    }
    
    # migration fixed
    #
    if(fixedM)
      A1$m <- rep(sdf$migrat[i],times=nrow(A1$interM))
    else {
      # 
      # Calculate migration probabilty using trophic level 
        b <- sdf$migrat[i]*max(tl)/mean(tl)
        A1$m <- tl/max(tl)*b
    }
    
    # initial population values 0
    #
    yini <- rep(0,times=nrow(A1$interM))

    # Add Fishing species
    #
    if(!is.null(pesca)){
      A1$interM[, pesca] <- 50*A1$interM[, pesca]
      A1$m[pesca] <- 0
      yini[pesca] <- 1/selfLim[3]
      A1$interM[pesca,pesca] <- -selfLim[3]
      A1$r[pesca] <- 1
    } 
    
    # Simulation
    #
    A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,yini,300,0.1)
    #plotGLVmodel(A2,meanTreshold=50)
    #plot_NetAssemblyModel_eqw(A2,20)
    if(eliminateLowMean==0)
      d <- (A2$STime[,300]>0)
    else
      d <- apply(A2$STime[,200:300],1,mean)>eliminateLowMean
    g <- fromGLVadjToIgraph(A1$interM,d)
    #dg <- components(g)
    #g <- induced_subgraph(g, which(dg$membership == which.max(dg$csize)))
    rS <- sum(d)
    A3 <- A2$A[d,d]
    rL <- sum(A3)
    rC <- rL/(rS*rS)
    list(S=rS,L=rL,C=rC,intensity=sdf$intens[i],migration=sdf$migrat[i],g=g,modelRes=A2,interM=A1$interM)
  }
  dfSC <- as_tibble(dfSC) %>% mutate_at(vars(S:migration), unlist) 
  
  return(dfSC)
}


#' Simulates ONE time the Generalized Lotka Volterra model from a metaweb varying intensity of interactions and migration
#'
#' @param adjM        GLV Lotka-Volterra adjacency matrix that acts as a meta-web
#' @param intens      vector of Intensities of interactions to simulate
#' @param migrat      vector of migration probabilities 
#' @param tl          Trophic level using the predation adjacency matrix e.g. `TrophInd( ifelse(adjM<0,1,0) )$TL`
#' 
#' @param fixedM      if TRUE the migration probability is equal for each species
#'                    if FALSE the migration is centered at the value but dependente on the trophic level
#' @param selfLim     numeric vector, maximal self limitation values for mutualistic,basal,predator species
#'                    this is the diagonal of the Lotka-Volterra matrix a lower value implies a greter 
#'                    carring capacity.  
#' @param eliminateLowMean numeric value, eliminate from the calculation of C and S the species with mean 
#'                         lower than this value
#'
#' @return 
#' @export
#'
#' @examples
#' 
simulateGLVIntMigr_one <- function(adjM, intens, migrat,tl=NULL,fixedM=TRUE,selfLim=c(0,0.001,0),intGrowRate=NULL,
                                   eliminateLowMean=4,pesca=NULL){
  require(igraph)
  require(EcoNetwork)
  require(MetaWebAssemblyModels)
  require(dplyr)

  stopifnot(length(selfLim)==3)
  
  if(!fixedM)
    if(is.null(tl)) {
      tl <- TrophInd( ifelse(adjM<0,1,0) )$TL
    }
  
  
  # Self limitation term c(m,b,p)
  # m = for Mutualistic 
  # b = for basal species
  # p = for antagonistic
  adjM <- adjM*intens 
  A1 <- generateGLVparmsFromAdj(adjM,0.1,intens,selfLim,migrMin = 0.0,preserveInt=TRUE)
  if(!is.null(intGrowRate))
  {
    rsign <- sign(A1$r)
    stopifnot(length(A1$r)==length(intGrowRate))
    A1$r <- intGrowRate*rsign
  }
  
  
  # migration fixed
  #
  if(fixedM)
    A1$m <- rep(migrat,times=nrow(A1$interM))
  else {
    # 
    # Calculate migration probabilty using trophic level 
    b <- migrat*max(tl)/mean(tl)
    A1$m <- tl/max(tl)*b
  }
  
  # initial population values 0
  #
  yini <- rep(0,times=nrow(A1$interM))

  if(!is.null(pesca)){
    A1$interM[, pesca] <- 50*A1$interM[, pesca]
    A1$m[pesca] <- 0
    yini[pesca] <- 1/selfLim[3]
    A1$interM[pesca,pesca] <- -selfLim[3]
    A1$r[pesca] <- 1
  } 


  # Simulation
  #
  A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,yini,300,0.1)
  #plotGLVmodel(A2,meanTreshold=50, type="r")
  #plot_NetAssemblyModel_eqw(A2,20)
  if(eliminateLowMean==0)
    d <- (A2$STime[,300]>0)
  else
    d <- apply(A2$STime[,200:300],1,mean)>eliminateLowMean

  g <- fromGLVadjToIgraph(A1$interM,d)
  #dg <- components(g)
  #g <- induced_subgraph(g, which(dg$membership == which.max(dg$csize)))
  rS <- sum(d)
  A3 <- A2$A[d,d]
  rL <- sum(A3)
  rC <- rL/(rS*rS)
  dfSC <-list(S=rS,L=rL,C=rC,intensity=intens,migration=migrat,g=g)

  return(dfSC)
}

# Estimacion de parametros usando Aproximate Bayesian computation 
# 
#
estima_ABC_GLVIntMigr <- function(adjM,dat,parMin,parMax,intGrowthRate=NULL,fixedM=FALSE,dlim=1,sim=1000,ncores=0,
                                  eliminateLowMean=4,pesca=NULL)
{
  require(doFuture)
  require(igraph)
  require(EcoNetwork)
  require(MetaWebAssemblyModels)
  require(dplyr)
  registerDoFuture()
  if(ncores) {
    cn <- future::availableCores()
    if(ncores>cn)
      ncores <- cn
    plan(multiprocess, workers=ncores)
    # cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug
    # cl <- parallel::makeCluster(ncores)
    # doParallel::registerDoParallel(cl)
    # on.exit(parallel::stopCluster(cl))
  } else {
    plan(sequential)
  }
  if(!fixedM)
    tropl <- TrophInd( ifelse(adjM<0,1,0) )$TL
  else
    tropl <- NULL
  
    
  intenInter <- runif(sim,parMin[1],parMax[1])
  migrat     <- runif(sim,parMin[2],parMax[2])
  selfLimBasal <- runif(sim,parMin[3],parMax[3])
  selfLimPred <- runif(sim,parMin[4],parMax[4])
  df <- data_frame()
  
  out <- foreach(i=1:sim,.combine='rbind',.inorder=FALSE) %dopar%
     {
       
       # Corre el modelo
       df <- simulateGLVIntMigr_one(adjM,intenInter[i],migrat[i],tl=tropl, fixedM=fixedM, selfLim = c(0,selfLimBasal[i],selfLimPred[i]),intGrowRate=intGrowthRate,eliminateLowMean,pesca)

       # selecciona la ultima parte
       cost <-   sqrt((df$S - dat$S)^2)/dat$S+sqrt((df$L-dat$L)^2)/dat$L
       
       data.frame(cost,intensity=intenInter[i],migration=migrat[i],selfLimBasal=selfLimBasal[i],selfLimPred=selfLimPred[i],S=df$S,L=df$L,C=df$C)
     }
  
  out <- out %>% filter(cost<dlim)
}

# Cost function for priceFit algorithm
#
cost_GLVIntMigr <- function(parm,adjM,data,tl,fixedM)
{
  p <- as.list(parm)
  df <- simulateGLVIntMigr_one(adjM,p$intensity ,p$migration,fixedM = fixedM,selfLim = c(0,p$selfLimBasal,p$selfLimPred),tl=tl)
  
  # sum( (df$S - data$S)^2+(df$L-data$L)^2 )
  sqrt((df$S - data$S)^2)/data$S+sqrt((df$L-data$L)^2)/data$L
}


par_pricefit <- function (par, minpar = rep(-1e+08, length(par)), maxpar = rep(1e+08, 
                                                               length(par)), func, npop = max(5 * length(par), 50), numiter = 10000, 
          centroid = 3, varleft = 1e-08, ncores=0, ...) 
{
  require(future)
  require(future.apply)
  if(ncores==0)
    plan(sequential)
  else{
    cn <- future::availableCores()
    if(ncores>cn)
      ncores <- cn
    plan(multiprocess, workers=ncores)
    
  }

  cost <- function(par) func(par, ...)
  npar <- length(par)
  tiny <- 1e-08
  varleft <- max(tiny, varleft)
  populationpar <- matrix(nrow = npop, ncol = npar, byrow = TRUE, 
                          data = minpar + runif(npar * npop) * rep((maxpar - minpar), 
                                                                   npop))
  colnames(populationpar) <- names(par)
  populationpar[1, ] <- par
  populationcost <- future_apply(populationpar, FUN = cost, MARGIN = 1)
  iworst <- which.max(populationcost)
  worstcost <- populationcost[iworst]
  iter <- 0
  while (iter < numiter & (max(populationcost) - min(populationcost)) > 
         (min(populationcost) * varleft)) {
    iter <- iter + 1
    selectpar <- sample(1:npop, size = centroid)
    mirrorpar <- sample(1:npop, size = 1)
    newpar <- colMeans(populationpar[selectpar, ])
    newpar <- 2 * newpar - populationpar[mirrorpar, ]
    newpar <- pmin(pmax(newpar, minpar), maxpar)
    newcost %<-% cost(newpar)
    if (newcost < worstcost) {
      populationcost[iworst] <- newcost
      populationpar[iworst, ] <- newpar
      iworst <- which.max(populationcost)
      worstcost <- populationcost[iworst]
    }
  }
  ibest <- which.min(populationcost)
  bestpar <- populationpar[ibest, ]
  bestcost <- populationcost[ibest]
  return(list(par = bestpar, cost = bestcost, poppar = populationpar, 
              popcost = populationcost))
}


# Estimacion de parametros usando Aproximate Bayesian computation 
# 
#
estima_ABC_metaWebNetAssembly <- function(adjM,dat,parMin,parMax,dlim,sim=1000,ncores=0)
{
  require(doFuture)
  require(igraph)
  require(EcoNetwork)
  require(MetaWebAssemblyModels)
  require(dplyr)
  registerDoFuture()
  if(ncores) {
    cn <- future::availableCores()
    if(ncores>cn)
      ncores <- cn
    plan(multiprocess, workers=ncores)
    # cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug
    # cl <- parallel::makeCluster(ncores)
    # doParallel::registerDoParallel(cl)
    # on.exit(parallel::stopCluster(cl))
  } else {
    plan(sequential)
  }

  mm <- runif(sim,parMin[1],parMax[1])
  ee <- runif(sim,parMin[2],parMax[2])
  df <- data_frame()
  tf <- 500
  
  out <- foreach(i=1:sim,.combine='rbind',.inorder=FALSE) %dopar%
  {
    
    # Corre el modelo
    AA <- metaWebNetAssembly(adjM,mm[i],ee[i],tf)
    dfA <- data.frame(S=AA$S[200:tf],L=as.numeric(AA$L[200:tf]),T=c(200:tf))
    df  <- data.frame(S=mean(dfA$S),L=mean(dfA$L),C=mean(dfA$L)/(mean(dfA$S)*mean(dfA$S)))

    # selecciona la ultima parte
    cost <-   sqrt((df$S - dat$S)^2)/dat$S+sqrt((df$L-dat$L)^2)/dat$L
    
    data.frame(cost,migrat=mm[i],extinc=ee[i],S=df$S,L=df$L,C=df$C)
  }
  
  out <- out %>% filter(cost<dlim)
}


# Cost function for priceFit algorithm
#
cost_metaWebNetAssembly <- function(parm,adjM,data,tf)
{
  p <- as.list(parm)
  
  # Corre el modelo
  AA <- metaWebNetAssembly(adjM,p$migrat,p$extin,tf)
  dfA <- data.frame(S=AA$S[200:tf],L=as.numeric(AA$L[200:tf]),T=c(200:tf))
  df  <- data.frame(S=mean(dfA$S),L=mean(dfA$L),C=mean(dfA$L)/(mean(dfA$S)*mean(dfA$S)))
  
  
  # sum( (df$S - data$S)^2+(df$L-data$L)^2 )
  sqrt((df$S - data$S)^2)/data$S+sqrt((df$L-data$L)^2)/data$L
}


# For GLV model Convert output of priceFit algorithm with cost<mCost to an ABC data.frame 
#
conv_PriceFit_ABCdataFrame <- function(pf,mCost,fixedM,tl,A,netw){
  require(doFuture)
  pp <- data.frame(cost=pf$popcost,intensity=pf$poppar[,"intensity"],migration=pf$poppar[,"migration"],selfLimBasal=pf$poppar[,"selfLimBasal"],selfLimPred=pf$poppar[,"selfLimPred"],Network=netw,fixedM=fixedM ) %>% filter(cost<mCost)
  # Add missing Parameters S,C,L
  #
  registerDoFuture()
  plan(multiprocess)
  out <- foreach(i=1:nrow(pp),.combine='rbind',.inorder=FALSE) %dopar%
  {
    df <- simulateGLVIntMigr_one(A,pp$intensity[i],pp$migration[i],tl=tl, fixedM=pp$fixedM, selfLim = c(0,pp$selfLimBasal[i],pp$selfLimPred[i]))
    data.frame(S=df$S,L=df$L,C=df$C)
  }
  pp <- pp %>% mutate(S=out$S,L=out$L,C=out$C)
}


# For metaWebAssemlby model Convert output of priceFit algorithm with cost<mCost to an ABC data.frame 
#
metaWeb_conv_PriceFit_ABCdataFrame <- function(pf,mCost,fixedM,A,netw){
  require(doFuture)
  p <- data.frame(cost=pf$popcost,migrat=pf$poppar[,"migrat"],extinc=pf$poppar[,"extinc"],Network=netw,fixedM=fixedM ) %>% filter(cost<mCost)
  # Add missing Parameters S,C,L
  #
  tf <- 500
  registerDoFuture()
  plan(multiprocess)
  out <- foreach(i=1:nrow(p),.combine='rbind',.inorder=FALSE) %dopar%
  {
    AA <- metaWebNetAssembly(A,p$migrat[i],p$extinc[i],tf)
    dfA <- data.frame(S=AA$S[200:tf],L=as.numeric(AA$L[200:tf]),T=c(200:tf))
    df  <- data.frame(S=mean(dfA$S),L=mean(dfA$L),C=mean(dfA$L)/(mean(dfA$S)*mean(dfA$S)))

  }
  p <- p %>% mutate(S=out$S,L=out$L,C=out$C)
}



#' Simulates to obtain a Time Series of the Generalized Lotka Volterra model from a metaweb varying intensity of interactions and migration
#'
#' @param adjM        GLV Lotka-Volterra adjacency matrix that acts as a meta-web
#' @param intens      vector of Intensities of interactions to simulate
#' @param migrat      vector of migration probabilities 
#' @param tl          Trophic level using the predation adjacency matrix e.g. `TrophInd( ifelse(adjM<0,1,0) )$TL`
#' 
#' @param fixedM      if TRUE the migration probability is equal for each species
#'                    if FALSE the migration is centered at the value but dependente on the trophic level
#' @param selfLim     numeric vector, maximal self limitation values for mutualistic,basal,predator species
#'                    this is the diagonal of the Lotka-Volterra matrix a lower value implies a greter 
#'                    carring capacity.  
#'
#' @param intGrowthRate vector of intrinsec growth rates, if NULL is generated randomly between 0 and 1 
#'                      with sign depending if it is a basal or not.
#' @param eliminateLowMean numeric, eliminate species with mean lower than this value
#'      
#' @return 
#' @export
#'
#' @examples
#' 
simulateGLVIntMigr_TS <- function(adjM, intens, migrat,tl=NULL,fixedM=TRUE,selfLim=c(0,0.001,0),intGrowthRate=NULL, eliminateLowMean=4,pesca=NULL,preserveInt=FALSE){
  require(igraph)
  require(EcoNetwork)
  require(MetaWebAssemblyModels)
  require(dplyr)
  
  stopifnot(length(selfLim)==3)
  if(!fixedM)
    stopifnot(!is.null(tl))
  
  # Self limitation term c(m,b,p)
  # m = for Mutualistic 
  # b = for basal species
  # p = for antagonistic
  
  adjM <- adjM*intens 
  A1 <- generateGLVparmsFromAdj(adjM,0.1,intens,selfLim,migrMin = 0.0,preserveInt=preserveInt)

  #
  # To replace A1$r we have to maintain the sign
  #
  if(!is.null(intGrowthRate))
  {
    if(length(intGrowthRate)!= length(A1$r))
      stop('The length of intGrowthRate must be the same that the number of species in adjM')
    
    rsign <- sign(A1$r)
    A1$r <- intGrowthRate*rsign
  }
  
  # migration fixed
  #
  if(fixedM)
    A1$m <- rep(migrat,times=nrow(A1$interM))
  else {
    # 
    # Calculate migration probabilty using trophic level 
    b <- migrat*max(tl)/mean(tl)
    A1$m <- tl/max(tl)*b
  }
  
  # initial population values 0
  #
  yini <- rep(0,times=nrow(A1$interM))

  # Add fishing species
  #
  if(!is.null(pesca)){
    A1$interM[, pesca] <- 50*A1$interM[, pesca]
    A1$m[pesca] <- 0
    yini[pesca] <- 1/selfLim[3]
    A1$interM[pesca,pesca] <- -selfLim[3]
    A1$r[pesca] <- 1
  }  
  # Simulation
  #
  A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,yini,300,0.1)
  #plotGLVmodel(A2,meanTreshold=4)
  #plot_NetAssemblyModel_eqw(A2,20)

  d <- apply(A2$STime[,200:300],1,mean)>eliminateLowMean
  g1 <- plotGLVmodel(A2,meanTreshold=eliminateLowMean)
  print(g1)
  if(eliminateLowMean==0)  
    g1a <- plotGLVmodel(A2,timeRange=250:300,type="r")
  else
    g1a <- plotGLVmodel(A2,timeRange=250:300,lowMeanTreshold=eliminateLowMean,type="r")
  print(g1a)
  
  #
  # Plot with a window and mean 
  #
  l1 <- plot_NetAssemblyModel_eqw(A2,20)
  
  #g <- fromGLVadjToIgraph(A1$interM,d)
  #dg <- components(g)
  #g <- induced_subgraph(g, which(dg$membership == which.max(dg$csize)))
  rS <- sum(d)
  A3 <- A2$A[d,d]
  rL <- sum(A3)
  rC <- rL/(rS*rS)
  dfSC <-list(S=rS,L=rL,C=rC,intensity=intens,migration=migrat,g1=g1,g1a=g1a,g2=l1$g1,g3=l1$g2)
  
  return(dfSC)
}


#' Calculates the distribution of abundances of species 
#'
#' @param modelRes vector of metaWebNetAssemblyGLV model results
#'
#' @return
#' @export
#'
#' @examples
calc_abund_distr <- function(modelRes,timeInterval=200:300){
  df <- lapply(seq_along(modelRes), function(i){
    A2 <- modelRes[[i]]
    d <- apply(A2$STime[,timeInterval],1,mean)
    data_frame( Dens=d,Species=1:length(d),Rep=i,Rank=dense_rank(desc(d)))
  })
  do.call(rbind,df)
}

#' Calculates the quantitative connectance using the interaction matrix and a vector of species abundances
#' 
#' The Quantitative connectance takes into account both the distribution of per capita interaction strengths
#' among species in the web and the distribution of species’ abundances and quantifies the diversity of network fluxes
#' if all the species have the same flux is equal to the directed connectance 
#' 
#' @references 
#' Fahimipour, A.K. & Hein, A.M. (2014). The dynamics of assembling food webs. Ecol. Lett., 17, 606–613
#' 
#' 5. Bersier, L.-F., Banašek-Richter, C. & Cattin, M.-F. (2002). Quantitative Descriptors Of Food-web Matrices. Ecology, 83, 2394–2407
#' 
#' 1. Ulanowicz, R.E. & Wolff, W.F. (1991). Ecosystem flow networks: Loaded dice? Math. Biosci., 103, 45–68
#' 
#' @param interM per capita interaction strength matrix  
#' @param d      species' abundances vector 
#'
#' @return 
#' @export
#'
#' @examples
calc_quantitative_connectance <- function(interM,d){
    
    nsp <- length(d) # number of species
    TT <- matrix(0,nsp,nsp)
    Hout <- numeric(nsp)
    Hin  <- numeric(nsp)
    
    for(i in seq_along(d)) {
      for(j in seq_along(d)) {
        TT[i,j] <- abs(interM[i,j])*d[i]*d[j]
      }
    }
    TTcol <- colSums(TT)    # Outflow
    TTrow <- rowSums(TT)    # Inflow
    
    for(k in seq_along(d)) {
      for(i in seq_along(d)) {
        if(TT[k,i]>0){
          Hin[k] <- Hin[k] - TT[k,i]/TTrow[k]*log(TT[k,i]/TTrow[k])
        }
        if(TT[i,k]>0){
          Hout[k] <- Hout[k] - TT[i,k]/TTcol[k]*log(TT[i,k]/TTcol[k])
        }
      }
    }
    nout <- exp(Hout)
    nout[TTrow==0] <- 0
    nin <- exp(Hout)
    nin[TTcol==0] <- 0

    Cq <- 1/(2*nsp^2)*(sum(exp(Hout))+sum(exp(Hin)))
}


#' Calculates cuantitative connectance for a dataframe with results of simulations done  
#' [meweasmo::metaWebNetAssemblyGLV()] obtained with the function simulateGLVIntMigr() 
#'
#' @param mdl data.frame with simulations output from simulateGLVIntMigr()
#'
#' @return data.frame with tS= total number of species, tL= number of links, tC= directed connectance,
#'         Cq=quantitative connectance 
#' @export
#'
#' @examples
calc_Cq_fromGLVsims <- function(mdl,ncores=0){
  require("future.apply")
  
  if(ncores) {
    cn <- future::availableCores()
    if(ncores>cn)
      ncores <- cn
    plan(multiprocess, workers=ncores)
  } else {
    plan(sequential)
  }
  
  df <- future_lapply(seq_len(nrow(mdl)), function(i){
    tA <- mdl$interM[[i]]
    #td <- apply(mdl$modelRes[[i]]$STime[,200:300],1,mean)
    finalT <- ncol(mdl$modelRes[[i]]$STime)
    td <- mdl$modelRes[[i]]$STime[,finalT]
    Cq <- calc_quantitative_connectance(tA,td)
    rS <- sum(td>0)
    rL <- sum(mdl$modelRes[[i]]$A)
    C <- rL/(rS*rS)
    data_frame( tS=rS,tL=rL,tC=C,Cq=Cq )
  })
  do.call(rbind,df)

}
    
#' Calculates variability for a dataframe with results of 
#' metaWebNetAssemblyGLV() obtained with the function simulateGLVIntMigr() 
#'
#' @param mdl data.frame with simulations output from simulateGLVIntMigr()
#'
#' @return data.frame with tS= total number of species, tL= number of links, tC= directed connectance,
#'         Cq=quantitative connectance 
#' @export
#'
#' @examples
calc_Variability_fromGLVsims <- function(mdl,ncores=0){
  
  require("future.apply")
  
  if(ncores) {
    cn <- future::availableCores()
    if(ncores>cn)
      ncores <- cn
    plan(multiprocess, workers=ncores)
  } else {
    plan(sequential)
  }
  
  df <- future_lapply(seq_len(nrow(mdl)), function(i){
    tA <- mdl$interM[[i]]
    #td <- apply(mdl$modelRes[[i]]$STime[,200:300],1,mean)
    finalT <- ncol(mdl$modelRes[[i]]$STime)
    initT <- finalT - 100
    td <- mdl$modelRes[[i]]$STime[,initT:finalT]
    totalN <- colSums(td)
    totalV <- sd(totalN)/mean(totalN)    
    totalSpV <- sd(mdl$modelRes[[i]]$S[initT:finalT])/nrow(td)
    data_frame( VDens=totalV,VSpecies=totalSpV )
  })
  do.call(rbind,df)
  
}



calc_modularity <- function(ig,ncores=0){
  
  if(inherits(ig,"igraph")) {
    ig <- list(ig)
  } else if(class(ig[[1]])!="igraph") {
    stop("parameter ig must be an igraph object")
  }
  
  registerDoFuture()
  if(ncores) {
    cn <- future::availableCores()
    if(ncores>cn)
      ncores <- cn
    future::plan(multiprocess, workers=ncores)
    on.exit(future::plan(sequential))
  } else {
    future::plan(sequential)
  }
  df <-  foreach(g=ig,.combine='rbind',.inorder=FALSE,.packages=c('igraph','NetIndices')) %dopar%
  {
    m<-cluster_spinglass(g,weights=NA)
    modl <- m$modularity
    data.frame(modularity=modl)
  }
return(df)  
}