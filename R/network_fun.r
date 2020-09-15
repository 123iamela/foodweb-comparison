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
