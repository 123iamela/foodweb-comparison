
# Comulative degree distributions
# Degree distribution: proportion of links per node as a function of links
# The degree distribution of a food web can tell us a lot about the amount of specialization and 
# generalization in the web (in degree), as well as vulnerability (out degree)
 
deg.B<-degree(B.FW.edge.list)
deg.A<-degree(Ab.FW.edge.list)
deg.Aa<-degree(Aa.FW.edge.list)
 
# Using the degree distribution gives a better way to visualize any differences
# Looking at the in degree tells us about how general the diets of consumers are
dd.B.in<-degree.distribution(B.FW.edge.list,mode="in",cumulative=T)
dd.A.in<-degree.distribution(Ab.FW.edge.list,mode="in",cumulative=T)
dd.Aa.in<-degree.distribution(Aa.FW.edge.list,mode="in",cumulative=T) 

# Out degree is a measure of the vulnerability of organisms, telling us how many consumers eat each species.
dd.B.out<-degree.distribution(B.FW.edge.list,mode="out",cumulative=T)
dd.A.out<-degree.distribution(Ab.FW.edge.list,mode="out",cumulative=T)
dd.Aa.out<-degree.distribution(Aa.FW.edge.list,mode="out",cumulative=T)


par(mfrow=c(1,1))
par(oma=c(2.7,1.3,0.7,1)) # change outer figure margins
par(mai=c(1,1,0.7,0.7)) # size of figure margins
op <- par(family = "serif")
par(op)

plot(dd.B.in,log="y", las=1, type="n",  axes=T ,tck=F, cex.lab=1, ann=F, cex.axis=1)
points(dd.B.in,log="y", xlab="Links", col="darkgrey", pch=16)
points(dd.A.in,log="y", xlab="Links", col="black", pch=16)
# points(dd.Aa.in,log="y", xlab="Links", col="black", pch=1)
#mtext(2, text="Log cumulative degree", line=4.5,font=2, cex=1)
#mtext(1, text="In-degree",line=3, font=2, cex=1)
axis(1,  lwd=1,las=2,lty=1, labels=F, cex.axis=1.5, cex=1.2)
axis(2, labels=FALSE, cex.axis=1.5)
legend("topright", legend=c("Boreal", "Arctic"), pch=c(16,16), col=c("darkgrey", "black"), bty="n", cex=1.2)
par(op)
# "Arctic II"
# tck=0.02 # tick marks inside

plot(dd.B.out,log="y", type="n", las=1, axes=T ,tck=F, cex.lab=1, ann=F, cex.axis=1, xlim=c(0,80))
points(dd.B.out,log="y", xlab="Links", col="darkgrey", pch=16)
points(dd.A.out,log="y", xlab="Links", col="black", pch=16)
# points(dd.Aa.out,log="y", xlab="Links", col="black", pch=1)
#mtext(2, text="Log cumulative degree", line=4.5,font=2, cex=1)
#mtext(1, text="Out-degree",line=3, font=2, cex=1)
axis(1,  lwd=1,las=2,lty=1, labels=F, cex.axis=1,xlim=c(0,80))
axis(2,labels=FALSE, cex.axis=1.2)
legend("topright", legend=c("Boreal", "Arctic"), pch=c(16,16), col=c("darkgrey", "black"), bty="n", cex=1.2)

legend 
