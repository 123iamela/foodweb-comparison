######################################################################
### Stability of ecosystems: Lotka-Volterra Multispecies model
######################################################################
 
###################################
# FUNCTION DEFINITIONS
###################################
 
###
# lvm(t,x,parms)
# Use:    Function to calculate derivative of multispecies Lotka-Volterra equations
# Input:
#     t: time (not used here, because there is no explicit time dependence)
#     x: vector containing the current abundance of all species
#     parms:     dummy variable, which is not used here (normally used to pass on
#        parameter values, but not needed here because a and r are defined globally)
# Output:
#    dx: derivative of Lotka-Volterra equations at point x
lvm <- function(t,x,parms){
dx <- (r - a %*% x) * x
list(dx)
}
 
###
# n.integrate(time,init.x,model)
# Use:    Numerical integration of model
# Input:
#     t:     list with elements time$start, time$end, and time$steps, giving start and
#        endpoint of integration and the number of time points in between
#     init.x: vector containing initial values (at time = time$start) of all species
#     model:     name of the function to integrate (here lvm)
# Output:
#    data frame with n+1 columns. The first column contains the time points at which
#    x is evaluated. The next n columns are the values of the n species at these
#    time points
# Description:
# Generates a vector of time points for the integration and uses function lsoda (from library
# odesolve) to integrate the model
n.integrate <- function(time=time,init.x= init.x,model=model){
t.out <- seq(time$start,time$end,length=time$steps)
as.data.frame(lsoda(init.x,t.out,model,parms=parms))
}
 
###
# cutoff(out,tol)
# Use: Get index of species (i) that are extinct (i.e. whose frequency is smaller than tol)
#      (ii) species that survive (i.e. whose freq is larger than tol) and (iii) return a vector
#      of the abundance of species, in which those with freq < tol are set to 0.
# Input:
#    out:     This is the output of n.integrate (i.e. a data frame with a column for the
#        time points and the abundance of the n species)
#    tol:    The threshold value below which a species is considered to be extinct
# Output:
#    list with three entries
#    (i) indices of species that are extinct
#       (ii) indices of species that survive
#       (iii) vector of species abundances where those with freq < tol are set to 0. (This vector
#    can be used as the new initial conditions to continue the integration of the model
# Description:
#    Takes output of n.integrate and determines which species
#    are more abundant than tol at the last time point.
 
cutoff <- function(out,tol){
n.end <- length(out[,1])
frequency <- out[n.end,2:(n+1)] #Note that species abundance is in columns 2 to n+1
extinct <- which(frequency<tol)
survive <- which(frequency>=tol)
frequency[extinct] <- 0
list(extinct=extinct,survive=survive,freq=as.numeric(frequency))
}
 
###
# generate.parameters(mat,rep,index,sparse)
# Use:      Gets reproduction rate and interaction matrix and generates new reproduction rates and interaction
#     coefficients of species with given set of indices
# Input:
#    mat:    current matrix of interaction coefficients
#    rep:    current reproduction rates of the species
#    index:     vector of indices of species for which new parameters have to be created (use for
#        example output of cutoff() as indices for extinct species)
#    sparse:    Sparsity of interaction matrix. The parameters sparse controls how many
#        interactions are set to 0. sparse has to between 0 and 1, where sparse = 1
#        corresponds to full connectivity.
# Output:
#    list containing updated reproduction rates and interaction matrix
generate.parameters <- function(mat,rep,index,sparse){
length.index <- length(index)
# generate random reproduction rates for species with indices in index vector
rep[index]=runif(length.index)
# generate uniformly distributed random interaction coefficients to replace rows of
# interaction matrix with corresponding indices
m1 <- matrix(runif(n*length.index),ncol=n)
# generate matrix of 0's and 1's as a function of sparse (using rbinom())
m2 <- matrix(rbinom(n*length.index,1,sparse),ncol=n)
# replace the corresponding rows in the interaction matrix
mat[index,] <- m1*m2
# Do the same for the columns
mat[,index] <- matrix(runif(n*length.index),nrow=n)*matrix(rbinom(n*length.index,1,sparse),nrow=n)
# make sure that diagonal elements are not zero (i.e. regenerate random coefficients for diag elements)
diag(mat)[index] <- runif(length.index)
list(matrix=mat,reproduction=rep)
}
 
### Plotting routines
###
# plot.lvm.time(out)
# Use:    Plots all species against time
# Input:
#    out: output of n.integrate
plot.lvm.time <- function(out){
n <- length(out[1,])-1
t.range <- range(out[,1])
y.range <- range(out[,2:(n+1)])
plot(t.range,y.range,xlab="time",ylab="abundance",type="n")
for (i in c(2:(n+1))) {lines(out[,1],out[,i],col=i)}
}
 
###
# plot.matrix(matrix,index)
# Use:    Plot interaction matrix (setting all species in index to zero)
# Input:
#    matrix: interaction matrix
#    index:    indices of species that are set to 0. (These would be the species
#         indices determined by function cutoff() that are set to zero because
#        the species abundance is smaller than tol).
plot.matrix <- function(matrix,index){
m <- matrix
m[index,] <-0
m[,index] <-0
image(t(m),axes=F,zlim=c(0,1),col=gray(c(32:0)/32),xlab="species affected",ylab="species affecting")
axis(1, at = seq(0,1,length=n),labels = 1:n)
axis(2, at = seq(0,1,length=n),labels = 1:n)
}
 
###
# plot.frequency(out)
# Use: Plots abundance of species at the end of simulation
# Input:
#     out: output of n.integrate()
# Plot frequency at end of simulation
plot.frequency <- function(out) {
plot(c(1:n), out[length(out[,1]),2:(n+1)],type="h",xlab="species index",ylab="abundance")
}
 
###########################
# MAIN PROGRAM
###########################
 
### LOAD LIBRARIES
#load R library for ordinary differential equation solvers
library(deSolve)
 
### INITIALIZE PARAMETER SETTINGS
# Integration window
time<- list(start=0,end=30,steps=100)
# dummy variable for lvm() function defined above
parms <- c(0) ### dummy variable (can have any numerical value)
# Number of Species
n<-10
# Sparsity of interaction matrix (0<sparse<1) [sparse = 1 ==> full connectivity]
sparse <- 0.2
# Generate parameters
# first generate reproduction vector and interaction matrix containing only zeros.
r <- numeric(n)
a <- matrix(0,nrow=n,ncol=n)
# Now generate random numbers for parameters according to function generate.parameters
tmp <- generate.parameters(a,r,c(1:n),sparse)
r <- tmp$reproduction
a <- tmp$matrix
# Initial conditions (i.e. frequencies of n species at time = time$start)
init.x <- runif(n)
# Threshold at which species are assumed to be extinct
tol <- 0.0001
 
### Example for a loop that integrates the LVM and allows new species to invade
outt <- init.x
for (i in c(1:20)){
# integrate lvm model
out <-n.integrate(time=time,init.x=init.x,model=lvm)
# plot time course
plot.lvm.time(out)
# save out and attach it to outt
outt <- rbind(outt,out)
# determine which species went extinct and which survived
tmp <- cutoff(out,tol)
# save index of extinct species (i.e. species with frequency smaller than tol)
index.extinct <- tmp$extinct
# generate new parameters (reproduction rates and interaction coefficients) to replace extinct species
tmp2 <- generate.parameters(a,r,index.extinct,sparse)
a <- tmp2$matrix
r <- tmp2$reproduction
# generate new init.x to continue integration
init.x <- tmp$freq
# introduce new species with a frequency of 10* tol
init.x[index.extinct] <- 10* tol
# generate new time window to continue integration
time <- list(start=time$end, end = time$end+time$end-time$start,
steps=100)
}
# plot entire time course
plot.lvm.time(outt)