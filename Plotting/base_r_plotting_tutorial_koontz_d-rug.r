### Title: Back to basics: High quality plots using base R graphics
###	An interactive tutorial for the Davis R Users Group meeting on April 24, 2015
###
### Date created: 20150418
### Last updated: 20150423
###
### Author: Michael Koontz
### Email: mikoontz@gmail.com
###	Twitter: @michaeljkoontz
###
### Purpose: Introduce basic to intermediate plotting capabilities of base R graphics
###
### Basic methods
###		1) Basic scatterplot and labeling a plot (Line 44)
###		2) Plotting groups in different ways on the same plot (Line 72)
###		3) Adding a legend (Line 120)
###		4) Adding a best fit line (Line 150)
###		5) Adding a 95% confidence interval (Line 150)
###		6) Shaded confidence intervals (Line 223)
###		7) Bar plots (Line 260)
###		8) Error bars (Line 274)
###
### Intermediate methods
###   	1) Using other graphics devices like pdf() (Line 324)
###   	2) Using par() for multipanel plots (Line 380) 
###		3) Using par() for margin adjustments (Line 438)
###		4) Using axis() and mtext() (Line 484)
###		5) Pretty print from plotmath (Line 601)


# We'll start with the very tractable 'trees' dataset, which is built into R. It describes the girth, height, and volume of 31 black cherry trees.

trees
dim(trees)
head(trees)

# Remember how we access the columns of a data.frame:
trees$Girth
trees$Volume

#Basic plot function takes x argument and a y argument
#Default plot type is points, but you can change it to lines or both points and lines by adding the 'type' argument

plot(x=trees$Girth, y=trees$Volume)
plot(x=trees$Girth, y=trees$Volume, type="l")
plot(x=trees$Girth, y=trees$Volume, type="b")

# pch: 'plotting character' changes the type of point that is used (default is an open circle); remember pch=19!
plot(x=trees$Girth, y=trees$Volume, pch=19)

# main: adds a title
plot(x=trees$Girth, y=trees$Volume, pch=19, main="Girth vs. Volume for Black Cherry Trees")

# xlab: adds an x axis label
plot(x=trees$Girth, y=trees$Volume, pch=19, main="Girth vs. Volume for Black Cherry Trees", xlab="Tree Girth (in)")

# ylab: adds a y axis label
plot(x=trees$Girth, y=trees$Volume, pch=19, main="Girth vs. Volume for Black Cherry Trees", xlab="Tree Girth (in)", ylab="Tree Volume (cu ft)")

# las: rotates axis labels; las=1 makes them all parallel to reading direction
plot(x=trees$Girth, y=trees$Volume, pch=19, main="Girth vs. Volume for Black Cherry Trees", xlab="Tree Girth (in)", ylab="Tree Volume (cu ft)", las=1)

# col: select a color for the plotting characters
plot(x=trees$Girth, y=trees$Volume, pch=19, main="Girth vs. Volume for Black Cherry Trees", xlab="Tree Girth (in)", ylab="Tree Volume (cu ft)", las=1, col="blue")

# We can use the c() function to make a vector and have several colors, plotting characters, etc. per plot.

plot(x=trees$Girth, y=trees$Volume, pch=19, main="Girth vs. Volume for Black Cherry Trees", xlab="Tree Girth (in)", ylab="Tree Volume (cu ft)", las=1, col=c("black", "blue"))

plot(x=trees$Girth, y=trees$Volume, pch=c(1,19), main="Girth vs. Volume for Black Cherry Trees", xlab="Tree Girth (in)", ylab="Tree Volume (cu ft)", las=1, col="blue")

#------------
# Plotting by group
#------------

### Those different colors and plotting characters that we just saw were arbitrary. The 2-element vector of colors or plotting characters just repeats for the whole data frame. What if we want to have more meaningful coloration, with a different color for each group? 

### We'll use the iris dataset to illustrate one way to do this. This dataframe describes the sepal length, sepal width, petal length, petal width, and species for 150 different irises.

# First look at the data:
iris
head(iris)
dim(iris)
str(iris)

# We can extend the idea of passing a vector of colors to the col= argument in the plot() function call.

# Let's cheat first, and see what the finished product will look like. First I define a new object with the three colors that I want to use.
plot.colors <- c("violet", "purple", "blue")
plot.colors

# Here's the cheating bit: I just looked at this dataframe and saw that there are exactly 50 observations for each species.
iris

# I use the repeat function, rep() and the each= argument, to create a new vector with each element of plot.colors repeated 50 times in turn.
color.vector <- rep(x=plot.colors, each=50)
color.vector

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, col=color.vector)

# Notice the lengths of the x-vector, the y-vector, and the color vector are all the same.
length(iris$Petal.Length)
length(iris$Sepal.Length)
length(color.vector)

# What if we want to automate the process? We can take advantage of the fact that the Species column is a factor. 
head(iris)
iris$Species
str(iris)
as.numeric(iris$Species)

plot.colors <- c("violet", "purple", "blue")

color.vector <- plot.colors[iris$Species]

dev.off() # Just clearing the present plots

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, col=color.vector, main="Iris sepal length vs. petal length", xlab="Petal length", ylab="Sepal length", las=1)

#-----------
# Let's add a legend
#-----------

# We use the legend() function to add a legend to an existing plot

legend("topleft", pch=19, col=plot.colors, legend=unique(iris$Species))

# You can customize the legend if you wish.
plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, col=color.vector, main="Iris sepal length vs. petal length", xlab="Petal length", ylab="Sepal length", las=1)

# Here I pass a character vector to the legend= argument so that I can include the first letter of the species name
# The bty="n" argument suppresses the border around the legend. (A personal preference)
legend("topleft", pch=19, col=plot.colors, legend=c("I. setosa", "I. versicolor", "I. virginica"), bty="n")

# Italicize the labels in the legend using text.font=3
plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, col=color.vector, main="Iris sepal length vs. petal length", xlab="Petal length", ylab="Sepal length", las=1)

legend("topleft", pch=19, col=plot.colors, legend=c("I. setosa", "I. versicolor", "I. virginica"), bty="n", text.font=3)


#------------------
# Test yourself. 
#------------------

#Using the ToothGrowth dataset built into R, plot the tooth length (the len column) as a function of the vitamin C dosage (the dose column). Use a different color for each method of administering the vitamin C (the supp column).

head(ToothGrowth)


#------------------
# Add a linear best fit line and confidence interval to a plot
#------------------

# We'll use a simple linear regression for this, but the general recipe is the same every time.
# The Recipe
#	1) Estimate the parameters of the best fit line
#	2) Make up your own x values that span the range of your data
#	3) Get your y values by applying your mathematical model (e.g. a straight line) with the best fit parameters to your fabricated x values
#	4) Plot these new y values against your fabricated x values.

# Save your model fit to an object. Here, we model Sepal.Length as a function of Petal.Length
model1 <- lm(Sepal.Length ~ Petal.Length, data=iris)

# Now we have the parameter estimates for our y=ax+b line. The estimate for (Intercept) is b, and the estimate for Petal.Length is a.
summary(model1)

# Make up our own x values; put them in a dataframe!
range(iris$Petal.Length)
xvals <- seq(from=1, to=7, by=0.1)
xvals
df <- data.frame(Petal.Length=xvals)
df

# Take advantage of the predict() function, which returns a matrix with one row for each of your new x values, and 3 columns: 'fit' is the expected y value, 'lwr' is the lower 95% confidence interval, and 'upr' is the upper 95% confidence interval. Note this only works this seamlessly using the predict() function on an lm() model

CI <- predict(model1, newdata=df, interval="confidence")
CI <- as.data.frame(CI) # Coerce the matrix to a dataframe, so we can access the column using the $ operator.

dim(df) # 61 x-values
head(df) # Here are the first 6 of them 
dim(CI) # 61 records, 3 columns
head(CI) # y-values, lower 95% CI bounds, and upper 95% CI bounds


# Plot the actual data (in the iris dataframe)
# This code copied from above
plot.colors <- c("violet", "purple", "blue") 
color.vector <- plot.colors[iris$Species]

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="Iris sepal length vs. petal length", xlab="Petal length", ylab="Sepal length", col=color.vector)

# Plot our best fit line. The x values are the Petal.Length column from the 'df' dataframe, and the y values are the 'fit' column from the CI dataframe. 
# Note that I use the lines() function, which just adds features to an existing plot.
# The lwd= argument changes the line width
lines(x=df$Petal.Length, y=CI$fit, lwd=2)

# Plot the confidence intervals
# The lty= argument changes the line type. There are 6 different line types, and you can just put a number 1 through 6 if you'd like. Default is "solid" (aka 1)
lines(x=df$Petal.Length, y=CI$lwr, lwd=2, lty="dashed", col="red")
lines(x=df$Petal.Length, y=CI$upr, lwd=2, lty="dashed", col="red")

#----------------------
# Going further... Prediction/forecast intervals
#----------------------

forecast <- predict(model1, newdata=df, interval="prediction")
forecast <- as.data.frame(forecast)
head(forecast)

# Plot the actual data (in the iris dataframe)
# This code copied from above
plot.colors <- c("violet", "purple", "blue") 
color.vector <- plot.colors[iris$Species]

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="Iris sepal length vs. petal length", xlab="Petal length", ylab="Sepal length", col=color.vector)

# New code using the forecast object
lines(x=df$Petal.Length, y=forecast$fit, lwd=2)
lines(x=df$Petal.Length, y=forecast$lwr, lwd=2, col="red", lty="dashed")
lines(x=df$Petal.Length, y=forecast$upr, lwd=2, col="red", lty="dashed")


#-------------
# Shaded region of forecast interval
#-------------

# This plot() call is copied from above
plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="Iris sepal length vs. petal length", xlab="Petal length", ylab="Sepal length", col=color.vector)

lines(x=df$Petal.Length, y=forecast$fit, lwd=2)

# New code
polygon.x <- c(df$Petal.Length, rev(df$Petal.Length))
polygon.y <- c(forecast$lwr, rev(forecast$upr))

# By default, R won't fill in the polygon
polygon(x=polygon.x, y=polygon.y)

# But we also may not want an opaque polygon
polygon(x=polygon.x, y=polygon.y, col='darkgrey')

# The adjustcolor() function is nice for 'turning down' the opacity. It takes a color and the opacity level as arguments. I also use border=NA to suppress the border of the polygon

# First recreate the plot
plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="Iris sepal length vs. petal length", xlab="Petal length", ylab="Sepal length", col=color.vector)

lines(x=df$Petal.Length, y=forecast$fit, lwd=2)

polygon(x=polygon.x, y=polygon.y, col=adjustcolor("black", alpha.f=0.4), border=NA)

# Add our legend back in
legend("topleft", pch=19, col=plot.colors, legend=c("I. setosa", "I. versicolor", "I. virginica"), bty="n", text.font=3)

#--------------
# Test yourself.
#--------------
# Layer the 95% confidence interval as a shaded region on top of the iris data, the best fit line for a Sepal.Length~Petal.Length model, and the forecast interval.


#------------
# Barplots
#------------

model2 <- lm(Sepal.Length ~ Species, data=iris)

bar.heights <- predict(model2, newdata=data.frame(Species=c("setosa", "versicolor", "virginica")))

# The basic barplot function
barplot(bar.heights)

# Let's add some flair
barplot(bar.heights, names.arg=c("I. setosa", "I. versicolor", "I. virginica"), las=1, col=adjustcolor(plot.colors, alpha.f=0.5), main="Sepal length for 3 Irises", ylab="Sepal length (cm)")

#---------------
# Error bars
#---------------
# Adding error bars to our barplot. These can be added to scatter plots in a similar way.

# We'll plot error bars representing 5 standard errors so you can see them more easily.

d <- summary(model2)
CI <- 5 * coef(d)[ ,'Std. Error']
lwr <- bar.heights - CI
upr <- bar.heights + CI

# I used the ylim= argument to pass a 2-element numeric vector specifying the y extent of the barplot. I added some extra room on the top to account for error bars.
# Importantly, assign the barplot to an object. I called it 'b' but you can call it whatever you like.

b <- barplot(bar.heights, names.arg=c("I. setosa", "I. versicolor", "I. virginica"), las=1, ylim=c(0,7.5), col=adjustcolor(plot.colors, alpha.f=0.5), main="Sepal length for 3 Irises", ylab="Sepal length (cm)")

# The object that you called your barplot is interpretted by R as the x values in the middle of each bar
b

# We'll use the arrows() function to add arrows to an existing plot. With some modifications, our arrows will have an arrowhead at each end (code=3), and the 'arrowhead' will actually be perpendicular to the arrow shaft (angle=90)
# Specify where each arrow starts (x0= and y=) and ends (x1= and y1=)
arrows(x0=b, x1=b, y0=lwr, y1=upr, code=3, angle=90, length=0.1)


#---------------
# Test yourself
#---------------
# These data represent survivorship of plant seedlings in 4 different treatments: ambient, watered, heated + watered, and heated. Make a bar plot with their 95% confidence intervals. Note these are asymmetric (more uncertainty above the mean than below) like what might come from a logistic regression model.

prop <- c(0.18, 0.25, 0.13, 0.05)
asympLCL <- c(0.14, 0.20, 0.11, 0.035)
asympUCL <- c(0.24, 0.33, 0.18, 0.09)

#---------------
# Test yourself. Error bars on scatter plots.
#---------------
# The randomly generated data below are measurements of the number of the number of angels who get their wings as a function of the number of bells that have been rung. There is some uncertainty in measuring wing acquisition (represented as the offset from the sampled mean). How would you add error bars to a scatter plot?

set.seed(13)
n <- 20 # Number of experimental trials
a <- 12
b <- 1.5
xvals <- round(runif(n)*50)
yvals <- round(a + b*xvals + rnorm(n, sd=5))
offset <- rpois(n, lambda=10)
lwr <- yvals - offset
upr <- yvals + offset


#-----------------
# Base R plotting skills: Other graphics devices
#-----------------

# 1) Use the pdf() graphics device (or png() or postscript()) to get a permanent record of your plot in the exact final format.

# First set your working directory. A working directory is where R assumes all input and output files that it interfaces with should be. When you send your code to a collaborator, they can just change the working directory once instead of changing the filepath for every input (e.g. reading data into R) or output (e.g. making a plot)

# This is how I set mine:
#	1) Make sure your script file is saved in a folder particular for the given project 
#	2) Run the file.choose() function
#	3) Navigate to your script file and open it
#	4) Copy the file path in the console up until the actual file name (i.e. just the folder path)
#	5) Paste that folder into the setwd() function in quotes
#	6) Delete the file.choose() command

file.choose()

pdf("iris plot.pdf")
plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, xlab="Petal length", ylab="Sepal length", col=plot.colors[iris$Species])

lines(x=df$Petal.Length, y=forecast$fit, lwd=2)

polygon.x <- c(df$Petal.Length, rev(df$Petal.Length))
polygon.y <- c(forecast$lwr, rev(forecast$upr))

polygon(x=polygon.x, y=polygon.y, col=adjustcolor("black", alpha.f=0.4), border=NA)

# Add our legend back in
legend("topleft", pch=19, col=plot.colors, legend=c("I. setosa", "I. versicolor", "I. virginica"), bty="n", text.font=3)

# Importantly, turn the device off at the end of your plotting block to complete the .pdf file.
dev.off()


#--------------
# Figure for a paper
#--------------
# To make a figure for a paper that meets particular size and style guidelines, add some more arguments to the pdf() function call.

pdf("iris plot for paper.pdf", width=3.3, height=3.3, pointsize=9, family="Times")

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, xlab="Petal length", ylab="Sepal length", col=plot.colors[iris$Species])
lines(x=df$Petal.Length, y=forecast$fit, lwd=2)

polygon.x <- c(df$Petal.Length, rev(df$Petal.Length))
polygon.y <- c(forecast$lwr, rev(forecast$upr))

polygon(x=polygon.x, y=polygon.y, col=adjustcolor("black", alpha.f=0.4), border=NA)

# Add our legend back in
legend("topleft", pch=19, col=plot.colors, legend=c("I. setosa", "I. versicolor", "I. virginica"), bty="n", text.font=3)

# Importantly, turn the device off at the end of your plotting block to complete the .pdf file.
dev.off()

#-----------------
# Base R plotting skills: Managing par()
#-----------------

# 2) You can change some fundamental plotting parameters by using par() before a plot.

# This is how you can:
#		Make multipanel plots
#		Give your plots more room in the inner margins
#		Give your plots some room in the outer margins

#---------------
# Multi-panel plots
#---------------

setosa <- subset(iris, subset=Species=="setosa")
versicolor <- subset(iris, subset=Species=="versicolor")
virginica <- subset(iris, subset=Species=="virginica")

# Using par(mfrow=c(number of rows, number of columns))
par(mfrow=c(2,2))

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="All Irises", ylab="Sepal Length", xlab="Petal Length")

plot(x=setosa$Petal.Length, y=setosa$Sepal.Length, pch=19, las=1, col=plot.colors[1], main="I. setosa", ylab="Sepal Length", xlab="Petal Length")

plot(x=versicolor$Petal.Length, y=versicolor$Sepal.Length, pch=19, las=1, col=plot.colors[2], main="I. versicolor", ylab="Sepal Length", xlab="Petal Length")

plot(x=virginica$Petal.Length, y=virginica$Sepal.Length, pch=19, las=1, col=plot.colors[3], main="I. virginica", ylab="Sepal Length", xlab="Petal Length")

dev.off()

# Using layout(matrix())
plot.matrix <- matrix(c(1,1,2,3), nrow=2, ncol=2)
plot.matrix

layout(mat=plot.matrix)

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="All Irises", ylab="Sepal Length", xlab="Petal Length")

plot(x=versicolor$Petal.Length, y=versicolor$Sepal.Length, pch=19, las=1, col=plot.colors[2], main="I. versicolor", ylab="Sepal Length", xlab="Petal Length")

plot(x=virginica$Petal.Length, y=virginica$Sepal.Length, pch=19, las=1, col=plot.colors[3], main="I. virginica", ylab="Sepal Length", xlab="Petal Length")

# Or we can adjust the widths of the plots (the heights, too). The widths= argument should be the same length as the number of columns while the heights= argument should be the same length as the number of rows. The values represent multipliers.

layout(mat=plot.matrix, widths=c(2,1))
plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="All Irises", ylab="Sepal Length", xlab="Petal Length")

plot(x=versicolor$Petal.Length, y=versicolor$Sepal.Length, pch=19, las=1, col=plot.colors[2], main="I. versicolor", ylab="Sepal Length", xlab="Petal Length")

plot(x=virginica$Petal.Length, y=virginica$Sepal.Length, pch=19, las=1, col=plot.colors[3], main="I. virginica", ylab="Sepal Length", xlab="Petal Length")

dev.off()

#----------------
# Specifying margins
#----------------
# We want to make sure not to waste white space, so we can specify the plotting margins using par(mar=c(bottom, left, top, right)). This can be important if your labels aren't fitting on the plot.
# The default is par(mar=c(5.1, 4.1, 4.1, 2.1)). 

# Let's say I want to make the text bigger because the plot is for a presentation. We can use the cex (character expander) family of arguments in the plot() call. 
#	cex= makes the plotted points larger or smaller
#	cex.main= makes the title larger
# 	cex.lab= makes the axis labels larger
#	cex.axis= makes the numbers of the axis larger

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="All Irises", ylab="Sepal Length", xlab="Petal Length")

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="All Irises", ylab="Sepal Length", xlab="Petal Length", cex=2, cex.main=4, cex.lab=3, cex.axis=1.5)

# A solution
par(mar=c(5,6,5,2))
plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="All Irises", ylab="Sepal Length", xlab="Petal Length", cex=2, cex.main=4, cex.lab=3, cex.axis=1.5)

# Note that the par(mgp=c(labels, tick marks, axis line)) method will adjust the location of those components away from the edge of the plot.
# Default is par(mgp=c(3,1,0))

par(mar=c(6,7,5,2), mgp=c(4, 1.5, 0))
plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="All Irises", ylab="Sepal Length", xlab="Petal Length", cex=2, cex.main=4, cex.lab=3, cex.axis=1.5)

#--------------
# Specifying outer margins
#--------------
# Outer margins are especially great for multi-panel plots. Let's recreate the one that we had earlier and add some outer margins with par(oma=c(bottom, left, top, right))
# The default is par(oma=c(0,0,0,0))

dev.off()
par(mfrow=c(2,2), oma=c(5,5,5,0), mar=c(2,2,4,2))

# Note here that I use the font.main= argument in each of these plot calls to make the title in normal font for the first panel and then italicized for the other panels (it is bold by default; equivalent to a font.main=2)

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="All Irises", ylab=NA, xlab=NA, cex.main=2, font.main=1)

plot(x=setosa$Petal.Length, y=setosa$Sepal.Length, pch=19, las=1, col=plot.colors[1], main="I. setosa", ylab=NA, xlab=NA, cex.main=2, font.main=3)

plot(x=versicolor$Petal.Length, y=versicolor$Sepal.Length, pch=19, las=1, col=plot.colors[2], main="I. versicolor", ylab=NA, xlab=NA, cex.main=2, font.main=3)

plot(x=virginica$Petal.Length, y=virginica$Sepal.Length, pch=19, las=1, col=plot.colors[3], main="I. virginica", ylab=NA, xlab=NA, cex.main=2, font.main=3)


#-----------------
# Base R plotting skills: Using mtext() and axis()
#-----------------
# The two add on functions that I find myself using most often are mtext(), which adds text to a margin, and axis(), which draws in a user-specified axis.

# mtext() lets me:
#	1) Make axis labels that are shared between more than one panel
#	2) rotate the y-axis label so it is easier to read for presentations

# axis() lets me:
#	1) Specify exactly how many tick marks and where they'll be
#	2) Specify exactly how those tick marks will be labelled

# We actually set up the previous 4-panel plot nicely to see the functionality of mtext() and axis(). Let's recreate that plot.

par(mfrow=c(2,2), oma=c(5,6,5,0), mar=c(2,2,3,2))

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="All Irises", ylab=NA, xlab=NA, cex.main=2, font.main=1)

plot(x=setosa$Petal.Length, y=setosa$Sepal.Length, pch=19, las=1, col=plot.colors[1], main="I. setosa", ylab=NA, xlab=NA, cex.main=2, font.main=3)

plot(x=versicolor$Petal.Length, y=versicolor$Sepal.Length, pch=19, las=1, col=plot.colors[2], main="I. versicolor", ylab=NA, xlab=NA, cex.main=2, font.main=3)

plot(x=virginica$Petal.Length, y=virginica$Sepal.Length, pch=19, las=1, col=plot.colors[3], main="I. virginica", ylab=NA, xlab=NA, cex.main=2, font.main=3)

# By default, an mtext() call will put an axis label on the specified side (1=bottom, 2=left, 3=top, 4=right) of the most recently created panel. It also puts it right at the margin. Use the line= argument to move it further from the margin however far you choose. These line units are the same units as those used by par(mar=c()).

mtext(side=1, text="Petal Length", line=3)

# We can specify outer=TRUE to put this text on the outer margins. This works only because we used par(oma=c()) to make the outer margins greater than 0 on some sides.

# Recreate the plot
par(mfrow=c(2,2), oma=c(5,6,5,0), mar=c(2,2,3,2))

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="All Irises", ylab=NA, xlab=NA, cex.main=2, font.main=1)

plot(x=setosa$Petal.Length, y=setosa$Sepal.Length, pch=19, las=1, col=plot.colors[1], main="I. setosa", ylab=NA, xlab=NA, cex.main=2, font.main=3)

plot(x=versicolor$Petal.Length, y=versicolor$Sepal.Length, pch=19, las=1, col=plot.colors[2], main="I. versicolor", ylab=NA, xlab=NA, cex.main=2, font.main=3)

plot(x=virginica$Petal.Length, y=virginica$Sepal.Length, pch=19, las=1, col=plot.colors[3], main="I. virginica", ylab=NA, xlab=NA, cex.main=2, font.main=3)

# We can just use cex= to make the text bigger here (rather than cex.lab=) because we're already in the mtext() function
mtext(side=1, text="Petal Length", outer=TRUE, line=3, cex=2)

# Do the same for the y axis
mtext(side=2, text="Sepal Length", outer=TRUE, line=3, cex=2)

# And the same for an overall title
mtext(side=3, text="Sepal length vs. petal length for Irises", outer=TRUE, line=2, cex=2.1)

# We can also use the las= argument to rotate axis text using mtext(). We can make the y-axis label in the reading direction this way.

# Recreate the plot. Add some extra outer margin space to the left.

par(mfrow=c(2,2), oma=c(5,9,5,0), mar=c(2,2,3,2))

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="All Irises", ylab=NA, xlab=NA, cex.main=2, font.main=1)

plot(x=setosa$Petal.Length, y=setosa$Sepal.Length, pch=19, las=1, col=plot.colors[1], main="I. setosa", ylab=NA, xlab=NA, cex.main=2, font.main=3)

plot(x=versicolor$Petal.Length, y=versicolor$Sepal.Length, pch=19, las=1, col=plot.colors[2], main="I. versicolor", ylab=NA, xlab=NA, cex.main=2, font.main=3)

plot(x=virginica$Petal.Length, y=virginica$Sepal.Length, pch=19, las=1, col=plot.colors[3], main="I. virginica", ylab=NA, xlab=NA, cex.main=2, font.main=3)

mtext(side=1, text="Petal Length", outer=TRUE, line=3, cex=2)

# Use las=1 here to rotate the text. Use \n to insert a carriage return.
mtext(side=2, text="Sepal\nLength", outer=TRUE, line=1.5, cex=2, las=1)

# I added an \n in the middle of the character string to insert a carriage return. I also changed the line= argument from 2 to 0.
mtext(side=3, text="Sepal length vs. petal length\nfor Irises", outer=TRUE, line=0, cex=2.1)

#-----------------
# Test yourself
#-----------------
# Add a best fit line and shaded 95% confidence interval to each panel of a 4-panel Iris plot similar to the one we just created.
 

#------------------
# Using axis()
#------------------
# Put the axis tick marks where you want them, and label them how you like

# First reset the plotting device.
dev.off()

# Here's the original plot
plot.colors <- c("violet", "purple", "blue")
color.vector <- plot.colors[iris$Species]

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="Iris sepal length vs. petal length", xlab="Petal length", ylab="Sepal length", col=color.vector)

# In the plot() call, suppress the x-axis with xaxt="n" and the y-axis with yaxt="n". Then we can draw them in how we'd like.

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="Iris sepal length vs. petal length", xlab="Petal length", ylab="Sepal length", col=color.vector, yaxt="n")

# Use axis() in a similar way to using mtext(). The side= argument specifies where you want the axis (1=bottom, 2=left, 3=top, 4=right)

axis(side=2, at=c(4.5, 6, 7.5), las=1)

# Or change how the tick marks are labeled
plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="Iris sepal length vs. petal length", xlab="Petal length", ylab="Sepal length", col=color.vector, yaxt="n")

axis(side=2, at=c(4.5, 5, 6, 7.5), labels=c("4.5", "Critical\nPoint", "6", "7.5"), las=1)

# Or meaningful points
s.tick <- mean(setosa$Sepal.Length)
ve.tick <- mean(versicolor$Sepal.Length)
vi.tick <- mean(virginica$Sepal.Length)

ticks <- c(s.tick, ve.tick, vi.tick)
ticks

plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="Iris sepal length vs. petal length", xlab="Petal length", ylab="Sepal length", col=color.vector, yaxt="n")

axis(side=2, at=ticks, labels=round(ticks,digits=2), las=1)

#-----------------
# Base R plotting skills: Using Pretty Print
#-----------------
# R can render some nice typesetting for axis labels, titles, etc.

# Here's the breakdown (I'll use the main= argument for illustration):
#	1) Want only text? Just use a text string. 
#			ex: main="text here"
#	2) Want evaluated statements and text? Use paste().
#			ex: main=paste("Mean sepal length =", mean(iris$Sepal.Length))
#	3) Want special symbols? Use expression().
#			ex: main=expression(alpha)
#	4) Want special symbols and text? Use expression(paste())
#			ex: main=expression(paste(mu, "Sepal length="))
#	5) Want text, symbols, and evaluated statements? Use bquote(). A tilda (~) goes between each component. Wrap statements you want evaluated in .()
#			ex: main=bquote(mu[sepal ~ length] ~ "=" .(mean(iris$Sepal.Length)))

# See ?plotmath for all of the fancy typesetting you can use.

# Original plot with just text as title
plot(x=iris$Petal.Length, y=iris$Sepal.Length, pch=19, las=1, main="Iris sepal length vs. petal length", xlab="Petal length", ylab="Sepal length", col=color.vector)

# Evaluated statement AND text using paste()
plot(	x=iris$Petal.Length, 
		y=iris$Sepal.Length, 
		pch=19, 
		las=1, 
		main=paste("Mean sepal length =", mean(iris$Sepal.Length)), 
		xlab="Petal length", 
		ylab="Sepal length", 
		col=color.vector)
		
# Just a special symbol? Use expression()
plot(	x=iris$Petal.Length, 
		y=iris$Sepal.Length, 
		pch=19, 
		las=1, 
		main=expression(mu), 
		xlab="Petal length", 
		ylab="Sepal length", 
		col=color.vector)
		
# Special symbol plus text? Use expression(paste())
plot(	x=iris$Petal.Length, 
		y=iris$Sepal.Length, 
		pch=19, 
		las=1, 
		main=expression(paste(mu["sepal length"])), 
		xlab="Petal length", 
		ylab="Sepal length", 
		col=color.vector)
		
# Special symbol plus text plus evaluated statements? Use bquote(). Put a ~ in between each component. Equal signs go in quotes. Wrap statements you want evaluated in .()
plot(	x=iris$Petal.Length, 
		y=iris$Sepal.Length, 
		pch=19, 
		las=1, 
		main=bquote(mu[sepal ~ length] ~ "=" ~ .(mean(iris$Sepal.Length))), 
		xlab="Petal length", 
		ylab="Sepal length", 
		col=color.vector)