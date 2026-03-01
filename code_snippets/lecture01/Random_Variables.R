
#########################
#Load R packages
#########################
library(extraDistr)
library(moments)
library(mvtnorm) 



############################
#A discrete random variable
###########################
#Note: Need to load package extraDistr first
set.seed(100) #fixes random number generator
#create a vector with 10000 random numbers X ~ discrete unif(a,b)
X <- rdunif(10000, min = 1, max = 6)
#create a vector that contains a sequence (seq) of numbers 
x <- seq(1,6,1) 
pdf <- ddunif(x, min = min(x), max = max(x)) #uniform pmf
cdf <- pdunif(x, min = min(x), max = max(x)) #uniform cdf


#Visualize probability mass function (PMF)
#png(file="pmf.png",width=1000, height=1000, res = 200)
plot(x, pdf, type = "h", 
     lwd = 6, xlim = c(0,7), ylim = c(0,0.25),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     main = "Probability mass function (pmf)", 
     xlab = "x", ylab = "Pr(X = x)", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
points(x, pdf, pch = 19, cex = 2, col ="darkgreen")
#dev.off()


#Visualize cumulative distribution function (CDF)
#png(file="cdf.png",width=1000, height=1000, res = 200)
plot(x, cdf, type = "p", 
     lwd = 6, xlim = c(0,7), ylim = c(-0.1,1.1),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     main = "Cumulative distribution function (cdf)", 
     xlab = "x", ylab = "1 - Pr(X > x)", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
segments(x0 = 0, y0 = 0, x1 = 1, y1 = 0, col = "darkgreen", lwd = 3)
segments(x0 = 1, y0 = 1/6, x1 = 2, y1 = 1/6, col = "darkgreen", lwd = 3)
segments(x0 = 2, y0 = 2/6, x1 = 3, y1 = 2/6, col = "darkgreen", lwd = 3)
segments(x0 = 3, y0 = 3/6, x1 = 4, y1 = 3/6, col = "darkgreen", lwd = 3)
segments(x0 = 4, y0 = 4/6, x1 = 5, y1 = 4/6, col = "darkgreen", lwd = 3)
segments(x0 = 5, y0 = 5/6, x1 = 6, y1 = 5/6, col = "darkgreen", lwd = 3)
segments(x0 = 6, y0 = 6/6, x1 = 7, y1 = 6/6, col = "darkgreen", lwd = 3)
#dev.off()



###############################
#A continuous random variable
##############################
set.seed(100)
#create a vector with 10000 random numbers X ~N(0,1)
X <- rnorm(10000,0,1)  
#create a sequence-vector 
x <- seq(-10,10,0.01)    
pdf <- dnorm(x)
cdf <- pnorm(x)


#Visualize probability density function (PDF)
#png(file="pdf.png",width=1000, height=1000, res = 200)
plot(x, pdf, type = "l", 
     lwd = 6, xlim = c(-4,4), ylim = c(0,0.5),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     main = "Probability density function (pdf)", 
     xlab = "x", ylab = "Density", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
xvec <-  subset(x, x >= 1 & x <= 2)
yvec <-  dnorm(xvec)
polygon(c(xvec, 2, 1), c(yvec, 0, 0), col="darkgrey")
#dev.off()


#Visualize cumulative distribution function (CDF)
#png(file="cdf.png",width=1000, height=1000, res = 200)
plot(x, cdf, type = "l", 
     lwd = 6, xlim = c(-4,4), ylim = c(0,1),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     main = "Cumulative distribution function (cdf)", 
     xlab = "x", ylab = "1 - Pr(X > x)", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
#dev.off()



############################
#Mean and median,
#Variance and standard deviation
############################

mean(X)	#compute the mean of X
median(X)	#compute the median of X
var(X)	#compute the variance of X
sd(X)	#compute the standard deviation of X


#Visualize expected value and standard deviation
#png(file="sd.png",width=1000, height=1000, res = 200)
plot(x, pdf, type = "l", 
     lwd = 6, xlim = c(-4,4), ylim = c(0,0.5),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     main = "Probability density function (pdf)", 
     xlab = "x", ylab = "Density", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)

abline(v=mean(X), lwd = 4, lty = "dashed") 

arrows(x0 = 0, y0 = dnorm(-sd(X)), 
       x1 = -sd(X), y1 = dnorm(-sd(X)), lwd = 4, lty = "solid")
arrows(x0 = 0, y0 = dnorm(+sd(X)), 
       x1 = +sd(X), y1 = dnorm(+sd(X)), lwd = 4, lty = "solid")
#dev.off()



############################
#Skewness and kurtosis
############################

#create a sequence vector 
x <- seq(-10,10,0.01)  
#compute corresponding pdf values for beta and lognormal distribution
pdf1 <- dbeta(x,5,2)
pdf2 <- dlnorm(x,-0.75,1.25)

#Graph with lognormal and beta distribution
#png(file="Skewness+Kurtosis.png",width=1000, height=1000, res = 200)
plot(x, pdf1, type = "l", 
     lwd = 6, xlim = c(0,1), ylim = c(0,3),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     main = "Skewness and kurtosis", 
     xlab = "x", ylab = "Density", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
lines(x, pdf2, lwd = 6, col ="black")
#dev.off()

#Draw realizations of beta-distributed random variable (a)
#Draw realizations of beta-distributed random variable (b)
#Draw realizations of lognormal random variable (c) 
a <- rbeta(10000,5,2)  #left-skewed beta
b <- rbeta(10000,2,5)  #right-skewed beta
c <- rlnorm(10000,0,1) #lognormal
#Then eyeball the histogram with hist() 
hist(b)

#Compute mean
mean(a)
mean(b)
mean(c)
#Compute sd
sd(a)
sd(b)
sd(c)
#Compute skewness
skewness(a)
skewness(b)
skewness(c)
#Compute kurtosis
kurtosis(a)
kurtosis(b)
kurtosis(c)



############################
#Covariance
############################

#Characterize the random variables X and Y
sd1 <- sd2 <- 1.25	#set standard deviations of X and Y
rho <- 0.6		      #set correlation between X and Y
cov <- rho*sd1*sd2	#compute covariance

#Create variance-covariance matrix and mean vector
S <- matrix(c(sd1^2,cov,cov,sd2^2),2,2)
M <- rep(0, nrow(S))

#generate data vectors
x <- seq(-4, 4, length = 50) #sequence with 50 steps between -4 and 4
y <- x #y is a clone of x 
z <- matrix(NA, nrow = length(x), ncol = length(y)) #z is matrix for the surface values
#Determine bivariate values of the bivariate normal pdf that correspond to (x,y)
for(i in 1:length(x)){
  for(j in 1:length(y)){
    z[i,j] <- dmvnorm(c(x[i],y[j]), M, S)}}

#Visualize bivariate normal pdf
#png(file="mvnormal.png",width=1000, height=1000, res = 200)
par(mar=c(3,3,3,3))
persp(x, y, z, theta = -25, phi = 25, expand = 0.75, col = "lightgrey", 
      ltheta = 120, shade = 0.25, ticktype = "simple", nticks = 1, 
      xlab = "x", ylab = "y", zlab = "z", box = TRUE, 
      main = "Joint distribution of two random variables", xaxs = "i", yaxs = "i", 
      cex.axis = 1, cex.lab = 1, cex.main = 1.25)
#dev.off()



############################
#Correlations
############################

#Characterize the random variables X and Y
sd1 <- sd2 <- 1.25	#set standard deviations of X and Y
rho <- 0.8		      #set correlation between X and Y
cov <- rho*sd1*sd2	#compute covariance

#Create variance-covariance matrix and mean vector
S <- matrix(c(sd1^2,cov,cov,sd2^2),2,2)
M <- rep(0, nrow(S))

#create a matrix with 10000 random numbers for X and Y ~N2(0,1)
N <- rmvnorm(n = 10000, mean = M, sigma = S, method=c("chol"))
X <- N[,1]
Y <- N[,2]

#create scatterplot of X and Y to illustrate the correlation
#png(file="correlation",width=1000, height=1000, res = 200)
#par(mar=c(5,5,3,5)+0.1) #sets margins: bottom, left, top, and right 
plot(X,Y, pch = 1, cex = 0.5, xlim = c(-4,4), ylim = c(-4,4),
     xaxs = "i", yaxs = "i", xlab = "Random variable X", 
     ylab = "Random variable Y", main = "Scatterplot (positive correlation)",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5) 
#dev.off()


