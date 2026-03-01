
############################
#Probability distributions
############################

###############################################
#Normal distribution
#############################################

#create a sequence-vector 
x <- seq(-10,10,0.01)    
#compute pdf-values corresponding to x
a <- dnorm(x, 0, 1)
b <- dnorm(x, 2, 2)
c <- dnorm(x, -1, 3)


#Visualize pdf 
#png(file="normal.png",width=1000, height=1000, res = 200)
plot(x, a, type = "l", 
     lwd = 6, xlim = c(-4,4), ylim = c(0,0.5),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     main = "Normal distributions (pdf)", 
     xlab = "x", ylab = "Density", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
lines(x, b, lwd = 6, col = "black")
lines(x, c, lwd = 6, col = "grey")
legend("topright",c("N(0,1)","N(2,2)", "N(-1,3)"), pch=c(NA,NA,NA), 
       lty = c(1,1,1), lwd = 6, col = c("darkgreen", "black", "grey"),
       bty = TRUE, bg=FALSE, cex=1.25, box.col=TRUE, horiz=FALSE)
#dev.off()



###############################################
#Standardizing variables
#############################################

#fix random number generator
set.seed(100)
#draw 10000 realizations from N(2,2)
x <- rnorm(10000, 2, 2)
mean(x)

#create variable y by demeaning x
y <- x - mean(x)
mean(y)

#create variable z by demeaning and scaling x
z <- (x-mean(x))/sd(x)
sd(z)

#Visualize x and z
#png(file="standardization.png",width=1000, height=1000, res = 200)
plot(x, pch = 19, col = "grey", cex = 0.5, xlim = c(0,length(x)), ylim = c(-4,8),     
     xaxs = "i", yaxs = "i", main = "Standardization of a random variable", 
     xlab = "Draw", ylab = "x", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
points(z, pch = 19, col = "darkgreen", cex = 0.5)
abline(h = sd(z), lwd = 3, col = "black", lty =2)
abline(h = mean(z), lwd = 3, col = "black", lty =1)
abline(h = -sd(z), lwd = 3, col = "black", lty =2)
#dev.off()



###############################################
#Lognormal distribution
#############################################

#create a sequence-vector 
x <- seq(-10,10,0.01)  
#compute pdf-values corresponding to x
b <- dlnorm(x,1.5,2)

#Visualize pdf
#png(file="lognormal.png",width=1000, height=1000, res = 200)
plot(x, b, type = "l", 
     lwd = 6, xlim = c(0,1), ylim = c(0,0.5),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     main = "Lognormal distribution (pdf)", 
     xlab = "x", ylab = "Density", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
#dev.off()



###############################################
#Poisson distribution
#############################################

#set the Poisson intensity
lambda <- 10      
#create a sequence-vector 
x <- seq(0,200,1)
#compute pmf-values corresponding to x
p <- dpois(x,lambda)

#Visualize pmf
#png(file="Poisson.png",width=1000, height=1000, res = 200)
plot(x, p, type = "h", lwd = 4,
     pch = 19, col = "darkgreen", xlim = c(0,25), ylim = c(0,0.2),
     main = "Poisson distribution (pmf)", 
     xlab = "k", ylab = "P(X = k)",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,
     xaxs = "i", yaxs = "i")
points(x, p, pch = 19, cex = 1.5, col ="darkgreen")
#dev.off()



########################################################
#Central limit theorem
#######################################################

#fix random number generator
#(create reproducible random numbers)
set.seed(100)    
#number of samples to draw  (variables in sum) 
variables <- 30   
#sample size for each variable
draws <- 10000    
#storage matrix and vector
data <- matrix(NA, draws, variables)
distr <- rep(NA, draws)
#draw 10000 realizations from an exponential distribution for each variable
for(i in 1:variables){data[,i] <- rexp(draws,0.75)}
#sum realizations across variables to generate distribution of their sum
for(i in 1:draws){distr[i] <- sum(data[i,1:variables])}

#create a histogram of the random numbers
#png(file="CLT.png",width=1000, height=1000, res = 200)
hist(distr, freq = FALSE,  breaks = 50, xlim = c(0,70), ylim = c(0,0.2),
     col = "grey", xlab = "x", ylab = "Density", 
     xaxs = "i", yaxs = "i", main = "number of variables n = 30",  
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
#overlay standard normal pdf
x <- seq(0,150,0.01)
y <- dnorm(x, mean = mean(distr), sd = sd(distr))
lines(x, y, lwd = 4, col = "darkgreen") 
legend("topright", "N(mean,sd)", lty=1, col = "darkgreen", lwd = 4, bty = TRUE)
#dev.off()



#######################################################
##Law of large numbers
#######################################################

#Illustration of the law of large numbers
#create reproducable random numbers
set.seed(10) 
#storage vector for the means
runs <- 1000
means <- seq(0,runs,1) 
#loop to create samples of growing size (1-1000)
for(i in 1:length(means)) { 
  means[i] <- mean(rnorm(i,0,1))} #store means

#Illustration of the law of large numbers
#png(file="LLN.png",width=1000, height=1000, res = 200)
plot(means, type = "l", lwd = 2, xlim = c(0,runs), 	
     ylim = c(-1,1), xaxs = "i", yaxs = "i", xlab = "No. of Trials", 
     ylab = "Sample Mean", main = 	"LLN visual example", col = "darkgreen", 
     cex.lab = 1.5, cex.axis = 1.5)
#Add thick dashed line (red) for population mean
abline(h = 0, lty = 2, lwd = 2, col = "red")
#Add legend
legend("topright", "Sampling from N(0,1)", 
       bty = "n", cex = 1.5)
#dev.off()
