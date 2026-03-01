
#######################
##Stochastic processes 
#######################

#Input for Wiener process
T <- 1                          # Time period in years 
N <- 365                        # Time steps per year (e.g. 365 days)
dt <- rep(1/N,T*N)              # Vector of constant time steps
t <- c(0,cumsum(dt))            # Time line including start date (t=0)


#fix random number generator
set.seed(100)                                    
#Simulate increments (dW is the vector of increments)
dW <- rnorm(length(dt), mean = 0, sd = 1)*sqrt(dt)


#Visualize distribution of increments
hist(dW, freq = FALSE,  breaks = 30, 
     col = "grey", xlab = "x", ylab = "Density", 
     xaxs = "i", yaxs = "i", main = "Histogram of dW",  
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
#overlay standard normal pdf
x <- seq(-1,1,0.001)
y <- dnorm(x, mean = 0, sd = sqrt(dt[1]))
lines(x, y, lwd = 4, col = "darkgreen") 
legend("topright", "N(0,sqrt(dt))", lty=1, col = "darkgreen", lwd = 4, bty = TRUE)


#Turn increments into a path
W0 <- 0                       #starting value is W0
Wt <- c(W0, W0 + cumsum(dW))  #Wt contains the path

#Add a drift and variance rate
#Convert the path to a generalized Wiener process
a <- 0.75                                        #constant drift
b <- 1                                           #constant variance
X0 <- 0                                          #starting value is X0
Xt <- c(X0, X0 + cumsum(a*dt) + b*cumsum(dW))    #Xt contains the path


#Plot sample paths of Wt and Xt
#png(file="Wiener.png",width=2000, height=1000, res = 200)
plot(t, Wt, type = "l", 
     lwd = 3, xlim = c(0,T), ylim = c(-0.4,1),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     #main = "Wiener processes", 
     xlab = "Time t", ylab = "X", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
lines(t, Xt, lwd = 3, col = "black")
legend("topleft",c("Generalized Wiener process", "Simple Wiener process"), pch=c(NA,NA), 
       lty = c(1,1), lwd = 6, col = c("black", "darkgreen"),
       bty = TRUE, bg=FALSE, cex=1.25, box.col=TRUE, horiz=FALSE)
abline(h=0, lty=2, lwd = 2)
lines(t,a*t, lty=2, lwd = 2)
#dev.off()




##########################
##Monte Carlo Simulations
##########################

#Input for simulation of Wiener process
paths <- 1000                   # number of paths to be simulated
T <- 1                          # Time period in years
N <- 365                        # Time steps per year (e.g. 365 days)
dt <- rep(1/N,T*N)              # Vector of constant time steps
t <- c(0,cumsum(dt))            # Time line including start date (t=0)
W0 <- 0                         # Starting value of the paths

#fix random number generator
#set.seed(100)                                    
#create storage matrix for the paths
Wt <- matrix(NA, paths, length(t))
#paths will be stored in rows
#columns represent the time steps


#Simulate increments for the chosen number of paths
for(i in 1:paths){
        #creates a vector of increments (fresh in each iteration of the loop)
        dW <- sqrt(dt)*rnorm(length(dt), mean = 0, sd = 1)  
        #turns the increments into a path and stores the path in row i of Wt
        Wt[i,] <- c(W0, W0 + cumsum(dW))}                

#Visualize all paths
#png(file="Monte Carlo.png",width=1000, height=1000, res = 200)
plot(t, Wt[1,], type = "l", 
     lwd = 3, xlim = c(0,T), ylim = c(-4,4),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     main = "Monte Carlo Simulation", 
     xlab = "Time t", ylab = "X", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
for(i in 2:length(Wt[,1])){lines(t, Wt[i,], lwd = 1, col = "darkgrey")}
lines(t,Wt[1,], lwd = 3, col = "darkgreen")
abline(h=0, lty=2, lwd = 2)
#dev.off()


