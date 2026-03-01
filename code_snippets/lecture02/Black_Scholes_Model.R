
################################
# Simulated Option Prices (GBM)
################################

# Black-Scholes input 
#------------------------------
S0 <- 6500                    # Current stock price (here: S&P500 index)
T <- 1                        # Time period (e.g. one year) and maturity of the option
s <-   365                    # Time steps per year (e.g. 365 days)
dt <- rep(1/s,T*s)            # Vector of constant time steps
t <- c(0,cumsum(dt))          # Time line including start date (t=0)
paths <- 100000               # Number of simulated paths
r <- 0.01                     # Risk-free rate (r is now also the GBM drift!)
X <- 7000                     # Strike price
sigma <- 0.2                  # Return volatility 
#------------------------------

# Simulate stock prices at maturity
#-------------------------------------------
# set.seed(13)
# Simulate Wiener process levels at maturity
WT <- sqrt(T)*rnorm(paths, mean = 0, sd = 1)             
# Simulate stock prices at maturity
ST <- S0*exp((r-sigma^2/2)*T + sigma*WT)
#-------------------------------------------

# Turn simulated stock prices into option payoffs
#-------------------------------------------
CT <- rep(NA, length(ST)) # Call option 
    for(i in 1:length(ST)){CT[i] <- max(ST[i] - X, 0)}
PT <- rep(NA, length(ST)) # Put option 
    for(i in 1:length(ST)){PT[i] <- max(X - ST[i], 0)}
#-------------------------------------------

# Plot call option payoffs
#-------------------------------------------
# png(file="CallOptionPayoffs.png",width=1000, height=1000, res = 200)
hist(CT, breaks = 200, col = "darkgrey", xaxs = "i", yaxs = "i", 
     main = "Distribution of a call option payoff", freq = FALSE,
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,  
     xlim = c(0,4000), ylim = c(0,0.001),
     xlab =   expression(paste(plain(C[T]))))
# dev.off()
#-------------------------------------------

# Simulated option prices
#------------------------------
C0 <- mean(CT)*exp(-r*T)
P0 <- mean(PT)*exp(-r*T)
C0
P0
#------------------------------



###############################
# Black-Scholes Option Prices
##############################

# Black-Scholes option prices
#------------------------------
d1 <- (log(S0/X)+(r+sigma^2/2)*T) / (sigma*sqrt(T))
d2 <- d1 - sigma*sqrt(T)

C0 <- S0*pnorm(d1) - X*exp(-r*T)*pnorm(d2)
P0 <- X*exp(-r*T)*pnorm(-d2) - S0*pnorm(-d1) 

C0
P0
#------------------------------

# Put-Call Parity
#------------------------------
B0 = X*exp(-r*T)
P0 - C0 
B0 - S0
#------------------------------



###############################
# Black-Scholes Greeks
##############################

# Vary initial stock price S0
#------------------------------
S0x <- seq(0,12000,1) 

# Compute d1, d2 and option prices that correspond to the various S0
d1 <- (log(S0x/X)+(r+sigma^2/2)*T) / (sigma*sqrt(T))
d2 <- d1 - sigma*sqrt(T)
C0 <- S0x*pnorm(d1) - X*exp(-r*T)*pnorm(d2)
P0 <- X*exp(-r*T)*pnorm(-d2) - S0x*pnorm(-d1) 

# Compute the option deltas
DeltaC <- pnorm(d1)
DeltaP <- pnorm(d1)-1

# Compute the option vegas
Vega <- S0x*dnorm(d1)*sqrt(T)

# Compute the option thetas
ThetaC <- -(S0x*dnorm(d1)*sigma)/(2*sqrt(T)) - r*X*exp(-r*T)*pnorm(d2)
ThetaP <- -(S0x*dnorm(d1)*sigma)/(2*sqrt(T)) + r*X*exp(-r*T)*pnorm(-d2)

#------------------------------



###############################
# Graphical illustrations
##############################

# Plot the deltas for different S0
#--------------------------------------------------------------------
# png(file="OptionDeltas.png",width=1000, height=1000, res = 200)
plot(S0x, DeltaC, type = "l",col = "darkgreen", lwd = 6, 
     xlab = expression(paste(plain(S[0]))), ylab = "Delta",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, main = "Option Deltas",
     xaxs = "i", yaxs = "i", xlim = c(3000,11000), ylim = c(-1,1))
lines(S0x, DeltaP, lwd = 6, type = "l", col = "black")
abline(h = 0, lwd = 2, lty = "dotted")
abline(v=X, lwd = 3, lty = "solid")
abline(v=S0, lwd = 3, lty = "dashed")
legend("topleft",c("Call", "Put"), pch=c(NA,NA), 
       lty = c(1,1), lwd = 6, col = c("darkgreen", "black"),
       bty = TRUE, bg=FALSE, cex=1.25, box.col=TRUE, horiz=FALSE)
# dev.off()
#--------------------------------------------------------------------

# png(file="OptionVegas.png",width=1000, height=1000, res = 200)
#--------------------------------------------------------------------
plot(S0x, Vega, type = "l",col = "darkgreen", lwd = 6, 
     xlab = expression(paste(plain(S[0]))), ylab = "Vega",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, main = "Option Vegas",
     xaxs = "i", yaxs = "i", xlim = c(3000,11000),  ylim = c(0,3000))
abline(v=X, lwd = 3, lty = "solid")
abline(v=S0, lwd = 3, lty = "dashed")
legend("topleft",c("Call = Put"), pch=c(NA,NA), 
       lty = c(1,1), lwd = 6, col = c("darkgreen"),
       bty = TRUE, bg=FALSE, cex=1.25, box.col=TRUE, horiz=FALSE)
# dev.off()
#--------------------------------------------------------------------

# png(file="OptionThetas.png",width=1000, height=1000, res = 200)
#--------------------------------------------------------------------
plot(S0x, ThetaC, type = "l",col = "darkgreen", lwd = 6, 
     xlab = expression(paste(plain(S[0]))), ylab = "Theta",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, main = "Option Thetas",
     xaxs = "i", yaxs = "i", xlim = c(3000,11000),  ylim = c(-400,200))
abline(v=X, lwd = 3, lty = "solid")
abline(v=S0, lwd = 3, lty = "dashed")
lines(S0x, ThetaP, lwd = 6, type = "l", col = "black")
legend("bottomleft",c("Call", "Put"), pch=c(NA,NA), 
       lty = c(1,1), lwd = 6, col = c("darkgreen", "black"),
       bty = TRUE, bg=FALSE, cex=1.25, box.col=TRUE, horiz=FALSE)
# dev.off()
#--------------------------------------------------------------------

