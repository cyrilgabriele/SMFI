
##################################
# Money-back guarantee product
# (with regular premium payments)
##################################

# Product parameters 
#----------------
premium <- 100                           # Annual premium payments
T <- 10                                  # Maturity of the contract (in years)
F0 <- 0                                  # Product begins with first premium
#----------------

# GBM input
#----------------
s <- 1                                   # Time steps per year (e.g. 365 days)
dt <- rep(1/s,T*s)                       # Vector of constant time steps
t <- c(0,cumsum(dt))                     # Time line including start date (t=0)
r <-  0.01                               # Risk-free rate for drift
sigma <- 0.15                            # Return standard deviation
paths <- 1000000                          # Number of simulated paths
#----------------

# Storage matrices
#------------------------------
Rt <- matrix(NA, paths, T*s)    # GBM returns
Ft <- Rt                        # GBM fund levels
#------------------------------



##############################################
# Run simulations with closed-form transitions
##############################################

# Wiener process and continuous GBM returns under Q
#-----------------------------------------
set.seed(13)                              # Use common random numbers
for(i in 1:paths){
  dW <- sqrt(dt)*rnorm(T*s, 0, 1)         # Wiener process increments
  Rt[i,] <- (r-sigma^2/2)*dt + sigma*dW } # GBM returns (continuous)
head(Rt)
#-----------------------------------------

# Turn GBM return paths into paths of the fund account
#----------------------------------------------
Ft[,1] <- (rep(F0,paths) + premium)*exp(Rt[,1])
for(i in 1:paths){
for(t in 2:T){
  Ft[i,t] <- (Ft[i,t-1] + premium)*exp(Rt[i,t]) }} 
head(Ft)
#----------------------------------------------

# Payoffs at maturity
#-----------------------------------------------------
FT <- Ft[,T] # Extract fund distribution at maturity T
GT <- FT     # Clone vector FT
LT <- FT     # Clone vector FT
for(i in 1:paths){ # Run loop to compute payoffs GT and LT
  GT[i] <- max(T*premium - FT[i], 0) # Fill vector GT with values
  LT[i] <- max(T*premium, FT[i])}    # Fill vector LT with values
#-----------------------------------------------------

# Price of the guarantee
#--------------------------------------
G01 <- exp(-r*T)*mean(GT)
G01
# Interpret guarantee price (in currency units) 
# relative to annual premium (85%) or to 
# total premiums paid (8.5%)
#--------------------------------------



##########################################
# Run simulations with discrete returns
##########################################

# Discretized GBM return paths under Q
# (Run in loop to save path after path)
#--------------------------------------
set.seed(13)                       # Use common random numbers
for(i in 1:paths) {
  dW <- sqrt(dt)*rnorm(s*T, 0, 1)  # Wiener process increments
  Rt[i,] <- r*dt + sigma*dW }      # GBM returns (discrete)
head(Rt)
#--------------------------------------

# Turn GBM return paths into paths of the fund account
#-----------------------------------------------------
Ft <- Rt # Clone Rt
Ft[,1] <- (rep(F0,paths) + premium)*(1 + Rt[,1])
for(i in 1:paths){
for(t in 2:T){
  Ft[i,t] <- (Ft[i,t-1] + premium)*(1 + Rt[i,t]) }} 
#-----------------------------------------------------

# Payoffs at maturity
#-----------------------------------------------------
FT <- Ft[,T] # Extract fund distribution at maturity T
GT <- FT     # Clone of FT
LT <- FT     # Clone of FT
for(i in 1:paths){ # Run loop to compute payoffs (GT and LT)
  GT[i] <- max(T*premium - FT[i], 0) # Fill vector GT with values
  LT[i] <- max(T*premium, FT[i])}    # Fill vector LT with values
#-----------------------------------------------------

# Price of the guarantee
#--------------------------------------
G02 <- exp(-r*T)*mean(GT)
G02
#--------------------------------------

# Compare exact transitions to Euler
#----------------------------------------------------------------------
D <- G02 - G01 # compute difference in guarantee price
D              # return approximation error (price difference)
D*1000000      # the difference matters in a portfolio of 1 mn contracts
               # Euler leads to overreserving of the guarantee 
#----------------------------------------------------------------------

# Critical reflection:
#-----------------------------------------------------------
# Increasing the no. of paths will not narrow the difference
# It will only make the simulation results more stable (lower standard error)
# The difference will only shrink for higher time grid density
# I.e., if we introduce intra-year time-steps (less time-step bias)
# However, that will complicate the simulation
#-----------------------------------------------------------




##############################################
# Graphically illustrate payoff distributions
##############################################

# Approximate FT (sum of lognormals) by matching the first two moments
# (Fenton-Wilkinson approximation, see Fenton, 1960)
#--------------------------------------
ml <- log(mean(FT)^2/sqrt(var(FT)+mean(FT)^2))
sl <- sqrt(log(var(FT)/mean(FT)^2 +1))
#--------------------------------------

# Plot fund value distribution at time T
#-----------------------------------------------------
# png(file="FundPayoff.png",width=1000, height=1000, res = 200)
hist(FT, breaks = 200, col = "darkgrey", xaxs = "i", yaxs = "i", 
     main = "Fund account distribution", freq = FALSE, 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
     xlim = c(0, 3000), ylim = c(0, 0.002),
     xlab =   expression(paste(plain(F[T]))))
abline(v = T*premium, lwd = 3, lty = 2)
x <- seq(0,4000,0.01)
# y <- dnorm(x, mean(FT), sd(FT))
y <- dlnorm(x, ml, sl)
lines(x,y, lwd = 4, col = "darkgreen")
# dev.off()    
#-----------------------------------------------------

# Plot guarantee payoff distribution at time T
#-----------------------------------------------------
# png(file="GuaranteePayoff.png",width=1000, height=1000, res = 200)
hist(GT, breaks = 30, col = "darkgrey", xaxs = "i", yaxs = "i", 
     main = "Guarantee payoff distribution", freq = FALSE, 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
     xlim = c(0, 3000), ylim = c(0, 0.002),
     xlab =   expression(paste(plain(G[T]))))
abline(v = T*premium, lwd = 3, lty = 2)
x <- seq(0,4000,0.01)
y <- dnorm(x, mean(FT), sd(FT))
y <- dlnorm(x, ml, sl)
lines(x,y, lwd = 4, col = "darkgreen")
# dev.off()    
#-----------------------------------------------------

# Plot combined product payoff distribution at time T
#-----------------------------------------------------
# png(file="ProductPayoff.png",width=1000, height=1000, res = 200)
hist(LT, breaks = 200, col = "darkgrey", xaxs = "i", yaxs = "i", 
     main = "Product payoff distribution", freq = FALSE, 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
     xlim = c(0, 3000), ylim = c(0, 0.002),
     xlab =   expression(paste(plain(L[T]))))
abline(v = T*premium, lwd = 3, lty = 2)
x <- seq(0,4000,0.01)
# y <- dnorm(x, mean(FT), sd(FT))
y <- dlnorm(x, ml, sl)
lines(x,y, lwd = 4, col = "darkgreen")
# dev.off()    
#-----------------------------------------------------

