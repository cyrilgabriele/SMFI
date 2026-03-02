
###############################
# Implied Volatilities
##############################

# Black-Scholes parameters
#------------------------------
S0 <- 21
X <- 20
T <- 0.25/
r <- 0.01
sigma <- 0.235
#------------------------------

# Create functions for Black-Scholes option prices
#--------------------------------------------------
C0 <- function(sigma, S0, X, T, r){                         
  d1 <- (log(S0/X)+(r+sigma^2/2)*T) / (sigma*sqrt(T))
  d2 <- d1 - sigma*sqrt(T)
  S0*pnorm(d1) - X*exp(-r*T)*pnorm(d2)}  

C0(sigma, S0, X, T, r)

P0 <- function(sigma, S0, X, T, r){
  d1 <- (log(S0/X)+(r+sigma^2/2)*T) / (sigma*sqrt(T))
  d2 <- d1 - sigma*sqrt(T)
  X*exp(-r*T)*pnorm(-d2) - S0*pnorm(-d1)} 

P0(sigma, S0, X, T, r)
#--------------------------------------------------

# Treat B-S prices as observed market prices
#--------------------------------------------------
C0m <- C0(sigma, S0, X, T, r)
P0m <- P0(sigma, S0, X, T, r)
#--------------------------------------------------

# Create function for root search
#--------------------------------------------------
fc <- function(sigma){abs(C0(sigma, S0, X, T, r) - C0m)}
fp <- function(sigma){abs(P0(sigma, S0, X, T, r) - P0m)}
#--------------------------------------------------

# Visualize difference function 
#--------------------------------------------------
# png(file="ImpliedVol.png", width=1000, height=1000, res = 200)
plot(fc, type = "l",col = "darkgreen", lwd = 6, 
     xlab = "Volatility Parameter", 
     ylab = "|BS Price - Market Price|",
     main = "Extracting implied volatility",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
     xaxs = "i", yaxs = "i")
# dev.off()
#--------------------------------------------------

# Iterative search procedure
#--------------------------------------------------
# Note: optimize searches the interval for a 
# minimum or maximum of the function f 
# with respect to its first argument.
optimize(f = fc, interval = c(0,1), tol = 0.0001)
optimize(f = fp, interval = c(0,1), tol = 0.0001)
#--------------------------------------------------



#################################
# Implied volatility smile/skew
###############################

# DAX Option data (March 3, 2022)
#--------------------------------------------------
S0 <- 13980
T <- 9/12
r <- 0.001
X <- c(12500, 13000, 13500, 14000, 14500, 15000, 
       15500, 16000, 16500, 17000, 17500, 18000) 
C0m <- c(22.14, 18.53, 15.13, 12.05, 9.31, 6.92, 
         4.93, 3.38, 2.23, 1.44, 0.93, 0.62)*100
#--------------------------------------------------

# Extract implied volatilities
#-----------------------------------------------------------------
IV <- rep(NA, length(X))   # Storage vector
for(i in 1:length(X)){     # Back out values across spectrum of strikes
  fc <- function(sigma){abs(C0(sigma, S0, X[i], T, r) - C0m[i])}        
  IV[i] <- optimize(f = fc, interval = c(0,1), tol = 0.0001)$minimum}
#-----------------------------------------------------------------

# Visualize skew
#-----------------------------------------------------------------
# png(file="VolatilitySkew.png", width=1000, height=1000, res = 200)
plot(X/S0, IV, type = "l",col = "darkgreen", lwd = 6, 
     xlab = expression(paste(plain(X/S[0]))), ylab = "Implied Vol.",
     main = "Volatility Skew (DAX Options)",
     xlim = c(0.85,1.3), ylim = c(0.1,0.35),
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
     xaxs = "i", yaxs = "i")
abline(v = 1, lwd = 4, lty = "dashed") 
# dev.off()
#-----------------------------------------------------------------

