
####################################################
# Prepare data for Geometric Brownian motion (GBM)
####################################################

# Import historical return time series 
#----------------------------------------------------------------
SPX <- read.csv2("GBMInput.csv", header = TRUE, stringsAsFactors=FALSE)
names(SPX) <- c("date","index")          # Set column names
SPX$date <- as.Date(SPX$date,"%d.%m.%Y") # Y = four-digit year, m = month, d = day
SPX$index <- as.numeric(SPX$index) 
head(SPX)
str(SPX)
#----------------------------------------------------------------

# Visually inspect index
#----------------------------------------------------------------
plot(SPX$date,SPX$index, type = "l", lwd = 3, 
     col = "darkgreen", xaxs = "i", yaxs = "i",
     cex.lab = 1.5, cex.axis = 1.5,
     ylim = c(0,7000), xlab = "Time", ylab = "S&P500 Index")
#----------------------------------------------------------------

# Compute discrete returns and log returns
#----------------------------------------------------------------
SPX$returns <- c(NA, diff(SPX$index)/SPX$index[-length(SPX$index)])
SPX$log.returns <- c(NA, log(SPX$index[-1]/SPX$index[-length(SPX$index)]))
head(SPX)
#----------------------------------------------------------------

# Use log returns for calibration
#----------------------------------------------------------------
# Convert monthly to annual figures
# Note: standard deviation converts by square-root of time
g <- mean(SPX$log.returns, na.rm = TRUE)*12
sigma <- sd(SPX$log.returns, na.rm = TRUE)*sqrt(12) 
mu <- g + sigma^2/2
#----------------------------------------------------------------



####################################
# Simulate Geometric Brownian Motion
####################################

# Input 
#---------------------------------
S0 <- SPX$index[length(SPX$index)]       # Current stock price (here: S&P500 index)
T <- 1                                   # Time period (e.g. one year)
s <- 365                                 # Time steps per year (e.g. 365 days)
dt <- rep(1/s,T*s)                       # Vector of constant time steps
t <- c(0,cumsum(dt))                     # Time line including start date (t=0)
paths <- 100000                          # Number of simulated paths
#---------------------------------

# Simulate Wiener process (increments and levels)
#-----------------------------------------------
# set.seed(13)
dW <- sqrt(dt)*rnorm(n = s*T, mean = 0, sd = 1)   # Vector of path increments
Wt <- c(0, cumsum(dW))                            # Vector of path levels (starts at W0=0)
WT <- sqrt(T)*rnorm(paths, mean = 0, sd = 1)      # Vector of simulated values at maturity
#-----------------------------------------------

# Simulate return and price paths with Euler discretization
#-----------------------------------------------
Rt <- mu*dt + sigma*dW          # Return path
St <- c(S0, S0*cumprod(1 + Rt)) # Convert into price path
# Notes: 
# 1) We use the Euler discretization as universal simulation method for didactic reasons.
# 2) Thus, we apply discrete compounding, which converges to continuous compounding for Δt -> 0
# 3) For the GBM, exact exponential transitions over any discrete time step are available.
# These exact transitions can be used to avoid time-step bias without shrinking Δt.
# We will use these exact transitions to model the money back guarantee under iv)
#-----------------------------------------------

# Plot price path
#-----------------------------------------------
# png(file="GBMpaths.png", width=2000, height=1000, res = 200)
plot(t,St, type = "l",col = "darkgreen", lwd = 3, 
     xlab = "Time t", ylab = "Stock Price",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
     xaxs = "i", yaxs = "i", ylim = c(4000,10000))
#lines(t, St2, lwd = 3, type = "l", col = "black")
#lines(t, St3, lwd = 3, type = "l", col = "darkgrey")
# dev.off()
#-----------------------------------------------

# Simulate stock price and return distribution at T
#-----------------------------------------------
ST <- S0*exp((mu-sigma^2/2)*T + sigma*WT)
RT <- (mu-sigma^2/2)*T + sigma*WT
#-----------------------------------------------

# Compute means and standard deviations 
#-----------------------------------------------
mean(ST)
sd(ST)
mean(RT)
sd(RT)
#-----------------------------------------------

# Compare to closed-form solution
#-----------------------------------------------
# Mean of the price distribution
meanST <- S0 * exp(mu*T) 
# Standard deviation of the price distribution
sdST <- sqrt(S0^2 * exp(2*mu*T) * (exp(sigma^2*T)-1)) 
# Mean of return distribution
meanRT <- (mu-sigma^2/2)*T 
# Standard deviation of the return distribution
sdRT <- sigma*sqrt(T) 
# Show results
meanST
sdST
meanRT
sdRT
#-----------------------------------------------



#############################
# Value-at-Risk
#############################

# Portfolio value and alpha level
#-----------------------------------------------
V0 <- 1000000
alpha <- 0.01
#-----------------------------------------------

# Simulated portfolio value and loss (currency)
#-----------------------------------------------
VT <- V0 * exp(RT)
LT <- V0 - VT
#-----------------------------------------------

# Analytical VaR (currency)
#-----------------------------------------------
z_alpha <- qnorm(alpha)
R_alpha_star <- (mu - sigma^2/2)*T + z_alpha * sigma*sqrt(T)
VaR_alpha <- V0 * (1 - exp(R_alpha_star))  #VaR in currency units
#-----------------------------------------------

# Convert to thousands for plotting
#-----------------------------------------------
LT_k  <- LT / 1000
VaR_k <- VaR_alpha / 1000
#-----------------------------------------------



#############################
# Graphical illustrations
#############################

# Plot price distribution
#--------------------------------------------------------------------
# png(file="GBMprices.png",width=1500, height=1000, res = 200)
hist(ST, breaks = 100, col = "darkgrey", xaxs = "i", yaxs = "i", 
     main = "Price distribution", freq = FALSE, 
     xlim = c(4000,15000), ylim = c(0,4*10^-4),
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
     xlab =   expression(paste(plain(S[T]))))
x <- seq(0,max(ST),1)
y <- dlnorm(x, meanlog = log(S0) + (mu-sigma^2/2)*T, sdlog = sigma*sqrt(T))
lines(x,y, lwd = 6, col = "darkgreen")
legend("topright",c("Simulated GBM", "Analytical GBM"), pch=c(NA,NA), 
       lty = c(1,1), lwd = 6, col = c("darkgrey", "darkgreen"),
       bty = TRUE, bg=FALSE, cex=1.25, box.col=TRUE, horiz=FALSE)
# dev.off()
#---------------------------------------------------------------------

# Plot return distribution
#--------------------------------------------------------------------
# png(file="GBMreturns.png",width=1500, height=1000, res = 200)
hist(RT, breaks = 100, col = "darkgrey", xaxs = "i", yaxs = "i", 
     main = "Return distribution", freq = FALSE, 
     xlim = c(-0.5,0.8), ylim = c(0,3),
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
     xlab =   expression(paste(plain(R["0,T"]))))
x <- seq(-1,1,0.01)
y <- dnorm(x, mean = (mu-sigma^2/2)*T, sd = sigma*sqrt(T))
lines(x,y, lwd = 6, col = "darkgreen")
legend("topright",c("Simulated GBM", "Analytical GBM"), pch=c(NA,NA), 
        lty = c(1,1), lwd = 6, col = c("darkgrey", "darkgreen"),
        bty = TRUE, bg=FALSE, cex=1.25, box.col=TRUE, horiz=FALSE)
# dev.off()
#--------------------------------------------------------------------

# Plot return distribution with VaR
#--------------------------------------------------------------------
#png(file="GBMReturnQuantile.png",width=1500, height=1000, res = 200)
hist(RT, breaks = 100, col = "darkgrey", xaxs = "i", yaxs = "i",
     main = "Return distribution with left-tail quantile",
     freq = FALSE,
     xlim = c(-0.5,0.8), ylim = c(0,3),
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,
     xlab = expression(paste(plain(R["0,T"]))))

# Analytical density
x <- seq(-1,1,0.01)
y <- dnorm(x, mean = meanRT, sd = sdRT)
lines(x,y, lwd = 6, col = "darkgreen")

# Mark quantile
abline(v = R_alpha_star, col = "red", lwd = 5)

# Optional: shade left tail
x_tail <- seq(-1, R_alpha_star, 0.01)
y_tail <- dnorm(x_tail, mean = meanRT, sd = sdRT)
polygon(c(x_tail, rev(x_tail)),
        c(y_tail, rep(0,length(y_tail))),
        col = rgb(1,0,0,0.3), border = NA)

legend("topright",
       c("Simulated GBM", "Analytical GBM", expression(R[alpha]^"*")),
       lty = c(1,1,1), lwd = c(6,6,5),
       col = c("darkgrey","darkgreen","red"),
       bty = TRUE, cex = 1.25)
#dev.off()
#--------------------------------------------------------------------

# Plot loss distribution with VaR
#--------------------------------------------------------------------
#png(file="GBMVaR.png",width=1500, height=1000, res = 200)
hist(LT_k, breaks = 100, col = "darkgrey", xaxs = "i", yaxs = "i",
     main = "Loss distribution with VaR",
     freq = FALSE,
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,
     xlab = "Loss (thousand $)")

# VaR marker
abline(v = VaR_k, col = "red", lwd = 5)

# Analytical loss density overlay:
xk <- seq(min(LT_k), min(max(LT_k), (V0 - 1) / 1000), length.out = 2000)
l  <- xk * 1000

densL <- dlnorm(V0 - l, meanlog = log(V0) + mu, sdlog = sigma)  # per-$
densL_k <- densL * 1000                                        # per-$thousand

lines(xk, densL_k, lwd = 6, col = "darkgreen")

# Shade right tail beyond VaR under analytical density
xk_tail <- xk[xk >= VaR_k]
l_tail  <- xk_tail * 1000

dens_tail <- dlnorm(V0 - l_tail, meanlog = log(V0) + mu, sdlog = sigma)  # per-$
dens_tail_k <- dens_tail * 1000                                         # per-$thousand

polygon(c(xk_tail, rev(xk_tail)),
        c(dens_tail_k, rep(0, length(dens_tail_k))),
        col = rgb(1,0,0,0.3),
        border = NA)

# Legend
legend("topright",
       c("Simulated losses", "Analytical GBM", expression(VaR[alpha])),
       lty = c(1,1,1), lwd = c(6,6,5),
       col = c("darkgrey","darkgreen","red"),
       bty = TRUE, cex = 1.25)
#dev.off()
#--------------------------------------------------------------------
