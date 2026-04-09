
#########################################
# Appendix A. Markowitz with Investment Limits
#########################################

# Packages and working directory
#-----------------------------------------------------------
library(quadprog)    # For short-sale constrained (SSC) frontier
library(LowRankQP)  # For frontier with investment limits (IL)
#getwd()             # Get working directory
#setwd("")           # Set working directory
#-----------------------------------------------------------


# Import data
#-----------------------------------------------------------
# Save time series data from .csv file into a data frame
data <- read.csv2("Asset_Classes.csv", header = TRUE, stringsAsFactors = FALSE)

# Convert the first column of the data frame to date format
data$Date <- as.Date(data$Date, "%Y-%m-%d")

# Convert remaining columns to numerical format
for(i in 2:length(data[1,])){data[,i] <- as.numeric(data[,i])}

# Display head and check structure (str) of the data frame
head(data)
str(data)
#-----------------------------------------------------------


# Estimate vector of mean returns
#-----------------------------------------------------------
# Number of assets (subtract date column)
N <- length(data[1,]) - 1

# Create vector to store annualized mean returns
MU <- rep(NA, N)

# Fill vector with annualized mean returns
for(i in 2:length(data[1,])){
  MU[i-1] <- round(mean(data[,i]) * 12, 4)}
#-----------------------------------------------------------


# Estimate variance-covariance matrix
#-----------------------------------------------------------
# Annualized variance-covariance matrix
SIGMA <- cov(data[,2:length(data[1,])]) * 12

# Variances found on main diagonal
VAR <- diag(SIGMA)

# Compute return standard deviations (volatilities)
SD <- sqrt(VAR)
#-----------------------------------------------------------


# Compute mean and standard deviation for an equal-weight portfolio
#-----------------------------------------------------------
w.eq <- rep(1/N, N)                                        # Equal-weight vector
muP.eq <- as.numeric(t(w.eq) %*% MU)                      # Portfolio mean return
sigmaP.eq <- as.numeric(sqrt(t(w.eq) %*% SIGMA %*% w.eq)) # Portfolio standard deviation
muP.eq
sigmaP.eq
#-----------------------------------------------------------



########################################
# Short-sale constrained (SSC) frontier
########################################

#### Prepare objects for solve.QP ####
#-----------------------------------------------------------
# Define number and type of constraints
n <- 2 + N    # Total number of constraints (2 equality + N short-sale)
meq <- 2      # First meq constraints treated as equalities

# dvec set to zero (we specify the target function in terms of the variance only)
dvec <- rep(0, N)

# Amat includes left-hand-side values of the constraints
Amat <- matrix(NA, N, n)
Amat[,1] <- MU              # (i)   w1*mu1 + ... + wN*muN = muP  (target return)
Amat[,2] <- rep(1, N)       # (ii)  SUM(w_i) = 1                 (budget constraint)
for(i in 1:N){              # (iii) w_i >= 0                      (short-sale constraints)
  constraint <- rep(0, N)
  constraint[i] <- 1
  Amat[,i+2] <- constraint}
#-----------------------------------------------------------

# Set range of portfolio means along which to minimize variance
#-----------------------------------------------------------
stepsize <- 0.00011
mu_bar <- seq(min(MU), max(MU), stepsize)
#-----------------------------------------------------------

# Create objects to store SSC efficient frontier
#-----------------------------------------------------------
w.SSC <- matrix(NA, N, length(mu_bar))       # Matrix for weights
muP.SSC <- rep(NA, length(mu_bar))           # Vector for means
sigmaP.SSC <- rep(NA, length(mu_bar))        # Vector for standard deviations
#-----------------------------------------------------------

#### Optimization for SSC Efficient Frontier ####
#-----------------------------------------------------------
for(i in 1:length(mu_bar)){
  w.SSC[,i] <- solve.QP(Dmat = SIGMA,
                         dvec = dvec,
                         Amat = Amat,
                         bvec = c(mu_bar[i], 1, rep(0, N)),
                         meq = meq)$solution
  muP.SSC[i] <- t(w.SSC[,i]) %*% MU
  sigmaP.SSC[i] <- sqrt(t(w.SSC[,i]) %*% SIGMA %*% w.SSC[,i])}
#-----------------------------------------------------------



##############################################
# Frontier with investment limits (IL)
##############################################

#### Prepare objects for LowRankQP ####
#-----------------------------------------------------------
# The LowRankQP function minimizes (?LowRankQP for details):
# min  1/2 * t(w) %*% Vmat %*% w + t(dvec) %*% w
# s.t. Amat %*% w = bvec  and  0 <= w <= uvec
#
# Note: short-sale constraints are included automatically (0 <= w).
# Amat only includes the equality constraints (transposed format).

Amat.IL <- matrix(NA, N, 2)
Amat.IL[,1] <- MU              # (i)  w1*mu1 + ... + wN*muN = muP  (target return)
Amat.IL[,2] <- rep(1, N)       # (ii) SUM(w_i) = 1                 (budget constraint)
Amat.IL <- t(Amat.IL)          # Transpose (required by LowRankQP)

# uvec contains upper bounds (investment limits) for each asset
uvec <- c(0.20, 1.00, 0.10, 0.25, 0.05, 1.00)
#-----------------------------------------------------------

# Determine achievable portfolio return range under investment limits
#-----------------------------------------------------------
# Maximum return: invest max allowed amounts in highest-yielding asset classes
Stocks       <- uvec[1]                            # 20% max in stocks
CorpBonds    <- uvec[3]                            # 10% max in corporate bonds
Alternatives <- uvec[5]                            # 5% max in alternatives
RealEstate   <- 0                                  # Zero in real estate
MoneyMarket  <- 0                                  # Zero in money market
Treasuries   <- 1 - (uvec[1] + uvec[3] + uvec[5]) # Residual in treasuries

w.max.r <- c(Stocks, Treasuries, CorpBonds, RealEstate, Alternatives, MoneyMarket)
max.r <- as.numeric(w.max.r %*% MU)

# Minimum return: invest everything in the lowest-yielding asset class
w.min.r <- c(0, 0, 0, 0, 0, uvec[6])
min.r <- as.numeric(w.min.r %*% MU)

max.r
min.r
#-----------------------------------------------------------

# Set range of portfolio means along which to minimize variance
#-----------------------------------------------------------
mu_bar.IL <- seq(min.r, max.r, stepsize)
#-----------------------------------------------------------

# Create objects to store IL efficient frontier
#-----------------------------------------------------------
w.IL <- matrix(NA, N, length(mu_bar.IL))       # Matrix for weights
muP.IL <- rep(NA, length(mu_bar.IL))            # Vector for means
sigmaP.IL <- rep(NA, length(mu_bar.IL))         # Vector for standard deviations
#-----------------------------------------------------------

#### Optimization for IL Efficient Frontier ####
#-----------------------------------------------------------
for(i in 1:length(mu_bar.IL)){
  w.IL[,i] <- LowRankQP(Vmat = SIGMA,
                         dvec = rep(0, N),
                         Amat = Amat.IL,
                         bvec = c(mu_bar.IL[i], 1),
                         uvec = uvec,
                         method = "CHOL")$alpha
  muP.IL[i] <- t(w.IL[,i]) %*% MU
  sigmaP.IL[i] <- sqrt(t(w.IL[,i]) %*% SIGMA %*% w.IL[,i])}
#-----------------------------------------------------------



#########################################
# Plot: SSC vs. IL efficient frontiers
#########################################

#-----------------------------------------------------------
# png(file="InvestmentLimits.png", width=1000, height=1000, res = 200)
plot(sigmaP.SSC, muP.SSC, type = "l",
     lwd = 6, xlim = c(0, 0.12), ylim = c(0, 0.12),
     xaxs = "i", yaxs = "i", col = "darkgreen",
     main = "Investment Limits",
     xlab = "Return Standard Deviation", ylab = "Expected Return",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
lines(sigmaP.IL, muP.IL, lwd = 6, col = "black")
points(SD, MU, pch = 21, cex = 2, lwd = 2, col = "black", bg = "darkgrey")
legend("bottomright", c("Short-Sale Constraints", "Investment Limits"),
       lty = 1, lwd = 6, col = c("darkgreen", "black"), cex = 1.2)
# dev.off()
#-----------------------------------------------------------



#########################################
# Appendix B. Two-Fund Theorem (Tobin)
#########################################

# Packages and working directory
#-----------------------------------------------------------
library(quadprog)   # Load quadprog package
#getwd()            # Get working directory
#setwd("")          # Set working directory
#-----------------------------------------------------------


# Import data
#-----------------------------------------------------------
# Save time series data from .csv file into a data frame
data <- read.csv2("SMI_Stocks.csv", header = TRUE, sep = ";",
                  stringsAsFactors = FALSE, strip.white = TRUE)

# Convert the first column of the data frame to date format
data$Date <- as.Date(data$Date, "%d.%m.%Y")

# Convert remaining columns to numerical format
for(i in 2:length(data[1,])){data[,i] <- as.numeric(data[,i])}

# Display head and check structure (str) of the data frame
head(data)
str(data)
#-----------------------------------------------------------


# Compute returns from prices
#-----------------------------------------------------------
# Clone data frame to store returns (instead of prices)
returns <- data

# Overwrite prices in first line with NAs
returns[1,2:length(returns[1,])] <- NA

# Compute and store returns
for(j in 2:length(data[1,])){                   # Outer loop: columns
  for(i in 2:length(data[,1])){                 # Inner loop: rows
    returns[i,j] <- data[i,j]/data[i-1,j]-1  }} # Save returns to dataframe

# Check return data frame
head(returns)
#-----------------------------------------------------------


# Estimate vector of mean returns
#-----------------------------------------------------------
# Number of assets (subtract date column)
N <- length(data[1,]) - 1

# Create vector to store mean returns
MU <- rep(NA, N)

# Fill vector with annualized mean returns
for(i in 2:length(returns[1,])){
  MU[i-1] <- round(mean(returns[,i], na.rm = TRUE) * 12, 4)}
#-----------------------------------------------------------


# Estimate variance-covariance matrix
#-----------------------------------------------------------
# Na.omit for case-wise exclusion
SIGMA <- cov(na.omit(returns[,-1])) * 12

# Variances found on main diagonal
VAR <- diag(SIGMA)

# Compute return standard deviations (volatilities)
SD <- sqrt(VAR)
#-----------------------------------------------------------



########################################
# SSC efficient frontier
########################################

#### Prepare objects for solve.QP ####
#-----------------------------------------------------------
n <- 2 + N     # Total number of constraints
meq <- 2       # First meq constraints treated as equalities
dvec <- rep(0, N)

Amat <- matrix(NA, N, n)
Amat[,1] <- MU              # (i)   w1*mu1 + ... + wN*muN = muP  (target return)
Amat[,2] <- rep(1, N)       # (ii)  SUM(w_i) = 1                 (budget constraint)
for(i in 1:N){              # (iii) w_i >= 0                      (short-sale constraints)
  constraint <- rep(0, N)
  constraint[i] <- 1
  Amat[,i+2] <- constraint}
#-----------------------------------------------------------

# Set range of portfolio means along which to minimize variance
#-----------------------------------------------------------
stepsize <- 0.00011
mu_bar <- seq(min(MU), max(MU), stepsize)
#-----------------------------------------------------------

# Create objects to store efficient frontier
#-----------------------------------------------------------
w <- matrix(NA, N, length(mu_bar))        # Matrix for weights
mu <- rep(NA, length(mu_bar))             # Vector for means
sigma <- rep(NA, length(mu_bar))          # Vector for standard deviations
#-----------------------------------------------------------

#### Optimization for Efficient Frontier ####
#-----------------------------------------------------------
for(i in 1:length(mu_bar)){
  w[,i] <- solve.QP(Dmat = SIGMA,
                     dvec = dvec,
                     Amat = Amat,
                     bvec = c(mu_bar[i], 1, rep(0, N)),
                     meq = meq)$solution
  mu[i] <- t(w[,i]) %*% MU
  sigma[i] <- sqrt(t(w[,i]) %*% SIGMA %*% w[,i])}
#-----------------------------------------------------------

# Find minimum variance portfolio
#-----------------------------------------------------------
min_var <- which.min(sigma)
#-----------------------------------------------------------



########################################
# Tobin's two-fund theorem
########################################

# Add a risk-free rate
#-----------------------------------------------------------
rf <- 0.02
#-----------------------------------------------------------

# Numerically identify the tangency portfolio
#-----------------------------------------------------------
# Step 1: Compute midpoints between consecutive frontier points
y <- (mu[1:(length(mu)-1)] + mu[2:length(mu)]) / 2
x <- (sigma[1:(length(sigma)-1)] + sigma[2:length(sigma)]) / 2

# Step 2: Numerical differentiation for the frontier slope at each midpoint
slopesCurve <- diff(mu) / diff(sigma)

# Step 3: Compute straight-line slopes from rf to each frontier midpoint
slopesLine <- (y - rf) / x

# Step 4: Restrict to upper branch of the hyperbola (positive slopes only)
negative_indices <- which(slopesCurve <= 0)
slopesCurve[negative_indices] <- NA
slopesLine[negative_indices] <- NA

# Step 5: Find the point where curve slope equals line slope (tangency condition)
differences <- (slopesLine - slopesCurve)^2
tangency <- which.min(differences)
#-----------------------------------------------------------



#########################################
# Plot: Two-Fund Theorem
#########################################

#-----------------------------------------------------------
# png(file="Tobin.png", width=1000, height=1000, res = 200)
plot(sigma, mu, type = "l",
     lwd = 6, xlim = c(0, 0.3), ylim = c(0, 0.3),
     xaxs = "i", yaxs = "i", col = "darkgreen",
     main = "Two-Fund Theorem",
     xlab = "Return Standard Deviation", ylab = "Expected Return",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
points(SD, MU, pch = 16, cex = 1, lwd = 2, col = "black")
abline(a = rf, b = slopesLine[tangency], lwd = 2)
points(sigma[tangency], mu[tangency], pch = 19, cex = 1.5, lwd = 2, col = "black", bg = "darkgrey")
points(sigma[min_var], mu[min_var], pch = 19, cex = 1.5, lwd = 2, col = "black", bg = "darkgrey")
# dev.off()
#-----------------------------------------------------------
