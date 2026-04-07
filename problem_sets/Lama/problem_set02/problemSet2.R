############################################################
# Problem Set 2: Adriano Lama
############################################################

#----------------------------------------------------------
#Exercise 1
#----------------------------------------------------------

# Save time series data from .csv file into a data frame
data <- read.csv2("Asset_Classes.csv", header = TRUE, stringsAsFactors = FALSE)

# Check structure
str(data)

# Convert the first column of the data frame to date format
data$Date <- as.Date(data$Date, "%d.%m.%Y")

# Convert the second column of the data frame to numerical format
for(i in 2:length(data[1,])){data[,i] <- as.numeric(data[,i])}

# Display head and check structure (str) of the data frame
head(data)  
str(data)


#----------------------------------------------------------
# Exercise 2
#----------------------------------------------------------


# Load package
#-----------------------------------------------------------
library(quadprog)


# Prepare asset return data
#-----------------------------------------------------------
# Exclude Date column and risk-free rate column
returns <- data[, 2:(ncol(data)-1)]

# Number of assets
N <- length(returns[1,])


# Estimate vector of mean returns
#-----------------------------------------------------------
# Create vector to store annualized mean returns
MU <- rep(NA, N)

# Fill vector with annualized mean returns
for(i in 1:N){
  MU[i] <- round(mean(returns[,i], na.rm = TRUE) * 12, 4)
}
MU


# Estimate variance-covariance matrix
#-----------------------------------------------------------
# Annualized variance-covariance matrix
SIGMA <- cov(returns, use = "complete.obs") * 12

# Variances found on main diagonal
VAR <- diag(SIGMA)

# Compute return standard deviations (volatilities)
SD <- sqrt(VAR)

# Display results
SIGMA
SD


##################
# Markowitz model
##################

#### Prepare objects for optimization function ####
#-----------------------------------------------------------
# Define number and type of constraints
n <- 2 + N     # Total number of constraints
meq <- 2       # First meq constraints treated as equalities

# dvec set to zero (we specify the target function in terms of the variance only)
dvec <- rep(0, N)

# Amat includes left-hand-side values of the constraints
Amat <- matrix(NA, N, n)
Amat[,1] <- MU              # (i)   w1*mu1 + ... + wN*muN = muP  (target return)
Amat[,2] <- rep(1, N)       # (ii)  SUM(w_i) = 1                 (budget constraint)
for(i in 1:N){              # (iii) w_i >= 0                     (short-sale constraints)
  constraint <- rep(0, N)
  constraint[i] <- 1
  Amat[,i+2] <- constraint
}


# Set range of portfolio means along which to minimize variance
#-----------------------------------------------------------
# Use slightly larger step size so the grid does not hit max(MU) exactly
# and solve.QP avoids the degenerate corner solution.
stepsize <- 0.0001
mu_bar <- seq(min(MU), max(MU), by = stepsize)
mu_bar <- mu_bar[mu_bar < max(MU)]


# Create objects to store efficient frontier
#-----------------------------------------------------------
w <- matrix(NA, N, length(mu_bar))        # Matrix for weights
mu <- rep(NA, length(mu_bar))             # Vector for means
sigma <- rep(NA, length(mu_bar))          # Vector for standard deviations


#### Optimization for Efficient Frontier ####
#-----------------------------------------------------------
for(i in 1:length(mu_bar)){
  w[,i] <- solve.QP(Dmat = SIGMA,
                    dvec = dvec,
                    Amat = Amat,
                    bvec = c(mu_bar[i], 1, rep(0, N)),
                    meq = meq)$solution
  mu[i] <- t(w[,i]) %*% MU
  sigma[i] <- sqrt(t(w[,i]) %*% SIGMA %*% w[,i])
}


# Find minimum variance portfolio
#-----------------------------------------------------------
min_var <- which.min(sigma)

# Weights of the minimum variance portfolio
w[, min_var]


# Show weights of all efficient portfolios
#-----------------------------------------------------------
colnames(w) <- paste0("Portfolio_", 1:ncol(w))
rownames(w) <- colnames(returns)


#########################################
# Plot the efficient frontier
#########################################

#-----------------------------------------------------------
plot(sigma, mu, type = "l",
     lwd = 6, xlim = c(0, max(sigma)*1.1), ylim = c(0, max(mu)*1.1),
     xaxs = "i", yaxs = "i", col = "darkgreen",
     main = "Efficient Frontier",
     xlab = "Return Standard Deviation", ylab = "Expected Return",
     cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.3)

points(SD, MU, pch = 16, cex = 1.2, lwd = 2, col = "black")
points(sigma[min_var], mu[min_var], pch = 19, cex = 1.5, lwd = 2,
       col = "black", bg = "darkgrey")

text(SD, MU, labels = colnames(returns), pos = 4, cex = 0.8)


#----------------------------------------------------------
# Exercise 3
# Optimal vs. Naive Diversification
#----------------------------------------------------------

# Split data into estimation sample and test sample
#-----------------------------------------------------------
data.est <- data[1:120, ]      # January 1993 to December 2002
data.test <- data[121:240, ]   # January 2003 to December 2012

# Use only the four risky asset classes
returns.est <- data.est[, 2:(ncol(data.est)-1)]
returns.test <- data.test[, 2:(ncol(data.test)-1)]


#----------------------------------------------------------
# Exercise 3a
# Reestimate input parameters and run Markowitz optimization
#----------------------------------------------------------

# Number of assets
N <- length(returns.est[1,])

# Estimate vector of mean returns
#-----------------------------------------------------------
MU <- rep(NA, N)

for(i in 1:length(returns.est[1,])){
  MU[i] <- round(mean(returns.est[,i], na.rm = TRUE) * 12, 4)
}
MU


# Estimate variance-covariance matrix
#-----------------------------------------------------------
SIGMA <- cov(returns.est, use = "complete.obs") * 12

VAR <- diag(SIGMA)
SD <- sqrt(VAR)

SIGMA
SD


##################
# Markowitz model
##################

#### Prepare objects for optimization function ####
#-----------------------------------------------------------
n <- 2 + N
meq <- 2
dvec <- rep(0, N)

Amat <- matrix(NA, N, n)
Amat[,1] <- MU
Amat[,2] <- rep(1, N)
for(i in 1:N){
  constraint <- rep(0, N)
  constraint[i] <- 1
  Amat[,i+2] <- constraint
}


# Set range of portfolio means along which to minimize variance
#-----------------------------------------------------------
stepsize <- 0.0001
mu_bar <- seq(min(MU), max(MU), by = stepsize)
mu_bar <- mu_bar[mu_bar < max(MU)]


# Create objects to store efficient frontier
#-----------------------------------------------------------
w <- matrix(NA, N, length(mu_bar))
mu <- rep(NA, length(mu_bar))
sigma <- rep(NA, length(mu_bar))


#### Optimization for Efficient Frontier ####
#-----------------------------------------------------------
for(i in 1:length(mu_bar)){
  w[,i] <- solve.QP(Dmat = SIGMA,
                    dvec = dvec,
                    Amat = Amat,
                    bvec = c(mu_bar[i], 1, rep(0, N)),
                    meq = meq)$solution
  mu[i] <- t(w[,i]) %*% MU
  sigma[i] <- sqrt(t(w[,i]) %*% SIGMA %*% w[,i])
}


# Find minimum variance portfolio
#-----------------------------------------------------------
min_var <- which.min(sigma)

# Weights of the minimum variance portfolio
w[, min_var]


#----------------------------------------------------------
# Exercise 3b
# 1/N strategy in the estimation sample
#----------------------------------------------------------

# Equal weights for all four asset classes
w.eq <- rep(1/N, N)

# Compute return time series of the 1/N portfolio
ret.eq <- rep(NA, nrow(returns.est))
for(i in 1:nrow(returns.est)){
  ret.eq[i] <- t(w.eq) %*% as.numeric(returns.est[i, ])
}

# Sample mean and sample standard deviation of historical returns
mean.eq <- mean(ret.eq)
sd.eq <- sd(ret.eq)

mean.eq
sd.eq


# Annualize 1/N portfolio statistics for the plot
#-----------------------------------------------------------
mu.eq <- mean.eq * 12
sigma.eq <- sd.eq * sqrt(12)

mu.eq
sigma.eq


#----------------------------------------------------------
# Exercise 3c
# Find efficient portfolio with same standard deviation as 1/N
#----------------------------------------------------------

same_sd <- which.min(abs(sigma.eq - sigma))

# Mean, standard deviation, and weights of the corresponding efficient portfolio
mu[same_sd]
sigma[same_sd]
w[, same_sd]


#########################################
# Plot the efficient frontier
#########################################

#-----------------------------------------------------------
plot(sigma, mu, type = "l",
     lwd = 6, xlim = c(0, max(c(sigma, SD, sigma.eq))*1.1),
     ylim = c(0, max(c(mu, MU, mu.eq))*1.1),
     xaxs = "i", yaxs = "i", col = "darkgreen",
     main = "Efficient Frontier",
     xlab = "Return Standard Deviation", ylab = "Expected Return",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)

# Asset class portfolios
points(SD, MU, pch = 16, cex = 1, lwd = 2, col = "black")

# Minimum variance portfolio
points(sigma[min_var], mu[min_var], pch = 19, cex = 1.5, lwd = 2,
       col = "black", bg = "darkgrey")

# 1/N portfolio
points(sigma.eq, mu.eq, pch = 17, cex = 1.5, lwd = 2, col = "blue")

# Efficient portfolio with same standard deviation as 1/N
points(sigma[same_sd], mu[same_sd], pch = 15, cex = 1.5, lwd = 2, col = "red")

# Labels for the four asset classes
text(SD, MU, labels = colnames(returns.est), pos = 4, cex = 0.8)

legend("bottomright",
       c("Efficient Frontier",
         "Minimum Variance Portfolio",
         "1/N Portfolio",
         "Efficient Portfolio with same SD as 1/N"),
       lty = c(1, NA, NA, NA),
       pch = c(NA, 19, 17, 15),
       lwd = c(6, NA, NA, NA),
       col = c("darkgreen", "black", "blue", "red"),
       cex = 1.1)


#----------------------------------------------------------
# Exercise 3d
# Out-of-sample Sharpe Ratios in the test period
#----------------------------------------------------------

# Risk-free rate in the test sample
#-----------------------------------------------------------
rf.test <- data.test[, ncol(data.test)]


# Compute realized returns of the 1/N portfolio in the test sample
#-----------------------------------------------------------
ret.eq.test <- rep(NA, nrow(returns.test))

for(i in 1:nrow(returns.test)){
  ret.eq.test[i] <- t(w.eq) %*% as.numeric(returns.test[i, ])
}

# Check realized returns
head(ret.eq.test)


# Compute realized returns of the efficient portfolio in the test sample
#-----------------------------------------------------------
ret.eff.test <- rep(NA, nrow(returns.test))

for(i in 1:nrow(returns.test)){
  ret.eff.test[i] <- t(w[, same_sd]) %*% as.numeric(returns.test[i, ])
}

# Check realized returns
head(ret.eff.test)


# Compute out-of-sample excess returns
#-----------------------------------------------------------
excess.eq <- ret.eq.test - rf.test
excess.eff <- ret.eff.test - rf.test

# Check excess returns
head(excess.eq)
head(excess.eff)


# Compute mean of out-of-sample excess returns
#-----------------------------------------------------------
mu_hat_e_eq <- mean(excess.eq)
mu_hat_e_eff <- mean(excess.eff)

mu_hat_e_eq
mu_hat_e_eff


# Compute standard deviation of out-of-sample excess returns
#-----------------------------------------------------------
sigma_hat_e_eq <- sd(excess.eq)
sigma_hat_e_eff <- sd(excess.eff)

sigma_hat_e_eq
sigma_hat_e_eff


# Compute out-of-sample Sharpe Ratios
#-----------------------------------------------------------
SR.eq <- mu_hat_e_eq / sigma_hat_e_eq
SR.eff <- mu_hat_e_eff / sigma_hat_e_eff

SR.eq
SR.eff


# Report: Out-of-sample Sharpe Ratios
# -----------------------------------------------------------
# 1/N Portfolio Sharpe Ratio:       0.1918
# Efficient Portfolio Sharpe Ratio:  0.1866
#
# Observation: The naïve 1/N portfolio achieves a slightly higher
# out-of-sample Sharpe Ratio (0.192 vs. 0.187) than the Markowitz-
# optimized efficient portfolio, despite the latter dominating the
# 1/N portfolio in the estimation sample (higher expected return
# for the same level of risk).
#
# Explanation: This result is driven by estimation risk. The
# Markowitz model treats estimated means, variances and
# covariances as if they were the true population parameters.
# In practice, especially expected returns are estimated with
# considerable error from historical data. The optimizer exploits
# these estimation errors, tilting weights heavily toward assets
# that happened to perform well in the estimation period (here:
# ~58% in corporate bonds, ~21% each in stocks and treasuries,
# 0% in real estate). When the market environment changes in
# the test period (2003-2012), these concentrated bets do not
# pay off as expected and the out-of-sample risk increases
# (sigma_hat = 0.0170 vs. 0.0144 for the 1/N portfolio).
# The 1/N strategy by contrast, is parameter-free and therefore
# immune to estimation error. Its equal weighting provides
# broader diversification, which proves more robust when the
# future deviates from the past. This is consistent with the 
# broader insight that, when the number of assets is small 
# and the estimation window is limited, naïve diversification 
# can match or even outperform mean-variance optimization 
# on a risk-adjusted basis out of sample.
