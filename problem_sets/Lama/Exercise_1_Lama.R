
####################################################
# Exercise 1
####################################################

################################################################################
# Exercise 1: Simulating a GBM
#-------------------------------------------------------------------------------
# adjust this accoridng to your working dir
#setwd("P:/HSG/SMFI/SMFI-main/problem_sets/Lama")

#----------------------------------------------------------------
# Exercise 1 a)
S0 <- 100                         # Current stock price
mu <- 0.08                        # Mean return rate 8%
sigma <- 0.2                      # Volatility of returns 20%
T <- 1                            # Price path over 1 year
s <- 250                          # Trading Days
delta_t <- rep(1/s,T*s)           # Daily time steps
t <- c(0,cumsum(delta_t))         # Time line including start date
paths <- 1                        # Paths
set.seed(13)                      # Random number generator

#----------------------------------------------------------------
delta_W <- sqrt(delta_t) * rnorm(n=T * s, mean = 0, sd=1) # Vector of path increments
Rt <- ((mu * delta_t)+(sigma*delta_W))                    # Return path

St <- c(S0, S0*cumprod(1+Rt))                             # Convert into price path



par(mfrow = c(2, 1), mar = c(4, 5, 2, 1))

plot(1:s, Rt, type = "l", col = "steelblue", lwd = 1,
     xlab = "Trading Day", ylab = "Daily Return",
     cex.axis = 1.5, cex.lab = 1.5,
     main = "Simulated Daily Returns")
abline(h = 0, lty = 2, col = "grey")

plot(t, St, type = "l", col = "darkgreen", lwd = 3,
     xlab = "Time t", ylab = "Stock Price",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,
     xaxs = "i", yaxs = "i", ylim = c(80, 120))

par(mfrow = c(1, 1))

################################################################################
# Exercise 1 b)
#-------------------------------------------------------------------------------
dar <- mean(Rt)               # Average daily returns
annualized_mu <- 250*dar      # Annulized Returns
print(annualized_mu)

# Report: The annualized sample mean differs from the true mu = 8% due to
# sampling variability. A single path consists of only 250 daily returns,
# which is a finite sample. The sample mean is an estimator of mu,
# but any single realization will deviate from the true value due to randomness.
# This deviation shrinks as the sample size increases.

################################################################################
# Exercise 1 c)
#-------------------------------------------------------------------------------
N <- 10000      # Number of paths
set.seed(13)

# Generate a Matrix
delta_W_matrix <- sqrt(1/s) * matrix(rnorm(N * s, mean = 0, sd = 1), nrow = N, ncol = s)

# Calculate Returns for all paths
Rt_matrix <- (mu * (1/s)) + (sigma * delta_W_matrix)

# Calculate the mean
dar_MC <- mean(Rt_matrix)

# Annulize the mean
annualized_mu_MC <- 250 * dar_MC

# Print the mean
print(annualized_mu_MC)

# Report: The Mean for 10'000 Paths is much closer to the true mu.
# It's because of the law of large numbers. It has to approx to 8%.

################################################################################
# Exercise 1 d)
#-------------------------------------------------------------------------------
T <- 3                            # Maturity T = 3 years
S0 <- 100                         # Current stock price
mu <- 0.08                        # Mean return rate 8%
sigma <- 0.2                      # Volatility of returns 20%
paths <- 10000
set.seed(13)

# Simulate W_T. It's normally distributed with mean 0 and variance T
W_T <- sqrt(T) * rnorm(paths, mean = 0, sd = 1)

# Compute 10'000 terminal prices using exact solution
S_T <- S0 * exp((mu - (sigma^2) / 2) * T + sigma * W_T)

# Plot Histogram
hist(S_T, breaks = 50, col = "steelblue", 
     main = "Distribution of Stock Price at T=3", 
     xlab = "Stock Price S_T", ylab = "Frequency")

# Compare sample average and expectation
sample_avg <- mean(S_T)
theo_exp <- S0 * exp(mu * T)

print(paste("Sample Average:", round(sample_avg, 4)))
print(paste("Theoretical Expectation:", round(theo_exp, 4)))

# Report:
# Shape: The histogram is right-skewed, which is expected because S_T follows
# a lognormal distribution. Prices cannot fall below zero but have no upper 
# bound, producing a right tail.
#
# Comparison: The sample average is close to the theoretical expectation
# E[S_T] = S0 * exp(mu * T) due to the law of large numbers. With N = 10'000
# simulations, the randomness averages out. The small remaining difference
# is finite-sample noise, which would shrink further as N increases.
