############################################################
# Problem Set 3
# Stochastic Modelling in Finance and Insurance
# Adriano Lama
############################################################

#clear workspace
rm(list = ls())

#readable numbers
options(scipen = 999)

############################################################
# Load package
############################################################

library(YieldCurve)


############################################################
# Exercise 1: Nelson-Siegel-Svensson Model
############################################################

############################################################
# 1a) Import yield curve data and estimate NSS parameters
############################################################

yc.data <- read.csv2(
  "YieldCurve.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)

yields <- as.numeric(yc.data[, 2])              # yields are in percent
maturity <- c(3/12, 6/12, 9/12, seq(1, 30, 1)) # maturities in years

# Fit Nelson-Siegel-Svensson model
NSSmodel <- Svensson(yields, maturity)

# Store parameters
b0 <- NSSmodel[1]
b1 <- NSSmodel[2]
b2 <- NSSmodel[3]
b3 <- NSSmodel[4]

# YieldCurve::Svensson estimates tau1 and tau2.
# In the slides/code we use lambda1 = 1/tau1 and lambda2 = 1/tau2.
lambda1 <- 1 / NSSmodel[5]
lambda2 <- 1 / NSSmodel[6]

params <- data.frame(
  Parameter = c("beta0", "beta1", "beta2", "beta3", "lambda1", "lambda2"),
  Estimate = c(b0, b1, b2, b3, lambda1, lambda2)
)

print(params)

############################################################
# 1b) Spot yield curve function y0(T)
############################################################

# Important:
# The fitted yields are in percent.
# For all later financial formulas, y0(T) returns DECIMAL form.
# Example: 0.025 means 2.5%.

y0 <- function(T) {
  F1 <- (1 - exp(-lambda1 * T)) / (lambda1 * T)
  F2 <- (1 - exp(-lambda1 * T)) / (lambda1 * T) - exp(-lambda1 * T)
  F3 <- (1 - exp(-lambda2 * T)) / (lambda2 * T) - exp(-lambda2 * T)
  
  (b0 + b1 * F1 + b2 * F2 + b3 * F3) / 100
}

m <- seq(0.25, 30, 0.01)

plot(
  maturity, yields,
  type = "p",
  pch = 21,
  col = "black",
  bg = "darkgreen",
  xlim = c(0, 31),
  ylim = c(0, 5),
  xlab = "Maturity (in years)",
  ylab = "Yield in % p.a.",
  main = "Exercise 1b: Spot Yield Curve"
)

lines(m, y0(m) * 100, col = "darkgreen", lwd = 2)


############################################################
# 1c) Arbitrage-free forward rate function f0(T1,T2)
############################################################

# Formula:
# f0(T1,T2) = [y0(T2)*T2 - y0(T1)*T1] / (T2 - T1)

f0 <- function(T1, T2) {
  (y0(T2) * T2 - y0(T1) * T1) / (T2 - T1)
}

T1 <- 1
T2.grid <- seq(1.01, 30, 0.01)

plot(
  T2.grid, f0(T1, T2.grid) * 100,
  type = "l",
  lwd = 2,
  col = "darkgreen",
  xlim = c(1, 30),
  ylim = c(0, 5),
  xlab = "Termination date T2",
  ylab = "Forward rate in % p.a.",
  main = "Exercise 1c: Forward Curve with T1 = 1"
)


############################################################
# 1d) Instantaneous forward curve function f0inst(T)
############################################################

f0inst <- function(T) {
  F1 <- exp(-lambda1 * T)
  F2 <- exp(-lambda1 * T) * lambda1 * T
  F3 <- exp(-lambda2 * T) * lambda2 * T
  
  (b0 + b1 * F1 + b2 * F2 + b3 * F3) / 100
}

m.inst <- seq(0, 30, 0.01)

plot(
  m.inst, f0inst(m.inst) * 100,
  type = "l",
  lwd = 2,
  col = "darkgreen",
  xlim = c(0, 30),
  ylim = c(0, 5),
  xlab = "Maturity (in years)",
  ylab = "Instantaneous forward rate in % p.a.",
  main = "Exercise 1d: Instantaneous Forward Curve"
)


############################################################
# Exercise 2: Hull-White Model
############################################################

############################################################
# 2a) Simulate 1000 short-rate paths
############################################################

alpha <- 0.20      # mean reversion speed
sigma <- 0.01      # short-rate volatility
paths <- 1000      # number of paths
Tsim <- 10         # simulation horizon in years
steps <- 250       # trading days per year
dt <- 1 / steps

t <- seq(0, Tsim, by = dt)

# Initial short rate, close to time 0
r0 <- y0(dt)

# Derivative of instantaneous forward curve
df0inst_dt <- function(T) {
  (
    -lambda1 * exp(-lambda1 * T) * (b1 - b2 + b2 * lambda1 * T) +
      lambda2 * exp(-lambda2 * T) * (b3 - b3 * lambda2 * T)
  ) / 100
}

# Simulate Hull-White paths by Euler discretization:
# r_t = r_{t-1} + (theta_{t-1} - alpha*r_{t-1})*dt + sigma*dW

set.seed(1337) #For reproductionability

rt <- matrix(NA, nrow = paths, ncol = length(t))
rt[, 1] <- r0

dW <- matrix(
  rnorm(paths * (length(t) - 1), mean = 0, sd = sqrt(dt)),
  nrow = paths
)

theta <- df0inst_dt(t[-length(t)]) + alpha * f0inst(t[-length(t)])

for (i in 2:length(t)) {
  rt[, i] <- rt[, i - 1] +
    (theta[i - 1] - alpha * rt[, i - 1]) * dt +
    sigma * dW[, i - 1]
}

# Plot paths
matplot(
  t, t(rt) * 100,
  type = "l",
  lty = 1,
  col = "grey",
  xlim = c(0, 10),
  xlab = "Time (in years)",
  ylab = "Short rate in % p.a.",
  main = "Exercise 2a: Hull-White Short Rate Paths"
)

lines(t, f0inst(t) * 100, col = "darkgreen", lwd = 3)

legend(
  "bottomright",
  legend = c("Short rate paths", "Instantaneous forward curve"),
  col = c("grey", "darkgreen"),
  lty = 1,
  lwd = c(1, 3),
  bty = "n"
)


############################################################
# 2b) Histograms at t = 1, 5, 10 with normal density
############################################################

get_index <- function(time) {
  round(time / dt) + 1
}

hist_with_normal <- function(time) {
  r <- rt[, get_index(time)] * 100
  
  hist(
    r,
    probability = TRUE,
    breaks = 30,
    col = "lightgrey",
    border = "white",
    main = paste("t =", time),
    xlab = "Short rate in % p.a."
  )
  
  x <- seq(min(r), max(r), length.out = 200)
  lines(x, dnorm(x, mean = mean(r), sd = sd(r)), col = "darkgreen", lwd = 2)
}

par(mfrow = c(1, 3))
hist_with_normal(1)
hist_with_normal(5)
hist_with_normal(10)
par(mfrow = c(1, 1))

summary.rates <- data.frame(
  Time = c(1, 5, 10),
  Mean_percent = c(
    mean(rt[, get_index(1)] * 100),
    mean(rt[, get_index(5)] * 100),
    mean(rt[, get_index(10)] * 100)
  ),
  SD_percent = c(
    sd(rt[, get_index(1)] * 100),
    sd(rt[, get_index(5)] * 100),
    sd(rt[, get_index(10)] * 100)
  )
)

print(summary.rates)

cat("\nExercise 2b interpretation:\n")
cat("The simulated short rates are approximately normally distributed.\n")
cat("This is expected because the Hull-White model is Gaussian: the short rate is driven by normal Wiener increments.\n")
cat("The standard deviation increases over time, but mean reversion prevents it from exploding.\n\n")


############################################################
# Hull-White affine term structure function
############################################################

# This function returns the future continuously-compounded spot rate
# y_time(T) for all simulated paths.
# T is a calendar maturity date, not time-to-maturity.

yHW <- function(time, T) {
  P0T <- exp(-y0(T) * T)
  P0t <- exp(-y0(time) * time)
  
  B <- (1 - exp(-alpha * (T - time))) / alpha
  
  lnA <- log(P0T / P0t) +
    B * f0inst(time) -
    sigma^2 / (4 * alpha) * (1 - exp(-2 * alpha * time)) * B^2
  
  idx <- get_index(time)
  
  -lnA / (T - time) + B / (T - time) * rt[, idx]
}


############################################################
# 2c) Expected price of a zero-coupon bond
############################################################

# One-year zero-coupon bond in three years:
# time = 3, maturity date = 4, notional = 100

y_3_4 <- yHW(3, 4)

zcb_price_t3 <- 100 * exp(-y_3_4 * (4 - 3))

expected_zcb_price_t3 <- mean(zcb_price_t3)

cat("Exercise 2c:\n")
cat("Expected price at t = 3 of a one-year zero-coupon bond with notional 100:\n")
cat(round(expected_zcb_price_t3, 4), "\n\n")

############################################################
# 2d) Caplet price under traditional risk-neutral measure Q*
############################################################

# Caplet data from the exercise
tk <- 2
tk1 <- 3
delta_tk <- tk1 - tk
N <- 1000000
yX <- 0.015

# Future one-year rate fixed at t = 2 for period [2,3]
y_2_3 <- yHW(tk, tk1)

# Payoff at t = 3
caplet_payoff_t3 <- N * delta_tk * pmax(y_2_3 - yX, 0)

# Traditional risk-neutral valuation:
# discount pathwise by exp(- integral_0^3 r_s ds)
discount_0_3 <- exp(-rowSums(rt[, 1:(get_index(tk1) - 1)] * dt))

caplet_price <- mean(discount_0_3 * caplet_payoff_t3)
caplet_se <- sd(discount_0_3 * caplet_payoff_t3) / sqrt(paths)

cat("Exercise 2d:\n")
cat("Caplet price under Q*:\n")
cat(round(caplet_price, 4), "\n")
cat("Monte Carlo standard error:\n")
cat(round(caplet_se, 4), "\n\n")


############################################################
# Exercise 3: GenAI in Financial Engineering
############################################################

# ----------------------------------------------------------
# 1. Prompt Used
# ----------------------------------------------------------
# "Act as an expert in financial engineering. I need to design an 
# innovative interest rate derivative product to hedge against 
# interest rate risk. Please provide the following:
# 1. Product Description: An intuitive explanation of a unique 
#    interest rate derivative. Explain the real-world use case 
#    and who would buy it.
# 2. Payoff Function: The formal mathematical payoff function 
#    of the product.
# 3. R Code: The complete R code to price this product. Assume 
#    the Hull-White short rate paths are already simulated and 
#    stored in a matrix called rt, with paths = 1000 and time 
#    steps dt. Discount the payoff pathwise to t=0 using the 
#    traditional risk-neutral measure Q*. Ensure the code 
#    prints the final expected price and the Monte Carlo 
#    standard error. Keep the code clean and well-commented."

# ----------------------------------------------------------
# 2. Product Description: Peak Rate Protection (Lookback Cap)
# ----------------------------------------------------------
# Standard interest rate caps consist of a series of caplets that 
# pay out if the floating rate exceeds a strike rate on specific, 
# predetermined fixing dates. However, interest rates can be highly 
# volatile between these discrete dates. The "Peak Rate Protection" 
# is a path-dependent derivative that continuously monitors the 
# short rate over a specified risk period. It pays out based on 
# the absolute highest rate (the peak) achieved during that entire 
# timeframe, rather than looking at discrete fixing dates.
#
# Real-World Use Case & Target Buyer:
# This product is highly attractive to real estate developers, private 
# equity firms, or infrastructure funds that rely on rolling, short-term 
# commercial paper or continuous revolving credit facilities. Because 
# their debt is constantly rolling over, they are exposed to interest 
# rate spikes at virtually any point in time. A corporate treasurer 
# would buy this to gain absolute certainty that the maximum financing 
# cost over the next T years will not exceed a catastrophic threshold.

# ----------------------------------------------------------
# 3. Formal Payoff Function
# ----------------------------------------------------------
# The formal payoff of the Peak Rate Protection at maturity T is 
# calculated by finding the maximum realized short rate r_t over 
# the continuous time interval [0, T], and comparing it to the 
# pre-agreed strike rate (cap rate) yX.
#
# Payoff_T = N * max( max_{0 <= t <= T} (r_t) - yX, 0 )
#
# Where:
# - N is the notional amount.
# - T is the maturity of the option (the end of the observation period).
# - r_t is the stochastic short rate at time t.
# - yX is the strike rate (cap rate).

# ----------------------------------------------------------
# 4. R Code and Pricing
# ----------------------------------------------------------

# 1. Define Option Parameters
N <- 1000000          # Notional amount (1,000,000)
T_mat <- 3            # Maturity of the option in years
yX <- 0.04            # Strike rate (cap rate) at 4.0% (decimal form)

# 2. Identify the time index for Maturity T
# This ensures we only look at the path up to the option's maturity
idx_T <- round(T_mat / dt) + 1

# 3. Calculate Path-Dependent Payoff
# Extract the simulated short rates from t=0 to t=T for all paths
rt_period <- rt[, 1:idx_T]

# Find the maximum short rate achieved on EACH path during the period
# apply(matrix, 1, max) applies the max function across the rows (paths)
r_max_per_path <- apply(rt_period, 1, max)

# Calculate the payoff at time T for each path
payoff_T <- N * pmax(r_max_per_path - yX, 0)

# 4. Pathwise Discounting under Q*
# Under traditional risk-neutral valuation, we discount each path 
# by exp(- integral_0^T r_s ds)
# We sum the short rates along the path and multiply by the step size dt
discount_factors <- exp(-rowSums(rt[, 1:(idx_T - 1)] * dt))

# 5. Calculate Present Value
# Discount the payoff of each path back to t=0
pv_payoffs <- discount_factors * payoff_T

# 6. Compute Expected Price and Standard Error
expected_price <- mean(pv_payoffs)
mc_standard_error <- sd(pv_payoffs) / sqrt(paths)

# 7. Output Results
cat("\n====================================================\n")
cat(" Exercise 3: Peak Rate Protection (Lookback Cap) Pricing \n")
cat("====================================================\n")
cat("Notional:               ", format(N, scientific = FALSE), "\n")
cat("Maturity:               ", T_mat, "Years\n")
cat("Strike Rate:            ", yX * 100, "%\n")
cat("----------------------------------------------------\n")
cat("Expected Price (Q*):    ", round(expected_price, 4), "\n")
cat("Monte Carlo SE:         ", round(mc_standard_error, 4), "\n")
cat("====================================================\n\n")

