############################################################
# Exercise 4
# Custom derivative product under GBM and risk-neutral pricing
#
# Product: Protected Target-Range Bull-Spread Note
#
# Payoff:
# X_T = G
#       + lambda * [ (S_T - K)^+ - (S_T - U)^+ ]
#       + B * 1_{L <= S_T <= U}
#
# Model under Q:
# dS_t = r S_t dt + sigma S_t dW_t^Q
#
# This script computes:
# 1) the analytical Black-Scholes price
# 2) a Monte Carlo verification under Q
############################################################

rm(list = ls())

############################
# 1. Parameters
############################

# Market parameters
S0    <- 100      # initial stock price
r     <- 0.02     # risk-free rate
sigma <- 0.22     # volatility
T     <- 2        # maturity in years

# Product parameters
G      <- 100     # guaranteed amount at maturity
lambda <- 0.70    # participation rate in the bull spread
K      <- 100     # lower strike
U      <- 130     # upper strike / cap
L      <- 95      # lower bound of target range
B      <- 6       # fixed bonus if ST in [L, U]

############################
# 2. Black-Scholes functions
############################

# Black-Scholes call price
bs_call <- function(S0, K, r, sigma, T) {
  d1 <- (log(S0 / K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  S0 * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
}

# Risk-neutral probability Q(ST > x)
# Under GBM with no dividends:
# Q(ST > x) = Phi(d2(x))
d2_fun <- function(S0, x, r, sigma, T) {
  (log(S0 / x) + (r - 0.5 * sigma^2) * T) / (sigma * sqrt(T))
}

# Present value of a digital corridor paying 1 if L <= ST <= U
range_digital_pv <- function(S0, L, U, r, sigma, T) {
  d2L <- d2_fun(S0, L, r, sigma, T)
  d2U <- d2_fun(S0, U, r, sigma, T)
  exp(-r * T) * (pnorm(d2L) - pnorm(d2U))
}

############################
# 3. Analytical price
############################

# Bond component
bond_pv <- G * exp(-r * T)

# Bull spread component
call_K <- bs_call(S0, K, r, sigma, T)
call_U <- bs_call(S0, U, r, sigma, T)
spread_pv <- lambda * (call_K - call_U)

# Bonus component
bonus_pv <- B * range_digital_pv(S0, L, U, r, sigma, T)

# Total analytical price
V0_analytic <- bond_pv + spread_pv + bonus_pv

############################
# 4. Output analytical results
############################

cat("========================================\n")
cat("ANALYTICAL PRICING RESULTS\n")
cat("========================================\n")
cat(sprintf("Bond component PV      : %.6f\n", bond_pv))
cat(sprintf("Call(K=%.2f)           : %.6f\n", K, call_K))
cat(sprintf("Call(U=%.2f)           : %.6f\n", U, call_U))
cat(sprintf("Bull spread PV         : %.6f\n", spread_pv))
cat(sprintf("Range bonus PV         : %.6f\n", bonus_pv))
cat(sprintf("Total analytical price : %.6f\n\n", V0_analytic))

############################
# 5. Monte Carlo verification under Q
############################

set.seed(13)

n_sim <- 1000000
Z <- rnorm(n_sim)

# Exact terminal stock price under Q
ST <- S0 * exp((r - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z)

# Payoff at maturity
payoff_T <- G +
  lambda * (pmax(ST - K, 0) - pmax(ST - U, 0)) +
  B * as.numeric(ST >= L & ST <= U)

# Discounted payoff
pv <- exp(-r * T) * payoff_T

# Monte Carlo estimator and standard error
V0_mc <- mean(pv)
SE_mc <- sd(pv) / sqrt(n_sim)
CI_low <- V0_mc - 1.96 * SE_mc
CI_high <- V0_mc + 1.96 * SE_mc

############################
# 6. Output Monte Carlo results
############################

cat("========================================\n")
cat("MONTE CARLO VERIFICATION\n")
cat("========================================\n")
cat(sprintf("Monte Carlo price      : %.6f\n", V0_mc))
cat(sprintf("Standard error         : %.6f\n", SE_mc))
cat(sprintf("95%% confidence interval: [%.6f, %.6f]\n", CI_low, CI_high))
cat(sprintf("Analytical - MC diff   : %.6f\n", V0_analytic - V0_mc))
cat("========================================\n")