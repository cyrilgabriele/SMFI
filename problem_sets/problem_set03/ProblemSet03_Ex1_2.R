############################
# Load packages
############################
library(YieldCurve)

set.seed(13)

############################
# Prepare data
############################

yc.data <- read.csv2("YieldCurve.csv", header = TRUE, stringsAsFactors = FALSE)

str(yc.data)
head(yc.data)

# Extract chosen yield curve column
# Adjust column index if your YieldCurve.csv uses another column.
yields <- as.numeric(yc.data[, 2])

# ================================================================
# EXERCISE 1
# Nelson-Siegel-Svensson Model
# ================================================================
# ---------------------------------------------------------------
# Exercise 1 a) Estimate NSS model
# ---------------------------------------------------------------
maturity <- c(3/12, 6/12, 9/12, seq(1, 30, 1))

NSSmodel <- Svensson(yields, maturity)

# Svensson() returns tau1 and tau2, which are reciprocals of lambda1/lambda2
b0 <- NSSmodel[1]
b1 <- NSSmodel[2]
b2 <- NSSmodel[3]
b3 <- NSSmodel[4]

lambda1 <- 1 / NSSmodel[5]
lambda2 <- 1 / NSSmodel[6]

cat("\nExercise 1 a) NSS parameter estimates\n")
cat(sprintf("b0      : %.6f\n", b0))
cat(sprintf("b1      : %.6f\n", b1))
cat(sprintf("b2      : %.6f\n", b2))
cat(sprintf("b3      : %.6f\n", b3))
cat(sprintf("lambda1 : %.6f\n", lambda1))
cat(sprintf("lambda2 : %.6f\n", lambda2))

# ---------------------------------------------------------------
# Exercise 1 b) Spot yield curve function y0(T)
# ---------------------------------------------------------------
# y0(T) returns DECIMAL rates, e.g. 0.025 for 2.5%.
y0 <- function(T){
  F1 <- (1 - exp(-lambda1 * T)) / (lambda1 * T)
  F2 <- (1 - exp(-lambda1 * T)) / (lambda1 * T) - exp(-lambda1 * T)
  F3 <- (1 - exp(-lambda2 * T)) / (lambda2 * T) - exp(-lambda2 * T)
  
  (b0 + b1 * F1 + b2 * F2 + b3 * F3) / 100
}

plot(maturity, yields,
     type = "p",
     cex = 1.25,
     xaxs = "i",
     yaxs = "i",
     pch = 21,
     col = "black",
     bg = "darkgreen",
     main = "Spot Yield Curve",
     xlim = c(0, 31),
     ylim = c(0, 5),
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main = 1.5,
     xlab = "Maturity in years",
     ylab = "Yield in % p.a.")

lines(maturity, y0(maturity) * 100,
      col = "darkgreen",
      lwd = 2)

# ---------------------------------------------------------------
# Exercise 1 c) Arbitrage-free forward rate function f0(T1,T2)
# ---------------------------------------------------------------
# f0(T1,T2) returns DECIMAL rates.
f0 <- function(T1, T2){
  (y0(T2) * T2 - y0(T1) * T1) / (T2 - T1)
}

T_1 <- 1
T_2 <- seq(T_1 + 0.01, 30, by = 0.01)

f0_curve <- f0(T_1, T_2)

plot(T_2, f0_curve * 100,
     type = "l",
     lwd = 2,
     col = "darkgreen",
     xaxs = "i",
     yaxs = "i",
     main = "Forward Rate Curve from T1 = 1",
     xlim = c(1, 30),
     ylim = c(0, 5),
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main = 1.5,
     xlab = "Termination date T2 in years",
     ylab = "Forward rate in % p.a.")

# ---------------------------------------------------------------
# Exercise 1 d) Instantaneous forward rate function f0inst(T)
# ---------------------------------------------------------------

# f0inst(T) returns DECIMAL rates.
f0inst <- function(T){
  F1 <- exp(-lambda1 * T)
  F2 <- exp(-lambda1 * T) * lambda1 * T
  F3 <- exp(-lambda2 * T) * lambda2 * T
  
  (b0 + b1 * F1 + b2 * F2 + b3 * F3) / 100
}

T_grid <- seq(0.25, 30, by = 0.01)

f0inst_curve <- f0inst(T_grid)

plot(T_grid, f0inst_curve * 100,
     type = "l",
     lwd = 2,
     col = "darkgreen",
     xaxs = "i",
     yaxs = "i",
     main = "Instantaneous Forward Rate Curve",
     xlim = c(0.25, 30),
     ylim = c(0, 5),
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main = 1.5,
     xlab = "Effective date T in years",
     ylab = "Instantaneous forward rate in % p.a.")

# ================================================================
# EXERCISE 2
# Hull-White Model
# ================================================================
# ---------------------------------------------------------------
# Exercise 2 a) Simulate Hull-White short-rate paths
# ---------------------------------------------------------------
paths <- 1000
alpha <- 0.20
sigma <- 0.01
T <- 10
s <- 250

dt <- rep(1 / s, T * s)
t <- c(0, cumsum(dt))

# Initial short rate: very short maturity spot rate
r0 <- y0(1 / s)

# Derivative of instantaneous forward curve
df0inst_dt <- function(T){
  (
    -lambda1 * exp(-lambda1 * T) * (b1 - b2 + b2 * lambda1 * T) +
      lambda2 * exp(-lambda2 * T) * (b3 - b3 * lambda2 * T)
  ) / 100
}

# Simulate short-rate paths -- LOOP version
rt <- matrix(NA, nrow = paths, ncol = length(t))
rt[, 1] <- r0

for(j in 1:paths){
  for(i in 2:length(t)){
    
    # Time-varying mean-reversion level
    theta <- df0inst_dt(t[i - 1]) + alpha * f0inst(t[i - 1])
    
    # Wiener increment: detla_W ~ N(0, delta_t)
    dW <- rnorm(n = 1, mean = 0, sd = sqrt(dt[1]))
    
    # Euler discretization of Hull-White:
    # r_t = r_{t-dt} + (theta_t - alpha*r_{t-dt})*dt + sigma*dW
    rt[j, i] <- rt[j, i - 1] +
      (theta - alpha * rt[j, i - 1]) * dt[1] +
      sigma * dW
  }
}

# Plot simulated short-rate paths with instantaneous forward curve
short_rate_ylim <- range(c(rt, f0inst(t)), na.rm = TRUE) * 100

plot(t, rt[1, ] * 100,
     type = "l",
     lwd = 1,
     xlim = range(t),
     ylim = short_rate_ylim,
     xaxs = "i",
     yaxs = "i",
     col = "darkgrey",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main = 1.5,
     main = "Hull-White Short-Rate Paths",
     xlab = "Time in years",
     ylab = "Short rate in % p.a.")

for(j in 2:paths){
  lines(t, rt[j, ] * 100,
        lwd = 1,
        col = "darkgrey")
}

lines(t, rt[1, ] * 100,
      lwd = 1,
      col = "black")

lines(t, f0inst(t) * 100,
      lwd = 3,
      col = "darkgreen")

legend("bottomright",
       legend = c("Short-rate paths", "Instantaneous forward curve"),
       lty = 1,
       lwd = c(1, 3),
       col = c("darkgrey", "darkgreen"),
       bty = "n",
       cex = 1.25)


# ---------------------------------------------------------------
# Exercise 2 b) Histograms after 1, 5, and 10 years
# ---------------------------------------------------------------
years <- c(1, 5, 10)

year.idx <- years * s + 1
rt.selected <- rt[, year.idx]

colnames(rt.selected) <- paste0("t = ", years)

par.old <- par(no.readonly = TRUE)
par(mfrow = c(1, 3))

for(i in seq_along(years)){
  
  rt_i <- rt.selected[, i]
  
  mu_i <- mean(rt_i)
  sd_i <- sd(rt_i)
  
  hist(rt_i * 100,
       breaks = 40,
       probability = TRUE,
       xaxs = "i",
       yaxs = "i",
       col = "lightgrey",
       border = "white",
       main = paste("Short rate after", years[i], "year(s)"),
       xlab = "Short rate in % p.a.",
       ylab = "Density",
       cex.axis = 1.2,
       cex.lab = 1.2,
       cex.main = 1.2)
  
  curve(dnorm(x, mean = mu_i * 100, sd = sd_i * 100),
        add = TRUE,
        col = "darkgreen",
        lwd = 3)
}

par(par.old)

rt.summary <- data.frame(
  year = years,
  mean_decimal = colMeans(rt.selected),
  sd_decimal = apply(rt.selected, 2, sd),
  mean_percent = colMeans(rt.selected) * 100,
  sd_percent = apply(rt.selected, 2, sd) * 100
)

cat("\nExercise 2 b) Short-rate summary\n")
print(rt.summary)

cat("
The simulated short rate is approximately normally distributed at all three horizons, 
which is consistent with the Gaussian Hull-White / Ornstein-Uhlenbeck structure. 
The normal density fits the histogram well. The dispersion increases with the horizon 
because uncertainty accumulates over time, although mean reversion prevents the variance 
from growing without bound. Negative short rates occur in some paths, which is a known 
feature of the Gaussian Hull-White model.
")
# ---------------------------------------------------------------
# Exercise 2 c) Expected future zero-coupon bond price
# ---------------------------------------------------------------
# Bond is observed in 3 years and has maturity of 1 year.
# Therefore: t = 3, T = 4.
t_bond <- 3
T_bond <- 4
notional <- 100

# Initial zero-coupon bond price from today's curve
P0 <- function(T){
  exp(-y0(T) * T)
}

# Hull-White B(t,T)
B <- function(t_bond, T_bond){
  (1 - exp(-alpha * (T_bond - t_bond))) / alpha
}

# Hull-White ln A(t,T)
ln_A <- function(t_bond, T_bond){
  log(P0(T_bond) / P0(t_bond)) +
    B(t_bond, T_bond) * f0inst(t_bond) -
    sigma^2 / (4 * alpha) *
    (1 - exp(-2 * alpha * t_bond)) *
    B(t_bond, T_bond)^2
}

# Hull-White log zero-coupon bond price
ln_P_tilde_t <- function(t_bond, T_bond, r_t){
  ln_A(t_bond, T_bond) -
    B(t_bond, T_bond) * r_t
}

idx_bond <- t_bond * s + 1

r_t_bond <- rt[, idx_bond]

ln_prices <- ln_P_tilde_t(t_bond, T_bond, r_t_bond)

bond_prices <- notional * exp(ln_prices)

expected_bond_price <- mean(bond_prices)

cat("\nExercise 2 c) Expected future zero-coupon bond price\n")
cat(sprintf("Expected bond price: %.6f\n", expected_bond_price))

# ---------------------------------------------------------------
# Exercise 2 d) Caplet price under traditional risk-neutral measure Q*
# ---------------------------------------------------------------
tk <- 2
tk_plus_1 <- 3

delta_tk <- tk_plus_1 - tk

N <- 1000000
y_x <- 0.015

idx_fix <- tk * s + 1
idx_pay <- tk_plus_1 * s + 1

# Simulated short rate at fixing date tk = 2
r_tk <- rt[, idx_fix]

# Compute log zero-coupon bond price ln P_tk(tk+1)
ln_P_tk_tk1 <- ln_P_tilde_t(tk, tk_plus_1, r_tk)

# Recover future zero-coupon yield y_tk(tk+1)
y_tilde_tk_tk1 <- -ln_P_tk_tk1 / delta_tk

# Caplet payoff at payment date tk+1
payoff_tk1 <- N * delta_tk * pmax(y_tilde_tk_tk1 - y_x, 0)

# Pathwise stochastic discount factor from 0 to tk+1 under Q*
DF_0_tk1 <- exp(
  -rowSums(rt[, 1:(idx_pay - 1), drop = FALSE]) * dt[1]
)

# Monte Carlo caplet price
caplet_price <- mean(DF_0_tk1 * payoff_tk1)

# Optional Monte Carlo standard error
caplet_price_se <- sd(DF_0_tk1 * payoff_tk1) / sqrt(paths)

cat("\nExercise 2 d) Caplet price under Q*\n")
cat(sprintf("Caplet price: %.6f\n", caplet_price))
cat(sprintf("Monte Carlo SE: %.6f\n", caplet_price_se))