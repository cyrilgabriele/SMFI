############################
# Load packages
############################
library(YieldCurve)
set.seed(13)
############################
# Prepare data
############################
# Import interest rate data 
#----------------------------------------------------------------
yc.data <- read.csv2("YieldCurve.csv", header = TRUE, stringsAsFactors=FALSE)
str(yc.data)   # Column names show the yield curves available in the file
head(yc.data)
yields <- as.numeric(yc.data[,2])         # Extract the chosen column
#----------------------------------------------------------------
#----------------------------------------------------------------
# EXERCISE 1
# a) BEGIN
###########################################
# Nelson-Siegel-Svensson yield curve model
###########################################
# Create maturity vector and fit model
#----------------------------------------------------------------
maturity <- c(3/12,6/12,9/12, seq(1,30,1)) # Create maturity vector for fitting
NSSmodel <- Svensson(yields, maturity)     # Fit Nelson-Siegel-Svensson model
# Note: the function estimates tau1 and tau2, 
# which are the reciprocals of lambda1 and lambda2
#----------------------------------------------------------------
# Store parameter values
#----------------------------------------------------------------
b0 <- NSSmodel[1]
b1 <- NSSmodel[2]
b2 <- NSSmodel[3]
b3 <- NSSmodel[4]
lambda1 <- 1/NSSmodel[5]
lambda2 <- 1/NSSmodel[6]

cat(sprintf("b0      : %.6f\n", b0))
cat(sprintf("b1      : %.6f\n", b1))
cat(sprintf("b2      : %.6f\n", b2))
cat(sprintf("b3      : %.6f\n", b3))
cat(sprintf("lambda1      : %.6f\n", lambda1))
cat(sprintf("lambda2      : %.6f\n", lambda2))
# a) END
#----------------------------------------------------------------
# b) BEGIN 
# Implement the Nelson-Siegel-Svensson model for the spot yield curve
#----------------------------------------------------------------
# Note: fit is done on yields in percent (numerical stability of Svensson),
# but y0() returns the yield in DECIMAL form (e.g. 0.025 for 2.5%) so that
# downstream formulas (Hull-White, Black) can use it directly without
# additional /100 conversions.
y0 <- function(T){
  F1 <- (1-exp(-lambda1*T))/(lambda1*T)
  F2 <- (1-exp(-lambda1*T))/(lambda1*T) - exp(-lambda1*T)
  F3 <- (1-exp(-lambda2*T))/(lambda2*T) - exp(-lambda2*T)
  (b0 + b1*F1 + b2*F2 + b3*F3)/100}


# Plot yield curve against raw yields
#----------------------------------------------------------------
# Raw yields are in percent (as stored in the CSV); y0() returns decimal,
# so we multiply by 100 to display on the same axis.
#png(file="EuroYieldCurve.png",width=2000, height=1000, res = 200)
plot(maturity, yields, type = "p", cex = 1.25,
     xaxs = "i", yaxs = "i", pch = 21, col = "black", bg = "darkgreen",
     main = "Euro Area Yield Curve",
     xlim = c(0,31), ylim = c(0,5),
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,
     xlab = "Maturity (in years)",
     ylab = "Yield in % p.a.")
lines(maturity, y0(maturity)*100, col = "darkgreen", lwd = 2)
#abline(h = 0, lwd = 3, lty = 2)
#dev.off()
#----------------------------------------------------------------
# b) END
#----------------------------------------------------------------
# c) BEGIN 
f0 <- function(T1,T2){(y0(T2)*T2-y0(T1)*T1)/(T2-T1)}

T_1 <- 1
T_2 <- 2:30
f0_curve <- f0(T_1, T_2)

plot(T_2, f0_curve*100, type = "l", lwd = 2, col = "darkgreen",
     xaxs = "i", yaxs = "i",
     main = "Forward Rate Curve from T1 = 1",
     xlim = c(1,30), ylim = c(0,5),
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,
     xlab = "T2 (in years)",
     ylab = "Forward rate in % p.a.")
points(T_2, f0_curve*100, pch = 21, col = "black", bg = "darkgreen")
# c) END 
#----------------------------------------------------------------
# d) BEGIN 
# Returns DECIMAL form (divide by 100 at the end), consistent with y0().
f0inst <- function(T){
  F1 <- exp(-lambda1*T)
  F2 <- exp(-lambda1*T)*lambda1*T
  F3 <- exp(-lambda2*T)*lambda2*T
  (b0 + b1*F1 + b2*F2 + b3*F3)/100}

T_2 <- 2:30

f0inst_curve <- f0inst(T_2)

plot(T_2, f0inst_curve*100, type = "l", lwd = 2, col = "darkgreen",
     xaxs = "i", yaxs = "i",
     main = "Instantaneous Forward Rate Curve",
     xlim = c(1,30), ylim = c(0,5),
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,
     xlab = "Maturity (in years)",
     ylab = "Instantaneous forward rate in % p.a.")
points(T_2, f0inst_curve*100, pch = 21, col = "black", bg = "darkgreen")
# d) END 
#----------------------------------------------------------------
#----------------------------------------------------------------
# EXERCISE 2
# What we reuse from Ex. 1: 
# y0()
# f0inst()
# df0inst_dt() => derived from f0inst()
############################
# The Hull-White model
############################
# Input
# ---------------------------------------------------------------
# Note on units (decimal convention throughout):
#   - all rates (y0, f0, f0inst, r_t) are in DECIMAL form,
#     consistent with the convention adopted in files (i) and (ii).
#   - sigma is in decimal per sqrt(year), so it can be read as a percent:
#     sigma = 0.01 means 1% per year (100 bp/year of short-rate vol).
#   - alpha is NOT a rate or percentage. It is a mean-reversion SPEED in
#     units of 1/year. Its economic interpretation is the half-life of a
#     deviation from the forward curve: half-life = ln(2)/alpha.
#     So alpha = 0.10 means deviations halve every ln(2)/0.10 ~ 7 years.
paths <- 1000       # Number of simulated short rate paths
alpha <- 0.20       # Mean reversion speed 
sigma <- 0.01       # Short-rate volatility 
T     <- 10         # Time span for simulation (in years)
s     <- 250        # Time steps per year (e.g. 365 days)
dt    <- rep(1/s,T*s)                 # Vector of constant time steps
t     <- c(0,cumsum(dt))              # Timeline including start date (t=0)
r0    <- y0(1/s)                      # Current value of the short rate (decimal)
# ---------------------------------------------------------------
# Partial derivative of the instantaneous forward rate
# (Nelson-Siegel-Svensson model)
# ---------------------------------------------------------------
# Divided by 100 so that the derivative is consistent with f0inst's
# decimal convention.
df0inst_dt <- function(T){
  ( - lambda1*exp(-lambda1*T)*(b1 - b2 + b2*lambda1*T) +
      lambda2*exp(-lambda2*T)*(b3 - b3*lambda2*T) )/100}
# ---------------------------------------------------------------
# Simulate short rate paths -- LOOP version (transparent, slower)
# ---------------------------------------------------------------
# Every step of the Euler scheme is explicit: for each path j,
# at each time i, we (1) compute the mean-reversion level theta,
# (2) draw a Wiener increment, (3) update the short rate.
rt <- matrix(NA, paths, length(t))
rt[,1] <- r0
rt[1:10,1:15] # View structure of the matrix

# Each Wiener increment is drawn individually at every step.
for(j in 1:length(rt[,1])){ # Outer loop: go through rows (from one path to the next)
  for(i in 2:length(rt[1,])){ # Inner loop: go through columns (simulate each path)
    # Time-varying mean reversion level
    theta <- df0inst_dt(t[i-1]) + alpha*f0inst(t[i-1])
    # Increment of Wiener process (for one step of length dt)
    dW <- rnorm(n = 1, mean = 0, sd = sqrt(dt[1]))
    # Levels of the short rate in the path
    rt[j,i] <- rt[j,i-1] + (theta - alpha*rt[j,i-1])*dt[1] + sigma*dW
  }
}
# Plot short rate paths
# ---------------------------------------------------------------
# Paths are in decimal; multiply by 100 for percent-scale display.
# png(file="HullWhitePaths.png",width=2000, height=1000, res = 200)
short_rate_ylim <- range(c(rt, f0inst(t)), na.rm = TRUE)*100

plot(t, rt[1,]*100, type = "l", lwd = 1, xlim = range(t), ylim = short_rate_ylim,
     xaxs = "i", yaxs = "i", col = "darkgrey",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,
     xlab = "Time (in Years)",
     ylab = "Short Rate in % p.a.")
for(i in 2:paths){
  lines(t, rt[i,]*100, lwd = 1, col = "darkgrey")}
lines(t, rt[1,]*100, lwd = 1, col = "black")
lines(t, f0inst(t)*100, lwd = 3, col = "darkgreen")
legend("bottomright", c("Short Rate Paths", "Inst. Forward Curve"),
       lty = 1, lwd = c(1, 3), col = c("darkgrey", "darkgreen"),
       bty = "n", cex = 1.25, horiz = FALSE)
# dev.off()
# ---------------------------------------------------------------
# b) BEGIN
# Histograms of the short rate after 1, 5, and 10 years
# ---------------------------------------------------------------
years <- c(1, 5, 10)
year.idx <- years*s + 1
rt.selected <- rt[, year.idx]
colnames(rt.selected) <- paste0("t = ", years)

# used to ssave the old default plot setting
par.old <- par(no.readonly = TRUE)
# used to set the column split into 3 plots
par(mfrow = c(1, 3))

for(i in seq_along(years)){
  rt_i <- rt.selected[, i]
  mu_i <- mean(rt_i)
  sd_i <- sd(rt_i)
  
  hist(rt_i*100, breaks = 40, probability = TRUE,
       xaxs = "i", yaxs = "i", col = "lightgrey", border = "white",
       main = paste("Short rate in", years[i], "year(s)"),
       xlab = "Short Rate in % p.a.",
       ylab = "Density",
       cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)
  curve(dnorm(x, mean = mu_i*100, sd = sd_i*100),
        add = TRUE, col = "darkgreen", lwd = 3)
}

# restore old settings
par(par.old)

rt.summary <- data.frame(
  year = years,
  mean = colMeans(rt.selected),
  sd = apply(rt.selected, 2, sd)
)
print(rt.summary)
# The histograms are close to normal, as expected from the Gaussian
# Hull-White model. The dispersion increases with the simulation horizon.
# b) END
# ---------------------------------------------------------------
# c) BEGIN 
t_bond <- 3
T_bond <- 4
notional <- 100

# Lecture notes page 31 with SPECIAL CASE: t = 0 => simplified verson used below
P0 <- function(T){
  exp(-y0(T) * T)
}

# Lecture notes page 30
B <- function(t_bond, T_bond) {
  (1 - exp(-alpha * (T_bond - t_bond))) / alpha
}

# Lecture notes page 30
ln_A <- function(t_bond, T_bond) {
  log(P0(T_bond) / P0(t_bond)) +
    (B(t_bond, T_bond) * f0inst(t_bond)) -
    ((sigma^2 / (4 * alpha)) *
       (1 - exp(-2*alpha*t_bond)) *
       B(t_bond, T_bond)^2)
}

# Lecture notes page 30
ln_P_tilde_t <- function(t_bond, T_bond, r_t) {
  ln_A(t_bond, T_bond) -
    B(t_bond, T_bond) * r_t
}

# t_bond * s = 3 * 250
# + 1 => since column 1 is t = 0 and first simulated step is column 2
idx <- t_bond * s + 1
r_t <- rt[, idx]
ln_prices <- ln_P_tilde_t(t_bond, T_bond, r_t)

bond_prices <- notional * exp(ln_prices)
expected_bond_price <- mean(bond_prices)

cat(sprintf("expected_bond_price      : %.6f\n", expected_bond_price))
# c) END
# ---------------------------------------------------------------
# d) BEGIN 
tk <- 2 
tk_plus_1 <- 3 
delta_tk <- tk_plus_1 - tk

N <- 1000000
y_x <- 0.015

# Indices
idx_fix <- tk * s + 1          # time t = 2 (= tk)
idx_pay <- tk_plus_1 * s + 1   # time T = 3 (= tk_plus_1)

# Short rate at fixing date t_k = 2
r_tk <- rt[, idx_fix]

# Use your existing Hull-White zero-bond log-price function
ln_P_tk_tk1 <- ln_P_tilde_t(tk, tk_plus_1, r_tk)

# Recover the future zero-coupon yield y_tilde_tk(tk+1)
y_tilde_tk_tk1 <- -ln_P_tk_tk1 / delta_tk

# Caplet payoff at time t_{k+1} = 3
payoff_tk1 <- N * delta_tk * pmax(y_tilde_tk_tk1 - y_x, 0)

# Pathwise discount factor from 0 to 3 under Q*
DF_0_tk1 <- exp(
  -rowSums(rt[, 1:(idx_pay - 1), drop = FALSE]) * dt[1]
)

# Caplet price under traditional risk-neutral measure Q*
caplet_price <- mean(DF_0_tk1 * payoff_tk1)

cat(sprintf("caplet_price: %.6f\n", caplet_price))
# d) END
