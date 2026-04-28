
############################
# The Hull-White model
############################
# Need to load content of files:
# SMFI R Code Chapter 4 - i. Yield Curve.R
# SMFI R Code Chapter 4 - ii. Forward Curve.R
source("Chapter_4_01_Yield_Curve.R")
source("Chapter_4_02_Forward_Curve.R")
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

paths <- 300        # Number of simulated short rate paths
alpha <- 0.30       # Mean reversion speed 
sigma <- 0.002      # Short-rate volatility 
T     <- 20         # Time span for simulation (in years)
s     <- 365        # Time steps per year (e.g. 365 days)
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

set.seed(13)
rt <- matrix(NA, paths, length(t))
rt[,1] <- r0
rt[1:10,1:15] # View structure of the matrix

# Wrap the simulation in system.time() to record how long it takes.
# loop_time will be compared with vec_time below to quantify the speed-up.
loop_time <- system.time({

# Both dimensions are explicit loops: the outer one over paths, the inner one
# over time. Each Wiener increment is drawn individually at every step.
for(j in 1:length(rt[,1])){ # Outer loop: go through rows (from one path to the next)
  for(i in 2:length(rt[1,])){ # Inner loop: go through columns (simulate each path)
      # Time-varying mean reversion level
      theta <- df0inst_dt(t[i-1]) + alpha*f0inst(t[i-1])
      # Increment of Wiener process (for one step of length dt)
      dW    <- rnorm(n = 1, mean = 0, sd = sqrt(dt[1]))
      # Levels of the short rate in the path
      rt[j,i] <- rt[j,i-1] + (theta - alpha*rt[j,i-1])*dt[1] + sigma*dW
    }
  }
})
# ---------------------------------------------------------------



# Simulate short rate paths -- VECTORISED version (same Euler scheme, fast)
# ---------------------------------------------------------------
# Path dimension is vectorised: all Wiener increments are drawn at once,
# theta is precomputed across the timeline, and each time step updates
# all paths simultaneously.
set.seed(13)
rt_vec <- matrix(NA, paths, length(t))        # Matrix to store short rate paths
rt_vec[,1] <- r0                              # Set starting values to r0
rt_vec[1:10,1:15]                             # View structure of the matrix

# All Wiener increments drawn in one call.
# Total count = paths * (length(t) - 1), i.e. one increment per path, per time step.
# Reshaped into a matrix with one column per time step and one row per path.
dW_all  <- matrix(rnorm(paths * (length(t)-1), mean = 0, sd = sqrt(dt[1])), nrow = paths)

# Time-varying mean-reversion level:
# Deterministic (depends only on t and today's forward curve)!
# Hence, precomputable outside the loop.
# The stochastic part of the drift, -alpha*r_t, must stay inside
# the time loop because it depends on the realised path.
theta_v <- df0inst_dt(t[-length(t)]) + alpha * f0inst(t[-length(t)])

# Wrap the simulation in system.time() to record how long it takes.
# vec_time will be compared with loop_time below to quantify the speed-up.
vec_time <- system.time({

# The time dimension stays as a loop (Markov-sequential: each step depends
# on the previous one), but all paths are updated in a single operation.
for(i in 2:length(t)){
    # Levels of the short rate: whole column (i.e., all paths) updated in one step
    rt_vec[,i] <- rt_vec[,i-1] + (theta_v[i-1] - alpha*rt_vec[,i-1])*dt[1] + sigma*dW_all[,i-1]
  }
})
# ---------------------------------------------------------------



# Compare timings
# ---------------------------------------------------------------
# Interpreting the output: each row shows three fields.
#   user    = CPU time spent in your code
#   system  = CPU time spent in OS-level calls (memory management, etc.)
#   elapsed = total time
# Compare the two versions using "elapsed" -- the actual wait time.
# ---------------------------------------------------------------
print(loop_time)        # loop-based timing
print(vec_time)         # vectorised timing
# ---------------------------------------------------------------
# Note: rt and rt_vec are statistically equivalent (same distribution)
# but path-by-path different -- the two versions consume random numbers
# in different orders, so identical seeds give different arrangements.



# Plot short rate paths
# ---------------------------------------------------------------
# Paths are in decimal; multiply by 100 for percent-scale display.
# png(file="HullWhitePaths.png",width=2000, height=1000, res = 200)
plot(t, rt[1,]*100, type = "l", lwd = 1, xlim = c(0,20), ylim = c(0,5),
     xaxs = "i", yaxs = "i", col = "darkgrey",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,
     xlab = "Time (in Years)",
     ylab = "Short Rate in % p.a.")
for(i in 2:length(rt[,1])){
lines(t, rt[i,]*100, lwd = 1, col = "darkgrey")}
lines(t, rt[1,]*100, lwd = 1, col = "black")
lines(c(0,m), f0inst(c(0,m))*100, lwd = 3, col = "darkgreen")
legend("bottomright", c("Short Rate Paths", "Inst. Forward Curve"),
       lty = 1, lwd = c(1, 3), col = c("darkgrey", "darkgreen"),
       bty = "n", cex = 1.25, horiz = FALSE)
# dev.off()
# ---------------------------------------------------------------



# Term structures associated with the short rates
# -------------------------------------------------------------
# Note: term structures differ across times t and paths j
# All inputs and outputs are in decimal (consistent with y0 and f0inst).
y <- function(t,T,j){
        # Initial zero bond prices (from current term structure, all in decimal)
        P0T <- exp(-y0(T)*T)
        P0t <- exp(-y0(t)*t)
        # B coefficient
        B <- (1 - exp(-alpha*(T-t)))/alpha
        # A coefficient
        lnA <- log(P0T/P0t) + B*f0inst(t) - sigma^2/(4*alpha) * (1 - exp(-2*alpha*t)) * B^2
        # Translate t into index to select adequate short rate from matrix rt
        x <- ceiling(t/dt[1])
        # Term structure at future time t: P(t,T) = exp(-y_t(T)*(T-t))
        ytT <- -lnA/(T-t) + B/(T-t)*rt[j,x]
        ytT}
# ---------------------------------------------------------------



# Select individual term structures at future times s in path j
# ---------------------------------------------------------------
# Important: input T needs to be greater than s
# (time = 1 and T = 2 for the 1-year spot rate in one year)
time <- 5     # Select future time (in years from now)
path <- 13    # Select path
# Display term structure of s on path j (decimal)
ys <- y(t = time, T = time+maturity, j = path)
# Corresponding zero bond prices (ys already in decimal, no /100 needed)
Ps = exp(-ys*(maturity))

# Plot the term structure at the future time on path j
# y() returns decimal; multiply by 100 for percent-scale display.
plot(maturity, y(t = time,T = time+maturity,j = path)*100, type = "l", lwd = 5,
     xlim = c(0,30), ylim = c(0,5), xaxs = "i", yaxs = "i",
     main = "Term structure at time s on path j", col = "darkgreen",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,
     xlab = "Maturity (in years)",
     ylab = "Yield in % p.a.")
# ---------------------------------------------------------------



# Compute mean term structure across all paths for a future time t
# -----------------------------------------------------------------
# Create matrix to store a term structure for a given time across all paths
term.structures <- matrix(NA, paths, length(m))
# Loop stores the future term structures for a given time across all paths
for(i in 1:paths){
        term.structures[i,] <- y(time, time+m, i)}
# Compute column means of matrix to obtain vector with average term structure
mean_y <- colMeans(term.structures)
# ---------------------------------------------------------------



# Plot average future term structure
# ---------------------------------------------------------------
# Multiply by 100 for percent-scale display.
#png(file="Hull-WhiteYC.png",width=1000, height=1000, res = 200)
plot(m, mean_y*100, type = "l", lwd = 5,
     xlim = c(0,30), ylim = c(0,5), xaxs = "i", yaxs = "i",
     main = "t = 5", col = "darkgreen",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,
     xlab = "Maturity (in years)",
     ylab = "Yield in % p.a.")
for(i in 1:paths){lines(m, term.structures[i,]*100, lwd = 1, col = "darkgrey")}
lines(m, mean_y*100, lwd = 5, col = "darkgreen")
#dev.off()
# ---------------------------------------------------------------
