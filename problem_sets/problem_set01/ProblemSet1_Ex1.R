# ============================================================
# Problem Set 1 — Exercise 1
# ============================================================
# ---------- Inputs ----------
S0    <- 100
mu    <- 0.08
sigma <- 0.20
T     <- 1
s     <- 250

dt    <- 1 / s
steps <- T * s
tgrid <- seq(0, T, by = dt)   # length = steps + 1

# ============================================================
# (a) One path, Euler returns + price path
# ============================================================
set.seed(13)

Z  <- rnorm(steps, mean = 0, sd = 1)     # Z_t ~ N(0,1)
dW <- sqrt(dt) * Z                       # dW_t ~ N(0, dt)
Wt <- c(0, cumsum(dW))                   # Brownian path levels W_t

Rt <- mu * dt + sigma * dW               # Euler discrete returns
St <- c(S0, S0 * cumprod(1 + Rt))        # Price path via discrete compounding

# Plot: returns (top) and prices (bottom)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

plot(tgrid[-1], Rt, type = "l", col = "steelblue", lwd = 2,
     xlab = "Time t", ylab = "Return R_t", main = "Euler-simulated discrete returns")
abline(h = mean(Rt), col = "black", lwd = 2)

plot(tgrid, St, type = "l", col = "darkgreen", lwd = 2,
     xlab = "Time t", ylab = "Price S_t", main = "Corresponding price path")
abline(h = mean(St), col = "black", lwd = 2)

par(mfrow = c(1, 1))   # reset so subsequent plots are standalone

# ============================================================
# (b) mu_hat from the single path
# ============================================================
R_bar  <- mean(Rt)
mu_hat <- (1 / dt) * R_bar   # = 250 * R_bar since dt = 1/250

mu_muHat_abs_diff <- abs(mu - mu_hat)

# ============================================================
# (c) N = 10,000 independent return paths 
# ============================================================
set.seed(13)   # required by the problem set

N_paths <- 10000

Z_mat  <- matrix(rnorm(steps * N_paths, mean = 0, sd = 1), nrow = steps, ncol = N_paths)
dW_mat <- sqrt(dt) * Z_mat
R_mat  <- mu * dt + sigma * dW_mat

mu_hat_mc <- (1 / dt) * mean(R_mat)      # = 250 * mean(R_mat)

mu_muHatmc_abs_diff     <- abs(mu - mu_hat_mc)
muHat_muHatmc_abs_diff  <- abs(mu_hat - mu_hat_mc)

# ============================================================
# (d) Exact solution, N = 10,000 terminal prices at T = 3
# ============================================================
set.seed(13)

S0      <- 100
mu      <- 0.08
sigma   <- 0.20
T       <- 3
N_paths <- 100000

W_T <- rnorm(N_paths, mean = 0, sd = 1)  
# formula *1) (see below)
S_T <- S0 * exp((mu - 0.5 * sigma^2) * T +     # GBM approximation formula
                sigma * sqrt(T) * W_T)

hist(S_T, breaks = 60, col = "lightblue", border = "white",
     freq = FALSE,
     main = "GBM terminal prices S_T (T = 3)",
     xlab = expression(S[T]))

mu_log    <- log(S0) + (mu - 0.5 * sigma^2) * T
sigma_log <- sigma * sqrt(T)
curve(dlnorm(x, meanlog = mu_log, sdlog = sigma_log), add = TRUE,
      col = "darkred", lwd = 2)

S_T_sample_mean <- mean(S_T)
S_T_theoretical <- S0 * exp(mu * T)
S_T_mean_diff      <- S_T_sample_mean - S_T_theoretical
S_T_abs_mean_diff  <- abs(S_T_mean_diff)
# Comment: 
# Shape: the histogram indicates a log normal distribution.
# This since W_T is sampled from a normal dist. and then taken exponentially in
# the formula *1). This results in a log-normal dist. See also the PDF
# as a black curve with the parameters implied by the model.
# Obviously we expect a (small) difference between the analytical closed-form 
# solution and the Monte Carlo approixmation approach. The abs diff of ~0.07 
# is in my opinion small w.r.t to the total value of ~127. 
