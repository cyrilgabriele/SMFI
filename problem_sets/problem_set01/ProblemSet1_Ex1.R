# ============================================================
# Problem Set 1 — Exercise 1 (a)–(c): GBM simulation (Euler)
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
# (a) One path (set.seed(13)), Euler returns + price path
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

# ============================================================
# (b) mu_hat from the single path
# ============================================================
R_bar  <- mean(Rt)
mu_hat <- (1 / dt) * R_bar   # = 250 * R_bar since dt = 1/250

mu_muHat_abs_diff <- abs(mu - mu_hat)

# ============================================================
# (c) N = 10,000 independent return paths (set.seed(13) again!)
# ============================================================
set.seed(13)   # required by the problem set

N_paths <- 10000

Z_mat  <- matrix(rnorm(steps * N_paths, mean = 0, sd = 1), nrow = steps, ncol = N_paths)
dW_mat <- sqrt(dt) * Z_mat
R_mat  <- mu * dt + sigma * dW_mat

mu_hat_mc <- (1 / dt) * mean(R_mat)      # = 250 * mean(R_mat)

mu_muHatmc_abs_diff     <- abs(mu - mu_hat_mc)
muHat_muHatmc_abs_diff  <- abs(mu_hat - mu_hat_mc)