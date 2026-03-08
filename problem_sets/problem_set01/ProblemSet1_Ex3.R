pi_0     <- 100000
I_0      <- 1000
mu      <- 0.05
sigma   <- 0.13
r       <- 0.01
g       <- 0.02
alpha   <- 0.30
T       <- 20
dt      <- 1
steps   <- 20
N_paths <- 100000
  
R_p_t <- function(I_t, I_t_minus_1, g, alpha) {
  pmax(g, alpha * ((I_t - I_t_minus_1) / I_t_minus_1))
}

P_t <- function(I_path, pi0, g, alpha) {
  R_p <- R_p_t(
    I_t         = I_path[-1],
    I_t_minus_1 = I_path[-length(I_path)],
    g           = g,
    alpha       = alpha
  )
  
  pi0 * prod(1 + R_p)
}

# a) simulation of I_T
set.seed(13)
I_mat_mu <- matrix(NA, nrow = N_paths, ncol = steps + 1)
I_mat_mu[, 1] <- I_0

for (t in 1:steps) {
  Z_t <- rnorm(N_paths)
  W_t <- sqrt(dt) * Z_t
  I_mat_mu[, t + 1] <- I_mat_mu[, t] * exp((mu - 0.5 * sigma^2) * dt +
                                             sigma * W_t)
}

I_T_mu <- I_mat_mu[, steps + 1]

hist(I_T_mu,
     main = expression(paste("Histogram of simulated ", tilde(I)[T])),
     xlab = expression(tilde(I)[T]),
     col = "lightblue",
     border = "white")
# ---
# b)
set.seed(13)
I_mat_mu <- matrix(NA, nrow = N_paths, ncol = steps + 1)
I_mat_mu[, 1] <- I_0

for (t in 1:steps) {
  Z_t <- rnorm(N_paths)
  W_t <- sqrt(dt) * Z_t
  I_mat_mu[, t + 1] <- I_mat_mu[, t] * exp((mu - 0.5 * sigma^2) * dt +
                                     sigma * W_t)
}

P_T_mu <- apply(I_mat_mu, 1, P_t, pi0 = pi_0, g = g, alpha = alpha)
hist(P_T_mu)
# ---
# c) 
set.seed(13)
I_mat_r <- matrix(NA, nrow = N_paths, ncol = steps + 1)
I_mat_r[, 1] <- I_0

for (t in 1:steps) {
  Z_t <- rnorm(N_paths)
  W_t <- sqrt(dt) * Z_t
  I_mat_r[, t + 1] <- I_mat_r[, t] * exp((r - 0.5 * sigma^2) * dt +
                                       sigma * W_t)
}

P_T_r <- apply(I_mat_r, 1, P_t, pi0 = pi_0, g = g, alpha = alpha)
hist(P_T_r)
expected_payoff <- mean(P_T_r)        # expected payoff at maturity under Q
V0 <- exp(-r * T) * expected_payoff   # product price at time 0

# economic explanation: 
# The product price differs from the initial investment because the contract
# promises a payoff stream that is richer than risk-free compounding:
# it guarantees 2% annually while pricing discounts at only 1%, and it
# additionally offers upside participation in the index. Therefore, its
# discounted expected payoff exceeds the initial investment.
