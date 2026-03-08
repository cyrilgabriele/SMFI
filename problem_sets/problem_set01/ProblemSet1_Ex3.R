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
Z <- rnorm(N_paths, mean=0, sd=1)
W_T <- sqrt(T) * Z
# usage here of the exact solution bc of hint but NOT risk-free
I_T_mu <- I_0 * exp((mu - 0.5 * sigma^2) * T +     # GBM exact formula
                  sigma * W_T)
# usage here of the exact solution bc of hint but AND risk-free
I_T_r <- I_0 * exp((r - 0.5 * sigma^2) * T +     # GBM exact formula
                      sigma * W_T)


hist(I_T_mu)
# ---
# b)
set.seed(13)
I_mat <- matrix(NA, nrow = N_paths, ncol = steps + 1)
I_mat[, 1] <- I_0

for (t in 1:steps) {
  Z_t <- rnorm(N_paths)
  W_t <- sqrt(dt) * Z_t
  I_mat[, t + 1] <- I_mat[, t] * exp((mu - 0.5 * sigma^2) * dt +
                                     sigma * W_t)
}

P_T <- apply(I_mat, 1, P_t, pi0 = pi_0, g = g, alpha = alpha)
hist(P_T)
# ---
# c) 

  
  
  