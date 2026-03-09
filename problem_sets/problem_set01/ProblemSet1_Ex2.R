# ============================================================
# Problem Set 1 — Exercise 2
# ============================================================
set.seed(13)

# ---------- Inputs ----------
S0      <- 100
sigma   <- 0.20
r       <- 0.015
T       <- 2
K1      <- 105
K2      <- 115
H       <- 20000
N_paths <- 100000

# a) simulation of S_T under Q (risk-free assumption)
Z <- rnorm(N_paths, mean=0, sd=1)
# usage here of the exact solution bc of hint
S_T <- S0 * exp((r - 0.5 * sigma^2) * T +     # GBM exact formula
                sigma * sqrt(T) * Z)

# b) payoff at maturity
payoff_indicator <- as.integer(K1 <= S_T & S_T <= K2)
D_T <- H * payoff_indicator

hist(payoff_indicator,
     breaks = seq(-0.5, 1.5, by = 1),
     col = "lightblue", border = "white",
     xaxt = "n",
     main = "Double digital payoff counter",
     xlab = "Pays payoff?")
axis(1, at = c(0, 1), labels = c("No", "Yes"))

# c) discount the expected payoff to get the price of the option today at t=0
probability_payofff = mean(payoff_indicator)
expected_payoff = mean(D_T)
V0 <- exp(-r * T) * expected_payoff
# V0 = 2352.942
cat(sprintf("V0      : %.6f\n", V0))
