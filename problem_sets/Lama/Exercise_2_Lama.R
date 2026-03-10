####################################################
# Exercise 2: Pricing a Double Digital Option
####################################################

# Exercise 2 a)
S0 <- 100                         # Current stock price
#mu <- 0.08                       # Mean return rate 8%
sigma <- 0.2                      # Volatility of returns 20%
r <- 0.015                        # Continuously compounded risk-free rate
T <- 2                            # Price path over 2 year
K1 <- 105                         # Lower Strike
K2 <- 115                         # Upper Strike
H <- 20000                        # Lump Sum Payoff
paths <- 100000                   # Number of Simulations
set.seed(13)                      # Set random number


# Simulate stock prices at maturity
#-------------------------------------------
# Simulate Wiener process levels at maturity
WT <- sqrt(T)*rnorm(paths, mean = 0, sd = 1)             
# Simulate stock prices at maturity
ST <- S0*exp((r-sigma^2/2)*T + sigma*WT)



# Compute option payoffs
#-------------------------------------------
# The indicator function: 1 if K1 <= S_T <= K2, else 0
DT <- H * as.numeric(ST >= K1 & ST <= K2)


# Discount expected payoff to get the option price
#-------------------------------------------
V0 <- exp(-r * T) * mean(DT)


# Report result
#-------------------------------------------
cat("Double Digital Option Price (V0):", round(V0, 2), "\n")
cat("Average payoff E^Q[D_T]:", round(mean(DT), 2), "\n")
cat("Fraction of paths in [K1, K2]:", round(mean(ST >= K1 & ST <= K2), 4), "\n")




# Report: The simulated price of the double digital option is V0 = 2352.94
# (Value depends on seed). This reflects the discounted expected payoff under
# the risk-neutral measure Q.
