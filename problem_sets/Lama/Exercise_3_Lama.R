################################################################################
# Exercise 3: Guarantee Product
################################################################################

#-------------------------------------------------------------------------------
# Setting the parameters
#-------------------------------------------------------------------------------
pi0 <- 100000         # Investment
I0 <- 1000            # Current index value
mu <- 0.05            # Drift
sigma <- 0.13         # Volatility
r <- 0.01             # Risk-free interest rate
g <- 0.02             # Annual guaranteed interest rate
alpha <- 0.3          # Participation rate
T <- 20               # Maturity
dt <- 1               # Annual time steps
paths <- 100000       # Paths
set.seed(13)          # Use common random numbers


# Storage matrices
#-------------------------------------------------------------------------------
It <- matrix(NA, nrow = paths, ncol = T)    # Index Levels
Rt <- matrix(NA, nrow = paths, ncol = T)    # Continuous returns

#-------------------------------------------------------------------------------
# Part 3. a):
#-------------------------------------------------------------------------------

# Generate all continuous returns under measure Q
for (i in 1:paths) {
  # Wiener process increments
  dW <- sqrt(dt)*rnorm(T, mean = 0, sd = 1)
  
  # Calculate return transition
  Rt[i,] <- (r-sigma^2/2)*dt + sigma*dW
}


# Turn returns into paths of the Index
It[,1] <- I0 * exp(Rt[, 1])

for(i in 1:paths){
  for(t in 2:T){
    It[i, t] <- It[i, t-1]*exp(Rt[i, t])
  }
}


# Report a histogram of the index values at time T
I_T <- It[, T] 

hist(I_T, breaks = 18, col = "darkgrey", xaxs = "i", yaxs = "i", 
     main = "Index distribution at T=20", freq = TRUE,
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
     xlab = expression(paste(plain(I[T]))))

#-------------------------------------------------------------------------------
# Part 3. b)
#-------------------------------------------------------------------------------

Rp <- matrix(NA, nrow = paths, ncol = T)    # Product annual returns
PT <- rep(NA, paths)                        # Final product payoff at maturity

# Calculate Rp and PT for each path
for(i in 1:paths){
  # 1. Calculate return for Year 1
  idx_return_1 <- (It[i, 1] - I0) / I0
  Rp[i, 1] <- max(g, alpha * idx_return_1)
  
  # 2. Calculate returns for Years 2 to T
  for(t in 2:T){
    idx_return_t <- (It[i, t] - It[i, t-1]) / It[i, t-1]
    Rp[i, t] <- max(g, alpha * idx_return_t)
  }
  
  # 3. Calculate the final payoff PT for this specific path
  PT[i] <- pi0 * prod(1 + Rp[i, ])
}

# Report: Histogram of the product payoffs
hist(PT, breaks = 100, col = "darkgreen", xaxs = "i", yaxs = "i", 
     main = "Product Payoff Distribution at Maturity T=20", freq = TRUE,
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
     xlab = expression(paste(plain(P[T]))))

#-------------------------------------------------------------------------------
# Part 3. c)
#-------------------------------------------------------------------------------
# The expected payoff under the risk-neutral measure Q is the mean of PT
expected_payoff <- mean(PT)

# Discount the expected payoff to present value (t=0)
V0_product <- exp(-r * T) * expected_payoff

# Report the final price
cat("Price of the Guarantee Product (V0):", round(V0_product, 2), "\n")
cat("Initial Investment (pi0):", pi0, "\n")
cat("Difference:", round(V0_product - pi0, 2), "\n")


# Report: The price of the guarantee product V0 145'184.3 exceeds the initial
# investment of pi0 = 100'000. This is because the product contains an floor 
# option: each year, the investor earns at least g = 2%, regardless of
# index performance. This downside protection has value. It is essentially a
# series of put options. Under Q, this optionality is priced in, making the fair
# value of the product higher than the initial investment.