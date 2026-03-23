
#######################
# Two-Asset Model 
#######################

# Parameter values
#-----------------------------------------------------------
muA <- 0.06      # Expected return stock A => historical data
muB <- 0.045     # Expected return stock B => historical data
sigmaA <- 0.1    # Return standard deviation stock A => historical data
sigmaB <- 0.06   # Return standard deviation stock B => historical data
rhoAB <- 1       # Correlation of the returns => derived from the historical data
#-----------------------------------------------------------
wA <- seq(0,1,0.005)  # Sequence of different weights for A
wB <- (1-wA)          # Corresponding weights for B
#-----------------------------------------------------------
  
# Compute combinations of expected return and 
# return standard deviation for the two-asset portfolio
#-----------------------------------------------------------
# weighted avarage of the mean returns for stock A and stock B
muP <- wA * muA + wB * muB              
sigmaP <- sqrt(wA^2 * sigmaA^2 + wB^2 * sigmaB^2 
               + 2*wA*wB*rhoAB*sigmaA*sigmaB)
#-----------------------------------------------------------

# Plot results
#-----------------------------------------------------------
# png(file="TwoAssetCase.png",width=1000, height=1000, res = 200)
plot(sigmaP, muP, type = "l", 
     lwd = 6, xlim = c(0,0.15), ylim = c(0.03,0.07),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     main = paste("Correlation =", rhoAB), 
     xlab = "Return Standard Deviation", ylab = "Expected Return", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
points(sigmaA, muA, pch = 21, cex = 2, lwd = 2, col = "black", bg = "darkgrey")
points(sigmaB, muB, pch = 21, cex = 2, lwd = 2, col = "black", bg = "darkgrey")
# dev.off()
#-----------------------------------------------------------



########################
# n-Asset Model 
########################

# Parameter values
#-----------------------------------------------------------
mu <- 0.06           # Expected return 
sigma <- 0.08        # Return standard deviation 
rho <- 0.5           # Pairwise correlation 
n <- seq(0,1000,1)   # Vector with varying numbers of assets
#-----------------------------------------------------------

# Compute combinations of expected return and 
# return standard deviation for the portfolio
# of n homogeneous assets
#-----------------------------------------------------------
muP <- mu
sigmaP <- sqrt(1/n^2*(n*sigma^2 + (n^2-n)*rho*sigma^2))
systematic <- sigma*sqrt(rho)
#-----------------------------------------------------------

# Plot results
#-----------------------------------------------------------
# png(file="nAssetCase.png",width=1000, height=1000, res = 200)
plot(n, sigmaP, type = "l", 
     lwd = 6, xlim = c(1,20), ylim = c(0,0.12),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     main = paste("Correlation =", rho), 
     xlab = "Number of Assets n", ylab = "Return Standard Deviation", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
abline(h = systematic, lty = 2, lwd = 4, col = "black")
# dev.off()
#-----------------------------------------------------------
