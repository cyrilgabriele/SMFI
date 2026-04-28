
############################
# Load packages
############################
#install.packages("YieldCurve")
library(YieldCurve)


############################
# Prepare data
############################

# Import interest rate data 
#----------------------------------------------------------------
yc.data <- read.csv2("EuroAreaYieldCurves.csv", header = TRUE, stringsAsFactors=FALSE)
str(yc.data)   # Column names show the yield curves available in the file
head(yc.data)
yields <- as.numeric(yc.data[,4])         # Extract the chosen column
#----------------------------------------------------------------



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
#----------------------------------------------------------------


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
#----------------------------------------------------------------


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



###############################################
# Archetypical yield curve shapes 
###############################################

m <- seq(0.25,30,0.1) # Define a fine-grained maturity grid

# Implement the Nelson-Siegel-Svensson model for the spot yield curve
#----------------------------------------------------------------
# Note: y0NSS is for pedagogical shape experiments with hand-chosen parameters
# and returns yields in PERCENT (unlike y0 above, which returns decimal).
# The parameters in the archetype plots below are tuned to give percent values.
y0NSS <- function(T, b0, b1, b2, b3, lambda1, lambda2){
  F1 <- (1-exp(-lambda1*T))/(lambda1*T)
  F2 <- (1-exp(-lambda1*T))/(lambda1*T) - exp(-lambda1*T)
  F3 <- (1-exp(-lambda2*T))/(lambda2*T) - exp(-lambda2*T)
  b0 + b1*F1 + b2*F2 + b3*F3}
#----------------------------------------------------------------


# Plot upward-sloping yield curve
#----------------------------------------------------------------
#png(file="upward-slopingYC.png",width=1000, height=1000, res = 200)
plot(m, y0NSS(m, 4.00, -3.00, 25.00, -25.00, 1.00, 1.00), type = "l", lwd = 5, 
     xlim = c(0,30), ylim = c(0,6), xaxs = "i", yaxs = "i",
     main = "upward-sloping", col = "darkgreen",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
     xlab = "Maturity (in years)",
     ylab = "Yield in % p.a.")
#dev.off()
#----------------------------------------------------------------


# Plot inverted yield curve
#----------------------------------------------------------------
#png(file="invertedYC.png",width=1000, height=1000, res = 200)
plot(m, y0NSS(m, 2.00, 3.00, 25.00, -25.00, 1.00, 1.00), type = "l", lwd = 5, 
     xlim = c(0,30), ylim = c(0,6), xaxs = "i", yaxs = "i",
     main = "inverted", col = "darkgreen",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
     xlab = "Maturity (in years)",
     ylab = "Yield in % p.a.")
#dev.off()
#----------------------------------------------------------------


# Plot hump-shaped yield curve
#----------------------------------------------------------------
#png(file="hump-shapedYC.png",width=1000, height=1000, res = 200)
plot(m, y0NSS(m, 2.00, -2.00, 10.00, 0, 1.00, 1.00), type = "l", lwd = 5, 
     xlim = c(0,30), ylim = c(0,6), xaxs = "i", yaxs = "i",
     main = "hump-shaped", col = "darkgreen",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, 
     xlab = "Maturity (in years)",
     ylab = "Yield in % p.a.")
#dev.off()
#----------------------------------------------------------------



###############################################
# Experiment with the NSS model
###############################################
# Baseline:  (b0, b1, b2, b3, lambda1, lambda2) = (3, -1, 0, 0, 1, 0.3)

# (1) beta_0 -- level (whole curve shifts up/down)
#----------------------------------------------------------------
plot(m,  y0NSS(m, 2, -1, 0, 0, 1, 0.3), type="l", lwd=3, col="darkred",
     xlim=c(0,30), ylim=c(0,6), xaxs="i", yaxs="i",
     main=expression(paste("Effect of ", beta[0])),
     xlab="Maturity (in years)", ylab="Yield in % p.a.")
lines(m, y0NSS(m, 3, -1, 0, 0, 1, 0.3), lwd=3, col="darkorange")
lines(m, y0NSS(m, 4, -1, 0, 0, 1, 0.3), lwd=3, col="darkgreen")
#----------------------------------------------------------------


# (2) beta_1 -- slope (short end swings; long end pinned at beta_0)
#----------------------------------------------------------------
plot(m,  y0NSS(m, 3, -3, 0, 0, 1, 0.3), type="l", lwd=3, col="darkred",
     xlim=c(0,30), ylim=c(0,6), xaxs="i", yaxs="i",
     main=expression(paste("Effect of ", beta[1])),
     xlab="Maturity (in years)", ylab="Yield in % p.a.")
lines(m, y0NSS(m, 3,  0, 0, 0, 1, 0.3), lwd=3, col="darkorange")
lines(m, y0NSS(m, 3,  2, 0, 0, 1, 0.3), lwd=3, col="darkgreen")
#----------------------------------------------------------------


# (3) beta_2 -- hump amplitude (+ hump, - trough at intermediate maturities)
#----------------------------------------------------------------
plot(m,  y0NSS(m, 3, -1, -5, 0, 1, 0.3), type="l", lwd=3, col="darkred",
     xlim=c(0,30), ylim=c(0,6), xaxs="i", yaxs="i",
     main=expression(paste("Effect of ", beta[2])),
     xlab="Maturity (in years)", ylab="Yield in % p.a.")
lines(m, y0NSS(m, 3, -1,  5, 0, 1, 0.3), lwd=3, col="darkorange")
lines(m, y0NSS(m, 3, -1, 10, 0, 1, 0.3), lwd=3, col="darkgreen")
#----------------------------------------------------------------


# (4) lambda_1 -- hump location (peak near 1.79/lambda_1)
#----------------------------------------------------------------
plot(m,  y0NSS(m, 3, -1, 5, 0, 0.2, 0.3), type="l", lwd=3, col="darkred",
     xlim=c(0,30), ylim=c(0,6), xaxs="i", yaxs="i",
     main=expression(paste("Effect of ", lambda[1])),
     xlab="Maturity (in years)", ylab="Yield in % p.a.")
lines(m, y0NSS(m, 3, -1, 5, 0, 1.0, 0.3), lwd=3, col="darkorange")
lines(m, y0NSS(m, 3, -1, 5, 0, 2.0, 0.3), lwd=3, col="darkgreen")
#----------------------------------------------------------------


# (5) beta_3 -- Svensson extension (NS single hump vs NSS double hump)
#----------------------------------------------------------------
plot(m,  y0NSS(m, 3, -1, 3,  0, 2, 0.3), type="l", lwd=3, col="darkred",
     xlim=c(0,30), ylim=c(0,6), xaxs="i", yaxs="i",
     main="Nelson-Siegel vs Svensson",
     xlab="Maturity (in years)", ylab="Yield in % p.a.")
lines(m, y0NSS(m, 3, -1, 3, -4, 2, 0.3), lwd=3, col="darkgreen")
#----------------------------------------------------------------


# (6) lambda_2 -- second-hump location (only matters when beta_3 != 0)
#----------------------------------------------------------------
plot(m,  y0NSS(m, 3, -1, 3, -4, 2, 0.1), type="l", lwd=3, col="darkred",
     xlim=c(0,30), ylim=c(0,6), xaxs="i", yaxs="i",
     main=expression(paste("Effect of ", lambda[2])),
     xlab="Maturity (in years)", ylab="Yield in % p.a.")
lines(m, y0NSS(m, 3, -1, 3, -4, 2, 0.3), lwd=3, col="darkorange")
lines(m, y0NSS(m, 3, -1, 3, -4, 2, 0.8), lwd=3, col="darkgreen")
#----------------------------------------------------------------

