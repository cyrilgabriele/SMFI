
############################
# Forward rates
############################
# Need to load content of files:
# SMFI R Code Chapter 4 - i. Yield Curve.R


# Important Note:
#----------------------------------------------------------------
# The forward rates computed here (both the finite-period f0(T1,T2)
# and the instantaneous f0inst(T)) are continuously compounded, consistent
# with the continuously-compounded spot rate y0(T) from the NSS model.
# All rates are returned in DECIMAL form (e.g. 0.025 for 2.5%).
#----------------------------------------------------------------


# Function for arbitrage-free forward rates (T2 > T1)
#----------------------------------------------------------------
# Inherits the decimal convention from y0(), so f0 also returns decimal.
f0 <- function(T1,T2){(y0(T2)*T2-y0(T1)*T1)/(T2-T1)}
# Forward rate examples
f0(0.5,1)       # Effective date: 6 months from now; Termination date: 1 year from now
f0(1,2)         # Effective date: 1 year from now; Termination date: 2 years from now
f0(10,12)       # Effective date: 10 years from now; Termination date: 12 years from now
#----------------------------------------------------------------



##########################################################
# Nelson-Siegel-Svensson instantaneous forward curve model
##########################################################

# Implement Nelson-Siegel-Svensson model for the instantaneous forward curve
#----------------------------------------------------------------
# Returns DECIMAL form (divide by 100 at the end), consistent with y0().
f0inst <- function(T){
        F1 <- exp(-lambda1*T)
        F2 <- exp(-lambda1*T)*lambda1*T
        F3 <- exp(-lambda2*T)*lambda2*T
        (b0 + b1*F1 + b2*F2 + b3*F3)/100}
#----------------------------------------------------------------


# Plot instantaneous forward curve
#----------------------------------------------------------------
# use "m" instead of "maturity" to obtain a smooth curve
# f0inst() returns decimal, so multiply by 100 for percent-scale display.
#png(file="EuroForwardCurve.png",width=2000, height=1000, res = 200)
plot(m, f0inst(m)*100, type = "l", cex = 1.25,
     xaxs = "i", yaxs = "i", lwd = 5, col = "darkgreen",
     main = "Euro Area Inst. Forward Curve",
     xlim = c(0,30), ylim = c(0,5),
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,
     xlab = "Maturity (in years)",
     ylab = "Yield in % p.a.")
#abline(h = 0, lwd = 3, lty = 2)
#dev.off()
#----------------------------------------------------------------
