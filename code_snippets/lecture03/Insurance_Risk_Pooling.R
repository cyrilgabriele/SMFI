
#########################################
# Insurance risk pooling 
#########################################

#------------------------------------------------------
n <- seq(1,3000,1) # Number of policies
m <- 1000          # Mean of individual claim distribution  
s <- 2000          # Standard deviation of individual claim distribution
c <- 0.1*m         # Constant risk loading c (10% of EL)
#------------------------------------------------------
EL <- n*m                    # Expected total portfolio loss
SD <- sqrt(n)*s              # Standard deviation of portfolio loss
PI <- n*(m+c)                # Aggregate portfolio premium
RP <- 1 - pnorm(c/s*sqrt(n)) # Ruin probability
#------------------------------------------------------

# Create graphs for EL and premium volume
#------------------------------------------------------
# png(file="PoolingEL.png",width=1000, height=1000, res = 200)
plot(n, PI, type = "l", 
     lwd = 6, xlim = c(1,length(n)), ylim = c(0,10^6.5),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     main = "Type II Pooling: EL and Premium",
     xlab = "No. of Risks", ylab = "Portfolio EL and Premium", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
lines(n, EL, lwd = 6)
legend("topleft", c(
       expression(paste("E(",plain(L[P]),")")),
       expression(paste(plain(Pi[P])))), 
       cex = 1.5, lty = c(1,1), lwd = 6, col = c("black", "darkgreen"))
# dev.off()
#------------------------------------------------------

# Create graphs for ruin probability
#------------------------------------------------------
# png(file="PoolingRP.png",width=1000, height=1000, res = 200)
plot(n, RP, type = "l", 
     lwd = 6, xlim = c(1,length(n)), ylim = c(0,0.5),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     main = "Type II Pooling: RP", 
     xlab = "No. of Risks", ylab = "Ruin Probability", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
# dev.off()
#------------------------------------------------------

# Create graph portfolio standard deviation 
#------------------------------------------------------
# png(file="PoolingSD.png",width=1000, height=1000, res = 200)
plot(n, SD, type = "l", 
     lwd = 6, xlim = c(1,length(n)), ylim = c(0,1.25*10^5),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     main = "Type II Pooling: SD", 
     xlab = "No. of Risks", ylab = "Portfolio Standard Deviation", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
# dev.off()
#------------------------------------------------------



# Illustrate distribution of realized portfolio losses
#-----------------------------------------------------------------------
n <- 1000                        # Number of policies in the portfolio
EL <- n*m                        # Expected portfolio loss
SD <- sqrt(n)*s                  # Standard deviation of portfolio loss
ELmn <- EL/1000000               # Express portfolio EL in millions
SDmn <- SD/1000000               # Express portfolio SD in millions
iter <- 10000                    # Number of simulations for realized portfolio losses
RL <- rnorm(iter,n*m,sqrt(n)*s)  # Simulate realized portfolio losses
RLmn <- RL/1000000               # Express RL in millions
#-----------------------------------------------------------------------
RAL <- RL/n                      # Realized average loss
mean(RAL)
sd(RAL)
#-----------------------------------------------------------------------

# Create a histogram of the realized total portfolio losses
#-------------------------------------------------------------------------
# png(file="PoolingDist.png",width=1000, height=1000, res = 200)
hist(RLmn, xlim = c(0.5, 4), ylim = c(0,7), freq = FALSE,  
     breaks = 30, col = "grey", xlab = "Total Portfolio Loss in Millions", ylab = "Density", 
     xaxs = "i", yaxs = "i", main = "portfolio size n = 1000",  
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
# Overlay standard normal pdf
x <- seq(0,4,0.001)
y <- dnorm(x, mean = ELmn, sd = SDmn)
lines(x, y, lwd = 4, col = "darkgreen") 
legend("topright", "N(EL,SD)", lty=1, col = "darkgreen", lwd = 4, bty = TRUE)
# dev.off()
#-------------------------------------------------------------------------

# Create a histogram of the realized average portfolio losses (per policy)
#-------------------------------------------------------------------------
# png(file="AverageLosses.png",width=1000, height=1000, res = 200)
hist(RAL/1000, xlim = c(0.5, 1.5), ylim = c(0,12), freq = FALSE,  
     breaks = 30, col = "grey", xlab = "Average Loss in Thousands", ylab = "Density", 
     xaxs = "i", yaxs = "i", main = "portfolio size n = 1000",  
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
# Overlay standard normal pdf
x <- seq(0,2000,0.001)
y <- dnorm(x, mean = m/1000, sd = s/sqrt(n)/1000)
lines(x, y, lwd = 4, col = "darkgreen") 
legend("topright", "N(m,s/sqrt(n))", lty=1, col = "darkgreen", lwd = 4, bty = TRUE)
# dev.off()
#-------------------------------------------------------------------------

