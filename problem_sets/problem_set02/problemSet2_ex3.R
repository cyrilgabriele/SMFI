library(quadprog)     # Load quadprog package

# Import & clean data (same as Ex. 1)
#-----------------------------------------------------------
# Save time series data from .csv file into a data frame
data <- read.csv2("Asset_Classes.csv", header = TRUE, stringsAsFactors = FALSE)

# Convert the first column of the data frame to date format
data$Date <- as.Date(data$Date, "%d.%m.%Y")

# Convert remaining columns to numerical format
for(i in 2:length(data[1,])){data[,i] <- as.numeric(data[,i])}

# Display head and check structure (str) of the data frame
head(data)
str(data)

# split now the data & use the hint i.e. by using the indices
# January 1993 - December 2002
train_data <- data[1:120, ]
# January 2003 - December 2012
test_data <- data[121:240, ]

#-----------------------------------------------------------
#-----------------------------------------------------------
# Markowitz model 

# expected returns, standard deviations, and the covariance matrix of the returns
# Number of risky assets: exclude Date and Rf
N <- ncol(train_data) - 2

# Create vector to store mean returns
MU <- rep(NA, N)

# Fill vector with monthly mean returns for each asset
for(i in 2:(ncol(train_data)-1)) {
  MU[i-1] <- mean(train_data[, i], na.rm = TRUE) * 12
}

# Covariance matrix of risky asset returns only
SIGMA <- cov(na.omit(train_data[, 2:(ncol(train_data)-1)])) * 12

# Variances and standard deviations
VAR <- diag(SIGMA)
SD <- sqrt(VAR)
#---------------------------
# preparation for the optimization
# Define number and type of constraints
n <- 2 + N   # Total number of constraints (2 basic + N short-sale constraints)
meq <- 2     # First "meq" constraints treated as equalities

# dvec set to zero (we specify the target function in terms of the variance only)
dvec <- rep(0,N)    

# Amat includes left-hand-side values of the constraints; see **slide 14**
Amat <- matrix(NA, N, n)    # Create matrix to store the constraints
Amat[,1] <- MU              # (i)   w1*mu1 + w2*mu2 = mup (find minimum variance for a given mean)
Amat[,2] <- rep(1, N)       # (ii)  SUM(w_i) = 1          (budget constraint) 
for(i in 1:N){              # (iii) w_i >= 0              (short-sale constraints)
  constraint <- rep(0,N)    # Create vector with N zeros
  constraint[i] <- 1        # Set 1 for the short-sale constraints (w*1), 0 for unconstrained frontier
  Amat[,i+2] <- constraint}

stepsize <- 0.0001                        # Stepsize for vector with fixed mean returns for optimization (first constraint)
mu_bar <- seq(min(MU),max(MU),stepsize)   # Fixed mean returns for optimization (first constraint) 
#mu_bar <- seq(0,1, stepsize)             # Use this for unconstrained frontier

#-----------------------------------------------------------
# Create objects to store efficient frontier
#-----------------------------------------------------------
w <- matrix(NA, N, length(mu_bar))        # Matrix for weights of efficient portfolios
# w[,1:10]                                # View first 10 columns of matrix w to see structure
mu <- rep(NA,length(mu_bar))              # Vector for means of efficient portfolios
sigma <- rep(NA,length(mu_bar))           # Vector for standard deviations of efficient portfolios
#-----------------------------------------------------------
#### Optimization for Efficient Frontier ####
# bvec = c(mu_bar[i], 1, rep(0,N))
# these contain the constrains from **page 14** of the lecture slides
#-----------------------------------------------------------
for(i in 1:length(mu_bar)) {
  w[,i] <- solve.QP(Dmat = SIGMA, 
                    dvec = dvec, 
                    Amat = Amat,
                    bvec = c(mu_bar[i], 1, rep(0,N)),
                    meq = meq)$solution
  mu[i] <- t(w[,i]) %*% MU
  sigma[i] <- sqrt(t(w[,i]) %*% SIGMA %*% w[,i])}


identical(mu,mu_bar)
#-----------------------------------------------------------^
# Find minimum variance portfolio
#-----------------------------------------------------------
# Look up the position of the smallest standard deviation 
# in the vector "sigma". This position index is stored as  
# scalar object "min_var"

min_var <- which.min(abs(min(sigma) - sigma)) 
#Alternative that finds closest value
#-----------------------------------------------------------
#########################################
# Plot the efficient frontier
#########################################
asset_names <- colnames(data)[2:(ncol(data)-1)]
eff <- mu >= mu[min_var]

x_all <- c(sigma[eff], SD)
y_all <- c(mu[eff], MU)

x_pad <- 0.05 * diff(range(x_all))
y_pad <- 0.05 * diff(range(y_all))

plot(sigma[eff], mu[eff], type = "l",
     lwd = 6,
     xlim = range(x_all) + c(-x_pad, x_pad),
     ylim = range(y_all) + c(-y_pad, y_pad),
     xaxs = "i", yaxs = "i", col = "darkgreen",
     main = "Efficient Frontier",
     xlab = "Return Standard Deviation", ylab = "Expected Return",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)

points(SD, MU, pch = 16, cex = 1.2, col = "black")
text(SD, MU, labels = asset_names, pos = 4, cex = 1)

points(sigma[min_var], mu[min_var],
       pch = 21, cex = 1.8, lwd = 2,
       col = "black", bg = "darkgrey")

#########################################
# Plot the 4 basic asset classes in the mu-sigma space
#########################################
asset_names <- colnames(train_data)[2:(ncol(train_data)-1)]

x_pad <- 0.05 * diff(range(SD))
y_pad <- 0.05 * diff(range(MU))

plot(SD, MU, type = "p",
     pch = 21, cex = 2, lwd = 2,
     col = "black", bg = "darkgrey",
     xlim = range(SD) + c(-x_pad, x_pad),
     ylim = range(MU) + c(-y_pad, y_pad),
     xaxs = "i", yaxs = "i",
     main = expression(paste("Basic Asset Class Portfolios in ", mu, "-", sigma, " Space")),
     xlab = "Return Standard Deviation", ylab = "Expected Return",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)

text(SD, MU, labels = asset_names, pos = 4, cex = 1)

#########################################
# 1/N Strategy
#########################################
# Monthly return vector of the equally weighted portfolio
ret_1N <- rowMeans(train_data[, 2:(ncol(train_data)-1)])

# Sample mean and sample standard deviation of historical monthly returns
mean_1N_monthly <- mean(ret_1N) 
mean_1N_annualized <- mean_1N_monthly * 12

sd_1N_monthly   <- sd(ret_1N)
sd_1N_annualized <- sqrt(12) * sd_1N_monthly 

idx_nearest_efficient_portfolio <- which.min(abs(sigma - sd_1N_annualized)) 
mu[idx_nearest_efficient_portfolio]
sigma[idx_nearest_efficient_portfolio]
w[, idx_nearest_efficient_portfolio]

#########################################
# Replot efficient frontier and add 1/N portfolio
#########################################

asset_names <- colnames(train_data)[2:(ncol(train_data)-1)]
eff <- mu >= mu[min_var]

x_all <- c(sigma[eff], SD, sd_1N_annualized)
y_all <- c(mu[eff], MU, mean_1N_annualized)

x_pad <- 0.05 * diff(range(x_all))
y_pad <- 0.05 * diff(range(y_all))

plot(sigma[eff], mu[eff], type = "l",
     lwd = 6,
     xlim = range(x_all) + c(-x_pad, x_pad),
     ylim = range(y_all) + c(-y_pad, y_pad),
     xaxs = "i", yaxs = "i", col = "darkgreen",
     main = "Efficient Frontier with 1/N Portfolio",
     xlab = "Return Standard Deviation", ylab = "Expected Return",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)

# Asset class portfolios
points(SD, MU, pch = 16, cex = 1.2, col = "black")
text(SD, MU, labels = asset_names, pos = 4, cex = 1)
# Minimum-variance portfolio
points(sigma[min_var], mu[min_var],
       pch = 21, cex = 1.8, lwd = 2,
       col = "black", bg = "darkgrey")
# --- --- ---
# 1/N portfolio
points(sd_1N_annualized, mean_1N_annualized,
       pch = 16, cex = 1.5, col = "blue")
text(sd_1N_annualized, mean_1N_annualized,
     labels = "1/N", pos = 4, cex = 1, col = "blue")
# --- --- --- 
# Efficient portfolio with same standard deviation as 1/N
points(sigma[idx_neareast_efficient_portfolio],
       mu[idx_neareast_efficient_portfolio],
       pch = 17, cex = 1.6, col = "red")
text(sigma[idx_neareast_efficient_portfolio],
     mu[idx_neareast_efficient_portfolio],
     labels = "Efficient portfolio", pos = 4, cex = 1, col = "red")
