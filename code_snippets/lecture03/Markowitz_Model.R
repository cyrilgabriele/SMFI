
#########################################
# Prepare input data
#########################################

# Packages and working directory 
#-----------------------------------------------------------
library(quadprog)   # Load quadprog package
#getwd()            # Get working directory
#setwd("")          # Set working directory 
#-----------------------------------------------------------


# Import data
#-----------------------------------------------------------
# Save time series data from .csv file into a data frame
data <- read.csv2("./SMI_Stocks.csv", header = TRUE, sep =";", stringsAsFactors = FALSE, strip.white = TRUE)

# Convert the first column of the data frame to date format
data$Date <- as.Date(data$Date, "%d.%m.%Y")  

# Convert the second column of the data frame to numerical format
for(i in 2:length(data[1,])){data[,i] <- as.numeric(data[,i])}

# Display head and check structure (str) of the data frame
head(data)  
str(data)
#-----------------------------------------------------------


# Compute returns from prices
#-----------------------------------------------------------
# Clone data frame to store returns (instead of prices)
returns <- data                                

# Overwrite prices in first line with NAs
returns[1,2:length(returns[1,])] <- NA         

# Compute and store returns (columns 2 to 19)
for(j in 2:length(data[1,])){                   # Outer loop: columns
  for(i in 2:length(data[,1])){                 # Inner loop: rows
    returns[i,j] <- data[i,j]/data[i-1,j]-1  }} # Save returns to dataframe

# Check return data frame
head(returns)                                  
#-----------------------------------------------------------


# Estimate vector of mean returns
#-----------------------------------------------------------
# Number of assets (subtract date column)
N <- length(data[1,])-1                        

# Create vector to store mean returns for 20 stocks
MU <- rep(NA, N)                               

# Fill vector with mean returns
for(i in 2:length(returns[1,])){               
  MU[i-1] <- round(mean(returns[,i], na.rm = TRUE)*12,4)}
#-----------------------------------------------------------


# Estimate variance-covariance matrix
#-----------------------------------------------------------
# Na.omit for case-wise exclusion
SIGMA <- cov(na.omit(returns[,-1]))*12         

# Variances found on main diagonal
VAR <- diag(SIGMA)

# Compute return standard deviations (volatilities)
SD <- sqrt(VAR)
#-----------------------------------------------------------



##################
# Markowitz model 
##################

#### Prepare objects for optimization function ####
#### Objects are: dvec, Dmat, bvec, Amat
#-----------------------------------------------------------
# Use ?solve.QP to open help page and see notation of the target function 
?solve.QP # portfolio weights vector is "b"

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

# bvec includes right-hand-side values of the constraints (b0); see slide 14
# (must be indexed and run inside the loop)
#-----------------------------------------------------------


# Set range of portfolio means along which to minimize variance
#-----------------------------------------------------------
stepsize <- 0.0001                        # Stepsize for vector with fixed mean returns for optimization (first constraint)
mu_bar <- seq(min(MU),max(MU),stepsize)   # Fixed mean returns for optimization (first constraint) 
#mu_bar <- seq(0,1, stepsize)             # Use this for unconstrained frontier
#-----------------------------------------------------------

# Create objects to store efficient frontier
#-----------------------------------------------------------
w <- matrix(NA, N, length(mu_bar))        # Matrix for weights of efficient portfolios
w[,1:10]                                  # View first 10 columns of matrix w to see structure
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
#-----------------------------------------------------------

# Find minimum variance portfolio
#-----------------------------------------------------------
min_var <- match(min(sigma, na.rm = TRUE), sigma)
# Look up the position of the smallest standard deviation 
# in the vector "sigma". This position index is stored as  
# scalar object "min_var"

min_var <- which.min(abs(min(sigma) - sigma)) 
#Alternative that finds closest value
#-----------------------------------------------------------



#########################################
# Plot the efficient frontier
#########################################

#png(file="placeholder.png",width=1000, height=1000, res = 200)
#-----------------------------------------------------------
plot(sigma, mu, type = "l", 
     lwd = 6, xlim = c(0,0.3), ylim = c(0,0.3),
     xaxs = "i", yaxs = "i", col = "darkgreen", 
     main = "Efficient Frontier", 
     xlab = "Return Standard Deviation", ylab = "Expected Return", 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
points(SD, MU, pch = 16, cex = 1, lwd = 2, col = "black")
points(sigma[min_var], mu[min_var], pch = 19, cex = 1.5, lwd = 2, col = "black", bg = "darkgrey")
#dev.off()
#-----------------------------------------------------------
