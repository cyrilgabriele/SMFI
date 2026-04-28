
#############################################
# The Black model for European bond options
#############################################
# Need to load content of files:
# SMFI R Code Chapter 4 - i. Yield Curve.R
# SMFI R Code Chapter 4 - ii. Forward Curve.R
#
# Note: all rates are in DECIMAL form, consistent with y0(), f0(), and
# f0inst() from files (i) and (ii). No /100 conversions needed.

# Input parameters
# ---------------------------------------------------------------
T <- 1               # Maturity of the bond option
X <- 100             # Strike price
sigmaF <- 0.05       # Volatility of the forward returns
F0 <- 100            # Current forward bond price for maturity T
DF <- exp(-y0(T)*T)  # Discount factor (y0 already decimal)
# ---------------------------------------------------------------


# Implementation of the model for call options
# ---------------------------------------------------------------
C0 <- function(sigmaF, F0, X, T){
        d1 <- (log(F0/X)+(sigmaF^2/2)*T) / (sigmaF*sqrt(T))
        d2 <- d1 - sigmaF*sqrt(T)
        DF*(F0*pnorm(d1) - X*pnorm(d2))}
# ---------------------------------------------------------------


# Implementation of the model for put options
# ---------------------------------------------------------------
P0 <- function(sigmaF, F0, X, T){
        d1 <- (log(F0/X)+(sigmaF^2/2)*T) / (sigmaF*sqrt(T))
        d2 <- d1 - sigmaF*sqrt(T)
        DF*(X*pnorm(-d2) - F0*pnorm(-d1))}
# ---------------------------------------------------------------


# Compute option prices
# ---------------------
C0(sigmaF, F0, X, T)
P0(sigmaF, F0, X, T)
# ---------------------



######################################
# The Black model for caps and floors
######################################

# Input parameters
# ---------------------------------------------------------------
N <- 1000000               # Notional of the caplet/floorlet
tk <- 3                    # Fixing date of the underlying interest rate
tk1 <- 4                   # Maturity of the caplet/floorlet (t_{k+1})
delta_tk <- tk1 - tk       # Length of period between fixing date and maturity of the option
yX <- f0(tk,tk1)           # Strike rate of the caplet/floorlet (f0 = atm; decimal)
sigmaf <- 0.2              # Volatility of the forward rate
DF <- exp(-y0(tk1)*tk1)    # Discount factor (y0 already decimal)
# ---------------------------------------------------------------
# Note: the Black model does not work for negative forward rates!
# (short end of the empirical yield curve was negative in 2022)


# Implementation of the model for caplets
# ---------------------------------------------------------------
# Note: forward rate is determined endogenously (from current yield curve)
c0 <- function(sigmaf, yX){
        d1 <- (log(f0(tk,tk1)/yX)+(sigmaf^2/2)*tk) / (sigmaf*sqrt(tk))
        d2 <- d1 - sigmaf*sqrt(tk)
        N*delta_tk*DF*(f0(tk,tk1)*pnorm(d1) - yX*pnorm(d2))}
# ---------------------------------------------------------------


# Implementation of the model for floorlets
# ---------------------------------------------------------------
# Note: forward rate is determined endogenously (from current yield curve)
p0 <- function(sigmaf, yX){
        d1 <- (log(f0(tk,tk1)/yX)+(sigmaf^2/2)*tk) / (sigmaf*sqrt(tk))
        d2 <- d1 - sigmaf*sqrt(tk)
        N*delta_tk*DF*(yX*pnorm(-d2) - f0(tk,tk1)*pnorm(-d1))}
# ---------------------------------------------------------------

# Compute option prices
# ---------------------
c0(sigmaf, yX)
p0(sigmaf, yX)
# ---------------------



######################################
# The Black model for swaptions
######################################

# Input parameters
# ---------------------------------------------------------------
N <- 1000000                    # Notional of the swaption
T <- 2                          # Maturity of the swaption (start date of the swap)
Ts <- 6                         # Maturity of the underlying swap (Ts > T)
spread <- 0.02                  # Swap spread
f0S <- f0(T,Ts) + spread        # Forward swap rate (f0 already decimal)
sX <- 0.02                      # Strike rate of the swaption (decimal)
sigmaS <- 0.2                   # Volatility of the forward swap rate
# ---------------------------------------------------------------
# Note: the Black model does not work for negative forward rates!
# (short end of the empirical yield curve was negative in 2022)


# Implementation of the model for payer swaption
# ---------------------------------------------------------------
# Note: forward rate is determined endogenously (from current yield curve).
# To simplify, we assume that the forward swap rate is the forward rate
# plus the swap spread.
# The sum runs over swap payment dates T+1, T+2, ..., Ts (after the swaption
# matures and the swap begins). Loop variable k avoids shadowing tk from the
# caplet section.
Pswaption0 <- function(sigmaS, sX, T, Ts){
        PV <- rep(NA, Ts-T)
        for(k in (T+1):Ts){
                DF <- exp(-y0(k)*k)
                d1 <- (log(f0S/sX)+(sigmaS^2/2)*T) / (sigmaS*sqrt(T))
                d2 <- d1 - sigmaS*sqrt(T)
                PV[k-T] <- N*DF*(f0S*pnorm(d1) - sX*pnorm(d2))}
                sum(PV)}
# ---------------------------------------------------------------


# Implementation of the model for receiver swaption
# ---------------------------------------------------------------
# Note: forward rate is determined endogenously (from current yield curve)
Rswaption0 <- function(sigmaS, sX, T, Ts){
        PV <- rep(NA, Ts-T)
        for(k in (T+1):Ts){
                DF <- exp(-y0(k)*k)
                d1 <- (log(f0S/sX)+(sigmaS^2/2)*T) / (sigmaS*sqrt(T))
                d2 <- d1 - sigmaS*sqrt(T)
                PV[k-T] <- N*DF*(sX*pnorm(-d2) - f0S*pnorm(-d1))}
                sum(PV)}
# ---------------------------------------------------------------


# Compute option prices
# ---------------------------
Pswaption0(sigmaS, sX, T, Ts)
Rswaption0(sigmaS, sX, T, Ts)
# ---------------------------
