########################
# GET LFSR PARAMS for SB
# C PETERSON - 2021
########################
# Example of getting low-fecundity stock-recruit params for Sandbar shark
#

# length of simulation
years <- 1000


N0 <- 1000 # arbitrarily chosen


################# define natural mortality
M_const <- c(
  0.1604, 0.1604, 0.1604, 0.1604, 0.1604, 0.1578, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168,
  0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 
  0.1168, 0.1168, 0.1168, 0.1168, 0.1168
)


#################### Run SImulation
# Create matrices to store results
Nay_M <- matrix(nrow = A, ncol = years)
Nay_F <- matrix(nrow = A, ncol = years)
M0 <- vector(length = years)
Npups <- vector(length = years)




####### LFSR - FOR STOCK SYNTHESIS ####################


# Inputs:
M <- M_const

# Define 1st recruitment class
Nay_M[1, 1] <- N0 / 2
Nay_F[1, 1] <- N0 / 2

# populate first row
for (i in 1:(A - 1)) {
  Nay_M[i + 1, 1] <- Nay_M[i, 1] * exp(-M[i])
  Nay_F[i + 1, 1] <- Nay_F[i, 1] * exp(-M[i])
}

Nay_M[A, 1] <- (Nay_M[A - 1, 1] * exp(-M[A - 1])) + (Nay_M[A, 1] * exp(-M[A]))
Nay_F[A, 1] <- (Nay_F[A - 1, 1] * exp(-M[A - 1])) + (Nay_F[A, 1] * exp(-M[A]))

####### FOR LFSR
#
# TO ESTIMATE Survival-based Recruitmetn parameters from (Taylor et al. 2012)
# 1. generate initial N at age vector (Na0), given M at age vector and arbitrary N0
# 2. calculate Npups0, given maturity & fecundity at age vectors and Na0 as estimated above
# 3. Z0 = -ln(N0/Npups0)
# 4. set Zmin as fixed M for age 0
# 5. zfrac = 1- (zmin/z0)
# 6. beta =  ( ln( 1 - (ln(h/0.2) / (z0 * zfrac) ) ) ) / (ln(0.2))


Npups0 <- sum((Nay_F[, 1] / 2.5) * mat_F * fec)
Npups0 # 2961.558
# Npups0 = sum(Nay_F[,1]*mat_F*fec)
z0 <- -log(N0 / Npups0)
z0
# z0 = 1.085716
(s0 <- exp(-z0)) # 0.3424836
zmin <- 0.1604
(smax <- exp(-zmin)) # 0.851803
(zfrac <- (log(smax) - log(s0)) / (-log(s0))) # check. zfrac = 0.8522633
(zfrac <- 1 - (zmin / z0))
h <- 0.3
beta <- (log(1 - (log(h / 0.2) / (z0 * zfrac)))) / log(0.2)
beta # 0.3582577

(S_frac <- (zmin - z0) / (0 - z0)) # 0.8522633

Z_max <- z0 + S_frac * (0.0 - z0) # 0.1604




Z0 <- 1.085716
Zmin <- 0.1604
zfrac <- 0.8522633
beta <- 0.3582577
Npups_eq <- 2961.558 # Npups0
