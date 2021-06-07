### SB DFA SIMULATION FUNCTIONS ###
### C PETERSON - 2021 ###


### Load packages ###
library(scales)
library(MARSS)
library(vioplot)

rbvn <- function(n, m1, s1, m2, s2, rho) {
  X1 <- rnorm(n, m1, s1)
  X2 <- rnorm(n, m2 + (s2 / s1) * rho *
    (X1 - m1), sqrt((1 - rho^2) * s2^2))
  data.frame(Linf = X1, K = X2)
}



# Missing data function #
DFA_Sim <- function(Tot_F, CV1, CV2, CV3, CV4, CV5, CV6, CV7,
                    q_survey1, q_survey2, q_survey3, q_survey4, q_survey5, q_survey6, q_survey7,
                    Iay_eq = Iay_eq,
                    Niters = 100, c = 5, plot = TRUE) {

  # Inputs:
  M <- M_const
  # GET F
  Prop_F <- c(0.2, 0.044, 0.754, 0.002)
  F_const1 <- Tot_F * Prop_F[1]
  F_const2 <- Tot_F * Prop_F[2]
  F_const3 <- Tot_F * Prop_F[3]
  F_const4 <- Tot_F * Prop_F[4]
  FM <- rbind(F_const1, F_const2, F_const3, F_const4)
  My_LFreqs <- function(data, multC = 100, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90), Linf = Linf_F, K = K_F, t0 = t0_F) {
    data_N <- cbind(ages, round(data * multC))
    LF <- data.frame(Lbins)
    for (l in Lyrs) {
      df <- as.data.frame(data_N[, c(1, l + 1)])
      colnames(df) <- c("ages", "count")
      df.expanded <- df[rep(row.names(df), df$count), 1:2]
      vonBert <- rbvn(nrow(df.expanded), m1 = Linf, s1 = Linf * 0.1, m2 = K, s2 = K * 0.1, rho = -0.9)
      df.expanded$L <- vonBert$Linf * (1 - exp(-(vonBert$K) * (jitter(df.expanded$ages, amount = 0.5) - t0)))
      hst <- hist(df.expanded$L, breaks = c(seq(30, 262, by = 2)), plot = F)
      LF <- cbind(LF, hst$counts)
    }
    colnames(LF) <- c("LBins", as.character(Lyrs))
    return(LF)
  }


  # Iay
  par(mfrow = c(1, 1))
  plot(1:years, apply(Iay_eq, 2, sum), type = "l", ylim = c(0, max(apply(Iay_eq, 2, sum), na.rm = T) * 3))
  abline(v = 51, lwd = 2, col = "red")
  abline(v = 90, lwd = 2, col = "red")
  # plot(1:years, apply(Iay, 2, sum), type='l')
  # Iy = apply(Iay, 2, sum)


  ####### ADD STOCHASTICITY ###########
  Iay1_F <- matrix(nrow = A, ncol = years)
  Iay2_F <- matrix(nrow = A, ncol = years)
  Iay3_F <- matrix(nrow = A, ncol = years)
  Iay4_F <- matrix(nrow = A, ncol = years)
  Iay5_F <- matrix(nrow = A, ncol = years)
  Iay6_F <- matrix(nrow = A, ncol = years)
  Iay7_F <- matrix(nrow = A, ncol = years)
  Iay1_M <- matrix(nrow = A, ncol = years)
  Iay2_M <- matrix(nrow = A, ncol = years)
  Iay3_M <- matrix(nrow = A, ncol = years)
  Iay4_M <- matrix(nrow = A, ncol = years)
  Iay5_M <- matrix(nrow = A, ncol = years)
  Iay6_M <- matrix(nrow = A, ncol = years)
  Iay7_M <- matrix(nrow = A, ncol = years)
  Iay1_CV <- matrix(nrow = A, ncol = years)
  Iay2_CV <- matrix(nrow = A, ncol = years)
  Iay3_CV <- matrix(nrow = A, ncol = years)
  Iay4_CV <- matrix(nrow = A, ncol = years)
  Iay5_CV <- matrix(nrow = A, ncol = years)
  Iay6_CV <- matrix(nrow = A, ncol = years)
  Iay7_CV <- matrix(nrow = A, ncol = years)
  Iy1 <- matrix(nrow = Niters, ncol = years)
  Iy2 <- matrix(nrow = Niters, ncol = years)
  Iy3 <- matrix(nrow = Niters, ncol = years)
  Iy4 <- matrix(nrow = Niters, ncol = years)
  Iy5 <- matrix(nrow = Niters, ncol = years)
  Iy6 <- matrix(nrow = Niters, ncol = years)
  Iy7 <- matrix(nrow = Niters, ncol = years)
  Iy1_CV <- matrix(nrow = Niters, ncol = years)
  Iy2_CV <- matrix(nrow = Niters, ncol = years)
  Iy3_CV <- matrix(nrow = Niters, ncol = years)
  Iy4_CV <- matrix(nrow = Niters, ncol = years)
  Iy5_CV <- matrix(nrow = Niters, ncol = years)
  Iy6_CV <- matrix(nrow = Niters, ncol = years)
  Iy7_CV <- matrix(nrow = Niters, ncol = years)
  Cy_1 <- matrix(nrow = Niters, ncol = years)
  Cy_2 <- matrix(nrow = Niters, ncol = years)
  Cy_3 <- matrix(nrow = Niters, ncol = years)
  Cy_4 <- matrix(nrow = Niters, ncol = years)
  Ny_F <- matrix(nrow = Niters, ncol = years)
  Ny_M <- matrix(nrow = Niters, ncol = years)
  Cy <- matrix(nrow = Niters, ncol = years)
  Ny <- matrix(nrow = Niters, ncol = years)


  M0 <- matrix(nrow = Niters, ncol = years)
  Npups <- matrix(nrow = Niters, ncol = years)
  FSS <- matrix(nrow = Niters, ncol = years) # Female spawning stock

  ResultsList <- list()


  Z0 <- 1.085716
  Zmin <- 0.1604
  zfrac <- 0.8522633
  beta <- 0.3582577
  Npups_eq <- 2961.558 # Npups0



  set.seed(430)

  # F_const = rep(0, yrs)

  for (k in 1:Niters) {
    # paste("Iteration",k,sep=" ")
    # Create matrices to store results
    Nay_F <- matrix(nrow = A, ncol = years)
    Nay_M <- matrix(nrow = A, ncol = years)
    Cay1_F <- matrix(nrow = A, ncol = years)
    Cay2_F <- matrix(nrow = A, ncol = years)
    Cay3_F <- matrix(nrow = A, ncol = years)
    Cay4_F <- matrix(nrow = A, ncol = years)
    Cay1_M <- matrix(nrow = A, ncol = years)
    Cay2_M <- matrix(nrow = A, ncol = years)
    Cay3_M <- matrix(nrow = A, ncol = years)
    Cay4_M <- matrix(nrow = A, ncol = years)



    # Set up equilibrium conditions
    Nay_F[, 1] <- Equilibrium_Nay_F
    Nay_M[, 1] <- Equilibrium_Nay_M

    for (i in 1:(years - 1)) {
      # i=i+1

      matjit <- rnorm(length(mat_F), mat_F, mat_F * 0.01)
      matjit <- ifelse(matjit > 1, 1, matjit)
      FSS[k, i] <- sum(Nay_F[, i] * matjit) # Fecund Stock size
      Npups[k, i] <- sum((Nay_F[, i] / 2.5) * matjit * rnorm(length(fec), fec, fec * 0.1)) #
      M0[k, i] <- rnorm(1, (((1 - (Npups[k, i] / Npups_eq)^beta) * (Zmin - Z0)) + Z0), 0.1)
      # M0[i] = ( (1- (Npups[i]/Npups_eq)^beta ) * (Zmin - Z0) ) + Z0


      Nay_F[1, i + 1] <- 0.5 * Npups[k, i] * exp(-M0[k, i]) # pups produced that survive to age 1; Recruits
      Nay_M[1, i + 1] <- 0.5 * Npups[k, i] * exp(-M0[k, i]) # pups produced that survive to age 1; Recruits


      Zay_F <- vector(length = length(M))
      Zay_M <- vector(length = length(M))
      for (j in 1:A) {
        Zay_F[j] <- M[j] + (Sel[j, 2] * FM[1, i]) + (Sel[j, 4] * FM[2, i]) + (Sel[j, 6] * FM[3, i]) + (Sel[j, 8] * FM[4, i])
        Zay_M[j] <- M[j] + (Sel[j, 3] * FM[1, i]) + (Sel[j, 5] * FM[2, i]) + (Sel[j, 7] * FM[3, i]) + (Sel[j, 8] * FM[4, i])
      } # END J LOOP

      for (j in 1:(A - 2)) {
        Mjit <- rnorm(length(M), M, M * 0.1)
        Nay_F[j + 1, i + 1] <- Nay_F[j, i] * exp(-Zay_F[j])
        Nay_M[j + 1, i + 1] <- Nay_M[j, i] * exp(-Zay_M[j])


        Cay1_F[j, i] <- Nay_F[j, i] * ((Sel[j, 2] * FM[1, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
        Cay2_F[j, i] <- Nay_F[j, i] * ((Sel[j, 4] * FM[2, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
        Cay3_F[j, i] <- Nay_F[j, i] * ((Sel[j, 6] * FM[3, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
        Cay4_F[j, i] <- Nay_F[j, i] * ((Sel[j, 8] * FM[4, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))

        Cay1_M[j, i] <- Nay_M[j, i] * ((Sel[j, 3] * FM[1, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
        Cay2_M[j, i] <- Nay_M[j, i] * ((Sel[j, 5] * FM[2, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
        Cay3_M[j, i] <- Nay_M[j, i] * ((Sel[j, 7] * FM[3, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
        Cay4_M[j, i] <- Nay_M[j, i] * ((Sel[j, 8] * FM[4, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
      } # END J LOOP

      Nay_F[A, i + 1] <- Nay_F[A - 1, i] * exp(-Zay_F[A - 1]) #+ Nay_F[A,i]*exp(-Zay[A])
      Nay_M[A, i + 1] <- Nay_M[A - 1, i] * exp(-Zay_M[A - 1]) #+ Nay_M[A,i]*exp(-Zay[A])

      Cay1_F[A - 1, i] <- Nay_F[A - 1, i] * ((Sel[A - 1, 2] * FM[1, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
      Cay2_F[A - 1, i] <- Nay_F[A - 1, i] * ((Sel[A - 1, 4] * FM[2, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
      Cay3_F[A - 1, i] <- Nay_F[A - 1, i] * ((Sel[A - 1, 6] * FM[3, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
      Cay4_F[A - 1, i] <- Nay_F[A - 1, i] * ((Sel[A - 1, 8] * FM[4, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))

      Cay1_M[A - 1, i] <- Nay_M[A - 1, i] * ((Sel[A - 1, 3] * FM[1, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
      Cay2_M[A - 1, i] <- Nay_M[A - 1, i] * ((Sel[A - 1, 5] * FM[2, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
      Cay3_M[A - 1, i] <- Nay_M[A - 1, i] * ((Sel[A - 1, 7] * FM[3, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
      Cay4_M[A - 1, i] <- Nay_M[A - 1, i] * ((Sel[A - 1, 8] * FM[4, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))



      Cay1_F[A, i] <- Nay_F[A, i] * ((Sel[A, 2] * FM[1, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
      Cay2_F[A, i] <- Nay_F[A, i] * ((Sel[A, 4] * FM[2, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
      Cay3_F[A, i] <- Nay_F[A, i] * ((Sel[A, 6] * FM[3, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
      Cay4_F[A, i] <- Nay_F[A, i] * ((Sel[A, 8] * FM[4, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))

      Cay1_M[A, i] <- Nay_M[A, i] * ((Sel[A, 3] * FM[1, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
      Cay2_M[A, i] <- Nay_M[A, i] * ((Sel[A, 5] * FM[2, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
      Cay3_M[A, i] <- Nay_M[A, i] * ((Sel[A, 7] * FM[3, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
      Cay4_M[A, i] <- Nay_M[A, i] * ((Sel[A, 8] * FM[4, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))


      for (h in 1:A) {
        # h=h+1
        Iay1_CV[h, i] <- runif(1, min = CV1 - 0.1, max = CV1 + 0.1)
        Iay1_F[h, i] <- max(q_survey1[i] * Sel[h, 17] * Nay_F[h, i] *
          exp(rnorm(1, 0, (q_survey1[i] * Sel[h, 17] * Nay_F[h, i]) * Iay1_CV[h, i]) -
            (((q_survey1[i] * Sel[h, 17] * Nay_F[h, i] * Iay1_CV[h, i])^2) / 2)), 0) # Uncertainty based off CV
        Iay1_M[h, i] <- max(q_survey1[i] * Sel[h, 18] * Nay_M[h, i] *
          exp(rnorm(1, 0, (q_survey1[i] * Sel[h, 18] * Nay_M[h, i]) * Iay1_CV[h, i]) -
            (((q_survey1[i] * Sel[h, 18] * Nay_M[h, i] * Iay1_CV[h, i])^2) / 2)), 0) # Uncertainty based off CV

        Iay2_CV[h, i] <- runif(1, min = CV2 - 0.1, max = CV2 + 0.1)
        Iay2_F[h, i] <- max(q_survey2[i] * Sel[h, 9] * Nay_F[h, i] *
          exp(rnorm(1, 0, (q_survey2[i] * Sel[h, 9] * Nay_F[h, i]) * Iay2_CV[h, i]) -
            (((q_survey2[i] * Sel[h, 9] * Nay_F[h, i] * Iay2_CV[h, i])^2) / 2)), 0) # Uncertainty based off CV
        Iay2_M[h, i] <- max(q_survey2[i] * Sel[h, 10] * Nay_M[h, i] *
          exp(rnorm(1, 0, (q_survey2[i] * Sel[h, 10] * Nay_M[h, i]) * Iay2_CV[h, i]) -
            (((q_survey2[i] * Sel[h, 10] * Nay_M[h, i] * Iay2_CV[h, i])^2) / 2)), 0) # Uncertainty based off CV

        Iay3_CV[h, i] <- runif(1, min = CV3 - 0.1, max = CV3 + 0.1)
        Iay3_F[h, i] <- max(q_survey3[i] * Sel[h, 28] * Nay_F[h, i] *
          exp(rnorm(1, 0, (q_survey3[i] * Sel[h, 28] * Nay_F[h, i]) * Iay3_CV[h, i]) -
            (((q_survey3[i] * Sel[h, 28] * Nay_F[h, i] * Iay3_CV[h, i])^2) / 2)), 0) # Uncertainty based off CV
        Iay3_M[h, i] <- max(q_survey3[i] * Sel[h, 29] * Nay_M[h, i] *
          exp(rnorm(1, 0, (q_survey3[i] * Sel[h, 29] * Nay_M[h, i]) * Iay3_CV[h, i]) -
            (((q_survey3[i] * Sel[h, 29] * Nay_M[h, i] * Iay3_CV[h, i])^2) / 2)), 0) # Uncertainty based off CV

        Iay4_CV[h, i] <- runif(1, min = CV4 - 0.1, max = CV4 + 0.1)
        Iay4_F[h, i] <- max(q_survey4[i] * Sel[h, 19] * Nay_F[h, i] *
          exp(rnorm(1, 0, (q_survey4[i] * Sel[h, 19] * Nay_F[h, i]) * Iay4_CV[h, i]) -
            (((q_survey4[i] * Sel[h, 19] * Nay_F[h, i] * Iay4_CV[h, i])^2) / 2)), 0)
        Iay4_M[h, i] <- max(q_survey4[i] * Sel[h, 21] * Nay_M[h, i] *
          exp(rnorm(1, 0, (q_survey4[i] * Sel[h, 21] * Nay_M[h, i]) * Iay4_CV[h, i]) -
            (((q_survey4[i] * Sel[h, 21] * Nay_M[h, i] * Iay4_CV[h, i])^2) / 2)), 0)

        Iay5_CV[h, i] <- runif(1, min = CV5 - 0.1, max = CV5 + 0.1)
        Iay5_F[h, i] <- max(q_survey5[i] * Sel[h, 23] * Nay_F[h, i] *
          exp(rnorm(1, 0, (q_survey5[i] * Sel[h, 23] * Nay_F[h, i]) * Iay5_CV[h, i]) -
            (((q_survey5[i] * Sel[h, 23] * Nay_F[h, i] * Iay5_CV[h, i])^2) / 2)), 0)
        Iay5_M[h, i] <- max(q_survey5[i] * Sel[h, 24] * Nay_M[h, i] *
          exp(rnorm(1, 0, (q_survey5[i] * Sel[h, 24] * Nay_M[h, i]) * Iay5_CV[h, i]) -
            (((q_survey5[i] * Sel[h, 24] * Nay_M[h, i] * Iay5_CV[h, i])^2) / 2)), 0) # LOTS OF UNCERTAINTY!based off CV;

        Iay6_CV[h, i] <- runif(1, min = CV6 - 0.1, max = CV6 + 0.1)
        Iay6_F[h, i] <- max(q_survey6[i] * Sel[h, 13] * Nay_F[h, i] *
          exp(rnorm(1, 0, (q_survey6[i] * Sel[h, 13] * Nay_F[h, i]) * Iay6_CV[h, i]) -
            (((q_survey6[i] * Sel[h, 13] * Nay_F[h, i] * Iay6_CV[h, i])^2) / 2)), 0) # LOTS OF UNCERTAINTY!based off CV;
        Iay6_M[h, i] <- max(q_survey6[i] * Sel[h, 12] * Nay_M[h, i] *
          exp(rnorm(1, 0, (q_survey6[i] * Sel[h, 12] * Nay_M[h, i]) * Iay6_CV[h, i]) -
            (((q_survey6[i] * Sel[h, 12] * Nay_M[h, i] * Iay6_CV[h, i])^2) / 2)), 0)

        Iay7_CV[h, i] <- runif(1, min = CV7 - 0.1, max = CV7 + 0.1)
        Iay7_F[h, i] <- max(q_survey7[i] * Sel[h, 27] * Nay_F[h, i] *
          exp(rnorm(1, 0, (q_survey7[i] * Sel[h, 27] * Nay_F[h, i]) * Iay7_CV[h, i]) -
            (((q_survey7[i] * Sel[h, 27] * Nay_F[h, i] * Iay7_CV[h, i])^2) / 2)), 0)
        Iay7_M[h, i] <- max(q_survey7[i] * Sel[h, 26] * Nay_M[h, i] *
          exp(rnorm(1, 0, (q_survey7[i] * Sel[h, 26] * Nay_M[h, i]) * Iay7_CV[h, i]) -
            (((q_survey7[i] * Sel[h, 26] * Nay_M[h, i] * Iay7_CV[h, i])^2) / 2)), 0)
      } # end h
    } # ends i loop



    #####  GET LENGTH FREQUENCIES #####
    # USE VON BERT CURVE TO ESTIMATE LENGTH COMP.
    # Lt = Linf * (1 - exp(-K*(t - t0)))

    Linf_M <- 175.5
    K_M <- 0.143
    t0_M <- -2.388
    Linf_F <- 183.3
    K_F <- 0.124
    t0_F <- -3.098
    # L0 = Linf * (1 - exp(-K*(0-t0)))



    LFreqList <- list()

    ## Female Fishery LFreqs ##
    LFreqList[[paste("LF_Cay1F", k, sep = "_")]] <- My_LFreqs(data = Cay1_F, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90))
    LFreqList[[paste("LF_Cay2F", k, sep = "_")]] <- My_LFreqs(data = Cay2_F, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90))
    LFreqList[[paste("LF_Cay3F", k, sep = "_")]] <- My_LFreqs(data = Cay3_F, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90))
    LFreqList[[paste("LF_Cay4F", k, sep = "_")]] <- My_LFreqs(data = Cay4_F, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90))

    ## Male Fishery LFreqs ##function(data, multC = 100, Lbins = seq(30, 260, by=2), Lyrs=c(51:90), Linf=Linf_F, K=K_F, t0=t0_F)
    LFreqList[[paste("LF_Cay1M", k, sep = "_")]] <- My_LFreqs(
      data = Cay1_M, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Cay2M", k, sep = "_")]] <- My_LFreqs(
      data = Cay2_M, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Cay3M", k, sep = "_")]] <- My_LFreqs(
      data = Cay3_M, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Cay4M", k, sep = "_")]] <- My_LFreqs(
      data = Cay4_M, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )

    ## Female Survey LFreqs ##
    LFreqList[[paste("LF_Iay1F", k, sep = "_")]] <- My_LFreqs(data = Iay1_F, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:53, 59:64, 67:90))
    LFreqList[[paste("LF_Iay2F", k, sep = "_")]] <- My_LFreqs(data = Iay2_F, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(58:90))
    LFreqList[[paste("LF_Iay3F", k, sep = "_")]] <- My_LFreqs(data = Iay3_F, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(64:90))
    LFreqList[[paste("LF_Iay4F", k, sep = "_")]] <- My_LFreqs(data = Iay4_F, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(67:90))
    LFreqList[[paste("LF_Iay5F", k, sep = "_")]] <- My_LFreqs(data = Iay5_F, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(73:90))
    LFreqList[[paste("LF_Iay6F", k, sep = "_")]] <- My_LFreqs(data = Iay6_F, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(66:79))
    LFreqList[[paste("LF_Iay7F", k, sep = "_")]] <- My_LFreqs(data = Iay7_F, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(seq(68, 89, by = 3)))

    ## Male Survey LFreqs ##
    LFreqList[[paste("LF_Iay1M", k, sep = "_")]] <- My_LFreqs(
      data = Iay1_M, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:53, 59:64, 67:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Iay2M", k, sep = "_")]] <- My_LFreqs(
      data = Iay2_M, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(58:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Iay3M", k, sep = "_")]] <- My_LFreqs(
      data = Iay3_M, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(64:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Iay4M", k, sep = "_")]] <- My_LFreqs(
      data = Iay4_M, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(67:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Iay5M", k, sep = "_")]] <- My_LFreqs(
      data = Iay5_M, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(73:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Iay6M", k, sep = "_")]] <- My_LFreqs(
      data = Iay6_M, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(66:79),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Iay7M", k, sep = "_")]] <- My_LFreqs(
      data = Iay7_M, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(seq(68, 89, by = 3)),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )


    save(LFreqList, file = paste("LFreqList_", k, ".RData", sep = ""))

    #### END LFREQS #####



    Iy1[k, ] <- apply(Iay1_F, 2, sum) + apply(Iay1_M, 2, sum)
    Iy2[k, ] <- apply(Iay2_F, 2, sum) + apply(Iay2_M, 2, sum)
    Iy3[k, ] <- apply(Iay3_F, 2, sum) + apply(Iay3_M, 2, sum)
    Iy4[k, ] <- apply(Iay4_F, 2, sum) + apply(Iay4_M, 2, sum)
    Iy5[k, ] <- apply(Iay5_F, 2, sum) + apply(Iay5_M, 2, sum)
    Iy6[k, ] <- apply(Iay6_F, 2, sum) + apply(Iay6_M, 2, sum)
    Iy7[k, ] <- apply(Iay7_F, 2, sum) + apply(Iay7_M, 2, sum)

    lines(1:years, Iy1[k, ], type = "l", col = "grey")
    lines(1:years, Iy2[k, ], type = "l", col = "blue")
    lines(1:years, Iy3[k, ], type = "l", col = "red")
    lines(1:years, Iy4[k, ], type = "l", col = "green")
    lines(1:years, Iy5[k, ], type = "l", col = "orange")
    lines(1:years, Iy6[k, ], type = "l", col = "mediumorchid2")
    lines(1:years, Iy7[k, ], type = "l", col = "cyan")
    Iy1_CV[k, ] <- apply(Iay1_CV, 2, mean, na.rm = T)
    Iy2_CV[k, ] <- apply(Iay2_CV, 2, mean, na.rm = T)
    Iy3_CV[k, ] <- apply(Iay3_CV, 2, mean, na.rm = T)
    Iy4_CV[k, ] <- apply(Iay4_CV, 2, mean, na.rm = T)
    Iy5_CV[k, ] <- apply(Iay5_CV, 2, mean, na.rm = T)
    Iy6_CV[k, ] <- apply(Iay6_CV, 2, mean, na.rm = T)
    Iy7_CV[k, ] <- apply(Iay7_CV, 2, mean, na.rm = T)

    Cy_1[k, ] <- apply(Cay1_F, 2, sum, na.rm = T) + apply(Cay1_M, 2, sum, na.rm = T)
    Cy_2[k, ] <- apply(Cay2_F, 2, sum, na.rm = T) + apply(Cay2_M, 2, sum, na.rm = T)
    Cy_3[k, ] <- apply(Cay3_F, 2, sum, na.rm = T) + apply(Cay3_M, 2, sum, na.rm = T)
    Cy_4[k, ] <- apply(Cay4_F, 2, sum, na.rm = T) + apply(Cay4_M, 2, sum, na.rm = T)
    Cy[k, ] <- Cy_1[k, ] + Cy_2[k, ] + Cy_3[k, ] + Cy_4[k, ]

    Ny_F[k, ] <- apply(Nay_F, 2, sum, na.rm = T)
    Ny_M[k, ] <- apply(Nay_M, 2, sum, na.rm = T)
    Ny[k, ] <- Ny_F[k, ] + Ny_M[k, ]
    # lines(1:years, Ny[k,]/100, type='l', lwd=2)


    # write
    AnnualNC <- list()
    assign(paste("NayF", k, sep = "_"), Nay_F)
    assign(paste("NayM", k, sep = "_"), Nay_M)
    AnnualNC[[paste("NayF", k, sep = "_")]] <- Nay_F
    AnnualNC[[paste("NayM", k, sep = "_")]] <- Nay_M

    assign(paste("Cay1F", k, sep = "_"), Cay1_F)
    assign(paste("Cay2F", k, sep = "_"), Cay2_F)
    assign(paste("Cay3F", k, sep = "_"), Cay3_F)
    assign(paste("Cay4F", k, sep = "_"), Cay4_F)
    assign(paste("Cay1M", k, sep = "_"), Cay1_M)
    assign(paste("Cay2M", k, sep = "_"), Cay2_M)
    assign(paste("Cay3M", k, sep = "_"), Cay3_M)
    assign(paste("Cay4M", k, sep = "_"), Cay4_M)
    AnnualNC[[paste("Cay1F", k, sep = "_")]] <- Cay1_F
    AnnualNC[[paste("Cay2F", k, sep = "_")]] <- Cay2_F
    AnnualNC[[paste("Cay3F", k, sep = "_")]] <- Cay3_F
    AnnualNC[[paste("Cay4F", k, sep = "_")]] <- Cay4_F
    AnnualNC[[paste("Cay1M", k, sep = "_")]] <- Cay1_M
    AnnualNC[[paste("Cay2M", k, sep = "_")]] <- Cay2_M
    AnnualNC[[paste("Cay3M", k, sep = "_")]] <- Cay3_M
    AnnualNC[[paste("Cay4M", k, sep = "_")]] <- Cay4_M

    save(AnnualNC, file = paste("AnnualNC_", k, ".RData", sep = ""))
  } # end k loop


  #

  ResultsList[["Iy1"]] <- Iy1
  ResultsList[["Iy2"]] <- Iy2
  ResultsList[["Iy3"]] <- Iy3
  ResultsList[["Iy4"]] <- Iy4
  ResultsList[["Iy5"]] <- Iy5
  ResultsList[["Iy6"]] <- Iy6
  ResultsList[["Iy7"]] <- Iy7
  ResultsList[["Iy1_CV"]] <- Iy1_CV
  ResultsList[["Iy2_CV"]] <- Iy2_CV
  ResultsList[["Iy3_CV"]] <- Iy3_CV
  ResultsList[["Iy4_CV"]] <- Iy4_CV
  ResultsList[["Iy5_CV"]] <- Iy5_CV
  ResultsList[["Iy6_CV"]] <- Iy6_CV
  ResultsList[["Iy7_CV"]] <- Iy7_CV

  ResultsList[["M0"]] <- M0
  ResultsList[["FSS"]] <- FSS
  ResultsList[["Npups"]] <- Npups

  ResultsList[["Cy"]] <- Cy
  ResultsList[["Cy_1"]] <- Cy_1
  ResultsList[["Cy_2"]] <- Cy_2
  ResultsList[["Cy_3"]] <- Cy_3
  ResultsList[["Cy_4"]] <- Cy_4

  ResultsList[["Ny"]] <- Ny

  save(ResultsList, file = "ResultsList.RData")




  ########## TEST DFA #############
  # 51 - 90;
  # library(MARSS)
  yrs <- 51:90
  biasAll <- matrix(nrow = Niters, ncol = length(yrs))
  bias <- vector(length = Niters)
  RMSE <- vector(length = Niters)
  biasAll.BT <- matrix(nrow = Niters, ncol = length(yrs))
  RMSE.BT <- vector(length = Niters)

  FitRatios <- matrix(nrow = Niters, ncol = 7)
  FactorLoadings <- vector()
  DFATrends <- vector()
  DFATrendsBT <- vector()
  DFATrendsSE <- vector()
  DFATrendsSEBT <- vector()
  upCI_DFATrends <- vector()
  lowCI_DFATrends <- vector()

  datz_SD <- matrix(nrow = Niters, ncol = 7)

  # i=1

  for (i in 1:Niters) {
    I1 <- c(Iy1[i, 51:53], rep(NA, 5), Iy1[i, 59:64], NA, NA, Iy1[i, 67:90])
    I2 <- c(rep(NA, 7), Iy2[i, 58:90])
    I3 <- c(rep(NA, 13), Iy3[i, 64:90])
    I4 <- c(rep(NA, 16), Iy4[i, 67:90])
    I5 <- c(rep(NA, 22), Iy5[i, 73:90])
    I6 <- c(rep(NA, 15), Iy6[i, 66:79], rep(NA, 11))
    I7 <- c(
      rep(NA, 17), Iy7[i, 68], NA, NA, Iy7[i, 71], NA, NA, Iy7[i, 74], NA, NA, Iy7[i, 77], NA, NA, Iy7[i, 80],
      NA, NA, Iy7[i, 83], NA, NA, Iy7[i, 86], NA, NA, Iy7[i, 89], NA
    )

    assign(paste("dat", i, sep = ""), rbind(I1, I2, I3, I4, I5, I6, I7))
    # NOTE: FOR MISSING DATA SCENARIO:
    # assign(paste("dat",i,sep=""), rbind(Iy1[i,51:90], Iy2[i,51:90], Iy3[i,51:90], Iy4[i,51:90], Iy5[i,51:90], Iy6[i,51:90], Iy7[i,51:90]))  # --> for complete data

    dat.a <- get(paste("dat", i, sep = ""))

    dat <- dat.a * c

    TT <- ncol(dat)
    N.ts <- nrow(dat)
    # Standardize data
    datL <- log(dat)
    y.bar <- apply(datL, 1, mean, na.rm = TRUE)
    dat.dm <- (datL - y.bar) / y.bar
    gsd <- sd(c(dat.dm), na.rm = T)
    dat.z <- dat.dm / gsd
    y.bar
    gsd
    datz_SD[i, ] <- apply(dat.z, 1, sd, na.rm = T)
    datz_SD[i, ]
    min(datL, na.rm = T)

    dat.z <- as.matrix(dat.z)
    rownames(dat.z) <- c("Survey1", "Survey2", "Survey3", "Survey4", "Survey5", "Survey6", "Survey7")




    ##### SB DELTA LOGNORMAL ####
    cntl.list <- list(minit = 200, maxit = 500000, abstol = 0.02, allow.degen = FALSE, conv.test.slope.tol = 0.05)
    R <- diag(c(CV1, CV2, CV3, CV4, CV5, CV6, CV7), nrow = 7, ncol = 7)
    # dfa = MARSS(dat.z, model = list(m=1,R="diagonal and equal"), control=cntl.list, form="dfa", z.score=TRUE)
    dfa <- MARSS(dat.z, model = list(m = 1, R = R), control = cntl.list, form = "dfa", z.score = FALSE)

    pars <- MARSSparamCIs(dfa, nboot = 10000)


    # get the inverse of the rotation matrix
    # NOTE: retained example code if more than one common trends were estimated OR where covariates were included
    # H.inv = varimax(coef(dfa, type="matrix")$Z)$rotmat # for when more than 1 common trends
    Z.rot <- coef(dfa, type = "matrix")$Z # %*% H.inv
    trends.rot <- dfa$states
    ts.trends <- t(trends.rot)
    par.mat <- coef(dfa, type = "matrix")
    fit.b <- par.mat$Z %*% dfa$states # + matrix(par.mat$D, nrow=N.ts) %*% covar # for when covariates included

    assign(paste("Z.rot", i, sep = "."), Z.rot)
    assign(paste("Z.upCI", i, sep = "."), pars$par.upCI$Z)
    assign(paste("Z.lowCI", i, sep = "."), pars$par.lowCI$Z)
    assign(paste("trends.rot", i, sep = "."), trends.rot)


    # plot factor loadings
    survey <- rownames(dat.z)
    minZ <- 0.00
    ylims <- c(-1.1 * max(abs(Z.rot)), 1.1 * max(abs(Z.rot)))
    par(mfrow = c(nrow(trends.rot), 1), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
    for (h in 1:nrow(trends.rot)) {
      plot(c(1:N.ts)[abs(Z.rot[, h]) > minZ], as.vector(Z.rot[abs(Z.rot[, h]) > minZ, h]),
        type = "h", lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylims, xlim = c(0, N.ts + 1)
      )
      for (j in 1:N.ts) {
        if (Z.rot[j, h] > minZ) {
          text(j, -0.05, survey[j], srt = 90, adj = 1, cex = 0.9)
        }
        if (Z.rot[j, h] < -minZ) {
          text(j, 0.05, survey[j], srt = 90, adj = 0, cex = 0.9)
        }
        abline(h = 0, lwd = 1, col = "gray")
        abline(h = 0.2, col = "gray", lty = 2)
        abline(h = -0.2, col = "gray", lty = 2)
      } # end j loop
      mtext(paste("Factor loadings on trend", h, sep = " "), side = 3, line = .5)
    } # end h loop


    # PLOT WITH CI's

    index <- dfa$states[1, ]
    indexSE <- dfa$states.se[1, ]
    lowerCI <- dfa$states[1, ] - 1.96 * dfa$states.se[1, ]
    upperCI <- dfa$states[1, ] + 1.96 * dfa$states.se[1, ]
    assign(paste("index", i, sep = ""), index)
    assign(paste("indexSE", i, sep = ""), indexSE)
    assign(paste("lowerCI", i, sep = ""), lowerCI)
    assign(paste("upperCI", i, sep = ""), upperCI)



    par(mfrow = c(4, 2))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type = "n", bty = "L", xlab = "Year", ylab = "Abundance", ylim = c(min(lowerCI), max(upperCI))))
    x <- as.numeric(df$yrs)
    polygon(c(x, rev(x)), c(df$lowerCI, rev(df$upperCI)), col = "grey", border = NA)
    with(df, lines(x, index, type = "l", lwd = 2, pch = 16))
    title(paste("Iteration", i, sep = " "))
    abline(h = 0)


    N <- Ny[i, 51:90]
    Nscale <- (N - mean(N)) / sd(N)
    NL <- log(N)
    NLscale <- (NL - mean(NL)) / sd(NL)
    index.z <- (index - mean(index)) / sd(index)
    index.z <- ifelse(is.na(index.z), 0, index.z)
    index.z <- ifelse(abs(index.z) == Inf, 0, index.z)


    indexSEBT <- indexSE * gsd
    assign(paste("indexSEBT", i, sep = ""), indexSEBT)

    indexBT <- exp(index * gsd + ((indexSEBT^2) / 2))
    assign(paste("indexBT", i, sep = ""), indexBT)



    lines(yrs, rescale(NL, to = c(min(index), max(index))), type = "l", col = "red", lwd = 2)
    assign(paste("N", i, sep = ""), N)
    assign(paste("Nscale", i, sep = ""), Nscale)
    assign(paste("NL", i, sep = ""), NL)
    assign(paste("NLscale", i, sep = ""), NLscale)


    biasAll[i, ] <- index.z - NLscale
    bias[i] <- mean(index.z - NLscale)
    RMSE[i] <- sqrt(sum((index.z - NLscale)^2) / length(index.z))

    indexBT.z <- (indexBT - mean(indexBT)) / sd(indexBT)
    biasAll.BT[i, ] <- indexBT.z - Nscale
    RMSE.BT[i] <- sqrt(sum((indexBT.z - Nscale)^2) / length(indexBT.z))


    # Plot fitted values
    survey <- rownames(dat.z)
    # par(mfcol=c(ceiling(N.ts/2),2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
    for (n in 1:length(survey)) {
      plot(yrs, dat.z[n, ], xlab = "", ylab = "abundance index", bty = "L", xaxt = "n", ylim = c(-4, 4), pch = 16, col = "blue")
      axis(1, labels = TRUE)
      lines(yrs, fit.b[n, ], lwd = 2)
      lines(yrs, rescale(NL, to = c(min(fit.b[n, ]), max(fit.b[n, ]))), type = "l", col = "red", lwd = 2)
      title(paste("Sandbar", survey[n], sep = " "))
    }


    #### ASSESS MODEL FITS
    sumResids <- rowSums((dat.z - fit.b)^2, na.rm = TRUE)
    sumObserved <- rowSums(dat.z^2, na.rm = TRUE)
    FitRatio <- sumResids / sumObserved
    FitRatio
    FitRatios[i, ] <- FitRatio

    # mean(FitRatio)
    # indexBT
    # dat.a


    FactorLoadings <- rbind(FactorLoadings, c(get(paste("Z.rot", i, sep = ".")), get(paste("Z.lowCI", i, sep = ".")), get(paste("Z.upCI", i, sep = "."))))
    DFATrends <- rbind(DFATrends, index)
    DFATrendsBT <- rbind(DFATrendsBT, indexBT)
    DFATrendsSE <- rbind(DFATrendsSE, indexSE)
    DFATrendsSEBT <- rbind(DFATrendsSEBT, indexSEBT)
    upCI_DFATrends <- rbind(upCI_DFATrends, upperCI)
    lowCI_DFATrends <- rbind(lowCI_DFATrends, lowerCI)


    # par(mfrow=c(1,1))
    # plot(fit.b, dat.z-fit.b)
  } # END LOOP

  colnames(FactorLoadings) <- c(
    "FL1", "FL2", "FL3", "FL4", "FL5", "FL6", "FL7",
    "lowCI_FL1", "lowCI_FL2", "lowCI_FL3", "lowCI_FL4", "lowCI_FL5", "lowCI_FL6", "lowCI_FL7",
    "upCI_FL1", "upCI_FL2", "upCI_FL3", "upCI_FL4", "upCI_FL5", "upCI_FL6", "upCI_FL7"
  )
  FR <- apply(FitRatios, 1, mean)
  FLoadings <- cbind(FactorLoadings, meanFactorLoading = FR)

  DFA_Results <- list()
  DFA_Results[["FitRatios"]] <- FitRatios
  DFA_Results[["DFATrends"]] <- DFATrends
  DFA_Results[["DFATrendsBT"]] <- DFATrendsBT
  DFA_Results[["DFATrendsSE"]] <- DFATrendsSE
  DFA_Results[["DFATrendsSEBT"]] <- DFATrendsSEBT
  DFA_Results[["upCI_DFATrends"]] <- upCI_DFATrends
  DFA_Results[["lowCI_DFATrends"]] <- lowCI_DFATrends
  DFA_Results[["FLoadings"]] <- FLoadings
  DFA_Results[["bias"]] <- bias
  DFA_Results[["biasAll"]] <- biasAll
  DFA_Results[["RMSE"]] <- RMSE
  DFA_Results[["biasAll.BT"]] <- biasAll.BT
  DFA_Results[["RMSE.BT"]] <- RMSE.BT
  DFA_Results[["GSD"]] <- gsd

  DFA_Results[["datz_SD"]] <- datz_SD

  save(DFA_Results, file = "DFA_Results.RData")




  ###################################  plotting ###################################
  if (plot == TRUE) {
    png(
      filename = "Trend.png",
      type = "cairo",
      units = "mm",
      width = 300,
      height = 300,
      pointsize = 12,
      res = 300
    )

    par(mfrow = c(10, 10))
    par(mar = c(0.3, 0.3, 0.1, 0.1), tcl = -0.1, mgp = c(2, 0.6, 0))
    for (l in 1:Niters) {
      index <- get(paste("index", l, sep = ""))
      lowerCI <- get(paste("lowerCI", l, sep = ""))
      upperCI <- get(paste("upperCI", l, sep = ""))
      df <- data.frame(yrs, index, lowerCI, upperCI)
      with(df, plot(index ~ yrs, type = "n", bty = "L", xlab = "", ylab = "", ylim = c(min(lowerCI), max(upperCI)), axes = FALSE))
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      x <- as.numeric(df$yrs)
      polygon(c(x, rev(x)), c(df$lowerCI, rev(df$upperCI)), col = "grey", border = NA)
      with(df, lines(x, index, type = "l", lwd = 2, pch = 16))
      # title(paste("Iteration",l,sep=" "))
      abline(h = 0)
      NLscale <- get(paste("NLscale", l, sep = ""))
      lines(yrs, NLscale, type = "l", col = "red", lwd = 2)
      text(x = min(yrs) + 4, y = 0.8 * max(min(lowerCI)), labels = format(FR[l], digits = 3), cex = 1)
    }

    dev.off()


    png(
      filename = "RawIndices.png",
      type = "cairo",
      units = "mm",
      width = 300,
      height = 300,
      pointsize = 12,
      res = 300
    )

    par(mfrow = c(10, 10))
    par(mar = c(0.3, 0.3, 0.1, 0.1), tcl = -0.1, mgp = c(2, 0.6, 0))
    for (l in 1:Niters) {
      plot(51:90, c(Iy1[l, 51:53], rep(NA, 5), Iy1[l, 59:64], NA, NA, Iy1[l, 67:90]),
        type = "l", col = "grey", xlim = c(51, 90), axes = F,
        ylim = c(0, max(c(Iy1[l, ], Iy2[l, ], Iy3[l, ], Iy4[l, ], Iy5[l, ], Iy6[l, ], Iy7[l, ]), na.rm = T) + 0.1)
      )
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      lines(51:90, c(rep(NA, 7), Iy2[l, 58:90]), type = "l", col = "blue")
      lines(51:90, c(rep(NA, 13), Iy3[l, 64:90]), type = "l", col = "red")
      lines(51:90, c(rep(NA, 16), Iy4[l, 67:90]), type = "l", col = "green")
      lines(51:90, c(rep(NA, 22), Iy5[l, 73:90]), type = "l", col = "orange")
      lines(51:90, c(rep(NA, 15), Iy6[l, 66:79], rep(NA, 11)), type = "l", col = "mediumorchid2")
      points(51:90, c(
        rep(NA, 17), Iy7[l, 68], NA, NA, Iy7[l, 71], NA, NA, Iy7[l, 74], NA, NA, Iy7[l, 77], NA, NA, Iy7[l, 80],
        NA, NA, Iy7[l, 83], NA, NA, Iy7[l, 86], NA, NA, Iy7[l, 89], NA
      ), type = "l", col = "cyan")
      # title(paste("Iteration",l, sep=" "))
    }

    dev.off()


    png(
      filename = "FactorLoadings.png",
      type = "cairo",
      units = "mm",
      width = 300,
      height = 300,
      pointsize = 12,
      res = 300
    )

    # Before plotting factor loadings, we'll need to assign Z.rot, trends.rot
    par(mfrow = c(10, 10))
    par(mar = c(0.3, 0.3, 0.1, 0.1), tcl = -0.1, mgp = c(2, 0.6, 0))
    survey <- rownames(dat.z)
    minZ <- 0.00
    for (l in 1:Niters) {
      Z.rot <- get(paste("Z.rot", l, sep = "."))
      trends.rot <- get(paste("trends.rot", l, sep = "."))
      ylims <- c(-1.1 * max(abs(Z.rot)), 1.1 * max(abs(Z.rot)))
      # par(mfrow=c(nrow(trends.rot),1), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
      for (h in 1:nrow(trends.rot)) {
        plot(c(1:N.ts)[abs(Z.rot[, h]) > minZ], as.vector(Z.rot[abs(Z.rot[, h]) > minZ, h]),
          type = "h", lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylims, xlim = c(0, N.ts + 1), axes = FALSE
        )
        axis(1, labels = FALSE)
        axis(2, labels = FALSE)
        abline(h = 0, lwd = 1, col = "gray")
        abline(h = 0.2, col = "gray", lty = 2)
        abline(h = -0.2, col = "gray", lty = 2)
      } # end h loop
    }

    dev.off()



    png(
      filename = "Trend_scale.png",
      type = "cairo",
      units = "mm",
      width = 300,
      height = 300,
      pointsize = 12,
      res = 300
    )


    par(mfrow = c(10, 10))
    par(mar = c(0.3, 0.3, 0.1, 0.1), tcl = -0.1, mgp = c(2, 0.6, 0))
    for (l in 1:Niters) {
      index <- get(paste("index", l, sep = ""))
      lowerCI <- get(paste("lowerCI", l, sep = ""))
      upperCI <- get(paste("upperCI", l, sep = ""))
      df <- data.frame(yrs, index, lowerCI, upperCI)
      with(df, plot(index ~ yrs, type = "n", bty = "L", xlab = "", ylab = "", ylim = c(min(lowerCI), max(upperCI)), axes = FALSE))
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      x <- as.numeric(df$yrs)
      polygon(c(x, rev(x)), c(df$lowerCI, rev(df$upperCI)), col = "grey", border = NA)
      with(df, lines(x, index, type = "l", lwd = 2, pch = 16))
      # title(paste("Iteration",l,sep=" "))
      abline(h = 0)
      NL <- get(paste("NL", l, sep = ""))
      lines(yrs, rescale(NL, to = c(min(index), max(index))), type = "l", col = "red", lwd = 2)
      text(x = min(yrs) + 4, y = 0.8 * min(lowerCI), labels = format(FR[l], digits = 3), cex = 1)
    }


    dev.off()




    png(
      filename = "Trend_N.png",
      type = "cairo",
      units = "mm",
      width = 300,
      height = 300,
      pointsize = 12,
      res = 300
    )

    par(mfrow = c(10, 10))
    par(mar = c(0.3, 0.3, 0.1, 0.1), tcl = -0.1, mgp = c(2, 0.6, 0))
    for (l in 1:Niters) {
      index <- get(paste("index", l, sep = ""))
      NL <- get(paste("NL", l, sep = ""))
      plot(yrs, N, type = "l", col = "red", axes = FALSE)
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      lines(yrs, rescale(index, to = c(min(N), max(N))))
      text(x = min(yrs) + 4, y = 1.05 * min(N), labels = format(FR[l], digits = 3), cex = 1)
    }

    dev.off()





    png(
      filename = "BackTransformed_Trend_scale.png",
      type = "cairo",
      units = "mm",
      width = 300,
      height = 300,
      pointsize = 12,
      res = 300
    )


    par(mfrow = c(10, 10))
    par(mar = c(0.3, 0.3, 0.1, 0.1), tcl = -0.1, mgp = c(2, 0.6, 0))
    for (l in 1:Niters) {
      indexBT <- get(paste("indexBT", l, sep = ""))
      indexSEBT <- get(paste("indexSEBT", l, sep = ""))
      lowerCI <- indexBT - (1.96 * indexSEBT)
      upperCI <- indexBT + (1.96 * indexSEBT)
      df <- data.frame(yrs, indexBT, lowerCI, upperCI)
      with(df, plot(index ~ yrs, type = "n", bty = "L", xlab = "", ylab = "", ylim = c(min(lowerCI), max(upperCI)), axes = FALSE))
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      x <- as.numeric(df$yrs)
      polygon(c(x, rev(x)), c(df$lowerCI, rev(df$upperCI)), col = "grey", border = NA)
      with(df, lines(x, indexBT, type = "l", lwd = 2, pch = 16))
      # title(paste("Iteration",l,sep=" "))
      # abline(h=0)
      N <- get(paste("N", l, sep = ""))
      lines(yrs, rescale(N, to = c(min(indexBT), max(indexBT))), type = "l", col = "red", lwd = 2)
    }


    dev.off()



    png(
      filename = "Trend_Standardize.png",
      type = "cairo",
      units = "mm",
      width = 300,
      height = 300,
      pointsize = 12,
      res = 300
    )

    par(mfrow = c(10, 10))
    par(mar = c(0.3, 0.3, 0.1, 0.1), tcl = -0.1, mgp = c(2, 0.6, 0))
    for (l in 1:Niters) {
      index <- get(paste("index", l, sep = ""))
      NL <- get(paste("NL", l, sep = ""))
      plot(yrs, (NL - mean(NL)) / sd(NL), type = "l", col = "red", axes = FALSE, lwd = 1.5)
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      lines(yrs, (index - mean(index)) / sd(index), type = "l", lwd = 1.5)
      text(x = min(yrs) + 4, y = 0.8 * min((index - mean(index)) / sd(index)), labels = format(FR[l], digits = 3), cex = 1)
    }

    dev.off()
  } # end if plot == TRUE


  ##################################################### END #################################################

  return(list(RMSE = RMSE, biasAll = biasAll, FitRatios = FitRatios, RMSE.BT = RMSE.BT, biasAll.BT = biasAll.BT))
} # END FUNCTION


# Complete data function #

DFA_Sim_Full <- function(Tot_F, CV1, CV2, CV3, CV4, CV5, CV6, CV7,
                         q_survey1, q_survey2, q_survey3, q_survey4, q_survey5, q_survey6, q_survey7,
                         Iay_eq = Iay_eq,
                         Niters = 100, c = 5, plot = TRUE) {

  # Inputs:
  M <- M_const
  # GET F
  Prop_F <- c(0.2, 0.044, 0.754, 0.002)
  F_const1 <- Tot_F * Prop_F[1]
  F_const2 <- Tot_F * Prop_F[2]
  F_const3 <- Tot_F * Prop_F[3]
  F_const4 <- Tot_F * Prop_F[4]
  FM <- rbind(F_const1, F_const2, F_const3, F_const4)
  My_LFreqs <- function(data, multC = 100, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90), Linf = Linf_F, K = K_F, t0 = t0_F) {
    data_N <- cbind(ages, round(data * multC))
    LF <- data.frame(Lbins)
    for (l in Lyrs) {
      df <- as.data.frame(data_N[, c(1, l + 1)])
      colnames(df) <- c("ages", "count")
      df.expanded <- df[rep(row.names(df), df$count), 1:2]
      vonBert <- rbvn(nrow(df.expanded), m1 = Linf, s1 = Linf * 0.1, m2 = K, s2 = K * 0.1, rho = -0.9)
      df.expanded$L <- vonBert$Linf * (1 - exp(-(vonBert$K) * (jitter(df.expanded$ages, amount = 0.5) - t0)))
      hst <- hist(df.expanded$L, breaks = c(seq(30, 262, by = 2)), plot = F)
      LF <- cbind(LF, hst$counts)
    }
    colnames(LF) <- c("LBins", as.character(Lyrs))
    return(LF)
  }


  # Iay
  par(mfrow = c(1, 1))
  plot(1:years, apply(Iay_eq, 2, sum), type = "l", ylim = c(0, max(apply(Iay_eq, 2, sum), na.rm = T) * 3))
  abline(v = 51, lwd = 2, col = "red")
  abline(v = 90, lwd = 2, col = "red")
  # plot(1:years, apply(Iay, 2, sum), type='l')
  # Iy = apply(Iay, 2, sum)


  ####### ADD STOCHASTICITY ###########
  Iay1_F <- matrix(nrow = A, ncol = years)
  Iay2_F <- matrix(nrow = A, ncol = years)
  Iay3_F <- matrix(nrow = A, ncol = years)
  Iay4_F <- matrix(nrow = A, ncol = years)
  Iay5_F <- matrix(nrow = A, ncol = years)
  Iay6_F <- matrix(nrow = A, ncol = years)
  Iay7_F <- matrix(nrow = A, ncol = years)
  Iay1_M <- matrix(nrow = A, ncol = years)
  Iay2_M <- matrix(nrow = A, ncol = years)
  Iay3_M <- matrix(nrow = A, ncol = years)
  Iay4_M <- matrix(nrow = A, ncol = years)
  Iay5_M <- matrix(nrow = A, ncol = years)
  Iay6_M <- matrix(nrow = A, ncol = years)
  Iay7_M <- matrix(nrow = A, ncol = years)
  Iay1_CV <- matrix(nrow = A, ncol = years)
  Iay2_CV <- matrix(nrow = A, ncol = years)
  Iay3_CV <- matrix(nrow = A, ncol = years)
  Iay4_CV <- matrix(nrow = A, ncol = years)
  Iay5_CV <- matrix(nrow = A, ncol = years)
  Iay6_CV <- matrix(nrow = A, ncol = years)
  Iay7_CV <- matrix(nrow = A, ncol = years)
  Iy1 <- matrix(nrow = Niters, ncol = years)
  Iy2 <- matrix(nrow = Niters, ncol = years)
  Iy3 <- matrix(nrow = Niters, ncol = years)
  Iy4 <- matrix(nrow = Niters, ncol = years)
  Iy5 <- matrix(nrow = Niters, ncol = years)
  Iy6 <- matrix(nrow = Niters, ncol = years)
  Iy7 <- matrix(nrow = Niters, ncol = years)
  Iy1_CV <- matrix(nrow = Niters, ncol = years)
  Iy2_CV <- matrix(nrow = Niters, ncol = years)
  Iy3_CV <- matrix(nrow = Niters, ncol = years)
  Iy4_CV <- matrix(nrow = Niters, ncol = years)
  Iy5_CV <- matrix(nrow = Niters, ncol = years)
  Iy6_CV <- matrix(nrow = Niters, ncol = years)
  Iy7_CV <- matrix(nrow = Niters, ncol = years)
  Cy_1 <- matrix(nrow = Niters, ncol = years)
  Cy_2 <- matrix(nrow = Niters, ncol = years)
  Cy_3 <- matrix(nrow = Niters, ncol = years)
  Cy_4 <- matrix(nrow = Niters, ncol = years)
  Ny_F <- matrix(nrow = Niters, ncol = years)
  Ny_M <- matrix(nrow = Niters, ncol = years)
  Cy <- matrix(nrow = Niters, ncol = years)
  Ny <- matrix(nrow = Niters, ncol = years)


  M0 <- matrix(nrow = Niters, ncol = years)
  Npups <- matrix(nrow = Niters, ncol = years)
  FSS <- matrix(nrow = Niters, ncol = years) # Female spawning stock

  ResultsList <- list()


  Z0 <- 1.085716
  Zmin <- 0.1604
  zfrac <- 0.8522633
  beta <- 0.3582577
  Npups_eq <- 2961.558 # Npups0


  set.seed(430)

  # F_const = rep(0, yrs)

  for (k in 1:Niters) {
    # paste("Iteration",k,sep=" ")
    # Create matrices to store results
    Nay_F <- matrix(nrow = A, ncol = years)
    Nay_M <- matrix(nrow = A, ncol = years)
    Cay1_F <- matrix(nrow = A, ncol = years)
    Cay2_F <- matrix(nrow = A, ncol = years)
    Cay3_F <- matrix(nrow = A, ncol = years)
    Cay4_F <- matrix(nrow = A, ncol = years)
    Cay1_M <- matrix(nrow = A, ncol = years)
    Cay2_M <- matrix(nrow = A, ncol = years)
    Cay3_M <- matrix(nrow = A, ncol = years)
    Cay4_M <- matrix(nrow = A, ncol = years)



    # Set up equilibrium conditions
    Nay_F[, 1] <- Equilibrium_Nay_F
    Nay_M[, 1] <- Equilibrium_Nay_M

    for (i in 1:(years - 1)) {
      # i=i+1

      matjit <- rnorm(length(mat_F), mat_F, mat_F * 0.01)
      matjit <- ifelse(matjit > 1, 1, matjit)
      FSS[k, i] <- sum(Nay_F[, i] * matjit) # Fecund Stock size
      Npups[k, i] <- sum((Nay_F[, i] / 2.5) * matjit * rnorm(length(fec), fec, fec * 0.1)) #
      M0[k, i] <- rnorm(1, (((1 - (Npups[k, i] / Npups_eq)^beta) * (Zmin - Z0)) + Z0), 0.1)
      # M0[i] = ( (1- (Npups[i]/Npups_eq)^beta ) * (Zmin - Z0) ) + Z0


      Nay_F[1, i + 1] <- 0.5 * Npups[k, i] * exp(-M0[k, i]) # pups produced that survive to age 1; Recruits
      Nay_M[1, i + 1] <- 0.5 * Npups[k, i] * exp(-M0[k, i]) # pups produced that survive to age 1; Recruits


      Zay_F <- vector(length = length(M))
      Zay_M <- vector(length = length(M))
      for (j in 1:A) {
        Zay_F[j] <- M[j] + (Sel[j, 2] * FM[1, i]) + (Sel[j, 4] * FM[2, i]) + (Sel[j, 6] * FM[3, i]) + (Sel[j, 8] * FM[4, i])
        Zay_M[j] <- M[j] + (Sel[j, 3] * FM[1, i]) + (Sel[j, 5] * FM[2, i]) + (Sel[j, 7] * FM[3, i]) + (Sel[j, 8] * FM[4, i])
      } # END J LOOP

      for (j in 1:(A - 2)) {
        Mjit <- rnorm(length(M), M, M * 0.1)
        Nay_F[j + 1, i + 1] <- Nay_F[j, i] * exp(-Zay_F[j])
        Nay_M[j + 1, i + 1] <- Nay_M[j, i] * exp(-Zay_M[j])


        Cay1_F[j, i] <- Nay_F[j, i] * ((Sel[j, 2] * FM[1, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
        Cay2_F[j, i] <- Nay_F[j, i] * ((Sel[j, 4] * FM[2, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
        Cay3_F[j, i] <- Nay_F[j, i] * ((Sel[j, 6] * FM[3, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
        Cay4_F[j, i] <- Nay_F[j, i] * ((Sel[j, 8] * FM[4, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))

        Cay1_M[j, i] <- Nay_M[j, i] * ((Sel[j, 3] * FM[1, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
        Cay2_M[j, i] <- Nay_M[j, i] * ((Sel[j, 5] * FM[2, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
        Cay3_M[j, i] <- Nay_M[j, i] * ((Sel[j, 7] * FM[3, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
        Cay4_M[j, i] <- Nay_M[j, i] * ((Sel[j, 8] * FM[4, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
      } # END J LOOP

      Nay_F[A, i + 1] <- Nay_F[A - 1, i] * exp(-Zay_F[A - 1]) #+ Nay_F[A,i]*exp(-Zay[A])
      Nay_M[A, i + 1] <- Nay_M[A - 1, i] * exp(-Zay_M[A - 1]) #+ Nay_M[A,i]*exp(-Zay[A])

      Cay1_F[A - 1, i] <- Nay_F[A - 1, i] * ((Sel[A - 1, 2] * FM[1, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
      Cay2_F[A - 1, i] <- Nay_F[A - 1, i] * ((Sel[A - 1, 4] * FM[2, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
      Cay3_F[A - 1, i] <- Nay_F[A - 1, i] * ((Sel[A - 1, 6] * FM[3, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
      Cay4_F[A - 1, i] <- Nay_F[A - 1, i] * ((Sel[A - 1, 8] * FM[4, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))

      Cay1_M[A - 1, i] <- Nay_M[A - 1, i] * ((Sel[A - 1, 3] * FM[1, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
      Cay2_M[A - 1, i] <- Nay_M[A - 1, i] * ((Sel[A - 1, 5] * FM[2, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
      Cay3_M[A - 1, i] <- Nay_M[A - 1, i] * ((Sel[A - 1, 7] * FM[3, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
      Cay4_M[A - 1, i] <- Nay_M[A - 1, i] * ((Sel[A - 1, 8] * FM[4, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))



      Cay1_F[A, i] <- Nay_F[A, i] * ((Sel[A, 2] * FM[1, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
      Cay2_F[A, i] <- Nay_F[A, i] * ((Sel[A, 4] * FM[2, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
      Cay3_F[A, i] <- Nay_F[A, i] * ((Sel[A, 6] * FM[3, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))
      Cay4_F[A, i] <- Nay_F[A, i] * ((Sel[A, 8] * FM[4, i]) / Zay_F[j]) * (1 - exp(-Zay_F[j]))

      Cay1_M[A, i] <- Nay_M[A, i] * ((Sel[A, 3] * FM[1, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
      Cay2_M[A, i] <- Nay_M[A, i] * ((Sel[A, 5] * FM[2, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
      Cay3_M[A, i] <- Nay_M[A, i] * ((Sel[A, 7] * FM[3, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))
      Cay4_M[A, i] <- Nay_M[A, i] * ((Sel[A, 8] * FM[4, i]) / Zay_M[j]) * (1 - exp(-Zay_M[j]))


      for (h in 1:A) {
        # h=h+1
        Iay1_CV[h, i] <- runif(1, min = CV1 - 0.1, max = CV1 + 0.1)
        Iay1_F[h, i] <- max(q_survey1[i] * Sel[h, 17] * Nay_F[h, i] *
          exp(rnorm(1, 0, (q_survey1[i] * Sel[h, 17] * Nay_F[h, i]) * Iay1_CV[h, i]) -
            (((q_survey1[i] * Sel[h, 17] * Nay_F[h, i] * Iay1_CV[h, i])^2) / 2)), 0) # uncertainty based off CV;
        Iay1_M[h, i] <- max(q_survey1[i] * Sel[h, 18] * Nay_M[h, i] *
          exp(rnorm(1, 0, (q_survey1[i] * Sel[h, 18] * Nay_M[h, i]) * Iay1_CV[h, i]) -
            (((q_survey1[i] * Sel[h, 18] * Nay_M[h, i] * Iay1_CV[h, i])^2) / 2)), 0)

        Iay2_CV[h, i] <- runif(1, min = CV2 - 0.1, max = CV2 + 0.1)
        Iay2_F[h, i] <- max(q_survey2[i] * Sel[h, 9] * Nay_F[h, i] *
          exp(rnorm(1, 0, (q_survey2[i] * Sel[h, 9] * Nay_F[h, i]) * Iay2_CV[h, i]) -
            (((q_survey2[i] * Sel[h, 9] * Nay_F[h, i] * Iay2_CV[h, i])^2) / 2)), 0)
        Iay2_M[h, i] <- max(q_survey2[i] * Sel[h, 10] * Nay_M[h, i] *
          exp(rnorm(1, 0, (q_survey2[i] * Sel[h, 10] * Nay_M[h, i]) * Iay2_CV[h, i]) -
            (((q_survey2[i] * Sel[h, 10] * Nay_M[h, i] * Iay2_CV[h, i])^2) / 2)), 0)

        Iay3_CV[h, i] <- runif(1, min = CV3 - 0.1, max = CV3 + 0.1)
        Iay3_F[h, i] <- max(q_survey3[i] * Sel[h, 28] * Nay_F[h, i] *
          exp(rnorm(1, 0, (q_survey3[i] * Sel[h, 28] * Nay_F[h, i]) * Iay3_CV[h, i]) -
            (((q_survey3[i] * Sel[h, 28] * Nay_F[h, i] * Iay3_CV[h, i])^2) / 2)), 0) #
        Iay3_M[h, i] <- max(q_survey3[i] * Sel[h, 29] * Nay_M[h, i] *
          exp(rnorm(1, 0, (q_survey3[i] * Sel[h, 29] * Nay_M[h, i]) * Iay3_CV[h, i]) -
            (((q_survey3[i] * Sel[h, 29] * Nay_M[h, i] * Iay3_CV[h, i])^2) / 2)), 0)

        Iay4_CV[h, i] <- runif(1, min = CV4 - 0.1, max = CV4 + 0.1)
        Iay4_F[h, i] <- max(q_survey4[i] * Sel[h, 19] * Nay_F[h, i] *
          exp(rnorm(1, 0, (q_survey4[i] * Sel[h, 19] * Nay_F[h, i]) * Iay4_CV[h, i]) -
            (((q_survey4[i] * Sel[h, 19] * Nay_F[h, i] * Iay4_CV[h, i])^2) / 2)), 0)
        Iay4_M[h, i] <- max(q_survey4[i] * Sel[h, 21] * Nay_M[h, i] *
          exp(rnorm(1, 0, (q_survey4[i] * Sel[h, 21] * Nay_M[h, i]) * Iay4_CV[h, i]) -
            (((q_survey4[i] * Sel[h, 21] * Nay_M[h, i] * Iay4_CV[h, i])^2) / 2)), 0)

        Iay5_CV[h, i] <- runif(1, min = CV5 - 0.1, max = CV5 + 0.1)
        Iay5_F[h, i] <- max(q_survey5[i] * Sel[h, 23] * Nay_F[h, i] *
          exp(rnorm(1, 0, (q_survey5[i] * Sel[h, 23] * Nay_F[h, i]) * Iay5_CV[h, i]) -
            (((q_survey5[i] * Sel[h, 23] * Nay_F[h, i] * Iay5_CV[h, i])^2) / 2)), 0)
        Iay5_M[h, i] <- max(q_survey5[i] * Sel[h, 24] * Nay_M[h, i] *
          exp(rnorm(1, 0, (q_survey5[i] * Sel[h, 24] * Nay_M[h, i]) * Iay5_CV[h, i]) -
            (((q_survey5[i] * Sel[h, 24] * Nay_M[h, i] * Iay5_CV[h, i])^2) / 2)), 0)

        Iay6_CV[h, i] <- runif(1, min = CV6 - 0.1, max = CV6 + 0.1)
        Iay6_F[h, i] <- max(q_survey6[i] * Sel[h, 13] * Nay_F[h, i] *
          exp(rnorm(1, 0, (q_survey6[i] * Sel[h, 13] * Nay_F[h, i]) * Iay6_CV[h, i]) -
            (((q_survey6[i] * Sel[h, 13] * Nay_F[h, i] * Iay6_CV[h, i])^2) / 2)), 0)
        Iay6_M[h, i] <- max(q_survey6[i] * Sel[h, 12] * Nay_M[h, i] *
          exp(rnorm(1, 0, (q_survey6[i] * Sel[h, 12] * Nay_M[h, i]) * Iay6_CV[h, i]) -
            (((q_survey6[i] * Sel[h, 12] * Nay_M[h, i] * Iay6_CV[h, i])^2) / 2)), 0)

        Iay7_CV[h, i] <- runif(1, min = CV7 - 0.1, max = CV7 + 0.1)
        Iay7_F[h, i] <- max(q_survey7[i] * Sel[h, 27] * Nay_F[h, i] *
          exp(rnorm(1, 0, (q_survey7[i] * Sel[h, 27] * Nay_F[h, i]) * Iay7_CV[h, i]) -
            (((q_survey7[i] * Sel[h, 27] * Nay_F[h, i] * Iay7_CV[h, i])^2) / 2)), 0)
        Iay7_M[h, i] <- max(q_survey7[i] * Sel[h, 26] * Nay_M[h, i] *
          exp(rnorm(1, 0, (q_survey7[i] * Sel[h, 26] * Nay_M[h, i]) * Iay7_CV[h, i]) -
            (((q_survey7[i] * Sel[h, 26] * Nay_M[h, i] * Iay7_CV[h, i])^2) / 2)), 0)
      } # end h
    } # ends i loop



    #####  GET LENGTH FREQUENCIES #####
    # USE VON BERT CURVE TO ESTIMATE LENGTH COMP.
    # Lt = Linf * (1 - exp(-K*(t - t0)))

    Linf_M <- 175.5
    K_M <- 0.143
    t0_M <- -2.388
    Linf_F <- 183.3
    K_F <- 0.124
    t0_F <- -3.098
    # L0 = Linf * (1 - exp(-K*(0-t0)))



    LFreqList <- list()

    ## Female Fishery LFreqs ##
    LFreqList[[paste("LF_Cay1F", k, sep = "_")]] <- My_LFreqs(data = Cay1_F, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90))
    LFreqList[[paste("LF_Cay2F", k, sep = "_")]] <- My_LFreqs(data = Cay2_F, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90))
    LFreqList[[paste("LF_Cay3F", k, sep = "_")]] <- My_LFreqs(data = Cay3_F, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90))
    LFreqList[[paste("LF_Cay4F", k, sep = "_")]] <- My_LFreqs(data = Cay4_F, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90))

    ## Male Fishery LFreqs ##
    LFreqList[[paste("LF_Cay1M", k, sep = "_")]] <- My_LFreqs(
      data = Cay1_M, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Cay2M", k, sep = "_")]] <- My_LFreqs(
      data = Cay2_M, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Cay3M", k, sep = "_")]] <- My_LFreqs(
      data = Cay3_M, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Cay4M", k, sep = "_")]] <- My_LFreqs(
      data = Cay4_M, multC = 10, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )

    ## Female Survey LFreqs ##
    LFreqList[[paste("LF_Iay1F", k, sep = "_")]] <- My_LFreqs(data = Iay1_F, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90))
    LFreqList[[paste("LF_Iay2F", k, sep = "_")]] <- My_LFreqs(data = Iay2_F, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90))
    LFreqList[[paste("LF_Iay3F", k, sep = "_")]] <- My_LFreqs(data = Iay3_F, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90))
    LFreqList[[paste("LF_Iay4F", k, sep = "_")]] <- My_LFreqs(data = Iay4_F, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90))
    LFreqList[[paste("LF_Iay5F", k, sep = "_")]] <- My_LFreqs(data = Iay5_F, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90))
    LFreqList[[paste("LF_Iay6F", k, sep = "_")]] <- My_LFreqs(data = Iay6_F, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90))
    LFreqList[[paste("LF_Iay7F", k, sep = "_")]] <- My_LFreqs(data = Iay7_F, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90))

    ## Male Survey LFreqs ##
    LFreqList[[paste("LF_Iay1M", k, sep = "_")]] <- My_LFreqs(
      data = Iay1_M, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Iay2M", k, sep = "_")]] <- My_LFreqs(
      data = Iay2_M, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Iay3M", k, sep = "_")]] <- My_LFreqs(
      data = Iay3_M, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Iay4M", k, sep = "_")]] <- My_LFreqs(
      data = Iay4_M, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Iay5M", k, sep = "_")]] <- My_LFreqs(
      data = Iay5_M, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Iay6M", k, sep = "_")]] <- My_LFreqs(
      data = Iay6_M, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )
    LFreqList[[paste("LF_Iay7M", k, sep = "_")]] <- My_LFreqs(
      data = Iay7_M, multC = 1, Lbins = seq(30, 260, by = 2), Lyrs = c(51:90),
      Linf = Linf_M, K = K_M, t0 = t0_F
    )


    save(LFreqList, file = paste("LFreqList_", k, ".RData", sep = ""))

    #### END LFREQS #####



    Iy1[k, ] <- apply(Iay1_F, 2, sum) + apply(Iay1_M, 2, sum)
    Iy2[k, ] <- apply(Iay2_F, 2, sum) + apply(Iay2_M, 2, sum)
    Iy3[k, ] <- apply(Iay3_F, 2, sum) + apply(Iay3_M, 2, sum)
    Iy4[k, ] <- apply(Iay4_F, 2, sum) + apply(Iay4_M, 2, sum)
    Iy5[k, ] <- apply(Iay5_F, 2, sum) + apply(Iay5_M, 2, sum)
    Iy6[k, ] <- apply(Iay6_F, 2, sum) + apply(Iay6_M, 2, sum)
    Iy7[k, ] <- apply(Iay7_F, 2, sum) + apply(Iay7_M, 2, sum)

    lines(1:years, Iy1[k, ], type = "l", col = "grey")
    lines(1:years, Iy2[k, ], type = "l", col = "blue")
    lines(1:years, Iy3[k, ], type = "l", col = "red")
    lines(1:years, Iy4[k, ], type = "l", col = "green")
    lines(1:years, Iy5[k, ], type = "l", col = "orange")
    lines(1:years, Iy6[k, ], type = "l", col = "mediumorchid2")
    lines(1:years, Iy7[k, ], type = "l", col = "cyan")
    Iy1_CV[k, ] <- apply(Iay1_CV, 2, mean, na.rm = T)
    Iy2_CV[k, ] <- apply(Iay2_CV, 2, mean, na.rm = T)
    Iy3_CV[k, ] <- apply(Iay3_CV, 2, mean, na.rm = T)
    Iy4_CV[k, ] <- apply(Iay4_CV, 2, mean, na.rm = T)
    Iy5_CV[k, ] <- apply(Iay5_CV, 2, mean, na.rm = T)
    Iy6_CV[k, ] <- apply(Iay6_CV, 2, mean, na.rm = T)
    Iy7_CV[k, ] <- apply(Iay7_CV, 2, mean, na.rm = T)

    Cy_1[k, ] <- apply(Cay1_F, 2, sum, na.rm = T) + apply(Cay1_M, 2, sum, na.rm = T)
    Cy_2[k, ] <- apply(Cay2_F, 2, sum, na.rm = T) + apply(Cay2_M, 2, sum, na.rm = T)
    Cy_3[k, ] <- apply(Cay3_F, 2, sum, na.rm = T) + apply(Cay3_M, 2, sum, na.rm = T)
    Cy_4[k, ] <- apply(Cay4_F, 2, sum, na.rm = T) + apply(Cay4_M, 2, sum, na.rm = T)
    Cy[k, ] <- Cy_1[k, ] + Cy_2[k, ] + Cy_3[k, ] + Cy_4[k, ]

    Ny_F[k, ] <- apply(Nay_F, 2, sum, na.rm = T)
    Ny_M[k, ] <- apply(Nay_M, 2, sum, na.rm = T)
    Ny[k, ] <- Ny_F[k, ] + Ny_M[k, ]
    # lines(1:years, Ny[k,]/100, type='l', lwd=2)


    # write
    AnnualNC <- list()
    assign(paste("NayF", k, sep = "_"), Nay_F)
    assign(paste("NayM", k, sep = "_"), Nay_M)
    AnnualNC[[paste("NayF", k, sep = "_")]] <- Nay_F
    AnnualNC[[paste("NayM", k, sep = "_")]] <- Nay_M

    assign(paste("Cay1F", k, sep = "_"), Cay1_F)
    assign(paste("Cay2F", k, sep = "_"), Cay2_F)
    assign(paste("Cay3F", k, sep = "_"), Cay3_F)
    assign(paste("Cay4F", k, sep = "_"), Cay4_F)
    assign(paste("Cay1M", k, sep = "_"), Cay1_M)
    assign(paste("Cay2M", k, sep = "_"), Cay2_M)
    assign(paste("Cay3M", k, sep = "_"), Cay3_M)
    assign(paste("Cay4M", k, sep = "_"), Cay4_M)
    AnnualNC[[paste("Cay1F", k, sep = "_")]] <- Cay1_F
    AnnualNC[[paste("Cay2F", k, sep = "_")]] <- Cay2_F
    AnnualNC[[paste("Cay3F", k, sep = "_")]] <- Cay3_F
    AnnualNC[[paste("Cay4F", k, sep = "_")]] <- Cay4_F
    AnnualNC[[paste("Cay1M", k, sep = "_")]] <- Cay1_M
    AnnualNC[[paste("Cay2M", k, sep = "_")]] <- Cay2_M
    AnnualNC[[paste("Cay3M", k, sep = "_")]] <- Cay3_M
    AnnualNC[[paste("Cay4M", k, sep = "_")]] <- Cay4_M

    save(AnnualNC, file = paste("AnnualNC_", k, ".RData", sep = ""))
  } # end k loop


  #

  ResultsList[["Iy1"]] <- Iy1
  ResultsList[["Iy2"]] <- Iy2
  ResultsList[["Iy3"]] <- Iy3
  ResultsList[["Iy4"]] <- Iy4
  ResultsList[["Iy5"]] <- Iy5
  ResultsList[["Iy6"]] <- Iy6
  ResultsList[["Iy7"]] <- Iy7
  ResultsList[["Iy1_CV"]] <- Iy1_CV
  ResultsList[["Iy2_CV"]] <- Iy2_CV
  ResultsList[["Iy3_CV"]] <- Iy3_CV
  ResultsList[["Iy4_CV"]] <- Iy4_CV
  ResultsList[["Iy5_CV"]] <- Iy5_CV
  ResultsList[["Iy6_CV"]] <- Iy6_CV
  ResultsList[["Iy7_CV"]] <- Iy7_CV

  ResultsList[["M0"]] <- M0
  ResultsList[["FSS"]] <- FSS
  ResultsList[["Npups"]] <- Npups

  ResultsList[["Cy"]] <- Cy
  ResultsList[["Cy_1"]] <- Cy_1
  ResultsList[["Cy_2"]] <- Cy_2
  ResultsList[["Cy_3"]] <- Cy_3
  ResultsList[["Cy_4"]] <- Cy_4

  ResultsList[["Ny"]] <- Ny

  save(ResultsList, file = "ResultsList.RData")





  ########## TEST DFA #############
  # 51 - 90;
  # library(MARSS)
  yrs <- 51:90
  biasAll <- matrix(nrow = Niters, ncol = length(yrs))
  bias <- vector(length = Niters)
  RMSE <- vector(length = Niters)
  biasAll.BT <- matrix(nrow = Niters, ncol = length(yrs))
  RMSE.BT <- vector(length = Niters)

  FitRatios <- matrix(nrow = Niters, ncol = 7)
  FactorLoadings <- vector()
  DFATrends <- vector()
  DFATrendsBT <- vector()
  DFATrendsSE <- vector()
  DFATrendsSEBT <- vector()
  upCI_DFATrends <- vector()
  lowCI_DFATrends <- vector()

  datz_SD <- matrix(nrow = Niters, ncol = 7)

  # i=1


  for (i in 1:Niters) {
    I1 <- Iy1[i, 51:90]
    I2 <- Iy2[i, 51:90]
    I3 <- Iy3[i, 51:90]
    I4 <- Iy4[i, 51:90]
    I5 <- Iy5[i, 51:90]
    I6 <- Iy6[i, 51:90]
    I7 <- Iy7[i, 51:90]
    # when missing data
    # I1 = c(Iy1[i,51:53], rep(NA, 5), Iy1[i,59:64], NA, NA, Iy1[i,67:90])
    # I2 = c(rep(NA, 7), Iy2[i,58:90])
    # I3 = c(rep(NA, 13), Iy3[i,64:90])
    # I4 = c(rep(NA, 16), Iy4[i,67:90])
    # I5 = c(rep(NA, 22), Iy5[i,73:90])
    # I6 = c(rep(NA, 15), Iy6[i,66:79], rep(NA, 11))
    # I7 = c(rep(NA, 17), Iy7[i,68], NA, NA, Iy7[i,71], NA, NA, Iy7[i,74], NA, NA, Iy7[i,77], NA, NA, Iy7[i,80],
    #        NA, NA, Iy7[i,83], NA, NA, Iy7[i,86], NA, NA, Iy7[i,89], NA)
    #
    assign(paste("dat", i, sep = ""), rbind(I1, I2, I3, I4, I5, I6, I7))
    # assign(paste("dat",i,sep=""), rbind(Iy1[i,51:90], Iy2[i,51:90], Iy3[i,51:90], Iy4[i,51:90], Iy5[i,51:90], Iy6[i,51:90], Iy7[i,51:90]))

    dat.a <- get(paste("dat", i, sep = ""))




    dat <- dat.a * c

    TT <- ncol(dat)
    N.ts <- nrow(dat)
    # Standardize data
    datL <- log(dat)
    y.bar <- apply(datL, 1, mean, na.rm = TRUE)
    dat.dm <- (datL - y.bar) / y.bar
    gsd <- sd(c(dat.dm), na.rm = T)
    dat.z <- dat.dm / gsd
    y.bar
    gsd
    datz_SD[i, ] <- apply(dat.z, 1, sd, na.rm = T)
    datz_SD[i, ]



    dat.z <- as.matrix(dat.z)
    rownames(dat.z) <- c("Survey1", "Survey2", "Survey3", "Survey4", "Survey5", "Survey6", "Survey7")




    ##### fit DFA ####
    cntl.list <- list(minit = 200, maxit = 500000, abstol = 0.02, allow.degen = FALSE, conv.test.slope.tol = 0.05)
    R <- diag(c(CV1, CV2, CV3, CV4, CV5, CV6, CV7), nrow = 7, ncol = 7)
    # dfa = MARSS(dat.z, model = list(m=1,R="diagonal and equal"), control=cntl.list, form="dfa", z.score=TRUE)
    dfa <- MARSS(dat.z, model = list(m = 1, R = R), control = cntl.list, form = "dfa", z.score = FALSE)


    pars <- MARSSparamCIs(dfa, nboot = 10000)

    # get the inverse of the rotation matrix
    # H.inv = varimax(coef(dfa, type="matrix")$Z)$rotmat # if more than 1 common trend
    Z.rot <- coef(dfa, type = "matrix")$Z # %*% H.inv # if more than 1 common trend
    trends.rot <- dfa$states
    ts.trends <- t(trends.rot)
    par.mat <- coef(dfa, type = "matrix")
    fit.b <- par.mat$Z %*% dfa$states # + matrix(par.mat$D, nrow=N.ts) %*% covar # if DFA included covariates

    assign(paste("Z.rot", i, sep = "."), Z.rot)
    assign(paste("Z.upCI", i, sep = "."), pars$par.upCI$Z)
    assign(paste("Z.lowCI", i, sep = "."), pars$par.lowCI$Z)
    assign(paste("trends.rot", i, sep = "."), trends.rot)


    # plot factor loadings
    survey <- rownames(dat.z)
    minZ <- 0.00
    ylims <- c(-1.1 * max(abs(Z.rot)), 1.1 * max(abs(Z.rot)))
    par(mfrow = c(nrow(trends.rot), 1), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
    for (h in 1:nrow(trends.rot)) {
      plot(c(1:N.ts)[abs(Z.rot[, h]) > minZ], as.vector(Z.rot[abs(Z.rot[, h]) > minZ, h]),
        type = "h", lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylims, xlim = c(0, N.ts + 1)
      )
      for (j in 1:N.ts) {
        if (Z.rot[j, h] > minZ) {
          text(j, -0.05, survey[j], srt = 90, adj = 1, cex = 0.9)
        }
        if (Z.rot[j, h] < -minZ) {
          text(j, 0.05, survey[j], srt = 90, adj = 0, cex = 0.9)
        }
        abline(h = 0, lwd = 1, col = "gray")
        abline(h = 0.2, col = "gray", lty = 2)
        abline(h = -0.2, col = "gray", lty = 2)
      } # end j loop
      mtext(paste("Factor loadings on trend", h, sep = " "), side = 3, line = .5)
    } # end h loop


    # PLOT WITH CI's

    index <- dfa$states[1, ]
    indexSE <- dfa$states.se[1, ]
    lowerCI <- dfa$states[1, ] - 1.96 * dfa$states.se[1, ]
    upperCI <- dfa$states[1, ] + 1.96 * dfa$states.se[1, ]

    assign(paste("index", i, sep = ""), index)
    assign(paste("indexSE", i, sep = ""), indexSE)
    assign(paste("lowerCI", i, sep = ""), lowerCI)
    assign(paste("upperCI", i, sep = ""), upperCI)



    par(mfrow = c(4, 2))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type = "n", bty = "L", xlab = "Year", ylab = "Abundance", ylim = c(min(lowerCI), max(upperCI))))
    x <- as.numeric(df$yrs)
    polygon(c(x, rev(x)), c(df$lowerCI, rev(df$upperCI)), col = "grey", border = NA)
    with(df, lines(x, index, type = "l", lwd = 2, pch = 16))
    title(paste("Iteration", i, sep = " "))
    abline(h = 0)


    N <- Ny[i, 51:90]
    Nscale <- (N - mean(N)) / sd(N)
    NL <- log(N)
    NLscale <- (NL - mean(NL)) / sd(NL)
    index.z <- (index - mean(index)) / sd(index)
    index.z <- ifelse(is.na(index.z), 0, index.z)
    index.z <- ifelse(abs(index.z) == Inf, 0, index.z)



    indexSEBT <- indexSE * gsd
    assign(paste("indexSEBT", i, sep = ""), indexSEBT)

    indexBT <- exp(index * gsd + ((indexSEBT^2) / 2))
    assign(paste("indexBT", i, sep = ""), indexBT)


    lines(yrs, rescale(NL, to = c(min(index), max(index))), type = "l", col = "red", lwd = 2)
    assign(paste("N", i, sep = ""), N)
    assign(paste("Nscale", i, sep = ""), Nscale)
    assign(paste("NL", i, sep = ""), NL)
    assign(paste("NLscale", i, sep = ""), NLscale)


    biasAll[i, ] <- index.z - NLscale
    bias[i] <- mean(index.z - NLscale)
    RMSE[i] <- sqrt(sum((index.z - NLscale)^2) / length(index.z))

    indexBT.z <- (indexBT - mean(indexBT)) / sd(indexBT)
    biasAll.BT[i, ] <- indexBT.z - Nscale
    RMSE.BT[i] <- sqrt(sum((indexBT.z - Nscale)^2) / length(indexBT.z))

    # Plot fitted values
    survey <- rownames(dat.z)
    # par(mfcol=c(ceiling(N.ts/2),2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
    for (n in 1:length(survey)) {
      plot(yrs, dat.z[n, ], xlab = "", ylab = "abundance index", bty = "L", xaxt = "n", ylim = c(-4, 4), pch = 16, col = "blue")
      axis(1, labels = TRUE)
      lines(yrs, fit.b[n, ], lwd = 2)
      lines(yrs, rescale(NL, to = c(min(fit.b[n, ]), max(fit.b[n, ]))), type = "l", col = "red", lwd = 2)
      title(paste("Sandbar", survey[n], sep = " "))
    }


    #### ASSESS MODEL FITS
    sumResids <- rowSums((dat.z - fit.b)^2, na.rm = TRUE)
    sumObserved <- rowSums(dat.z^2, na.rm = TRUE)
    FitRatio <- sumResids / sumObserved
    FitRatio
    FitRatios[i, ] <- FitRatio
    # mean(FitRatio)
    # indexBT
    # dat.a




    FactorLoadings <- rbind(FactorLoadings, c(get(paste("Z.rot", i, sep = ".")), get(paste("Z.lowCI", i, sep = ".")), get(paste("Z.upCI", i, sep = "."))))
    DFATrends <- rbind(DFATrends, index)
    DFATrendsBT <- rbind(DFATrendsBT, indexBT)
    DFATrendsSE <- rbind(DFATrendsSE, indexSE)
    DFATrendsSEBT <- rbind(DFATrendsSEBT, indexSEBT)
    upCI_DFATrends <- rbind(upCI_DFATrends, upperCI)
    lowCI_DFATrends <- rbind(lowCI_DFATrends, lowerCI)
  } # END LOOP

  colnames(FactorLoadings) <- c(
    "FL1", "FL2", "FL3", "FL4", "FL5", "FL6", "FL7",
    "lowCI_FL1", "lowCI_FL2", "lowCI_FL3", "lowCI_FL4", "lowCI_FL5", "lowCI_FL6", "lowCI_FL7",
    "upCI_FL1", "upCI_FL2", "upCI_FL3", "upCI_FL4", "upCI_FL5", "upCI_FL6", "upCI_FL7"
  )
  FR <- apply(FitRatios, 1, mean)
  FLoadings <- cbind(FactorLoadings, meanFactorLoading = FR)

  DFA_Results <- list()
  DFA_Results[["FitRatios"]] <- FitRatios
  DFA_Results[["DFATrends"]] <- DFATrends
  DFA_Results[["DFATrendsSE"]] <- DFATrendsSE
  DFA_Results[["DFATrendsBT"]] <- DFATrendsBT
  DFA_Results[["DFATrendsSEBT"]] <- DFATrendsSEBT
  DFA_Results[["upCI_DFATrends"]] <- upCI_DFATrends
  DFA_Results[["lowCI_DFATrends"]] <- lowCI_DFATrends
  DFA_Results[["FLoadings"]] <- FLoadings
  DFA_Results[["bias"]] <- bias
  DFA_Results[["biasAll"]] <- biasAll
  DFA_Results[["RMSE"]] <- RMSE
  DFA_Results[["biasAll.BT"]] <- biasAll.BT
  DFA_Results[["RMSE.BT"]] <- RMSE.BT
  DFA_Results[["GSD"]] <- gsd

  DFA_Results[["datz_SD"]] <- datz_SD

  save(DFA_Results, file = "DFA_Results.RData")




  ###################################  plotting ###################################

  if (plot == TRUE) {
    png(
      filename = "Trend.png",
      type = "cairo",
      units = "mm",
      width = 300,
      height = 300,
      pointsize = 12,
      res = 300
    )

    par(mfrow = c(10, 10))
    par(mar = c(0.3, 0.3, 0.1, 0.1), tcl = -0.1, mgp = c(2, 0.6, 0))
    for (l in 1:Niters) {
      index <- get(paste("index", l, sep = ""))
      lowerCI <- get(paste("lowerCI", l, sep = ""))
      upperCI <- get(paste("upperCI", l, sep = ""))
      df <- data.frame(yrs, index, lowerCI, upperCI)
      with(df, plot(index ~ yrs, type = "n", bty = "L", xlab = "", ylab = "", ylim = c(min(lowerCI), max(upperCI)), axes = FALSE))
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      x <- as.numeric(df$yrs)
      polygon(c(x, rev(x)), c(df$lowerCI, rev(df$upperCI)), col = "grey", border = NA)
      with(df, lines(x, index, type = "l", lwd = 2, pch = 16))
      # title(paste("Iteration",l,sep=" "))
      abline(h = 0)
      NLscale <- get(paste("NLscale", l, sep = ""))
      lines(yrs, NLscale, type = "l", col = "red", lwd = 2)
      text(x = min(yrs) + 4, y = 0.8 * max(min(lowerCI)), labels = format(FR[l], digits = 3), cex = 1)
    }

    dev.off()


    png(
      filename = "RawIndices.png",
      type = "cairo",
      units = "mm",
      width = 300,
      height = 300,
      pointsize = 12,
      res = 300
    )

    par(mfrow = c(10, 10))
    par(mar = c(0.3, 0.3, 0.1, 0.1), tcl = -0.1, mgp = c(2, 0.6, 0))
    for (l in 1:Niters) {
      plot(51:90, Iy1[l, 51:90],
        type = "l", col = "grey", xlim = c(51, 90), axes = F,
        ylim = c(0, max(c(Iy1[l, ], Iy2[l, ], Iy3[l, ], Iy4[l, ], Iy5[l, ], Iy6[l, ], Iy7[l, ]), na.rm = T) + 0.1)
      )
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      lines(51:90, Iy2[l, 51:90], type = "l", col = "blue")
      lines(51:90, Iy3[l, 51:90], type = "l", col = "red")
      lines(51:90, Iy4[l, 51:90], type = "l", col = "green")
      lines(51:90, Iy5[l, 51:90], type = "l", col = "orange")
      lines(51:90, Iy6[l, 51:90], type = "l", col = "mediumorchid2")
      lines(51:90, Iy7[l, 51:90], type = "l", col = "cyan")
      # title(paste("Iteration",l, sep=" "))
    }

    dev.off()


    png(
      filename = "FactorLoadings.png",
      type = "cairo",
      units = "mm",
      width = 300,
      height = 300,
      pointsize = 12,
      res = 300
    )

    # Before plotting factor loadings, we'll need to assign Z.rot, trends.rot
    par(mfrow = c(10, 10))
    par(mar = c(0.3, 0.3, 0.1, 0.1), tcl = -0.1, mgp = c(2, 0.6, 0))
    survey <- rownames(dat.z)
    minZ <- 0.00
    for (l in 1:Niters) {
      Z.rot <- get(paste("Z.rot", l, sep = "."))
      trends.rot <- get(paste("trends.rot", l, sep = "."))
      ylims <- c(-1.1 * max(abs(Z.rot)), 1.1 * max(abs(Z.rot)))
      # par(mfrow=c(nrow(trends.rot),1), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
      for (h in 1:nrow(trends.rot)) {
        plot(c(1:N.ts)[abs(Z.rot[, h]) > minZ], as.vector(Z.rot[abs(Z.rot[, h]) > minZ, h]),
          type = "h", lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylims, xlim = c(0, N.ts + 1), axes = FALSE
        )
        axis(1, labels = FALSE)
        axis(2, labels = FALSE)
        abline(h = 0, lwd = 1, col = "gray")
        abline(h = 0.2, col = "gray", lty = 2)
        abline(h = -0.2, col = "gray", lty = 2)
      } # end h loop
    }

    dev.off()



    png(
      filename = "Trend_scale.png",
      type = "cairo",
      units = "mm",
      width = 300,
      height = 300,
      pointsize = 12,
      res = 300
    )


    par(mfrow = c(10, 10))
    par(mar = c(0.3, 0.3, 0.1, 0.1), tcl = -0.1, mgp = c(2, 0.6, 0))
    for (l in 1:Niters) {
      index <- get(paste("index", l, sep = ""))
      lowerCI <- get(paste("lowerCI", l, sep = ""))
      upperCI <- get(paste("upperCI", l, sep = ""))
      df <- data.frame(yrs, index, lowerCI, upperCI)
      with(df, plot(index ~ yrs, type = "n", bty = "L", xlab = "", ylab = "", ylim = c(min(lowerCI), max(upperCI)), axes = FALSE))
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      x <- as.numeric(df$yrs)
      polygon(c(x, rev(x)), c(df$lowerCI, rev(df$upperCI)), col = "grey", border = NA)
      with(df, lines(x, index, type = "l", lwd = 2, pch = 16))
      # title(paste("Iteration",l,sep=" "))
      abline(h = 0)
      NL <- get(paste("NL", l, sep = ""))
      lines(yrs, rescale(NL, to = c(min(index), max(index))), type = "l", col = "red", lwd = 2)
      text(x = min(yrs) + 4, y = 0.8 * min(lowerCI), labels = format(FR[l], digits = 3), cex = 1)
    }


    dev.off()




    png(
      filename = "Trend_N.png",
      type = "cairo",
      units = "mm",
      width = 300,
      height = 300,
      pointsize = 12,
      res = 300
    )

    par(mfrow = c(10, 10))
    par(mar = c(0.3, 0.3, 0.1, 0.1), tcl = -0.1, mgp = c(2, 0.6, 0))
    for (l in 1:Niters) {
      index <- get(paste("index", l, sep = ""))
      NL <- get(paste("NL", l, sep = ""))
      plot(yrs, N, type = "l", col = "red", axes = FALSE)
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      lines(yrs, rescale(index, to = c(min(N), max(N))))
      text(x = min(yrs) + 4, y = 1.05 * min(N), labels = format(FR[l], digits = 3), cex = 1)
    }

    dev.off()





    png(
      filename = "BackTransformed_Trend_scale.png",
      type = "cairo",
      units = "mm",
      width = 300,
      height = 300,
      pointsize = 12,
      res = 300
    )

    par(mfrow = c(10, 10))
    par(mar = c(0.3, 0.3, 0.1, 0.1), tcl = -0.1, mgp = c(2, 0.6, 0))
    for (l in 1:Niters) {
      indexBT <- get(paste("indexBT", l, sep = ""))
      indexSEBT <- get(paste("indexSEBT", l, sep = ""))
      lowerCI <- indexBT - (1.96 * indexSEBT)
      upperCI <- indexBT + (1.96 * indexSEBT)
      df <- data.frame(yrs, indexBT, lowerCI, upperCI)
      with(df, plot(index ~ yrs, type = "n", bty = "L", xlab = "", ylab = "", ylim = c(min(lowerCI), max(upperCI)), axes = FALSE))
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      x <- as.numeric(df$yrs)
      polygon(c(x, rev(x)), c(df$lowerCI, rev(df$upperCI)), col = "grey", border = NA)
      with(df, lines(x, indexBT, type = "l", lwd = 2, pch = 16))
      # title(paste("Iteration",l,sep=" "))
      # abline(h=0)
      N <- get(paste("N", l, sep = ""))
      lines(yrs, rescale(N, to = c(min(indexBT), max(indexBT))), type = "l", col = "red", lwd = 2)
    }

    dev.off()



    png(
      filename = "Trend_Standardize.png",
      type = "cairo",
      units = "mm",
      width = 300,
      height = 300,
      pointsize = 12,
      res = 300
    )

    par(mfrow = c(10, 10))
    par(mar = c(0.3, 0.3, 0.1, 0.1), tcl = -0.1, mgp = c(2, 0.6, 0))
    for (l in 1:Niters) {
      index <- get(paste("index", l, sep = ""))
      NL <- get(paste("NL", l, sep = ""))
      plot(yrs, (NL - mean(NL)) / sd(NL), type = "l", col = "red", axes = FALSE, lwd = 1.5)
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      lines(yrs, (index - mean(index)) / sd(index), type = "l", lwd = 1.5)
      text(x = min(yrs) + 4, y = 0.8 * min((index - mean(index)) / sd(index)), labels = format(FR[l], digits = 3), cex = 1)
    }

    dev.off()
  } # end plot==T





  ##################################################### END #################################################

  return(list(RMSE = RMSE, biasAll = biasAll, FitRatios = FitRatios, RMSE.BT = RMSE.BT, biasAll.BT = biasAll.BT))
} # END FUNCTION



#### Example Application ####
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial27")
SB27 = DFA_Sim(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q0,
               q_survey4 = q0,
               q_survey5 = q10,
               q_survey6 = q0,
               q_survey7 = q3,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c = c(12, 3, 1, 0.4, 200, 0.8, 5))

save(SB27, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial27\\SB27.RData")


### Read in available life history / fishery info: ####

## Life History Information
Linf_M <- 175.5
K_M <- 0.143
t0_M <- -2.388
Linf_F <- 183.3
K_F <- 0.124
t0_F <- -3.098
# vonBert = rbvn(nrow(df.expanded),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
ages <- 1:31


# VON BERT:  Lt = Linf *(1-exp(-K*(t-t0)))
# SCHNUTE: L1 + (L2 - L1) * ((1-exp(-K*(t-t1)))/(1-exp(-K*(t2-t1))))
(L0_F <- Linf_F * (1 - exp(-K_F * (0 - t0_F)))) # L0 = 58.46759
(L1_F <- Linf_F * (1 - exp(-K_F * (1 - t0_F)))) # L1 = 73.02556
(L0_M <- Linf_M * (1 - exp(-K_M * (0 - t0_M)))) # L0 = 50.76955
(L1_M <- Linf_M * (1 - exp(-K_M * (1 - t0_M)))) # L1 = 67.38937

## Selectivity read in from .csv file
# Sel = read.csv("D:\\vspace1\\DFA_Simulation\\SB\\Selectivity_new.csv")
# FullSel = Sel
# Sel = Sel[-1,]

## equilibrium numbers at age
Equilibrium_Nay_F <- c(
  500.00065, 425.90205, 362.78465, 309.02105, 263.22506, 224.21589, 191.48498, 170.37628, 151.59453,
  134.88322, 120.01412, 106.78413, 95.01259, 84.53869, 75.21941, 66.92746, 59.54959, 52.98503,
  47.14412, 41.94710, 37.32298, 33.20861, 29.54780, 26.29054, 23.39236, 20.81366, 18.51922,
  16.47772, 14.66127, 13.04506, 11.60701
)

Equilibrium_Nay_M <- c(
  500.00065, 425.90205, 362.78465, 309.02105, 263.22506, 224.21589, 191.48498, 170.37628, 151.59453,
  134.88322, 120.01412, 106.78413, 95.01259, 84.53869, 75.21941, 66.92746, 59.54959, 52.98503,
  47.14412, 41.94710, 37.32298, 33.20861, 29.54780, 26.29054, 23.39236, 20.81366, 18.51922,
  16.47772, 14.66127, 13.04506, 11.60701
)

ages <- 1:31
A <- max(ages)
years <- 100

# Define maturity and fecundity
mat_F <- as.vector(c(0, 0, 0, 0, 0, 0.01, 0.02, 0.03, 0.06, 0.12, 0.21, 0.33, 0.49, 0.65, 0.78, 0.88, 0.93, 0.96, 0.98, 0.99, 0.99, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
fec <- c(
  6.61072819, 7.027399402, 7.39547835, 7.720631873, 8.00786594, 8.261602725, 8.485748685, 8.683754708, 8.858669237,
  9.013185206, 9.149681498, 9.27025957, 9.376775809, 9.470870107, 9.553991113, 9.627418534, 9.692282837, 9.749582655, 9.80020016,
  9.844914642, 9.884414515, 9.919307906, 9.950132025, 9.97736143, 10.00141534, 10.02266407, 10.04143478, 10.05801644, 10.07266435,
  10.08560401, 10.09703465
)

M_const <- c(
  0.1604, 0.1604, 0.1604, 0.1604, 0.1604, 0.1578, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168,
  0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168, 0.1168,
  0.1168, 0.1168, 0.1168, 0.1168
)


F1 <- as.vector(c(rep(0, 45), rep(0.1, 10), rep(0.3, 10), rep(0.2, 10), rep(0.05, 25))) # F1
load(file = "D:\\vspace1\\DFA_Simulation\\SB\\Iay_eq_F1.RData")

# Define q timeseries

q0 <- rep(0.01, 100) # q0 = constant
q1 <- rep(c(0.01, 0.02), 50) # q1 = period 2
q2 <- rep(c(0.01, 0.02, 0.015), length = 100) # q2 = period 3
q3 <- c(rep(0.01, 50), rep(seq(0.01, 0.045, by = 0.005), each = 5), rep(0.045, 10)) # q3 = step increase (step=5)
q4 <- c(rep(0.05, 50), rep(seq(0.05, 0.0325, by = -0.0025), each = 5), rep(0.0325, 10)) # q4 = step decrease (step=5)
q5 <- c(rep(seq(0.03, 0.05, by = 0.005), seq(0.05, 0.03, by = -0.005), length = 100)) # q5 = increase then decrease
q6 <- c(rep(0.01, 50), seq(0.01, 0.049, by = 0.001), rep(0.049, 10)) # q6 = increase
q7 <- c(rep(0.049, 50), seq(0.049, 0.01, by = -0.001), rep(0.01, 10)) # q7 = decrease
set.seed(430)
q8 <- runif(100, 0.03, 0.05) # q8 = random between 0.01 and 0.05
set.seed(430)
q9 <- jitter(q6, amount = 0.005) # q9 = jitter increase
set.seed(431)
q10 <- jitter(q7, amount = 0.005) # q10 = jitter decrease
set.seed(432)
q11 <- jitter(q1, amount = 0.005) # q11 = jitter period 2
set.seed(433)
q12 <- jitter(q2, amount = 0.005) # q12 = jitter period 3                                                         # jitter decrease
set.seed(434)
q13 <- jitter(q0, amount = 0.005) # q13 = jitter constant
set.seed(435)
q14 <- jitter(q0, amount = 0.003) # q14 = jitter constant 2
set.seed(436)
q15 <- jitter(q0, amount = 0.002) # q15 = jitter constant 3
set.seed(437)
q16 <- jitter(q0, amount = 0.004) # q16 = jitter constant 3
set.seed(438)
q17 <- jitter(q0, amount = 0.003) # q17 = jitter constant 4
set.seed(439)
q18 <- jitter(q0, amount = 0.002) # q18 = jitter constant 5
set.seed(440)
q19 <- jitter(q0, amount = 0.004) # q19 = jitter constant 6
set.seed(441)
q20 <- jitter(q0, amount = 0.003) # q20 = jitter constant 7
set.seed(441)
q21 <- jitter(q0, amount = 0.002) # q21 = jitter constant 8
