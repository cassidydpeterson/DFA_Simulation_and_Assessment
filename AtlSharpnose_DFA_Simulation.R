#####################################
# Atlantic sharpnose shark
# DFA Simulation
# C PETERSON
#####################################

### load packages
library(scales)
library(MARSS)
library(vioplot)

rbvn <- function(n, m1, s1, m2, s2, rho) {
  X1 <- rnorm(n, m1, s1)
  X2 <- rnorm(n, m2 + (s2 / s1) * rho *
    (X1 - m1), sqrt((1 - rho^2) * s2^2))
  data.frame(Linf = X1, K = X2)
}


######### SIMULATION FUNCTIONs ###########
# 3 survey indices
DFA_Sim <- function(F_var, CV1, CV2, CV3, q_survey1, q_survey2, q_survey3, Iay_eq = Iay_eq, c = 5, plot = T) {
  # Iay
  par(mfrow = c(1, 1))
  plot(1:years, apply(Iay_eq, 2, sum), type = "l", ylim = c(0, max(apply(Iay_eq, 2, sum), na.rm = T) + max(apply(Iay_eq, 2, sum), na.rm = T) * 0.1))
  abline(v = 40, lwd = 2, col = "red")
  abline(v = 65, lwd = 2, col = "red")

  F_var <- F_var
  CV1 <- CV1
  CV2 <- CV2
  CV3 <- CV3
  q_survey1 <- q_survey1
  q_survey2 <- q_survey2
  q_survey3 <- q_survey3
  c <- c
  c1 <- c

  ####### ADD STOCHASTICITY ###########
  Niters <- 100
  Iay1 <- matrix(nrow = A, ncol = years)
  Iay2 <- matrix(nrow = A, ncol = years)
  Iay3 <- matrix(nrow = A, ncol = years)
  Iay1_CV <- matrix(nrow = A, ncol = years)
  Iay2_CV <- matrix(nrow = A, ncol = years)
  Iay3_CV <- matrix(nrow = A, ncol = years)
  Cy <- matrix(nrow = Niters, ncol = years)
  Ny <- matrix(nrow = Niters, ncol = years)
  Iy1 <- matrix(nrow = Niters, ncol = years)
  Iy2 <- matrix(nrow = Niters, ncol = years)
  Iy3 <- matrix(nrow = Niters, ncol = years)
  Iy1_CV <- matrix(nrow = Niters, ncol = years)
  Iy2_CV <- matrix(nrow = Niters, ncol = years)
  Iy3_CV <- matrix(nrow = Niters, ncol = years)
  M0 <- matrix(nrow = Niters, ncol = years)
  Npups <- matrix(nrow = Niters, ncol = years)
  FSS <- matrix(nrow = Niters, ncol = years) # Female spawning stock
  Z0 <- 1.904416
  Zmin <- 0.209
  Npups_eq <- 3320.671
  beta <- 0.5807613
  ResultsList <- list()

  set.seed(430)


  for (k in 1:Niters) {
    # paste("Iteration",k,sep=" ")
    # Create matrices to store results
    Nay <- matrix(nrow = A, ncol = years)
    Cay <- matrix(nrow = A, ncol = years)

    # Inputs:
    M <- M_const
    FM <- F_var

    # Set up equilibrium conditions
    Nay[, 1] <- Equilibrium_Nay

    for (i in 1:(years - 1)) {
      matjit <- rnorm(length(mat), mat, mat * 0.01)
      matjit <- ifelse(matjit > 1, 1, matjit)
      FSS[k, i] <- sum(Nay[, i] * matjit) # Fecund Stock size
      Npups[k, i] <- sum(Nay[, i] * matjit * rnorm(length(fec), fec, fec * 0.1)) # equilibrium Npups = 3320.671
      M0[k, i] <- rnorm(1, (((1 - (Npups[k, i] / Npups_eq)^beta) * (Zmin - Z0)) + Z0), 0.1) # equilibrium M0 = 1.9044160
      Nay[1, i + 1] <- Npups[k, i] * exp(-M0[k, i])

      Mjit <- rnorm(length(M), M, M * 0.1)

      Zay <- vector(length = length(M))
      for (j in 1:A) {
        Zay[j] <- Mjit[j] + (Sel_catch[j, 2] * FM[i] / 4) + (Sel_catch[j, 3] * FM[i] / 4) + (Sel_catch[j, 4] * FM[i] / 4) + (Sel_catch[j, 5] * FM[i] / 4)
      } # end j loop Z


      for (j in 1:(A - 2)) {
        Nay[j + 1, i + 1] <- Nay[j, i] * exp(-Zay[j])
        Cay[j, i] <- Nay[j, i] * ((Sel_catch[j, 2] * FM[i] / 4) / Zay[j]) * (1 - exp(-Zay[j])) +
          Nay[j, i] * ((Sel_catch[j, 3] * FM[i] / 4) / Zay[j]) * (1 - exp(-Zay[j])) +
          Nay[j, i] * ((Sel_catch[j, 4] * FM[i] / 4) / Zay[j]) * (1 - exp(-Zay[j])) +
          Nay[j, i] * ((Sel_catch[j, 5] * FM[i] / 4) / Zay[j]) * (1 - exp(-Zay[j]))
      } # END J LOOP

      Nay[A, i + 1] <- Nay[A - 1, i] * exp(-Zay[A - 1]) +
        Nay[A, i] * exp(-Zay[A])

      Cay[A - 1, i] <- Nay[A - 1, i] * ((Sel_catch[A - 1, 2] * FM[i] / 4) / Zay[A - 1]) * (1 - exp(-Zay[A - 1])) +
        Nay[A - 1, i] * ((Sel_catch[A - 1, 3] * FM[i] / 4) / Zay[A - 1]) * (1 - exp(-Zay[A - 1])) +
        Nay[A - 1, i] * ((Sel_catch[A - 1, 4] * FM[i] / 4) / Zay[A - 1]) * (1 - exp(-Zay[A - 1])) +
        Nay[A - 1, i] * ((Sel_catch[A - 1, 5] * FM[i] / 4) / Zay[A - 1]) * (1 - exp(-Zay[A - 1]))

      Cay[A, i] <- Nay[A, i] * ((Sel_catch[A, 2] * FM[i] / 4) / Zay[A]) * (1 - exp(-Zay[A])) +
        Nay[A, i] * ((Sel_catch[A, 3] * FM[i] / 4) / Zay[A]) * (1 - exp(-Zay[A])) +
        Nay[A, i] * ((Sel_catch[A, 4] * FM[i] / 4) / Zay[A]) * (1 - exp(-Zay[A])) +
        Nay[A, i] * ((Sel_catch[A, 5] * FM[i] / 4) / Zay[A]) * (1 - exp(-Zay[A]))



      for (h in 1:A) {
        Iay1_CV[h, i] <- runif(1, min = CV1 - 0.1, max = CV1 + 0.1)
        Iay1[h, i] <- max(q_survey1[i] * Sel_survey[h, 2] * Nay[h, i] *
          exp(rnorm(1, 0, (q_survey1[i] * Sel_survey[h, 2] * Nay[h, i]) * Iay1_CV[h, i]) -
            (((q_survey1[i] * Sel_survey[h, 2] * Nay[h, i] * Iay1_CV[h, i])^2) / 2)), 0) # uncertainty based off CV;
        Iay2_CV[h, i] <- runif(1, min = CV2 - 0.1, max = CV2 + 0.1)
        Iay2[h, i] <- max(q_survey2[i] * Sel_survey[h, 3] * Nay[h, i] *
          exp(rnorm(1, 0, (q_survey2[i] * Sel_survey[h, 3] * Nay[h, i]) * Iay2_CV[h, i]) -
            (((q_survey2[i] * Sel_survey[h, 3] * Nay[h, i] * Iay2_CV[h, i])^2) / 2)), 0)


        Iay3_CV[h, i] <- runif(1, min = CV3 - 0.1, max = CV3 + 0.1)
        Iay3[h, i] <- max(q_survey3[i] * Sel_survey[h, 2] * Nay[h, i] *
          exp(rnorm(1, 0, (q_survey3[i] * Sel_survey[h, 2] * Nay[h, i]) * Iay3_CV[h, i]) -
            (((q_survey3[i] * Sel_survey[h, 2] * Nay[h, i] * Iay3_CV[h, i])^2) / 2)), 0)
      } # end h
    } # ends i loop

    # USE VON BERT CURVE TO ESTIMATE LENGTH COMP.
    # Lt = Linf * (1 - exp(-K*(t - t0)))
    # For sexes combined.
    Linf <- 81.6
    K <- 0.49
    t0 <- -0.97

    # L0 = Linf * (1 - exp(-K*(0-t0)))

    # GET LENGTH FREQUENCIES
    Cay_N <- cbind(ages, round(Cay))
    assign(paste("LF_Cay", k, sep = "_"), data.frame(LBins = seq(20, 122, by = 2)))
    for (f in 40:65) {
      df <- as.data.frame(Cay_N[, c(1, f + 1)])
      colnames(df) <- c("ages", "count")
      df.expanded <- df[rep(row.names(df), df$count), 1:2]
      vonBert <- rbvn(nrow(df.expanded), m1 = Linf, s1 = Linf * 0.1, m2 = K, s2 = K * 0.1, rho = -0.9)
      df.expanded$L <- vonBert$Linf * (1 - exp(-(vonBert$K) * (jitter(df.expanded$ages, amount = 0.5) - t0)))
      hst <- hist(df.expanded$L, breaks = c(seq(20, 124, by = 2)), plot = F)
      assign(paste("LF_Cay", k, sep = "_"), cbind(get(paste("LF_Cay", k, sep = "_")), hst$counts))
    } # end f loop

    Iay1_N <- cbind(ages, round(Iay1 * 100))
    assign(paste("LF_Iay1", k, sep = "_"), data.frame(LBins = seq(20, 122, by = 2)))
    for (l in 40:65) {
      df <- as.data.frame(Iay1_N[, c(1, l + 1)])
      colnames(df) <- c("ages", "count")
      df.expanded <- df[rep(row.names(df), df$count), 1:2]
      vonBert <- rbvn(nrow(df.expanded), m1 = Linf, s1 = Linf * 0.1, m2 = K, s2 = K * 0.1, rho = -0.9)
      df.expanded$L <- vonBert$Linf * (1 - exp(-(vonBert$K) * (jitter(df.expanded$ages, amount = 0.5) - t0)))
      hst <- hist(df.expanded$L, breaks = c(seq(20, 124, by = 2)), plot = F)
      assign(paste("LF_Iay1", k, sep = "_"), cbind(get(paste("LF_Iay1", k, sep = "_")), hst$counts))
    } # end l loop

    Iay2_N <- cbind(ages, round(Iay2 * 100))
    assign(paste("LF_Iay2", k, sep = "_"), data.frame(LBins = seq(20, 122, by = 2)))
    for (m in 40:65) {
      df <- as.data.frame(Iay2_N[, c(1, m + 1)])
      colnames(df) <- c("ages", "count")
      df.expanded <- df[rep(row.names(df), df$count), 1:2]
      vonBert <- rbvn(nrow(df.expanded), m1 = Linf, s1 = Linf * 0.1, m2 = K, s2 = K * 0.1, rho = -0.9)
      df.expanded$L <- vonBert$Linf * (1 - exp(-(vonBert$K) * (jitter(df.expanded$ages, amount = 0.5) - t0)))
      hst <- hist(df.expanded$L, breaks = c(seq(20, 124, by = 2)), plot = F)
      assign(paste("LF_Iay2", k, sep = "_"), cbind(get(paste("LF_Iay2", k, sep = "_")), hst$counts))
    } # end m loop


    Iay3_N <- cbind(ages, round(Iay3 * 100))
    assign(paste("LF_Iay3", k, sep = "_"), data.frame(LBins = seq(20, 122, by = 2)))
    for (n in 40:65) {
      df <- as.data.frame(Iay3_N[, c(1, n + 1)])
      colnames(df) <- c("ages", "count")
      df.expanded <- df[rep(row.names(df), df$count), 1:2]
      vonBert <- rbvn(nrow(df.expanded), m1 = Linf, s1 = Linf * 0.1, m2 = K, s2 = K * 0.1, rho = -0.9)
      df.expanded$L <- vonBert$Linf * (1 - exp(-(vonBert$K) * (jitter(df.expanded$ages, amount = 0.5) - t0)))
      hst <- hist(df.expanded$L, breaks = c(seq(20, 124, by = 2)), plot = F)
      assign(paste("LF_Iay3", k, sep = "_"), cbind(get(paste("LF_Iay3", k, sep = "_")), hst$counts))
    } # end n loop

    # save length-frequencies
    LFreqList <- list()
    LFreqList[[paste("LF_Cay", k, sep = "_")]] <- get(paste("LF_Cay", k, sep = "_"))
    LFreqList[[paste("LF_Iay1", k, sep = "_")]] <- get(paste("LF_Iay1", k, sep = "_"))
    LFreqList[[paste("LF_Iay2", k, sep = "_")]] <- get(paste("LF_Iay2", k, sep = "_"))
    LFreqList[[paste("LF_Iay3", k, sep = "_")]] <- get(paste("LF_Iay3", k, sep = "_"))

    save(LFreqList, file = paste("LFreqList_", k, ".RData", sep = ""))




    Iy1[k, ] <- apply(Iay1, 2, sum)
    Iy2[k, ] <- apply(Iay2, 2, sum)
    Iy3[k, ] <- apply(Iay3, 2, sum)
    lines(1:years, Iy1[k, ], type = "l", col = "grey")
    lines(1:years, Iy2[k, ], type = "l", col = "blue")
    lines(1:years, Iy3[k, ], type = "l", col = "red")
    Iy1_CV[k, ] <- apply(Iay1_CV, 2, mean, na.rm = T)
    Iy2_CV[k, ] <- apply(Iay2_CV, 2, mean, na.rm = T)
    Iy3_CV[k, ] <- apply(Iay3_CV, 2, mean, na.rm = T)

    Cy[k, ] <- apply(Cay, 2, sum)
    Ny[k, ] <- apply(Nay, 2, sum)



    # Save
    assign(paste("Nay", k, sep = "_"), Nay)
    assign(paste("Cay", k, sep = "_"), Cay)

    # Save
    AnnualNC <- list()
    AnnualNC[[paste("Nay", k, sep = "_")]] <- Nay
    AnnualNC[[paste("Cay", k, sep = "_")]] <- Cay
    AnnualNC[[paste("Iay1", k, sep = "_")]] <- Iay1
    AnnualNC[[paste("Iay2", k, sep = "_")]] <- Iay2
    AnnualNC[[paste("Iay3", k, sep = "_")]] <- Iay3

    save(AnnualNC, file = paste("AnnualNC_", k, ".RData", sep = ""))
  } # end k loop



  ResultsList[["Iy1"]] <- Iy1
  ResultsList[["Iy2"]] <- Iy2
  ResultsList[["Iy3"]] <- Iy3

  ResultsList[["Iy1_CV"]] <- Iy1_CV
  ResultsList[["Iy2_CV"]] <- Iy2_CV
  ResultsList[["Iy3_CV"]] <- Iy3_CV


  ResultsList[["M0"]] <- M0
  ResultsList[["FSS"]] <- FSS
  ResultsList[["Npups"]] <- Npups

  ResultsList[["Cy"]] <- Cy
  ResultsList[["Ny"]] <- Ny

  save(ResultsList, file = "ResultsList.RData")



  ########## TEST DFA #############

  # library(MARSS)
  yrs <- 40:(years - 35)
  biasAll <- matrix(nrow = Niters, ncol = length(yrs))
  bias <- vector(length = Niters)
  RMSE <- vector(length = Niters)
  biasAll.BT <- matrix(nrow = Niters, ncol = length(yrs))
  RMSE.BT <- vector(length = Niters)

  FitRatios <- matrix(nrow = Niters, ncol = 3)
  FactorLoadings <- vector()
  DFATrends <- vector()
  DFATrendsBT <- vector()
  DFATrendsSE <- vector()
  DFATrendsSEBT <- vector()
  upCI_DFATrends <- vector()
  lowCI_DFATrends <- vector()

  datz_SD <- matrix(nrow = Niters, ncol = 3)
  # i=1


  for (i in 1:Niters) {
    assign(paste("dat", i, sep = ""), rbind(Iy2[i, 40:(years - 35)], Iy1[i, 40:(years - 35)], Iy3[i, 40:(years - 35)]))


    dat.a <- get(paste("dat", i, sep = ""))

    dat <- dat.a * c

    TT <- ncol(dat)
    N.ts <- nrow(dat)
    # Standardize data
    datL <- log(dat)
    y.bar <- apply(datL, 1, mean, na.rm = TRUE)

    # y.bar = mean(datL)
    dat.dm <- (datL - y.bar) / y.bar
    gsd <- sd(c(dat.dm))
    dat.z <- dat.dm / gsd
    # y.bar
    # gsd
    # apply(dat.z, 1, sd)

    datz_SD[i, ] <- apply(dat.z, 1, sd, na.rm = T)

    # dat.dm = (datL - y.bar)
    dat.z <- as.matrix(dat.z)
    rownames(dat.z) <- c("Survey2", "Survey1", "Survey3")



    ######## SN DFA ########
    cntl.list <- list(minit = 200, maxit = 500000, abstol = 0.02, allow.degen = FALSE, conv.test.slope.tol = 0.05)
    R <- diag(c(CV2, CV1, CV3), nrow = 3, ncol = 3)

    dfa <- MARSS(dat.z, model = list(m = 1, R = R), control = cntl.list, form = "dfa", z.score = FALSE)


    pars <- MARSSparamCIs(dfa, nboot = 10000)

    # get the inverse of the rotation matrix
    # H.inv = varimax(coef(dfa, type="matrix")$Z)$rotmat # if n common trends > 1
    Z.rot <- coef(dfa, type = "matrix")$Z # %*% H.inv
    trends.rot <- dfa$states
    ts.trends <- t(trends.rot)
    par.mat <- coef(dfa, type = "matrix")
    fit.b <- par.mat$Z %*% dfa$states # + matrix(par.mat$D, nrow=N.ts) %*% covar # if covariates included

    assign(paste("Z.rot", i, sep = "."), Z.rot)
    assign(paste("Z.upCI", i, sep = "."), pars$par.upCI$Z)
    assign(paste("Z.lowCI", i, sep = "."), pars$par.lowCI$Z)
    assign(paste("trends.rot", i, sep = "."), trends.rot)


    # # plot factor loadings
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
    yrs <- 40:(years - 35)

    index <- dfa$states[1, ]
    indexSE <- dfa$states.se[1, ]
    lowerCI <- dfa$states[1, ] - 1.96 * dfa$states.se[1, ]
    upperCI <- dfa$states[1, ] + 1.96 * dfa$states.se[1, ]
    assign(paste("index", i, sep = ""), index)
    assign(paste("indexSE", i, sep = ""), indexSE)
    assign(paste("lowerCI", i, sep = ""), lowerCI)
    assign(paste("upperCI", i, sep = ""), upperCI)



    par(mfrow = c(2, 2))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type = "n", bty = "L", xlab = "Year", ylab = "Abundance", ylim = c(min(lowerCI), max(upperCI))))
    x <- as.numeric(df$yrs)
    polygon(c(x, rev(x)), c(df$lowerCI, rev(df$upperCI)), col = "grey", border = NA)
    with(df, lines(x, index, type = "l", lwd = 2, pch = 16))
    title(paste("Iteration", i, sep = " "))
    abline(h = 0)

    N <- Ny[i, 40:65]
    NL <- log(N)
    Nscale <- (N - mean(N)) / sd(N)
    NLscale <- (N - mean(N)) / sd(N)
    index.z <- (index - mean(index)) / sd(index)
    index.z <- ifelse(is.na(index.z), 0, index.z)
    index.z <- ifelse(abs(index.z) == Inf, 0, index.z)

    indexSEBT <- indexSE * gsd
    indexBT <- exp(index * gsd + ((indexSEBT^2) / 2))
    assign(paste("indexBT", i, sep = ""), indexBT)
    assign(paste("indexSEBT", i, sep = ""), indexSEBT)


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
      title(paste("Atlantic sharpnose", survey[n], sep = " "))
    }


    #### ASSESS MODEL FITS
    sumResids <- rowSums((dat.z - fit.b)^2, na.rm = TRUE)
    sumObserved <- rowSums(dat.z^2, na.rm = TRUE)
    FitRatio <- sumResids / sumObserved
    FitRatio
    # mean(FitRatio)
    # indexBT
    # dat.a

    FitRatios[i, ] <- FitRatio



    FactorLoadings <- rbind(FactorLoadings, c(get(paste("Z.rot", i, sep = ".")), get(paste("Z.lowCI", i, sep = ".")), get(paste("Z.upCI", i, sep = "."))))
    DFATrends <- rbind(DFATrends, index)
    DFATrendsBT <- rbind(DFATrendsBT, indexBT)
    DFATrendsSE <- rbind(DFATrendsSE, indexSE)
    DFATrendsSEBT <- rbind(DFATrendsSEBT, indexSEBT)
    upCI_DFATrends <- rbind(upCI_DFATrends, upperCI)
    lowCI_DFATrends <- rbind(lowCI_DFATrends, lowerCI)
  } # END LOOP



  colnames(FactorLoadings) <- c(
    "FL1", "FL2", "FL3",
    "lowCI_FL1", "lowCI_FL2", "lowCI_FL3",
    "upCI_FL1", "upCI_FL2", "upCI_FL3"
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





  ############  plotting ###################
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
      plot(1:years, Iy1[l, ],
        type = "l", col = "black", xlim = c(min(yrs), max(yrs)), axes = F,
        ylim = c(0, max(c(Iy1[l, yrs], Iy2[l, yrs], Iy3[l, yrs])) + 0.1)
      )
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      lines(1:years, Iy2[l, ], type = "l", col = "blue")
      lines(1:years, Iy3[l, ], type = "l", col = "red")
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
      plot(yrs, NL, type = "l", col = "red", axes = FALSE)
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      lines(yrs, rescale(index, to = c(min(NL), max(NL))))
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
    }

    dev.off()
  } # end if plot==T





  ############ END #####################

  return(list(RMSE = RMSE, biasAll = biasAll, FitRatios = FitRatios, RMSE.BT = RMSE.BT, biasAll.BT = biasAll.BT))
} # END FUNCTION



# 4 survey indices
DFA_Sim4 <- function(F_var, CV1, CV2, CV3, CV4, q_survey1, q_survey2, q_survey3, q_survey4, Iay_eq = Iay_eq, c = 5, plot = T) {
  # Iay
  par(mfrow = c(1, 1))
  plot(1:years, apply(Iay_eq, 2, sum), type = "l", ylim = c(0, max(apply(Iay_eq, 2, sum), na.rm = T) + max(apply(Iay_eq, 2, sum), na.rm = T) * 0.1))
  abline(v = 40, lwd = 2, col = "red")
  abline(v = 65, lwd = 2, col = "red")
  # plot(1:years, apply(Iay, 2, sum), type='l')
  # Iy = apply(Iay, 2, sum)

  F_var <- F_var
  CV1 <- CV1
  CV2 <- CV2
  CV3 <- CV3
  q_survey1 <- q_survey1
  q_survey2 <- q_survey2
  q_survey3 <- q_survey3
  c <- c
  c1 <- c

  ####### ADD STOCHASTICITY ###########
  Niters <- 100
  Iay1 <- matrix(nrow = A, ncol = years)
  Iay2 <- matrix(nrow = A, ncol = years)
  Iay3 <- matrix(nrow = A, ncol = years)
  Iay4 <- matrix(nrow = A, ncol = years)
  Iay1_CV <- matrix(nrow = A, ncol = years)
  Iay2_CV <- matrix(nrow = A, ncol = years)
  Iay3_CV <- matrix(nrow = A, ncol = years)
  Iay4_CV <- matrix(nrow = A, ncol = years)
  Cy <- matrix(nrow = Niters, ncol = years)
  Ny <- matrix(nrow = Niters, ncol = years)
  Iy1 <- matrix(nrow = Niters, ncol = years)
  Iy2 <- matrix(nrow = Niters, ncol = years)
  Iy3 <- matrix(nrow = Niters, ncol = years)
  Iy4 <- matrix(nrow = Niters, ncol = years)
  Iy1_CV <- matrix(nrow = Niters, ncol = years)
  Iy2_CV <- matrix(nrow = Niters, ncol = years)
  Iy3_CV <- matrix(nrow = Niters, ncol = years)
  Iy4_CV <- matrix(nrow = Niters, ncol = years)
  M0 <- matrix(nrow = Niters, ncol = years)
  Npups <- matrix(nrow = Niters, ncol = years)
  FSS <- matrix(nrow = Niters, ncol = years) # Female spawning stock
  Z0 <- 1.904416
  Zmin <- 0.209
  Npups_eq <- 3320.671
  beta <- 0.5807613
  ResultsList <- list()


  set.seed(430)

  # F_const = rep(0, yrs)

  for (k in 1:Niters) {
    # Create matrices to store results
    Nay <- matrix(nrow = A, ncol = years)
    Cay <- matrix(nrow = A, ncol = years)
    # M0 = vector(length=years)

    # Inputs:
    M <- M_const
    FM <- F_var

    # Set up equilibrium conditions
    Nay[, 1] <- Equilibrium_Nay

    for (i in 1:(years - 1)) {
      matjit <- rnorm(length(mat), mat, mat * 0.01)
      matjit <- ifelse(matjit > 1, 1, matjit)
      FSS[k, i] <- sum(Nay[, i] * matjit) # Fecund Stock size
      Npups[k, i] <- sum(Nay[, i] * matjit * rnorm(length(fec), fec, fec * 0.1)) # equilibrium Npups = 3320.671
      M0[k, i] <- rnorm(1, (((1 - (Npups[k, i] / Npups_eq)^beta) * (Zmin - Z0)) + Z0), 0.1) # Equilibrium M0 = 1.9044160
      Nay[1, i + 1] <- Npups[k, i] * exp(-M0[k, i])

      Mjit <- rnorm(length(M), M, M * 0.1)

      Zay <- vector(length = length(M))
      for (j in 1:A) {
        Zay[j] <- Mjit[j] + (Sel_catch[j, 2] * FM[i] / 4) + (Sel_catch[j, 3] * FM[i] / 4) + (Sel_catch[j, 4] * FM[i] / 4) + (Sel_catch[j, 5] * FM[i] / 4)
      } # end j loop Z


      for (j in 1:(A - 2)) {
        Nay[j + 1, i + 1] <- Nay[j, i] * exp(-Zay[j])
        Cay[j, i] <- Nay[j, i] * ((Sel_catch[j, 2] * FM[i] / 4) / Zay[j]) * (1 - exp(-Zay[j])) +
          Nay[j, i] * ((Sel_catch[j, 3] * FM[i] / 4) / Zay[j]) * (1 - exp(-Zay[j])) +
          Nay[j, i] * ((Sel_catch[j, 4] * FM[i] / 4) / Zay[j]) * (1 - exp(-Zay[j])) +
          Nay[j, i] * ((Sel_catch[j, 5] * FM[i] / 4) / Zay[j]) * (1 - exp(-Zay[j]))
      } # END J LOOP
      Nay[A, i + 1] <- Nay[A - 1, i] * exp(-Zay[A - 1]) +
        Nay[A, i] * exp(-Zay[A])

      Cay[A - 1, i] <- Nay[A - 1, i] * ((Sel_catch[A - 1, 2] * FM[i] / 4) / Zay[A - 1]) * (1 - exp(-Zay[A - 1])) +
        Nay[A - 1, i] * ((Sel_catch[A - 1, 3] * FM[i] / 4) / Zay[A - 1]) * (1 - exp(-Zay[A - 1])) +
        Nay[A - 1, i] * ((Sel_catch[A - 1, 4] * FM[i] / 4) / Zay[A - 1]) * (1 - exp(-Zay[A - 1])) +
        Nay[A - 1, i] * ((Sel_catch[A - 1, 5] * FM[i] / 4) / Zay[A - 1]) * (1 - exp(-Zay[A - 1]))

      Cay[A, i] <- Nay[A, i] * ((Sel_catch[A, 2] * FM[i] / 4) / Zay[A]) * (1 - exp(-Zay[A])) +
        Nay[A, i] * ((Sel_catch[A, 3] * FM[i] / 4) / Zay[A]) * (1 - exp(-Zay[A])) +
        Nay[A, i] * ((Sel_catch[A, 4] * FM[i] / 4) / Zay[A]) * (1 - exp(-Zay[A])) +
        Nay[A, i] * ((Sel_catch[A, 5] * FM[i] / 4) / Zay[A]) * (1 - exp(-Zay[A]))




      for (h in 1:A) {
        Iay1_CV[h, i] <- runif(1, min = CV1 - 0.1, max = CV1 + 0.1)
        Iay1[h, i] <- max(q_survey1[i] * Sel_survey[h, 2] * Nay[h, i] *
          exp(rnorm(1, 0, (q_survey1[i] * Sel_survey[h, 2] * Nay[h, i]) * Iay1_CV[h, i]) -
            (((q_survey1[i] * Sel_survey[h, 2] * Nay[h, i] * Iay1_CV[h, i])^2) / 2)), 0) # uncertainty based off CV
        Iay2_CV[h, i] <- runif(1, min = CV2 - 0.1, max = CV2 + 0.1)
        Iay2[h, i] <- max(q_survey2[i] * Sel_survey[h, 3] * Nay[h, i] *
          exp(rnorm(1, 0, (q_survey2[i] * Sel_survey[h, 3] * Nay[h, i]) * Iay2_CV[h, i]) -
            (((q_survey2[i] * Sel_survey[h, 3] * Nay[h, i] * Iay2_CV[h, i])^2) / 2)), 0)


        Iay3_CV[h, i] <- runif(1, min = CV3 - 0.1, max = CV3 + 0.1)
        Iay3[h, i] <- max(q_survey3[i] * Sel_survey[h, 2] * Nay[h, i] *
          exp(rnorm(1, 0, (q_survey3[i] * Sel_survey[h, 2] * Nay[h, i]) * Iay3_CV[h, i]) -
            (((q_survey3[i] * Sel_survey[h, 2] * Nay[h, i] * Iay3_CV[h, i])^2) / 2)), 0)

        Iay4_CV[h, i] <- runif(1, min = CV4 - 0.1, max = CV4 + 0.1)
        Iay4[h, i] <- max(q_survey4[i] * Sel_survey[h, 2] * Nay[h, i] *
          exp(rnorm(1, 0, (q_survey4[i] * Sel_survey[h, 2] * Nay[h, i]) * Iay4_CV[h, i]) -
            (((q_survey4[i] * Sel_survey[h, 2] * Nay[h, i] * Iay4_CV[h, i])^2) / 2)), 0)
      } # end h
    } # ends i loop

    # USE VON BERT CURVE TO ESTIMATE LENGTH COMP.
    # Lt = Linf * (1 - exp(-K*(t - t0)))
    Linf <- 81.6
    K <- 0.49
    t0 <- -0.97

    # L0 = Linf * (1 - exp(-K*(0-t0)))


    # GET LENGTH FREQUENCIES
    Cay_N <- cbind(ages, round(Cay))
    assign(paste("LF_Cay", k, sep = "_"), data.frame(LBins = seq(20, 122, by = 2)))
    for (f in 40:65) {
      dfC <- as.data.frame(Cay_N[, c(1, f + 1)])
      colnames(dfC) <- c("ages", "count")
      df.expandedC <- dfC[rep(row.names(dfC), dfC$count), 1:2]
      vonBertC <- rbvn(nrow(df.expandedC), m1 = Linf, s1 = Linf * 0.1, m2 = K, s2 = K * 0.1, rho = -0.9)
      df.expandedC$L <- vonBertC$Linf * (1 - exp(-(vonBertC$K) * (jitter(df.expandedC$ages, amount = 0.5) - t0)))
      hstC <- hist(df.expandedC$L, breaks = c(seq(20, 124, by = 2)), plot = F)
      assign(paste("LF_Cay", k, sep = "_"), cbind(get(paste("LF_Cay", k, sep = "_")), hstC$counts))
    } # end f loop

    Iay1_N <- cbind(ages, round(Iay1 * 100))
    assign(paste("LF_Iay1", k, sep = "_"), data.frame(LBins = seq(20, 122, by = 2)))
    for (l in 40:65) {
      df1 <- as.data.frame(Iay1_N[, c(1, l + 1)])
      colnames(df1) <- c("ages", "count")
      df.expanded1 <- df1[rep(row.names(df1), df1$count), 1:2]
      vonBert1 <- rbvn(nrow(df.expanded1), m1 = Linf, s1 = Linf * 0.1, m2 = K, s2 = K * 0.1, rho = -0.9)
      df.expanded1$L <- vonBert1$Linf * (1 - exp(-(vonBert1$K) * (jitter(df.expanded1$ages, amount = 0.5) - t0)))
      hst1 <- hist(df.expanded1$L, breaks = c(seq(20, 124, by = 2)), plot = F)
      assign(paste("LF_Iay1", k, sep = "_"), cbind(get(paste("LF_Iay1", k, sep = "_")), hst1$counts))
    } # end l loop

    Iay2_N <- cbind(ages, round(Iay2 * 100))
    assign(paste("LF_Iay2", k, sep = "_"), data.frame(LBins = seq(20, 122, by = 2)))
    for (m in 40:65) {
      df2 <- as.data.frame(Iay2_N[, c(1, m + 1)])
      colnames(df2) <- c("ages", "count")
      df.expanded2 <- df2[rep(row.names(df2), df2$count), 1:2]
      vonBert2 <- rbvn(nrow(df.expanded2), m1 = Linf, s1 = Linf * 0.1, m2 = K, s2 = K * 0.1, rho = -0.9)
      df.expanded2$L <- vonBert2$Linf * (1 - exp(-(vonBert2$K) * (jitter(df.expanded2$ages, amount = 0.5) - t0)))
      hst2 <- hist(df.expanded2$L, breaks = c(seq(20, 124, by = 2)), plot = F)
      assign(paste("LF_Iay2", k, sep = "_"), cbind(get(paste("LF_Iay2", k, sep = "_")), hst2$counts))
    } # end m loop


    Iay3_N <- cbind(ages, round(Iay3 * 100))
    assign(paste("LF_Iay3", k, sep = "_"), data.frame(LBins = seq(20, 122, by = 2)))
    for (n in 40:65) {
      df3 <- as.data.frame(Iay3_N[, c(1, n + 1)])
      colnames(df3) <- c("ages", "count")
      df.expanded3 <- df3[rep(row.names(df3), df3$count), 1:2]
      vonBert3 <- rbvn(nrow(df.expanded3), m1 = Linf, s1 = Linf * 0.1, m2 = K, s2 = K * 0.1, rho = -0.9)
      df.expanded3$L <- vonBert3$Linf * (1 - exp(-(vonBert3$K) * (jitter(df.expanded3$ages, amount = 0.5) - t0)))
      hst3 <- hist(df.expanded3$L, breaks = c(seq(20, 124, by = 2)), plot = F)
      assign(paste("LF_Iay3", k, sep = "_"), cbind(get(paste("LF_Iay3", k, sep = "_")), hst3$counts))
    } # end n loop

    Iay4_N <- cbind(ages, round(Iay4 * 100))
    assign(paste("LF_Iay4", k, sep = "_"), data.frame(LBins = seq(20, 122, by = 2)))
    for (u in 40:65) {
      df4 <- as.data.frame(Iay4_N[, c(1, u + 1)])
      colnames(df4) <- c("ages", "count")
      df.expanded4 <- df4[rep(row.names(df4), df4$count), 1:2]
      vonBert4 <- rbvn(nrow(df.expanded4), m1 = Linf, s1 = Linf * 0.1, m2 = K, s2 = K * 0.1, rho = -0.9)
      df.expanded4$L <- vonBert4$Linf * (1 - exp(-(vonBert4$K) * (jitter(df.expanded4$ages, amount = 0.5) - t0)))
      hst4 <- hist(df.expanded4$L, breaks = c(seq(20, 124, by = 2)), plot = F)
      assign(paste("LF_Iay4", k, sep = "_"), cbind(get(paste("LF_Iay4", k, sep = "_")), hst4$counts))
    } # end u loop



    # save length-frequencies
    LFreqList <- list()
    LFreqList[[paste("LF_Cay", k, sep = "_")]] <- get(paste("LF_Cay", k, sep = "_"))
    LFreqList[[paste("LF_Iay1", k, sep = "_")]] <- get(paste("LF_Iay1", k, sep = "_"))
    LFreqList[[paste("LF_Iay2", k, sep = "_")]] <- get(paste("LF_Iay2", k, sep = "_"))
    LFreqList[[paste("LF_Iay3", k, sep = "_")]] <- get(paste("LF_Iay3", k, sep = "_"))
    LFreqList[[paste("LF_Iay4", k, sep = "_")]] <- get(paste("LF_Iay4", k, sep = "_"))

    save(LFreqList, file = paste("LFreqList_", k, ".RData", sep = ""))




    Iy1[k, ] <- apply(Iay1, 2, sum)
    Iy2[k, ] <- apply(Iay2, 2, sum)
    Iy3[k, ] <- apply(Iay3, 2, sum)
    Iy4[k, ] <- apply(Iay4, 2, sum)
    lines(1:years, Iy1[k, ], type = "l", col = "grey")
    lines(1:years, Iy2[k, ], type = "l", col = "blue")
    lines(1:years, Iy3[k, ], type = "l", col = "red")
    lines(1:years, Iy4[k, ], type = "l", col = "green")
    Iy1_CV[k, ] <- apply(Iay1_CV, 2, mean, na.rm = T)
    Iy2_CV[k, ] <- apply(Iay2_CV, 2, mean, na.rm = T)
    Iy3_CV[k, ] <- apply(Iay3_CV, 2, mean, na.rm = T)
    Iy4_CV[k, ] <- apply(Iay4_CV, 2, mean, na.rm = T)

    Cy[k, ] <- apply(Cay, 2, sum)
    Ny[k, ] <- apply(Nay, 2, sum)


    # Save
    assign(paste("Nay", k, sep = "_"), Nay)
    assign(paste("Cay", k, sep = "_"), Cay)


    AnnualNC <- list()
    AnnualNC[[paste("Nay", k, sep = "_")]] <- Nay
    AnnualNC[[paste("Cay", k, sep = "_")]] <- Cay
    AnnualNC[[paste("Iay1", k, sep = "_")]] <- Iay1
    AnnualNC[[paste("Iay2", k, sep = "_")]] <- Iay2
    AnnualNC[[paste("Iay3", k, sep = "_")]] <- Iay3
    AnnualNC[[paste("Iay4", k, sep = "_")]] <- Iay4

    save(AnnualNC, file = paste("AnnualNC_", k, ".RData", sep = ""))
  } # end k loop



  ResultsList[["Iy1"]] <- Iy1
  ResultsList[["Iy2"]] <- Iy2
  ResultsList[["Iy3"]] <- Iy3
  ResultsList[["Iy4"]] <- Iy4

  ResultsList[["Iy1_CV"]] <- Iy1_CV
  ResultsList[["Iy2_CV"]] <- Iy2_CV
  ResultsList[["Iy3_CV"]] <- Iy3_CV
  ResultsList[["Iy4_CV"]] <- Iy4_CV


  ResultsList[["M0"]] <- M0
  ResultsList[["FSS"]] <- FSS
  ResultsList[["Npups"]] <- Npups

  ResultsList[["Cy"]] <- Cy
  ResultsList[["Ny"]] <- Ny

  save(ResultsList, file = "ResultsList.RData")



  ########## TEST DFA #############

  # library(MARSS)
  yrs <- 40:(years - 35)
  biasAll <- matrix(nrow = Niters, ncol = length(yrs))
  bias <- vector(length = Niters)
  biasAll.BT <- matrix(nrow = Niters, ncol = length(yrs))
  RMSE.BT <- vector(length = Niters)

  RMSE <- vector(length = Niters)
  FitRatios <- matrix(nrow = Niters, ncol = 4)
  FactorLoadings <- vector()
  DFATrends <- vector()
  DFATrendsBT <- vector()
  DFATrendsSE <- vector()
  DFATrendsSEBT <- vector()
  upCI_DFATrends <- vector()
  lowCI_DFATrends <- vector()

  datz_SD <- matrix(nrow = Niters, ncol = 4)
  # i=1

  for (i in 1:Niters) {
    assign(paste("dat", i, sep = ""), rbind(Iy1[i, 40:(years - 35)], Iy2[i, 40:(years - 35)], Iy3[i, 40:(years - 35)], c(rep(NA, 15), Iy4[i, 55:(years - 35)])))



    dat.a <- get(paste("dat", i, sep = ""))
    dat <- dat.a * c


    TT <- ncol(dat)
    N.ts <- nrow(dat)

    datL <- log(dat)
    y.bar <- apply(datL, 1, mean, na.rm = TRUE)

    dat.dm <- (datL - y.bar) / y.bar
    gsd <- sd(c(dat.dm), na.rm = TRUE)
    dat.z <- dat.dm / gsd
    # gsd
    # y.bar
    # apply(dat.z, 1, mean, na.rm=T)
    # apply(dat.z, 1, sd, na.rm=T)
    datz_SD[i, ] <- apply(dat.z, 1, sd, na.rm = T)

    dat.z <- as.matrix(dat.z)
    rownames(dat.z) <- c("Survey1", "Survey2", "Survey3", "Survey4")



    ##### SN DFA ####
    cntl.list <- list(minit = 200, maxit = 500000, abstol = 0.02, allow.degen = FALSE, conv.test.slope.tol = 0.05)
    R <- diag(c(CV1, CV2, CV3, CV4), nrow = 4, ncol = 4)

    dfa <- MARSS(dat.z, model = list(m = 1, R = R), control = cntl.list, form = "dfa", z.score = FALSE)

    pars <- MARSSparamCIs(dfa, nboot = 10000)

    # get the inverse of the rotation matrix
    # H.inv = varimax(coef(dfa, type="matrix")$Z)$rotmat # if n common trends > 1
    Z.rot <- coef(dfa, type = "matrix")$Z # %*% H.inv
    trends.rot <- dfa$states
    ts.trends <- t(trends.rot)
    par.mat <- coef(dfa, type = "matrix")
    fit.b <- par.mat$Z %*% dfa$states # + matrix(par.mat$D, nrow=N.ts) %*% covar # If DFA mod includes covariate

    assign(paste("Z.rot", i, sep = "."), Z.rot)
    assign(paste("Z.upCI", i, sep = "."), pars$par.upCI$Z)
    assign(paste("Z.lowCI", i, sep = "."), pars$par.lowCI$Z)
    assign(paste("trends.rot", i, sep = "."), trends.rot)


    #
    # # plot factor loadings
    par(mfrow = c(3, 2))
    survey <- rownames(dat.z)
    minZ <- 0.00
    ylims <- c(-1.1 * max(abs(Z.rot)), 1.1 * max(abs(Z.rot)))
    # par(mfrow=c(nrow(trends.rot),1), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
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
    yrs <- 40:(years - 35)

    index <- dfa$states[1, ]
    indexSE <- dfa$states.se[1, ]
    lowerCI <- dfa$states[1, ] - 1.96 * dfa$states.se[1, ]
    upperCI <- dfa$states[1, ] + 1.96 * dfa$states.se[1, ]
    assign(paste("index", i, sep = ""), index)
    assign(paste("indexSE", i, sep = ""), indexSE)
    assign(paste("lowerCI", i, sep = ""), lowerCI)
    assign(paste("upperCI", i, sep = ""), upperCI)

    # par(mfrow=c(2,2))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type = "n", bty = "L", xlab = "Year", ylab = "Abundance", ylim = c(min(lowerCI), max(upperCI))))
    x <- as.numeric(df$yrs)
    polygon(c(x, rev(x)), c(df$lowerCI, rev(df$upperCI)), col = "grey", border = NA)
    with(df, lines(x, index, type = "l", lwd = 2, pch = 16))
    title(paste("Iteration", i, sep = " "))
    abline(h = 0)

    N <- Ny[i, 40:65]
    NL <- log(N)
    Nscale <- (N - mean(N)) / sd(N)
    NLscale <- (N - mean(N)) / sd(N)
    index.z <- (index - mean(index)) / sd(index)
    index.z <- ifelse(is.na(index.z), 0, index.z)
    index.z <- ifelse(abs(index.z) == Inf, 0, index.z)


    indexSEBT <- indexSE * gsd
    indexBT <- exp(index * gsd + ((indexSEBT^2) / 2))
    assign(paste("indexBT", i, sep = ""), indexBT)
    assign(paste("indexSEBT", i, sep = ""), indexSEBT)

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
      title(paste("Atlantic sharpnose", survey[n], sep = " "))
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
    "FL1", "FL2", "FL3", "FL4",
    "lowCI_FL1", "lowCI_FL2", "lowCI_FL3", "lowCI_FL4",
    "upCI_FL1", "upCI_FL2", "upCI_FL3", "upCI_FL4"
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
      plot(1:years, Iy1[l, ],
        type = "l", col = "black", xlim = c(min(yrs), max(yrs)), axes = F,
        ylim = c(0, max(c(Iy1[l, yrs], Iy2[l, yrs], Iy3[l, yrs])) + 0.1)
      )
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      lines(1:years, Iy2[l, ], type = "l", col = "blue")
      lines(1:years, Iy3[l, ], type = "l", col = "red")
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
        plot(c(1:N.ts)[abs(Z.rot[, h]) >= minZ], as.vector(Z.rot[abs(Z.rot[, h]) >= minZ, h]),
          type = "h", lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylims, xlim = c(0, N.ts + 1), axes = FALSE
        )
        axis(1, labels = FALSE)
        axis(2, labels = FALSE)
        #   for(j in 1:N.ts) {
        #     if(Z.rot[j,h] > minZ) {text(j, -0.05, survey[j], srt=90, adj=1, cex=0.9)}
        #     if(Z.rot[j,h] < -minZ) {text(j, 0.05, survey[j], srt=90, adj=0, cex=0.9)}
        #   } # end j loop
        #   mtext(paste("Factor loadings on trend",h,sep=" "),side=3,line=.5)
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
      plot(yrs, NL, type = "l", col = "red", axes = FALSE)
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)
      lines(yrs, rescale(index, to = c(min(NL), max(NL))))
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
    }

    dev.off()
  } # end if plot=T





  ########## END ########

  return(list(RMSE = RMSE, biasAll = biasAll, FitRatios = FitRatios, RMSE.BT = RMSE.BT, biasAll.BT = biasAll.BT))
} # END FUNCTION





########### EXAMPLE APPLICATION ###########

# F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))) # Trials 19-24
setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial19")
ATrial19 <- DFA_Sim(
  F_var = as.vector(c(rep(0, 40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))),
  CV1 = 0.5,
  CV2 = 0.5,
  CV3 = 0.5,
  q_survey1 = rep(0.001, length = years),
  q_survey2 = rep(0.001, length = years),
  q_survey3 = rep(0.001, length = years),
  Iay_eq = Iay_eq_F4,
  c = c(15, 10, 10)
)
save(ATrial19, file = "ATrial19.RData")





###### Load data / biology & fishery info ####

# Age classes: 18
# Longevity: 23
# Atlantic coast: age at 50% maturity; A50 = 3
# GOM: A50 = 1.6 -- used in 2013 assessment
# restricted migratory patterns around florida; two genetically singular stock
# M
# Virgin conditions: 1950
# Consider stock-recruitment relationship following a functional form...

# Reference Porch 2003: Preliminary assessment of atlantic white marlin (Tetrapturus albidus) using a state-space
#   implementation of an age-structured production model

# SEDAR 34 estimates steepness of 0.57.


# 3 catch streams included:
#   3 commercial: bottom longline, gillnets, lines
#   recreatioanl catches
#   shrimp bycatch

# Von Bert params
Linf <- 81.6
K <- 0.49
t0 <- -0.97

# VON BERT:  Lt = Linf *(1-exp(-K*(t-t0)))
# SCHNUTE: L1 + (L2 - L1) * ((1-exp(-K*(t-t1)))/(1-exp(-K*(t2-t1))))
L0 <- Linf * (1 - exp(-K * (0 - t0))) # L0 = 30.8694
L1 <- Linf * (1 - exp(-K * (1 - t0))) # L1 = 50.5211


# Read in selectivity
# AGE-BASED SELECTIVITY:
# Sel_catch = read.csv("D:\\vspace1\\DFA_Simulation\\Atl_SN\\CatchSelectivity.csv")
# Sel_survey = read.csv("D:\\vspace1\\DFA_Simulation\\Atl_SN\\SurveySelectivity.csv")



# ages
ages <- 1:18
A <- max(ages)

# Define maturity function
mat <- as.vector(c(0.185, 0.953, 0.999, rep(1, 15)))

# Define fecundity (female pups)
fec <- c(0.401, 0.762, 1.133, 1.448, 1.686, 1.852, 1.963, 2.035, 2.081, 2.11, 2.128, 2.139, 2.146, 2.15, 2.153, 2.155, 2.156, 2.156)




# define natural mortality
M_const <- as.vector(rep(0.209, A))

# Equilibrium numbers at age
Equilibrium_Nay <- c(494.47970, 401.21848, 325.54676, 264.14709, 214.32769, 173.90447, 141.10526, 114.49213, 92.89837, 75.37730, 61.16078, 49.62556, 40.26595, 32.67160, 26.50958, 21.50975, 17.45290, 75.08402)



# GET LFSR PARAMs

z0 <- 1.904416
s0 <- exp(-z0)
s0 # 0.1489096
# zmin = 0.06157842
zmin <- 0.209
(smax <- exp(-zmin)) # 0.8113952
(zfrac <- (log(smax) - log(s0)) / (-log(s0))) # 0.8902551
h <- 0.56
beta <- (log(1 - (log(h / 0.2) / (z0 * zfrac)))) / log(0.2)
beta
# zmax = 2.0
(S_frac <- (zmin - z0) / (0 - z0)) # 0.8902551

Z_max <- z0 + S_frac * (0.0 - z0) # 0.209




# length of simulation
years <- 100


# Fishing mortality will vary over time. For now, set it to a constant 0.1

# define natural mortality
M_const <- as.vector(rep(0.209, A))

# Assume constant fishing mortality
F_var <- as.vector(rep(0.0, years))
F_var <- as.vector(rep(0.2, years)) # Trials 1-6
F_var <- as.vector(c(rep(0.4, 50), rep(0, 50))) # Trials 13-18
F_var <- as.vector(c(rep(0.0, 50), rep(0.4, 50))) # Trials 7-12
F_var <- as.vector(c(rep(0, 40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))) # Trials 19-24


## LOAD EQUILIBRIUM SURVEYS from saved .Rdata
# load("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\I_eq_list.RData")
# load("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\I_eq_list.RData")
#
# Iay_eq_F1 = I_eq_list$Iay_eq_F1
# Iay_eq_F2 = I_eq_list$Iay_eq_F2
# Iay_eq_F3 = I_eq_list$Iay_eq_F3
# Iay_eq_F4 = I_eq_list$Iay_eq_F4
