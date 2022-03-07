################### START HERE ################
# DFA TUTORIAL FROM ICES MANUSCRIPT: DOI: 10.1093/icesjms/fsab051

## This is Trial SB 3 from ICES manuscript
# all catchabilities constant except Survey 2 (I2) experiences decrease in q over time

library(MARSS)
library(scales)
# setwd("O:/Users/cassidy.peterson/My Documents/DFA_Simulation/SB/Trial3")
load('dat.RData') # get dataset called dat.a
dat.a
load('Ny.RData') # get dataset called dat.a
Ny
CV1=0.38; CV2=0.48; CV3=0.65; CV4=0.24; CV5=0.30; CV6=0.36; CV7=0.40;
yrs = 51:90



  # recall that we have a unique way of "standardizing" the indices to allow for back-transformation out of log-space.
  # 1. multiply each index by a constant  -- this needs to be iteratively determined to meet "c requirements"
  # 2. de-mean
  # 3. calculate GSD (global standard deviation)
  # 4. standardize (demeaned indices / GSD)


  # DEFINE c VECTOR -- this will have to be determined iteratively, start with something basic
# c = c(5, 5, 5, 5, 5, 5, 5) #
c = c(18, 360, 1, 0.45, 10.5, 0.79, 0.37) # ending point


  dat = dat.a * c # multiply dat.a by vector of constants (c)

  TT=ncol(dat)
  N.ts = nrow(dat)

  # Standardize data
  datL=log(dat)
  y.bar = apply(datL, 1, mean, na.rm=TRUE) # get mean of each index; ENSURE y.bar >0
  dat.dm = (datL-y.bar) / y.bar  # get demeaned index
  gsd = sd(c(dat.dm), na.rm=T) # get GSD or global standard deviation for all indices
  # ensure gsd<<1; if gsd>1, try new vector c
  dat.z = dat.dm/gsd # calculate "z-scored" indices (using demeaned indices and GSD)

  # REQUIREMENTS FOR C -- change c vector 
  y.bar # means must be >0 (otherwise demeaning will change sign/direction of index); means should be >1; ideally >~2 (which isn't always possible)
  gsd # gsd should be small! <<~0.5
  apply(dat.z, 1, sd, na.rm=T) # MOST IMPORTANT DIAGNOSTIC: sd should be 1.
  # this rescalign approach is to approximate a z-score. 
  # if sd for an index is <1, then increase corresponding element in c matrix; vice versa || index specific c values
  # helpful to tune 1 at a time
  min(datL, na.rm=T) # ensure that minimum datL >0

  #### STOP ###################################################
  # IF REQUIREMENTS FOR c ARE NOT MET, THEN ITERATIVELY MODIFY c VECTOR UNTIL REQUIREMENTS ARE MET!




  ##### RUN DFA
  datz_SD = apply(dat.z, 1, sd, na.rm=T)
  dat.z = as.matrix(dat.z)
  rownames(dat.z) = c("Survey1","Survey2","Survey3","Survey4","Survey5","Survey6","Survey7")

  cntl.list = list(minit=200, maxit=500000, abstol=0.02, allow.degen=FALSE, conv.test.slope.tol = 0.05)
  R = diag(c(CV1, CV2, CV3, CV4, CV5, CV6, CV7), nrow=7, ncol=7)

  dfa = MARSS(dat.z, model = list(m=1,R=R), control=cntl.list, form="dfa", z.score=FALSE) # z.score=T would re-z-score for you and we don't want that
  # we do this to get a single global sd -- z.score=T would give us an index-specific z-score. not a global z-score. 
  

  # Get results
  pars = MARSSparamCIs(dfa, nboot=10000)

  # IF calculating more than 1 common trend/common factor
  ## get the inverse of the rotation matrix
  ##H.inv = varimax(coef(dfa, type="matrix")$Z)$rotmat
  Z.rot = coef(dfa, type="matrix")$Z #%*% H.inv # H.inv not necessary if 1 common trend
  trends.rot = dfa$states
  ts.trends = t(trends.rot)
  par.mat=coef(dfa, type="matrix")
  fit.b = par.mat$Z %*% dfa$states # + matrix(par.mat$D, nrow=N.ts) %*% COVARIATE # only necessary if includign covariates (e.g., regional SST or NAO, etc.) to the model

  # Calculate CI's
  index = dfa$states[1,]
  indexSE = dfa$states.se[1,]
  lowerCI = dfa$states[1,]-1.96*dfa$states.se[1,]
  upperCI = dfa$states[1,]+1.96*dfa$states.se[1,]

  # Calculate simulated N and rescale
  N = Ny[51:90]
  Nscale = (N-mean(N))/sd(N)
  NL = log(N)
  NLscale = (NL - mean(NL))/sd(NL)
  index.z = (index - mean(index))/sd(index)
  index.z = ifelse(is.na(index.z), 0, index.z)
  index.z = ifelse(abs(index.z)==Inf, 0, index.z)

  # back-transform index and SE
  indexSEBT = indexSE * gsd
  indexBT = exp(index*gsd + ((indexSEBT^2)/2) )



 ### PLOTTING ###
  # PLOT factor loadings
  survey = rownames(dat.z)
  minZ = 0.00
  ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
  par(mfrow=c(nrow(trends.rot),1), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
  for(h in 1:nrow(trends.rot)) {
    plot(c(1:N.ts)[abs(Z.rot[,h])>minZ], as.vector(Z.rot[abs(Z.rot[,h])>minZ,h]),
         type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1))
    for(j in 1:N.ts) {
      if(Z.rot[j,h] > minZ) {text(j, -0.05, survey[j], srt=90, adj=1, cex=0.9)}
      if(Z.rot[j,h] < -minZ) {text(j, 0.05, survey[j], srt=90, adj=0, cex=0.9)}
      abline(h=0, lwd=1, col="gray")
      abline(h=0.2, col='gray', lty=2)
      abline(h=-0.2, col='gray', lty=2)
    } # end j loop
    mtext(paste("Factor loadings on trend",h,sep=" "),side=3,line=.5)
  } # end h loop




  # PLOT COMMON TREND w/ CIs
  par(mfrow=c(4,2))
  df <- data.frame(yrs, index, lowerCI, upperCI)
  with(df, plot(index ~ yrs, type='n', bty='L',xlab='Year', ylab='Abundance', ylim=c(min(lowerCI), max(upperCI))))
  x <- as.numeric(df$yrs)
  polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
  with(df, lines(x, index, type='l', lwd=2, pch=16))
  title(paste("Common Trend"))
  abline(h=0)

  # superimpose simulated index
  lines(yrs,rescale(NL,to=c(min(index),max(index))), type='l', col='red', lwd=2) ## Note -- this is rescaling to the min/max of the common trend, consider that rescaling to other values may be more appropriate
  lines(yrs,rescale(NL,to=c(-2.5,15)), type='l', col='darkred', lwd=1, lty=2) # for example...




  # Plot fitted values
  survey = rownames(dat.z)
  # par(mfcol=c(ceiling(N.ts/2),2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
  for(n in 1:length(survey)){
    plot(yrs,dat.z[n,],xlab="",ylab="abundance index",bty="L", xaxt="n", ylim=c(-4,4), pch=16, col="blue")
    axis(1,labels=TRUE)
    lines(yrs,fit.b[n,], lwd=2)
    lines(yrs,rescale(NL, to=c(min(fit.b[n,]), max(fit.b[n,]))), type='l', col='red', lwd=2)
    title(paste("Sandbar",survey[n], sep=" "))
  }



  #### ASSESS DFA MODEL FITS using post-hoc test
  sumResids = rowSums((dat.z-fit.b)^2, na.rm=TRUE)
  sumObserved = rowSums(dat.z^2, na.rm=TRUE)
  FitRatio = sumResids / sumObserved ; FitRatio
  FitRatio
  mean(FitRatio)

  rbind(indexBT, dat.a)
  # indexBT
  # dat.a


  # CALCULATE BIAS AND RMSE
  biasAll = index.z - NLscale
  bias = mean(index.z - NLscale)
  RMSE = sqrt(sum((index.z-NLscale)^2)/length(index.z))

  indexBT.z = (indexBT-mean(indexBT))/sd(indexBT)
  biasAll.BT = indexBT.z - Nscale
  RMSE.BT = sqrt(sum((indexBT.z-Nscale)^2)/length(indexBT.z))






############################ IGNORE BELOW ########################
  # # these are notes that I used to generate example data
  # # setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial2")
  # setwd("O:\\Users\\cassidy.peterson\\My Documents\\DFA_Simulation\\SB\\Trial3")
  # load("DFA_Results.RData")
  # load("ResultsList.RData")
  #
  # Iy1<-ResultsList$Iy1
  # Iy2<-ResultsList$Iy2
  # Iy3<-ResultsList$Iy3
  # Iy4<-ResultsList$Iy4
  # Iy5<-ResultsList$Iy5
  # Iy6<-ResultsList$Iy6
  # Iy7<-ResultsList$Iy7
  #
  # Niters=1
  # i=1
  # yrs = 51:90
  #
  #
  # I1 = c(Iy1[i,51:53], rep(NA, 5), Iy1[i,59:64], NA, NA, Iy1[i,67:90])
  # I2 = c(rep(NA, 7), Iy2[i,58:90])
  # I3 = c(rep(NA, 13), Iy3[i,64:90])
  # I4 = c(rep(NA, 16), Iy4[i,67:90])
  # I5 = c(rep(NA, 22), Iy5[i,73:90])
  # I6 = c(rep(NA, 15), Iy6[i,66:79], rep(NA, 11))
  # I7 = c(rep(NA, 17), Iy7[i,68], NA, NA, Iy7[i,71], NA, NA, Iy7[i,74], NA, NA, Iy7[i,77], NA, NA, Iy7[i,80], NA, NA, Iy7[i,83], NA, NA, Iy7[i,86], NA, NA, Iy7[i,89], NA)
  #
  # assign(paste("dat",i,sep=""), rbind(I1, I2, I3, I4, I5, I6, I7))
  #
  # dat1
  #
  # dat.a = get(paste("dat",i,sep="")) #
  # save(dat.a, file="dat.RData")
  #
  #
  # load("ResultsList.RData")
  # Ny<-ResultsList$Ny[1,]
  # save(Ny, file="Ny.RData")
