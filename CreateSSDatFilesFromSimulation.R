
# C Peterson
# July 2018
# Automating Atl_SN simulations to SS data input files
# how to change source file?

# read in packages
library(r4ss)
library(scales)

#############################

base_dat = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")


# base_dat = SS_readdat(file="C:/Users/cpeterson/Documents/DFA_Simulation/Atl_SN_324/snA71.dat")
# base_starter = SS_readstarter(file="C:/Users/cpeterson/Documents/DFA_Simulation/Atl_SN_324/starter.ss")
# base_ctl = SS_readctl(file="C:/Users/cpeterson/Documents/DFA_Simulation/Atl_SN_324/snA71.ctl")
# # REQUIRES LONG BAYESIAN PRIOR LINES-- Doesn't recognize short lines
# SS_parlines(ctlfile="C:/Users/cpeterson/Documents/DFA_Simulation/Atl_SN_324/control.ss_new")
# # SS_changepars(ctlfile="C:/Users/cpeterson/Documents/DFA_Simulation/Atl_SN_324/control.ss_new")

# names(base_dat)
# base_dat$CPUE
# class(base_dat$CPUE)
# 
# x = data.frame(x=c(1993:2000), seas=c(rep(1,8)), index=c(rep(2,8)), 
#                obs=c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8), se_log=c(rep(1.5,8)))
# 
# SN_new_data = base_dat
# SN_new_data$CPUE = x
# SN_new_data
# 
# SS_writedat(datlist = SN_new_data, outfile="C:/Users/cpeterson/Documents/DFA_Simulation/NEW_SS_dat_file.dat")


### NEED TO FIX THIS... 
SS_changepars(ctlfile="DFA_Simulation/Atl_SN_324/control.ss_new", 
              newctlfile=paste("C:/Users/cpeterson/Documents/DFA_Simulation/Atl_SN/", trial, "/SS/SN_SS_", i, "/SN.ctl", sep=""),
              linenums=NULL)


Niters = 100

# FOR WHEN Nindices = 3

# base_dat = SS_readdat(file="C:/Users/cpeterson/Documents/DFA_Simulation/Atl_SN_324/snA71.dat")
# # For DFA index: 
# base_DFA_dat = SS_readdat(file="C:/Users/cpeterson/Documents/DFA_Simulation/Atl_SN_324/snA71.dat")


setwd("C:/Users/cpeterson/Documents/DFA_Simulation/Atl_SN/Base/Trial1")
paste("C:/Users/cpeterson/Documents/DFA_Simulation/Atl_SN/", trial, "/SS/SN_SS_", i, "/SN.dat", sep="")


base_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")

base_dat_10 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial10\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_10 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial10\\SS_DFA\\SN_SS_1\\SN.dat")

base_dat_22 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial22\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_22 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial22\\SS_DFA\\SN_SS_1\\SN.dat")


# TEST
trial = "Base/Trial22"
base_dat = base_dat_22
base_DFA_dat = base_DFA_dat_22
Niters = 1
age = FALSE
length = TRUE 
Nindices=3
dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/"
EqCatch=NULL



####################################################################################################################################
# ### WRITE FUNCTION #### NEW 2019
####################################################################################################################################

datfile = function(dir = dir, trial = "Base/Trial1", base_dat = base_dat, base_DFA_dat = base_DFA_dat, 
                   Niters = 100,  length = TRUE, Nindices=3, EqCatch=NULL, ...){
  
  setwd(paste(dir,trial, sep=""))
  
  new_dat = base_dat
  new_DFA_dat = base_DFA_dat
  
  
  load("ResultsList.RData")
  load("DFA_Results.RData")
  
  
  if(Nindices ==3){
    
    Cy = ResultsList$Cy
    Iy1 = ResultsList$Iy1
    Iy2 = ResultsList$Iy2
    Iy3 = ResultsList$Iy3
    Iy4 = ResultsList$Iy4
    Iy1_CV = ResultsList$Iy1_CV
    Iy2_CV = ResultsList$Iy2_CV
    Iy3_CV = ResultsList$Iy3_CV
    Iy4_CV = ResultsList$Iy4_CV
    DFATrends = DFA_Results$DFATrends
    DFATrendsSE = DFA_Results$DFATrendsSE
    FLoadings = DFA_Results$FLoadings
    
    
    for(i in 1:Niters){
      
      # catch
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[40:65], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV= unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 3), seas=rep(1, (2018-1992)*3), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992)), 
                              obs = c(I1, I2, I3), se_log = c(se_log1, se_log2, se_log3))
      
      new_dat$CPUE = index.data
      
      
      #for DFA data
      # I_DFA = unlist(DFATrends[i,] - min(DFATrends[i,]) + 0.1 )
      CV_DFA = unlist(DFATrendsSE[i, ]) / unlist(DFATrends[i,])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # mean_raw = mean(c(I1, I2, I3))
      # min_raw = min(c(I1, I2, I3))
      # trend = unlist(DFATrends[i,])
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I2), max(I2))) # BASE
      I_DFA = unlist(DFATrends[i,] + (-2 * min(DFATrends[i,])) ) # TRY THIS!!! MUCH BETTER!
      
      # mean2 = mean(I_DFA)
      # I_DFA2 = I_DFA + (mean_raw - mean2)
      # I_DFA = I_DFA2
      #mean(I_DFA2)
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      
      
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
      new_DFA_dat$CPUE = index_DFA.data
      
      
      ### Length Comp ####
      if(length == TRUE){
        
        load(paste("LFreqList_",i,".RData",sep=""))
        
        # From commercial catch
        lencomp_C=LFreqList[[1]]
        lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,2:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(1, 2018-1992)), 
                           Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        # From Surveys
        lencomp_I1=LFreqList[[2]]
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        
        lencomp_I2=LFreqList[[3]]
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        
        lencomp_I3 = LFreqList[[4]]
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        
        # for DFA assessment
        FL = FLoadings[i,1:3]
        lambda = data.frame(lambda=ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda$lambda[1])*LF_I1) + round((as.numeric(lambda$lambda[2])*LF_I2)) + round(as.numeric(lambda$lambda[3])*LF_I3)
        # LF_DFA = LF_I1+LF_I2+LF_I3
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      # SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      # SS_writedat(new_DFA_dat, outfile = "N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/Base/Trial10a/SS_DFA/SN_SS_1/SN.dat", overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices==3)
  
  
  
  ### 4 survey indices
  if(Nindices == 4){
   
    Cy = ResultsList$Cy
    Iy1 = ResultsList$Iy1
    Iy2 = ResultsList$Iy2
    Iy3 = ResultsList$Iy3
    Iy4 = ResultsList$Iy4
    Iy1_CV = ResultsList$Iy1_CV
    Iy2_CV = ResultsList$Iy2_CV
    Iy3_CV = ResultsList$Iy3_CV
    Iy4_CV = ResultsList$Iy4_CV
    DFATrends = DFA_Results$DFATrends
    DFATrendsSE = DFA_Results$DFATrendsSE
    FLoadings = DFA_Results$FLoadings
    
    for(i in 1:Niters){
      
      # catch
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[40:65], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      I4 = unlist(Iy4[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV = unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      I4_CV = unlist(Iy4_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      se_log4 = sqrt(log(1+(I4_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 4), seas=rep(1, (2018-1992)*4), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992), rep(5, 2018-1992)), 
                              obs = c(I1, I2, I3, I4), se_log = c(se_log1, se_log2, se_log3, se_log4))
      
      new_dat$CPUE = index.data
      new_dat$N_cpue = nrow(index.data)
      new_dat$Nsurveys = Nindices
      new_dat$fleetnames = c("Comm_Fish", "Survey1", "Survey2", "Survey3", "Survey4")
      new_dat$surveytiming = c(-1, 0.5, 0.5, 0.5, 0.5)
      new_dat$areas = c(rep(1,5))
      
      #for DFA data
      # I_DFA = unlist(DFATrends[i,41:66])
      if(DFATrends[i,2] > 1e-100) {#unlist(DFATrends[i,2]) < 1e-100
        CV_DFA = unlist(DFATrendsSE[i, ]) / unlist(DFATrends[i,])
        SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
        # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I2), max(I2)))
        I_DFA = unlist(DFATrends[i,] + (-2 * min(DFATrends[i,])) ) # TRY THIS!!! MUCH BETTER!
      }
      if(DFATrends[i,2] < 1e-100) {#unlist(DFATrends[i,2]) < 1e-100
        CV_DFA = unlist(DFATrendsSE[i, ]) / 1
        SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
        # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I2), max(I2)))
        I_DFA = rep(1, ncol(DFATrends)-1) 
      }
      # SE_DFA = unlist(DFATrendsSE[i, 41:66]) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # #for DFA data
      # DFATrends = read.csv("DFATrends.csv")
      # DFATrendsSE = read.csv("DFATrendsSE.csv")
      # # I_DFA = unlist(DFATrends[i,] - min(DFATrends[i,]) + 0.1 )
      # CV_DFA = DFATrendsSE[i, ] / DFATrends[i,]
      # SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I2), max(I2)))
      # # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
      new_DFA_dat$CPUE = index_DFA.data
      new_DFA_dat$N_cpue = nrow(index_DFA.data)
      
      
      ### Length Comp ####
      if(length == TRUE){
        
        load(paste("LFreqList_",i,".RData",sep=""))
        
        # From commercial catch
        lencomp_C=LFreqList[[1]]
        lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,2:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(1, 2018-1992)), 
                           Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        # From Surveys
        lencomp_I1=LFreqList[[2]]
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        
        lencomp_I2=LFreqList[[3]]
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        
        lencomp_I3 = LFreqList[[4]]
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        
        lencomp_I4 = LFreqList[[5]]
        colnames(lencomp_I4) = NULL
        LF_I4 = as.data.frame(t(lencomp_I4[,2:ncol(lencomp_I4)]))
        Lfreq_I4 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(5, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I4, 1, sum)), LF_I4)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3, Lfreq_I4)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        # for DFA assessment
        FL = FLoadings[i,1:4]
        lambda = data.frame(lambda=ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda$lambda[1])*LF_I1) + round((as.numeric(lambda$lambda[2])*LF_I2)) + round(as.numeric(lambda$lambda[3])*LF_I3) + round(as.numeric(lambda$lambda[4])*LF_I4)
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      # SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
} # end datfile function

datfileF2 = function(dir = dir, trial = "Base/Trial1", base_dat = base_dat, base_DFA_dat = base_DFA_dat, 
                     Niters = 100, length = TRUE, Nindices=3, EqCatch=NULL, ...){
  
  
  setwd(paste(dir,trial, sep=""))
  
  new_dat = base_dat
  new_DFA_dat = base_DFA_dat
  
  
  load("ResultsList.RData")
  load("DFA_Results.RData")
  
  
  if(Nindices ==3){
   
    Cy = ResultsList$Cy
    Iy1 = ResultsList$Iy1
    Iy2 = ResultsList$Iy2
    Iy3 = ResultsList$Iy3
    Iy4 = ResultsList$Iy4
    Iy1_CV = ResultsList$Iy1_CV
    Iy2_CV = ResultsList$Iy2_CV
    Iy3_CV = ResultsList$Iy3_CV
    Iy4_CV = ResultsList$Iy4_CV
    DFATrends = DFA_Results$DFATrends
    DFATrendsSE = DFA_Results$DFATrendsSE
    FLoadings = DFA_Results$FLoadings
    
    
    for(i in 1:Niters){
      
      # catch
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[40:65], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV= unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 3), seas=rep(1, (2018-1992)*3), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992)), 
                              obs = c(I1, I2, I3), se_log = c(se_log1, se_log2, se_log3))
      
      new_dat$CPUE = index.data
      
      
      #for DFA data
      # I_DFA = unlist(DFATrends[i,] - min(DFATrends[i,]) + 0.1 )
      CV_DFA = unlist(DFATrendsSE[i, ]) / unlist(DFATrends[i,])
      # SE_DFA = unlist(DFATrendsSE[i, ])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) 
      # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # mean_raw = mean(c(I1, I2, I3))
      # min_raw = min(c(I1, I2, I3))
      # trend = unlist(DFATrends[i,])
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I2), max(I2))) ## BASE USE 
      I_DFA = unlist(DFATrends[i,] + -2* min(DFATrends[i,])  ) # THIS ONE!
      # I_DFA = I2
      # SE_DFA = se_log2
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(0,1)) + 0.1
      # mean2 = mean(I_DFA)
      # I_DFA2 = trend + max(c(I1, I2, I3))
      # I_DFA2 = trend - (min(trend) - min_raw)
      # I_DFA2 = I_DFA + (mean_raw - mean2)
      # I_DFA = I_DFA2
      # mean(I_DFA2)
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I2)+(mean_raw - mean2), max(I2)))
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      
      
      
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)
      new_DFA_dat$CPUE = index_DFA.data
      
      
      ### Length Comp ####
      if(length == TRUE){
        
        load(paste("LFreqList_",i,".RData",sep=""))
        
        # From commercial catch
        lencomp_C=LFreqList[[1]]
        lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,13:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=2004:2018, seas = (rep(1, 2018-2003)), FltSry = (rep(1, 2018-2003)), 
                           Sex = (rep(0, 2018-2003)), Part = (rep(0, 2018-2003)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        # From Surveys
        lencomp_I1=LFreqList[[2]]
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        
        lencomp_I2=LFreqList[[3]]
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
       
        lencomp_I3 = LFreqList[[4]]
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3)
        # tail(lencomp.dat, 20)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        # for DFA assessment
        FL = FLoadings[i,1:3]
        # FL = as.vector(c(1,1,1))
        # lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
        lambda = as.vector(lambda=ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda$lambda[1])*LF_I1) + round((as.numeric(lambda$lambda[2])*LF_I2)) + round(as.numeric(lambda$lambda[3])*LF_I3)
        # LF_DFA = (LF_I1) + (LF_I2) + (LF_I3)
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        
      } # end if length == TRUE
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      # SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      # SS_writedat(new_DFA_dat, outfile = "N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/Base/Trial10h/SS_DFA/SN_SS_1/SN.dat", overwrite=TRUE)
      
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
  ### 4 survey indices
  if(Nindices == 4){
    
    Cy = ResultsList$Cy
    Iy1 = ResultsList$Iy1
    Iy2 = ResultsList$Iy2
    Iy3 = ResultsList$Iy3
    Iy4 = ResultsList$Iy4
    Iy1_CV = ResultsList$Iy1_CV
    Iy2_CV = ResultsList$Iy2_CV
    Iy3_CV = ResultsList$Iy3_CV
    Iy4_CV = ResultsList$Iy4_CV
    DFATrends = DFA_Results$DFATrends
    DFATrendsSE = DFA_Results$DFATrendsSE
    FLoadings = DFA_Results$FLoadings
    
    
    for(i in 1:Niters){
      
      # catch
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[40:65], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      I4 = unlist(Iy4[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV = unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      I4_CV = unlist(Iy4_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      se_log4 = sqrt(log(1+(I4_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 4), seas=rep(1, (2018-1992)*4), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992), rep(5, 2018-1992)), 
                              obs = c(I1, I2, I3, I4), se_log = c(se_log1, se_log2, se_log3, se_log4))
      
      new_dat$CPUE = index.data
      new_dat$N_cpue = nrow(index.data)
      new_dat$Nsurveys = Nindices
      new_dat$fleetnames = c("Comm_Fish", "Survey1", "Survey2", "Survey3", "Survey4")
      new_dat$surveytiming = c(-1, 0.5, 0.5, 0.5, 0.5)
      new_dat$areas = c(rep(1,5))
      
      
      #for DFA data
      CV_DFA = unlist(DFATrendsSE[i, ]) / unlist(DFATrends[i,])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) 
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I2), max(I2)))
      I_DFA = unlist(DFATrends[i,] + (-2 * min(DFATrends[i,])) ) # TRY THIS!!! MUCH BETTER!
      
      # I_DFA = unlist(DFATrends[i,41:66])
      # SE_DFA = unlist(DFATrendsSE[i, 41:66]) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      
      
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
      new_DFA_dat$CPUE = index_DFA.data
      
      
      ### Length Comp ####
      if(length == TRUE){
        
        load(paste("LFreqList_",i,".RData",sep=""))
        
        # From commercial catch
        lencomp_C=LFreqList[[1]]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,13:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=2004:2018, seas = (rep(1, 2018-2003)), FltSry = (rep(1, 2018-2003)), 
                           Sex = (rep(0, 2018-2003)), Part = (rep(0, 2018-2003)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        
        # From Surveys
        lencomp_I1=LFreqList[[2]]
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        
        lencomp_I2=LFreqList[[3]]
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        
        lencomp_I3 = LFreqList[[4]]
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        
        lencomp_I4 = LFreqList[[5]]
        colnames(lencomp_I4) = NULL
        LF_I4 = as.data.frame(t(lencomp_I4[,2:ncol(lencomp_I4)]))
        Lfreq_I4 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(5, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I4, 1, sum)), LF_I4)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3, Lfreq_I4)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        # for DFA assessment
        FL = FLoadings[i,1:4]
        lambda = data.frame(lambda=ifelse(FL<0, 0, FL/max(FL)))
        LF_DFA = round(as.numeric(lambda$lambda[1])*LF_I1) + round(as.numeric(lambda$lambda[2])*LF_I2) + 
          round(as.numeric(lambda$lambda[3])*LF_I3) + round(as.numeric(lambda$lambda[4])*LF_I4)
        # LF_DFA = LF_I1+LF_I2+LF_I3+LF_I4
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      # SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
} # end datfile function

datfileF3 = function(dir = dir, trial = "Base/Trial1", base_dat = base_dat, base_DFA_dat = base_DFA_dat, 
                     Niters = 100,  length = TRUE, Nindices=3, EqCatch=NULL, ...){
  
  setwd(paste(dir,trial, sep=""))
  
  new_dat = base_dat
  new_DFA_dat = base_DFA_dat
  
  
  load("ResultsList.RData")
  load("DFA_Results.RData")
  
  
  if(Nindices ==3){
    Cy = ResultsList$Cy
    Iy1 = ResultsList$Iy1
    Iy2 = ResultsList$Iy2
    Iy3 = ResultsList$Iy3
    Iy4 = ResultsList$Iy4
    Iy1_CV = ResultsList$Iy1_CV
    Iy2_CV = ResultsList$Iy2_CV
    Iy3_CV = ResultsList$Iy3_CV
    Iy4_CV = ResultsList$Iy4_CV
    DFATrends = DFA_Results$DFATrends
    DFATrendsSE = DFA_Results$DFATrendsSE
    FLoadings = DFA_Results$FLoadings
    
    for(i in 1:Niters){
      
      # catch
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[40:65], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV= unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 3), seas=rep(1, (2018-1992)*3), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992)), 
                              obs = c(I1, I2, I3), se_log = c(se_log1, se_log2, se_log3))
      
      new_dat$CPUE = index.data
      
      
      #for DFA data
      # I_DFA = unlist(DFATrends[i,] - min(DFATrends[i,]) + 0.1 )
      CV_DFA = unlist(DFATrendsSE[i, ]) / unlist(DFATrends[i,])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I2), max(I2))) # BASE
      I_DFA = unlist(DFATrends[i,] + (-2 * min(DFATrends[i,])) ) # TRY THIS!!! MUCH BETTER!
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      
      
      
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
      new_DFA_dat$CPUE = index_DFA.data
      
      
      ### Length Comp ####
      if(length == TRUE){
        
        load(paste("LFreqList_",i,".RData",sep=""))
        
        # From commercial catch
        lencomp_C=LFreqList[[1]]
        lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,2:ncol(lencomp_C)]))[1:11,]
        Lfreq_Comm = cbind(year=1993:2003, seas = (rep(1, 2003-1992)), FltSry = (rep(1, 2003-1992)), 
                           Sex = (rep(0, 2003-1992)), Part = (rep(0, 2003-1992)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        
        # From Surveys
        lencomp_I1=LFreqList[[2]]
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        
        lencomp_I2=LFreqList[[3]]
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        
        lencomp_I3 = LFreqList[[4]]
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        
        # for DFA assessment
        FL = FLoadings[i,1:3]
        lambda = data.frame(lambda=ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda$lambda[1])*LF_I1) + round((as.numeric(lambda$lambda[2])*LF_I2)) + round(as.numeric(lambda$lambda[3])*LF_I3)
        # LF_DFA = LF_I1+LF_I2+LF_I3
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
     
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      # SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
  ### 4 survey indices
  if(Nindices == 4){
   
    Cy = ResultsList$Cy
    Iy1 = ResultsList$Iy1
    Iy2 = ResultsList$Iy2
    Iy3 = ResultsList$Iy3
    Iy4 = ResultsList$Iy4
    Iy1_CV = ResultsList$Iy1_CV
    Iy2_CV = ResultsList$Iy2_CV
    Iy3_CV = ResultsList$Iy3_CV
    Iy4_CV = ResultsList$Iy4_CV
    DFATrends = DFA_Results$DFATrends
    DFATrendsSE = DFA_Results$DFATrendsSE
    FLoadings = DFA_Results$FLoadings
    
    
    for(i in 1:Niters){
      
      # catch
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[40:65], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      I4 = unlist(Iy4[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV = unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      I4_CV = unlist(Iy4_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      se_log4 = sqrt(log(1+(I4_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 4), seas=rep(1, (2018-1992)*4), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992), rep(5, 2018-1992)), 
                              obs = c(I1, I2, I3, I4), se_log = c(se_log1, se_log2, se_log3, se_log4))
      
      new_dat$CPUE = index.data
      new_dat$N_cpue = nrow(index.data)
      new_dat$Nsurveys = Nindices
      new_dat$fleetnames = c("Comm_Fish", "Survey1", "Survey2", "Survey3", "Survey4")
      new_dat$surveytiming = c(-1, 0.5, 0.5, 0.5, 0.5)
      new_dat$areas = c(rep(1,5))
      
      
      #for DFA data
      CV_DFA = unlist(DFATrendsSE[i, ]) / unlist(DFATrends[i,])
      SE_DFA = sqrt(log(1+(CV_DFA^2)))
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I2), max(I2)))
      I_DFA = unlist(DFATrends[i,] + (-2 * min(DFATrends[i,])) ) # TRY THIS!!! MUCH BETTER!
      # I_DFA = unlist(DFATrends[i,41:66])
      # SE_DFA = unlist(DFATrendsSE[i, 41:66]) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # # #for DFA data
      # DFATrends = read.csv("DFATrends.csv")
      # DFATrendsSE = read.csv("DFATrendsSE.csv")
      # # I_DFA = unlist(DFATrends[i,] - min(DFATrends[i,]) + 0.1 )
      # CV_DFA = DFATrendsSE[i, ] / DFATrends[i,]
      # SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I2), max(I2)))
      # # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)
      new_DFA_dat$CPUE = index_DFA.data
      new_DFA_dat$N_cpue = nrow(index_DFA.data)
      
      
      ### Length Comp ####
      if(length == TRUE){
        
        load(paste("LFreqList_",i,".RData",sep=""))
        
        # From commercial catch
        lencomp_C=LFreqList[[1]]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,2:ncol(lencomp_C)]))[1:11,]
        Lfreq_Comm = cbind(year=1993:2003, seas = (rep(1, 2003-1992)), FltSry = (rep(1, 2003-1992)), 
                           Sex = (rep(0, 2003-1992)), Part = (rep(0, 2003-1992)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        
        # From Surveys
        lencomp_I1=LFreqList[[2]]
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        
        lencomp_I2=LFreqList[[3]]
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        
        lencomp_I3 = LFreqList[[4]]
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        
        lencomp_I4 = LFreqList[[5]]
        colnames(lencomp_I4) = NULL
        LF_I4 = as.data.frame(t(lencomp_I4[,2:ncol(lencomp_I4)]))
        Lfreq_I4 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(5, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I4, 1, sum)), LF_I4)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3, Lfreq_I4)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        # for DFA assessment
        FL = FLoadings[i,1:4]
        lambda = data.frame(lambda=ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda$lambda[1])*LF_I1) + round((as.numeric(lambda$lambda[2])*LF_I2)) + round(as.numeric(lambda$lambda[3])*LF_I3) + round(as.numeric(lambda$lambda[4])*LF_I4)
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE

      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      # SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
} # end datfile function

datfileF4 = function(dir = dir, trial = "Base/Trial1", base_dat = base_dat, base_DFA_dat = base_DFA_dat, 
                     Niters = 100, length = TRUE, Nindices=3, EqCatch=NULL, ...){
  
  setwd(paste(dir,trial, sep=""))
  
  new_dat = base_dat
  new_DFA_dat = base_DFA_dat
  
  load("ResultsList.RData")
  load("DFA_Results.RData")
  
  
  if(Nindices ==3){
    
    Cy = ResultsList$Cy
    Iy1 = ResultsList$Iy1
    Iy2 = ResultsList$Iy2
    Iy3 = ResultsList$Iy3
    Iy1_CV = ResultsList$Iy1_CV
    Iy2_CV = ResultsList$Iy2_CV
    Iy3_CV = ResultsList$Iy3_CV
    DFATrends = DFA_Results$DFATrends
    DFATrendsSE = DFA_Results$DFATrendsSE
    FLoadings = DFA_Results$FLoadings
    
    
    for(i in 1:Niters){
      
      # catch
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[40:65], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV= unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 3), seas=rep(1, (2018-1992)*3), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992)), 
                              obs = c(I1, I2, I3), se_log = c(se_log1, se_log2, se_log3))
      
      new_dat$CPUE = index.data
      # par(mfrow=c(3,1))
      # s1=index.data[index.data$index==2,]
      # s2=index.data[index.data$index==3,]
      # s3=index.data[index.data$index==4,]
      # b1=base_dat$CPUE[base_dat$CPUE$index==2,]
      # b2=base_dat$CPUE[base_dat$CPUE$index==3,]
      # b3=base_dat$CPUE[base_dat$CPUE$index==4,]
      # plot(s1$year, s1$obs, type='l', lwd=2)
      # lines(b1$year, b1$obs, type='l', lwd=2, col='red')
      # plot(s2$year, s2$obs, type='l', lwd=2)
      # lines(b2$year, b2$obs, type='l', lwd=2, col='red')
      # plot(s3$year, s3$obs, type='l', lwd=2)
      # lines(b3$year, b3$obs, type='l', lwd=2, col='red')
      
      
      #for DFA data
      # I_DFA = unlist(DFATrends[i,] - min(DFATrends[i,]) + 0.1 )
      CV_DFA = unlist(DFATrendsSE[i, ]) / unlist(DFATrends[i,])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I2), max(I2)))
      I_DFA = unlist(DFATrends[i,] + (-2 * min(DFATrends[i,])) ) # TRY THIS!!! MUCH BETTER!
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      
      
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
      new_DFA_dat$CPUE = index_DFA.data
      
      
      ### Length Comp ####
      if(length == TRUE){
        # From commercial catch
        load(paste("LFreqList_",i,".RData",sep=""))
        lencomp_C=LFreqList[[1]]
        # lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        # lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,3:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=1994:2018, seas = (rep(1, 2018-1993)), FltSry = (rep(1, 2018-1993)), 
                           Sex = (rep(0, 2018-1993)), Part = (rep(0, 2018-1993)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        
        # From Surveys
        lencomp_I1=LFreqList[[2]]
        # lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        
        lencomp_I2 = LFreqList[[3]]
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        
        lencomp_I3 = LFreqList[[4]]
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        
        # for DFA assessment
        FL = FLoadings[i,1:3]
        lambda = data.frame(lambda=ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda$lambda[1])*LF_I1) + round((as.numeric(lambda$lambda[2])*LF_I2)) + round(as.numeric(lambda$lambda[3])*LF_I3)
        # LF_DFA = LF_I1+LF_I2+LF_I3
        
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      # SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
  ### 4 survey indices
  if(Nindices == 4){

    Cy = ResultsList$Cy
    Iy1 = ResultsList$Iy1
    Iy2 = ResultsList$Iy2
    Iy3 = ResultsList$Iy3
    Iy4 = ResultsList$Iy4
    Iy1_CV = ResultsList$Iy1_CV
    Iy2_CV = ResultsList$Iy2_CV
    Iy3_CV = ResultsList$Iy3_CV
    Iy4_CV = ResultsList$Iy4_CV
    DFATrends = DFA_Results$DFATrends
    DFATrendsSE = DFA_Results$DFATrendsSE
    FLoadings = DFA_Results$FLoadings
    
    
    for(i in 1:Niters){
      
      # catch
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[40:65], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      I4 = unlist(Iy4[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV = unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      I4_CV = unlist(Iy4_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      se_log4 = sqrt(log(1+(I4_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 4), seas=rep(1, (2018-1992)*4), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992), rep(5, 2018-1992)), 
                              obs = c(I1, I2, I3, I4), se_log = c(se_log1, se_log2, se_log3, se_log4))
      
      new_dat$CPUE = index.data
      new_dat$N_cpue = nrow(index.data)
      new_dat$Nsurveys = Nindices
      new_dat$fleetnames = c("Comm_Fish", "Survey1", "Survey2", "Survey3", "Survey4")
      new_dat$surveytiming = c(-1, 0.5, 0.5, 0.5, 0.5)
      new_dat$areas = c(rep(1,5))
      
      
      #for DFA data
      CV_DFA = unlist(DFATrendsSE[i, ]) / unlist(DFATrends[i,])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) 
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I2), max(I2)))
      I_DFA = unlist(DFATrends[i,] + (-2 * min(DFATrends[i,])) ) # TRY THIS!!! MUCH BETTER!
      
      # I_DFA = unlist(DFATrends[i,41:66])
      # SE_DFA = unlist(DFATrendsSE[i, 41:66]) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # #for DFA data
      # DFATrends = read.csv("DFATrends.csv")
      # DFATrendsSE = read.csv("DFATrendsSE.csv")
      # # I_DFA = unlist(DFATrends[i,] - min(DFATrends[i,]) + 0.1 )
      # CV_DFA = DFATrendsSE[i, ] / DFATrends[i,]
      # SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I2), max(I2)))
      # # I_DFA = rescale(unlist(DFATrends[i,]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)
      new_DFA_dat$CPUE = index_DFA.data
      new_DFA_dat$N_cpue = nrow(index_DFA.data)
      
      
      ### Length Comp ####
      if(length == TRUE){
        
        load(paste("LFreqList_",i,".RData",sep=""))
        
        # From commercial catch
        lencomp_C = LFreqList[[1]]
        # lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        # lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,3:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=1994:2018, seas = (rep(1, 2018-1993)), FltSry = (rep(1, 2018-1993)), 
                           Sex = (rep(0, 2018-1993)), Part = (rep(0, 2018-1993)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        
        # From Surveys
        lencomp_I1 = LFreqList[[2]]
        # lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        
        
        lencomp_I2 = LFreqList[[3]]
        # lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        
        
        lencomp_I3 = LFreqList[[4]]
        # lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        
        
        lencomp_I4 = LFreqList[[5]]
        # lencomp_I4 = read.csv(paste("LF_Iay4_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I4) = NULL
        LF_I4 = as.data.frame(t(lencomp_I4[,2:ncol(lencomp_I4)]))
        Lfreq_I4 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(5, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I4, 1, sum)), LF_I4)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3, Lfreq_I4)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        # for DFA assessment
        FL = FLoadings[i,1:4]
        lambda = data.frame(lambda=ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda$lambda[1])*LF_I1) + round((as.numeric(lambda$lambda[2])*LF_I2)) + round(as.numeric(lambda$lambda[3])*LF_I3) + round(as.numeric(lambda$lambda[4])*LF_I4)
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      # SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices==4)
  
  
  
} # end datfile function

########## Write data files for conflicting indices only ##############
# NOTE: save for dfa data file is commented out
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/", trial="Base/Trial22", 
          base_dat=base_dat_22, base_DFA_dat=base_DFA_dat_22, Niters=100,  length=TRUE, Nindices=3)



####################################################################################################################################
# ### WRITE FUNCTION #### OLD FORMAT
####################################################################################################################################

datfile = function(dir = dir, trial = "Base/Trial1", base_dat = base_dat, base_DFA_dat = base_DFA_dat, 
                   Niters = 100, age = FALSE, length = TRUE, Nindices=3, EqCatch=NULL, ...){
  
  setwd(paste(dir,trial, sep=""))
  
  new_dat = base_dat
  new_DFA_dat = base_DFA_dat
  
  if(Nindices ==3){
    Cy = read.csv("Cy.csv")
    Iy1 = read.csv("Iy1.csv")
    Iy2 = read.csv("Iy2.csv")
    Iy3 = read.csv("Iy3.csv")
    Iy1_CV = read.csv("Iy1_CV.csv")
    Iy2_CV = read.csv("Iy2_CV.csv")
    Iy3_CV = read.csv("Iy3_CV.csv")
    DFATrends = read.csv("DFATrends.csv")
    DFATrendsSE = read.csv("DFATrendsSE.csv")
    FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
    
    
    for(i in 1:Niters){
 
      # catch
      Cy = read.csv("Cy.csv")
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      Iy1 = read.csv("Iy1.csv")
      Iy2 = read.csv("Iy2.csv")
      Iy3 = read.csv("Iy3.csv")
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      # CV
      Iy1_CV = read.csv("Iy1_CV.csv")
      Iy2_CV = read.csv("Iy2_CV.csv")
      Iy3_CV = read.csv("Iy3_CV.csv")
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV= unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 3), seas=rep(1, (2018-1992)*3), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992)), 
                              obs = c(I1, I2, I3), se_log = c(se_log1, se_log2, se_log3))
      
      new_dat$CPUE = index.data
      
      
      #for DFA data
      DFATrends = read.csv("DFATrends.csv")
      DFATrendsSE = read.csv("DFATrendsSE.csv")
      # I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] - min(DFATrends[i,2:ncol(DFATrends)]) + 0.1 )
      CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # mean_raw = mean(c(I1, I2, I3))
      # min_raw = min(c(I1, I2, I3))
      # trend = unlist(DFATrends[i,2:ncol(DFATrends)])
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2))) # BASE
      I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + (-2 * min(DFATrends[i,2:ncol(DFATrends)])) ) # TRY THIS!!! MUCH BETTER!
      
      # mean2 = mean(I_DFA)
      # I_DFA2 = I_DFA + (mean_raw - mean2)
      # I_DFA = I_DFA2
      #mean(I_DFA2)
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      
      
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
      new_DFA_dat$CPUE = index_DFA.data
      
      
      ### Length Comp ####
      if(length == TRUE){
        # From commercial catch
        lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,2:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(1, 2018-1992)), 
                           Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        # From Surveys
        lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)

        
        
        # for DFA assessment
        FL = FLoadings[i,1:3]
        lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round((as.numeric(lambda[2])*LF_I2)) + round(as.numeric(lambda[3])*LF_I3)
        # LF_DFA = LF_I1+LF_I2+LF_I3
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      
      ### Age Comp? ####
      if(age == TRUE){
        
      } # end if age == TRUE
      
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      # SS_writedat(new_DFA_dat, outfile = "N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/Base/Trial10a/SS_DFA/SN_SS_1/SN.dat", overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices==3)
  
  
  
  ### 4 survey indices
  if(Nindices == 4){
    Cy = read.csv("Cy.csv")
    Iy1 = read.csv("Iy1.csv")
    Iy2 = read.csv("Iy2.csv")
    Iy3 = read.csv("Iy3.csv")
    Iy4 = read.csv("Iy4.csv")
    Iy1_CV = read.csv("Iy1_CV.csv")
    Iy2_CV = read.csv("Iy2_CV.csv")
    Iy3_CV = read.csv("Iy3_CV.csv")
    Iy4_CV = read.csv("Iy4_CV.csv")
    DFATrends = read.csv("DFATrends.csv")
    DFATrendsSE = read.csv("DFATrendsSE.csv")
    FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
    
    for(i in 1:Niters){
      
      # catch
      Cy = read.csv("Cy.csv")
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      Iy1 = read.csv("Iy1.csv")
      Iy2 = read.csv("Iy2.csv")
      Iy3 = read.csv("Iy3.csv")
      Iy4 = read.csv("Iy4.csv")
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      I4 = unlist(Iy4[i,41:66])
      # CV
      Iy1_CV = read.csv("Iy1_CV.csv")
      Iy2_CV = read.csv("Iy2_CV.csv")
      Iy3_CV = read.csv("Iy3_CV.csv")
      Iy4_CV = read.csv("Iy4_CV.csv")
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV = unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      I4_CV = unlist(Iy4_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      se_log4 = sqrt(log(1+(I4_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 4), seas=rep(1, (2018-1992)*4), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992), rep(5, 2018-1992)), 
                              obs = c(I1, I2, I3, I4), se_log = c(se_log1, se_log2, se_log3, se_log4))
      
      new_dat$CPUE = index.data
      new_dat$N_cpue = nrow(index.data)
      new_dat$Nsurveys = Nindices
      new_dat$fleetnames = c("Comm_Fish", "Survey1", "Survey2", "Survey3", "Survey4")
      new_dat$surveytiming = c(-1, 0.5, 0.5, 0.5, 0.5)
      new_dat$areas = c(rep(1,5))
      
      #for DFA data
      DFATrends = read.csv("DFATrends.csv")
      DFATrendsSE = read.csv("DFATrendsSE.csv")
      # I_DFA = unlist(DFATrends[i,41:66])
      if(DFATrends[i,2] > 1e-100) {#unlist(DFATrends[i,2]) < 1e-100
        CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
        SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
        # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
        I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + (-2 * min(DFATrends[i,2:ncol(DFATrends)])) ) # TRY THIS!!! MUCH BETTER!
      }
      if(DFATrends[i,2] < 1e-100) {#unlist(DFATrends[i,2]) < 1e-100
        CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / 1
        SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
        # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
        I_DFA = rep(1, ncol(DFATrends)-1) 
      }
      # SE_DFA = unlist(DFATrendsSE[i, 41:66]) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # #for DFA data
      # DFATrends = read.csv("DFATrends.csv")
      # DFATrendsSE = read.csv("DFATrendsSE.csv")
      # # I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] - min(DFATrends[i,2:ncol(DFATrends)]) + 0.1 )
      # CV_DFA = DFATrendsSE[i, 2:ncol(DFATrendsSE)] / DFATrends[i,2:ncol(DFATrends)]
      # SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
      # # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
      new_DFA_dat$CPUE = index_DFA.data
      new_DFA_dat$N_cpue = nrow(index_DFA.data)
      
      
      ### Length Comp ####
      if(length == TRUE){
        # From commercial catch
        lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,2:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(1, 2018-1992)), 
                           Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        # From Surveys
        lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        lencomp_I4 = read.csv(paste("LF_Iay4_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I4) = NULL
        LF_I4 = as.data.frame(t(lencomp_I4[,2:ncol(lencomp_I4)]))
        Lfreq_I4 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(5, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I4, 1, sum)), LF_I4)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3, Lfreq_I4)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        # for DFA assessment
        FL = FLoadings[i,1:4]
        lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round((as.numeric(lambda[2])*LF_I2)) + round(as.numeric(lambda[3])*LF_I3) + round(as.numeric(lambda[4])*LF_I4)
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      
      ### Age Comp? ####
      if(age == TRUE){
        
      } # end if age == TRUE
      
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
} # end datfile function

datfileF2 = function(dir = dir, trial = "Base/Trial1", base_dat = base_dat, base_DFA_dat = base_DFA_dat, 
                     Niters = 100, age = FALSE, length = TRUE, Nindices=3, EqCatch=NULL, ...){
  
  
  setwd(paste(dir,trial, sep=""))
  
  new_dat = base_dat
  new_DFA_dat = base_DFA_dat
  
  
  
  if(Nindices ==3){
    Cy = read.csv("Cy.csv")
    Iy1 = read.csv("Iy1.csv")
    Iy2 = read.csv("Iy2.csv")
    Iy3 = read.csv("Iy3.csv")
    Iy1_CV = read.csv("Iy1_CV.csv")
    Iy2_CV = read.csv("Iy2_CV.csv")
    Iy3_CV = read.csv("Iy3_CV.csv")
    DFATrends = read.csv("DFATrends.csv")
    DFATrendsSE = read.csv("DFATrendsSE.csv")
    FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
    
    
    for(i in 1:Niters){
      
      # catch
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV= unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 3), seas=rep(1, (2018-1992)*3), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992)), 
                              obs = c(I1, I2, I3), se_log = c(se_log1, se_log2, se_log3))
      
      new_dat$CPUE = index.data
      
      
      #for DFA data
      # I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] - min(DFATrends[i,2:ncol(DFATrends)]) + 0.1 )
      CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
      # SE_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) 
      # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # mean_raw = mean(c(I1, I2, I3))
      # min_raw = min(c(I1, I2, I3))
      # trend = unlist(DFATrends[i,2:ncol(DFATrends)])
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2))) ## BASE USE 
      I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + -2* min(DFATrends[i,2:ncol(DFATrends)])  ) # THIS ONE!
      # I_DFA = I2
      # SE_DFA = se_log2
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(0,1)) + 0.1
      # mean2 = mean(I_DFA)
      # I_DFA2 = trend + max(c(I1, I2, I3))
      # I_DFA2 = trend - (min(trend) - min_raw)
      # I_DFA2 = I_DFA + (mean_raw - mean2)
      # I_DFA = I_DFA2
      # mean(I_DFA2)
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2)+(mean_raw - mean2), max(I2)))
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      
      
      
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)
      new_DFA_dat$CPUE = index_DFA.data
      
      
      ### Length Comp ####
      if(length == TRUE){
        # From commercial catch
        lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,13:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=2004:2018, seas = (rep(1, 2018-2003)), FltSry = (rep(1, 2018-2003)), 
                           Sex = (rep(0, 2018-2003)), Part = (rep(0, 2018-2003)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        # From Surveys
        lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3)
        # tail(lencomp.dat, 20)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        # for DFA assessment
        FL = FLoadings[i,1:3]
        # FL = as.vector(c(1,1,1))
        # lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
        lambda = as.vector(ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round((as.numeric(lambda[2])*LF_I2)) + round(as.numeric(lambda[3])*LF_I3)
        # LF_DFA = (LF_I1) + (LF_I2) + (LF_I3)
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        
      } # end if length == TRUE
      
      
      
      ### Age Comp? ####
      if(age == TRUE){
        
      } # end if age == TRUE
      
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      # SS_writedat(new_DFA_dat, outfile = "N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/Base/Trial10h/SS_DFA/SN_SS_1/SN.dat", overwrite=TRUE)
      
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
  ### 4 survey indices
  if(Nindices == 4){
    Cy = read.csv("Cy.csv")
    Iy1 = read.csv("Iy1.csv")
    Iy2 = read.csv("Iy2.csv")
    Iy3 = read.csv("Iy3.csv")
    Iy4 = read.csv("Iy4.csv")
    Iy1_CV = read.csv("Iy1_CV.csv")
    Iy2_CV = read.csv("Iy2_CV.csv")
    Iy3_CV = read.csv("Iy3_CV.csv")
    Iy4_CV = read.csv("Iy4_CV.csv")
    DFATrends = read.csv("DFATrends.csv")
    DFATrendsSE = read.csv("DFATrendsSE.csv")
    FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
    
    
    for(i in 1:Niters){
      
      # catch
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      I4 = unlist(Iy4[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV = unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      I4_CV = unlist(Iy4_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      se_log4 = sqrt(log(1+(I4_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 4), seas=rep(1, (2018-1992)*4), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992), rep(5, 2018-1992)), 
                              obs = c(I1, I2, I3, I4), se_log = c(se_log1, se_log2, se_log3, se_log4))
      
      new_dat$CPUE = index.data
      new_dat$N_cpue = nrow(index.data)
      new_dat$Nsurveys = Nindices
      new_dat$fleetnames = c("Comm_Fish", "Survey1", "Survey2", "Survey3", "Survey4")
      new_dat$surveytiming = c(-1, 0.5, 0.5, 0.5, 0.5)
      new_dat$areas = c(rep(1,5))
      
      
      #for DFA data
      CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) 
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
      I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + (-2 * min(DFATrends[i,2:ncol(DFATrends)])) ) # TRY THIS!!! MUCH BETTER!
      
      # I_DFA = unlist(DFATrends[i,41:66])
      # SE_DFA = unlist(DFATrendsSE[i, 41:66]) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      
      
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
      new_DFA_dat$CPUE = index_DFA.data
      
      
      ### Length Comp ####
      if(length == TRUE){
        # From commercial catch
        lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,13:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=2004:2018, seas = (rep(1, 2018-2003)), FltSry = (rep(1, 2018-2003)), 
                           Sex = (rep(0, 2018-2003)), Part = (rep(0, 2018-2003)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        # From Surveys
        lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        lencomp_I4 = read.csv(paste("LF_Iay4_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I4) = NULL
        LF_I4 = as.data.frame(t(lencomp_I4[,2:ncol(lencomp_I4)]))
        Lfreq_I4 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(5, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I4, 1, sum)), LF_I4)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3, Lfreq_I4)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        # for DFA assessment
        FL = FLoadings[i,1:4]
        lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
        LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round(as.numeric(lambda[2])*LF_I2) + 
          round(as.numeric(lambda[3])*LF_I3) + round(as.numeric(lambda[4])*LF_I4)
        # LF_DFA = LF_I1+LF_I2+LF_I3+LF_I4
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      
      ### Age Comp? ####
      if(age == TRUE){
        
      } # end if age == TRUE
      
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
} # end datfile function

datfileF3 = function(dir = dir, trial = "Base/Trial1", base_dat = base_dat, base_DFA_dat = base_DFA_dat, 
                     Niters = 100, age = FALSE, length = TRUE, Nindices=3, EqCatch=NULL, ...){
  
  setwd(paste(dir,trial, sep=""))
  
  new_dat = base_dat
  new_DFA_dat = base_DFA_dat
  
  if(Nindices ==3){
    Cy = read.csv("Cy.csv")
    Iy1 = read.csv("Iy1.csv")
    Iy2 = read.csv("Iy2.csv")
    Iy3 = read.csv("Iy3.csv")
    Iy1_CV = read.csv("Iy1_CV.csv")
    Iy2_CV = read.csv("Iy2_CV.csv")
    Iy3_CV = read.csv("Iy3_CV.csv")
    DFATrends = read.csv("DFATrends.csv")
    DFATrendsSE = read.csv("DFATrendsSE.csv")
    FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
    
    
    for(i in 1:Niters){
      
      # catch
      Cy = read.csv("Cy.csv")
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      Iy1 = read.csv("Iy1.csv")
      Iy2 = read.csv("Iy2.csv")
      Iy3 = read.csv("Iy3.csv")
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      # CV
      Iy1_CV = read.csv("Iy1_CV.csv")
      Iy2_CV = read.csv("Iy2_CV.csv")
      Iy3_CV = read.csv("Iy3_CV.csv")
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV= unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 3), seas=rep(1, (2018-1992)*3), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992)), 
                              obs = c(I1, I2, I3), se_log = c(se_log1, se_log2, se_log3))
      
      new_dat$CPUE = index.data
      
      
      #for DFA data
      DFATrends = read.csv("DFATrends.csv")
      DFATrendsSE = read.csv("DFATrendsSE.csv")
      # I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] - min(DFATrends[i,2:ncol(DFATrends)]) + 0.1 )
      CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2))) # BASE
      I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + (-2 * min(DFATrends[i,2:ncol(DFATrends)])) ) # TRY THIS!!! MUCH BETTER!
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      
      
      
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
      new_DFA_dat$CPUE = index_DFA.data
      
      
      ### Length Comp ####
      if(length == TRUE){
        # From commercial catch
        lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,2:ncol(lencomp_C)]))[1:11,]
        Lfreq_Comm = cbind(year=1993:2003, seas = (rep(1, 2003-1992)), FltSry = (rep(1, 2003-1992)), 
                           Sex = (rep(0, 2003-1992)), Part = (rep(0, 2003-1992)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        # From Surveys
        lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        
        # for DFA assessment
        FL = FLoadings[i,1:3]
        lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round((as.numeric(lambda[2])*LF_I2)) + round(as.numeric(lambda[3])*LF_I3)
        # LF_DFA = LF_I1+LF_I2+LF_I3
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      
      ### Age Comp? ####
      if(age == TRUE){
        
      } # end if age == TRUE
      
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
  ### 4 survey indices
  if(Nindices == 4){
    Cy = read.csv("Cy.csv")
    Iy1 = read.csv("Iy1.csv")
    Iy2 = read.csv("Iy2.csv")
    Iy3 = read.csv("Iy3.csv")
    Iy4 = read.csv("Iy4.csv")
    Iy1_CV = read.csv("Iy1_CV.csv")
    Iy2_CV = read.csv("Iy2_CV.csv")
    Iy3_CV = read.csv("Iy3_CV.csv")
    Iy4_CV = read.csv("Iy4_CV.csv")
    DFATrends = read.csv("DFATrends.csv")
    DFATrendsSE = read.csv("DFATrendsSE.csv")
    FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
    
    for(i in 1:Niters){
      
      # catch
      Cy = read.csv("Cy.csv")
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      I4 = unlist(Iy4[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV = unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      I4_CV = unlist(Iy4_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      se_log4 = sqrt(log(1+(I4_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 4), seas=rep(1, (2018-1992)*4), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992), rep(5, 2018-1992)), 
                              obs = c(I1, I2, I3, I4), se_log = c(se_log1, se_log2, se_log3, se_log4))
      
      new_dat$CPUE = index.data
      new_dat$N_cpue = nrow(index.data)
      new_dat$Nsurveys = Nindices
      new_dat$fleetnames = c("Comm_Fish", "Survey1", "Survey2", "Survey3", "Survey4")
      new_dat$surveytiming = c(-1, 0.5, 0.5, 0.5, 0.5)
      new_dat$areas = c(rep(1,5))
      
      
      #for DFA data
      CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
      SE_DFA = sqrt(log(1+(CV_DFA^2)))
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
      I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + (-2 * min(DFATrends[i,2:ncol(DFATrends)])) ) # TRY THIS!!! MUCH BETTER!
      # I_DFA = unlist(DFATrends[i,41:66])
      # SE_DFA = unlist(DFATrendsSE[i, 41:66]) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # # #for DFA data
      # DFATrends = read.csv("DFATrends.csv")
      # DFATrendsSE = read.csv("DFATrendsSE.csv")
      # # I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] - min(DFATrends[i,2:ncol(DFATrends)]) + 0.1 )
      # CV_DFA = DFATrendsSE[i, 2:ncol(DFATrendsSE)] / DFATrends[i,2:ncol(DFATrends)]
      # SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
      # # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)
      new_DFA_dat$CPUE = index_DFA.data
      new_DFA_dat$N_cpue = nrow(index_DFA.data)
      
      
      ### Length Comp ####
      if(length == TRUE){
        # From commercial catch
        lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,2:ncol(lencomp_C)]))[1:11,]
        Lfreq_Comm = cbind(year=1993:2003, seas = (rep(1, 2003-1992)), FltSry = (rep(1, 2003-1992)), 
                           Sex = (rep(0, 2003-1992)), Part = (rep(0, 2003-1992)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        
        # From Surveys
        lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        lencomp_I4 = read.csv(paste("LF_Iay4_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I4) = NULL
        LF_I4 = as.data.frame(t(lencomp_I4[,2:ncol(lencomp_I4)]))
        Lfreq_I4 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(5, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I4, 1, sum)), LF_I4)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3, Lfreq_I4)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        # for DFA assessment
        FL = FLoadings[i,1:4]
        lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round((as.numeric(lambda[2])*LF_I2)) + round(as.numeric(lambda[3])*LF_I3) + round(as.numeric(lambda[4])*LF_I4)
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      
      ### Age Comp? ####
      if(age == TRUE){
        
      } # end if age == TRUE
      
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
} # end datfile function

datfileF4 = function(dir = dir, trial = "Base/Trial1", base_dat = base_dat, base_DFA_dat = base_DFA_dat, 
                     Niters = 100, age = FALSE, length = TRUE, Nindices=3, EqCatch=NULL, ...){
  
  setwd(paste(dir,trial, sep=""))
  
  new_dat = base_dat
  new_DFA_dat = base_DFA_dat
  
  if(Nindices ==3){
    Cy = read.csv("Cy.csv")
    Iy1 = read.csv("Iy1.csv")
    Iy2 = read.csv("Iy2.csv")
    Iy3 = read.csv("Iy3.csv")
    Iy1_CV = read.csv("Iy1_CV.csv")
    Iy2_CV = read.csv("Iy2_CV.csv")
    Iy3_CV = read.csv("Iy3_CV.csv")
    DFATrends = read.csv("DFATrends.csv")
    DFATrendsSE = read.csv("DFATrendsSE.csv")
    FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
    
    
    for(i in 1:Niters){
      
      # catch
      Cy = read.csv("Cy.csv")
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      Iy1 = read.csv("Iy1.csv")
      Iy2 = read.csv("Iy2.csv")
      Iy3 = read.csv("Iy3.csv")
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      # CV
      Iy1_CV = read.csv("Iy1_CV.csv")
      Iy2_CV = read.csv("Iy2_CV.csv")
      Iy3_CV = read.csv("Iy3_CV.csv")
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV= unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 3), seas=rep(1, (2018-1992)*3), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992)), 
                              obs = c(I1, I2, I3), se_log = c(se_log1, se_log2, se_log3))
      
      new_dat$CPUE = index.data
      
      
      #for DFA data
      DFATrends = read.csv("DFATrends.csv")
      DFATrendsSE = read.csv("DFATrendsSE.csv")
      # I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] - min(DFATrends[i,2:ncol(DFATrends)]) + 0.1 )
      CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
      I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + (-2 * min(DFATrends[i,2:ncol(DFATrends)])) ) # TRY THIS!!! MUCH BETTER!
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      
      
      
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
      new_DFA_dat$CPUE = index_DFA.data
      
      
      ### Length Comp ####
      if(length == TRUE){
        # From commercial catch
        lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,3:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=1994:2018, seas = (rep(1, 2018-1993)), FltSry = (rep(1, 2018-1993)), 
                           Sex = (rep(0, 2018-1993)), Part = (rep(0, 2018-1993)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        # From Surveys
        lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        
        # for DFA assessment
        FL = FLoadings[i,1:3]
        lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round((as.numeric(lambda[2])*LF_I2)) + round(as.numeric(lambda[3])*LF_I3)
        # LF_DFA = LF_I1+LF_I2+LF_I3
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      
      ### Age Comp? ####
      if(age == TRUE){
        
      } # end if age == TRUE
      
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
  ### 4 survey indices
  if(Nindices == 4){
    Cy = read.csv("Cy.csv")
    Iy1 = read.csv("Iy1.csv")
    Iy2 = read.csv("Iy2.csv")
    Iy3 = read.csv("Iy3.csv")
    Iy4 = read.csv("Iy4.csv")
    Iy1_CV = read.csv("Iy1_CV.csv")
    Iy2_CV = read.csv("Iy2_CV.csv")
    Iy3_CV = read.csv("Iy3_CV.csv")
    Iy4_CV = read.csv("Iy4_CV.csv")
    DFATrends = read.csv("DFATrends.csv")
    DFATrendsSE = read.csv("DFATrendsSE.csv")
    FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
    
    for(i in 1:Niters){
      
      # catch
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      I4 = unlist(Iy4[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV = unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      I4_CV = unlist(Iy4_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      se_log4 = sqrt(log(1+(I4_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 4), seas=rep(1, (2018-1992)*4), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992), rep(5, 2018-1992)), 
                              obs = c(I1, I2, I3, I4), se_log = c(se_log1, se_log2, se_log3, se_log4))
      
      new_dat$CPUE = index.data
      new_dat$N_cpue = nrow(index.data)
      new_dat$Nsurveys = Nindices
      new_dat$fleetnames = c("Comm_Fish", "Survey1", "Survey2", "Survey3", "Survey4")
      new_dat$surveytiming = c(-1, 0.5, 0.5, 0.5, 0.5)
      new_dat$areas = c(rep(1,5))
      
      
      #for DFA data
      CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) 
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
      I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + (-2 * min(DFATrends[i,2:ncol(DFATrends)])) ) # TRY THIS!!! MUCH BETTER!
      
      # I_DFA = unlist(DFATrends[i,41:66])
      # SE_DFA = unlist(DFATrendsSE[i, 41:66]) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # #for DFA data
      # DFATrends = read.csv("DFATrends.csv")
      # DFATrendsSE = read.csv("DFATrendsSE.csv")
      # # I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] - min(DFATrends[i,2:ncol(DFATrends)]) + 0.1 )
      # CV_DFA = DFATrendsSE[i, 2:ncol(DFATrendsSE)] / DFATrends[i,2:ncol(DFATrends)]
      # SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
      # # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)
      new_DFA_dat$CPUE = index_DFA.data
      new_DFA_dat$N_cpue = nrow(index_DFA.data)
      
      
      ### Length Comp ####
      if(length == TRUE){
        # From commercial catch
        lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,3:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=1994:2018, seas = (rep(1, 2018-1993)), FltSry = (rep(1, 2018-1993)), 
                           Sex = (rep(0, 2018-1993)), Part = (rep(0, 2018-1993)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        
        # From Surveys
        lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        lencomp_I4 = read.csv(paste("LF_Iay4_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I4) = NULL
        LF_I4 = as.data.frame(t(lencomp_I4[,2:ncol(lencomp_I4)]))
        Lfreq_I4 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(5, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I4, 1, sum)), LF_I4)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3, Lfreq_I4)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        # for DFA assessment
        FL = FLoadings[i,1:4]
        lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round((as.numeric(lambda[2])*LF_I2)) + round(as.numeric(lambda[3])*LF_I3) + round(as.numeric(lambda[4])*LF_I4)
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      
      ### Age Comp? ####
      if(age == TRUE){
        
      } # end if age == TRUE
      
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
} # end datfile function



### FIT BASE ####
base_dat$CPUE

datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial1a", 
        base_dat=base_dat, base_DFA_dat=base_DFA_dat, Niters=100, age=FALSE, length=TRUE, Nindices=3)

datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial1", 
        base_dat=base_dat, base_DFA_dat=base_DFA_dat, Niters=10, age=FALSE, length=TRUE, Nindices=3)


base_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")
# base_dat_1 = SS_readdat(file="Y:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
# base_DFA_dat_1 = SS_readdat(file="Y:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")

datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial1", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial2", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial3", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial4", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial5", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial6", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)


datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Gradual/Trial4", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Gradual/Trial5", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Gradual/Trial6", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)


base_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial7", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial8", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial9", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial10", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial11", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial12", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)

datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Gradual/Trial10", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Gradual/Trial11", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Gradual/Trial12", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)



base_dat_13 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial13\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_13 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial13\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial13", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial14", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial15", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial16", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial17", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial18", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=3)

datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Gradual/Trial16", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Gradual/Trial17", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Gradual/Trial18", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=3)




base_dat_19 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial19\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_19 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial19\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial19", 
          base_dat=base_dat_19, base_DFA_dat=base_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial20", 
          base_dat=base_dat_19, base_DFA_dat=base_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial21", 
          base_dat=base_dat_19, base_DFA_dat=base_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial22", 
          base_dat=base_dat_19, base_DFA_dat=base_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial23", 
          base_dat=base_dat_19, base_DFA_dat=base_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial24", 
          base_dat=base_dat_19, base_DFA_dat=base_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=3)


datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Gradual/Trial22", 
          base_dat=baseG_dat_19, base_DFA_dat=baseG_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Gradual/Trial23", 
          base_dat=baseG_dat_19, base_DFA_dat=baseG_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Gradual/Trial24", 
          base_dat=base_dat_19, base_DFA_dat=base_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=3)


## FOUR INDICES ####
Four_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial1\\SS\\SN_SS_1\\SN.dat")
Four_DFA_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")

datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial1", 
        base_dat=Four_dat_1, base_DFA_dat=Four_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial2", 
        base_dat=Four_dat_1, base_DFA_dat=Four_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial3", 
        base_dat=Four_dat_1, base_DFA_dat=Four_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial4", 
        base_dat=Four_dat_1, base_DFA_dat=Four_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial5", 
        base_dat=Four_dat_1, base_DFA_dat=Four_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial6", 
        base_dat=Four_dat_1, base_DFA_dat=Four_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=4)

datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndicesGradual/Trial4", 
        base_dat=Four_dat_1, base_DFA_dat=Four_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndicesGradual/Trial5", 
        base_dat=Four_dat_1, base_DFA_dat=Four_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndicesGradual/Trial6", 
        base_dat=Four_dat_1, base_DFA_dat=Four_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=4)



Four_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial7\\SS\\SN_SS_1\\SN.dat")
Four_DFA_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial7\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial7", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial8", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial9", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial10", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial11", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial12", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)

datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndicesGradual/Trial10", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndicesGradual/Trial11", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndicesGradual/Trial12", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)



Four_dat_13 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial13\\SS\\SN_SS_1\\SN.dat")
Four_DFA_dat_13 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial13\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial13", 
          base_dat=Four_dat_13, base_DFA_dat=Four_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial14", 
          base_dat=Four_dat_13, base_DFA_dat=Four_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial15", 
          base_dat=Four_dat_13, base_DFA_dat=Four_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial16", 
          base_dat=Four_dat_13, base_DFA_dat=Four_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial17", 
          base_dat=Four_dat_13, base_DFA_dat=Four_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial18", 
          base_dat=Four_dat_13, base_DFA_dat=Four_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=4)

datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndicesGradual/Trial16", 
          base_dat=Four_dat_13, base_DFA_dat=Four_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndicesGradual/Trial17", 
          base_dat=Four_dat_13, base_DFA_dat=Four_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndicesGradual/Trial18", 
          base_dat=Four_dat_13, base_DFA_dat=Four_DFA_dat_13, Niters=100, age=FALSE, length=TRUE, Nindices=4)



Four_dat_19 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial19\\SS\\SN_SS_1\\SN.dat")
Four_DFA_dat_19 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial19\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial19", 
          base_dat=Four_dat_19, base_DFA_dat=Four_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial20", 
          base_dat=Four_dat_19, base_DFA_dat=Four_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial21", 
          base_dat=Four_dat_19, base_DFA_dat=Four_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial22", 
          base_dat=Four_dat_19, base_DFA_dat=Four_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial23", 
          base_dat=Four_dat_19, base_DFA_dat=Four_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial24", 
          base_dat=Four_dat_19, base_DFA_dat=Four_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=4)

datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndicesGradual/Trial22", 
          base_dat=Four_dat_19, base_DFA_dat=Four_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndicesGradual/Trial23", 
          base_dat=Four_dat_19, base_DFA_dat=Four_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndicesGradual/Trial24", 
          base_dat=Four_dat_19, base_DFA_dat=Four_DFA_dat_19, Niters=100, age=FALSE, length=TRUE, Nindices=4)










################# FIT DATA FILES ###########
Four_dat_2 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial2\\SS\\SN_SS_1\\SN.dat")
Four_DFA_dat_2 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial2\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial1", 
          base_dat=Four_dat_2, base_DFA_dat=Four_DFA_dat_2, Niters=100, age=FALSE, length=TRUE, Nindices=4)

Four_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial1\\SS\\SN_SS_1\\SN.dat")
Four_DFA_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial19", 
          base_dat=Four_dat_1, base_DFA_dat=Four_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=4)

base_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")
# base_dat_1 = SS_readdat(file="Y:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
# base_DFA_dat_1 = SS_readdat(file="Y:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")

datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial1", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=3)



base_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial7", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial8", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial9", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial10", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial11", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial12", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=1, age=FALSE, length=TRUE, Nindices=3)


base_dat_19 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial19\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_19 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial19\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial19", 
          base_dat=base_dat_19, base_DFA_dat=base_DFA_dat_19, Niters=1, age=FALSE, length=TRUE, Nindices=3)



Base_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS\\SN_SS_1\\SN.dat")
Base_DFA_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial7", 
          base_dat=Base_dat_7, base_DFA_dat=Base_DFA_dat_7, Niters=1, age=FALSE, length=TRUE, Nindices=4)

Four_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial7\\SS\\SN_SS_1\\SN.dat")
Four_DFA_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial7\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial7", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial8", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial9", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial10", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial11", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial12", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)

datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndicesGradual/Trial10", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndicesGradual/Trial11", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndicesGradual/Trial12", 
          base_dat=Four_dat_7, base_DFA_dat=Four_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=4)



Base_dat_13 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial13\\SS\\SN_SS_1\\SN.dat")
Base_DFA_dat_13 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial13\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial13", 
        base_dat=Base_dat_13, base_DFA_dat=Base_DFA_dat_13, Niters=1, age=FALSE, length=TRUE, Nindices=4)


Base_dat_19 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial19\\SS\\SN_SS_1\\SN.dat")
Base_DFA_dat_19 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial19\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial19", 
          base_dat=Base_dat_19, base_DFA_dat=Base_DFA_dat_19, Niters=1, age=FALSE, length=TRUE, Nindices=4)




######### RUNS ###############


Four_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial1\\SS\\SN_SS_1\\SN.dat")
Four_DFA_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndices\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="FourIndices/Trial13", 
          base_dat=Four_dat_1, base_DFA_dat=Four_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=4)



base_dat_14 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial14\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_14 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial14\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial13", 
          base_dat=base_dat_14, base_DFA_dat=base_DFA_dat_14, Niters=100, age=FALSE, length=TRUE, Nindices=3)


datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial24", 
          base_dat=base_dat_19, base_DFA_dat=base_DFA_dat_19, Niters=1, age=FALSE, length=TRUE, Nindices=3)


base_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial1", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)

base_dat_2 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_2 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial2", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)

base_dat_3 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_3 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial3", 
        base_dat=base_dat_3, base_DFA_dat=base_DFA_dat_3, Niters=100, age=FALSE, length=TRUE, Nindices=3)

base_dat_4 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_4 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial4", 
        base_dat=base_dat_4, base_DFA_dat=base_DFA_dat_4, Niters=100, age=FALSE, length=TRUE, Nindices=3)

base_dat_5 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_5 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial5", 
        base_dat=base_dat_5, base_DFA_dat=base_DFA_dat_5, Niters=100, age=FALSE, length=TRUE, Nindices=3)

base_dat_6 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_6 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial6", 
        base_dat=base_dat_6, base_DFA_dat=base_DFA_dat_6, Niters=100, age=FALSE, length=TRUE, Nindices=3)




datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial19", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial20", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial21", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial22", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial23", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial24", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=3)



base_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial7", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)

datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial8", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)

datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial9", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)

datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial10", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)

datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial11", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)

datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial12", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)





# TEST
base_dat_10 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial10\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_10 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial10\\SS_DFA\\SN_SS_1\\SN.dat")
base_dat = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial13\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial13\\SS_DFA\\SN_SS_1\\SN.dat")
trial = "Base/Trial13"
base_dat = base_dat
base_DFA_dat = base_DFA_dat
Niters = 1 
age = FALSE
length = TRUE 
Nindices=3
dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/"

#### F2 #####
# For iterations BASE (A) 7 - 12  &  GRADUAL BASE (C) 10-12

plot(y = Iy1[1,41:66]/max(Iy1[1,41:66]), x=1993:2018, lwd=2, type='b', ylim=c(0,1.5), col='grey10')
lines(y = Iy2[1,41:66]/max(Iy2[1,41:66]), x=1993:2018, lwd=2, type='b', lty=2, col='grey30')
lines(y = Iy3[1,41:66]/max(Iy3[1,41:66]), x=1993:2018, lwd=2, type='b', lty=3, col='grey50')
lines(y=I_DFA/max(I_DFA), x=1993:2018, lwd=2, type='b', lty=1, col='blue')



base_dat_10 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\BaseTrial10_ModelDevelopment\\Update10.9.2018_1\\Iteration50\\SN.dat")
base_DFA_dat_10 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\BaseTrial10_ModelDevelopment\\Update10.9.2018_1\\Iteration50_DFA\\SN.dat")

base_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS_DFA\\SN_SS_1\\SN.dat")


datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial10", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)


datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial7", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial8", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial9", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial11", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial12", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)



base_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS_DFA\\SN_SS_1\\SN.dat")

datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial7", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)


base_dat_11 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial11\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_11 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial11\\SS_DFA\\SN_SS_1\\SN.dat")

datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial11", 
          base_dat=base_dat_11, base_DFA_dat=base_DFA_dat_11, Niters=100, age=FALSE, length=TRUE, Nindices=3)



base_DFA_dat_10 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial10\\SS_DFA\\SN_SS_1\\SN.dat")


datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial11", 
          base_dat=base_dat_10, base_DFA_dat=base_DFA_dat_10, Niters=100, age=FALSE, length=TRUE, Nindices=3)
# datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial11", 
          # base_dat=base_dat_10, base_DFA_dat=base_DFA_dat_10, Niters=1, age=FALSE, length=TRUE, Nindices=3)


## For Trial 12...
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial12", 
          base_dat=base_dat_11, base_DFA_dat=base_DFA_dat_11, Niters=1, age=FALSE, length=TRUE, Nindices=3)


## For Trial 7...
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial7", 
          base_dat=base_dat_11, base_DFA_dat=base_DFA_dat_11, Niters=100, age=FALSE, length=TRUE, Nindices=3)


## For Trial 8...
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial8", 
          base_dat=base_dat_11, base_DFA_dat=base_DFA_dat_11, Niters=100, age=FALSE, length=TRUE, Nindices=3)




## For Trial 8...
datfileF2(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial9", 
          base_dat=base_dat_11, base_DFA_dat=base_DFA_dat_11, Niters=100, age=FALSE, length=TRUE, Nindices=3)









######################
# For TRIAL 1
# base_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
# base_DFA_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")
base_dat = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS_DFA\\SN_SS_1\\SN.dat")

datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial1", 
          base_dat=base_dat, base_DFA_dat=base_DFA_dat, Niters=1, age=FALSE, length=TRUE, Nindices=3)

base_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")
datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial1", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)

datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial2", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)

datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial3", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)

datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial4", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)

datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial5", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)

datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial6", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=100, age=FALSE, length=TRUE, Nindices=3)





#### F3 #####
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial13", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial14", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial15", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial16", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial17", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF3(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial18", 
        base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=3)




#### F4 ####
base_dat_13 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial13\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_13 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial13\\SS_DFA\\SN_SS_1\\SN.dat")
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial19", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial20", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial21", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial22", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial23", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial24", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=1, age=FALSE, length=TRUE, Nindices=3)

datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Gradual/Trial22", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Gradual/Trial23", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=1, age=FALSE, length=TRUE, Nindices=3)
datfileF4(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Gradual/Trial24", 
          base_dat=base_dat_13, base_DFA_dat=base_DFA_dat_13, Niters=1, age=FALSE, length=TRUE, Nindices=3)


#### NO UNCERTAINTY RUN #####


base_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS_DFA\\SN_SS_1\\SN.dat")

datfileF2(dir="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base_F2_NOUNCERTAINTY", trial="Base/Trial7", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=1, age=FALSE, length=TRUE, Nindices=3)
setwd("N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base_F2_NOUNCERTAINTY")
SS_writedat(new_dat, outfile = paste("SS\\SN.dat", sep=""), overwrite=TRUE)
SS_writedat(new_DFA_dat, outfile = paste("SS_DFA\\SN.dat", sep=""), overwrite=TRUE)

## No Uncertainty to uncertainty
base_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_7 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\SS_DFA\\SN_SS_1\\SN.dat")

datfileF2(dir="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\", trial="Base/Trial7", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)



### USE TRIAL 7 TO POPULATE TRIALS: 8-10
datfileF2(dir="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base_F2_NOUNCERTAINTY", trial="Base/Trial7", 
          base_dat=base_dat_7, base_DFA_dat=base_DFA_dat_7, Niters=100, age=FALSE, length=TRUE, Nindices=3)


################## NO UNCERTAINTY TRIAL 1 ###### 
base_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")

datfileF(dir="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base_F1_NOUNCERTAINTY", trial="Base/Trial1", 
          base_dat=base_dat_1, base_DFA_dat=base_DFA_dat_1, Niters=1, age=FALSE, length=TRUE, Nindices=3)
setwd("N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base_F1_NOUNCERTAINTY")
SS_writedat(new_dat, outfile = paste("SS\\SN.dat", sep=""), overwrite=TRUE)
SS_writedat(new_DFA_dat, outfile = paste("SS_DFA\\SN.dat", sep=""), overwrite=TRUE)





# # For TRIAL 2
# base_dat_2 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial2\\SS\\SN_SS_1\\SN.dat")
# base_DFA_dat_2 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial2\\SS_DFA\\SN_SS_1\\SN.dat")
# 
# datfile(dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/", trial="Base/Trial2", 
#         base_dat=base_dat_2, base_DFA_dat=base_DFA_dat_2, Niters=100, age=FALSE, length=TRUE, Nindices=3)


#####################################################################################################################################
#### WRITE FUNCTION FOR UPDATING FITTED DAT FILES
#####################################################################################################################################

base_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS\\SN_SS_1\\SN.dat")
base_DFA_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\SS_DFA\\SN_SS_1\\SN.dat")


datF1 = function(dir = dir, trial = "Base/Trial1", 
                   Niters = 100, age = FALSE, length = TRUE, Nindices=3, EqCatch=NULL, ...){
  
  
  # for(t in 1:Niters){
    setwd(paste(dir,trial, sep=""))
    
    
    if(Nindices ==3){
      Cy = read.csv("Cy.csv")
      Iy1 = read.csv("Iy1.csv")
      Iy2 = read.csv("Iy2.csv")
      Iy3 = read.csv("Iy3.csv")
      Iy1_CV = read.csv("Iy1_CV.csv")
      Iy2_CV = read.csv("Iy2_CV.csv")
      Iy3_CV = read.csv("Iy3_CV.csv")
      DFATrends = read.csv("DFATrends.csv")
      DFATrendsSE = read.csv("DFATrendsSE.csv")
      FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
      
      
      for(i in 1:Niters){
        
        # setwd(paste(dir,trial, sep=""))
        
        new_dat = SS_readdat(file=paste(dir,trial,"\\SS\\SN_SS_", i, "\\SN.dat", sep="") )
        new_DFA_dat = SS_readdat(file=paste(dir,trial,'\\SS_DFA\\SN_SS_',i,"\\SN.dat", sep="") )
        
        
        
        # catch
        Cy = read.csv("Cy.csv")
        catch = Cy[i,] #100 years for iteration i
        catch = unlist(catch)
        catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
        new_dat$catch = catch.data
        new_DFA_dat$catch = catch.data
        new_dat$N_catch = nrow(catch.data)
        new_DFA_dat$N_catch = nrow(catch.data)
        new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
        new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
        
        
        # CPUE [year, seas, index, obs, se_log]
        # Index
        Iy1 = read.csv("Iy1.csv")
        Iy2 = read.csv("Iy2.csv")
        Iy3 = read.csv("Iy3.csv")
        I1 = unlist(Iy1[i,41:66])
        I2 = unlist(Iy2[i,41:66])
        I3 = unlist(Iy3[i,41:66])
        # CV
        Iy1_CV = read.csv("Iy1_CV.csv")
        Iy2_CV = read.csv("Iy2_CV.csv")
        Iy3_CV = read.csv("Iy3_CV.csv")
        I1_CV = unlist(Iy1_CV[i, 41:66])
        I2_CV= unlist(Iy2_CV[i, 41:66])
        I3_CV = unlist(Iy3_CV[i, 41:66])
        se_log1 = sqrt(log(1+(I1_CV^2)))
        se_log2 = sqrt(log(1+(I2_CV^2)))
        se_log3 = sqrt(log(1+(I3_CV^2)))
        
        index.data = data.frame(year = rep(1993:2018, 3), seas=rep(1, (2018-1992)*3), 
                                index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992)), 
                                obs = c(I1, I2, I3), se_log = c(se_log1, se_log2, se_log3))
        
        new_dat$CPUE = index.data
        
        
        #for DFA data
        DFATrends = read.csv("DFATrends.csv")
        DFATrendsSE = read.csv("DFATrendsSE.csv")
        # I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] - min(DFATrends[i,2:ncol(DFATrends)]) + 0.1 )
        CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
        SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
        # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
        # mean_raw = mean(c(I1, I2, I3))
        # min_raw = min(c(I1, I2, I3))
        # trend = unlist(DFATrends[i,2:ncol(DFATrends)])
        # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2))) # BASE
        I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + (-2 * min(DFATrends[i,2:ncol(DFATrends)])) ) # TRY THIS!!! MUCH BETTER!
        
        # mean2 = mean(I_DFA)
        # I_DFA2 = I_DFA + (mean_raw - mean2)
        # I_DFA = I_DFA2
        #mean(I_DFA2)
        # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
        
        
        index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                    obs = I_DFA, se=SE_DFA)############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
        new_DFA_dat$CPUE = index_DFA.data
        
        
        ### Length Comp ####
        if(length == TRUE){
          # From commercial catch
          lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
          lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
          colnames(lencomp_C) = NULL
          LF_C = as.data.frame(t(lencomp_C[,2:ncol(lencomp_C)]))
          Lfreq_Comm = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(1, 2018-1992)), 
                             Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                             Nsamp = (apply(LF_C, 1, sum)), LF_C)
          # From Surveys
          lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
          colnames(lencomp_I1) = NULL
          LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
          Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                           Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                           Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
          lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
          colnames(lencomp_I2) = NULL
          LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
          Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                           Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                           Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
          lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
          colnames(lencomp_I3) = NULL
          LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
          Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                           Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                           Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
          
          lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3)
          
          new_dat$lencomp = lencomp.dat
          new_dat$N_lencomp = nrow(lencomp.dat)
          
          
          
          # for DFA assessment
          FL = FLoadings[i,1:3]
          lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
          
          LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round((as.numeric(lambda[2])*LF_I2)) + round(as.numeric(lambda[3])*LF_I3)
          # LF_DFA = LF_I1+LF_I2+LF_I3
          Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                            Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                            Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
          lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
          
          new_DFA_dat$lencomp = lencomp_DFA.dat
          new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
          
        } # end if length == TRUE
        
        
        
        
        SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
        SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
        # SS_writedat(new_DFA_dat, outfile = "N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/Base/Trial10a/SS_DFA/SN_SS_1/SN.dat", overwrite=TRUE)
      }# End i loop
      
    } # end if(Nindices==3)
    
    
    
    ### 4 survey indices
    if(Nindices == 4){
      
      Cy = read.csv("Cy.csv")
      Iy1 = read.csv("Iy1.csv")
      Iy2 = read.csv("Iy2.csv")
      Iy3 = read.csv("Iy3.csv")
      Iy4 = read.csv("Iy4.csv")
      Iy1_CV = read.csv("Iy1_CV.csv")
      Iy2_CV = read.csv("Iy2_CV.csv")
      Iy3_CV = read.csv("Iy3_CV.csv")
      Iy4_CV = read.csv("Iy4_CV.csv")
      DFATrends = read.csv("DFATrends.csv")
      DFATrendsSE = read.csv("DFATrendsSE.csv")
      FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
      
      for(i in 1:Niters){
        
        new_dat = SS_readdat(file=paste(dir,trial,"\\SS\\SN_SS_", i, "\\SN.dat", sep="") )
        new_DFA_dat = SS_readdat(file=paste(dir,trial,'\\SS_DFA\\SN_SS_',i,"\\SN.dat", sep="") )
        
        
        # catch
        Cy = read.csv("Cy.csv")
        catch = Cy[i,] #100 years for iteration i
        catch = unlist(catch)
        catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
        new_dat$catch = catch.data
        new_DFA_dat$catch = catch.data
        new_dat$N_catch = nrow(catch.data)
        new_DFA_dat$N_catch = nrow(catch.data)
        new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
        new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
        
        # CPUE [year, seas, index, obs, se_log]
        # Index
        Iy1 = read.csv("Iy1.csv")
        Iy2 = read.csv("Iy2.csv")
        Iy3 = read.csv("Iy3.csv")
        Iy4 = read.csv("Iy4.csv")
        I1 = unlist(Iy1[i,41:66])
        I2 = unlist(Iy2[i,41:66])
        I3 = unlist(Iy3[i,41:66])
        I4 = unlist(Iy4[i,41:66])
        # CV
        Iy1_CV = read.csv("Iy1_CV.csv")
        Iy2_CV = read.csv("Iy2_CV.csv")
        Iy3_CV = read.csv("Iy3_CV.csv")
        Iy4_CV = read.csv("Iy4_CV.csv")
        I1_CV = unlist(Iy1_CV[i, 41:66])
        I2_CV = unlist(Iy2_CV[i, 41:66])
        I3_CV = unlist(Iy3_CV[i, 41:66])
        I4_CV = unlist(Iy4_CV[i, 41:66])
        se_log1 = sqrt(log(1+(I1_CV^2)))
        se_log2 = sqrt(log(1+(I2_CV^2)))
        se_log3 = sqrt(log(1+(I3_CV^2)))
        se_log4 = sqrt(log(1+(I4_CV^2)))
        
        index.data = data.frame(year = rep(1993:2018, 4), seas=rep(1, (2018-1992)*4), 
                                index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992), rep(5, 2018-1992)), 
                                obs = c(I1, I2, I3, I4), se_log = c(se_log1, se_log2, se_log3, se_log4))
        
        new_dat$CPUE = index.data
        new_dat$N_cpue = nrow(index.data)
        new_dat$Nsurveys = Nindices
        new_dat$fleetnames = c("Comm_Fish", "Survey1", "Survey2", "Survey3", "Survey4")
        new_dat$surveytiming = c(-1, 0.5, 0.5, 0.5, 0.5)
        new_dat$areas = c(rep(1,5))
        
        #for DFA data
        DFATrends = read.csv("DFATrends.csv")
        DFATrendsSE = read.csv("DFATrendsSE.csv")
        # I_DFA = unlist(DFATrends[i,41:66])
        if(DFATrends[i,2] > 1e-100) {#unlist(DFATrends[i,2]) < 1e-100
          CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
          SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
          # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
          I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + (-2 * min(DFATrends[i,2:ncol(DFATrends)])) ) # TRY THIS!!! MUCH BETTER!
        }
        if(DFATrends[i,2] < 1e-100) {#unlist(DFATrends[i,2]) < 1e-100
          CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / 1
          SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
          # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
          I_DFA = rep(1, ncol(DFATrends)-1) 
        }
        # SE_DFA = unlist(DFATrendsSE[i, 41:66]) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
        # #for DFA data
        # DFATrends = read.csv("DFATrends.csv")
        # DFATrendsSE = read.csv("DFATrendsSE.csv")
        # # I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] - min(DFATrends[i,2:ncol(DFATrends)]) + 0.1 )
        # CV_DFA = DFATrendsSE[i, 2:ncol(DFATrendsSE)] / DFATrends[i,2:ncol(DFATrends)]
        # SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
        # # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
        # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
        # # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
        index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                    obs = I_DFA, se=SE_DFA)############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
        new_DFA_dat$CPUE = index_DFA.data
        new_DFA_dat$N_cpue = nrow(index_DFA.data)
        
        
        ### Length Comp ####
        if(length == TRUE){
          # From commercial catch
          lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
          lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
          colnames(lencomp_C) = NULL
          LF_C = as.data.frame(t(lencomp_C[,2:ncol(lencomp_C)]))
          Lfreq_Comm = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(1, 2018-1992)), 
                             Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                             Nsamp = (apply(LF_C, 1, sum)), LF_C)
          # From Surveys
          lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
          colnames(lencomp_I1) = NULL
          LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
          Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                           Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                           Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
          lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
          colnames(lencomp_I2) = NULL
          LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
          Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                           Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                           Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
          lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
          colnames(lencomp_I3) = NULL
          LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
          Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                           Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                           Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
          lencomp_I4 = read.csv(paste("LF_Iay4_", i, ".csv", sep=""), header=TRUE, row.names=2)
          colnames(lencomp_I4) = NULL
          LF_I4 = as.data.frame(t(lencomp_I4[,2:ncol(lencomp_I4)]))
          Lfreq_I4 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(5, 2018-1992)), 
                           Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                           Nsamp = (apply(LF_I4, 1, sum)), LF_I4)
          
          lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3, Lfreq_I4)
          
          new_dat$lencomp = lencomp.dat
          new_dat$N_lencomp = nrow(lencomp.dat)
          
          
          # for DFA assessment
          FL = FLoadings[i,1:4]
          lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
          
          LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round((as.numeric(lambda[2])*LF_I2)) + round(as.numeric(lambda[3])*LF_I3) + round(as.numeric(lambda[4])*LF_I4)
          Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                            Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                            Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
          lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
          
          new_DFA_dat$lencomp = lencomp_DFA.dat
          new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
          
        } # end if length == TRUE
        
        
        
        SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
        SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      }# End i loop
      
    } # end if(Nindices=TRUE)
  
  
} # end datfile function
datF2 = function(dir = dir, trial = "Base/Trial1", 
                     Niters = 100, age = FALSE, length = TRUE, Nindices=3, EqCatch=NULL, ...){
  
  
  setwd(paste(dir,trial, sep=""))
  
  # new_dat = base_dat
  # new_DFA_dat = base_DFA_dat
  
  
  
  if(Nindices ==3){
    Cy = read.csv("Cy.csv")
    Iy1 = read.csv("Iy1.csv")
    Iy2 = read.csv("Iy2.csv")
    Iy3 = read.csv("Iy3.csv")
    Iy1_CV = read.csv("Iy1_CV.csv")
    Iy2_CV = read.csv("Iy2_CV.csv")
    Iy3_CV = read.csv("Iy3_CV.csv")
    DFATrends = read.csv("DFATrends.csv")
    DFATrendsSE = read.csv("DFATrendsSE.csv")
    FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
    
    
    for(i in 1:Niters){
      
      new_dat = SS_readdat(file=paste(dir,trial,"\\SS\\SN_SS_", i, "\\SN.dat", sep="") )
      new_DFA_dat = SS_readdat(file=paste(dir,trial,'\\SS_DFA\\SN_SS_',i,"\\SN.dat", sep="") )
      
      
      # catch
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV= unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 3), seas=rep(1, (2018-1992)*3), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992)), 
                              obs = c(I1, I2, I3), se_log = c(se_log1, se_log2, se_log3))
      
      new_dat$CPUE = index.data
      
      
      #for DFA data
      # I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] - min(DFATrends[i,2:ncol(DFATrends)]) + 0.1 )
      CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
      # SE_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) 
      # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # mean_raw = mean(c(I1, I2, I3))
      # min_raw = min(c(I1, I2, I3))
      # trend = unlist(DFATrends[i,2:ncol(DFATrends)])
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2))) ## BASE USE 
      I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + -2* min(DFATrends[i,2:ncol(DFATrends)])  ) # THIS ONE!
      # I_DFA = I2
      # SE_DFA = se_log2
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(0,1)) + 0.1
      # mean2 = mean(I_DFA)
      # I_DFA2 = trend + max(c(I1, I2, I3))
      # I_DFA2 = trend - (min(trend) - min_raw)
      # I_DFA2 = I_DFA + (mean_raw - mean2)
      # I_DFA = I_DFA2
      # mean(I_DFA2)
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2)+(mean_raw - mean2), max(I2)))
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      
      
      
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)
      new_DFA_dat$CPUE = index_DFA.data
      
      
      ### Length Comp ####
      if(length == TRUE){
        # From commercial catch
        lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,13:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=2004:2018, seas = (rep(1, 2018-2003)), FltSry = (rep(1, 2018-2003)), 
                           Sex = (rep(0, 2018-2003)), Part = (rep(0, 2018-2003)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        # From Surveys
        lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3)
        # tail(lencomp.dat, 20)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        # for DFA assessment
        FL = FLoadings[i,1:3]
        # FL = as.vector(c(1,1,1))
        # lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
        lambda = as.vector(ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round((as.numeric(lambda[2])*LF_I2)) + round(as.numeric(lambda[3])*LF_I3)
        # LF_DFA = (LF_I1) + (LF_I2) + (LF_I3)
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        
      } # end if length == TRUE
      
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      # SS_writedat(new_DFA_dat, outfile = "N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/Base/Trial10h/SS_DFA/SN_SS_1/SN.dat", overwrite=TRUE)
      
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
  ### 4 survey indices
  if(Nindices == 4){
    
    Cy = read.csv("Cy.csv")
    Iy1 = read.csv("Iy1.csv")
    Iy2 = read.csv("Iy2.csv")
    Iy3 = read.csv("Iy3.csv")
    Iy4 = read.csv("Iy4.csv")
    Iy1_CV = read.csv("Iy1_CV.csv")
    Iy2_CV = read.csv("Iy2_CV.csv")
    Iy3_CV = read.csv("Iy3_CV.csv")
    Iy4_CV = read.csv("Iy4_CV.csv")
    DFATrends = read.csv("DFATrends.csv")
    DFATrendsSE = read.csv("DFATrendsSE.csv")
    FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
    
    
    for(i in 1:Niters){
      
      new_dat = SS_readdat(file=paste(dir,trial,"\\SS\\SN_SS_", i, "\\SN.dat", sep="") )
      new_DFA_dat = SS_readdat(file=paste(dir,trial,'\\SS_DFA\\SN_SS_',i,"\\SN.dat", sep="") )
      
      
      # catch
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      I4 = unlist(Iy4[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV = unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      I4_CV = unlist(Iy4_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      se_log4 = sqrt(log(1+(I4_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 4), seas=rep(1, (2018-1992)*4), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992), rep(5, 2018-1992)), 
                              obs = c(I1, I2, I3, I4), se_log = c(se_log1, se_log2, se_log3, se_log4))
      
      new_dat$CPUE = index.data
      new_dat$N_cpue = nrow(index.data)
      new_dat$Nsurveys = Nindices
      new_dat$fleetnames = c("Comm_Fish", "Survey1", "Survey2", "Survey3", "Survey4")
      new_dat$surveytiming = c(-1, 0.5, 0.5, 0.5, 0.5)
      new_dat$areas = c(rep(1,5))
      
      
      #for DFA data
      CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) 
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
      I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + (-2 * min(DFATrends[i,2:ncol(DFATrends)])) ) # TRY THIS!!! MUCH BETTER!
      
      # I_DFA = unlist(DFATrends[i,41:66])
      # SE_DFA = unlist(DFATrendsSE[i, 41:66]) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      
      
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
      new_DFA_dat$CPUE = index_DFA.data
      
      
      ### Length Comp ####
      if(length == TRUE){
        # From commercial catch
        lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,13:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=2004:2018, seas = (rep(1, 2018-2003)), FltSry = (rep(1, 2018-2003)), 
                           Sex = (rep(0, 2018-2003)), Part = (rep(0, 2018-2003)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        # From Surveys
        lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        lencomp_I4 = read.csv(paste("LF_Iay4_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I4) = NULL
        LF_I4 = as.data.frame(t(lencomp_I4[,2:ncol(lencomp_I4)]))
        Lfreq_I4 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(5, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I4, 1, sum)), LF_I4)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3, Lfreq_I4)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        # for DFA assessment
        FL = FLoadings[i,1:4]
        lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
        LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round(as.numeric(lambda[2])*LF_I2) + 
          round(as.numeric(lambda[3])*LF_I3) + round(as.numeric(lambda[4])*LF_I4)
        # LF_DFA = LF_I1+LF_I2+LF_I3+LF_I4
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
} # end datfile function

datF3 = function(dir = dir, trial = "Base/Trial1", 
                     Niters = 100, age = FALSE, length = TRUE, Nindices=3, EqCatch=NULL, ...){
  
  setwd(paste(dir,trial, sep=""))
  
  # new_dat = base_dat
  # new_DFA_dat = base_DFA_dat
  
  if(Nindices ==3){
    Cy = read.csv("Cy.csv")
    Iy1 = read.csv("Iy1.csv")
    Iy2 = read.csv("Iy2.csv")
    Iy3 = read.csv("Iy3.csv")
    Iy1_CV = read.csv("Iy1_CV.csv")
    Iy2_CV = read.csv("Iy2_CV.csv")
    Iy3_CV = read.csv("Iy3_CV.csv")
    DFATrends = read.csv("DFATrends.csv")
    DFATrendsSE = read.csv("DFATrendsSE.csv")
    FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
    
    
    for(i in 1:Niters){
      
      new_dat = SS_readdat(file=paste(dir,trial,"\\SS\\SN_SS_", i, "\\SN.dat", sep="") )
      new_DFA_dat = SS_readdat(file=paste(dir,trial,'\\SS_DFA\\SN_SS_',i,"\\SN.dat", sep="") )
      
      
      # catch
      Cy = read.csv("Cy.csv")
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      Iy1 = read.csv("Iy1.csv")
      Iy2 = read.csv("Iy2.csv")
      Iy3 = read.csv("Iy3.csv")
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      # CV
      Iy1_CV = read.csv("Iy1_CV.csv")
      Iy2_CV = read.csv("Iy2_CV.csv")
      Iy3_CV = read.csv("Iy3_CV.csv")
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV= unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 3), seas=rep(1, (2018-1992)*3), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992)), 
                              obs = c(I1, I2, I3), se_log = c(se_log1, se_log2, se_log3))
      
      new_dat$CPUE = index.data
      
      
      #for DFA data
      DFATrends = read.csv("DFATrends.csv")
      DFATrendsSE = read.csv("DFATrendsSE.csv")
      # I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] - min(DFATrends[i,2:ncol(DFATrends)]) + 0.1 )
      CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2))) # BASE
      I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + (-2 * min(DFATrends[i,2:ncol(DFATrends)])) ) # TRY THIS!!! MUCH BETTER!
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      
      
      
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
      new_DFA_dat$CPUE = index_DFA.data
      
      
      ### Length Comp ####
      if(length == TRUE){
        # From commercial catch
        lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,2:ncol(lencomp_C)]))[1:11,]
        Lfreq_Comm = cbind(year=1993:2003, seas = (rep(1, 2003-1992)), FltSry = (rep(1, 2003-1992)), 
                           Sex = (rep(0, 2003-1992)), Part = (rep(0, 2003-1992)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        # From Surveys
        lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        
        # for DFA assessment
        FL = FLoadings[i,1:3]
        lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round((as.numeric(lambda[2])*LF_I2)) + round(as.numeric(lambda[3])*LF_I3)
        # LF_DFA = LF_I1+LF_I2+LF_I3
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
  ### 4 survey indices
  if(Nindices == 4){
    Cy = read.csv("Cy.csv")
    Iy1 = read.csv("Iy1.csv")
    Iy2 = read.csv("Iy2.csv")
    Iy3 = read.csv("Iy3.csv")
    Iy4 = read.csv("Iy4.csv")
    Iy1_CV = read.csv("Iy1_CV.csv")
    Iy2_CV = read.csv("Iy2_CV.csv")
    Iy3_CV = read.csv("Iy3_CV.csv")
    Iy4_CV = read.csv("Iy4_CV.csv")
    DFATrends = read.csv("DFATrends.csv")
    DFATrendsSE = read.csv("DFATrendsSE.csv")
    FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
    
    for(i in 1:Niters){
      
      new_dat = SS_readdat(file=paste(dir,trial,"\\SS\\SN_SS_", i, "\\SN.dat", sep="") )
      new_DFA_dat = SS_readdat(file=paste(dir,trial,'\\SS_DFA\\SN_SS_',i,"\\SN.dat", sep="") )
      
      
      # catch
      Cy = read.csv("Cy.csv")
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      I4 = unlist(Iy4[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV = unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      I4_CV = unlist(Iy4_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      se_log4 = sqrt(log(1+(I4_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 4), seas=rep(1, (2018-1992)*4), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992), rep(5, 2018-1992)), 
                              obs = c(I1, I2, I3, I4), se_log = c(se_log1, se_log2, se_log3, se_log4))
      
      new_dat$CPUE = index.data
      new_dat$N_cpue = nrow(index.data)
      new_dat$Nsurveys = Nindices
      new_dat$fleetnames = c("Comm_Fish", "Survey1", "Survey2", "Survey3", "Survey4")
      new_dat$surveytiming = c(-1, 0.5, 0.5, 0.5, 0.5)
      new_dat$areas = c(rep(1,5))
      
      
      #for DFA data
      CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
      SE_DFA = sqrt(log(1+(CV_DFA^2)))
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
      I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + (-2 * min(DFATrends[i,2:ncol(DFATrends)])) ) # TRY THIS!!! MUCH BETTER!
      # I_DFA = unlist(DFATrends[i,41:66])
      # SE_DFA = unlist(DFATrendsSE[i, 41:66]) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # # #for DFA data
      # DFATrends = read.csv("DFATrends.csv")
      # DFATrendsSE = read.csv("DFATrendsSE.csv")
      # # I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] - min(DFATrends[i,2:ncol(DFATrends)]) + 0.1 )
      # CV_DFA = DFATrendsSE[i, 2:ncol(DFATrendsSE)] / DFATrends[i,2:ncol(DFATrends)]
      # SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
      # # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)
      new_DFA_dat$CPUE = index_DFA.data
      new_DFA_dat$N_cpue = nrow(index_DFA.data)
      
      
      ### Length Comp ####
      if(length == TRUE){
        # From commercial catch
        lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,2:ncol(lencomp_C)]))[1:11,]
        Lfreq_Comm = cbind(year=1993:2003, seas = (rep(1, 2003-1992)), FltSry = (rep(1, 2003-1992)), 
                           Sex = (rep(0, 2003-1992)), Part = (rep(0, 2003-1992)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        
        # From Surveys
        lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        lencomp_I4 = read.csv(paste("LF_Iay4_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I4) = NULL
        LF_I4 = as.data.frame(t(lencomp_I4[,2:ncol(lencomp_I4)]))
        Lfreq_I4 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(5, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I4, 1, sum)), LF_I4)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3, Lfreq_I4)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        # for DFA assessment
        FL = FLoadings[i,1:4]
        lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round((as.numeric(lambda[2])*LF_I2)) + round(as.numeric(lambda[3])*LF_I3) + round(as.numeric(lambda[4])*LF_I4)
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
} # end datfile function

datF4 = function(dir = dir, trial = "Base/Trial1",  
                     Niters = 100, age = FALSE, length = TRUE, Nindices=3, EqCatch=NULL, ...){
  
  setwd(paste(dir,trial, sep=""))
  
  # new_dat = base_dat
  # new_DFA_dat = base_DFA_dat
  
  if(Nindices ==3){
    Cy = read.csv("Cy.csv")
    Iy1 = read.csv("Iy1.csv")
    Iy2 = read.csv("Iy2.csv")
    Iy3 = read.csv("Iy3.csv")
    Iy1_CV = read.csv("Iy1_CV.csv")
    Iy2_CV = read.csv("Iy2_CV.csv")
    Iy3_CV = read.csv("Iy3_CV.csv")
    DFATrends = read.csv("DFATrends.csv")
    DFATrendsSE = read.csv("DFATrendsSE.csv")
    FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
    
    
    for(i in 1:Niters){
      
      new_dat = SS_readdat(file=paste(dir,trial,"\\SS\\SN_SS_", i, "\\SN.dat", sep="") )
      new_DFA_dat = SS_readdat(file=paste(dir,trial,'\\SS_DFA\\SN_SS_',i,"\\SN.dat", sep="") )
      
      
      # catch
      Cy = read.csv("Cy.csv")
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      Iy1 = read.csv("Iy1.csv")
      Iy2 = read.csv("Iy2.csv")
      Iy3 = read.csv("Iy3.csv")
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      # CV
      Iy1_CV = read.csv("Iy1_CV.csv")
      Iy2_CV = read.csv("Iy2_CV.csv")
      Iy3_CV = read.csv("Iy3_CV.csv")
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV= unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 3), seas=rep(1, (2018-1992)*3), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992)), 
                              obs = c(I1, I2, I3), se_log = c(se_log1, se_log2, se_log3))
      
      new_dat$CPUE = index.data
      
      
      #for DFA data
      DFATrends = read.csv("DFATrends.csv")
      DFATrendsSE = read.csv("DFATrendsSE.csv")
      # I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] - min(DFATrends[i,2:ncol(DFATrends)]) + 0.1 )
      CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
      I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + (-2 * min(DFATrends[i,2:ncol(DFATrends)])) ) # TRY THIS!!! MUCH BETTER!
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      
      
      
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
      new_DFA_dat$CPUE = index_DFA.data
      
      
      ### Length Comp ####
      if(length == TRUE){
        # From commercial catch
        lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,3:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=1994:2018, seas = (rep(1, 2018-1993)), FltSry = (rep(1, 2018-1993)), 
                           Sex = (rep(0, 2018-1993)), Part = (rep(0, 2018-1993)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        # From Surveys
        lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        
        # for DFA assessment
        FL = FLoadings[i,1:3]
        lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round((as.numeric(lambda[2])*LF_I2)) + round(as.numeric(lambda[3])*LF_I3)
        # LF_DFA = LF_I1+LF_I2+LF_I3
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
  ### 4 survey indices
  if(Nindices == 4){
    Cy = read.csv("Cy.csv")
    Iy1 = read.csv("Iy1.csv")
    Iy2 = read.csv("Iy2.csv")
    Iy3 = read.csv("Iy3.csv")
    Iy4 = read.csv("Iy4.csv")
    Iy1_CV = read.csv("Iy1_CV.csv")
    Iy2_CV = read.csv("Iy2_CV.csv")
    Iy3_CV = read.csv("Iy3_CV.csv")
    Iy4_CV = read.csv("Iy4_CV.csv")
    DFATrends = read.csv("DFATrends.csv")
    DFATrendsSE = read.csv("DFATrendsSE.csv")
    FLoadings = read.csv(paste("FactorLoadings_Sum.csv", sep=""), header=TRUE, row.names=1)
    
    for(i in 1:Niters){
      
      new_dat = SS_readdat(file=paste(dir,trial,"\\SS\\SN_SS_", i, "\\SN.dat", sep="") )
      new_DFA_dat = SS_readdat(file=paste(dir,trial,'\\SS_DFA\\SN_SS_',i,"\\SN.dat", sep="") )
      
      
      # catch
      catch = Cy[i,] #100 years for iteration i
      catch = unlist(catch)
      catch.data = data.frame(Comm_Fish = catch[41:66], year = 1993:2018, seas = rep(1, length(1993:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      new_DFA_dat$init_equil = ifelse(is.null(EqCatch), catch.data[1,1], EqCatch)
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = unlist(Iy1[i,41:66])
      I2 = unlist(Iy2[i,41:66])
      I3 = unlist(Iy3[i,41:66])
      I4 = unlist(Iy4[i,41:66])
      # CV
      I1_CV = unlist(Iy1_CV[i, 41:66])
      I2_CV = unlist(Iy2_CV[i, 41:66])
      I3_CV = unlist(Iy3_CV[i, 41:66])
      I4_CV = unlist(Iy4_CV[i, 41:66])
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      se_log4 = sqrt(log(1+(I4_CV^2)))
      
      index.data = data.frame(year = rep(1993:2018, 4), seas=rep(1, (2018-1992)*4), 
                              index = c(rep(2, 2018-1992), rep(3, 2018-1992), rep(4, 2018-1992), rep(5, 2018-1992)), 
                              obs = c(I1, I2, I3, I4), se_log = c(se_log1, se_log2, se_log3, se_log4))
      
      new_dat$CPUE = index.data
      new_dat$N_cpue = nrow(index.data)
      new_dat$Nsurveys = Nindices
      new_dat$fleetnames = c("Comm_Fish", "Survey1", "Survey2", "Survey3", "Survey4")
      new_dat$surveytiming = c(-1, 0.5, 0.5, 0.5, 0.5)
      new_dat$areas = c(rep(1,5))
      
      
      #for DFA data
      CV_DFA = unlist(DFATrendsSE[i, 2:ncol(DFATrendsSE)]) / unlist(DFATrends[i,2:ncol(DFATrends)])
      SE_DFA = sqrt(log(1+(CV_DFA^2))) 
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
      I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] + (-2 * min(DFATrends[i,2:ncol(DFATrends)])) ) # TRY THIS!!! MUCH BETTER!
      
      # I_DFA = unlist(DFATrends[i,41:66])
      # SE_DFA = unlist(DFATrendsSE[i, 41:66]) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # #for DFA data
      # DFATrends = read.csv("DFATrends.csv")
      # DFATrendsSE = read.csv("DFATrendsSE.csv")
      # # I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] - min(DFATrends[i,2:ncol(DFATrends)]) + 0.1 )
      # CV_DFA = DFATrendsSE[i, 2:ncol(DFATrendsSE)] / DFATrends[i,2:ncol(DFATrends)]
      # SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      # # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2)))
      # # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      index_DFA.data = data.frame(year = 1993:2018, seas=rep(1, 2018-1992), index=rep(2, 2018-1992), 
                                  obs = I_DFA, se=SE_DFA)
      new_DFA_dat$CPUE = index_DFA.data
      new_DFA_dat$N_cpue = nrow(index_DFA.data)
      
      
      ### Length Comp ####
      if(length == TRUE){
        # From commercial catch
        lencomp_C = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE, row.names=2)
        lencomp_C_labels = read.csv(paste("LF_Cay_", i, ".csv", sep=""), header=TRUE)[,2]
        colnames(lencomp_C) = NULL
        LF_C = as.data.frame(t(lencomp_C[,3:ncol(lencomp_C)]))
        Lfreq_Comm = cbind(year=1994:2018, seas = (rep(1, 2018-1993)), FltSry = (rep(1, 2018-1993)), 
                           Sex = (rep(0, 2018-1993)), Part = (rep(0, 2018-1993)), 
                           Nsamp = (apply(LF_C, 1, sum)), LF_C)
        
        # From Surveys
        lencomp_I1 = read.csv(paste("LF_Iay1_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I1) = NULL
        LF_I1 = as.data.frame(t(lencomp_I1[,2:ncol(lencomp_I1)]))
        Lfreq_I1 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I1, 1, sum)), LF_I1)
        lencomp_I2 = read.csv(paste("LF_Iay2_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I2) = NULL
        LF_I2 = as.data.frame(t(lencomp_I2[,2:ncol(lencomp_I2)]))
        Lfreq_I2 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(3, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I2, 1, sum)), LF_I2)
        lencomp_I3 = read.csv(paste("LF_Iay3_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I3) = NULL
        LF_I3 = as.data.frame(t(lencomp_I3[,2:ncol(lencomp_I3)]))
        Lfreq_I3 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(4, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I3, 1, sum)), LF_I3)
        lencomp_I4 = read.csv(paste("LF_Iay4_", i, ".csv", sep=""), header=TRUE, row.names=2)
        colnames(lencomp_I4) = NULL
        LF_I4 = as.data.frame(t(lencomp_I4[,2:ncol(lencomp_I4)]))
        Lfreq_I4 = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(5, 2018-1992)), 
                         Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                         Nsamp = (apply(LF_I4, 1, sum)), LF_I4)
        
        lencomp.dat = rbind(Lfreq_Comm, Lfreq_I1, Lfreq_I2, Lfreq_I3, Lfreq_I4)
        
        new_dat$lencomp = lencomp.dat
        new_dat$N_lencomp = nrow(lencomp.dat)
        
        
        # for DFA assessment
        FL = FLoadings[i,1:4]
        lambda = data.frame(ifelse(FL<0, 0, FL/max(FL)))
        
        LF_DFA = round(as.numeric(lambda[1])*LF_I1) + round((as.numeric(lambda[2])*LF_I2)) + round(as.numeric(lambda[3])*LF_I3) + round(as.numeric(lambda[4])*LF_I4)
        Lfreq_DFA = cbind(year=1993:2018, seas = (rep(1, 2018-1992)), FltSry = (rep(2, 2018-1992)), 
                          Sex = (rep(0, 2018-1992)), Part = (rep(0, 2018-1992)), 
                          Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
        lencomp_DFA.dat = rbind(Lfreq_Comm, Lfreq_DFA)
        
        new_DFA_dat$lencomp = lencomp_DFA.dat
        new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
        
      } # end if length == TRUE
      
      
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
      SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SN_SS_", i, "/SN.dat", sep=""), overwrite=TRUE)
    }# End i loop
    
  } # end if(Nindices=TRUE)
  
  
  
} # end datfile function



####################################################################################################################################