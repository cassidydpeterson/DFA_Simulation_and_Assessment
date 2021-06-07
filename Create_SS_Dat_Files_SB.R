#### TURN SIMULATION RESULTS INTO DATA FILES
#### SANDBAR SHARK 
#### March 2019



# read in packages
library(r4ss)
# library(scales)

#############################
### READ IN ###

base_dat_1 = SS_readdat(file="D:\\vspace1\\DFA_Simulation\\SB\\Trial101\\SS\\SB_SS_1\\sandbar.dat")
base_dat_DFA_1 = SS_readdat(file="D:\\vspace1\\DFA_Simulation\\SB\\Trial101\\SS_DFA\\SB_SS_1\\sandbar.dat")

base_dat_M_1 = SS_readdat(file="D:\\vspace1\\DFA_Simulation\\SB\\Trial1\\SS\\SB_SS_1\\sandbar.dat")
base_dat_DFAM_1 = SS_readdat(file="D:\\vspace1\\DFA_Simulation\\SB\\Trial1\\SS_DFA\\SB_SS_1\\sandbar.dat")

dir1="D:/vspace1/DFA_Simulation/SB/"


# base_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\SB\\Build_SS_mod\\sandbar.dat")
# # base_DFA_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\SB\\Build_SS_mod\\DFA_SS\\SB.dat")
# base_dat_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\SB\\Build_SS_mod\\SS_NewDat_1\\sandbar.dat")
# base_dat_DFA_1 = SS_readdat(file="N:\\Documents\\DFA_Simulation\\SB\\Build_SS_mod\\SS_DFA_NewDat_1\\sandbar.dat")

# TEST
trial = "Trial1"
trial = "Trial101"
base_dat = base_dat_1
base_DFA_dat = base_dat_DFA_1
Niters = 100 
# dir="N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/"
dir1="D:/vspace1/DFA_Simulation/SB/"
EqCatch=NULL
missing=TRUE

####################################################################################################################################
# ### WRITE FUNCTION ####
####################################################################################################################################
#
datfile_SS = function(dir = dir, trial = "Trial1", base_dat = base_dat, base_DFA_dat = base_DFA_dat, 
                   Niters = 100, missing=TRUE, EqCatch=c(NULL,NULL,NULL,NULL), ...){
  
  setwd(paste(dir,trial, sep=""))
  
  new_dat = base_dat
  new_DFA_dat = base_DFA_dat
  
  load("ResultsList.RData")
  load("DFA_Results.RData")
  # load("N:/Documents/DFA_Simulation/SB/Trial1/ResultsList.RData")
  # load("N:/Documents/DFA_Simulation/SB/Trial1/DFA_Results.RData")
  # > names(ResultsList)
  # "Iy1"    "Iy2"    "Iy3"    "Iy4"    "Iy5"    "Iy6"    "Iy7"    "Iy1_CV" "Iy2_CV" "Iy3_CV" "Iy4_CV" "Iy5_CV" "Iy6_CV" "Iy7_CV" "M0"     
  # "FSS"    "Npups"  "Cy"  "Cy_1"  "Cy_2"  "Cy_3"  "Cy_4"   "Ny" 
  # > names(DFA_Results)
  # [1] "FitRatios"       "DFATrends"       "DFATrendsSE"     "upCI_DFATrends"  "lowCI_DFATrends" "FLoadings"      
  # "bias"            "biasAll"         "RMSE"
  
  if(missing==FALSE){
    
    for(i in 1:Niters){
      # catch
      catch = unlist(ResultsList$Cy[i,]) #100 years for iteration i
      catch.data = data.frame(Comm_Fish1 = ResultsList$Cy_1[i,51:90], Comm_Fish2 = ResultsList$Cy_2[i,51:90], 
                              Comm_Fish3 = ResultsList$Cy_3[i,51:90], Comm_Fish4 = ResultsList$Cy_4[i,51:90],
                              year = 1979:2018, seas = rep(1, length(1979:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = unlist(ifelse(c(is.null(EqCatch),is.null(EqCatch),is.null(EqCatch),is.null(EqCatch)), catch.data[1,1:4], EqCatch))
      new_DFA_dat$init_equil = unlist(ifelse(c(is.null(EqCatch),is.null(EqCatch),is.null(EqCatch),is.null(EqCatch)), catch.data[1,1:4], EqCatch))
      
      
      # CPUE [year, seas, index, obs, se_log]
      # Index
      I1 = ResultsList$Iy1[i,51:90]
      I2 = ResultsList$Iy2[i,51:90]
      I3 = ResultsList$Iy3[i,51:90]
      I4 = ResultsList$Iy4[i,51:90]
      I5 = ResultsList$Iy5[i,51:90]
      I6 = ResultsList$Iy6[i,51:90]
      I7 = ResultsList$Iy7[i,51:90]
      
      # CV
      I1_CV = ResultsList$Iy1_CV[i, 51:90]
      I2_CV = ResultsList$Iy2_CV[i, 51:90]
      I3_CV = ResultsList$Iy3_CV[i, 51:90]
      I4_CV = ResultsList$Iy4_CV[i, 51:90]
      I5_CV = ResultsList$Iy5_CV[i, 51:90]
      I6_CV = ResultsList$Iy6_CV[i, 51:90]
      I7_CV = ResultsList$Iy7_CV[i, 51:90]
      
      se_log1 = sqrt(log(1+(I1_CV^2)))
      se_log2 = sqrt(log(1+(I2_CV^2)))
      se_log3 = sqrt(log(1+(I3_CV^2)))
      se_log4 = sqrt(log(1+(I4_CV^2)))
      se_log5 = sqrt(log(1+(I5_CV^2)))
      se_log6 = sqrt(log(1+(I6_CV^2)))
      se_log7 = sqrt(log(1+(I7_CV^2)))
      
      
      index.data = data.frame(year = rep(1979:2018, 7), seas=rep(1, length(1979:2018)*7), 
                              index = rep(5:11, each=length(1979:2018)), obs = c(I1, I2, I3, I4, I5, I6, I7), 
                              se_log = c(se_log1, se_log2, se_log3, se_log4, se_log5, se_log6, se_log7))
      
      new_dat$CPUE = index.data
      
      
      
      #for DFA data
      I_DFA = DFA_Results$DFATrendsBT[i,] 
      SE_DFA = DFA_Results$DFATrendsSEBT[i,]
      # I_DFA = unlist(DFATrends[i,2:ncol(DFATrends)] - min(DFATrends[i,2:ncol(DFATrends)]) + 0.1 )
      # I_DFA = DFA_Results$DFATrends[i,] + (-2 * min(DFA_Results$DFATrends[i,]))  # TRY THIS!!! MUCH BETTER!
      # up_CI = DFA_Results$upCI_DFATrends[i,] + (-2 * min(DFA_Results$DFATrends[i,]))
      # SD = (up_CI - I_DFA )/1.96
      # CV_DFA = SD/ I_DFA
      # CV_DFA = DFA_Results$DFATrendsSE[i,] / DFA_Results$DFATrends[i,]
      # SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
      
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I2), max(I2))) # BASE
      # I_DFA = unlist(DFATrends[i, 2:ncol(DFATrends)] + 10)
      # mean_raw = mean(c(I1, I2, I3))
      # min_raw = min(c(I1, I2, I3))
      # trend = unlist(DFATrends[i,2:ncol(DFATrends)])
      # mean2 = mean(I_DFA)
      # I_DFA2 = I_DFA + (mean_raw - mean2)
      # I_DFA = I_DFA2
      #mean(I_DFA2)
      # I_DFA = rescale(unlist(DFATrends[i,2:ncol(DFATrends)]), to = c(min(I1, I2, I3), max(I1, I2, I3)))
      
      
      index_DFA.data = data.frame(year = 1979:2018, seas=rep(1, length(1979:2018)), index=rep(5, length(1979:2018)), 
                                  obs = I_DFA, se=SE_DFA)
      ############## CHECK FORMAT OF SEs FOR THIS !!!!!!!!!!!!!!!
      new_DFA_dat$CPUE = index_DFA.data
      
      
      ### Length Comp ####
      load(paste("LFreqList_",i,".RData", sep=""))
      # load("N:/Documents/DFA_Simulation/SB/Trial1/LFreqList_1.RData")
      # load(paste("N:/Documents/DFA_Simulation/SB/",trial,"/LFreqList_", i, ".RData", sep=""))
      # names(LFreqList)
      # 1.  "LF_Cay1F_2"  2. "LF_Cay2F_2"  3. "LF_Cay3F_2"  4. "LF_Cay4F_2" 
      # 5.  "LF_Cay1M_2"  6. "LF_Cay2M_2"  7. "LF_Cay3M_2"  8. "LF_Cay4M_2" 
      # 9.  "LF_Iay1F_2" 10. "LF_Iay2F_2" 11. "LF_Iay3F_2" 12. "LF_Iay4F_2" 13. "LF_Iay5F_2" 14. "LF_Iay6F_2" 15. "LF_Iay7F_2" 
      # 16. "LF_Iay1M_2" 17. "LF_Iay2M_2" 18. "LF_Iay3M_2" 19. "LF_Iay4M_2" 20. "LF_Iay5M_2" 21. "LF_Iay6M_2" 22. "LF_Iay7M_2"
      
      
      
      
      
      
      ########## ADD A CAVEAT IN CODE THAT REMOVES LINE IF Nsamp=0 ####################################
      
      
      
      
      # From commercial catch ------------------
      # C1
      LF_C1_F = as.data.frame(t(LFreqList[[1]]))
      colnames(LF_C1_F) = LF_C1_F[1,]
      LF_C1_F = LF_C1_F[-1,]
      LF_C1_M = as.data.frame(t(LFreqList[[5]]))
      colnames(LF_C1_M) = LF_C1_M[1,]
      LF_C1_M = LF_C1_M[-1,]
      LF_C1 = cbind(LF_C1_F, LF_C1_M)
      Lfreq_C1 = cbind(year=1979:2018, seas = (rep(1,length(1979:2018))), FltSry = (rep(1, length(1979:2018))), 
                       Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                       Nsamp = (apply(LF_C1, 1, sum)), LF_C1)
      
      # C2
      LF_C2_F = as.data.frame(t(LFreqList[[2]]))
      colnames(LF_C2_F) = LF_C2_F[1,]
      LF_C2_F = LF_C2_F[-1,]
      LF_C2_M = as.data.frame(t(LFreqList[[6]]))
      colnames(LF_C2_M) = LF_C2_M[1,]
      LF_C2_M = LF_C2_M[-1,]
      LF_C2 = cbind(LF_C2_F, LF_C2_M)
      Lfreq_C2 = cbind(year=1979:2018, seas = (rep(1,length(1979:2018))), FltSry = (rep(2, length(1979:2018))), 
                       Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                       Nsamp = (apply(LF_C2, 1, sum)), LF_C2)
      
      # C3
      LF_C3_F = as.data.frame(t(LFreqList[[3]]))
      colnames(LF_C3_F) = LF_C3_F[1,]
      LF_C3_F = LF_C3_F[-1,]
      LF_C3_M = as.data.frame(t(LFreqList[[7]]))
      colnames(LF_C3_M) = LF_C3_M[1,]
      LF_C3_M = LF_C3_M[-1,]
      LF_C3 = cbind(LF_C3_F, LF_C3_M)
      Lfreq_C3 = cbind(year=1979:2018, seas = (rep(1,length(1979:2018))), FltSry = (rep(3, length(1979:2018))), 
                       Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                       Nsamp = (apply(LF_C3, 1, sum)), LF_C3)
      
      # C4
      LF_C4_F = as.data.frame(t(LFreqList[[4]]))
      colnames(LF_C4_F) = LF_C4_F[1,]
      LF_C4_F = LF_C4_F[-1,]
      LF_C4_M = as.data.frame(t(LFreqList[[8]]))
      colnames(LF_C4_M) = LF_C4_M[1,]
      LF_C4_M = LF_C4_M[-1,]
      LF_C4 = cbind(LF_C4_F, LF_C4_M)
      Lfreq_C4_2 = cbind(year=1979:2018, seas = (rep(1,length(1979:2018))), FltSry = (rep(4, length(1979:2018))), 
                         Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                         Nsamp = (apply(LF_C4, 1, sum)), LF_C4)
      Lfreq_C4 = Lfreq_C4_2[Lfreq_C4_2[,6] > 0,]
      
      
      # From Surveys -----------------------
      # S1
      LF_S1_F = as.data.frame(t(LFreqList[[9]]))
      colnames(LF_S1_F) = LF_S1_F[1,]
      LF_S1_F = LF_S1_F[-1,]
      LF_S1_M = as.data.frame(t(LFreqList[[16]]))
      colnames(LF_S1_M) = LF_S1_M[1,]
      LF_S1_M = LF_S1_M[-1,]
      LF_S1 = cbind(LF_S1_F, LF_S1_M)
      Lfreq_S1_2 = cbind(year=1979:2018, seas = (rep(1,length(1979:2018))), FltSry = (rep(5, length(1979:2018))), 
                       Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                       Nsamp = (apply(LF_S1, 1, sum)), LF_S1)
      Lfreq_S1 = Lfreq_S1_2[Lfreq_S1_2[,6] > 0,]
      
      
      # S2
      LF_S2_F = as.data.frame(t(LFreqList[[10]]))
      colnames(LF_S2_F) = LF_S2_F[1,]
      LF_S2_F = LF_S2_F[-1,]
      LF_S2_M = as.data.frame(t(LFreqList[[17]]))
      colnames(LF_S2_M) = LF_S2_M[1,]
      LF_S2_M = LF_S2_M[-1,]
      LF_S2 = cbind(LF_S2_F, LF_S2_M)
      Lfreq_S2_2 = cbind(year=1979:2018, seas = (rep(1,length(1979:2018))), FltSry = (rep(6, length(1979:2018))), 
                       Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                       Nsamp = (apply(LF_S2, 1, sum)), LF_S2)
      Lfreq_S2 = Lfreq_S2_2[Lfreq_S2_2[,6] > 0,]
      
      
      # S3
      LF_S3_F = as.data.frame(t(LFreqList[[11]]))
      colnames(LF_S3_F) = LF_S3_F[1,]
      LF_S3_F = LF_S3_F[-1,]
      LF_S3_M = as.data.frame(t(LFreqList[[18]]))
      colnames(LF_S3_M) = LF_S3_M[1,]
      LF_S3_M = LF_S3_M[-1,]
      LF_S3 = cbind(LF_S3_F, LF_S3_M)
      Lfreq_S3_2 = cbind(year=1979:2018, seas = (rep(1,length(1979:2018))), FltSry = (rep(7, length(1979:2018))), 
                       Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                       Nsamp = (apply(LF_S3, 1, sum)), LF_S3)
      Lfreq_S3 = Lfreq_S3_2[Lfreq_S3_2[,6] > 0,]
      
      
      # S4
      LF_S4_F = as.data.frame(t(LFreqList[[12]]))
      colnames(LF_S4_F) = LF_S4_F[1,]
      LF_S4_F = LF_S4_F[-1,]
      LF_S4_M = as.data.frame(t(LFreqList[[19]]))
      colnames(LF_S4_M) = LF_S4_M[1,]
      LF_S4_M = LF_S4_M[-1,]
      LF_S4 = cbind(LF_S4_F, LF_S4_M)
      Lfreq_S4_2 = cbind(year=1979:2018, seas = (rep(1,length(1979:2018))), FltSry = (rep(8, length(1979:2018))), 
                       Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                       Nsamp = (apply(LF_S4, 1, sum)), LF_S4)
      Lfreq_S4 = Lfreq_S4_2[Lfreq_S4_2[,6] > 0,]
      
      
      # S5
      LF_S5_F = as.data.frame(t(LFreqList[[13]]))
      colnames(LF_S5_F) = LF_S5_F[1,]
      LF_S5_F = LF_S5_F[-1,]
      LF_S5_M = as.data.frame(t(LFreqList[[20]]))
      colnames(LF_S5_M) = LF_S5_M[1,]
      LF_S5_M = LF_S5_M[-1,]
      LF_S5 = cbind(LF_S5_F, LF_S5_M)
      Lfreq_S5_2 = cbind(year=1979:2018, seas = (rep(1,length(1979:2018))), FltSry = (rep(9, length(1979:2018))), 
                       Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                       Nsamp = (apply(LF_S5, 1, sum)), LF_S5)
      Lfreq_S5 = Lfreq_S5_2[Lfreq_S5_2[,6] > 0,]
      
      
      # S6
      LF_S6_F = as.data.frame(t(LFreqList[[14]]))
      colnames(LF_S6_F) = LF_S6_F[1,]
      LF_S6_F = LF_S6_F[-1,]
      LF_S6_M = as.data.frame(t(LFreqList[[21]]))
      colnames(LF_S6_M) = LF_S6_M[1,]
      LF_S6_M = LF_S6_M[-1,]
      LF_S6 = cbind(LF_S6_F, LF_S6_M)
      Lfreq_S6_2 = cbind(year=1979:2018, seas = (rep(1,length(1979:2018))), FltSry = (rep(10, length(1979:2018))), 
                       Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                       Nsamp = (apply(LF_S6, 1, sum)), LF_S6)
      Lfreq_S6 = Lfreq_S6_2[Lfreq_S6_2[,6] > 0,]
      
      
      # S7
      LF_S7_F = as.data.frame(t(LFreqList[[15]]))
      colnames(LF_S7_F) = LF_S7_F[1,]
      LF_S7_F = LF_S7_F[-1,]
      LF_S7_M = as.data.frame(t(LFreqList[[22]]))
      colnames(LF_S7_M) = LF_S7_M[1,]
      LF_S7_M = LF_S7_M[-1,]
      LF_S7 = cbind(LF_S7_F, LF_S7_M)
      Lfreq_S7_2 = cbind(year=1979:2018, seas = (rep(1,length(1979:2018))), FltSry = (rep(11, length(1979:2018))), 
                       Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                       Nsamp = (apply(LF_S7, 1, sum)), LF_S7)
      Lfreq_S7 = Lfreq_S7_2[Lfreq_S7_2[,6] > 0,]
      
      
      
      #####
      lencomp.dat = rbind(Lfreq_C1, Lfreq_C2, Lfreq_C3, Lfreq_C4, Lfreq_S1, Lfreq_S2, Lfreq_S3, Lfreq_S4, Lfreq_S5, Lfreq_S6, Lfreq_S7)
      
      new_dat$lencomp = lencomp.dat
      new_dat$N_lencomp = nrow(lencomp.dat)
      
      
      
      # for DFA assessment
      FL = DFA_Results$FLoadings[i,1:7]
      lambda = ifelse(FL<0, 0, FL/max(FL))
      
      LF_DFA = round(as.numeric(lambda[1])*LF_S1) + round((as.numeric(lambda[2])*LF_S2)) + round(as.numeric(lambda[3])*LF_S3) + 
        round((as.numeric(lambda[4])*LF_S4)) + round(as.numeric(lambda[5])*LF_S5) + round((as.numeric(lambda[6])*LF_S6)) + 
        round(as.numeric(lambda[7])*LF_S7)
      # LF_DFA = LF_I1+LF_I2+LF_I3
      Lfreq_DFA1 = cbind(year=1979:2018, seas = (rep(1, length(1979:2018))), FltSry = (rep(5, length(1979:2018))), 
                        Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                        Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
      Lfreq_DFA = Lfreq_DFA1[Lfreq_DFA1[,6] > 0,]
      
      
      lencomp_DFA.dat = rbind(Lfreq_C1, Lfreq_C2, Lfreq_C3, Lfreq_C4, Lfreq_DFA)
      
      new_DFA_dat$lencomp = lencomp_DFA.dat
      new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
      
      
      
      SS_writedat(new_dat, outfile = "D:\\vspace1\\DFA_Simulation\\SB\\METHODS TEST\\SB_SS_1\\sandbar.dat", overwrite=TRUE, version="3.24")
      SS_writedat(new_DFA_dat, outfile = "D:\\vspace1\\DFA_Simulation\\SB\\METHODS TEST\\SB_SS_1_DFA\\sandbar.dat", overwrite=TRUE, version="3.24")
      
      SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SB_SS_", i, "/sandbar.dat", sep=""), overwrite=TRUE, version="3.24")
      SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SB_SS_", i, "/sandbar.dat", sep=""), overwrite=TRUE, version="3.24")
      
      # SS_writedat(new_DFA_dat, outfile = "N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/Base/Trial10a/SS_DFA/SN_SS_1/SN.dat", overwrite=TRUE)
    }# End i loop
  }# END IF MISSING = FALSE
    
  
  if(missing==TRUE){
    for(i in 1:Niters){
      # catch
      catch = unlist(ResultsList$Cy[i,]) #100 years for iteration i
      catch.data = data.frame(Comm_Fish1 = ResultsList$Cy_1[i,51:90], Comm_Fish2 = ResultsList$Cy_2[i,51:90], 
                              Comm_Fish3 = ResultsList$Cy_3[i,51:90], Comm_Fish4 = ResultsList$Cy_4[i,51:90],
                              year = 1979:2018, seas = rep(1, length(1979:2018)))
      new_dat$catch = catch.data
      new_DFA_dat$catch = catch.data
      new_dat$N_catch = nrow(catch.data)
      new_DFA_dat$N_catch = nrow(catch.data)
      new_dat$init_equil = unlist(ifelse(c(is.null(EqCatch),is.null(EqCatch),is.null(EqCatch),is.null(EqCatch)), catch.data[1,1:4], EqCatch))
      new_DFA_dat$init_equil = unlist(ifelse(c(is.null(EqCatch),is.null(EqCatch),is.null(EqCatch),is.null(EqCatch)), catch.data[1,1:4], EqCatch))
      
      
    # FOR MISSING DATA ! 
    I1 = ResultsList$Iy1[i,c(51:53,59:64,67:90)]
    Y1 = c(1979:1981,1987:1992,1995:2018)
    I2 = ResultsList$Iy2[i,58:90]
    Y2 = 1986:2018
    I3 = ResultsList$Iy3[i,64:90]
    Y3 = 1992:2018
    I4 = ResultsList$Iy4[i,67:90]
    Y4 = 1995:2018
    I5 = ResultsList$Iy5[i,73:90]
    Y5 = 2001:2018
    I6 = ResultsList$Iy6[i,66:79]
    Y6 = 1994:2007
    I7 = ResultsList$Iy7[i,c(68,71,74,77,80,83,86,89)]
    Y7 = c(1996, 1999, 2002, 2005, 2008, 2011, 2014, 2017)
    
    I1_CV = ResultsList$Iy1_CV[i,c(51:53,59:64,67:90)]
    I2_CV = ResultsList$Iy2_CV[i,58:90]
    I3_CV = ResultsList$Iy3_CV[i,64:90]
    I4_CV = ResultsList$Iy4_CV[i,67:90]
    I5_CV = ResultsList$Iy5_CV[i,73:90]
    I6_CV = ResultsList$Iy6_CV[i,66:79]
    I7_CV = ResultsList$Iy7_CV[i,c(68,71,74,77,80,83,86,89)]
    
    se_log1 = sqrt(log(1+(I1_CV^2)))
    se_log2 = sqrt(log(1+(I2_CV^2)))
    se_log3 = sqrt(log(1+(I3_CV^2)))
    se_log4 = sqrt(log(1+(I4_CV^2)))
    se_log5 = sqrt(log(1+(I5_CV^2)))
    se_log6 = sqrt(log(1+(I6_CV^2)))
    se_log7 = sqrt(log(1+(I7_CV^2)))
    
    Ys = c(Y1, Y2, Y3, Y4, Y5, Y6, Y7)
    index.data = data.frame(year = Ys, seas=rep(1, length(Ys)), 
                            index = c(rep(5, length(Y1)), rep(6, length(Y2)), rep(7, length(Y3)), rep(8, length(Y4)), 
                                      rep(9, length(Y5)), rep(10, length(Y6)), rep(11, length(Y7)) ), 
                            obs = c(I1, I2, I3, I4, I5, I6, I7), 
                            se_log = c(se_log1, se_log2, se_log3, se_log4, se_log5, se_log6, se_log7))
    
    new_dat$CPUE = index.data
    new_dat$N_cpue = nrow(index.data)
    
    
    #for DFA data
    I_DFA = DFA_Results$DFATrendsBT[i,] 
    SE_DFA = DFA_Results$DFATrendsSEBT[i,]
    # I_DFA = DFA_Results$DFATrends[i,] + (-2 * min(DFA_Results$DFATrends[i,]))  # BASE
    # # I_DFA2 = DFA_Results$DFATrends[i,] + (-3 * min(DFA_Results$DFATrends[i,]))  # TRY 
    # # I_DFA3 = DFA_Results$DFATrends[i,] + (-4 * min(DFA_Results$DFATrends[i,]))  # TRY
    # up_CI = DFA_Results$upCI_DFATrends[i,] + (-2 * min(DFA_Results$DFATrends[i,]))
    # SD = (up_CI - I_DFA )/1.96
    # CV_DFA = SD/ I_DFA
    # # CV_DFA = DFA_Results$DFATrendsSE[i,] / DFA_Results$DFATrends[i,]
    # SE_DFA = sqrt(log(1+(CV_DFA^2))) ############ CHECK FOR FORMAT OF SEs FOR THIS@!!!!!!!!!!!!!!!!!
 
    index_DFA.data = data.frame(year = 1979:2018, seas=rep(1, length(1979:2018)), index=rep(5, length(1979:2018)), 
                                obs = I_DFA, se=SE_DFA)
    new_DFA_dat$CPUE = index_DFA.data
    new_DFA_dat$N_cpue = nrow(index_DFA.data)
    
    ### Length Comp ####
    load(paste("LFreqList_",i,".RData", sep=""))
    # load("N:/Documents/DFA_Simulation/SB/Trial1/LFreqList_1.RData")
    # load(paste("N:/Documents/DFA_Simulation/SB/Trial1/LFreqList_", i, ".RData", sep=""))
    # names(LFreqList)
    # 1.  "LF_Cay1F_2"  2. "LF_Cay2F_2"  3. "LF_Cay3F_2"  4. "LF_Cay4F_2" 
    # 5.  "LF_Cay1M_2"  6. "LF_Cay2M_2"  7. "LF_Cay3M_2"  8. "LF_Cay4M_2" 
    # 9.  "LF_Iay1F_2" 10. "LF_Iay2F_2" 11. "LF_Iay3F_2" 12. "LF_Iay4F_2" 13. "LF_Iay5F_2" 14. "LF_Iay6F_2" 15. "LF_Iay7F_2" 
    # 16. "LF_Iay1M_2" 17. "LF_Iay2M_2" 18. "LF_Iay3M_2" 19. "LF_Iay4M_2" 20. "LF_Iay5M_2" 21. "LF_Iay6M_2" 22. "LF_Iay7M_2"
    
    # LF_C1_F = t(LFreqList[[1]])
    # colnames(LF_C1_F) = LF_C1_F[1,]
    # LF_C1_F = LF_C1_F[-1,]
    # 
    # LF_C1_M = t(LFreqList[[5]])
    # colnames(LF_C1_M) = LF_C1_M[1,]
    # LF_C1_M = LF_C1_M[-1,]
    # LF_C1 = cbind(LF_C1_F, LF_C1_M)
    # 
    # nrow(LF_C1_M)
    
    # From commercial catch ------------------
    # C1
    LF_C1_F = as.data.frame(t(LFreqList[[1]]))
    colnames(LF_C1_F) = LF_C1_F[1,]
    LF_C1_F = LF_C1_F[-1,]
    LF_C1_M = as.data.frame(t(LFreqList[[5]]))
    colnames(LF_C1_M) = LF_C1_M[1,]
    LF_C1_M = LF_C1_M[-1,]
    LF_C1 = cbind(LF_C1_F, LF_C1_M)
    Lfreq_C1 = cbind(year=1979:2018, seas = (rep(1,length(1979:2018))), FltSry = (rep(1, length(1979:2018))), 
                     Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                     Nsamp = (apply(LF_C1, 1, sum)), LF_C1)
    
    # C2
    LF_C2_F = as.data.frame(t(LFreqList[[2]]))
    colnames(LF_C2_F) = LF_C2_F[1,]
    LF_C2_F = LF_C2_F[-1,]
    LF_C2_M = as.data.frame(t(LFreqList[[6]]))
    colnames(LF_C2_M) = LF_C2_M[1,]
    LF_C2_M = LF_C2_M[-1,]
    LF_C2 = cbind(LF_C2_F, LF_C2_M)
    Lfreq_C2 = cbind(year=1979:2018, seas = (rep(1,length(1979:2018))), FltSry = (rep(2, length(1979:2018))), 
                     Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                     Nsamp = (apply(LF_C2, 1, sum)), LF_C2)
    
    # C3
    LF_C3_F = as.data.frame(t(LFreqList[[3]]))
    colnames(LF_C3_F) = LF_C3_F[1,]
    LF_C3_F = LF_C3_F[-1,]
    LF_C3_M = as.data.frame(t(LFreqList[[7]]))
    colnames(LF_C3_M) = LF_C3_M[1,]
    LF_C3_M = LF_C3_M[-1,]
    LF_C3 = cbind(LF_C3_F, LF_C3_M)
    Lfreq_C3 = cbind(year=1979:2018, seas = (rep(1,length(1979:2018))), FltSry = (rep(3, length(1979:2018))), 
                     Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                     Nsamp = (apply(LF_C3, 1, sum)), LF_C3)
    
    # C4
    LF_C4_F = as.data.frame(t(LFreqList[[4]]))
    colnames(LF_C4_F) = LF_C4_F[1,]
    LF_C4_F = LF_C4_F[-1,]
    LF_C4_M = as.data.frame(t(LFreqList[[8]]))
    colnames(LF_C4_M) = LF_C4_M[1,]
    LF_C4_M = LF_C4_M[-1,]
    LF_C4 = cbind(LF_C4_F, LF_C4_M)
    Lfreq_C4_2 = cbind(year=1979:2018, seas = (rep(1,length(1979:2018))), FltSry = (rep(4, length(1979:2018))), 
                     Sex = (rep(3, length(1979:2018))), Part = (rep(0, length(1979:2018))), 
                     Nsamp = (apply(LF_C4, 1, sum)), LF_C4)
    Lfreq_C4 = Lfreq_C4_2[Lfreq_C4_2[,6] > 0,]
    
    
    
    # From Surveys -----------------------
    # S1
    LF_S1_F = as.data.frame(t(LFreqList[[9]]))
    colnames(LF_S1_F) = LF_S1_F[1,]
    LF_S1_F = LF_S1_F[-1,]
    LF_S1_M = as.data.frame(t(LFreqList[[16]]))
    colnames(LF_S1_M) = LF_S1_M[1,]
    LF_S1_M = LF_S1_M[-1,]
    LF_S1 = cbind(LF_S1_F, LF_S1_M)
    Lfreq_S1_2 = cbind(year=Y1, seas = (rep(1,length(Y1))), FltSry = (rep(5, length(Y1))), 
                     Sex = (rep(3, length(Y1))), Part = (rep(0, length(Y1))), 
                     Nsamp = (apply(LF_S1, 1, sum)), LF_S1)
    Lfreq_S1 = Lfreq_S1_2[Lfreq_S1_2[,6] > 0,]
    
    # S2
    LF_S2_F = as.data.frame(t(LFreqList[[10]]))
    colnames(LF_S2_F) = LF_S2_F[1,]
    LF_S2_F = LF_S2_F[-1,]
    LF_S2_M = as.data.frame(t(LFreqList[[17]]))
    colnames(LF_S2_M) = LF_S2_M[1,]
    LF_S2_M = LF_S2_M[-1,]
    LF_S2 = cbind(LF_S2_F, LF_S2_M)
    Lfreq_S2_2 = cbind(year=Y2, seas = (rep(1,length(Y2))), FltSry = (rep(6, length(Y2))), 
                     Sex = (rep(3, length(Y2))), Part = (rep(0, length(Y2))), 
                     Nsamp = (apply(LF_S2, 1, sum)), LF_S2)
    Lfreq_S2 = Lfreq_S2_2[Lfreq_S2_2[,6] > 0,]
    
    # S3
    LF_S3_F = as.data.frame(t(LFreqList[[11]]))
    colnames(LF_S3_F) = LF_S3_F[1,]
    LF_S3_F = LF_S3_F[-1,]
    LF_S3_M = as.data.frame(t(LFreqList[[18]]))
    colnames(LF_S3_M) = LF_S3_M[1,]
    LF_S3_M = LF_S3_M[-1,]
    LF_S3 = cbind(LF_S3_F, LF_S3_M)
    Lfreq_S3_2 = cbind(year=Y3, seas = (rep(1,length(Y3))), FltSry = (rep(7, length(Y3))), 
                     Sex = (rep(3, length(Y3))), Part = (rep(0, length(Y3))), 
                     Nsamp = (apply(LF_S3, 1, sum)), LF_S3)
    Lfreq_S3 = Lfreq_S3_2[Lfreq_S3_2[,6] > 0,]
    
    # S4
    LF_S4_F = as.data.frame(t(LFreqList[[12]]))
    colnames(LF_S4_F) = LF_S4_F[1,]
    LF_S4_F = LF_S4_F[-1,]
    LF_S4_M = as.data.frame(t(LFreqList[[19]]))
    colnames(LF_S4_M) = LF_S4_M[1,]
    LF_S4_M = LF_S4_M[-1,]
    LF_S4 = cbind(LF_S4_F, LF_S4_M)
    Lfreq_S4_2 = cbind(year=Y4, seas = rep(1,length(Y4)), FltSry = rep(8, length(Y4)), 
                     Sex = rep(3, length(Y4)), Part = rep(0, length(Y4)), 
                     Nsamp = apply(LF_S4, 1, sum), LF_S4)
    Lfreq_S4 = Lfreq_S4_2[Lfreq_S4_2[,6] > 0,]
    
    # S5
    LF_S5_F = as.data.frame(t(LFreqList[[13]]))
    colnames(LF_S5_F) = LF_S5_F[1,]
    LF_S5_F = LF_S5_F[-1,]
    LF_S5_M = as.data.frame(t(LFreqList[[20]]))
    colnames(LF_S5_M) = LF_S5_M[1,]
    LF_S5_M = LF_S5_M[-1,]
    LF_S5 = cbind(LF_S5_F, LF_S5_M)
    Lfreq_S5_2 = cbind(year=Y5, seas = (rep(1,length(Y5))), FltSry = (rep(9, length(Y5))), 
                     Sex = (rep(3, length(Y5))), Part = (rep(0, length(Y5))), 
                     Nsamp = (apply(LF_S5, 1, sum)), LF_S5)
    Lfreq_S5 = Lfreq_S5_2[Lfreq_S5_2[,6] > 0,]
    
    # S6
    LF_S6_F = as.data.frame(t(LFreqList[[14]]))
    colnames(LF_S6_F) = LF_S6_F[1,]
    LF_S6_F = LF_S6_F[-1,]
    LF_S6_M = as.data.frame(t(LFreqList[[21]]))
    colnames(LF_S6_M) = LF_S6_M[1,]
    LF_S6_M = LF_S6_M[-1,]
    LF_S6 = cbind(LF_S6_F, LF_S6_M)
    Lfreq_S6_2 = cbind(year=Y6, seas = (rep(1,length(Y6))), FltSry = (rep(10, length(Y6))), 
                     Sex = (rep(3, length(Y6))), Part = (rep(0, length(Y6))), 
                     Nsamp = (apply(LF_S6, 1, sum)), LF_S6)
    Lfreq_S6 = Lfreq_S6_2[Lfreq_S6_2[,6] > 0,]
    
    # S7
    LF_S7_F = as.data.frame(t(LFreqList[[15]]))
    colnames(LF_S7_F) = LF_S7_F[1,]
    LF_S7_F = LF_S7_F[-1,]
    LF_S7_M = as.data.frame(t(LFreqList[[22]]))
    colnames(LF_S7_M) = LF_S7_M[1,]
    LF_S7_M = LF_S7_M[-1,]
    LF_S7 = cbind(LF_S7_F, LF_S7_M)
    Lfreq_S7_2 = cbind(year=Y7, seas = (rep(1,length(Y7))), FltSry = (rep(11, length(Y7))), 
                     Sex = (rep(3, length(Y7))), Part = (rep(0, length(Y7))), 
                     Nsamp = (apply(LF_S7, 1, sum)), LF_S7)
    Lfreq_S7 = Lfreq_S7_2[Lfreq_S7_2[,6] > 0,]
    
    
    #####
    lencomp.dat = rbind(Lfreq_C1, Lfreq_C2, Lfreq_C3, Lfreq_C4, Lfreq_S1, Lfreq_S2, Lfreq_S3, Lfreq_S4, Lfreq_S5, Lfreq_S6, Lfreq_S7)
    # write.csv(lencomp.dat, file="N:\\Documents\\DFA_Simulation\\SB\\Build_SS_mod\\lencomp.csv")
    new_dat$lencomp = lencomp.dat
    new_dat$N_lencomp = nrow(lencomp.dat)
    
    
    # for DFA assessment
    FL = DFA_Results$FLoadings[i,1:7]
    lambda = ifelse(FL<0, 0, FL/max(FL))
    
    
    # L1 = round(lambda[1]*LF_S1)
    # L2 = round(lambda[2]*LF_S2)
    # L3 = round(lambda[3]*LF_S3)
    # L4 = round(lambda[4]*LF_S4)
    # L5 = round(lambda[5]*LF_S5)
    # L6 = round(lambda[6]*LF_S6)
    # L7 = round(lambda[7]*LF_S7) 
    # 
    # L1 = cbind(Yr=Y1, L1)
    # L2 = cbind(Yr=Y2, L2)
    # L3 = cbind(Yr=Y3, L3)
    # L4 = cbind(Yr=Y4, L4)
    # L5 = cbind(Yr=Y5, L5)
    # L6 = cbind(Yr=Y6, L6)
    # L7 = cbind(Yr=Y7, L7)
    # 
    # rbind(L1, L2, L3, L4, L5, L6, L7)
    # d = rbind(L1, L2)
    # plot(as.vector(summaryBy(d~Yr, d, FUN=sum, keep.names=TRUE)[5,-2]) , as.vector(L2[2,-1] + L1[4,-1]) )
    # rbind(L1, L2, L3, L4, L5, L6, L7)
    # 
    # rownames(LF_S1) = Y1
    # rownames(LF_S2) = Y2
    # rownames(LF_S3) = Y3
    # rownames(LF_S4) = Y4
    # rownames(LF_S5) = Y5
    # rownames(LF_S6) = Y6
    # rownames(LF_S7) = Y7
    # xy = rbind(x,y)
    # z = rowsum(xy,group=rownames(xy))
    # z[order(as.numeric(row.names(z))),]
    LF1 = rbind(round(lambda[1]*LF_S1),round(lambda[2]*LF_S2),round(lambda[3]*LF_S3),round(lambda[4]*LF_S4),round(lambda[5]*LF_S5),
                round(lambda[6]*LF_S6),round(lambda[7]*LF_S7))
    # LF1 = cbind( Yr = c(Y1, Y2, Y3, Y4, Y5, Y6, Y7), LF1  )
    LF_DFA = rowsum(LF1, group=c(Y1, Y2, Y3, Y4, Y5, Y6, Y7))
    
    
    # LF_DFA = round(as.numeric(lambda[1])*LF_S1) + round((as.numeric(lambda[2])*LF_S2)) + round(as.numeric(lambda[3])*LF_S3) + 
    #   round((as.numeric(lambda[4])*LF_S4)) + round(as.numeric(lambda[5])*LF_S5) + round((as.numeric(lambda[6])*LF_S6)) + 
    #   round(as.numeric(lambda[7])*LF_S7)
    # LF_DFA = LF_I1+LF_I2+LF_I3
    Lfreq_DFA1 = cbind(year=as.numeric(rownames(LF_DFA)), seas = (rep(1, length(rownames(LF_DFA)))), FltSry = (rep(5, length(rownames(LF_DFA)))), 
                      Sex = (rep(3, length(rownames(LF_DFA)))), Part = (rep(0, length(rownames(LF_DFA)))), 
                      Nsamp = (apply(LF_DFA, 1, sum)), LF_DFA)
    Lfreq_DFA = Lfreq_DFA1[Lfreq_DFA1[,6] > 0,]
    
    lencomp_DFA.dat = rbind(Lfreq_C1, Lfreq_C2, Lfreq_C3, Lfreq_C4, Lfreq_DFA)
    # write.csv(lencomp_DFA.dat, file="N:\\Documents\\DFA_Simulation\\SB\\Build_SS_mod\\DFA_SS\\lencomp_DFA.dat.csv")
    
    new_DFA_dat$lencomp = lencomp_DFA.dat
    new_DFA_dat$N_lencomp = nrow(lencomp_DFA.dat)
  
    
    
    SS_writedat(new_dat, outfile = paste(dir, trial, "/SS/SB_SS_", i, "/sandbar.dat", sep=""), overwrite=TRUE, version="3.24")
    SS_writedat(new_DFA_dat, outfile = paste(dir, trial, "/SS_DFA/SB_SS_", i, "/sandbar.dat", sep=""), overwrite=TRUE, version="3.24")
    
    
    } # END i loop
    
  } # END if missing==TRUE
  
  
  
} # end datfile function



dir1="D:/vspace1/DFA_Simulation/SB/"
base_dat_1 = SS_readdat(file="D:\\vspace1\\DFA_Simulation\\SB\\Trial138\\SS\\SB_SS_1\\sandbar.dat")
base_dat_DFA_1 = SS_readdat(file="D:\\vspace1\\DFA_Simulation\\SB\\Trial138\\SS_DFA\\SB_SS_1\\sandbar.dat")
datfile_SS(dir=dir1, trial="Trial138", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)


datfile_SS(dir=dir1, trial="Trial1", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial2", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial3", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial4", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial5", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial6", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial7", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial8", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial9", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial10", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial11", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial12", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial13", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial14", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial15", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial16", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial17", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial18", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial19", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial20", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial21", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial22", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial23", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial24", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial25", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial26", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial27", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial28", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial29", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial30", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial31", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial32", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial33", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial34", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial35", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial36", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial37", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial38", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)
datfile_SS(dir=dir1, trial="Trial39", base_dat=base_dat_M_1, base_DFA_dat=base_dat_DFAM_1,
           Niters=100, missing=TRUE)



datfile_SS(dir=dir1, trial="Trial101", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial102", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial103", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial104", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial105", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial106", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial107", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial108", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial109", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial110", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial111", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial112", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial113", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial114", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial115", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial116", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial117", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial118", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial119", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial120", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial121", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial122", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial123", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial124", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial125", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial126", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial127", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial128", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial129", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial130", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial131", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial132", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial133", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial134", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial135", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial136", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial137", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial138", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
datfile_SS(dir=dir1, trial="Trial139", base_dat=base_dat_1, base_DFA_dat=base_dat_DFA_1,
           Niters=100, missing=FALSE)
