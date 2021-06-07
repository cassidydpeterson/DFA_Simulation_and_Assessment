#####################################
# Atlantic sharpnose shark simulation
# DFA Simulation/Management Test
# WORKING VERSION: Updated Mar 2020
#####################################


# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/ATL_SN_FEMALE_SurvivalSR_LogErr_DFA_Simulation_workspace.R.RData")

# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/ATL_SN_DFA_Simulation_FEMALEONLY_SurvivalSR_workspace.R.RData")
# load("Y:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/ATL_SN_DFA_Simulation_FEMALEONLY_SurvivalSR_workspace.R.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/SN_BothSexes_Workspace.RData")
# load("C:/Users/cpeterson/Documents/DFA_Simulation/Atl_SN/SN_BothSexes_Workspace.RData")

library(scales)
library(MARSS)
library(vioplot)

###########################
# from: http://blog.revolutionanalytics.com/2016/08/simulating-form-the-bivariate-normal-distribution-in-r-1.html
rbvn<-function (n, m1, s1, m2, s2, rho) {
  X1 <- rnorm(n, m1, s1)
  X2 <- rnorm(n, m2 + (s2/s1) * rho * 
                (X1 - m1), sqrt((1 - rho^2)*s2^2))
  data.frame(Linf = X1, K = X2)
}

Linf=81.6
K=0.49
t0=-0.97
# vonBert = rbvn(nrow(df.expanded),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
ages = 1:23
# vonBert$L0 = 

# VON BERT:  Lt = Linf *(1-exp(-K*(t-t0)))
# SCHNUTE: L1 + (L2 - L1) * ((1-exp(-K*(t-t1)))/(1-exp(-K*(t2-t1))))
L0 = Linf * (1-exp(-K*(0-t0))) # L0 = 30.8694
L1 = Linf * (1-exp(-K*(1-t0))) # L1 = 50.5211



# AGE-BASED SELECTIVITY:
Sel_catch = read.csv("D:\\vspace1\\DFA_Simulation\\Atl_SN\\CatchSelectivity.csv")
Sel_survey = read.csv("D:\\vspace1\\DFA_Simulation\\Atl_SN\\SurveySelectivity.csv")
# Sel_catch = read.csv("N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\CatchSelectivity.csv")
# Sel_survey = read.csv("N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\SurveySelectivity.csv")
# Sel_catch = read.csv("C:\\Users\\cpeterson\\Documents\\DFA_Simulation\\Atl_SN\\CatchSelectivity.csv")
# Sel_survey = read.csv("C:\\Users\\cpeterson\\Documents\\DFA_Simulation\\Atl_SN\\SurveySelectivity.csv")

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


#age structure
ages=1:18
A=max(ages)

#initial population
No=1000

#length of simulation
yrs=1000


# #density dependent pup mortality 
# 
# a=0.012
# b=0.006
# 
# a/b
# 
# Npups=0:5000
# Mo=(a*Npups)/(1+b*Npups)
# Npups_eq = 3320.671
# 
# par(mfrow=c(1,1))
# plot(Npups,Mo,type='l',lwd=2,xlab="Number of 1-yr-old Pups (nos.)",ylab="Instantaneous Mortality of Newborn Pups")
# plot(Npups, exp(-Mo), type='l',lwd=2,xlab="Number of 1-yr-old Pups (nos.)",ylab="Survival of Newborn Pups")
# 
# plot(Npups/Npups_eq, Mo,type='l',lwd=2,xlab="Number of 1-yr-old Pups (nos.)",ylab="Instantaneous Mortality of Newborn Pups")
# plot(Npups/Npups_eq, exp(-Mo), type='l',lwd=2,xlab="Number of 1-yr-old Pups (nos.)",ylab="Survival of Newborn Pups")
# 
# # equilibrium Npups = 3320.671

#F= fishing mortality (estimated Fmsy=0.06)
FM=0.00

iterations=1000




#ages
ages = 1:18
A = max(ages)
# consider plus group to age 23
# Define maturity function 
mat = as.vector(c(0.185, 0.953, 0.999, rep(1,15)))

# Define fecundity (female pups)
fec = c(0.401, 0.762, 1.133, 1.448, 1.686, 1.852, 1.963, 2.035, 2.081, 2.11, 2.128, 2.139, 2.146, 2.15, 2.153, 2.155, 2.156, 2.156)
#fec=fec*2



# length of simulation
years = 2000


# Fishing mortality will vary over time. For now, set it to a constant 0.1
# FM = 0.3

N0 = 505

################## CHANGE INPUTS ##################
#__________________________________________________


################# define natural mortality
M_const = as.vector(rep(0.209, A))

################# Assume constant fishing mortality
F_const = as.vector(rep(0.0,years))

# Store results



#################### Run SImulation  ##########
# Create matrices to store results
Nay = matrix(nrow=A, ncol=years)
Cay = matrix(nrow=A, ncol=years)
M0 = vector(length=years)
Npups = vector(length=years)

# Inputs: 
M = M_const
FM = F_const


Nay[1,1] = N0

# populate first row
for (i in 1:(A-1)){
  Nay[i+1,1] = Nay[i,1]*exp(-M[i])
}

# a=0.012 # implement as Monte Carlo appraoch.
# b=0.006
Z0 = 1.904416
Zmin = 0.209
beta = 0.58
Npups_eq = 3320.671

for (i in 1:(years-1)){

  Npups[i] = sum(Nay[,i]*mat*fec) # equilibrium Npups = 3320.671
  
  M0[i] = ( (1- (Npups[i]/Npups_eq)^beta ) * (Zmin - Z0) ) + Z0
  
  # M0[i] = (a*Npups[i])/(1+b*Npups[i])
  # M0[i] = 1.9044160 # THIS IS EQUILIBRIUM!!!
  Nay[1,i+1] = Npups[i]*exp(-M0[i])  # pups produced that survive to age 1; Nay[1,] equilibrium = 494.4797
  
  # Cay[j,1]=Nay[j,1]*((Sel_catch[j,2]*FM[1]/4)/(Sel_catch[j,2]*FM[1]/4+M[j]))*(1-exp(-(Sel_catch[j,2]*FM[1]/4+M[j]))) + 
  #   Nay[j,1]*((Sel_catch[j,3]*FM[1]/4)/(Sel_catch[j,3]*FM[1]/4+M[j]))*(1-exp(-(Sel_catch[j,3]*FM[1]/4+M[j]))) +
  #   Nay[j,1]*((Sel_catch[j,4]*FM[1]/4)/(Sel_catch[j,4]*FM[1]/4+M[j]))*(1-exp(-(Sel_catch[j,4]*FM[1]/4+M[j]))) +
  #   Nay[j,1]*((Sel_catch[j,5]*FM[1]/4)/(Sel_catch[j,5]*FM[1]/4+M[j]))*(1-exp(-(Sel_catch[j,5]*FM[1]/4+M[j])))
  
  for (j in 1:(A-2)){
    Nay[j+1,i+1] = Nay[j,i]*exp(-M[j] - (Sel_catch[j,2]*FM[i]/4) - (Sel_catch[j,3]*FM[i]/4) - (Sel_catch[j,4]*FM[i]/4) - (Sel_catch[j,5]*FM[i]/4)) 
    # Say[j,i]=Nay[j,i]*mat[j]
    # Pay[j,i]=Say[j,i]*fec[j]
    Cay[j,i+1]=Nay[j,i]*((Sel_catch[j,2]*FM[i]/4)/(Sel_catch[j,2]*FM[i]/4+M[j]))*(1-exp(-(Sel_catch[j,2]*FM[i]/4+M[j]))) + 
      Nay[j,i]*((Sel_catch[j,3]*FM[i]/4)/(Sel_catch[j,3]*FM[i]/4+M[j]))*(1-exp(-(Sel_catch[j,3]*FM[i]/4+M[j]))) +
      Nay[j,i]*((Sel_catch[j,4]*FM[i]/4)/(Sel_catch[j,4]*FM[i]/4+M[j]))*(1-exp(-(Sel_catch[j,4]*FM[i]/4+M[j]))) +
      Nay[j,i]*((Sel_catch[j,5]*FM[i]/4)/(Sel_catch[j,5]*FM[i]/4+M[j]))*(1-exp(-(Sel_catch[j,5]*FM[i]/4+M[j])))
    
  }
  Nay[A,i+1] = Nay[A-1,i]*exp(-M[A-1] - (Sel_catch[A-1,2]*FM[i]/4) - (Sel_catch[A-1,3]*FM[i]/4) - (Sel_catch[A-1,4]*FM[i]/4) - (Sel_catch[A-1,5]*FM[i]/4)) + 
    Nay[A,i]*exp(-M[A] - (Sel_catch[A,2]*FM[i]/4) - (Sel_catch[A,3]*FM[i]/4) - (Sel_catch[A,4]*FM[i]/4) - (Sel_catch[A,5]*FM[i]/4))
  
  Cay[A-1,i]=Nay[A,i]*((Sel_catch[A-1,2]*FM[i]/4)/(Sel_catch[A-1,2]*FM[i]/4+M[j]))*(1-exp(-(Sel_catch[A-1,2]*FM[i]/4+M[A-1]))) +
    Nay[A-1,i]*((Sel_catch[A-1,3]*FM[i]/4)/(Sel_catch[A-1,3]*FM[i]/4+M[A-1]))*(1-exp(-(Sel_catch[A-1,3]*FM[i]/4+M[A-1]))) +
    Nay[A-1,i]*((Sel_catch[A-1,4]*FM[i]/4)/(Sel_catch[A-1,4]*FM[i]/4+M[A-1]))*(1-exp(-(Sel_catch[A-1,4]*FM[i]/4+M[A-1]))) +
    Nay[A-1,i]*((Sel_catch[A-1,5]*FM[i]/4)/(Sel_catch[A-1,5]*FM[i]/4+M[A-1]))*(1-exp(-(Sel_catch[A-1,5]*FM[i]/4+M[A-1])))
  
  Cay[A,i]=Nay[A,i]*((Sel_catch[A,2]*FM[i]/4)/(Sel_catch[A,2]*FM[i]/4+M[j]))*(1-exp(-(Sel_catch[A,2]*FM[i]/4+M[A]))) +
    Nay[A,i]*((Sel_catch[A,3]*FM[i]/4)/(Sel_catch[A,3]*FM[i]/4+M[A]))*(1-exp(-(Sel_catch[A,3]*FM[i]/4+M[A]))) +
    Nay[A,i]*((Sel_catch[A,4]*FM[i]/4)/(Sel_catch[A,4]*FM[i]/4+M[A]))*(1-exp(-(Sel_catch[A,4]*FM[i]/4+M[A]))) +
    Nay[A,i]*((Sel_catch[A,5]*FM[i]/4)/(Sel_catch[A,5]*FM[i]/4+M[A]))*(1-exp(-(Sel_catch[A,5]*FM[i]/4+M[A])))
}



plot(1:years, colSums(Nay),type='l', ylim=c(0, 2750))#, ylim=c(0,max(colSums(Nay))+100))
plot(1:years, Npups,type='l')#, ylim=c(0,max(colSums(Nay))+100))
tail(colSums(Nay))


# Equilibrium conditions:
Nay[,1990:2000]
Equilibrium_Nay = Nay[,2000] 
# Equilibrium_Nay = c(494.47970, 401.21848, 325.54676, 264.14709, 214.32769, 173.90447, 141.10526, 114.49213,  92.89837,  75.37730,  
#                   61.16078,  49.62556,  40.26595,  32.67160, 26.50958,  21.50975,  17.45290,  75.08402)


####### FOR STOCK SYNTHESIS #####

z0 = 1.904416
s0 = exp(-z0); s0 #0.1489096
# zmin = 0.06157842
zmin = 0.209
(smax = exp(-zmin)) #0.8113952
( zfrac = (log(smax) - log(s0))/(-log(s0)) ) #0.8902551
h=0.56
beta = (log(1-(log(h/0.2)/(z0*zfrac))))/log(0.2)
beta
# zmax = 2.0
( S_frac = (zmin - z0)/(0-z0) ) #0.8902551

Z_max = z0 + S_frac*(0.0-z0) #0.209



### ALTER INPUTS #### 
# length of simulation
years = 100


# Fishing mortality will vary over time. For now, set it to a constant 0.1

# define natural mortality
M_const = as.vector(rep(0.209, A))

# Assume constant fishing mortality
F_var = as.vector(rep(0.0,years))
F_var = as.vector(rep(0.2,years)) # Trials 1-6
# F_var = as.vector(c(rep(0.1, 20), rep(0.2, 20), rep(0.3, 20), rep(0.4, 20), rep(0.2, 20)))
# F_var = as.vector(c(rep(0.1, 50), rep(0.2, 50)))
# F_var = as.vector(c(rep(0.4, 50), rep(0.2, 50)))
F_var = as.vector(c(rep(0.4, 50), rep(0, 50))) # I like this one. Trials 13-18
F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))) # I like this one. Trials 7-12
F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))) # Trials 19-24
# F_var = as.vector(c(rep(0,30), rep(0.4, 20), rep(0.2, 5), rep(0.05, 45))) # 
# F_var = as.vector(c(rep(0.0, 50), rep(0.3, 50)))
# F_var = as.vector(c(rep(0.0, 50), rep(0.2, 50)))
# F_var = as.vector(c(rep(0.6, 25), rep(0.2, 50), rep(0.6, 25)))
# F_var = as.vector(c(rep(0.1, 50), rep(0.2, 50)))
# F_var = as.vector(c(rep(0.05, 50), rep(0.1, 50))) # BASE CASE
# F_var = as.vector(c(rep(0.05, 50), rep(0.0, 50)))

# F_var = as.vector(rep(1,years)) # Test
# F_var = as.vector(c(rep(0.4, 100)))
#### INITIALIZE EQUILIBRIUM CONDITIONS ####
# Create matrices to store results
Nay = matrix(nrow=A, ncol=years)
Cay = matrix(nrow=A, ncol=years)
Iay = matrix(nrow=A, ncol=years)
IayN = matrix(nrow=A, ncol=years)
Iay1 = matrix(nrow=A, ncol=years)
Iay2 = matrix(nrow=A, ncol=years)
vulnerability_a = c(0.5, rep(0.99, length=(A-1)))
vulnerability_a1 = c(0.5, rep(0.99, length=(A-1)))
vulnerability_a2 = c(0.2, 0.6, 0.9, rep(1, length=(A-6)), 0.9, 0.8, 0.6)
q_survey=0.01
CV=0.5
CV1=0.5
CV2=0.5
CV3=0.5

q_survey1 = rep(0.01, length=years)
q_survey2 = rep(0.01, length=years)
q_survey3 = rep(0.01, length=years)
q_survey = c(rep(0.00025, 49), seq(0.0004, 0.0014, by=0.0001), rep(0.0015, 40))
q_survey = c(rep(0.002, 49), seq(0.0018, 0.0006, by=-0.0002), rep(0.0005, 44))


M0 = vector(length=years)
Npups = vector(length=years)
# REMEMBER: S (annual survival) is related to exp(-Z), where Z is total mortality. So, M0 ~= Z0 : S0 ~= exp(-M0). 

# Inputs: 
M = M_const
# FM = F_const
FM = F_var

# Set up equilibrium conditions
Nay[,1] = Equilibrium_Nay
Z0 = 1.904416
Zmin = 0.209
Npups_eq = 3320.671
beta = 0.5807613

  
for (i in 1:(years-1)){
  
  Npups[i] = sum(Nay[,i]*mat*fec) # equilibrium Npups = 3320.671
  M0[i] = ( (1- (Npups[i]/Npups_eq)^beta ) * (Zmin - Z0) ) + Z0 
  # M0[i] = 1.9044160 # THIS IS EQUILIBRIUM!!!
  Nay[1,i+1] = Npups[i]*exp(-M0[i])  # pups produced that survive to age 1
  
  
  for (j in 1:(A-2)){
    Nay[j+1,i+1] = Nay[j,i]*exp(-M[j] - (Sel_catch[j,2]*FM[i]/4) - (Sel_catch[j,3]*FM[i]/4) - (Sel_catch[j,4]*FM[i]/4) - (Sel_catch[j,5]*FM[i]/4)) 
    # Say[j,i]=Nay[j,i]*mat[j]
    # Pay[j,i]=Say[j,i]*fec[j]
    Cay[j,i]=Nay[j,i]*((Sel_catch[j,2]*FM[i]/4)/(Sel_catch[j,2]*FM[i]/4+M[j]))*(1-exp(-(Sel_catch[j,2]*FM[i]/4+M[j]))) + 
      Nay[j,i]*((Sel_catch[j,3]*FM[i]/4)/(Sel_catch[j,3]*FM[i]/4+M[j]))*(1-exp(-(Sel_catch[j,3]*FM[i]/4+M[j]))) +
      Nay[j,i]*((Sel_catch[j,4]*FM[i]/4)/(Sel_catch[j,4]*FM[i]/4+M[j]))*(1-exp(-(Sel_catch[j,4]*FM[i]/4+M[j]))) +
      Nay[j,i]*((Sel_catch[j,5]*FM[i]/4)/(Sel_catch[j,5]*FM[i]/4+M[j]))*(1-exp(-(Sel_catch[j,5]*FM[i]/4+M[j])))
    
    # Iay[j,i] = q_survey * vulnerability_a[j] * Nay[j,i] +  rnorm(1, 0, q_survey * vulnerability_a[j] * Nay[j,i]/CV) #LOTS OF UNCERTAINTY!based off CV; maybe randomly sample CV from a normal distribution to add more uncertainty?
    
  }
  Nay[A,i+1] = Nay[A-1,i]*exp(-M[A-1] - (Sel_catch[A-1,2]*FM[i]/4) - (Sel_catch[A-1,3]*FM[i]/4) - (Sel_catch[A-1,4]*FM[i]/4) - (Sel_catch[A-1,5]*FM[i]/4)) + 
    Nay[A,i]*exp(-M[A] - (Sel_catch[A,2]*FM[i]/4) - (Sel_catch[A,3]*FM[i]/4) - (Sel_catch[A,4]*FM[i]/4) - (Sel_catch[A,5]*FM[i]/4))
  
  Cay[A-1,i]=Nay[A,i]*((Sel_catch[A-1,2]*FM[i]/4)/(Sel_catch[A-1,2]*FM[i]/4+M[j]))*(1-exp(-(Sel_catch[A-1,2]*FM[i]/4+M[A-1]))) +
    Nay[A-1,i]*((Sel_catch[A-1,3]*FM[i]/4)/(Sel_catch[A-1,3]*FM[i]/4+M[A-1]))*(1-exp(-(Sel_catch[A-1,3]*FM[i]/4+M[A-1]))) +
    Nay[A-1,i]*((Sel_catch[A-1,4]*FM[i]/4)/(Sel_catch[A-1,4]*FM[i]/4+M[A-1]))*(1-exp(-(Sel_catch[A-1,4]*FM[i]/4+M[A-1]))) +
    Nay[A-1,i]*((Sel_catch[A-1,5]*FM[i]/4)/(Sel_catch[A-1,5]*FM[i]/4+M[A-1]))*(1-exp(-(Sel_catch[A-1,5]*FM[i]/4+M[A-1])))
  
  Cay[A,i]=Nay[A,i]*((Sel_catch[A,2]*FM[i]/4)/(Sel_catch[A,2]*FM[i]/4+M[j]))*(1-exp(-(Sel_catch[A,2]*FM[i]/4+M[A]))) +
    Nay[A,i]*((Sel_catch[A,3]*FM[i]/4)/(Sel_catch[A,3]*FM[i]/4+M[A]))*(1-exp(-(Sel_catch[A,3]*FM[i]/4+M[A]))) +
    Nay[A,i]*((Sel_catch[A,4]*FM[i]/4)/(Sel_catch[A,4]*FM[i]/4+M[A]))*(1-exp(-(Sel_catch[A,4]*FM[i]/4+M[A]))) +
    Nay[A,i]*((Sel_catch[A,5]*FM[i]/4)/(Sel_catch[A,5]*FM[i]/4+M[A]))*(1-exp(-(Sel_catch[A,5]*FM[i]/4+M[A])))
  
  for(j in 1:A){
    # Iay[j,i] =  q_survey1[i]*Sel_survey[h,2] * Nay[h,i]  
    
    IayN[j,i] = max(q_survey * Nay[j,i] +  rnorm(1, 0, (q_survey *  Nay[j,i])*CV),0) #LOTS OF UNCERTAINTY!based off CV;
    Iay[j,i] = max(q_survey * Nay[j,i] * exp(rnorm(1, 0, (q_survey *  Nay[j,i])*CV) - ((q_survey * Nay[j,i]*CV)^2/2) ), 0) #LOTS OF UNCERTAINTY!based off CV;
    
  }
  
  
}

Iay_eq = Iay


par(mfrow=c(1,1))
EquilibriumPop = colSums(Nay)
plot(1:years, EquilibriumPop,type='l', ylim=c(0,max(colSums(Nay))+1500))
lines(1:years, EquilibriumPop, type='l', col='blue')
abline(v=40, lwd=2, col="red")
abline(v=65, lwd=2, col="red")
# lines(1:years, colSums(Nay),type='l',col='black')

# Iay
par(mfrow=c(1,1))
plot(1:years, apply(Iay_eq, 2, sum), type='l',ylim=c(0,10), lwd=2)
lines(1:years, apply(IayN, 2, sum), type='l', col='red')
abline(v=40, lwd=2, col="red")
abline(v=65, lwd=2, col="red")
lines(1:years, EquilibriumPop/500, type='l', col='blue')
# plot(1:years, apply(Iay, 2, sum), type='l')
# Iy = apply(Iay, 2, sum)

# 
# plot(1:years, EquilibriumPop/150,type='l', ylim=c(0,20))
# lines(1:years, apply(Iay_eq, 2, sum), type='l')
# abline(v=51, lwd=2, col='red')

# Iay_eq_F1 = Iay_eq
# Iay_eq_F2 = Iay_eq
# Iay_eq_F3 = Iay_eq
# Iay_eq_F4 = Iay_eq
# 
# I_eq_list = list()
# I_eq_list[["Iay_eq_F1"]] = Iay_eq_F1
# I_eq_list[["Iay_eq_F2"]] = Iay_eq_F2
# I_eq_list[["Iay_eq_F3"]] = Iay_eq_F3
# I_eq_list[["Iay_eq_F4"]] = Iay_eq_F4
# save(I_eq_list, file="D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\I_eq_list.RData")
load("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\I_eq_list.RData")
load("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\I_eq_list.RData")

Iay_eq_F1 = I_eq_list$Iay_eq_F1
Iay_eq_F2 = I_eq_list$Iay_eq_F2
Iay_eq_F3 = I_eq_list$Iay_eq_F3
Iay_eq_F4 = I_eq_list$Iay_eq_F4


##### TRIALS ################
CV1=0.5
CV2=0.5
CV3=0.5

q_survey1 = rep(0.01, length=years)
q_survey2 = rep(0.01, length=years)
q_survey3 = rep(0.01, length=years)
q_survey = c(rep(0.00025, 49), seq(0.0004, 0.0014, by=0.0001), rep(0.0015, 40))
q_survey = c(rep(0.002, 49), seq(0.0018, 0.0006, by=-0.0002), rep(0.0005, 44))
# q_survey1 = c(rep(0.001, length=years/5), rep(0.0015, length=years/5),rep(0.002, length=years/5),rep(0.0025, length=years/5),rep(0.003, length=years/5))
# q_survey2 = c(rep(0.00005, length=years/2), rep(0.003, length=years/2))
# q_survey3 = c(rep(0.003, length=years/2), rep(0.00001, length=years/2))
# q_survey3 = c(rep(0.003, length=years))

# CHANGE WD
# setwd("Y:\\Documents\\DFA_Simulation\\Atl_SN\\Trial1")
# setwd("N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\Trial1")
# setwd(paste("N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\", paste("Trial",f,sep=""), sep=""))
##################################


############################################
# # Read in data
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial1\\ATrial1.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial2\\ATrial2.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial3\\ATrial3.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial4\\ATrial4.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial5\\ATrial5.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial6\\ATrial6.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial7\\ATrial7.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial8\\ATrial8.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial9\\ATrial9.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial10\\ATrial10.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial11\\ATrial11.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial12\\ATrial12.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial13\\ATrial13.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial14\\ATrial14.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial15\\ATrial15.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial16\\ATrial16.RData")
# load( file = "D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial17\\ATrial17.RData")


#############################



# w/ changing q
# when F4 --> c=c(7, 50, 2.5)
# when F1 --> c=c(1.2, 125, 1.5)
# when F2 --> c=c(10.5, 10, 10)
# when F3 --> c=c(20, 275, 3.5)


# w/ constant q
# when F1 -->   c=c(1.1, 7, 3)
# when F2 --> c=c(10.5, 10, 10)
# when F3 --> c=c(20.5, 19.5, 19.5)
# when F4 --> c=c(10, 10, 10)

# w/ gradual change q
# when F1 --> c(1.8, 50000, 20)
# when F2 --> c=c(500, 2.5, 800)
# when F3 --> c=c(40, 1000, 5)
# when F4 --> c=c(7, 80, 2.5)

#-------------------------------------------------------------
years = 100
c=5


F_var = as.vector(rep(0.2,years))
CV1=0.5
CV2=0.5
CV3=0.5
q_survey1 = rep(0.01, length=years)
q_survey2 = rep(0.01, length=years)
q_survey3 = rep(0.01, length=years)
Iay_eq = Iay_eq_F1

F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50)))
CV1=0.3
CV2=0.5
CV3=0.7
q_survey1 = rep(0.001, length=years)
q_survey2 = rep(0.001, length=years)
q_survey3 = rep(0.001, length=years)
Iay_eq = Iay_eq_F2

F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45)))
CV1=0.3
CV2=0.5
CV3=0.7
q_survey1 = rep(0.001, length=years)
q_survey2 = rep(0.001, length=years)
q_survey3 = rep(0.001, length=years)
Iay_eq = Iay_eq_F4



F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50)))
CV1=0.7
CV2=0.5
CV3=0.3
q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2))
q_survey2 = rep(0.001, length=years)
q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2))
Iay_eq = Iay_eq_F2


F_var = as.vector(c(rep(0.4, 50), rep(0, 50))) # I like this one. Trials 13-18
Iay_eq = Iay_eq_F3

#A10
F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50)))
CV1=0.5
CV2=0.5
CV3=0.5
q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2))
q_survey2 = rep(0.001, length=years)
q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2))
Iay_eq = Iay_eq_F2
c=100

q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35))
q_survey2 = rep(0.001, length=years)
q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40))

F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50)))
CV1=0.3
CV2=0.5
CV3=0.7
q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35))
q_survey2 = rep(0.001, length=years)
q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40))
Iay_eq = Iay_eq_F2
c=c(40, 1, 50)
F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50)))
CV1=0.7
CV2=0.5
CV3=0.3
q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35))
q_survey2 = rep(0.001, length=years)
q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40))
Iay_eq = Iay_eq_F2
c=c(500, 2.5, 800)
F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50)))
CV1=0.5
CV2=0.5
CV3=0.5
q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35))
q_survey2 = rep(0.001, length=years)
q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40))
Iay_eq = Iay_eq_F2
c=c(500, 2.5, 800)

F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45)))
CV1=0.3
CV2=0.5
CV3=0.7
q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2))
q_survey2 = rep(0.001, length=years)
q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2))
Iay_eq = Iay_eq_F4
c=c(7, 50, 2.5)



F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))) 
CV1=0.3
CV2=0.5
CV3=0.7
q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35))
q_survey2 = rep(0.001, length=years)
q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40))
Iay_eq = Iay_eq_F4
c=c(15, 250, 5)
reorg=T 



setwd("C:\\~")






########################## SIMULATION FUNCTION #################################################################################################
DFA_Sim = function(F_var, CV1, CV2, CV3, q_survey1, q_survey2, q_survey3, Iay_eq = Iay_eq, c=5, reorg=F){
  # Iay
  par(mfrow=c(1,1))
  plot(1:years, apply(Iay_eq, 2, sum), type='l',ylim=c(0,max(apply(Iay_eq, 2, sum),na.rm=T)+max(apply(Iay_eq, 2, sum), na.rm=T)*0.1))
  abline(v=40, lwd=2, col="red")
  abline(v=65, lwd=2, col="red")
  # plot(1:years, apply(Iay, 2, sum), type='l')
  # Iy = apply(Iay, 2, sum)
  
  F_var=F_var
  CV1=CV1
  CV2=CV2
  CV3=CV3
  q_survey1=q_survey1
  q_survey2=q_survey2
  q_survey3=q_survey3
  c=c
  c1=c
  
  ####### ADD STOCHASTICITY ###########
  Niters = 100
  Iay1 = matrix(nrow=A, ncol=years)
  Iay2 = matrix(nrow=A, ncol=years)
  Iay3 = matrix(nrow=A, ncol=years)
  Iay1_CV = matrix(nrow=A, ncol=years)
  Iay2_CV = matrix(nrow=A, ncol=years)
  Iay3_CV = matrix(nrow=A, ncol=years)
  Cy = matrix(nrow=Niters, ncol=years)
  Ny = matrix(nrow=Niters, ncol=years)
  Iy1 = matrix(nrow=Niters, ncol=years)
  Iy2 = matrix(nrow=Niters, ncol=years)
  Iy3 = matrix(nrow=Niters, ncol=years)
  Iy1_CV = matrix(nrow=Niters, ncol=years)
  Iy2_CV = matrix(nrow=Niters, ncol=years)
  Iy3_CV = matrix(nrow=Niters, ncol=years)
  M0 = matrix(nrow=Niters, ncol=years)
  Npups = matrix(nrow=Niters, ncol=years)
  FSS = matrix(nrow=Niters, ncol=years) # Female spawning stock 
  Z0 = 1.904416
  Zmin = 0.209
  Npups_eq = 3320.671
  beta=0.5807613
  ResultsList = list()
  
  set.seed(430)
  
  # F_const = rep(0, yrs)
  
  for(k in 1:Niters){
    # paste("Iteration",k,sep=" ")
    # Create matrices to store results
    Nay = matrix(nrow=A, ncol=years)
    Cay = matrix(nrow=A, ncol=years)
    
    # Inputs: 
    M = M_const
    # FM = F_const
    FM = F_var
    
    # Set up equilibrium conditions
    Nay[,1] = Equilibrium_Nay
    
    for (i in 1:(years-1)){
      
      matjit = rnorm(length(mat),mat,mat*0.01)
      matjit = ifelse(matjit>1, 1, matjit)
      FSS[k,i] = sum(Nay[,i]*matjit) #Fecund Stock size
      Npups[k,i] = sum(Nay[,i]*matjit*rnorm(length(fec), fec, fec*0.1)) # equilibrium Npups = 3320.671
      M0[k,i] = rnorm(1, ( ( (1- (Npups[k,i]/Npups_eq)^beta ) * (Zmin - Z0) ) + Z0 ) , 0.1)
      # M0[i] = 1.9044160 # THIS IS EQUILIBRIUM!!!
      Nay[1,i+1] = Npups[k,i]*exp(-M0[k,i])  # pups produced that survive to age 1 = Recruits!
      
      Mjit = rnorm(length(M), M, M*0.1) 
      
      Zay = vector(length=length(M))
      for (j in 1:A){
        Zay[j] = Mjit[j] + (Sel_catch[j,2]*FM[i]/4) + (Sel_catch[j,3]*FM[i]/4) + (Sel_catch[j,4]*FM[i]/4) + (Sel_catch[j,5]*FM[i]/4)
      } # end j loop Z
    
    
    for (j in 1:(A-2)){
      Nay[j+1,i+1] = Nay[j,i]*exp(-Zay[j])
      # Say[j,i]=Nay[j,i]*mat[j]
      # Pay[j,i]=Say[j,i]*fec[j]
      Cay[j,i]=Nay[j,i]*((Sel_catch[j,2]*FM[i]/4)/Zay[j])*(1-exp(-Zay[j])) +
        Nay[j,i]*((Sel_catch[j,3]*FM[i]/4)/Zay[j])*(1-exp(-Zay[j])) +
        Nay[j,i]*((Sel_catch[j,4]*FM[i]/4)/Zay[j])*(1-exp(-Zay[j])) +
        Nay[j,i]*((Sel_catch[j,5]*FM[i]/4)/Zay[j])*(1-exp(-Zay[j]))
      
    } # END J LOOP
    
    Nay[A,i+1] = Nay[A-1,i]*exp(-Zay[A-1]) +
      Nay[A,i]*exp(-Zay[A])
    
    Cay[A-1,i]=Nay[A-1,i]*((Sel_catch[A-1,2]*FM[i]/4)/Zay[A-1])*(1-exp(-Zay[A-1])) +
      Nay[A-1,i]*((Sel_catch[A-1,3]*FM[i]/4)/Zay[A-1])*(1-exp(-Zay[A-1])) +
      Nay[A-1,i]*((Sel_catch[A-1,4]*FM[i]/4)/Zay[A-1])*(1-exp(-Zay[A-1])) +
      Nay[A-1,i]*((Sel_catch[A-1,5]*FM[i]/4)/Zay[A-1])*(1-exp(-Zay[A-1]))
    
    Cay[A,i]=Nay[A,i]*((Sel_catch[A,2]*FM[i]/4)/Zay[A])*(1-exp(-Zay[A])) +
      Nay[A,i]*((Sel_catch[A,3]*FM[i]/4)/Zay[A])*(1-exp(-Zay[A])) +
      Nay[A,i]*((Sel_catch[A,4]*FM[i]/4)/Zay[A])*(1-exp(-Zay[A])) +
      Nay[A,i]*((Sel_catch[A,5]*FM[i]/4)/Zay[A])*(1-exp(-Zay[A]))
    
    
      
      for(h in 1:A){
        Iay1_CV[h,i] = runif(1, min=CV1-0.1, max=CV1+0.1)
        Iay1[h,i] = max(q_survey1[i]*Sel_survey[h,2] * Nay[h,i] * 
                          exp( rnorm(1, 0, (q_survey1[i]*Sel_survey[h,2]*Nay[h,i])*Iay1_CV[h,i]) - 
                                 ( ((q_survey1[i]*Sel_survey[h,2]*Nay[h,i]*Iay1_CV[h,i])^2)/2 ) ), 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay2_CV[h,i] = runif(1, min=CV2-0.1, max=CV2+0.1)
        # Iay2[h,i] = max(q_survey2[i]*Sel_survey[h,3] * Nay[h,i] +  rnorm(1, 0, (q_survey2[i]*Sel_survey[h,3]*Nay[h,i])*Iay2_CV[h,i] ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay2[h,i] = max(q_survey2[i]*Sel_survey[h,3] * Nay[h,i] * 
                          exp( rnorm(1, 0, (q_survey2[i]*Sel_survey[h,3]*Nay[h,i])*Iay2_CV[h,i]) - 
                                 ( ((q_survey2[i]*Sel_survey[h,3]*Nay[h,i]*Iay2_CV[h,i])^2)/2 ) ), 0) #LOTS OF UNCERTAINTY!based off CV;
        
        
        Iay3_CV[h,i] = runif(1, min=CV3-0.1, max=CV3+0.1)
        # Iay3[h,i] = max(q_survey3[i]*Sel_survey[h,2] * Nay[h,i] +  rnorm(1, 0, (q_survey3[i]*Sel_survey[h,2]*Nay[h,i])*Iay3_CV[h,i] ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay3[h,i] = max(q_survey3[i]*Sel_survey[h,2] * Nay[h,i] * 
                          exp( rnorm(1, 0, (q_survey3[i]*Sel_survey[h,2]*Nay[h,i])*Iay3_CV[h,i]) - 
                                 ( ((q_survey3[i]*Sel_survey[h,2]*Nay[h,i]*Iay3_CV[h,i])^2)/2 ) ), 0) #LOTS OF UNCERTAINTY!based off CV;
        
      } # end h
      
      # LCay1 = von bert curve to calculate length composition
      # random q: runif(1,max=0.005, min=0.0005)
    }# ends i loop
    
    # USE VON BERT CURVE TO ESTIMATE LENGTH COMP. Let von Bert parameters vary by Niter...???? 
    #   year (implement in i loop)... 
    #   Or should they vary by year class? Could create a matrix of von bert parameters based on age class 
    #       and pull out appropriate values as needed
    # Lt = Linf * (1 - exp(-K*(t - t0)))
    # For sexes combined.
    Linf=81.6
    K=0.49
    t0=-0.97
    
    # L0 = Linf * (1 - exp(-K*(0-t0)))

    # GET LENGTH FREQUENCIES
    Cay_N = cbind(ages, round(Cay))
    assign(paste("LF_Cay",k,sep="_"), data.frame(LBins = seq(20, 122, by=2)))
    for(f in 40:65){
      df = as.data.frame(Cay_N[,c(1,f+1)])
      colnames(df) = c("ages","count") 
      df.expanded = df[rep(row.names(df), df$count),1:2]
      vonBert = rbvn(nrow(df.expanded),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expanded$L = vonBert$Linf * ( 1 - exp(-(vonBert$K) * (jitter(df.expanded$ages,amount=0.5) - t0) ))
      hst=hist(df.expanded$L, breaks=c(seq(20, 124, by=2)), plot=F)
      assign(paste("LF_Cay",k,sep="_"), cbind(get(paste("LF_Cay",k,sep="_")), hst$counts))
    } # end f loop
    # write.csv(get(paste("LF_Cay",k,sep="_")), paste(paste("LF_Cay",k,sep="_"),".csv",sep=""))
    
    Iay1_N = cbind(ages, round(Iay1*100))
    assign(paste("LF_Iay1",k,sep="_"), data.frame(LBins = seq(20, 122, by=2)))
    for(l in 40:65){
      df = as.data.frame(Iay1_N[,c(1,l+1)])
      colnames(df) = c("ages","count") 
      df.expanded = df[rep(row.names(df), df$count),1:2]
      vonBert = rbvn(nrow(df.expanded),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expanded$L = vonBert$Linf * ( 1 - exp(-(vonBert$K) * (jitter(df.expanded$ages,amount=0.5) - t0) ))
      hst=hist(df.expanded$L, breaks=c(seq(20, 124, by=2)), plot=F)
      assign(paste("LF_Iay1",k,sep="_"), cbind(get(paste("LF_Iay1",k,sep="_")), hst$counts))
    } # end l loop
    # write.csv(get(paste("LF_Iay1",k,sep="_")), paste(paste("LF_Iay1",k,sep="_"),".csv",sep=""))
    
    Iay2_N = cbind(ages, round(Iay2*100))
    assign(paste("LF_Iay2",k,sep="_"), data.frame(LBins = seq(20, 122, by=2)))
    for(m in 40:65){
      df = as.data.frame(Iay2_N[,c(1,m+1)])
      colnames(df) = c("ages","count") 
      df.expanded = df[rep(row.names(df), df$count),1:2]
      vonBert = rbvn(nrow(df.expanded),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expanded$L = vonBert$Linf * ( 1 - exp(-(vonBert$K) * (jitter(df.expanded$ages,amount=0.5) - t0) ))
      hst=hist(df.expanded$L, breaks=c(seq(20, 124, by=2)), plot=F)
      assign(paste("LF_Iay2",k,sep="_"), cbind(get(paste("LF_Iay2",k,sep="_")), hst$counts))
    } # end m loop
    # write.csv(get(paste("LF_Iay2",k,sep="_")), paste(paste("LF_Iay2",k,sep="_"),".csv",sep=""))
    
    
    Iay3_N = cbind(ages, round(Iay3*100))
    assign(paste("LF_Iay3",k,sep="_"), data.frame(LBins = seq(20, 122, by=2)))
    for(n in 40:65){
      df = as.data.frame(Iay3_N[,c(1,n+1)])
      colnames(df) = c("ages","count") 
      df.expanded = df[rep(row.names(df), df$count),1:2]
      vonBert = rbvn(nrow(df.expanded),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expanded$L = vonBert$Linf * ( 1 - exp(-(vonBert$K) * (jitter(df.expanded$ages,amount=0.5) - t0) ))
      hst=hist(df.expanded$L, breaks=c(seq(20, 124, by=2)), plot=F)
      assign(paste("LF_Iay3",k,sep="_"), cbind(get(paste("LF_Iay3",k,sep="_")), hst$counts))
    } # end n loop
    # write.csv(get(paste("LF_Iay3",k,sep="_")), paste(paste("LF_Iay3",k,sep="_"),".csv",sep=""))
    
    #save length-frequencies
    LFreqList = list()
    LFreqList[[paste("LF_Cay",k,sep="_")]] = get(paste("LF_Cay",k,sep="_"))
    LFreqList[[paste("LF_Iay1",k,sep="_")]] = get(paste("LF_Iay1",k,sep="_"))
    LFreqList[[paste("LF_Iay2",k,sep="_")]] = get(paste("LF_Iay2",k,sep="_"))
    LFreqList[[paste("LF_Iay3",k,sep="_")]] = get(paste("LF_Iay3",k,sep="_"))
    
    save(LFreqList, file=paste("LFreqList_",k,".RData",sep=""))
    
    
    # Linf = rnorm(years, 81.9, 81.9*0.1) # CV = 0.2; 0.2 = sd / 81.9
    # K = rnorm(years, 0.48, 0.48*0.1)
    # t0 = rnorm(years, -0.99, 0.99*0.1)
    # LC = Linf * ( 1 - exp(-K*(c(0,ages) - t0)))
    # assign(paste("LC",k,sep="_"), cbind(c(0,ages), LC))
    # write.csv(LC, paste(paste("LC",k,sep="_"),".csv",sep=""))
    # assign(paste("M0",k,sep="_"), M0)
    
    
    
    Iy1[k,] = apply(Iay1, 2, sum)
    Iy2[k,] = apply(Iay2, 2, sum)
    Iy3[k,] = apply(Iay3, 2, sum)
    lines(1:years, Iy1[k,], type='l', col='grey')
    lines(1:years, Iy2[k,], type='l', col='blue')
    lines(1:years, Iy3[k,], type='l', col='red')
    Iy1_CV[k,] = apply(Iay1_CV, 2, mean, na.rm=T)
    Iy2_CV[k,] = apply(Iay2_CV, 2, mean, na.rm=T)
    Iy3_CV[k,] = apply(Iay3_CV, 2, mean, na.rm=T)
    
    Cy[k,] = apply(Cay, 2, sum)
    Ny[k,] = apply(Nay, 2, sum)
    

    
    #write 
    # write.csv(Nay, paste(paste("Nay", k, sep="_"),".csv",sep=""))
    assign(paste("Nay", k, sep="_"), Nay)
    # write.csv(Cay, paste(paste("Cay", k, sep="_"),".csv",sep=""))
    assign(paste("Cay", k, sep="_"), Cay)
    # write.csv() ########### WRITE CSVs to SAVE Nay AND Cay!!!
    
    #write 
    AnnualNC = list()
    AnnualNC[[paste("Nay",k,sep="_")]] = Nay
    AnnualNC[[paste("Cay",k,sep="_")]] = Cay
    AnnualNC[[paste("Iay1",k,sep="_")]] = Iay1
    AnnualNC[[paste("Iay2",k,sep="_")]] = Iay2
    AnnualNC[[paste("Iay3",k,sep="_")]] = Iay3
    
    save(AnnualNC, file=paste("AnnualNC_",k,".RData",sep=""))
    
    
    
  } # end k loop
  
  
  
  ResultsList[["Iy1"]] = Iy1
  ResultsList[["Iy2"]] = Iy2
  ResultsList[["Iy3"]] = Iy3
  
  ResultsList[["Iy1_CV"]] = Iy1_CV
  ResultsList[["Iy2_CV"]] = Iy2_CV
  ResultsList[["Iy3_CV"]] = Iy3_CV
  
  
  ResultsList[["M0"]] = M0
  ResultsList[["FSS"]] = FSS
  ResultsList[["Npups"]] = Npups
  
  ResultsList[["Cy"]] = Cy
  ResultsList[["Ny"]] = Ny
  
  save(ResultsList, file="ResultsList.RData")
  
  
  
  ########## TEST DFA #############
  
  # library(MARSS)
  yrs = 40:(years-35)
  biasAll = matrix(nrow=Niters, ncol=length(yrs))
  bias = vector(length=Niters)
  RMSE = vector(length=Niters)
  biasAll.BT = matrix(nrow=Niters, ncol=length(yrs))
  RMSE.BT = vector(length=Niters)
  
  FitRatios = matrix(nrow=Niters, ncol=3)
  FactorLoadings=vector()
  DFATrends = vector()
  DFATrendsBT = vector()
  DFATrendsSE = vector()
  DFATrendsSEBT = vector()
  upCI_DFATrends = vector()
  lowCI_DFATrends = vector()
  
  datz_SD = matrix(nrow=Niters, ncol=3)
  # i=1
  
  
  for(i in 1:Niters){
    # assign(paste("dat",i,sep=""), rbind(Iy1[i,40:(years-35)], Iy2[i,40:(years-35)], Iy3[i,40:(years-35)]))
    # assign(paste("dat",i,sep=""), rbind(Iy2[i,40:(years-35)], Iy1[i,40:(years-35)], Iy3[i,40:(years-35)]))
    # assign(paste("dat",i,sep=""), rbind(Iy2[i,40:(years-35)], Iy3[i,40:(years-35)], Iy1[i,40:(years-35)]))
    
    if(reorg == T){
      assign(paste("dat",i,sep=""), rbind(Iy2[i,40:(years-35)], Iy3[i,40:(years-35)], Iy1[i,40:(years-35)]))
      c=c(  c1[1], c1[3], c1[2])
      
    }
    
    if(reorg==F){
      assign(paste("dat",i,sep=""), rbind(Iy2[i,40:(years-35)], Iy1[i,40:(years-35)], Iy3[i,40:(years-35)]))
    }
    
    dat.a = get(paste("dat",i,sep=""))
    # c=c(10, 110, 7)
    
    dat=dat.a*c
    # dat=dat.a
    
    TT=ncol(dat)
    N.ts = nrow(dat)
    # Standardize data
    # Sigma = sqrt(apply(dat, 1, var, na.rm=TRUE))
    datL = log(dat)
    y.bar = apply(datL, 1, mean, na.rm=TRUE)
    
    # y.bar = mean(datL)
    # dat.z = (dat - y.bar) * (1/Sigma)
    dat.dm = (datL - y.bar) / y.bar
    gsd = sd(c(dat.dm))
    dat.z = dat.dm/gsd
    # y.bar
    # gsd
    # apply(dat.z, 1, sd)
    
    datz_SD[i,] = apply(dat.z, 1, sd, na.rm=T)
    
    # dat.dm = (datL - y.bar) 
    dat.z = as.matrix(dat.z)
    rownames(dat.z) = c("Survey2","Survey1","Survey3")
    if(reorg==T){ rownames(dat.z) = c('Survey2','Survey3','Survey1')}
   
    
    
    ##### SN DELTA LOGNORMAL ####
    cntl.list = list(minit=200, maxit=500000, abstol=0.02, allow.degen=FALSE, conv.test.slope.tol = 0.05)
    if(reorg==T) { R = diag(c( CV2, CV3, CV1), nrow=3, ncol=3) }
    if(reorg==F) { R = diag(c( CV2, CV1, CV3), nrow=3, ncol=3) }
    # R = diag(c( CV2, CV1, CV3), nrow=3, ncol=3)
    # R = diag(c( CV1, CV2, CV3), nrow=3, ncol=3)
    
    # dfa = MARSS(dat.z, model = list(m=1,R="diagonal and equal"), control=cntl.list, form="dfa", z.score=FALSE)
    # dfa = MARSS(dat.z, model = list(m=1,R=R), control=cntl.list, form="dfa", z.score=TRUE)
    dfa = MARSS(dat.z, model = list(m=1,R=R), control=cntl.list, form="dfa", z.score=FALSE)
    
    # xx=exp(dfa$states*gsd)
    # xx2=exp(dfa$states*gsd + ( (dfa$states.se^2) /2) )
    # # dfa100 = dfa
    # # xx100 = xx
    # # xx5 = xx
    # # xx10 = xx
    # N = Ny[1, yrs]
    # # xxx=xx/c
    # plot(yrs,xx, type='l')
    # lines(yrs, rescale(N, to=c(min(xx),max(xx))), col='red')
    # plot(yrs,xx2, type='l')
    # lines(yrs, rescale(N, to=c(min(xx2),max(xx2))), col='red')
    # lines(yrs, rescale(xx, to=c(min(xx2),max(xx2))), col='red')
    # plot(yrs,xxx, type='l')
    # lines(yrs, rescale(N, to=c(min(xxx),max(xxx))), col='red')
    
    
    
    # dfa=dfa2
    pars = MARSSparamCIs(dfa, nboot=10000)
    # pars2 = MARSSparamCIs(dfa2, nboot=10000)
    # MARSSparamCIs(Mdglm.FINAL, nboot=10000, alpha=0.1)
    # MARSSparamCIs(Mdglm.FINAL, nboot=10000, alpha=0.15)
    # dfa = dfa2
    
    # get the inverse of the rotation matrix
    #H.cMdglm.inv = varimax(coef(Mdglm.FINAL, type="matrix")$Z)$rotmat
    Z.rot = coef(dfa, type="matrix")$Z #%*% H.cMdglm.inv
    trends.rot = dfa$states
    ts.trends = t(trends.rot)
    par.mat=coef(dfa, type="matrix")
    fit.b = par.mat$Z %*% dfa$states # + matrix(par.cMdglm.mat$D, nrow=N.ts) %*% SBnao
    
    assign(paste("Z.rot", i, sep="."), Z.rot)
    assign(paste("Z.upCI", i, sep="."), pars$par.upCI$Z)
    assign(paste("Z.lowCI", i, sep="."), pars$par.lowCI$Z)
    assign(paste("trends.rot", i, sep="."), trends.rot)
    
    
    # 
    # # plot factor loadings
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
    
    
    # PLOT WITH CI's
    yrs = 40:(years-35)
    
    index = dfa$states[1,]
    indexSE = dfa$states.se[1,]
    lowerCI = dfa$states[1,]-1.96*dfa$states.se[1,]
    upperCI = dfa$states[1,]+1.96*dfa$states.se[1,]
    assign(paste("index",i,sep=""),index)
    assign(paste("indexSE",i,sep=""),indexSE)
    assign(paste("lowerCI",i,sep=""),lowerCI)
    assign(paste("upperCI",i,sep=""),upperCI)
    
    
    # INDEX = (index - mean(index))/sd(index)
    # cis = c(lowerCI, upperCI)
    # CIS = (cis - mean(cis))/sd(cis)
    # LOWERCI
    
    par(mfrow=c(2,2))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='Year', ylab='Abundance', ylim=c(min(lowerCI), max(upperCI))))
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, index, type='l', lwd=2, pch=16))
    title(paste("Iteration",i,sep=" "))
    abline(h=0)
    
    N = Ny[i,40:65]
    NL=log(N)
    Nscale = (N-mean(N))/sd(N)
    NLscale = (N-mean(N))/sd(N)
    index.z = (index - mean(index))/sd(index)
    index.z = ifelse(is.na(index.z), 0, index.z)
    index.z = ifelse(abs(index.z)==Inf, 0, index.z)
    
    indexSEBT = indexSE*gsd
    indexBT = exp(index*gsd+ ((indexSEBT^2)/2))
    # indexBT = exp(index*gsd * (1/c))
    assign(paste("indexBT",i,sep=""),indexBT)
    assign(paste("indexSEBT",i,sep=""),indexSEBT)
    
    
    lines(yrs,rescale(NL,to=c(min(index),max(index))), type='l', col='red', lwd=2)
    # lines(yrs,rescale(N,to=c(min(index),max(index))), type='l', col='red', lwd=2)
    assign(paste("N",i,sep=""), N)
    assign(paste("Nscale",i,sep=""),Nscale)
    assign(paste("NL",i,sep=""), NL)
    assign(paste("NLscale",i,sep=""),NLscale)
    # write.csv(N,paste("N",i,".csv",sep=""))
    
    biasAll[i,] = index.z - NLscale
    bias[i] = mean(index.z - NLscale)
    RMSE[i] = sqrt(sum((index.z-NLscale)^2)/length(index.z))
    
    indexBT.z = (indexBT-mean(indexBT))/sd(indexBT)
    biasAll.BT[i,] = indexBT.z - Nscale
    RMSE.BT[i] = sqrt(sum((indexBT.z-Nscale)^2)/length(indexBT.z))
    
    
    # Plot fitted values
    survey = rownames(dat.z)
    # par(mfcol=c(ceiling(N.ts/2),2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
    for(n in 1:length(survey)){
      plot(yrs,dat.z[n,],xlab="",ylab="abundance index",bty="L", xaxt="n", ylim=c(-4,4), pch=16, col="blue")
      axis(1,labels=TRUE)
      lines(yrs,fit.b[n,], lwd=2)
      lines(yrs,rescale(NL, to=c(min(fit.b[n,]), max(fit.b[n,]))), type='l', col='red', lwd=2)
      title(paste("Atlantic sharpnose",survey[n], sep=" "))
    }
    
    
    #### ASSESS MODEL FITS
    sumResids = rowSums((dat.z-fit.b)^2, na.rm=TRUE)
    sumObserved = rowSums(dat.z^2, na.rm=TRUE)
    FitRatio = sumResids / sumObserved ; FitRatio
    # mean(FitRatio)
    # indexBT
    # dat.a
    
    FitRatios[i,] = FitRatio
    
    # mean(FitRatio)
    
    
    FactorLoadings = rbind(FactorLoadings, c(get(paste("Z.rot",i,sep=".")),get(paste("Z.lowCI",i,sep=".")),get(paste("Z.upCI",i,sep="."))))
    DFATrends = rbind(DFATrends, index)
    DFATrendsBT = rbind(DFATrendsBT, indexBT)
    DFATrendsSE = rbind(DFATrendsSE, indexSE)
    DFATrendsSEBT = rbind(DFATrendsSEBT, indexSEBT)
    upCI_DFATrends = rbind(upCI_DFATrends, upperCI)
    lowCI_DFATrends = rbind(lowCI_DFATrends, lowerCI)
    
    
    # par(mfrow=c(1,1))
    # plot(fit.b, dat.z-fit.b)
  } # END LOOP
  
  
  
  colnames(FactorLoadings) = c("FL1","FL2","FL3",
                               "lowCI_FL1","lowCI_FL2","lowCI_FL3",
                               "upCI_FL1","upCI_FL2","upCI_FL3")
  FR = apply(FitRatios,1,mean)
  FLoadings = cbind(FactorLoadings, meanFactorLoading = FR)
  
  DFA_Results = list()
  DFA_Results[["FitRatios"]] = FitRatios
  DFA_Results[["DFATrends"]] = DFATrends
  DFA_Results[["DFATrendsBT"]] = DFATrendsBT
  DFA_Results[["DFATrendsSE"]] = DFATrendsSE
  DFA_Results[["DFATrendsSEBT"]] = DFATrendsSEBT
  DFA_Results[["upCI_DFATrends"]] = upCI_DFATrends
  DFA_Results[["lowCI_DFATrends"]] = lowCI_DFATrends
  DFA_Results[["FLoadings"]] = FLoadings
  DFA_Results[["bias"]] = bias
  DFA_Results[["biasAll"]] = biasAll
  DFA_Results[["RMSE"]] = RMSE
  DFA_Results[["biasAll.BT"]] = biasAll.BT
  DFA_Results[["RMSE.BT"]] = RMSE.BT
  DFA_Results[["GSD"]] = gsd
  DFA_Results[["datz_SD"]] = datz_SD
  
  save(DFA_Results, file="DFA_Results.RData")
  
  
  # write.csv(FitRatios, "FitRatios.csv")
  # 
  # colnames(FactorLoadings) = c("FL1","FL2","FL3","lowCI_FL1","lowCI_FL2","lowCI_FL3","upCI_FL1","upCI_FL2","upCI_FL3")
  # write.csv(DFATrends, "DFATrends.csv")
  # write.csv(DFATrendsSE, "DFATrendsSE.csv")
  # write.csv(upCI_DFATrends, "upCI_DFATrends.csv")
  # write.csv(lowCI_DFATrends, "lowCI_DFATrends.csv")
  # 
  # apply(FitRatios,1,mean)
  # FLoadings = cbind(FactorLoadings, SUM = apply(FactorLoadings[,1:3], 1, sum))
  # write.csv(FLoadings, "FactorLoadings_sum.csv")
  # 
  # write.csv(bias, "bias.csv")
  # write.csv(biasAll, "biasAll.csv")
  # write.csv(RMSE, "RMSE.csv")
  # 
  
  ###################################  plotting ###################################
  png(filename="Trend.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    index = get(paste("index",l,sep=""))
    lowerCI = get(paste("lowerCI",l,sep=""))
    upperCI = get(paste("upperCI",l,sep=""))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, index, type='l', lwd=2, pch=16))
    # title(paste("Iteration",l,sep=" "))
    abline(h=0)
    NLscale = get(paste("NLscale",l,sep=""))
    lines(yrs,NLscale, type='l', col='red', lwd=2)
  }
  
  dev.off()
  
  
  png(filename="RawIndices.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  par(mfrow=c(10,10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    plot(1:years, Iy1[l,], type='l', col='black', xlim=c(min(yrs),max(yrs)), axes=F,
         ylim=c(0, max(c(Iy1[l,yrs], Iy2[l,yrs], Iy3[l,yrs]))+0.1))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    lines(1:years, Iy2[l,], type='l', col='blue')
    lines(1:years, Iy3[l,], type='l', col='red')
    # title(paste("Iteration",l, sep=" "))
  }
  
  dev.off()
  

  
  png(filename="FactorLoadings.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  # Before plotting factor loadings, we'll need to assign Z.rot, trends.rot
  par(mfrow=c(10,10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  survey = rownames(dat.z)
  minZ = 0.00
  for(l in 1:Niters){
    Z.rot = get(paste("Z.rot",l,sep="."))
    trends.rot = get(paste("trends.rot",l,sep="."))
    ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
    # par(mfrow=c(nrow(trends.rot),1), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
    for(h in 1:nrow(trends.rot)) {
      plot(c(1:N.ts)[abs(Z.rot[,h])>minZ], as.vector(Z.rot[abs(Z.rot[,h])>minZ,h]),
           type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1), axes=FALSE)
      axis(1, labels=FALSE)
      axis(2, labels=FALSE)
      #   for(j in 1:N.ts) {
      #     if(Z.rot[j,h] > minZ) {text(j, -0.05, survey[j], srt=90, adj=1, cex=0.9)}
      #     if(Z.rot[j,h] < -minZ) {text(j, 0.05, survey[j], srt=90, adj=0, cex=0.9)}
      #   } # end j loop
      #   mtext(paste("Factor loadings on trend",h,sep=" "),side=3,line=.5)
      abline(h=0, lwd=1, col="gray")
      abline(h=0.2, col='gray', lty=2)
      abline(h=-0.2, col='gray', lty=2)
    } # end h loop
  }
  
  dev.off()
  
  
  
  png(filename="Trend_scale.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    index = get(paste("index",l,sep=""))
    lowerCI = get(paste("lowerCI",l,sep=""))
    upperCI = get(paste("upperCI",l,sep=""))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, index, type='l', lwd=2, pch=16))
    # title(paste("Iteration",l,sep=" "))
    abline(h=0)
    NL = get(paste("NL",l,sep=""))
    lines(yrs, rescale(NL, to=c(min(index), max(index))), type='l',col="red", lwd=2)
  }
  
  
  dev.off()
  
  
  
  
  png(filename="BackTransformed_Trend_scale.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    indexBT = get(paste("indexBT",l,sep=""))
    indexSEBT = get(paste("indexSEBT",l,sep=""))
    lowerCI = indexBT-(1.96*indexSEBT)
    upperCI = indexBT+(1.96*indexSEBT)
    df <- data.frame(yrs, indexBT, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, indexBT, type='l', lwd=2, pch=16))
    # title(paste("Iteration",l,sep=" "))
    # abline(h=0)
    N = get(paste("N",l,sep=""))
    lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT))), type='l',col="red", lwd=2)
  }
  
  
  dev.off()
  
  
  png(filename="Trend_N.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    
    index = get(paste("index",l,sep=""))
    NL = get(paste("NL",l,sep=""))
    plot(yrs, NL, type='l',col="red", axes=FALSE)
    axis(1, labels=FALSE)
    axis(2, labels=FALSE)
    lines(yrs, rescale(index, to=c(min(NL),max(NL))))
    
  }
  
  dev.off()
  
  
  
  
  
  png(filename="Trend_Standardize.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    
    index = get(paste("index",l,sep=""))
    NL = get(paste("NL",l,sep=""))
    plot(yrs, (NL-mean(NL))/sd(NL), type='l',col="red", axes=FALSE, lwd=1.5)
    axis(1, labels=FALSE)
    axis(2, labels=FALSE)
    lines(yrs, (index-mean(index))/sd(index), type='l', lwd=1.5)
    
  }
  
  dev.off()
  
  
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    index = get(paste("index",l,sep=""))
    lowerCI = get(paste("lowerCI",l,sep=""))
    upperCI = get(paste("upperCI",l,sep=""))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, index, type='l', lwd=2, pch=16))
    # title(paste("Iteration",l,sep=" "))
    abline(h=0)
    NL = get(paste("NL",l,sep=""))
    lines(yrs, rescale(NL, to=c(min(index), max(index))), type='l',col="red", lwd=2)
  }
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    indexBT = get(paste("indexBT",l,sep=""))
    indexSEBT = get(paste("indexSEBT",l,sep=""))
    lowerCI = indexBT-(1.96*indexSEBT)
    upperCI = indexBT+(1.96*indexSEBT)
    df <- data.frame(yrs, indexBT, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, indexBT, type='l', lwd=2, pch=16))
    # title(paste("Iteration",l,sep=" "))
    # abline(h=0)
    N = get(paste("N",l,sep=""))
    lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT))), type='l',col="red", lwd=2)
  }
  
  # 
  # png(filename="TwoAxesCompare.png", 
  #     type="cairo",
  #     units="mm", 
  #     width=300, 
  #     height=300, 
  #     pointsize=12, 
  #     res=600)
  # 
  # par(mfrow=c(10, 10))
  # par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  # for(l in 1:Niters){
  #   index = get(paste("index",l,sep=""))
  #   lowerCI = get(paste("lowerCI",l,sep=""))
  #   upperCI = get(paste("upperCI",l,sep=""))
  #   df <- data.frame(yrs, index, lowerCI, upperCI)
  #   with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
  #   axis(1,labels=FALSE)
  #   axis(2,labels=FALSE)
  #   x <- as.numeric(df$yrs)
  #   polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
  #   with(df, lines(x, index, type='l', lwd=2, pch=16))
  #   # title(paste("Iteration",l,sep=" "))
  #   abline(h=0)
  #   par(new=T)
  #   NL = get(paste("NL",l,sep=""))
  #   plot(yrs,NL, type='l', col='red', axes=FALSE, lwd=1.5)
  #   axis(4, labels=FALSE)
  # }
  # 
  # dev.off()
  # 
  
  
  
  
  
  
  
  png(filename="TwoAxesCompare_TryCorrected.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.3),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    
    index = get(paste("index",l,sep=""))
    NL = get(paste("NL",l,sep=""))
    plot(yrs, NL, type='l',col="red", axes=FALSE, lwd=1.5)
    axis(1, labels=FALSE)
    axis(2, labels=FALSE)
    par(new=T)
    if(FLoadings[l,10] >= 0){
      plot(yrs,index, type='l', col='black', axes=FALSE, lwd=1.5)
    }
    if(FLoadings[l,10] < 0){
      plot(yrs,-1*index, type='l', col='black', axes=FALSE, lwd=1.5)
    }
    axis(4, labels=FALSE)
    
  }
  
  dev.off()
  
  
  
  
  ##################################################### END #################################################
  
  return(list(RMSE=RMSE, biasAll=biasAll, FitRatios=FitRatios, RMSE.BT=RMSE.BT, biasAll.BT=biasAll.BT))} # END FUNCTION


########################



########################## FOUR INDICES SIMULATION FUNCTION ##############################################################################

# when change in q
# IF F1 --> c=c(1000, 2, 100, 3)
# IF F2 --> c=c(15, 40, 2000, 4) 
# IF F3 --> c=c(8000, 100, 3, 1.7) 
# IF F4 --> c=c(600, 10, 100, 2)  

# when constant q
# IF F1 -->   c=c(10, 10, 8, 250)
# IF F2 --> c=c(40, 60, 40, 5)
# IF F3 --> c=c(25, 30, 25, 1)
# IF F4 --> c=c(25, 30, 25, 5)

# when gradual change in q
# IF F1 --> c(5000, 2.5, 3, 6)
# IF F2 --> c(3.8, 130, 200, 9)
# IF F3 --> c=c(2000, 80, 9, 1.5)
# IF F4 --> c=c(500, 15, 5, 2.5)


#------------------------------------------------------------
# Data for run
F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50)))
CV1=0.7
CV2=0.5
CV3=0.3
CV4=0.5
q_survey1 = rep(0.001, length=years)
q_survey2 = rep(0.001, length=years)
q_survey3 = rep(0.001, length=years)
q_survey4 = rep(0.003, length=years)
Iay_eq = Iay_eq_F2
# 
# CV4 = 0.5
# q_survey4 = rep(0.003, length=years)
F_var = as.vector(rep(0.2,years))
CV1=0.5
CV2=0.5
CV3=0.5
CV4=0.5
q_survey1 = c(rep(0.00025, 49), seq(0.0004, 0.0014, by=0.0001), rep(0.0015, 40))
q_survey2 = rep(0.001, length=years)
q_survey3 = c(rep(0.002, 49), seq(0.0018, 0.0006, by=-0.0002), rep(0.0005, 44))
q_survey4 = rep(0.003, length=years)
Iay_eq = Iay_eq_F1

F_var = as.vector(rep(0.2,years))
CV1=0.5
CV2=0.5
CV3=0.5
CV4=0.5
q_survey1 = rep(0.001, length=years)
q_survey2 = rep(0.001, length=years)
q_survey3 = rep(0.001, length=years)
q_survey4 = rep(0.003, length=years)
Iay_eq = Iay_eq_F1


F_var = as.vector(c(rep(0.4, 50), rep(0, 50)))
CV1=0.5
CV2=0.5
CV3=0.5
CV4=0.5
q_survey1 = rep(0.001, length=years)
q_survey2 = rep(0.001, length=years)
q_survey3 = rep(0.001, length=years)
q_survey4 = rep(0.003, length=years)
Iay_eq = Iay_eq_F3


F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45)))
CV1=0.5
CV2=0.5
CV3=0.5
CV4=0.5
q_survey1 = rep(0.001, length=years)
q_survey2 = rep(0.001, length=years)
q_survey3 = rep(0.001, length=years)
q_survey4 = rep(0.003, length=years)
Iay_eq = Iay_eq_F4

q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35))
q_survey2 = rep(0.001, length=years)
q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40))
q_survey4 = rep(0.003, length=years)

F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50)))
CV1=0.5
CV2=0.5
CV3=0.5
CV4=0.5
q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2))
q_survey2 = rep(0.001, length=years)
q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2))
q_survey4 = rep(0.003, length=years)
Iay_eq = Iay_eq_F2
c=c(30, 40, 40, 4) 

F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50)))
CV1=0.3
CV2=0.5
CV3=0.7
CV4=0.5
q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2))
q_survey2 = rep(0.001, length=years)
q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2))
q_survey4 = rep(0.003, length=years)
Iay_eq = Iay_eq_F2
c=c(15, 40, 2000, 4) 


setwd("C:\\~")
###############


DFA_Sim4 = function(F_var, CV1, CV2, CV3, CV4, q_survey1, q_survey2, q_survey3, q_survey4, Iay_eq=Iay_eq, c=5, reorg=F){
  # Iay
  par(mfrow=c(1,1))
  plot(1:years, apply(Iay_eq, 2, sum), type='l',ylim=c(0,max(apply(Iay_eq, 2, sum),na.rm=T)+max(apply(Iay_eq, 2, sum), na.rm=T)*0.1))
  abline(v=40, lwd=2, col="red")
  abline(v=65, lwd=2, col="red")
  # plot(1:years, apply(Iay, 2, sum), type='l')
  # Iy = apply(Iay, 2, sum)
  
  F_var=F_var
  CV1=CV1
  CV2=CV2
  CV3=CV3
  q_survey1=q_survey1
  q_survey2=q_survey2
  q_survey3=q_survey3
  c=c
  c1=c
  
  ####### ADD STOCHASTICITY ###########
  Niters = 100
  Iay1 = matrix(nrow=A, ncol=years)
  Iay2 = matrix(nrow=A, ncol=years)
  Iay3 = matrix(nrow=A, ncol=years)
  Iay4 = matrix(nrow=A, ncol=years)
  Iay1_CV = matrix(nrow=A, ncol=years)
  Iay2_CV = matrix(nrow=A, ncol=years)
  Iay3_CV = matrix(nrow=A, ncol=years)
  Iay4_CV = matrix(nrow=A, ncol=years)
  Cy = matrix(nrow=Niters, ncol=years)
  Ny = matrix(nrow=Niters, ncol=years)
  Iy1 = matrix(nrow=Niters, ncol=years)
  Iy2 = matrix(nrow=Niters, ncol=years)
  Iy3 = matrix(nrow=Niters, ncol=years)
  Iy4 = matrix(nrow=Niters, ncol=years)
  Iy1_CV = matrix(nrow=Niters, ncol=years)
  Iy2_CV = matrix(nrow=Niters, ncol=years)
  Iy3_CV = matrix(nrow=Niters, ncol=years)
  Iy4_CV = matrix(nrow=Niters, ncol=years)
  M0 = matrix(nrow=Niters, ncol=years)
  Npups = matrix(nrow=Niters, ncol=years)
  FSS = matrix(nrow=Niters, ncol=years) # Female spawning stock 
  Z0 = 1.904416
  Zmin = 0.209
  Npups_eq = 3320.671
  beta=0.5807613
  ResultsList = list()
  
  # N_All = vector()
  # N_STD = vector()
  
  set.seed(430)
  
  # F_const = rep(0, yrs)
  
  for(k in 1:Niters){
    # paste("Iteration",k,sep=" ")
    # Create matrices to store results
    Nay = matrix(nrow=A, ncol=years)
    Cay = matrix(nrow=A, ncol=years)
    # M0 = vector(length=years)
    
    # Inputs: 
    M = M_const
    # FM = F_const
    FM = F_var
    
    # Set up equilibrium conditions
    Nay[,1] = Equilibrium_Nay
    
    for (i in 1:(years-1)){
      
      matjit = rnorm(length(mat),mat,mat*0.01)
      matjit = ifelse(matjit>1, 1, matjit)
      FSS[k,i] = sum(Nay[,i]*matjit) #Fecund Stock size
      Npups[k,i] = sum(Nay[,i]*matjit*rnorm(length(fec), fec, fec*0.1)) # equilibrium Npups = 3320.671
      M0[k,i] = rnorm(1, ( ( (1- (Npups[k,i]/Npups_eq)^beta ) * (Zmin - Z0) ) + Z0 ) , 0.1)
      # M0[i] = 1.9044160 # THIS IS EQUILIBRIUM!!!
      Nay[1,i+1] = Npups[k,i]*exp(-M0[k,i])  # pups produced that survive to age 1 = Recruits!
      
      Mjit = rnorm(length(M), M, M*0.1) 
      
      Zay = vector(length=length(M))
      for (j in 1:A){
        Zay[j] = Mjit[j] + (Sel_catch[j,2]*FM[i]/4) + (Sel_catch[j,3]*FM[i]/4) + (Sel_catch[j,4]*FM[i]/4) + (Sel_catch[j,5]*FM[i]/4)
      } #end j loop Z
      
      
      for (j in 1:(A-2)){
        Nay[j+1,i+1] = Nay[j,i]*exp(-Zay[j])
        # Say[j,i]=Nay[j,i]*mat[j]
        # Pay[j,i]=Say[j,i]*fec[j]
        Cay[j,i]=Nay[j,i]*((Sel_catch[j,2]*FM[i]/4)/Zay[j])*(1-exp(-Zay[j])) +
          Nay[j,i]*((Sel_catch[j,3]*FM[i]/4)/Zay[j])*(1-exp(-Zay[j])) +
          Nay[j,i]*((Sel_catch[j,4]*FM[i]/4)/Zay[j])*(1-exp(-Zay[j])) +
          Nay[j,i]*((Sel_catch[j,5]*FM[i]/4)/Zay[j])*(1-exp(-Zay[j]))
        
      } # END J LOOP
      Nay[A,i+1] = Nay[A-1,i]*exp(-Zay[A-1]) +
        Nay[A,i]*exp(-Zay[A])
      
      Cay[A-1,i]=Nay[A-1,i]*((Sel_catch[A-1,2]*FM[i]/4)/Zay[A-1])*(1-exp(-Zay[A-1])) +
        Nay[A-1,i]*((Sel_catch[A-1,3]*FM[i]/4)/Zay[A-1])*(1-exp(-Zay[A-1])) +
        Nay[A-1,i]*((Sel_catch[A-1,4]*FM[i]/4)/Zay[A-1])*(1-exp(-Zay[A-1])) +
        Nay[A-1,i]*((Sel_catch[A-1,5]*FM[i]/4)/Zay[A-1])*(1-exp(-Zay[A-1]))
      
      Cay[A,i]=Nay[A,i]*((Sel_catch[A,2]*FM[i]/4)/Zay[A])*(1-exp(-Zay[A])) +
        Nay[A,i]*((Sel_catch[A,3]*FM[i]/4)/Zay[A])*(1-exp(-Zay[A])) +
        Nay[A,i]*((Sel_catch[A,4]*FM[i]/4)/Zay[A])*(1-exp(-Zay[A])) +
        Nay[A,i]*((Sel_catch[A,5]*FM[i]/4)/Zay[A])*(1-exp(-Zay[A]))
      
      
      
      
      for(h in 1:A){
        Iay1_CV[h,i] = runif(1, min=CV1-0.1, max=CV1+0.1)
        Iay1[h,i] = max(q_survey1[i]*Sel_survey[h,2] * Nay[h,i] * 
                          exp( rnorm(1, 0, (q_survey1[i]*Sel_survey[h,2]*Nay[h,i])*Iay1_CV[h,i]) - 
                                 ( ((q_survey1[i]*Sel_survey[h,2]*Nay[h,i]*Iay1_CV[h,i])^2)/2 ) ), 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay2_CV[h,i] = runif(1, min=CV2-0.1, max=CV2+0.1)
        # Iay2[h,i] = max(q_survey2[i]*Sel_survey[h,3] * Nay[h,i] +  rnorm(1, 0, (q_survey2[i]*Sel_survey[h,3]*Nay[h,i])*Iay2_CV[h,i] ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay2[h,i] = max(q_survey2[i]*Sel_survey[h,3] * Nay[h,i] * 
                          exp( rnorm(1, 0, (q_survey2[i]*Sel_survey[h,3]*Nay[h,i])*Iay2_CV[h,i]) - 
                                 ( ((q_survey2[i]*Sel_survey[h,3]*Nay[h,i]*Iay2_CV[h,i])^2)/2 ) ), 0) #LOTS OF UNCERTAINTY!based off CV;
        
        
        Iay3_CV[h,i] = runif(1, min=CV3-0.1, max=CV3+0.1)
        # Iay3[h,i] = max(q_survey3[i]*Sel_survey[h,2] * Nay[h,i] +  rnorm(1, 0, (q_survey3[i]*Sel_survey[h,2]*Nay[h,i])*Iay3_CV[h,i] ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay3[h,i] = max(q_survey3[i]*Sel_survey[h,2] * Nay[h,i] * 
                          exp( rnorm(1, 0, (q_survey3[i]*Sel_survey[h,2]*Nay[h,i])*Iay3_CV[h,i]) - 
                                 ( ((q_survey3[i]*Sel_survey[h,2]*Nay[h,i]*Iay3_CV[h,i])^2)/2 ) ), 0) #LOTS OF UNCERTAINTY!based off CV;
        
        Iay4_CV[h,i] = runif(1, min=CV4-0.1, max=CV4+0.1)
        # Iay4[h,i] = max(q_survey4[i]*Sel_survey[h,2] * Nay[h,i] +  rnorm(1, 0, (q_survey4[i]*Sel_survey[h,2]*Nay[h,i])*Iay4_CV[h,i] ), 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay4[h,i] = max(q_survey4[i]*Sel_survey[h,2] * Nay[h,i] * 
                          exp( rnorm(1, 0, (q_survey4[i]*Sel_survey[h,2]*Nay[h,i])*Iay4_CV[h,i]) - 
                                 ( ((q_survey4[i]*Sel_survey[h,2]*Nay[h,i]*Iay4_CV[h,i])^2)/2 ) ), 0) #LOTS OF UNCERTAINTY!based off CV;
        
      } # end h
      
      
      # LCay1 = von bert curve to calculate length composition
      # random q: runif(1,max=0.005, min=0.0005)
    }# ends i loop
    
    # USE VON BERT CURVE TO ESTIMATE LENGTH COMP. Let von Bert parameters vary by Niter...???? 
    #   year (implement in i loop)... 
    #   Or should they vary by year class? Could create a matrix of von bert parameters based on age class 
    #       and pull out appropriate values as needed
    # Lt = Linf * (1 - exp(-K*(t - t0)))
    Linf=81.6
    K=0.49
    t0=-0.97
    
    # L0 = Linf * (1 - exp(-K*(0-t0)))
    
    
    # GET LENGTH FREQUENCIES
    Cay_N = cbind(ages, round(Cay))
    assign(paste("LF_Cay",k,sep="_"), data.frame(LBins = seq(20, 122, by=2)))
    for(f in 40:65){
      dfC = as.data.frame(Cay_N[,c(1,f+1)])
      colnames(dfC) = c("ages","count") 
      df.expandedC = dfC[rep(row.names(dfC), dfC$count),1:2]
      vonBertC = rbvn(nrow(df.expandedC),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expandedC$L = vonBertC$Linf * ( 1 - exp(-(vonBertC$K) * (jitter(df.expandedC$ages,amount=0.5) - t0) ))
      hstC=hist(df.expandedC$L, breaks=c(seq(20, 124, by=2)), plot=F)
      assign(paste("LF_Cay",k,sep="_"), cbind(get(paste("LF_Cay",k,sep="_")), hstC$counts))
    } # end f loop
    # write.csv(get(paste("LF_Cay",k,sep="_")), paste(paste("LF_Cay",k,sep="_"),".csv",sep=""))
    
    Iay1_N = cbind(ages, round(Iay1*100))
    assign(paste("LF_Iay1",k,sep="_"), data.frame(LBins = seq(20, 122, by=2)))
    for(l in 40:65){
      df1 = as.data.frame(Iay1_N[,c(1,l+1)])
      colnames(df1) = c("ages","count") 
      df.expanded1 = df1[rep(row.names(df1), df1$count),1:2]
      vonBert1 = rbvn(nrow(df.expanded1),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expanded1$L = vonBert1$Linf * ( 1 - exp(-(vonBert1$K) * (jitter(df.expanded1$ages,amount=0.5) - t0) ))
      hst1=hist(df.expanded1$L, breaks=c(seq(20, 124, by=2)), plot=F)
      assign(paste("LF_Iay1",k,sep="_"), cbind(get(paste("LF_Iay1",k,sep="_")), hst1$counts))
    } # end l loop
    # write.csv(get(paste("LF_Iay1",k,sep="_")), paste(paste("LF_Iay1",k,sep="_"),".csv",sep=""))
    
    Iay2_N = cbind(ages, round(Iay2*100))
    assign(paste("LF_Iay2",k,sep="_"), data.frame(LBins = seq(20, 122, by=2)))
    for(m in 40:65){
      df2 = as.data.frame(Iay2_N[,c(1,m+1)])
      colnames(df2) = c("ages","count") 
      df.expanded2 = df2[rep(row.names(df2), df2$count),1:2]
      vonBert2 = rbvn(nrow(df.expanded2),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expanded2$L = vonBert2$Linf * ( 1 - exp(-(vonBert2$K) * (jitter(df.expanded2$ages,amount=0.5) - t0) ))
      hst2=hist(df.expanded2$L, breaks=c(seq(20, 124, by=2)), plot=F)
      assign(paste("LF_Iay2",k,sep="_"), cbind(get(paste("LF_Iay2",k,sep="_")), hst2$counts))
    } # end m loop
    # write.csv(get(paste("LF_Iay2",k,sep="_")), paste(paste("LF_Iay2",k,sep="_"),".csv",sep=""))
    
    
    Iay3_N = cbind(ages, round(Iay3*100))
    assign(paste("LF_Iay3",k,sep="_"), data.frame(LBins = seq(20, 122, by=2)))
    for(n in 40:65){
      df3 = as.data.frame(Iay3_N[,c(1,n+1)])
      colnames(df3) = c("ages","count") 
      df.expanded3 = df3[rep(row.names(df3), df3$count),1:2]
      vonBert3 = rbvn(nrow(df.expanded3),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expanded3$L = vonBert3$Linf * ( 1 - exp(-(vonBert3$K) * (jitter(df.expanded3$ages,amount=0.5) - t0) ))
      hst3=hist(df.expanded3$L, breaks=c(seq(20, 124, by=2)), plot=F)
      assign(paste("LF_Iay3",k,sep="_"), cbind(get(paste("LF_Iay3",k,sep="_")), hst3$counts))
    } # end n loop
    # write.csv(get(paste("LF_Iay3",k,sep="_")), paste(paste("LF_Iay3",k,sep="_"),".csv",sep=""))
    
    Iay4_N = cbind(ages, round(Iay4*100))
    assign(paste("LF_Iay4",k,sep="_"), data.frame(LBins = seq(20, 122, by=2)))
    for(u in 40:65){
      df4 = as.data.frame(Iay4_N[,c(1,u+1)])
      colnames(df4) = c("ages","count") 
      df.expanded4 = df4[rep(row.names(df4), df4$count),1:2]
      vonBert4 = rbvn(nrow(df.expanded4),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expanded4$L = vonBert4$Linf * ( 1 - exp(-(vonBert4$K) * (jitter(df.expanded4$ages,amount=0.5) - t0) ))
      hst4=hist(df.expanded4$L, breaks=c(seq(20, 124, by=2)), plot=F)
      assign(paste("LF_Iay4",k,sep="_"), cbind(get(paste("LF_Iay4",k,sep="_")), hst4$counts))
    } # end u loop
    # write.csv(get(paste("LF_Iay4",k,sep="_")), paste(paste("LF_Iay4",k,sep="_"),".csv",sep=""))
    
    
    
    #save length-frequencies
    LFreqList = list()
    LFreqList[[paste("LF_Cay",k,sep="_")]] = get(paste("LF_Cay",k,sep="_"))
    LFreqList[[paste("LF_Iay1",k,sep="_")]] = get(paste("LF_Iay1",k,sep="_"))
    LFreqList[[paste("LF_Iay2",k,sep="_")]] = get(paste("LF_Iay2",k,sep="_"))
    LFreqList[[paste("LF_Iay3",k,sep="_")]] = get(paste("LF_Iay3",k,sep="_"))
    LFreqList[[paste("LF_Iay4",k,sep="_")]] = get(paste("LF_Iay4",k,sep="_"))
    
    save(LFreqList, file=paste("LFreqList_",k,".RData",sep=""))
    
    
    # Linf = rnorm(years, 81.9, 81.9*0.1) # CV = 0.2; 0.2 = sd / 81.9
    # K = rnorm(years, 0.48, 0.48*0.1)
    # t0 = rnorm(years, -0.99, 0.99*0.1)
    # LC = Linf * ( 1 - exp(-K*(c(0,ages) - t0)))
    # assign(paste("LC",k,sep="_"), cbind(c(0,ages), LC))
    # write.csv(LC, paste(paste("LC",k,sep="_"),".csv",sep=""))
    # assign(paste("M0",k,sep="_"), M0)
    
    
    
    # write.csv(Iay1, paste(paste("Iay1",k,sep="_"),".csv",sep=""))
    # write.csv(Iay2, paste(paste("Iay2",k,sep="_"),".csv",sep=""))
    # write.csv(Iay3, paste(paste("Iay3",k,sep="_"),".csv",sep=""))
    # write.csv(Iay4, paste(paste("Iay4",k,sep="_"),".csv",sep=""))
    
    Iy1[k,] = apply(Iay1, 2, sum)
    Iy2[k,] = apply(Iay2, 2, sum)
    Iy3[k,] = apply(Iay3, 2, sum)
    Iy4[k,] = apply(Iay4, 2, sum)
    lines(1:years, Iy1[k,], type='l', col='grey')
    lines(1:years, Iy2[k,], type='l', col='blue')
    lines(1:years, Iy3[k,], type='l', col='red')
    lines(1:years, Iy4[k,], type='l', col='green')
    Iy1_CV[k,] = apply(Iay1_CV, 2, mean, na.rm=T)
    Iy2_CV[k,] = apply(Iay2_CV, 2, mean, na.rm=T)
    Iy3_CV[k,] = apply(Iay3_CV, 2, mean, na.rm=T)
    Iy4_CV[k,] = apply(Iay4_CV, 2, mean, na.rm=T)
    
    Cy[k,] = apply(Cay, 2, sum)
    Ny[k,] = apply(Nay, 2, sum)
    
    
    # plot(1:years, colSums(Nay),type='l', ylim=c(0,max(colSums(Nay))+20))
    # lines(1:years, colSums(Nay),type='l',col='black')
    # paste("Nay_",k, sep='')
    
    #write 
    # write.csv(Nay, paste(paste("Nay", k, sep="_"),".csv",sep=""))
    assign(paste("Nay", k, sep="_"), Nay)
    # write.csv(Cay, paste(paste("Cay", k, sep="_"),".csv",sep=""))
    assign(paste("Cay", k, sep="_"), Cay)
    
    
    #write 
    AnnualNC = list()
    AnnualNC[[paste("Nay",k,sep="_")]] = Nay
    AnnualNC[[paste("Cay",k,sep="_")]] = Cay
    AnnualNC[[paste("Iay1",k,sep="_")]] = Iay1
    AnnualNC[[paste("Iay2",k,sep="_")]] = Iay2
    AnnualNC[[paste("Iay3",k,sep="_")]] = Iay3
    AnnualNC[[paste("Iay4",k,sep="_")]] = Iay4
    
    save(AnnualNC, file=paste("AnnualNC_",k,".RData",sep=""))
    
    
    
  } # end k loop
  
  
  
  ResultsList[["Iy1"]] = Iy1
  ResultsList[["Iy2"]] = Iy2
  ResultsList[["Iy3"]] = Iy3
  ResultsList[["Iy4"]] = Iy4
  
  ResultsList[["Iy1_CV"]] = Iy1_CV
  ResultsList[["Iy2_CV"]] = Iy2_CV
  ResultsList[["Iy3_CV"]] = Iy3_CV
  ResultsList[["Iy4_CV"]] = Iy4_CV
  
  
  ResultsList[["M0"]] = M0
  ResultsList[["FSS"]] = FSS
  ResultsList[["Npups"]] = Npups
  
  ResultsList[["Cy"]] = Cy
  ResultsList[["Ny"]] = Ny
  
  save(ResultsList, file="ResultsList.RData")
  
  # write.csv(Iy1, "Iy1.csv")
  # write.csv(Iy2, "Iy2.csv")
  # write.csv(Iy3, "Iy3.csv")
  # write.csv(Iy4, "Iy4.csv")
  # write.csv(M0, "M0.csv") # populates where each row is an iteration and each column is one year
  # write.csv(FSS, "FSS.csv") # populates where each row is an iteration and each column is one year
  # write.csv(Iy1_CV, "Iy1_CV.csv")
  # write.csv(Iy2_CV, "Iy2_CV.csv")
  # write.csv(Iy3_CV, "Iy3_CV.csv")
  # write.csv(Iy4_CV, "Iy4_CV.csv")
  # write.csv(Cy, "Cy.csv")
  
  
  ########## TEST DFA #############
  
  # library(MARSS)
  yrs = 40:(years-35)
  biasAll = matrix(nrow=Niters, ncol=length(yrs))
  bias = vector(length=Niters)
  biasAll.BT = matrix(nrow=Niters, ncol=length(yrs))
  RMSE.BT = vector(length=Niters)
  
  RMSE = vector(length=Niters)
  FitRatios = matrix(nrow=Niters, ncol=4)
  FactorLoadings=vector()
  DFATrends = vector()
  DFATrendsBT = vector()
  DFATrendsSE = vector()
  DFATrendsSEBT = vector()
  upCI_DFATrends = vector()
  lowCI_DFATrends = vector()
  
  datz_SD = matrix(nrow=Niters, ncol=4)
  # i=1
  
  for(i in 1:Niters){
    
    if(reorg == T){
      assign(paste("dat",i,sep=""), rbind(Iy2[i,40:(years-35)], Iy1[i,40:(years-35)], Iy3[i,40:(years-35)], c(rep(NA,15), Iy4[i,55:(years-35)])))
      c=c( c1[2], c1[1], c1[3], c1[4])
      
    }
    
    if(reorg==F){
      assign(paste("dat",i,sep=""), rbind(Iy1[i,40:(years-35)], Iy2[i,40:(years-35)], Iy3[i,40:(years-35)], c(rep(NA,15), Iy4[i,55:(years-35)])))
    }
    
    
    dat.a = get(paste("dat",i,sep=""))
    # c=c(40, 45, 60, 4) 
    dat=dat.a*c
    
    
    TT=ncol(dat)
    N.ts = nrow(dat)
    
    datL = log(dat)
    y.bar = apply(datL, 1, mean, na.rm=TRUE)

    dat.dm = (datL - y.bar) / y.bar
    gsd = sd(c(dat.dm),na.rm=TRUE)
    dat.z = dat.dm/gsd
    # gsd
    # y.bar
    # apply(dat.z, 1, mean, na.rm=T)
    # apply(dat.z, 1, sd, na.rm=T)
    datz_SD[i,] = apply(dat.z, 1, sd, na.rm=T)
   
    # dat.dm = (datL - y.bar) 
    dat.z = as.matrix(dat.z)
    rownames(dat.z) = c("Survey1","Survey2","Survey3","Survey4")
    if(reorg==T) { rownames(dat.z) = c("Survey2","Survey1","Survey3","Survey4") } 
    
    
    
    ##### SN DELTA LOGNORMAL ####
    cntl.list = list(minit=200, maxit=500000, abstol=0.02, allow.degen=FALSE, conv.test.slope.tol = 0.05)
    R = diag(c(CV1, CV2, CV3, CV4), nrow=4, ncol=4)
    if(reorg==T) { R = diag(c(CV2, CV1, CV3, CV4), nrow=4, ncol=4) }
    
    # dfa = MARSS(dat.z, model = list(m=1,R="diagonal and equal"), control=cntl.list, form="dfa", z.score=TRUE)
    # dfa = MARSS(dat.z, model = list(m=1,R="diagonal and equal"), control=cntl.list, form="dfa", z.score=FALSE)
    dfa = MARSS(dat.z, model = list(m=1,R=R), control=cntl.list, form="dfa", z.score=FALSE)
    
    pars = MARSSparamCIs(dfa, nboot=10000)
    # MARSSparamCIs(Mdglm.FINAL, nboot=10000, alpha=0.1)
    # MARSSparamCIs(Mdglm.FINAL, nboot=10000, alpha=0.15)
    
    # get the inverse of the rotation matrix
    #H.cMdglm.inv = varimax(coef(Mdglm.FINAL, type="matrix")$Z)$rotmat
    Z.rot = coef(dfa, type="matrix")$Z #%*% H.cMdglm.inv
    trends.rot = dfa$states
    ts.trends = t(trends.rot)
    par.mat=coef(dfa, type="matrix")
    fit.b = par.mat$Z %*% dfa$states # + matrix(par.cMdglm.mat$D, nrow=N.ts) %*% SBnao
    
    assign(paste("Z.rot", i, sep="."), Z.rot)
    assign(paste("Z.upCI", i, sep="."), pars$par.upCI$Z)
    assign(paste("Z.lowCI", i, sep="."), pars$par.lowCI$Z)
    assign(paste("trends.rot", i, sep="."), trends.rot)
    
    
    # 
    # # plot factor loadings
    par(mfrow=c(3,2))
    survey = rownames(dat.z)
    minZ = 0.00
    ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
    # par(mfrow=c(nrow(trends.rot),1), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
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
    
    
    # PLOT WITH CI's
    yrs = 40:(years-35)
    
    index = dfa$states[1,]
    indexSE = dfa$states.se[1,]
    lowerCI = dfa$states[1,]-1.96*dfa$states.se[1,]
    upperCI = dfa$states[1,]+1.96*dfa$states.se[1,]
    assign(paste("index",i,sep=""),index)
    assign(paste("indexSE",i,sep=""),indexSE)
    assign(paste("lowerCI",i,sep=""),lowerCI)
    assign(paste("upperCI",i,sep=""),upperCI)
    
    # par(mfrow=c(2,2))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='Year', ylab='Abundance', ylim=c(min(lowerCI), max(upperCI))))
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, index, type='l', lwd=2, pch=16))
    title(paste("Iteration",i,sep=" "))
    abline(h=0)
    
    N = Ny[i,40:65]
    NL=log(N)
    Nscale = (N-mean(N))/sd(N)
    NLscale = (N-mean(N))/sd(N)
    index.z = (index - mean(index))/sd(index)
    index.z = ifelse(is.na(index.z), 0, index.z)
    index.z = ifelse(abs(index.z)==Inf, 0, index.z)
    
    
    indexSEBT = indexSE*gsd
    indexBT = exp(index*gsd+ ((indexSEBT^2)/2))
    assign(paste("indexBT",i,sep=""),indexBT)
    assign(paste("indexSEBT",i,sep=""),indexSEBT)
    
    lines(yrs,rescale(NL,to=c(min(index),max(index))), type='l', col='red', lwd=2)

    assign(paste("N",i,sep=""), N)
    assign(paste("Nscale",i,sep=""),Nscale)
    assign(paste("NL",i,sep=""), NL)
    assign(paste("NLscale",i,sep=""),NLscale)
    
    # N_All = rbind(N_All, N)
    # N_STD = rbind(N_STD, Nscale)
    
    biasAll[i,] = index.z - NLscale
    bias[i] = mean(index.z - NLscale)
    RMSE[i] = sqrt(sum((index.z-NLscale)^2)/length(index.z))
    
    indexBT.z = (indexBT-mean(indexBT))/sd(indexBT)
    biasAll.BT[i,] = indexBT.z - Nscale
    RMSE.BT[i] = sqrt(sum((indexBT.z-Nscale)^2)/length(indexBT.z))
    
    
    # RMSE_yr[i] = sqrt(sum((Nscale - indexscale)^2)/length(indexscale))
    
    # Plot fitted values
    survey = rownames(dat.z)
    # par(mfcol=c(ceiling(N.ts/2),2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
    for(n in 1:length(survey)){
      plot(yrs,dat.z[n,],xlab="",ylab="abundance index",bty="L", xaxt="n", ylim=c(-4,4), pch=16, col="blue")
      axis(1,labels=TRUE)
      lines(yrs,fit.b[n,], lwd=2)
      lines(yrs,rescale(NL, to=c(min(fit.b[n,]), max(fit.b[n,]))), type='l', col='red', lwd=2)
      title(paste("Atlantic sharpnose",survey[n], sep=" "))
    }
    
    
    #### ASSESS MODEL FITS
    sumResids = rowSums((dat.z-fit.b)^2, na.rm=TRUE)
    sumObserved = rowSums(dat.z^2, na.rm=TRUE)
    FitRatio = sumResids / sumObserved ; FitRatio
    FitRatios[i,] = FitRatio
    # mean(FitRatio)
    # indexBT
    # dat.a
    
    FactorLoadings = rbind(FactorLoadings, c(get(paste("Z.rot",i,sep=".")),get(paste("Z.lowCI",i,sep=".")),get(paste("Z.upCI",i,sep="."))))
    DFATrends = rbind(DFATrends, index)
    DFATrendsBT = rbind(DFATrendsBT, indexBT)
    DFATrendsSE = rbind(DFATrendsSE, indexSE)
    DFATrendsSEBT = rbind(DFATrendsSEBT, indexSEBT)
    upCI_DFATrends = rbind(upCI_DFATrends, upperCI)
    lowCI_DFATrends = rbind(lowCI_DFATrends, lowerCI)
    
    
    # par(mfrow=c(1,1))
    # plot(fit.b, dat.z-fit.b)
  } # END LOOP
  
  
  
  
  colnames(FactorLoadings) = c("FL1","FL2","FL3","FL4",
                               "lowCI_FL1","lowCI_FL2","lowCI_FL3","lowCI_FL4",
                               "upCI_FL1","upCI_FL2","upCI_FL3","upCI_FL4")
  FR = apply(FitRatios,1,mean)
  FLoadings = cbind(FactorLoadings, meanFactorLoading = FR)
  
  DFA_Results = list()
  DFA_Results[["FitRatios"]] = FitRatios
  DFA_Results[["DFATrends"]] = DFATrends
  DFA_Results[["DFATrendsBT"]] = DFATrendsBT
  DFA_Results[["DFATrendsSE"]] = DFATrendsSE
  DFA_Results[["DFATrendsSEBT"]] = DFATrendsSEBT
  DFA_Results[["upCI_DFATrends"]] = upCI_DFATrends
  DFA_Results[["lowCI_DFATrends"]] = lowCI_DFATrends
  DFA_Results[["FLoadings"]] = FLoadings
  DFA_Results[["bias"]] = bias
  DFA_Results[["biasAll"]] = biasAll
  DFA_Results[["RMSE"]] = RMSE
  DFA_Results[["biasAll.BT"]] = biasAll.BT
  DFA_Results[["RMSE.BT"]] = RMSE.BT
  DFA_Results[["GSD"]] = gsd
  DFA_Results[["datz_SD"]] = datz_SD
  
  save(DFA_Results, file="DFA_Results.RData")
  
  # write.csv(FitRatios, "FitRatios.csv")
  # 
  # colnames(FactorLoadings) = c("FL1","FL2","FL3","FL4","lowCI_FL1","lowCI_FL2","lowCI_FL3","lowCI_FL4","upCI_FL1","upCI_FL2","upCI_FL3","upCI_FL4")
  # write.csv(DFATrends, "DFATrends.csv")
  # write.csv(DFATrendsSE, "DFATrendsSE.csv")
  # write.csv(upCI_DFATrends, "upCI_DFATrends.csv")
  # write.csv(lowCI_DFATrends, "lowCI_DFATrends.csv")
  # 
  # apply(FitRatios,1,mean)
  # FLoadings = cbind(FactorLoadings, SUM = apply(FactorLoadings[,1:4], 1, sum))
  # write.csv(FLoadings, "FactorLoadings_sum.csv")
  # 
  # write.csv(bias, "bias.csv")
  # write.csv(biasAll, "biasAll.csv")
  # write.csv(RMSE, "RMSE.csv")
  # # write.csv(N_All, "N_All.csv")
  # # write.csv(N_STD, "N_STD.csv")
  
  
  
  ###################################  plotting ###################################
  png(filename="Trend.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    index = get(paste("index",l,sep=""))
    lowerCI = get(paste("lowerCI",l,sep=""))
    upperCI = get(paste("upperCI",l,sep=""))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, index, type='l', lwd=2, pch=16))
    # title(paste("Iteration",l,sep=" "))
    abline(h=0)
    NLscale = get(paste("NLscale",l,sep=""))
    lines(yrs,NLscale, type='l', col='red', lwd=2)
  }
  
  dev.off()
  
  
  png(filename="RawIndices.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  par(mfrow=c(10,10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    plot(1:years, Iy1[l,], type='l', col='black', xlim=c(min(yrs),max(yrs)), axes=F,
         ylim=c(0, max(c(Iy1[l,yrs], Iy2[l,yrs], Iy3[l,yrs]))+0.1))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    lines(1:years, Iy2[l,], type='l', col='blue')
    lines(1:years, Iy3[l,], type='l', col='red')
    # title(paste("Iteration",l, sep=" "))
  }
  
  dev.off()
  
  
  
  png(filename="FactorLoadings.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  # Before plotting factor loadings, we'll need to assign Z.rot, trends.rot
  par(mfrow=c(10,10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  survey = rownames(dat.z)
  minZ = 0.00
  for(l in 1:Niters){
    Z.rot = get(paste("Z.rot",l,sep="."))
    trends.rot = get(paste("trends.rot",l,sep="."))
    ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
    # par(mfrow=c(nrow(trends.rot),1), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
    for(h in 1:nrow(trends.rot)) {
      plot(c(1:N.ts)[abs(Z.rot[,h])>=minZ], as.vector(Z.rot[abs(Z.rot[,h])>=minZ,h]),
           type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1), axes=FALSE)
      axis(1, labels=FALSE)
      axis(2, labels=FALSE)
      #   for(j in 1:N.ts) {
      #     if(Z.rot[j,h] > minZ) {text(j, -0.05, survey[j], srt=90, adj=1, cex=0.9)}
      #     if(Z.rot[j,h] < -minZ) {text(j, 0.05, survey[j], srt=90, adj=0, cex=0.9)}
      #   } # end j loop
      #   mtext(paste("Factor loadings on trend",h,sep=" "),side=3,line=.5)
      abline(h=0, lwd=1, col="gray")
      abline(h=0.2, col='gray', lty=2)
      abline(h=-0.2, col='gray', lty=2)
    } # end h loop
  }
  
  dev.off()
  
  
  
  png(filename="Trend_scale.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    index = get(paste("index",l,sep=""))
    lowerCI = get(paste("lowerCI",l,sep=""))
    upperCI = get(paste("upperCI",l,sep=""))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, index, type='l', lwd=2, pch=16))
    # title(paste("Iteration",l,sep=" "))
    abline(h=0)
    NL = get(paste("NL",l,sep=""))
    lines(yrs, rescale(NL, to=c(min(index), max(index))), type='l',col="red", lwd=2)
  }
  
  
  dev.off()
  
  
  
  
  png(filename="BackTransformed_Trend_scale.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    indexBT = get(paste("indexBT",l,sep=""))
    indexSEBT = get(paste("indexSEBT",l,sep=""))
    lowerCI = indexBT-(1.96*indexSEBT)
    upperCI = indexBT+(1.96*indexSEBT)
    df <- data.frame(yrs, indexBT, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, indexBT, type='l', lwd=2, pch=16))
    # title(paste("Iteration",l,sep=" "))
    # abline(h=0)
    N = get(paste("N",l,sep=""))
    lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT))), type='l',col="red", lwd=2)
  }
  
  
  dev.off()
  
  
  png(filename="Trend_N.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    
    index = get(paste("index",l,sep=""))
    NL = get(paste("NL",l,sep=""))
    plot(yrs, NL, type='l',col="red", axes=FALSE)
    axis(1, labels=FALSE)
    axis(2, labels=FALSE)
    lines(yrs, rescale(index, to=c(min(NL),max(NL))))
    
  }
  
  dev.off()
  
  
  
  
  
  png(filename="Trend_Standardize.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    
    index = get(paste("index",l,sep=""))
    NL = get(paste("NL",l,sep=""))
    plot(yrs, (NL-mean(NL))/sd(NL), type='l',col="red", axes=FALSE, lwd=1.5)
    axis(1, labels=FALSE)
    axis(2, labels=FALSE)
    lines(yrs, (index-mean(index))/sd(index), type='l', lwd=1.5)
    
  }
  
  dev.off()
  
  
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    index = get(paste("index",l,sep=""))
    lowerCI = get(paste("lowerCI",l,sep=""))
    upperCI = get(paste("upperCI",l,sep=""))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, index, type='l', lwd=2, pch=16))
    # title(paste("Iteration",l,sep=" "))
    abline(h=0)
    NL = get(paste("NL",l,sep=""))
    lines(yrs, rescale(NL, to=c(min(index), max(index))), type='l',col="red", lwd=2)
  }
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    indexBT = get(paste("indexBT",l,sep=""))
    indexSEBT = get(paste("indexSEBT",l,sep=""))
    lowerCI = indexBT-(1.96*indexSEBT)
    upperCI = indexBT+(1.96*indexSEBT)
    df <- data.frame(yrs, indexBT, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, indexBT, type='l', lwd=2, pch=16))
    # title(paste("Iteration",l,sep=" "))
    # abline(h=0)
    N = get(paste("N",l,sep=""))
    lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT))), type='l',col="red", lwd=2)
  }
  
  # 
  # png(filename="TwoAxesCompare.png", 
  #     type="cairo",
  #     units="mm", 
  #     width=300, 
  #     height=300, 
  #     pointsize=12, 
  #     res=600)
  # 
  # par(mfrow=c(10, 10))
  # par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  # for(l in 1:Niters){
  #   index = get(paste("index",l,sep=""))
  #   lowerCI = get(paste("lowerCI",l,sep=""))
  #   upperCI = get(paste("upperCI",l,sep=""))
  #   df <- data.frame(yrs, index, lowerCI, upperCI)
  #   with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
  #   axis(1,labels=FALSE)
  #   axis(2,labels=FALSE)
  #   x <- as.numeric(df$yrs)
  #   polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
  #   with(df, lines(x, index, type='l', lwd=2, pch=16))
  #   # title(paste("Iteration",l,sep=" "))
  #   abline(h=0)
  #   par(new=T)
  #   NL = get(paste("NL",l,sep=""))
  #   plot(yrs,NL, type='l', col='red', axes=FALSE, lwd=1.5)
  #   axis(4, labels=FALSE)
  # }
  # 
  # dev.off()
  # 
  
  
  
  
  
  
  FLoadings2 = ifelse(is.na(FLoadings)==TRUE, 0, FLoadings)
  png(filename="TwoAxesCompare_TryCorrected.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.3),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    
    index = get(paste("index",l,sep=""))
    NL = get(paste("NL",l,sep=""))
    plot(yrs, NL, type='l',col="red", axes=FALSE, lwd=1.5)
    axis(1, labels=FALSE)
    axis(2, labels=FALSE)
    par(new=T)
    if(FLoadings2[l,10] >= 0){
      plot(yrs,index, type='l', col='black', axes=FALSE, lwd=1.5)
    }
    if(FLoadings2[l,10] < 0){
      plot(yrs,-1*index, type='l', col='black', axes=FALSE, lwd=1.5)
    }
    axis(4, labels=FALSE)
    
  }
  
  dev.off()
  
  
  
  
  ##################################################### END #################################################
  
  return(list(RMSE=RMSE, biasAll=biasAll, FitRatios=FitRatios, RMSE.BT=RMSE.BT, biasAll.BT=biasAll.BT))} # END FUNCTION
########################



# BTrial1 = list(RMSE=RMSE, biasAll=biasAll, FitRatios=FitRatios, RMSE.BT=RMSE.BT, biasAll.BT=biasAll.BT)


###################### LOAD SAVED TRIALS ##################


load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial7/ATrial7.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial7/DFA_Results.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial7/ResultsList.RData")

# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial1/ATrial1.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial2/ATrial2.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial3/ATrial3.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial4/ATrial4.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial5/ATrial5.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial6/ATrial6.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial7/ATrial7.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial8/ATrial8.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial9/ATrial9.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial10/ATrial10.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial11/ATrial11.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial12/ATrial12.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial13/ATrial13.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial14/ATrial14.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial15/ATrial15.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial16/ATrial16.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial17/ATrial17.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial18/ATrial18.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial19/ATrial19.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial20/ATrial20.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial21/ATrial21.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial22/ATrial22.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial23/ATrial23.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial24/ATrial24.RData")
# 
# 
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial1/BTrial1.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial2/BTrial2.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial3/BTrial3.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial4/BTrial4.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial5/BTrial5.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial6/BTrial6.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial7/BTrial7.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial8/BTrial8.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial9/BTrial9.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial10/BTrial10.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial11/BTrial11.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial12/BTrial12.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial13/BTrial13.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial14/BTrial14.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial15/BTrial15.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial16/BTrial16.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial17/BTrial17.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial18/BTrial18.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial19/BTrial19.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial20/BTrial20.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial21/BTrial21.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial22/BTrial22.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial23/BTrial23.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial24/BTrial24.RData")
# 
# 
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial4/CTrial4.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial5/CTrial5.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial6/CTrial6.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial10/CTrial10.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial11/CTrial11.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial12/CTrial12.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial16/CTrial16.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial17/CTrial17.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial18/CTrial18.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial22/CTrial22.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial23/CTrial23.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial24/CTrial24.RData")
# 
# 
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial4/DTrial4.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial5/DTrial5.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial6/DTrial6.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial10/DTrial10.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial11/DTrial11.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial12/DTrial12.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial16/DTrial16.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial17/DTrial17.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial18/DTrial18.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial22/DTrial22.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial23/DTrial23.RData")
# load("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial24/DTrial24.RData")




load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial1/ATrial1.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial2/ATrial2.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial3/ATrial3.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial4/ATrial4.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial5/ATrial5.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial6/ATrial6.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial7/ATrial7.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial8/ATrial8.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial9/ATrial9.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial10/ATrial10.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial11/ATrial11.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial12/ATrial12.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial13/ATrial13.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial14/ATrial14.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial15/ATrial15.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial16/ATrial16.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial17/ATrial17.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial18/ATrial18.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial19/ATrial19.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial20/ATrial20.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial21/ATrial21.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial22/ATrial22.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial23/ATrial23.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Base/Trial24/ATrial24.RData")


load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial1/BTrial1.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial2/BTrial2.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial3/BTrial3.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial4/BTrial4.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial5/BTrial5.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial6/BTrial6.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial7/BTrial7.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial8/BTrial8.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial9/BTrial9.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial10/BTrial10.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial11/BTrial11.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial12/BTrial12.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial13/BTrial13.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial14/BTrial14.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial15/BTrial15.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial16/BTrial16.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial17/BTrial17.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial18/BTrial18.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial19/BTrial19.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial20/BTrial20.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial21/BTrial21.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial22/BTrial22.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial23/BTrial23.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndices/Trial24/BTrial24.RData")


load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial4/CTrial4.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial5/CTrial5.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial6/CTrial6.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial10/CTrial10.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial11/CTrial11.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial12/CTrial12.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial16/CTrial16.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial17/CTrial17.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial18/CTrial18.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial22/CTrial22.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial23/CTrial23.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/Gradual/Trial24/CTrial24.RData")


load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial4/DTrial4.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial5/DTrial5.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial6/DTrial6.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial10/DTrial10.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial11/DTrial11.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial12/DTrial12.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial16/DTrial16.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial17/DTrial17.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial18/DTrial18.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial22/DTrial22.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial23/DTrial23.RData")
load("D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/FourIndicesGradual/Trial24/DTrial24.RData")



# D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019
# D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019
################ RUN TRIALS ###########################
setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial1")
# setwd("Y:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1")
# setwd("C:\\Users\\cpeterson\\Documents\\DFA_Simulation\\Atl_SN\\Base\\Trial1")
ATrial1 = DFA_Sim(F_var = as.vector(rep(0.2,years)), 
                  CV1=0.5,
                  CV2=0.5,
                  CV3=0.5,
                  q_survey1 = rep(0.001, length=years),
                  q_survey2 = rep(0.001, length=years),
                  q_survey3 = rep(0.001, length=years),
                  Iay_eq = Iay_eq_F1,
                  c=c(10, 10, 10))
save(ATrial1, file="ATrial1.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial2")
ATrial2 = DFA_Sim(F_var = as.vector(rep(0.2,years)), 
                  CV1=0.3, 
                  CV2=0.5,
                  CV3=0.7,
                  q_survey1 = rep(0.001, length=years),
                  q_survey2 = rep(0.001, length=years),
                  q_survey3 = rep(0.001, length=years),
                  Iay_eq = Iay_eq_F1,
                  c=c(25, 10, 70) )
save(ATrial2, file="ATrial2.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial3")
ATrial3 = DFA_Sim(F_var = as.vector(rep(0.2,years)), 
                  CV1=0.7,
                  CV2=0.5,
                  CV3=0.3,
                  q_survey1 = rep(0.001, length=years),
                  q_survey2 = rep(0.001, length=years),
                  q_survey3 = rep(0.001, length=years),
                  Iay_eq = Iay_eq_F1,
                  c=c(20, 50, 10) )
save(ATrial3, file="ATrial3.RData")




setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial4")
ATrial4 = DFA_Sim(F_var = as.vector(rep(0.2,years)), 
                  CV1=0.5,
                  CV2=0.5,
                  CV3=0.5,
                  q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                  q_survey2 = rep(0.001, length=years),
                  q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                  Iay_eq = Iay_eq_F1,
                  c=c(3, 10000000, 100) )
save(ATrial4, file="ATrial4.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial5")
ATrial5 = DFA_Sim(F_var = as.vector(rep(0.2,years)), 
                  CV1=0.3,
                  CV2=0.5,
                  CV3=0.7,
                  q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                  q_survey2 = rep(0.001, length=years),
                  q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                  Iay_eq = Iay_eq_F1,
                  c=c(3, 10000000, 200) )
save(ATrial5, file="ATrial5.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial6")
ATrial6 = DFA_Sim(F_var = as.vector(rep(0.2,years)), 
                  CV1=0.7,
                  CV2=0.5,
                  CV3=0.3,
                  q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                  q_survey2 = rep(0.001, length=years),
                  q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                  Iay_eq = Iay_eq_F1,
                  c=c(3, 10000000, 50) )
save(ATrial6, file="ATrial6.RData")


# save.image("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/ATL_SN_FEMALE_SurvivalSR_LogErr_DFA_Simulation_workspace.R.RData")





############### NEW BASE F; F2 #########################



setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial7")
ATrial7 = DFA_Sim(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                  CV1=0.5,
                  CV2=0.5,
                  CV3=0.5,
                  q_survey1 = rep(0.001, length=years),
                  q_survey2 = rep(0.001, length=years),
                  q_survey3 = rep(0.001, length=years),
                  Iay_eq = Iay_eq_F2,
                  c=c(12, 10, 10))
save(ATrial7, file="ATrial7.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial8")
ATrial8 = DFA_Sim(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                  CV1=0.3,
                  CV2=0.5,
                  CV3=0.7,
                  q_survey1 = rep(0.001, length=years),
                  q_survey2 = rep(0.001, length=years),
                  q_survey3 = rep(0.001, length=years),
                  Iay_eq = Iay_eq_F2,
                  c=c(12, 9, 10))
save(ATrial8, file="ATrial8.RData")



setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial9")
ATrial9 = DFA_Sim(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                  CV1=0.7,
                  CV2=0.5,
                  CV3=0.3,
                  q_survey1 = rep(0.001, length=years),
                  q_survey2 = rep(0.001, length=years),
                  q_survey3 = rep(0.001, length=years),
                  Iay_eq = Iay_eq_F2,
                  c=c(12, 9, 10) )
save(ATrial9, file="ATrial9.RData")







setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial10")
ATrial10 = DFA_Sim(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                   CV1=0.5,
                   CV2=0.5,
                   CV3=0.5,
                   q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                   Iay_eq = Iay_eq_F2,
                   c=c(12, 9, 10) )
save(ATrial10, file="ATrial10.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial11")
ATrial11 = DFA_Sim(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                   CV1=0.3,
                   CV2=0.5,
                   CV3=0.7,
                   q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                   Iay_eq = Iay_eq_F2,
                   c=c(12, 9, 10) ,
                   reorg=T)
save(ATrial11, file="ATrial11.RData")



setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial12")
ATrial12 = DFA_Sim(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                   CV1=0.7,
                   CV2=0.5,
                   CV3=0.3,
                   q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                   Iay_eq = Iay_eq_F2,
                   c=c(12, 9, 10))
save(ATrial12, file="ATrial12.RData")


# save.image("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/ATL_SN_FEMALE_SurvivalSR_LogErr_DFA_Simulation_workspace.R.RData")












################################# NEW BASE F3 #################################




setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial13")
ATrial13 = DFA_Sim(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                   CV1=0.5,
                   CV2=0.5,
                   CV3=0.5,
                   q_survey1 = rep(0.001, length=years),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = rep(0.001, length=years),
                   Iay_eq = Iay_eq_F3, 
                   c=c(25, 20, 20))
save(ATrial13, file="ATrial13.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial14")
ATrial14 = DFA_Sim(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                   CV1=0.3,
                   CV2=0.5,
                   CV3=0.7,
                   q_survey1 = rep(0.001, length=years),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = rep(0.001, length=years),
                   Iay_eq = Iay_eq_F3, 
                   c=c(25, 20, 20) )
save(ATrial14, file="ATrial14.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial15")
ATrial15 = DFA_Sim(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                   CV1=0.7,
                   CV2=0.5,
                   CV3=0.3,
                   q_survey1 = rep(0.001, length=years),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = rep(0.001, length=years),
                   Iay_eq = Iay_eq_F3, 
                   c=c(25, 20, 20) )
save(ATrial15, file="ATrial15.RData")




setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial16")
ATrial16 = DFA_Sim(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                   CV1=0.5,
                   CV2=0.5,
                   CV3=0.5,
                   q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                   Iay_eq = Iay_eq_F3,
                   c=c(50, 1500, 7))
save(ATrial16, file="ATrial16.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial17")
ATrial17 = DFA_Sim(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                   CV1=0.3,
                   CV2=0.5,
                   CV3=0.7,
                   q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                   Iay_eq = Iay_eq_F3,
                   c=c(50, 1500, 7))
save(ATrial17, file="ATrial17.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial18")
ATrial18 = DFA_Sim(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                   CV1=0.7,
                   CV2=0.5,
                   CV3=0.3,
                   q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                   Iay_eq = Iay_eq_F3,
                   c=c(50, 1500, 7) )
save(ATrial18, file="ATrial18.RData")



# save.image("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/ATL_SN_FEMALE_SurvivalSR_LogErr_DFA_Simulation_workspace.R.RData")









################################# NEW BASE F4  #################################


# F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))) # Trials 19-24

setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial19")
ATrial19 = DFA_Sim(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                   CV1=0.5,
                   CV2=0.5,
                   CV3=0.5,
                   q_survey1 = rep(0.001, length=years),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = rep(0.001, length=years),
                   Iay_eq = Iay_eq_F4,
                   c=c(15, 10, 10))
save(ATrial19, file="ATrial19.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial20")
ATrial20 = DFA_Sim(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                   CV1=0.3,
                   CV2=0.5,
                   CV3=0.7,
                   q_survey1 = rep(0.001, length=years),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = rep(0.001, length=years),
                   Iay_eq = Iay_eq_F4,
                   c=c(15, 10, 10))
save(ATrial20, file="ATrial20.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial21")
ATrial21 = DFA_Sim(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                   CV1=0.7,
                   CV2=0.5,
                   CV3=0.3,
                   q_survey1 = rep(0.001, length=years),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = rep(0.001, length=years),
                   Iay_eq = Iay_eq_F4,
                   c=c(15, 10, 10))
save(ATrial21, file="ATrial21.RData")



setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial22")
ATrial22 = DFA_Sim(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                   CV1=0.5,
                   CV2=0.5,
                   CV3=0.5,
                   q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                   Iay_eq = Iay_eq_F4,
                   c=c(7, 60, 4) ,
                   reorg=F)
save(ATrial22, file="ATrial22.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial23")
ATrial23 = DFA_Sim(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                   CV1=0.3,
                   CV2=0.5,
                   CV3=0.7,
                   q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                   Iay_eq = Iay_eq_F4,
                   c=c(7, 60, 4)  ,
                   reorg=T) #c=c(7, 50, 2.5))
save(ATrial23, file="ATrial23.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial24")
ATrial24 = DFA_Sim(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                   CV1=0.7,
                   CV2=0.5,
                   CV3=0.3,
                   q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                   Iay_eq = Iay_eq_F4,
                   c=c(7, 60, 4) ,
                   reorg=F)

save(ATrial24, file="ATrial24.RData")






















# save.image("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/ATL_SN_FEMALE_SurvivalSR_LogErr_DFA_Simulation_workspace.R.RData")
















###########################################################################################
########### 4 SURVEYS
###########################################################################################
# TRIALS ###############
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\\\FourIndices\\Trial24")
# Trial4_24 = DFA_Sim4(F_var = as.vector(c(rep(0,30), rep(0.4, 15), rep(0.2, 10), rep(0.05, 45))), 
#                   CV1=0.7,
#                   CV2=0.5,
#                   CV3=0.3,
#                   CV4=0.5,
#                   q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
#                   q_survey2 = rep(0.001, length=years),
#                   q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
#                   q_survey4 = rep(0.003, length=years))

# remove(ATrial1, ATrial2, ATrial3, ATrial4, ATrial5, ATrial6, ATrial7, ATrial8, ATrial9, ATrial10, ATrial11, ATrial12, ATrial13, ATrial14, ATrial15, ATrial16)

################ RUN TRIALS ###########################
setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial1")
BTrial1 = DFA_Sim4(F_var = as.vector(rep(0.2,years)), 
                   CV1=0.5,
                   CV2=0.5,
                   CV3=0.5,
                   CV4=0.5,
                   q_survey1 = rep(0.001, length=years),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = rep(0.001, length=years),
                   q_survey4 = rep(0.003, length=years),
                   Iay_eq = Iay_eq_F1,
                   c=c(10, 7, 9, 250) )
save(BTrial1, file="BTrial1.RData")




setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial2")
BTrial2 = DFA_Sim4(F_var = as.vector(rep(0.2,years)), 
                   CV1=0.3,
                   CV2=0.5,
                   CV3=0.7,
                   CV4=0.5,
                   q_survey1 = rep(0.001, length=years),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = rep(0.001, length=years),
                   q_survey4 = rep(0.003, length=years),
                   Iay_eq = Iay_eq_F1,
                   c=c(6, 7, 15, 250) )
save(BTrial2, file="BTrial2.RData")
names(BTrial2)

setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial3")
BTrial3 = DFA_Sim4(F_var = as.vector(rep(0.2,years)), 
                   CV1=0.7,
                   CV2=0.5,
                   CV3=0.3,
                   CV4=0.5,
                   q_survey1 = rep(0.001, length=years),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = rep(0.001, length=years),
                   q_survey4 = rep(0.003, length=years),
                   Iay_eq = Iay_eq_F1,
                   c=c(20, 7, 5, 200) )
save(BTrial3, file="BTrial3.RData")






setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial4")
BTrial4 = DFA_Sim4(F_var = as.vector(rep(0.2,years)), 
                   CV1=0.5,
                   CV2=0.5,
                   CV3=0.5,
                   CV4=0.5,
                   q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                   q_survey4 = rep(0.003, length=years),
                   Iay_eq = Iay_eq_F1,
                   c=c(1000000, 3, 15, 10) )
save(BTrial4, file="BTrial4.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial5")
BTrial5 = DFA_Sim4(F_var = as.vector(rep(0.2,years)), 
                   CV1=0.3,
                   CV2=0.5,
                   CV3=0.7,
                   CV4=0.5,
                   q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                   q_survey4 = rep(0.003, length=years),
                   Iay_eq = Iay_eq_F1,
                   c=c(1000000, 3, 15, 10) )
save(BTrial5, file="BTrial5.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial6")
BTrial6 = DFA_Sim4(F_var = as.vector(rep(0.2,years)), 
                   CV1=0.7,
                   CV2=0.5,
                   CV3=0.3,
                   CV4=0.5,
                   q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                   q_survey4 = rep(0.003, length=years),
                   Iay_eq = Iay_eq_F1,
                   c=c(1000000, 3, 15, 10) )
save(BTrial6, file="BTrial6.RData")



# save.image("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/ATL_SN_FEMALE_SurvivalSR_LogErr_DFA_Simulation_workspace.R.RData")







############### NEW BASE F; F2 #########################



setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial7")
BTrial7 = DFA_Sim4(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                   CV1=0.5,
                   CV2=0.5,
                   CV3=0.5,
                   CV4=0.5,
                   q_survey1 = rep(0.001, length=years),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = rep(0.001, length=years),
                   q_survey4 = rep(0.003, length=years),
                   Iay_eq = Iay_eq_F2,
                   c=c(20, 35, 20, 1.5))
save(BTrial7, file="BTrial7.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial8")
BTrial8 = DFA_Sim4(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                   CV1=0.3,
                   CV2=0.5,
                   CV3=0.7,
                   CV4=0.5,
                   q_survey1 = rep(0.001, length=years),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = rep(0.001, length=years),
                   q_survey4 = rep(0.003, length=years),
                   Iay_eq = Iay_eq_F2,
                   c=c(20, 35, 20, 1.5))
save(BTrial8, file="BTrial8.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial9")
BTrial9 = DFA_Sim4(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                   CV1=0.7,
                   CV2=0.5,
                   CV3=0.3,
                   CV4=0.5,
                   q_survey1 = rep(0.001, length=years),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = rep(0.001, length=years),
                   q_survey4 = rep(0.003, length=years),
                   Iay_eq = Iay_eq_F2,
                   c=c(20, 35, 20, 1.5))
save(BTrial9, file="BTrial9.RData")



setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial10")
BTrial10 = DFA_Sim4(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                    CV1=0.5,
                    CV2=0.5,
                    CV3=0.5,
                    CV4=0.5,
                    q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F2,
                    c=c(20, 35, 20, 1.5)  )
save(BTrial10, file="BTrial10.RData")



setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial11")
BTrial11 = DFA_Sim4(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                    CV1=0.3,
                    CV2=0.5,
                    CV3=0.7,
                    CV4=0.5,
                    q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F2,
                    c=c(20, 35, 20, 1.5)  ,
                    reorg=T)
save(BTrial11, file="BTrial11.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial12")
BTrial12 = DFA_Sim4(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                    CV1=0.7,
                    CV2=0.5,
                    CV3=0.3,
                    CV4=0.5,
                    q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F2,
                    c=c(20, 35, 20, 1.5)  )

save(BTrial12, file="BTrial12.RData")



# save.image("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/ATL_SN_FEMALE_SurvivalSR_LogErr_DFA_Simulation_workspace.R.RData")









################################# NEW BASE F3 #################################




setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial13")
BTrial13 = DFA_Sim4(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                    CV1=0.5,
                    CV2=0.5,
                    CV3=0.5,
                    CV4=0.5,
                    q_survey1 = rep(0.001, length=years),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = rep(0.001, length=years),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F3,
                    c=c(100, 150, 100, 3.5) )
save(BTrial13, file="BTrial13.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial14")
BTrial14 = DFA_Sim4(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                    CV1=0.3,
                    CV2=0.5,
                    CV3=0.7,
                    CV4=0.5,
                    q_survey1 = rep(0.001, length=years),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = rep(0.001, length=years),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F3,
                    c=c(100, 150, 100, 3.5))
save(BTrial14, file="BTrial14.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial15")
BTrial15 = DFA_Sim4(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                    CV1=0.7,
                    CV2=0.5,
                    CV3=0.3,
                    CV4=0.5,
                    q_survey1 = rep(0.001, length=years),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = rep(0.001, length=years),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F3,
                    c=c(100, 150, 100, 3.5))

save(BTrial15, file="BTrial15.RData")



setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial16")
BTrial16 = DFA_Sim4(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                    CV1=0.5,
                    CV2=0.5,
                    CV3=0.5,
                    CV4=0.5,
                    q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F3, 
                    c=c(30000, 150, 20, 3.5) )
save(BTrial16, file="BTrial16.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial17")
BTrial17 = DFA_Sim4(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                    CV1=0.3,
                    CV2=0.5,
                    CV3=0.7,
                    CV4=0.5,
                    q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F3, 
                    c=c(30000, 150, 20, 3.5) )
save(BTrial17, file="BTrial17.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial18")
BTrial18 = DFA_Sim4(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                    CV1=0.7,
                    CV2=0.5,
                    CV3=0.3,
                    CV4=0.5,
                    q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F3, 
                    c=c(30000, 150, 20, 3.5) )

save(BTrial18, file="BTrial18.RData")



# save.image("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/ATL_SN_FEMALE_SurvivalSR_LogErr_DFA_Simulation_workspace.R.RData")








################################# NEW BASE F4  #################################




setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial19")
BTrial19 = DFA_Sim4(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                    CV1=0.5,
                    CV2=0.5,
                    CV3=0.5,
                    CV4=0.5,
                    q_survey1 = rep(0.001, length=years),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = rep(0.001, length=years),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F4,
                    c=c(25, 30, 20, 10)  )
save(BTrial19, file="BTrial19.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial20")
BTrial20 = DFA_Sim4(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                    CV1=0.3,
                    CV2=0.5,
                    CV3=0.7,
                    CV4=0.5,
                    q_survey1 = rep(0.001, length=years),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = rep(0.001, length=years),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F4,
                    c=c(25, 30, 20, 10)  )
save(BTrial20, file="BTrial20.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial21")
BTrial21 = DFA_Sim4(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                    CV1=0.7,
                    CV2=0.5,
                    CV3=0.3,
                    CV4=0.5,
                    q_survey1 = rep(0.001, length=years),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = rep(0.001, length=years),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F4,
                    c=c(25, 30, 20, 10)  )

save(BTrial21, file="BTrial21.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial22")
BTrial22 = DFA_Sim4(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                    CV1=0.5,
                    CV2=0.5,
                    CV3=0.5,
                    CV4=0.5,
                    q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F4,
                    c=c(3000, 30, 20, 10)  ,
                    reorg=F)
save(BTrial22, file="BTrial22.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial23")
BTrial23 = DFA_Sim4(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                    CV1=0.3,
                    CV2=0.5,
                    CV3=0.7,
                    CV4=0.5,
                    q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F4,
                    c=c(3000, 30, 20, 10),
                    reorg=T )
save(BTrial23, file="BTrial23.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndices\\Trial24")
BTrial24 = DFA_Sim4(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                    CV1=0.7,
                    CV2=0.5,
                    CV3=0.3,
                    CV4=0.5,
                    q_survey1 =  c(rep(0.00025, length=years/2), rep(0.001, length=years/2)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, length=years/2), rep(0.002, length=years/2)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F4,
                    c=c(3000, 30, 20, 10),
                    reorg=F )
save(BTrial24, file="BTrial24.RData")



# save.image("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/ATL_SN_FEMALE_SurvivalSR_LogErr_DFA_Simulation_workspace.R.RData")










##################################################################################### 
##  GRADUAL SHIFT IN Q 
##################################################################################### 



################ RUN TRIALS ###########################
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial1")
# CTrial1 = DFA_Sim(F_var = as.vector(rep(0.2,years)), 
#                  CV1=0.5,
#                  CV2=0.5,
#                  CV3=0.5,
#                  q_survey1 = rep(0.001, length=years),
#                  q_survey2 = rep(0.001, length=years),
#                  q_survey3 = rep(0.001, length=years))
# 
# 
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial2")
# CTrial2 = DFA_Sim(F_var = as.vector(rep(0.2,years)), 
#                  CV1=0.3,
#                  CV2=0.5,
#                  CV3=0.7,
#                  q_survey1 = rep(0.001, length=years),
#                  q_survey2 = rep(0.001, length=years),
#                  q_survey3 = rep(0.001, length=years))
# 
# 
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial3")
# CTrial3 = DFA_Sim(F_var = as.vector(rep(0.2,years)), 
#                  CV1=0.7,
#                  CV2=0.5,
#                  CV3=0.3,
#                  q_survey1 = rep(0.001, length=years),
#                  q_survey2 = rep(0.001, length=years),
#                  q_survey3 = rep(0.001, length=years))



setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial4")
CTrial4 = DFA_Sim(F_var = as.vector(rep(0.2,years)), 
                  CV1=0.5,
                  CV2=0.5,
                  CV3=0.5,
                  q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                  q_survey2 = rep(0.001, length=years),
                  q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                  Iay_eq = Iay_eq_F1,
                  c=c(4, 5000000, 100) )
save(CTrial4, file="CTrial4.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial5")
CTrial5 = DFA_Sim(F_var = as.vector(rep(0.2,years)), 
                  CV1=0.3,
                  CV2=0.5,
                  CV3=0.7,
                  q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                  q_survey2 = rep(0.001, length=years),
                  q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                  Iay_eq = Iay_eq_F1,
                  c=c(4, 10000000, 500) )
save(CTrial5, file="CTrial5.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial6")
CTrial6 = DFA_Sim(F_var = as.vector(rep(0.2,years)), 
                  CV1=0.7,
                  CV2=0.5,
                  CV3=0.3,
                  q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                  q_survey2 = rep(0.001, length=years),
                  q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                  Iay_eq = Iay_eq_F1,
                  c=c(3.5, 5000000, 50) )

save(CTrial6, file="CTrial6.RData")









############### NEW BASE F; F2 #########################



# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial7")
# CTrial7 = DFA_Sim(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
#                  CV1=0.5,
#                  CV2=0.5,
#                  CV3=0.5,
#                  q_survey1 = rep(0.001, length=years),
#                  q_survey2 = rep(0.001, length=years),
#                  q_survey3 = rep(0.001, length=years))
# 
# 
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial8")
# CTrial8 = DFA_Sim(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
#                  CV1=0.3,
#                  CV2=0.5,
#                  CV3=0.7,
#                  q_survey1 = rep(0.001, length=years),
#                  q_survey2 = rep(0.001, length=years),
#                  q_survey3 = rep(0.001, length=years))
# 
# 
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial9")
# CTrial9 = DFA_Sim(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
#                  CV1=0.7,
#                  CV2=0.5,
#                  CV3=0.3,
#                  q_survey1 = rep(0.001, length=years),
#                  q_survey2 = rep(0.001, length=years),
#                  q_survey3 = rep(0.001, length=years))






setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial10")
CTrial10 = DFA_Sim(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                   CV1=0.5,
                   CV2=0.5,
                   CV3=0.5,
                   q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                   Iay_eq = Iay_eq_F2,
                   c=c(25, 2.25, 40) ) #c=c(500, 2.5, 800))
save(CTrial10, file="CTrial10.RData")

#*** CHALLENGING  b/c c's change! 
setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial11")
CTrial11 = DFA_Sim(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                   CV1=0.3,
                   CV2=0.5,
                   CV3=0.7,
                   q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                   Iay_eq = Iay_eq_F2,
                   c=c(25, 2.25, 40) ,
                   reorg=T) # c=c(500, 2.5, 800)
save(CTrial11, file="CTrial11.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial12")
CTrial12 = DFA_Sim(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                   CV1=0.7,
                   CV2=0.5,
                   CV3=0.3,
                   q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                   Iay_eq = Iay_eq_F2,
                   c=c(25, 2.25, 40) )
save(CTrial12, file="CTrial12.RData")


# save.image("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/ATL_SN_FEMALE_SurvivalSR_LogErr_DFA_Simulation_workspace.R.RData")











################################# NEW BASE F3 #################################




# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial13")
# CTrial13 = DFA_Sim(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
#                   CV1=0.5,
#                   CV2=0.5,
#                   CV3=0.5,
#                   q_survey1 = rep(0.001, length=years),
#                   q_survey2 = rep(0.001, length=years),
#                   q_survey3 = rep(0.001, length=years))
# 
# 
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial14")
# CTrial14 = DFA_Sim(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
#                   CV1=0.3,
#                   CV2=0.5,
#                   CV3=0.7,
#                   q_survey1 = rep(0.001, length=years),
#                   q_survey2 = rep(0.001, length=years),
#                   q_survey3 = rep(0.001, length=years))
# 
# 
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial15")
# CTrial15 = DFA_Sim(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
#                   CV1=0.7,
#                   CV2=0.5,
#                   CV3=0.3,
#                   q_survey1 = rep(0.001, length=years),
#                   q_survey2 = rep(0.001, length=years),
#                   q_survey3 = rep(0.001, length=years))






setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial16")
CTrial16 = DFA_Sim(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                   CV1=0.5,
                   CV2=0.5,
                   CV3=0.5,
                   q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                   Iay_eq = Iay_eq_F3, 
                   c=c(50, 1000, 5) )
save(CTrial16, file="CTrial16.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial17")
CTrial17 = DFA_Sim(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                   CV1=0.3,
                   CV2=0.5,
                   CV3=0.7,
                   q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                   Iay_eq = Iay_eq_F3, 
                   c=c(50, 1000, 5) )
save(CTrial17, file="CTrial17.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial18")
CTrial18 = DFA_Sim(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                   CV1=0.7,
                   CV2=0.5,
                   CV3=0.3,
                   q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                   Iay_eq = Iay_eq_F3, 
                   c=c(45, 1000, 5) )
save(CTrial18, file="CTrial18.RData")











################################# NEW BASE F4  #################################




# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial19")
# CTrial19 = DFA_Sim(F_var = as.vector(c(rep(0,30), rep(0.4, 15), rep(0.2, 10), rep(0.05, 45))), 
#                   CV1=0.5,
#                   CV2=0.5,
#                   CV3=0.5,
#                   q_survey1 = rep(0.001, length=years),
#                   q_survey2 = rep(0.001, length=years),
#                   q_survey3 = rep(0.001, length=years))
# 
# 
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial20")
# CTrial20 = DFA_Sim(F_var = as.vector(c(rep(0,30), rep(0.4, 15), rep(0.2, 10), rep(0.05, 45))), 
#                   CV1=0.3,
#                   CV2=0.5,
#                   CV3=0.7,
#                   q_survey1 = rep(0.001, length=years),
#                   q_survey2 = rep(0.001, length=years),
#                   q_survey3 = rep(0.001, length=years))
# 
# 
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial21")
# CTrial21 = DFA_Sim(F_var = as.vector(c(rep(0,30), rep(0.4, 15), rep(0.2, 10), rep(0.05, 45))), 
#                   CV1=0.7,
#                   CV2=0.5,
#                   CV3=0.3,
#                   q_survey1 = rep(0.001, length=years),
#                   q_survey2 = rep(0.001, length=years),
#                   q_survey3 = rep(0.001, length=years))




# F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))) # Trials 19-24
#


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial22")
CTrial22 = DFA_Sim(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                   CV1=0.5,
                   CV2=0.5,
                   CV3=0.5,
                   q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                   Iay_eq = Iay_eq_F4,
                   c=c(15, 250, 5),
                   reorg=T )
save(CTrial22, file="CTrial22.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial23")
CTrial23 = DFA_Sim(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                   CV1=0.3,
                   CV2=0.5,
                   CV3=0.7,
                   q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                   Iay_eq = Iay_eq_F4,
                   c=c(15, 250, 5),
                   reorg=T )
save(CTrial23, file="CTrial23.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Gradual\\Trial24")
CTrial24 = DFA_Sim(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                   CV1=0.7,
                   CV2=0.5,
                   CV3=0.3,
                   q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                   Iay_eq = Iay_eq_F4,
                   c=c(15, 250, 5),
                   reorg=F )
save(CTrial24, file="CTrial24.RData")



# save.image("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/ATL_SN_FEMALE_SurvivalSR_LogErr_DFA_Simulation_workspace.R.RData")




















################ RUN TRIALS ###########################
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial1")
# DTrial1 = DFA_Sim4(F_var = as.vector(rep(0.2,years)), 
#                   CV1=0.5,
#                   CV2=0.5,
#                   CV3=0.5,
#                   CV4=0.5,
#                   q_survey1 = rep(0.001, length=years),
#                   q_survey2 = rep(0.001, length=years),
#                   q_survey3 = rep(0.001, length=years),
#                   q_survey4 = rep(0.003, length=years))
# 
# 
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial2")
# DTrial2 = DFA_Sim4(F_var = as.vector(rep(0.2,years)), 
#                   CV1=0.3,
#                   CV2=0.5,
#                   CV3=0.7,
#                   CV4=0.5,
#                   q_survey1 = rep(0.001, length=years),
#                   q_survey2 = rep(0.001, length=years),
#                   q_survey3 = rep(0.001, length=years),
#                   q_survey4 = rep(0.003, length=years))
# 
# 
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial3")
# DTrial3 = DFA_Sim4(F_var = as.vector(rep(0.2,years)), 
#                   CV1=0.7,
#                   CV2=0.5,
#                   CV3=0.3,
#                   CV4=0.5,
#                   q_survey1 = rep(0.001, length=years),
#                   q_survey2 = rep(0.001, length=years),
#                   q_survey3 = rep(0.001, length=years),
#                   q_survey4 = rep(0.003, length=years))






setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial4")
DTrial4 = DFA_Sim4(F_var = as.vector(rep(0.2,years)), 
                   CV1=0.5,
                   CV2=0.5,
                   CV3=0.5,
                   CV4=0.5,
                   q_survey1 = c(rep(0.00025, 49), seq(0.0004, 0.0014, by=0.0001), rep(0.0015, 40)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.002, 49), seq(0.0018, 0.0006, by=-0.0002), rep(0.0005, 44)),
                   q_survey4 = rep(0.003, length=years),
                   Iay_eq = Iay_eq_F1,
                   c=c(5000000, 3, 50000, 10))
save(DTrial4, file="DTrial4.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial5")
DTrial5 = DFA_Sim4(F_var = as.vector(rep(0.2,years)), 
                   CV1=0.3,
                   CV2=0.5,
                   CV3=0.7,
                   CV4=0.5,
                   q_survey1 = c(rep(0.00025, 49), seq(0.0004, 0.0014, by=0.0001), rep(0.0015, 40)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.002, 49), seq(0.0018, 0.0006, by=-0.0002), rep(0.0005, 44)),
                   q_survey4 = rep(0.003, length=years),
                   Iay_eq = Iay_eq_F1,
                   c=c(5000000, 3, 50000, 10))
save(DTrial5, file="DTrial5.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial6")
DTrial6 = DFA_Sim4(F_var = as.vector(rep(0.2,years)), 
                   CV1=0.7,
                   CV2=0.5,
                   CV3=0.3,
                   CV4=0.5,
                   q_survey1 = c(rep(0.00025, 49), seq(0.0004, 0.0014, by=0.0001), rep(0.0015, 40)),
                   q_survey2 = rep(0.001, length=years),
                   q_survey3 = c(rep(0.002, 49), seq(0.0018, 0.0006, by=-0.0002), rep(0.0005, 44)),
                   q_survey4 = rep(0.003, length=years),
                   Iay_eq = Iay_eq_F1,
                   c=c(5000000, 3, 50000, 10))
save(DTrial6, file="DTrial6.RData")









############### NEW BASE F; F2 #########################



# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial7")
# DTrial7 = DFA_Sim4(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
#                   CV1=0.5,
#                   CV2=0.5,
#                   CV3=0.5,
#                   CV4=0.5,
#                   q_survey1 = rep(0.001, length=years),
#                   q_survey2 = rep(0.001, length=years),
#                   q_survey3 = rep(0.001, length=years),
#                   q_survey4 = rep(0.003, length=years))
# 
# 
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial8")
# DTrial8 = DFA_Sim4(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
#                   CV1=0.3,
#                   CV2=0.5,
#                   CV3=0.7,
#                   CV4=0.5,
#                   q_survey1 = rep(0.001, length=years),
#                   q_survey2 = rep(0.001, length=years),
#                   q_survey3 = rep(0.001, length=years),
#                   q_survey4 = rep(0.003, length=years))
# 
# 
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial9")
# DTrial9 = DFA_Sim4(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
#                   CV1=0.7,
#                   CV2=0.5,
#                   CV3=0.3,
#                   CV4=0.5,
#                   q_survey1 = rep(0.001, length=years),
#                   q_survey2 = rep(0.001, length=years),
#                   q_survey3 = rep(0.001, length=years),
#                   q_survey4 = rep(0.003, length=years))






setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial10")
DTrial10 = DFA_Sim4(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                    CV1=0.5,
                    CV2=0.5,
                    CV3=0.5,
                    CV4=0.5,
                    q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F2,
                    c=c(2.5, 40, 30, 1.5) )
save(DTrial10, file="DTrial10.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial11")
DTrial11 = DFA_Sim4(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                    CV1=0.3,
                    CV2=0.5,
                    CV3=0.7,
                    CV4=0.5,
                    q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F2,
                    c=c(2.5, 40, 30, 1.5) ,
                    reorg=T)
save(DTrial11, file="DTrial11.RData")




DFA_Sim4_d12 = function(F_var, CV1, CV2, CV3, CV4, q_survey1, q_survey2, q_survey3, q_survey4, Iay_eq=Iay_eq, c=5, reorg=F){
  # Iay
  par(mfrow=c(1,1))
  plot(1:years, apply(Iay_eq, 2, sum), type='l',ylim=c(0,max(apply(Iay_eq, 2, sum),na.rm=T)+max(apply(Iay_eq, 2, sum), na.rm=T)*0.1))
  abline(v=40, lwd=2, col="red")
  abline(v=65, lwd=2, col="red")
  # plot(1:years, apply(Iay, 2, sum), type='l')
  # Iy = apply(Iay, 2, sum)
  
  F_var=F_var
  CV1=CV1
  CV2=CV2
  CV3=CV3
  q_survey1=q_survey1
  q_survey2=q_survey2
  q_survey3=q_survey3
  c=c
  c1=c
  
  ####### ADD STOCHASTICITY ###########
  Niters = 100
  Iay1 = matrix(nrow=A, ncol=years)
  Iay2 = matrix(nrow=A, ncol=years)
  Iay3 = matrix(nrow=A, ncol=years)
  Iay4 = matrix(nrow=A, ncol=years)
  Iay1_CV = matrix(nrow=A, ncol=years)
  Iay2_CV = matrix(nrow=A, ncol=years)
  Iay3_CV = matrix(nrow=A, ncol=years)
  Iay4_CV = matrix(nrow=A, ncol=years)
  Cy = matrix(nrow=Niters, ncol=years)
  Ny = matrix(nrow=Niters, ncol=years)
  Iy1 = matrix(nrow=Niters, ncol=years)
  Iy2 = matrix(nrow=Niters, ncol=years)
  Iy3 = matrix(nrow=Niters, ncol=years)
  Iy4 = matrix(nrow=Niters, ncol=years)
  Iy1_CV = matrix(nrow=Niters, ncol=years)
  Iy2_CV = matrix(nrow=Niters, ncol=years)
  Iy3_CV = matrix(nrow=Niters, ncol=years)
  Iy4_CV = matrix(nrow=Niters, ncol=years)
  M0 = matrix(nrow=Niters, ncol=years)
  Npups = matrix(nrow=Niters, ncol=years)
  FSS = matrix(nrow=Niters, ncol=years) # Female spawning stock 
  Z0 = 1.904416
  Zmin = 0.209
  Npups_eq = 3320.671
  beta=0.5807613
  ResultsList = list()
  
  # N_All = vector()
  # N_STD = vector()
  
  set.seed(430)
  
  # F_const = rep(0, yrs)
  
  for(k in 1:Niters){
    # paste("Iteration",k,sep=" ")
    # Create matrices to store results
    Nay = matrix(nrow=A, ncol=years)
    Cay = matrix(nrow=A, ncol=years)
    # M0 = vector(length=years)
    
    # Inputs: 
    M = M_const
    # FM = F_const
    FM = F_var
    
    # Set up equilibrium conditions
    Nay[,1] = Equilibrium_Nay
    
    for (i in 1:(years-1)){
      
      matjit = rnorm(length(mat),mat,mat*0.01)
      matjit = ifelse(matjit>1, 1, matjit)
      FSS[k,i] = sum(Nay[,i]*matjit) #Fecund Stock size
      Npups[k,i] = sum(Nay[,i]*matjit*rnorm(length(fec), fec, fec*0.1)) # equilibrium Npups = 3320.671
      M0[k,i] = rnorm(1, ( ( (1- (Npups[k,i]/Npups_eq)^beta ) * (Zmin - Z0) ) + Z0 ) , 0.1)
      # M0[i] = 1.9044160 # THIS IS EQUILIBRIUM!!!
      Nay[1,i+1] = Npups[k,i]*exp(-M0[k,i])  # pups produced that survive to age 1 = Recruits!
      
      Mjit = rnorm(length(M), M, M*0.1) 
      
      Zay = vector(length=length(M))
      for (j in 1:A){
        Zay[j] = Mjit[j] + (Sel_catch[j,2]*FM[i]/4) + (Sel_catch[j,3]*FM[i]/4) + (Sel_catch[j,4]*FM[i]/4) + (Sel_catch[j,5]*FM[i]/4)
      } #end j loop Z
      
      
      for (j in 1:(A-2)){
        Nay[j+1,i+1] = Nay[j,i]*exp(-Zay[j])
        # Say[j,i]=Nay[j,i]*mat[j]
        # Pay[j,i]=Say[j,i]*fec[j]
        Cay[j,i]=Nay[j,i]*((Sel_catch[j,2]*FM[i]/4)/Zay[j])*(1-exp(-Zay[j])) +
          Nay[j,i]*((Sel_catch[j,3]*FM[i]/4)/Zay[j])*(1-exp(-Zay[j])) +
          Nay[j,i]*((Sel_catch[j,4]*FM[i]/4)/Zay[j])*(1-exp(-Zay[j])) +
          Nay[j,i]*((Sel_catch[j,5]*FM[i]/4)/Zay[j])*(1-exp(-Zay[j]))
        
      } # END J LOOP
      Nay[A,i+1] = Nay[A-1,i]*exp(-Zay[A-1]) +
        Nay[A,i]*exp(-Zay[A])
      
      Cay[A-1,i]=Nay[A-1,i]*((Sel_catch[A-1,2]*FM[i]/4)/Zay[A-1])*(1-exp(-Zay[A-1])) +
        Nay[A-1,i]*((Sel_catch[A-1,3]*FM[i]/4)/Zay[A-1])*(1-exp(-Zay[A-1])) +
        Nay[A-1,i]*((Sel_catch[A-1,4]*FM[i]/4)/Zay[A-1])*(1-exp(-Zay[A-1])) +
        Nay[A-1,i]*((Sel_catch[A-1,5]*FM[i]/4)/Zay[A-1])*(1-exp(-Zay[A-1]))
      
      Cay[A,i]=Nay[A,i]*((Sel_catch[A,2]*FM[i]/4)/Zay[A])*(1-exp(-Zay[A])) +
        Nay[A,i]*((Sel_catch[A,3]*FM[i]/4)/Zay[A])*(1-exp(-Zay[A])) +
        Nay[A,i]*((Sel_catch[A,4]*FM[i]/4)/Zay[A])*(1-exp(-Zay[A])) +
        Nay[A,i]*((Sel_catch[A,5]*FM[i]/4)/Zay[A])*(1-exp(-Zay[A]))
      
      
      
      
      for(h in 1:A){
        Iay1_CV[h,i] = runif(1, min=CV1-0.1, max=CV1+0.1)
        Iay1[h,i] = max(q_survey1[i]*Sel_survey[h,2] * Nay[h,i] * 
                          exp( rnorm(1, 0, (q_survey1[i]*Sel_survey[h,2]*Nay[h,i])*Iay1_CV[h,i]) - 
                                 ( ((q_survey1[i]*Sel_survey[h,2]*Nay[h,i]*Iay1_CV[h,i])^2)/2 ) ), 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay2_CV[h,i] = runif(1, min=CV2-0.1, max=CV2+0.1)
        # Iay2[h,i] = max(q_survey2[i]*Sel_survey[h,3] * Nay[h,i] +  rnorm(1, 0, (q_survey2[i]*Sel_survey[h,3]*Nay[h,i])*Iay2_CV[h,i] ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay2[h,i] = max(q_survey2[i]*Sel_survey[h,3] * Nay[h,i] * 
                          exp( rnorm(1, 0, (q_survey2[i]*Sel_survey[h,3]*Nay[h,i])*Iay2_CV[h,i]) - 
                                 ( ((q_survey2[i]*Sel_survey[h,3]*Nay[h,i]*Iay2_CV[h,i])^2)/2 ) ), 0) #LOTS OF UNCERTAINTY!based off CV;
        
        
        Iay3_CV[h,i] = runif(1, min=CV3-0.1, max=CV3+0.1)
        # Iay3[h,i] = max(q_survey3[i]*Sel_survey[h,2] * Nay[h,i] +  rnorm(1, 0, (q_survey3[i]*Sel_survey[h,2]*Nay[h,i])*Iay3_CV[h,i] ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay3[h,i] = max(q_survey3[i]*Sel_survey[h,2] * Nay[h,i] * 
                          exp( rnorm(1, 0, (q_survey3[i]*Sel_survey[h,2]*Nay[h,i])*Iay3_CV[h,i]) - 
                                 ( ((q_survey3[i]*Sel_survey[h,2]*Nay[h,i]*Iay3_CV[h,i])^2)/2 ) ), 0) #LOTS OF UNCERTAINTY!based off CV;
        
        Iay4_CV[h,i] = runif(1, min=CV4-0.1, max=CV4+0.1)
        # Iay4[h,i] = max(q_survey4[i]*Sel_survey[h,2] * Nay[h,i] +  rnorm(1, 0, (q_survey4[i]*Sel_survey[h,2]*Nay[h,i])*Iay4_CV[h,i] ), 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay4[h,i] = max(q_survey4[i]*Sel_survey[h,2] * Nay[h,i] * 
                          exp( rnorm(1, 0, (q_survey4[i]*Sel_survey[h,2]*Nay[h,i])*Iay4_CV[h,i]) - 
                                 ( ((q_survey4[i]*Sel_survey[h,2]*Nay[h,i]*Iay4_CV[h,i])^2)/2 ) ), 0) #LOTS OF UNCERTAINTY!based off CV;
        
      } # end h
      
      
      # LCay1 = von bert curve to calculate length composition
      # random q: runif(1,max=0.005, min=0.0005)
    }# ends i loop
    
    # USE VON BERT CURVE TO ESTIMATE LENGTH COMP. Let von Bert parameters vary by Niter...???? 
    #   year (implement in i loop)... 
    #   Or should they vary by year class? Could create a matrix of von bert parameters based on age class 
    #       and pull out appropriate values as needed
    # Lt = Linf * (1 - exp(-K*(t - t0)))
    Linf=81.6
    K=0.49
    t0=-0.97
    
    # L0 = Linf * (1 - exp(-K*(0-t0)))
    
    
    # GET LENGTH FREQUENCIES
    Cay_N = cbind(ages, round(Cay))
    assign(paste("LF_Cay",k,sep="_"), data.frame(LBins = seq(20, 122, by=2)))
    for(f in 40:65){
      dfC = as.data.frame(Cay_N[,c(1,f+1)])
      colnames(dfC) = c("ages","count") 
      df.expandedC = dfC[rep(row.names(dfC), dfC$count),1:2]
      vonBertC = rbvn(nrow(df.expandedC),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expandedC$L = vonBertC$Linf * ( 1 - exp(-(vonBertC$K) * (jitter(df.expandedC$ages,amount=0.5) - t0) ))
      hstC=hist(df.expandedC$L, breaks=c(seq(20, 124, by=2)), plot=F)
      assign(paste("LF_Cay",k,sep="_"), cbind(get(paste("LF_Cay",k,sep="_")), hstC$counts))
    } # end f loop
    # write.csv(get(paste("LF_Cay",k,sep="_")), paste(paste("LF_Cay",k,sep="_"),".csv",sep=""))
    
    if(k==27){rnorm(1)}
    Iay1_N = cbind(ages, round(Iay1*100))
    assign(paste("LF_Iay1",k,sep="_"), data.frame(LBins = seq(20, 122, by=2)))
    for(l in 40:65){
      df1 = as.data.frame(Iay1_N[,c(1,l+1)])
      colnames(df1) = c("ages","count") 
      df.expanded1 = df1[rep(row.names(df1), df1$count),1:2]
      vonBert1 = rbvn(nrow(df.expanded1),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expanded1$L = vonBert1$Linf * ( 1 - exp(-(vonBert1$K) * (jitter(df.expanded1$ages,amount=0.5) - t0) ))
      hst1=hist(df.expanded1$L, breaks=c(seq(20, 124, by=2)), plot=F)
      assign(paste("LF_Iay1",k,sep="_"), cbind(get(paste("LF_Iay1",k,sep="_")), hst1$counts))
    } # end l loop
    # write.csv(get(paste("LF_Iay1",k,sep="_")), paste(paste("LF_Iay1",k,sep="_"),".csv",sep=""))
    
    Iay2_N = cbind(ages, round(Iay2*100))
    assign(paste("LF_Iay2",k,sep="_"), data.frame(LBins = seq(20, 122, by=2)))
    for(m in 40:65){
      df2 = as.data.frame(Iay2_N[,c(1,m+1)])
      colnames(df2) = c("ages","count") 
      df.expanded2 = df2[rep(row.names(df2), df2$count),1:2]
      vonBert2 = rbvn(nrow(df.expanded2),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expanded2$L = vonBert2$Linf * ( 1 - exp(-(vonBert2$K) * (jitter(df.expanded2$ages,amount=0.5) - t0) ))
      hst2=hist(df.expanded2$L, breaks=c(seq(20, 124, by=2)), plot=F)
      assign(paste("LF_Iay2",k,sep="_"), cbind(get(paste("LF_Iay2",k,sep="_")), hst2$counts))
    } # end m loop
    # write.csv(get(paste("LF_Iay2",k,sep="_")), paste(paste("LF_Iay2",k,sep="_"),".csv",sep=""))
    
    
    Iay3_N = cbind(ages, round(Iay3*100))
    assign(paste("LF_Iay3",k,sep="_"), data.frame(LBins = seq(20, 122, by=2)))
    for(n in 40:65){
      df3 = as.data.frame(Iay3_N[,c(1,n+1)])
      colnames(df3) = c("ages","count") 
      df.expanded3 = df3[rep(row.names(df3), df3$count),1:2]
      vonBert3 = rbvn(nrow(df.expanded3),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expanded3$L = vonBert3$Linf * ( 1 - exp(-(vonBert3$K) * (jitter(df.expanded3$ages,amount=0.5) - t0) ))
      hst3=hist(df.expanded3$L, breaks=c(seq(20, 124, by=2)), plot=F)
      assign(paste("LF_Iay3",k,sep="_"), cbind(get(paste("LF_Iay3",k,sep="_")), hst3$counts))
    } # end n loop
    # write.csv(get(paste("LF_Iay3",k,sep="_")), paste(paste("LF_Iay3",k,sep="_"),".csv",sep=""))
    
    Iay4_N = cbind(ages, round(Iay4*100))
    assign(paste("LF_Iay4",k,sep="_"), data.frame(LBins = seq(20, 122, by=2)))
    for(u in 40:65){
      df4 = as.data.frame(Iay4_N[,c(1,u+1)])
      colnames(df4) = c("ages","count") 
      df.expanded4 = df4[rep(row.names(df4), df4$count),1:2]
      vonBert4 = rbvn(nrow(df.expanded4),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expanded4$L = vonBert4$Linf * ( 1 - exp(-(vonBert4$K) * (jitter(df.expanded4$ages,amount=0.5) - t0) ))
      hst4=hist(df.expanded4$L, breaks=c(seq(20, 124, by=2)), plot=F)
      assign(paste("LF_Iay4",k,sep="_"), cbind(get(paste("LF_Iay4",k,sep="_")), hst4$counts))
    } # end u loop
    # write.csv(get(paste("LF_Iay4",k,sep="_")), paste(paste("LF_Iay4",k,sep="_"),".csv",sep=""))
    
    
    
    #save length-frequencies
    LFreqList = list()
    LFreqList[[paste("LF_Cay",k,sep="_")]] = get(paste("LF_Cay",k,sep="_"))
    LFreqList[[paste("LF_Iay1",k,sep="_")]] = get(paste("LF_Iay1",k,sep="_"))
    LFreqList[[paste("LF_Iay2",k,sep="_")]] = get(paste("LF_Iay2",k,sep="_"))
    LFreqList[[paste("LF_Iay3",k,sep="_")]] = get(paste("LF_Iay3",k,sep="_"))
    LFreqList[[paste("LF_Iay4",k,sep="_")]] = get(paste("LF_Iay4",k,sep="_"))
    
    save(LFreqList, file=paste("LFreqList_",k,".RData",sep=""))
    
    
    # Linf = rnorm(years, 81.9, 81.9*0.1) # CV = 0.2; 0.2 = sd / 81.9
    # K = rnorm(years, 0.48, 0.48*0.1)
    # t0 = rnorm(years, -0.99, 0.99*0.1)
    # LC = Linf * ( 1 - exp(-K*(c(0,ages) - t0)))
    # assign(paste("LC",k,sep="_"), cbind(c(0,ages), LC))
    # write.csv(LC, paste(paste("LC",k,sep="_"),".csv",sep=""))
    # assign(paste("M0",k,sep="_"), M0)
    
    
    
    # write.csv(Iay1, paste(paste("Iay1",k,sep="_"),".csv",sep=""))
    # write.csv(Iay2, paste(paste("Iay2",k,sep="_"),".csv",sep=""))
    # write.csv(Iay3, paste(paste("Iay3",k,sep="_"),".csv",sep=""))
    # write.csv(Iay4, paste(paste("Iay4",k,sep="_"),".csv",sep=""))
    
    Iy1[k,] = apply(Iay1, 2, sum)
    Iy2[k,] = apply(Iay2, 2, sum)
    Iy3[k,] = apply(Iay3, 2, sum)
    Iy4[k,] = apply(Iay4, 2, sum)
    lines(1:years, Iy1[k,], type='l', col='grey')
    lines(1:years, Iy2[k,], type='l', col='blue')
    lines(1:years, Iy3[k,], type='l', col='red')
    lines(1:years, Iy4[k,], type='l', col='green')
    Iy1_CV[k,] = apply(Iay1_CV, 2, mean, na.rm=T)
    Iy2_CV[k,] = apply(Iay2_CV, 2, mean, na.rm=T)
    Iy3_CV[k,] = apply(Iay3_CV, 2, mean, na.rm=T)
    Iy4_CV[k,] = apply(Iay4_CV, 2, mean, na.rm=T)
    
    Cy[k,] = apply(Cay, 2, sum)
    Ny[k,] = apply(Nay, 2, sum)
    
    
    # plot(1:years, colSums(Nay),type='l', ylim=c(0,max(colSums(Nay))+20))
    # lines(1:years, colSums(Nay),type='l',col='black')
    # paste("Nay_",k, sep='')
    
    #write 
    # write.csv(Nay, paste(paste("Nay", k, sep="_"),".csv",sep=""))
    assign(paste("Nay", k, sep="_"), Nay)
    # write.csv(Cay, paste(paste("Cay", k, sep="_"),".csv",sep=""))
    assign(paste("Cay", k, sep="_"), Cay)
    
    
    #write 
    AnnualNC = list()
    AnnualNC[[paste("Nay",k,sep="_")]] = Nay
    AnnualNC[[paste("Cay",k,sep="_")]] = Cay
    AnnualNC[[paste("Iay1",k,sep="_")]] = Iay1
    AnnualNC[[paste("Iay2",k,sep="_")]] = Iay2
    AnnualNC[[paste("Iay3",k,sep="_")]] = Iay3
    AnnualNC[[paste("Iay4",k,sep="_")]] = Iay4
    
    save(AnnualNC, file=paste("AnnualNC_",k,".RData",sep=""))
    
    
    
  } # end k loop
  
  
  
  ResultsList[["Iy1"]] = Iy1
  ResultsList[["Iy2"]] = Iy2
  ResultsList[["Iy3"]] = Iy3
  ResultsList[["Iy4"]] = Iy4
  
  ResultsList[["Iy1_CV"]] = Iy1_CV
  ResultsList[["Iy2_CV"]] = Iy2_CV
  ResultsList[["Iy3_CV"]] = Iy3_CV
  ResultsList[["Iy4_CV"]] = Iy4_CV
  
  
  ResultsList[["M0"]] = M0
  ResultsList[["FSS"]] = FSS
  ResultsList[["Npups"]] = Npups
  
  ResultsList[["Cy"]] = Cy
  ResultsList[["Ny"]] = Ny
  
  save(ResultsList, file="ResultsList.RData")
  
  # write.csv(Iy1, "Iy1.csv")
  # write.csv(Iy2, "Iy2.csv")
  # write.csv(Iy3, "Iy3.csv")
  # write.csv(Iy4, "Iy4.csv")
  # write.csv(M0, "M0.csv") # populates where each row is an iteration and each column is one year
  # write.csv(FSS, "FSS.csv") # populates where each row is an iteration and each column is one year
  # write.csv(Iy1_CV, "Iy1_CV.csv")
  # write.csv(Iy2_CV, "Iy2_CV.csv")
  # write.csv(Iy3_CV, "Iy3_CV.csv")
  # write.csv(Iy4_CV, "Iy4_CV.csv")
  # write.csv(Cy, "Cy.csv")
  
  
  ########## TEST DFA #############
  
  # library(MARSS)
  yrs = 40:(years-35)
  biasAll = matrix(nrow=Niters, ncol=length(yrs))
  bias = vector(length=Niters)
  biasAll.BT = matrix(nrow=Niters, ncol=length(yrs))
  RMSE.BT = vector(length=Niters)
  
  RMSE = vector(length=Niters)
  FitRatios = matrix(nrow=Niters, ncol=4)
  FactorLoadings=vector()
  DFATrends = vector()
  DFATrendsBT = vector()
  DFATrendsSE = vector()
  DFATrendsSEBT = vector()
  upCI_DFATrends = vector()
  lowCI_DFATrends = vector()
  
  datz_SD = matrix(nrow=Niters, ncol=4)
  # i=1
  
  for(i in 1:Niters){
    
    if(reorg == T){
      assign(paste("dat",i,sep=""), rbind(Iy2[i,40:(years-35)], Iy1[i,40:(years-35)], Iy3[i,40:(years-35)], c(rep(NA,15), Iy4[i,55:(years-35)])))
      c=c( c1[2], c1[1], c1[3], c1[4])
      
    }
    
    if(reorg==F){
      assign(paste("dat",i,sep=""), rbind(Iy1[i,40:(years-35)], Iy2[i,40:(years-35)], Iy3[i,40:(years-35)], c(rep(NA,15), Iy4[i,55:(years-35)])))
    }
    
    
    dat.a = get(paste("dat",i,sep=""))
    # c=c(40, 45, 60, 4) 
    dat=dat.a*c
    
    
    TT=ncol(dat)
    N.ts = nrow(dat)
    
    datL = log(dat)
    y.bar = apply(datL, 1, mean, na.rm=TRUE)
    
    dat.dm = (datL - y.bar) / y.bar
    gsd = sd(c(dat.dm),na.rm=TRUE)
    dat.z = dat.dm/gsd
    # gsd
    # y.bar
    # apply(dat.z, 1, mean, na.rm=T)
    # apply(dat.z, 1, sd, na.rm=T)
    datz_SD[i,] = apply(dat.z, 1, sd, na.rm=T)
    
    # dat.dm = (datL - y.bar) 
    dat.z = as.matrix(dat.z)
    rownames(dat.z) = c("Survey1","Survey2","Survey3","Survey4")
    
    
    
    ##### SN DELTA LOGNORMAL ####
    cntl.list = list(minit=200, maxit=500000, abstol=0.02, allow.degen=FALSE, conv.test.slope.tol = 0.05)
    R = diag(c(CV1, CV2, CV3, CV4), nrow=4, ncol=4)
    # dfa = MARSS(dat.z, model = list(m=1,R="diagonal and equal"), control=cntl.list, form="dfa", z.score=TRUE)
    # dfa = MARSS(dat.z, model = list(m=1,R="diagonal and equal"), control=cntl.list, form="dfa", z.score=FALSE)
    dfa = MARSS(dat.z, model = list(m=1,R=R), control=cntl.list, form="dfa", z.score=FALSE)
    
    pars = MARSSparamCIs(dfa, nboot=10000)
    # MARSSparamCIs(Mdglm.FINAL, nboot=10000, alpha=0.1)
    # MARSSparamCIs(Mdglm.FINAL, nboot=10000, alpha=0.15)
    
    # get the inverse of the rotation matrix
    #H.cMdglm.inv = varimax(coef(Mdglm.FINAL, type="matrix")$Z)$rotmat
    Z.rot = coef(dfa, type="matrix")$Z #%*% H.cMdglm.inv
    trends.rot = dfa$states
    ts.trends = t(trends.rot)
    par.mat=coef(dfa, type="matrix")
    fit.b = par.mat$Z %*% dfa$states # + matrix(par.cMdglm.mat$D, nrow=N.ts) %*% SBnao
    
    assign(paste("Z.rot", i, sep="."), Z.rot)
    assign(paste("Z.upCI", i, sep="."), pars$par.upCI$Z)
    assign(paste("Z.lowCI", i, sep="."), pars$par.lowCI$Z)
    assign(paste("trends.rot", i, sep="."), trends.rot)
    
    
    # 
    # # plot factor loadings
    par(mfrow=c(3,2))
    survey = rownames(dat.z)
    minZ = 0.00
    ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
    # par(mfrow=c(nrow(trends.rot),1), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
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
    
    
    # PLOT WITH CI's
    yrs = 40:(years-35)
    
    index = dfa$states[1,]
    indexSE = dfa$states.se[1,]
    lowerCI = dfa$states[1,]-1.96*dfa$states.se[1,]
    upperCI = dfa$states[1,]+1.96*dfa$states.se[1,]
    assign(paste("index",i,sep=""),index)
    assign(paste("indexSE",i,sep=""),indexSE)
    assign(paste("lowerCI",i,sep=""),lowerCI)
    assign(paste("upperCI",i,sep=""),upperCI)
    
    # par(mfrow=c(2,2))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='Year', ylab='Abundance', ylim=c(min(lowerCI), max(upperCI))))
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, index, type='l', lwd=2, pch=16))
    title(paste("Iteration",i,sep=" "))
    abline(h=0)
    
    N = Ny[i,40:65]
    NL=log(N)
    Nscale = (N-mean(N))/sd(N)
    NLscale = (N-mean(N))/sd(N)
    index.z = (index - mean(index))/sd(index)
    index.z = ifelse(is.na(index.z), 0, index.z)
    index.z = ifelse(abs(index.z)==Inf, 0, index.z)
    
    
    indexSEBT = indexSE*gsd
    indexBT = exp(index*gsd+ ((indexSEBT^2)/2))
    assign(paste("indexBT",i,sep=""),indexBT)
    assign(paste("indexSEBT",i,sep=""),indexSEBT)
    
    lines(yrs,rescale(NL,to=c(min(index),max(index))), type='l', col='red', lwd=2)
    
    assign(paste("N",i,sep=""), N)
    assign(paste("Nscale",i,sep=""),Nscale)
    assign(paste("NL",i,sep=""), NL)
    assign(paste("NLscale",i,sep=""),NLscale)
    
    # N_All = rbind(N_All, N)
    # N_STD = rbind(N_STD, Nscale)
    
    biasAll[i,] = index.z - NLscale
    bias[i] = mean(index.z - NLscale)
    RMSE[i] = sqrt(sum((index.z-NLscale)^2)/length(index.z))
    
    indexBT.z = (indexBT-mean(indexBT))/sd(indexBT)
    biasAll.BT[i,] = indexBT.z - Nscale
    RMSE.BT[i] = sqrt(sum((indexBT.z-Nscale)^2)/length(indexBT.z))
    
    
    # RMSE_yr[i] = sqrt(sum((Nscale - indexscale)^2)/length(indexscale))
    
    # Plot fitted values
    survey = rownames(dat.z)
    # par(mfcol=c(ceiling(N.ts/2),2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
    for(n in 1:length(survey)){
      plot(yrs,dat.z[n,],xlab="",ylab="abundance index",bty="L", xaxt="n", ylim=c(-4,4), pch=16, col="blue")
      axis(1,labels=TRUE)
      lines(yrs,fit.b[n,], lwd=2)
      lines(yrs,rescale(NL, to=c(min(fit.b[n,]), max(fit.b[n,]))), type='l', col='red', lwd=2)
      title(paste("Atlantic sharpnose",survey[n], sep=" "))
    }
    
    
    #### ASSESS MODEL FITS
    sumResids = rowSums((dat.z-fit.b)^2, na.rm=TRUE)
    sumObserved = rowSums(dat.z^2, na.rm=TRUE)
    FitRatio = sumResids / sumObserved ; FitRatio
    FitRatios[i,] = FitRatio
    # mean(FitRatio)
    # indexBT
    # dat.a
    
    FactorLoadings = rbind(FactorLoadings, c(get(paste("Z.rot",i,sep=".")),get(paste("Z.lowCI",i,sep=".")),get(paste("Z.upCI",i,sep="."))))
    DFATrends = rbind(DFATrends, index)
    DFATrendsBT = rbind(DFATrendsBT, indexBT)
    DFATrendsSE = rbind(DFATrendsSE, indexSE)
    DFATrendsSEBT = rbind(DFATrendsSEBT, indexSEBT)
    upCI_DFATrends = rbind(upCI_DFATrends, upperCI)
    lowCI_DFATrends = rbind(lowCI_DFATrends, lowerCI)
    
    
    # par(mfrow=c(1,1))
    # plot(fit.b, dat.z-fit.b)
  } # END LOOP
  
  
  
  
  colnames(FactorLoadings) = c("FL1","FL2","FL3","FL4",
                               "lowCI_FL1","lowCI_FL2","lowCI_FL3","lowCI_FL4",
                               "upCI_FL1","upCI_FL2","upCI_FL3","upCI_FL4")
  FR = apply(FitRatios,1,mean)
  FLoadings = cbind(FactorLoadings, meanFactorLoading = FR)
  
  DFA_Results = list()
  DFA_Results[["FitRatios"]] = FitRatios
  DFA_Results[["DFATrends"]] = DFATrends
  DFA_Results[["DFATrendsBT"]] = DFATrendsBT
  DFA_Results[["DFATrendsSE"]] = DFATrendsSE
  DFA_Results[["DFATrendsSEBT"]] = DFATrendsSEBT
  DFA_Results[["upCI_DFATrends"]] = upCI_DFATrends
  DFA_Results[["lowCI_DFATrends"]] = lowCI_DFATrends
  DFA_Results[["FLoadings"]] = FLoadings
  DFA_Results[["bias"]] = bias
  DFA_Results[["biasAll"]] = biasAll
  DFA_Results[["RMSE"]] = RMSE
  DFA_Results[["biasAll.BT"]] = biasAll.BT
  DFA_Results[["RMSE.BT"]] = RMSE.BT
  DFA_Results[["GSD"]] = gsd
  DFA_Results[["datz_SD"]] = datz_SD
  
  save(DFA_Results, file="DFA_Results.RData")
  
  # write.csv(FitRatios, "FitRatios.csv")
  # 
  # colnames(FactorLoadings) = c("FL1","FL2","FL3","FL4","lowCI_FL1","lowCI_FL2","lowCI_FL3","lowCI_FL4","upCI_FL1","upCI_FL2","upCI_FL3","upCI_FL4")
  # write.csv(DFATrends, "DFATrends.csv")
  # write.csv(DFATrendsSE, "DFATrendsSE.csv")
  # write.csv(upCI_DFATrends, "upCI_DFATrends.csv")
  # write.csv(lowCI_DFATrends, "lowCI_DFATrends.csv")
  # 
  # apply(FitRatios,1,mean)
  # FLoadings = cbind(FactorLoadings, SUM = apply(FactorLoadings[,1:4], 1, sum))
  # write.csv(FLoadings, "FactorLoadings_sum.csv")
  # 
  # write.csv(bias, "bias.csv")
  # write.csv(biasAll, "biasAll.csv")
  # write.csv(RMSE, "RMSE.csv")
  # # write.csv(N_All, "N_All.csv")
  # # write.csv(N_STD, "N_STD.csv")
  
  
  
  ###################################  plotting ###################################
  png(filename="Trend.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    index = get(paste("index",l,sep=""))
    lowerCI = get(paste("lowerCI",l,sep=""))
    upperCI = get(paste("upperCI",l,sep=""))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, index, type='l', lwd=2, pch=16))
    # title(paste("Iteration",l,sep=" "))
    abline(h=0)
    NLscale = get(paste("NLscale",l,sep=""))
    lines(yrs,NLscale, type='l', col='red', lwd=2)
  }
  
  dev.off()
  
  
  png(filename="RawIndices.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  par(mfrow=c(10,10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    plot(1:years, Iy1[l,], type='l', col='black', xlim=c(min(yrs),max(yrs)), axes=F,
         ylim=c(0, max(c(Iy1[l,yrs], Iy2[l,yrs], Iy3[l,yrs]))+0.1))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    lines(1:years, Iy2[l,], type='l', col='blue')
    lines(1:years, Iy3[l,], type='l', col='red')
    # title(paste("Iteration",l, sep=" "))
  }
  
  dev.off()
  
  
  
  png(filename="FactorLoadings.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  # Before plotting factor loadings, we'll need to assign Z.rot, trends.rot
  par(mfrow=c(10,10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  survey = rownames(dat.z)
  minZ = 0.00
  for(l in 1:Niters){
    Z.rot = get(paste("Z.rot",l,sep="."))
    trends.rot = get(paste("trends.rot",l,sep="."))
    ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
    # par(mfrow=c(nrow(trends.rot),1), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
    for(h in 1:nrow(trends.rot)) {
      plot(c(1:N.ts)[abs(Z.rot[,h])>=minZ], as.vector(Z.rot[abs(Z.rot[,h])>=minZ,h]),
           type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1), axes=FALSE)
      axis(1, labels=FALSE)
      axis(2, labels=FALSE)
      #   for(j in 1:N.ts) {
      #     if(Z.rot[j,h] > minZ) {text(j, -0.05, survey[j], srt=90, adj=1, cex=0.9)}
      #     if(Z.rot[j,h] < -minZ) {text(j, 0.05, survey[j], srt=90, adj=0, cex=0.9)}
      #   } # end j loop
      #   mtext(paste("Factor loadings on trend",h,sep=" "),side=3,line=.5)
      abline(h=0, lwd=1, col="gray")
      abline(h=0.2, col='gray', lty=2)
      abline(h=-0.2, col='gray', lty=2)
    } # end h loop
  }
  
  dev.off()
  
  
  
  png(filename="Trend_scale.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    index = get(paste("index",l,sep=""))
    lowerCI = get(paste("lowerCI",l,sep=""))
    upperCI = get(paste("upperCI",l,sep=""))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, index, type='l', lwd=2, pch=16))
    # title(paste("Iteration",l,sep=" "))
    abline(h=0)
    NL = get(paste("NL",l,sep=""))
    lines(yrs, rescale(NL, to=c(min(index), max(index))), type='l',col="red", lwd=2)
  }
  
  
  dev.off()
  
  
  
  
  png(filename="BackTransformed_Trend_scale.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    indexBT = get(paste("indexBT",l,sep=""))
    indexSEBT = get(paste("indexSEBT",l,sep=""))
    lowerCI = indexBT-(1.96*indexSEBT)
    upperCI = indexBT+(1.96*indexSEBT)
    df <- data.frame(yrs, indexBT, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, indexBT, type='l', lwd=2, pch=16))
    # title(paste("Iteration",l,sep=" "))
    # abline(h=0)
    N = get(paste("N",l,sep=""))
    lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT))), type='l',col="red", lwd=2)
  }
  
  
  dev.off()
  
  
  png(filename="Trend_N.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    
    index = get(paste("index",l,sep=""))
    NL = get(paste("NL",l,sep=""))
    plot(yrs, NL, type='l',col="red", axes=FALSE)
    axis(1, labels=FALSE)
    axis(2, labels=FALSE)
    lines(yrs, rescale(index, to=c(min(NL),max(NL))))
    
  }
  
  dev.off()
  
  
  
  
  
  png(filename="Trend_Standardize.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    
    index = get(paste("index",l,sep=""))
    NL = get(paste("NL",l,sep=""))
    plot(yrs, (NL-mean(NL))/sd(NL), type='l',col="red", axes=FALSE, lwd=1.5)
    axis(1, labels=FALSE)
    axis(2, labels=FALSE)
    lines(yrs, (index-mean(index))/sd(index), type='l', lwd=1.5)
    
  }
  
  dev.off()
  
  
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    index = get(paste("index",l,sep=""))
    lowerCI = get(paste("lowerCI",l,sep=""))
    upperCI = get(paste("upperCI",l,sep=""))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, index, type='l', lwd=2, pch=16))
    # title(paste("Iteration",l,sep=" "))
    abline(h=0)
    NL = get(paste("NL",l,sep=""))
    lines(yrs, rescale(NL, to=c(min(index), max(index))), type='l',col="red", lwd=2)
  }
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    indexBT = get(paste("indexBT",l,sep=""))
    indexSEBT = get(paste("indexSEBT",l,sep=""))
    lowerCI = indexBT-(1.96*indexSEBT)
    upperCI = indexBT+(1.96*indexSEBT)
    df <- data.frame(yrs, indexBT, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, indexBT, type='l', lwd=2, pch=16))
    # title(paste("Iteration",l,sep=" "))
    # abline(h=0)
    N = get(paste("N",l,sep=""))
    lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT))), type='l',col="red", lwd=2)
  }
  
  # 
  # png(filename="TwoAxesCompare.png", 
  #     type="cairo",
  #     units="mm", 
  #     width=300, 
  #     height=300, 
  #     pointsize=12, 
  #     res=600)
  # 
  # par(mfrow=c(10, 10))
  # par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  # for(l in 1:Niters){
  #   index = get(paste("index",l,sep=""))
  #   lowerCI = get(paste("lowerCI",l,sep=""))
  #   upperCI = get(paste("upperCI",l,sep=""))
  #   df <- data.frame(yrs, index, lowerCI, upperCI)
  #   with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
  #   axis(1,labels=FALSE)
  #   axis(2,labels=FALSE)
  #   x <- as.numeric(df$yrs)
  #   polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
  #   with(df, lines(x, index, type='l', lwd=2, pch=16))
  #   # title(paste("Iteration",l,sep=" "))
  #   abline(h=0)
  #   par(new=T)
  #   NL = get(paste("NL",l,sep=""))
  #   plot(yrs,NL, type='l', col='red', axes=FALSE, lwd=1.5)
  #   axis(4, labels=FALSE)
  # }
  # 
  # dev.off()
  # 
  
  
  
  
  
  
  FLoadings2 = ifelse(is.na(FLoadings)==TRUE, 0, FLoadings)
  png(filename="TwoAxesCompare_TryCorrected.png", 
      type="cairo",
      units="mm", 
      width=300, 
      height=300, 
      pointsize=12, 
      res=600)
  
  par(mfrow=c(10, 10))
  par(mar=c(0.3, 0.3, 0.1, 0.3),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:Niters){
    
    index = get(paste("index",l,sep=""))
    NL = get(paste("NL",l,sep=""))
    plot(yrs, NL, type='l',col="red", axes=FALSE, lwd=1.5)
    axis(1, labels=FALSE)
    axis(2, labels=FALSE)
    par(new=T)
    if(FLoadings2[l,10] >= 0){
      plot(yrs,index, type='l', col='black', axes=FALSE, lwd=1.5)
    }
    if(FLoadings2[l,10] < 0){
      plot(yrs,-1*index, type='l', col='black', axes=FALSE, lwd=1.5)
    }
    axis(4, labels=FALSE)
    
  }
  
  dev.off()
  
  
  
  
  ##################################################### END #################################################
  
  return(list(RMSE=RMSE, biasAll=biasAll, FitRatios=FitRatios, RMSE.BT=RMSE.BT, biasAll.BT=biasAll.BT))} # END FUNCTION


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial12")
DTrial12 = DFA_Sim4_d12(F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))), 
                        CV1=0.7,
                        CV2=0.5,
                        CV3=0.3,
                        CV4=0.5,
                        q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                        q_survey2 = rep(0.001, length=years),
                        q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                        q_survey4 = rep(0.003, length=years),
                        Iay_eq = Iay_eq_F2,
                        c=c(2.5, 40, 30, 1.5) )
save(DTrial12, file="DTrial12.RData")

# save.image("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/ATL_SN_FEMALE_SurvivalSR_LogErr_DFA_Simulation_workspace.R.RData")









################################# NEW BASE F3 #################################




# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial13")
# DTrial13 = DFA_Sim4(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
#                    CV1=0.5,
#                    CV2=0.5,
#                    CV3=0.5,
#                    CV4=0.5,
#                    q_survey1 = rep(0.001, length=years),
#                    q_survey2 = rep(0.001, length=years),
#                    q_survey3 = rep(0.001, length=years),
#                    q_survey4 = rep(0.003, length=years))
# 
# 
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial14")
# DTrial14 = DFA_Sim4(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
#                    CV1=0.3,
#                    CV2=0.5,
#                    CV3=0.7,
#                    CV4=0.5,
#                    q_survey1 = rep(0.001, length=years),
#                    q_survey2 = rep(0.001, length=years),
#                    q_survey3 = rep(0.001, length=years),
#                    q_survey4 = rep(0.003, length=years))
# 
# 
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial15")
# DTrial15 = DFA_Sim4(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
#                    CV1=0.7,
#                    CV2=0.5,
#                    CV3=0.3,
#                    CV4=0.5,
#                    q_survey1 = rep(0.001, length=years),
#                    q_survey2 = rep(0.001, length=years),
#                    q_survey3 = rep(0.001, length=years),
#                    q_survey4 = rep(0.003, length=years))






setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial16")
DTrial16 = DFA_Sim4(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                    CV1=0.5,
                    CV2=0.5,
                    CV3=0.5,
                    CV4=0.5,
                    q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F3, 
                    c=c(10000, 150, 15, 3.5))
save(DTrial16, file="DTrial16.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial17")
DTrial17 = DFA_Sim4(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                    CV1=0.3,
                    CV2=0.5,
                    CV3=0.7,
                    CV4=0.5,
                    q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F3, 
                    c=c(10000, 150, 15, 3.5))
save(DTrial17, file="DTrial17.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial18")
DTrial18 = DFA_Sim4(F_var = as.vector(c(rep(0.4, 50), rep(0, 50))), 
                    CV1=0.7,
                    CV2=0.5,
                    CV3=0.3,
                    CV4=0.5,
                    q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F3, 
                    c=c(10000, 150, 15, 3.5) )

save(DTrial18, file="DTrial18.RData")











################################# NEW BASE F4  #################################




# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial19")
# DTrial19 = DFA_Sim4(F_var = as.vector(c(rep(0,30), rep(0.4, 15), rep(0.2, 10), rep(0.05, 45))), 
#                    CV1=0.5,
#                    CV2=0.5,
#                    CV3=0.5,
#                    CV4=0.5,
#                    q_survey1 = rep(0.001, length=years),
#                    q_survey2 = rep(0.001, length=years),
#                    q_survey3 = rep(0.001, length=years),
#                    q_survey4 = rep(0.003, length=years))
# biasAll = read.csv("biasAll.csv")
# head(biasAll)
# biasAll = biasAll[,2:27]
# apply(biasAll, 1, runs.test) ## USE RUNS TEST to check for randomness/patterns in residuals
# 
# biasNormal = vector()
# for(i in 1:nrow(biasAll)){
#   b = as.numeric(biasAll[i,])
#   assign(paste('nb',i,sep=""), (b - mean(b))/sd(b) )
#   biasNormal = rbind(biasNormal, get(paste('nb',i,sep="")))
# }
# apply(biasNormal, 1, sd) ## USE standard deviation of normalized (standardized) residuals: SDNR??? 
# 
# biasAll
# N_STD = read.csv("N_STD.csv")
# N_STD = N_STD[,2:27]
# 
# 4.890716e-01/1.743064
# 
# RE = biasAll/N_STD
# MRE_yr = apply(RE, 2, median)
# MARE_yr = apply(abs(RE), 2, median)
# MRE_iter = apply(RE, 1, median)
# MARE_iter = apply(abs(RE), 1, median)
# sum(MARE_iter < 0.16)
# sum(MARE_iter < 0.35)
# sum(abs(MRE_iter) < 0.1)
# median(as.numeric(RE))


# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial20")
# DTrial20 = DFA_Sim4(F_var = as.vector(c(rep(0,30), rep(0.4, 15), rep(0.2, 10), rep(0.05, 45))), 
#                    CV1=0.3,
#                    CV2=0.5,
#                    CV3=0.7,
#                    CV4=0.5,
#                    q_survey1 = rep(0.001, length=years),
#                    q_survey2 = rep(0.001, length=years),
#                    q_survey3 = rep(0.001, length=years),
#                    q_survey4 = rep(0.003, length=years))
# 
# 
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial21")
# DTrial21 = DFA_Sim4(F_var = as.vector(c(rep(0,30), rep(0.4, 15), rep(0.2, 10), rep(0.05, 45))), 
#                    CV1=0.7,
#                    CV2=0.5,
#                    CV3=0.3,
#                    CV4=0.5,
#                    q_survey1 = rep(0.001, length=years),
#                    q_survey2 = rep(0.001, length=years),
#                    q_survey3 = rep(0.001, length=years),
#                    q_survey4 = rep(0.003, length=years))






setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial22")
DTrial22 = DFA_Sim4(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                    CV1=0.5,
                    CV2=0.5,
                    CV3=0.5,
                    CV4=0.5,
                    q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F4,
                    c=c(3000, 30, 7, 10) ,
                    reorg=F)
save(DTrial22, file="DTrial22.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial23")
DTrial23 = DFA_Sim4(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                    CV1=0.3,
                    CV2=0.5,
                    CV3=0.7,
                    CV4=0.5,
                    q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F4,
                    c=c(3000, 30, 7, 10) ,
                    reorg=F)
save(DTrial23, file="DTrial23.RData")


setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial24")
DTrial24 = DFA_Sim4(F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))), 
                    CV1=0.7,
                    CV2=0.5,
                    CV3=0.3,
                    CV4=0.5,
                    q_survey1 =   c(rep(0.00025, 49), seq(0.00025, 0.001, by=0.00005), rep(0.001, 35)),
                    q_survey2 = rep(0.001, length=years),
                    q_survey3 = c(rep(0.003, 49), seq(0.003, 0.002, by=-0.0001), rep(0.002, 40)),
                    q_survey4 = rep(0.003, length=years),
                    Iay_eq = Iay_eq_F4,
                    c=c(3000, 30, 7, 10) ,
                    reorg=F)
save(DTrial24, file="DTrial24.RData")


# save.image("N:/Documents/DFA_Simulation/DFA_Simulation/Atl_SN/ATL_SN_FEMALE_SurvivalSR_LogErr_DFA_Simulation_workspace.R.RData")




######## SAVE RESULTS ################

BigSNSimulationResults=list()
for(i in 1:24){
  BigSNSimulationResults[[paste("ATrial",i,sep="")]] = get(paste("ATrial",i,sep=""))
}
for(i in 1:24){
  BigSNSimulationResults[[paste("BTrial",i,sep="")]] = get(paste("BTrial",i,sep=""))
}
for(i in c(4:6, 10:12, 16:18, 22:24)){
  BigSNSimulationResults[[paste("CTrial",i,sep="")]] = get(paste("CTrial",i,sep=""))
}
for(i in c(4:6, 10:12, 16:18, 22:24)){
  BigSNSimulationResults[[paste("DTrial",i,sep="")]] = get(paste("DTrial",i,sep=""))
}

save(BigSNSimulationResults, file="D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/BigSNSimulationResults.RData")
load(file="D:/vspace1/DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/BigSNSimulationResults.RData")


####################### PLOT ########################################
# names(Trial24)

png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\ALL.png", 
    type="cairo",
    units="mm", 
    width=600, 
    height=250, 
    pointsize=16, 
    res=600)
############
par(mfrow=c(2,1), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))

vioplot(c( ATrial1$biasAll), c(ATrial2$biasAll), c( ATrial3$biasAll), 
        c( ATrial4$biasAll), c( ATrial5$biasAll), c( ATrial6$biasAll), 
        c( CTrial4$biasAll), c( CTrial5$biasAll), c( CTrial6$biasAll), 
        
        c( BTrial1$biasAll), c( BTrial2$biasAll), c( BTrial3$biasAll), 
        c( BTrial4$biasAll), c( BTrial5$biasAll), c( BTrial6$biasAll), 
        c( DTrial4$biasAll), c( DTrial5$biasAll), c( DTrial6$biasAll), 
        
        c( ATrial7$biasAll), c( ATrial8$biasAll), c(ATrial9$biasAll), 
        c( ATrial10$biasAll), c( ATrial11$biasAll), c( ATrial12$biasAll),
        c( CTrial10$biasAll), c( CTrial11$biasAll), c( CTrial12$biasAll),
        
        c( BTrial7$biasAll), c( BTrial8$biasAll), c( BTrial9$biasAll), 
        c( BTrial10$biasAll), c( BTrial11$biasAll), c( BTrial12$biasAll), 
        c( DTrial10$biasAll), c( DTrial11$biasAll), c( DTrial12$biasAll), 
        
        
        
        c( ATrial13$biasAll), c( ATrial14$biasAll), c( ATrial15$biasAll), 
        c( ATrial16$biasAll), c(ATrial17$biasAll), c( ATrial18$biasAll), 
        c( CTrial16$biasAll), c(CTrial17$biasAll), c( CTrial18$biasAll), 
        
        c( BTrial13$biasAll), c( BTrial14$biasAll), c( BTrial15$biasAll), 
        c( BTrial16$biasAll), c( BTrial17$biasAll), c( BTrial18$biasAll), 
        c( DTrial16$biasAll), c( DTrial17$biasAll), c( DTrial18$biasAll), 
        
        
        c( ATrial19$biasAll), c( ATrial20$biasAll), c( ATrial21$biasAll), 
        c( ATrial22$biasAll), c( ATrial23$biasAll), c( ATrial24$biasAll), 
        c( CTrial22$biasAll), c( CTrial23$biasAll), c( CTrial24$biasAll), 
        
        c( BTrial19$biasAll), c( BTrial20$biasAll), c( BTrial21$biasAll), 
        c( BTrial22$biasAll), c( BTrial23$biasAll), c( BTrial24$biasAll),
        c( DTrial22$biasAll), c( DTrial23$biasAll), c( DTrial24$biasAll),
        names=c(rep("",72)),
        
        ylim=c(-4.5,4.5), col=c(rep('darkolivegreen1',3), rep('darkolivegreen3',6), rep('darkolivegreen4',3), rep('darkolivegreen',6), rep(
          'skyblue',3), rep('deepskyblue2',6), rep('deepskyblue3',3), rep('deepskyblue4',6), rep(
            'orange',3), rep('darkorange',6), rep('darkorange2',3), rep('darkorange3',6), rep(
              'mediumorchid',3), rep('darkorchid1',6), rep('darkorchid',3), rep('darkorchid4',6)) )

title('Relative Error', line=-1.4)
abline(h=0)
vioplot(c( ATrial1$RMSE), c(ATrial2$RMSE), c( ATrial3$RMSE), 
        c( ATrial4$RMSE), c( ATrial5$RMSE), c( ATrial6$RMSE), 
        c( CTrial4$RMSE), c( CTrial5$RMSE), c( CTrial6$RMSE), 
        
        c( BTrial1$RMSE), c( BTrial2$RMSE), c( BTrial3$RMSE), 
        c( BTrial4$RMSE), c( BTrial5$RMSE), c( BTrial6$RMSE), 
        c( DTrial4$RMSE), c( DTrial5$RMSE), c( DTrial6$RMSE), 
        
        c( ATrial7$RMSE), c( ATrial8$RMSE), c(ATrial9$RMSE), 
        c( ATrial10$RMSE), c( ATrial11$RMSE), c( ATrial12$RMSE),
        c( CTrial10$RMSE), c( CTrial11$RMSE), c( CTrial12$RMSE),
        
        c( BTrial7$RMSE), c( BTrial8$RMSE), c( BTrial9$RMSE), 
        c( BTrial10$RMSE), c( BTrial11$RMSE), c( BTrial12$RMSE), 
        c( DTrial10$RMSE), c( DTrial11$RMSE), c( DTrial12$RMSE), 
        
        
        c( ATrial13$RMSE), c( ATrial14$RMSE), c( ATrial15$RMSE), 
        c( ATrial16$RMSE), c(ATrial17$RMSE), c( ATrial18$RMSE), 
        c( CTrial16$RMSE), c(CTrial17$RMSE), c( CTrial18$RMSE), 
        
        c( BTrial13$RMSE), c( BTrial14$RMSE), c( BTrial15$RMSE), 
        c( BTrial16$RMSE), c( BTrial17$RMSE), c( BTrial18$RMSE), 
        c( DTrial16$RMSE), c( DTrial17$RMSE), c( DTrial18$RMSE), 
        
        
        c( ATrial19$RMSE), c( ATrial20$RMSE), c( ATrial21$RMSE), 
        c( ATrial22$RMSE), c( ATrial23$RMSE), c( ATrial24$RMSE), 
        c( CTrial22$RMSE), c( CTrial23$RMSE), c( CTrial24$RMSE), 
        
        c( BTrial19$RMSE), c( BTrial20$RMSE), c( BTrial21$RMSE), 
        c( BTrial22$RMSE), c( BTrial23$RMSE), c( BTrial24$RMSE),
        c( DTrial22$RMSE), c( DTrial23$RMSE), c( DTrial24$RMSE),
        names=c(rep("",72)),
        ylim=c(0,2), col=c(rep('darkolivegreen1',3), rep('darkolivegreen3',6), rep('darkolivegreen4',3), rep('darkolivegreen',6), rep(
          'skyblue',3), rep('deepskyblue2',6), rep('deepskyblue3',3), rep('deepskyblue4',6), rep(
            'orange',3), rep('darkorange',6), rep('darkorange2',3), rep('darkorange3',6), rep(
              'mediumorchid',3), rep('darkorchid1',6), rep('darkorchid',3), rep('darkorchid4',6)) )

title('RMSE', line=-1.4)
##################
dev.off() 


png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\Backtransformed_ALL.png", 
    type="cairo",
    units="mm", 
    width=600, 
    height=250, 
    pointsize=16, 
    res=600)
############################
par(mfrow=c(2,1), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))

vioplot(c( ATrial1$biasAll.BT), c(ATrial2$biasAll.BT), c( ATrial3$biasAll.BT), 
        c( ATrial4$biasAll.BT), c( ATrial5$biasAll.BT), c( ATrial6$biasAll.BT), 
        c( CTrial4$biasAll.BT), c( CTrial5$biasAll.BT), c( CTrial6$biasAll.BT), 
        
        c( BTrial1$biasAll.BT), c( BTrial2$biasAll.BT), c( BTrial3$biasAll.BT), 
        c( BTrial4$biasAll.BT), c( BTrial5$biasAll.BT), c( BTrial6$biasAll.BT), 
        c( DTrial4$biasAll.BT), c( DTrial5$biasAll.BT), c( DTrial6$biasAll.BT), 
        
        c( ATrial7$biasAll.BT), c( ATrial8$biasAll.BT), c(ATrial9$biasAll.BT), 
        c( ATrial10$biasAll.BT), c( ATrial11$biasAll.BT), c( ATrial12$biasAll.BT),
        c( CTrial10$biasAll.BT), c( CTrial11$biasAll.BT), c( CTrial12$biasAll.BT),
        
        c( BTrial7$biasAll.BT), c( BTrial8$biasAll.BT), c( BTrial9$biasAll.BT), 
        c( BTrial10$biasAll.BT), c( BTrial11$biasAll.BT), c( BTrial12$biasAll.BT), 
        c( DTrial10$biasAll.BT), c( DTrial11$biasAll.BT), c( DTrial12$biasAll.BT), 
        
        
        
        c( ATrial13$biasAll.BT), c( ATrial14$biasAll.BT), c( ATrial15$biasAll.BT), 
        c( ATrial16$biasAll.BT), c(ATrial17$biasAll.BT), c( ATrial18$biasAll.BT), 
        c( CTrial16$biasAll.BT), c(CTrial17$biasAll.BT), c( CTrial18$biasAll.BT), 
        
        c( BTrial13$biasAll.BT), c( BTrial14$biasAll.BT), c( BTrial15$biasAll.BT), 
        c( BTrial16$biasAll.BT), c( BTrial17$biasAll.BT), c( BTrial18$biasAll.BT), 
        c( DTrial16$biasAll.BT), c( DTrial17$biasAll.BT), c( DTrial18$biasAll.BT), 
        
        
        c( ATrial19$biasAll.BT), c( ATrial20$biasAll.BT), c( ATrial21$biasAll.BT), 
        c( ATrial22$biasAll.BT), c( ATrial23$biasAll.BT), c( ATrial24$biasAll.BT), 
        c( CTrial22$biasAll.BT), c( CTrial23$biasAll.BT), c( CTrial24$biasAll.BT), 
        
        c( BTrial19$biasAll.BT), c( BTrial20$biasAll.BT), c( BTrial21$biasAll.BT), 
        c( BTrial22$biasAll.BT), c( BTrial23$biasAll.BT), c( BTrial24$biasAll.BT),
        c( DTrial22$biasAll.BT), c( DTrial23$biasAll.BT), c( DTrial24$biasAll.BT),
        names=c(rep("",72)),
        
        ylim=c(-4.5,4.5), col=c(rep('darkolivegreen1',3), rep('darkolivegreen3',6), rep('darkolivegreen4',3), rep('darkolivegreen',6), rep(
          'skyblue',3), rep('deepskyblue2',6), rep('deepskyblue3',3), rep('deepskyblue4',6), rep(
            'orange',3), rep('darkorange',6), rep('darkorange2',3), rep('darkorange3',6), rep(
              'mediumorchid',3), rep('darkorchid1',6), rep('darkorchid',3), rep('darkorchid4',6)) )
abline(h=0)
title('Back-transformed Relative Error', line=-1.4)
vioplot(c( ATrial1$RMSE.BT), c(ATrial2$RMSE.BT), c( ATrial3$RMSE.BT), 
        c( ATrial4$RMSE.BT), c( ATrial5$RMSE.BT), c( ATrial6$RMSE.BT), 
        c( CTrial4$RMSE.BT), c( CTrial5$RMSE.BT), c( CTrial6$RMSE.BT), 
        
        c( BTrial1$RMSE.BT), c( BTrial2$RMSE.BT), c( BTrial3$RMSE.BT), 
        c( BTrial4$RMSE.BT), c( BTrial5$RMSE.BT), c( BTrial6$RMSE.BT), 
        c( DTrial4$RMSE.BT), c( DTrial5$RMSE.BT), c( DTrial6$RMSE.BT), 
        
        c( ATrial7$RMSE.BT), c( ATrial8$RMSE.BT), c(ATrial9$RMSE.BT), 
        c( ATrial10$RMSE.BT), c( ATrial11$RMSE.BT), c( ATrial12$RMSE.BT),
        c( CTrial10$RMSE.BT), c( CTrial11$RMSE.BT), c( CTrial12$RMSE.BT),
        
        c( BTrial7$RMSE.BT), c( BTrial8$RMSE.BT), c( BTrial9$RMSE.BT), 
        c( BTrial10$RMSE.BT), c( BTrial11$RMSE.BT), c( BTrial12$RMSE.BT), 
        c( DTrial10$RMSE.BT), c( DTrial11$RMSE.BT), c( DTrial12$RMSE.BT), 
        
        
        c( ATrial13$RMSE.BT), c( ATrial14$RMSE.BT), c( ATrial15$RMSE.BT), 
        c( ATrial16$RMSE.BT), c(ATrial17$RMSE.BT), c( ATrial18$RMSE.BT), 
        c( CTrial16$RMSE.BT), c(CTrial17$RMSE.BT), c( CTrial18$RMSE.BT), 
        
        c( BTrial13$RMSE.BT), c( BTrial14$RMSE.BT), c( BTrial15$RMSE.BT), 
        c( BTrial16$RMSE.BT), c( BTrial17$RMSE.BT), c( BTrial18$RMSE.BT), 
        c( DTrial16$RMSE.BT), c( DTrial17$RMSE.BT), c( DTrial18$RMSE.BT), 
        
        
        c( ATrial19$RMSE.BT), c( ATrial20$RMSE.BT), c( ATrial21$RMSE.BT), 
        c( ATrial22$RMSE.BT), c( ATrial23$RMSE.BT), c( ATrial24$RMSE.BT), 
        c( CTrial22$RMSE.BT), c( CTrial23$RMSE.BT), c( CTrial24$RMSE.BT), 
        
        c( BTrial19$RMSE.BT), c( BTrial20$RMSE.BT), c( BTrial21$RMSE.BT), 
        c( BTrial22$RMSE.BT), c( BTrial23$RMSE.BT), c( BTrial24$RMSE.BT),
        c( DTrial22$RMSE.BT), c( DTrial23$RMSE.BT), c( DTrial24$RMSE.BT),
        names=c(rep("",72)),
        ylim=c(0,2), col=c(rep('darkolivegreen1',3), rep('darkolivegreen3',6), rep('darkolivegreen4',3), rep('darkolivegreen',6), rep(
          'skyblue',3), rep('deepskyblue2',6), rep('deepskyblue3',3), rep('deepskyblue4',6), rep(
            'orange',3), rep('darkorange',6), rep('darkorange2',3), rep('darkorange3',6), rep(
              'mediumorchid',3), rep('darkorchid1',6), rep('darkorchid',3), rep('darkorchid4',6)) )

title('Back-transformed RMSE', line=-1.4)
#####################
dev.off() 





png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\Backtransformed_3Surveys.png", 
    type="cairo",
    units="mm", 
    width=400, 
    height=170, 
    pointsize=10, 
    res=600)
################
par(mfrow=c(2,1), mar=c(0.6, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))

vioplot(c( ATrial1$biasAll.BT), c(ATrial2$biasAll.BT), c( ATrial3$biasAll.BT), 
        c( ATrial4$biasAll.BT), c( ATrial5$biasAll.BT), c( ATrial6$biasAll.BT), 
        c( CTrial4$biasAll.BT), c( CTrial5$biasAll.BT), c( CTrial6$biasAll.BT), 
        
        c( ATrial7$biasAll.BT), c( ATrial8$biasAll.BT), c(ATrial9$biasAll.BT), 
        c( ATrial10$biasAll.BT), c( ATrial11$biasAll.BT), c( ATrial12$biasAll.BT),
        c( CTrial10$biasAll.BT), c( CTrial11$biasAll.BT), c( CTrial12$biasAll.BT),
        
        c( ATrial13$biasAll.BT), c( ATrial14$biasAll.BT), c( ATrial15$biasAll.BT), 
        c( ATrial16$biasAll.BT), c(ATrial17$biasAll.BT), c( ATrial18$biasAll.BT), 
        c( CTrial16$biasAll.BT), c(CTrial17$biasAll.BT), c( CTrial18$biasAll.BT), 
        
        c( ATrial19$biasAll.BT), c( ATrial20$biasAll.BT), c( ATrial21$biasAll.BT), 
        c( ATrial22$biasAll.BT), c( ATrial23$biasAll.BT), c( ATrial24$biasAll.BT), 
        c( CTrial22$biasAll.BT), c( CTrial23$biasAll.BT), c( CTrial24$biasAll.BT), 
        
        names=c(rep("",36)),cex=1,pchMed=19, 
        
        ylim=c(-4.25,4.5), col=c(rep('darkolivegreen1',3), rep('darkolivegreen3',3), rep('darkolivegreen',3), rep(
          'skyblue',3), rep('deepskyblue3',3), rep('deepskyblue4',3), rep(
            'orange',3), rep('darkorange2',3), rep('darkorange3',3), rep(
              'mediumorchid',3), rep('darkorchid',3), rep('darkorchid4',3)) )
abline(h=0)
mtext('Back-transformed Relative Error', 2, line=1, cex=1.2)
mtext('Atlantic sharpnose: 3 surveys', 3, line=-1.4, cex=1.2)

par( mar=c(1.3, 2.1, 0.1, 0.3))
vioplot(c( ATrial1$RMSE.BT), c(ATrial2$RMSE.BT), c( ATrial3$RMSE.BT), 
        c( ATrial4$RMSE.BT), c( ATrial5$RMSE.BT), c( ATrial6$RMSE.BT), 
        c( CTrial4$RMSE.BT), c( CTrial5$RMSE.BT), c( CTrial6$RMSE.BT), 
        
        c( ATrial7$RMSE.BT), c( ATrial8$RMSE.BT), c(ATrial9$RMSE.BT), 
        c( ATrial10$RMSE.BT), c( ATrial11$RMSE.BT), c( ATrial12$RMSE.BT),
        c( CTrial10$RMSE.BT), c( CTrial11$RMSE.BT), c( CTrial12$RMSE.BT),
        
        c( ATrial13$RMSE.BT), c( ATrial14$RMSE.BT), c( ATrial15$RMSE.BT), 
        c( ATrial16$RMSE.BT), c(ATrial17$RMSE.BT), c( ATrial18$RMSE.BT), 
        c( CTrial16$RMSE.BT), c(CTrial17$RMSE.BT), c( CTrial18$RMSE.BT), 
        
        c( ATrial19$RMSE.BT), c( ATrial20$RMSE.BT), c( ATrial21$RMSE.BT), 
        c( ATrial22$RMSE.BT), c( ATrial23$RMSE.BT), c( ATrial24$RMSE.BT), 
        c( CTrial22$RMSE.BT), c( CTrial23$RMSE.BT), c( CTrial24$RMSE.BT), 
        
        names=c("A1","A2","A3","A4","A5","A6",
                "C4","C5","C6",
                "A7","A8","A9","A10","A11","A12",
                "C10","C11","C12",
                "A13","A14","A15","A16","A17","A18",
                "C16","C17","C18",
                "A19","A20","A21","A22","A23","A24",
                "C22","C23","C24"), cex=1,pchMed=19,
        ylim=c(0,2),  col=c(rep('darkolivegreen1',3), rep('darkolivegreen3',3), rep('darkolivegreen',3), rep(
          'skyblue',3), rep('deepskyblue3',3), rep('deepskyblue4',3), rep(
            'orange',3), rep('darkorange2',3), rep('darkorange3',3), rep(
              'mediumorchid',3), rep('darkorchid',3), rep('darkorchid4',3)) )

mtext('Back-transformed RMSE', 2, line=1, cex=1.2)
######################
dev.off() 


png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\Backtransformed_4Surveys.png", 
    type="cairo",
    units="mm", 
    width=400, 
    height=170, 
    pointsize=10, 
    res=600)
########################
par(mfrow=c(2,1), mar=c(0.6, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))

vioplot(c( BTrial1$biasAll.BT), c( BTrial2$biasAll.BT), c( BTrial3$biasAll.BT), 
        c( BTrial4$biasAll.BT), c( BTrial5$biasAll.BT), c( BTrial6$biasAll.BT), 
        c( DTrial4$biasAll.BT), c( DTrial5$biasAll.BT), c( DTrial6$biasAll.BT), 
        
        c( BTrial7$biasAll.BT), c( BTrial8$biasAll.BT), c( BTrial9$biasAll.BT), 
        c( BTrial10$biasAll.BT), c( BTrial11$biasAll.BT), c( BTrial12$biasAll.BT), 
        c( DTrial10$biasAll.BT), c( DTrial11$biasAll.BT), c( DTrial12$biasAll.BT), 
        
        c( BTrial13$biasAll.BT), c( BTrial14$biasAll.BT), c( BTrial15$biasAll.BT), 
        c( BTrial16$biasAll.BT), c( BTrial17$biasAll.BT), c( BTrial18$biasAll.BT), 
        c( DTrial16$biasAll.BT), c( DTrial17$biasAll.BT), c( DTrial18$biasAll.BT), 
        
        c( BTrial19$biasAll.BT), c( BTrial20$biasAll.BT), c( BTrial21$biasAll.BT), 
        c( BTrial22$biasAll.BT), c( BTrial23$biasAll.BT), c( BTrial24$biasAll.BT),
        c( DTrial22$biasAll.BT), c( DTrial23$biasAll.BT), c( DTrial24$biasAll.BT),
        
        names=c(rep("",36)),cex=1,pchMed=19, 
        
        ylim=c(-4.25,4.5), col=c(rep('darkolivegreen1',3), rep('darkolivegreen3',3), rep('darkolivegreen',3), rep(
          'skyblue',3), rep('deepskyblue3',3), rep('deepskyblue4',3), rep(
            'orange',3), rep('darkorange2',3), rep('darkorange3',3), rep(
              'mediumorchid',3), rep('darkorchid',3), rep('darkorchid4',3)) )
abline(h=0)
mtext('Back-transformed Relative Error', 2, line=1, cex=1.2)
mtext('Atlantic sharpnose: 4 surveys', 3, line=-1.4, cex=1.2)

par( mar=c(1.3, 2.1, 0.1, 0.3))
vioplot(c( BTrial1$RMSE.BT), c( BTrial2$RMSE.BT), c( BTrial3$RMSE.BT), 
        c( BTrial4$RMSE.BT), c( BTrial5$RMSE.BT), c( BTrial6$RMSE.BT), 
        c( DTrial4$RMSE.BT), c( DTrial5$RMSE.BT), c( DTrial6$RMSE.BT), 
        
        c( BTrial7$RMSE.BT), c( BTrial8$RMSE.BT), c( BTrial9$RMSE.BT), 
        c( BTrial10$RMSE.BT), c( BTrial11$RMSE.BT), c( BTrial12$RMSE.BT), 
        c( DTrial10$RMSE.BT), c( DTrial11$RMSE.BT), c( DTrial12$RMSE.BT), 
        
        c( BTrial13$RMSE.BT), c( BTrial14$RMSE.BT), c( BTrial15$RMSE.BT), 
        c( BTrial16$RMSE.BT), c( BTrial17$RMSE.BT), c( BTrial18$RMSE.BT), 
        c( DTrial16$RMSE.BT), c( DTrial17$RMSE.BT), c( DTrial18$RMSE.BT), 
        
        c( BTrial19$RMSE.BT), c( BTrial20$RMSE.BT), c( BTrial21$RMSE.BT), 
        c( BTrial22$RMSE.BT), c( BTrial23$RMSE.BT), c( BTrial24$RMSE.BT),
        c( DTrial22$RMSE.BT), c( DTrial23$RMSE.BT), c( DTrial24$RMSE.BT),
        
        names=c("B1","B2","B3","B4","B5","B6",
                "D4","D5","D6",
                "B7","B8","B9","B10","B11","B12",
                "D10","D11","D12",
                "B13","B14","B15","B16","B17","B18",
                "D16","D17","D18",
                "B19","B20","B21","B22","B23","B24",
                "D22","D23","D24"), cex=1,pchMed=19,
        ylim=c(0,2),  col=c(rep('darkolivegreen1',3), rep('darkolivegreen3',3), rep('darkolivegreen',3), rep(
          'skyblue',3), rep('deepskyblue3',3), rep('deepskyblue4',3), rep(
            'orange',3), rep('darkorange2',3), rep('darkorange3',3), rep(
              'mediumorchid',3), rep('darkorchid',3), rep('darkorchid4',3)) )

mtext('Back-transformed RMSE', 2, line=1, cex=1.2)
########################
dev.off()





png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\Backtransformed_3Surveys_grid.png", 
    type="cairo",
    units="mm", 
    width=400, 
    height=170, 
    pointsize=10, 
    res=600)
################
par(mfrow=c(2,1), mar=c(0.6, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0), oma=c(2, 0.5, 0.5, 0.5))

vioplot(c( ATrial1$biasAll.BT), c(ATrial2$biasAll.BT), c( ATrial3$biasAll.BT), 
        c( ATrial4$biasAll.BT), c( ATrial5$biasAll.BT), c( ATrial6$biasAll.BT), 
        c( CTrial4$biasAll.BT), c( CTrial5$biasAll.BT), c( CTrial6$biasAll.BT), 
        
        c( ATrial7$biasAll.BT), c( ATrial8$biasAll.BT), c(ATrial9$biasAll.BT), 
        c( ATrial10$biasAll.BT), c( ATrial11$biasAll.BT), c( ATrial12$biasAll.BT),
        c( CTrial10$biasAll.BT), c( CTrial11$biasAll.BT), c( CTrial12$biasAll.BT),
        
        c( ATrial13$biasAll.BT), c( ATrial14$biasAll.BT), c( ATrial15$biasAll.BT), 
        c( ATrial16$biasAll.BT), c(ATrial17$biasAll.BT), c( ATrial18$biasAll.BT), 
        c( CTrial16$biasAll.BT), c(CTrial17$biasAll.BT), c( CTrial18$biasAll.BT), 
        
        c( ATrial19$biasAll.BT), c( ATrial20$biasAll.BT), c( ATrial21$biasAll.BT), 
        c( ATrial22$biasAll.BT), c( ATrial23$biasAll.BT), c( ATrial24$biasAll.BT), 
        c( CTrial22$biasAll.BT), c( CTrial23$biasAll.BT), c( CTrial24$biasAll.BT), 
        
        names=c(rep("",36)),cex=1,pchMed=19,  panel.first=grid(NA,NULL, lty=1),
        
        ylim=c(-4.25,4.5), col=c(rep('darkolivegreen1',3), rep('darkolivegreen3',3), rep('darkolivegreen',3), rep(
          'skyblue',3), rep('deepskyblue3',3), rep('deepskyblue4',3), rep(
            'orange',3), rep('darkorange2',3), rep('darkorange3',3), rep(
              'mediumorchid',3), rep('darkorchid',3), rep('darkorchid4',3)) )
abline(h=0)
mtext('Back-transformed Relative Error', 2, line=1, cex=1.2)
mtext('Atlantic sharpnose: 3 surveys', 3, line=-1.4, cex=1.2)

par( mar=c(1.3, 2.1, 0.1, 0.3))
vioplot(c( ATrial1$RMSE.BT), c(ATrial2$RMSE.BT), c( ATrial3$RMSE.BT), 
        c( ATrial4$RMSE.BT), c( ATrial5$RMSE.BT), c( ATrial6$RMSE.BT), 
        c( CTrial4$RMSE.BT), c( CTrial5$RMSE.BT), c( CTrial6$RMSE.BT), 
        
        c( ATrial7$RMSE.BT), c( ATrial8$RMSE.BT), c(ATrial9$RMSE.BT), 
        c( ATrial10$RMSE.BT), c( ATrial11$RMSE.BT), c( ATrial12$RMSE.BT),
        c( CTrial10$RMSE.BT), c( CTrial11$RMSE.BT), c( CTrial12$RMSE.BT),
        
        c( ATrial13$RMSE.BT), c( ATrial14$RMSE.BT), c( ATrial15$RMSE.BT), 
        c( ATrial16$RMSE.BT), c(ATrial17$RMSE.BT), c( ATrial18$RMSE.BT), 
        c( CTrial16$RMSE.BT), c(CTrial17$RMSE.BT), c( CTrial18$RMSE.BT), 
        
        c( ATrial19$RMSE.BT), c( ATrial20$RMSE.BT), c( ATrial21$RMSE.BT), 
        c( ATrial22$RMSE.BT), c( ATrial23$RMSE.BT), c( ATrial24$RMSE.BT), 
        c( CTrial22$RMSE.BT), c( CTrial23$RMSE.BT), c( CTrial24$RMSE.BT), 
        
        # names=c("A1","A2","A3","A4","A5","A6",
        #         "C4","C5","C6",
        #         "A7","A8","A9","A10","A11","A12",
        #         "C10","C11","C12",
        #         "A13","A14","A15","A16","A17","A18",
        #         "C16","C17","C18",
        #         "A19","A20","A21","A22","A23","A24",
        #         "C22","C23","C24"),
        names=FALSE,
        cex=1,pchMed=19,
        ylim=c(0,2),  
        panel.first=grid(NA,NULL, lty=1),
        col=c(rep('darkolivegreen1',3), rep('darkolivegreen3',3), rep('darkolivegreen',3), rep(
          'skyblue',3), rep('deepskyblue3',3), rep('deepskyblue4',3), rep(
            'orange',3), rep('darkorange2',3), rep('darkorange3',3), rep(
              'mediumorchid',3), rep('darkorchid',3), rep('darkorchid4',3)) )

#names=FALSE , lwd=0.8)
axis(1, labels=FALSE, tck=0)
text(seq(1, 36, by=1), par("usr")[3] - 0.03, 
     labels = c("CF_c_1","CF_c_2","CF_c_3","CF_k_1","CF_k_2","CF_k_3",
                "CF_g_1","CF_g_2","CF_g_3",
                
                "IF_c_1","IF_c_2","IF_c_3","IF_k_1","IF_k_2","IF_k_3",
                "IF_g_1","IF_g_2","IF_g_3",
                
                "DF_c_1","DF_c_2","DF_c_3","DF_k_1","DF_k_2","DF_k_3",
                "DF_g_1","DF_g_2","DF_g_3",
                
                "UF_c_1","UF_c_2","UF_c_3","UF_k_1","UF_k_2","UF_k_3",
                "UF_g_1","UF_g_2","UF_g_3"), 
     srt = 75, xpd = NA, adj=1, cex=0.9)

mtext('Back-transformed RMSE', 2, line=1, cex=1.2)
######################
dev.off() 


png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\Backtransformed_4Surveys_grid.png", 
    type="cairo",
    units="mm", 
    width=400, 
    height=170, 
    pointsize=10, 
    res=600)
########################
par(mfrow=c(2,1), mar=c(0.6, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0), oma=c(2, 0.5, 0.5, 0.5))

vioplot(c( BTrial1$biasAll.BT), c( BTrial2$biasAll.BT), c( BTrial3$biasAll.BT), 
        c( BTrial4$biasAll.BT), c( BTrial5$biasAll.BT), c( BTrial6$biasAll.BT), 
        c( DTrial4$biasAll.BT), c( DTrial5$biasAll.BT), c( DTrial6$biasAll.BT), 
        
        c( BTrial7$biasAll.BT), c( BTrial8$biasAll.BT), c( BTrial9$biasAll.BT), 
        c( BTrial10$biasAll.BT), c( BTrial11$biasAll.BT), c( BTrial12$biasAll.BT), 
        c( DTrial10$biasAll.BT), c( DTrial11$biasAll.BT), c( DTrial12$biasAll.BT), 
        
        c( BTrial13$biasAll.BT), c( BTrial14$biasAll.BT), c( BTrial15$biasAll.BT), 
        c( BTrial16$biasAll.BT), c( BTrial17$biasAll.BT), c( BTrial18$biasAll.BT), 
        c( DTrial16$biasAll.BT), c( DTrial17$biasAll.BT), c( DTrial18$biasAll.BT), 
        
        c( BTrial19$biasAll.BT), c( BTrial20$biasAll.BT), c( BTrial21$biasAll.BT), 
        c( BTrial22$biasAll.BT), c( BTrial23$biasAll.BT), c( BTrial24$biasAll.BT),
        c( DTrial22$biasAll.BT), c( DTrial23$biasAll.BT), c( DTrial24$biasAll.BT),
        
        names=c(rep("",36)),cex=1,pchMed=19,  panel.first=grid(NA,NULL, lty=1),
        
        ylim=c(-4.25,4.5), col=c(rep('darkolivegreen1',3), rep('darkolivegreen3',3), rep('darkolivegreen',3), rep(
          'skyblue',3), rep('deepskyblue3',3), rep('deepskyblue4',3), rep(
            'orange',3), rep('darkorange2',3), rep('darkorange3',3), rep(
              'mediumorchid',3), rep('darkorchid',3), rep('darkorchid4',3)) )
abline(h=0)
mtext('Back-transformed Relative Error', 2, line=1, cex=1.2)
mtext('Atlantic sharpnose: 4 surveys', 3, line=-1.4, cex=1.2)

par( mar=c(1.3, 2.1, 0.1, 0.3))
vioplot(c( BTrial1$RMSE.BT), c( BTrial2$RMSE.BT), c( BTrial3$RMSE.BT), 
        c( BTrial4$RMSE.BT), c( BTrial5$RMSE.BT), c( BTrial6$RMSE.BT), 
        c( DTrial4$RMSE.BT), c( DTrial5$RMSE.BT), c( DTrial6$RMSE.BT), 
        
        c( BTrial7$RMSE.BT), c( BTrial8$RMSE.BT), c( BTrial9$RMSE.BT), 
        c( BTrial10$RMSE.BT), c( BTrial11$RMSE.BT), c( BTrial12$RMSE.BT), 
        c( DTrial10$RMSE.BT), c( DTrial11$RMSE.BT), c( DTrial12$RMSE.BT), 
        
        c( BTrial13$RMSE.BT), c( BTrial14$RMSE.BT), c( BTrial15$RMSE.BT), 
        c( BTrial16$RMSE.BT), c( BTrial17$RMSE.BT), c( BTrial18$RMSE.BT), 
        c( DTrial16$RMSE.BT), c( DTrial17$RMSE.BT), c( DTrial18$RMSE.BT), 
        
        c( BTrial19$RMSE.BT), c( BTrial20$RMSE.BT), c( BTrial21$RMSE.BT), 
        c( BTrial22$RMSE.BT), c( BTrial23$RMSE.BT), c( BTrial24$RMSE.BT),
        c( DTrial22$RMSE.BT), c( DTrial23$RMSE.BT), c( DTrial24$RMSE.BT),
        
        # names=c("B1","B2","B3","B4","B5","B6",
        #         "D4","D5","D6",
        #         "B7","B8","B9","B10","B11","B12",
        #         "D10","D11","D12",
        #         "B13","B14","B15","B16","B17","B18",
        #         "D16","D17","D18",
        #         "B19","B20","B21","B22","B23","B24",
        #         "D22","D23","D24"), 
        names=FALSE,
        cex=1,pchMed=19,
        ylim=c(0,2),   panel.first=grid(NA,NULL, lty=1),
        col=c(rep('darkolivegreen1',3), rep('darkolivegreen3',3), rep('darkolivegreen',3), rep(
          'skyblue',3), rep('deepskyblue3',3), rep('deepskyblue4',3), rep(
            'orange',3), rep('darkorange2',3), rep('darkorange3',3), rep(
              'mediumorchid',3), rep('darkorchid',3), rep('darkorchid4',3)) )
axis(1, labels=FALSE, tck=0)
text(seq(1, 36, by=1), par("usr")[3] - 0.03, 
     labels = c("CF_c_1","CF_c_2","CF_c_3","CF_k_1","CF_k_2","CF_k_3",
                "CF_g_1","CF_g_2","CF_g_3",
                
                "IF_c_1","IF_c_2","IF_c_3","IF_k_1","IF_k_2","IF_k_3",
                "IF_g_1","IF_g_2","IF_g_3",
                
                "DF_c_1","DF_c_2","DF_c_3","DF_k_1","DF_k_2","DF_k_3",
                "DF_g_1","DF_g_2","DF_g_3",
                
                "UF_c_1","UF_c_2","UF_c_3","UF_k_1","UF_k_2","UF_k_3",
                "UF_g_1","UF_g_2","UF_g_3"), 
     srt = 75, xpd = NA, adj=1, cex=0.9)

mtext('Back-transformed RMSE', 2, line=1, cex=1.2)
########################
dev.off()
###############################################################






legend('top', c('Constant q, 3 surveys', 'Knife-edge change in q, 3 surveys', 'Gradual change in q, 3 surveys',
                'Constant q, 4 surveys', 'Knife-edge change in q, 4 surveys', 'Gradual change in q, 4 surveys'), 
       fill=c("deepskyblue",'lightskyblue','cyan','darkolivegreen3','darkolivegreen1','green'),ncol=3)





#### BY F PATTERN ##########
png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\ByFPattern.png", 
    type="cairo",
    units="mm", 
    width=500, 
    height=167, 
    pointsize=16, 
    res=600)
##################
par(mfrow=c(1,3), mar=c(1.6, 1.6, 1.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))

vioplot(c(ATrial1$biasAll, ATrial2$biasAll, ATrial3$biasAll, ATrial4$biasAll, ATrial5$biasAll, ATrial6$biasAll, 
          CTrial4$biasAll, CTrial5$biasAll, CTrial6$biasAll, 
          BTrial1$biasAll, BTrial2$biasAll, BTrial3$biasAll, BTrial4$biasAll, BTrial5$biasAll, BTrial6$biasAll,
          DTrial4$biasAll, DTrial5$biasAll, DTrial6$biasAll),
        c(ATrial7$biasAll, ATrial8$biasAll, ATrial9$biasAll, ATrial10$biasAll, ATrial11$biasAll, ATrial12$biasAll, 
          CTrial10$biasAll, CTrial11$biasAll, CTrial12$biasAll,  
          BTrial7$biasAll, BTrial8$biasAll, BTrial9$biasAll, BTrial10$biasAll, BTrial11$biasAll, BTrial12$biasAll, 
          DTrial10$biasAll, DTrial11$biasAll, DTrial12$biasAll),
        c(ATrial13$biasAll, ATrial14$biasAll, ATrial15$biasAll, ATrial16$biasAll, ATrial17$biasAll, ATrial18$biasAll, 
          CTrial16$biasAll, CTrial17$biasAll, CTrial18$biasAll, 
          BTrial13$biasAll, BTrial14$biasAll, BTrial15$biasAll, BTrial16$biasAll, BTrial17$biasAll, BTrial18$biasAll, 
          DTrial16$biasAll, DTrial17$biasAll, DTrial18$biasAll),
        c(ATrial19$biasAll, ATrial20$biasAll, ATrial21$biasAll, ATrial22$biasAll, ATrial23$biasAll, ATrial24$biasAll, 
          CTrial22$biasAll, CTrial23$biasAll, CTrial24$biasAll,  
          BTrial19$biasAll, BTrial20$biasAll, BTrial21$biasAll, BTrial22$biasAll, BTrial23$biasAll, BTrial24$biasAll, 
          DTrial22$biasAll, DTrial23$biasAll, DTrial24$biasAll),
        ylim=c(-4.5,4.5), col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
title('Relative Error', line=-1.4)
# abline(h=0)

vioplot(c(ATrial1$RMSE, ATrial2$RMSE, ATrial3$RMSE, ATrial4$RMSE, ATrial5$RMSE, ATrial6$RMSE, 
          CTrial4$RMSE, CTrial5$RMSE, CTrial6$RMSE, 
          BTrial1$RMSE, BTrial2$RMSE, BTrial3$RMSE, BTrial4$RMSE, BTrial5$RMSE, BTrial6$RMSE,
          DTrial4$RMSE, DTrial5$RMSE, DTrial6$RMSE), 
        c(ATrial7$RMSE, ATrial8$RMSE, ATrial9$RMSE, ATrial10$RMSE, ATrial11$RMSE, ATrial12$RMSE, 
          CTrial10$RMSE, CTrial11$RMSE, CTrial12$RMSE,  
          BTrial7$RMSE, BTrial8$RMSE, BTrial9$RMSE, BTrial10$RMSE, BTrial11$RMSE, BTrial12$RMSE, 
          DTrial10$RMSE, DTrial11$RMSE, DTrial12$RMSE),
        c(ATrial13$RMSE, ATrial14$RMSE, ATrial15$RMSE, ATrial16$RMSE, ATrial17$RMSE, ATrial18$RMSE, 
          CTrial16$RMSE, CTrial17$RMSE, CTrial18$RMSE, 
          BTrial13$RMSE, BTrial14$RMSE, BTrial15$RMSE, BTrial16$RMSE, BTrial17$RMSE, BTrial18$RMSE, 
          DTrial16$RMSE, DTrial17$RMSE, DTrial18$RMSE),
        c(ATrial19$RMSE, ATrial20$RMSE, ATrial21$RMSE, ATrial22$RMSE, ATrial23$RMSE, ATrial24$RMSE, 
          CTrial22$RMSE, CTrial23$RMSE, CTrial24$RMSE,  
          BTrial19$RMSE, BTrial20$RMSE, BTrial21$RMSE, BTrial22$RMSE, BTrial23$RMSE, BTrial24$RMSE, 
          DTrial22$RMSE, DTrial23$RMSE, DTrial24$RMSE),
        ylim=c(0,2), col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
title('RMSE', line=-1.4)



vioplot(c(apply(ATrial1$FitRatios, 1, FUN=mean), apply(ATrial2$FitRatios, 1, FUN=mean), 
          apply(ATrial3$FitRatios, 1, FUN=mean), apply(ATrial4$FitRatios, 1, FUN=mean),
          apply(ATrial5$FitRatios, 1, FUN=mean), apply(ATrial6$FitRatios, 1, FUN=mean), 
          apply(CTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean),
          apply(CTrial6$FitRatios, 1, FUN=mean),
          apply(BTrial1$FitRatios, 1, FUN=mean), apply(BTrial2$FitRatios, 1, FUN=mean), 
          apply(BTrial3$FitRatios, 1, FUN=mean), apply(BTrial4$FitRatios, 1, FUN=mean),
          apply(BTrial5$FitRatios, 1, FUN=mean), apply(BTrial6$FitRatios, 1, FUN=mean), 
          apply(DTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean),
          apply(DTrial6$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial7$FitRatios, 1, FUN=mean), apply(ATrial8$FitRatios, 1, FUN=mean), 
          apply(ATrial9$FitRatios, 1, FUN=mean), apply(ATrial10$FitRatios, 1, FUN=mean),
          apply(ATrial11$FitRatios, 1, FUN=mean), apply(ATrial12$FitRatios, 1, FUN=mean), 
          apply(CTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean),
          apply(CTrial12$FitRatios, 1, FUN=mean),
          apply(BTrial7$FitRatios, 1, FUN=mean), apply(BTrial8$FitRatios, 1, FUN=mean), 
          apply(BTrial9$FitRatios, 1, FUN=mean), apply(BTrial10$FitRatios, 1, FUN=mean),
          apply(BTrial11$FitRatios, 1, FUN=mean), apply(BTrial12$FitRatios, 1, FUN=mean), 
          apply(DTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean),
          apply(DTrial12$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial13$FitRatios, 1, FUN=mean), apply(ATrial14$FitRatios, 1, FUN=mean), 
          apply(ATrial15$FitRatios, 1, FUN=mean), apply(ATrial16$FitRatios, 1, FUN=mean),
          apply(ATrial17$FitRatios, 1, FUN=mean), apply(ATrial18$FitRatios, 1, FUN=mean), 
          apply(CTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean),
          apply(CTrial18$FitRatios, 1, FUN=mean),
          apply(BTrial13$FitRatios, 1, FUN=mean), apply(BTrial14$FitRatios, 1, FUN=mean), 
          apply(BTrial15$FitRatios, 1, FUN=mean), apply(BTrial16$FitRatios, 1, FUN=mean),
          apply(BTrial17$FitRatios, 1, FUN=mean), apply(BTrial18$FitRatios, 1, FUN=mean), 
          apply(DTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean),
          apply(DTrial18$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial19$FitRatios, 1, FUN=mean), apply(ATrial20$FitRatios, 1, FUN=mean), 
          apply(ATrial21$FitRatios, 1, FUN=mean), apply(ATrial22$FitRatios, 1, FUN=mean),
          apply(ATrial23$FitRatios, 1, FUN=mean), apply(ATrial24$FitRatios, 1, FUN=mean), 
          apply(CTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean),
          apply(CTrial24$FitRatios, 1, FUN=mean),
          apply(BTrial19$FitRatios, 1, FUN=mean), apply(BTrial20$FitRatios, 1, FUN=mean), 
          apply(BTrial21$FitRatios, 1, FUN=mean), apply(BTrial22$FitRatios, 1, FUN=mean),
          apply(BTrial23$FitRatios, 1, FUN=mean), apply(BTrial24$FitRatios, 1, FUN=mean), 
          apply(DTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean),
          apply(DTrial24$FitRatios, 1, FUN=mean)),
        
        col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F')) #ylim=c(0,2), 
title('Fit Ratio', line=-1.4)
################
dev.off()


#### BY F PATTERN ##########
png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\ByFPattern_LgFont.png", 
    type="cairo",
    units="mm", 
    width=250, 
    height=175, 
    pointsize=16, 
    res=600)
#############
par(mfrow=c(2,1), mar=c(0.6, 2.1, 0.3, 0.1),tcl = -0.1, mgp = c(0.8, 0.2, 0), cex=1)

vioplot(c(ATrial1$biasAll, ATrial2$biasAll, ATrial3$biasAll, ATrial4$biasAll, ATrial5$biasAll, ATrial6$biasAll, 
          CTrial4$biasAll, CTrial5$biasAll, CTrial6$biasAll, 
          BTrial1$biasAll, BTrial2$biasAll, BTrial3$biasAll, BTrial4$biasAll, BTrial5$biasAll, BTrial6$biasAll,
          DTrial4$biasAll, DTrial5$biasAll, DTrial6$biasAll),
        c(ATrial7$biasAll, ATrial8$biasAll, ATrial9$biasAll, ATrial10$biasAll, ATrial11$biasAll, ATrial12$biasAll, 
          CTrial10$biasAll, CTrial11$biasAll, CTrial12$biasAll,  
          BTrial7$biasAll, BTrial8$biasAll, BTrial9$biasAll, BTrial10$biasAll, BTrial11$biasAll, BTrial12$biasAll, 
          DTrial10$biasAll, DTrial11$biasAll, DTrial12$biasAll),
        c(ATrial13$biasAll, ATrial14$biasAll, ATrial15$biasAll, ATrial16$biasAll, ATrial17$biasAll, ATrial18$biasAll, 
          CTrial16$biasAll, CTrial17$biasAll, CTrial18$biasAll, 
          BTrial13$biasAll, BTrial14$biasAll, BTrial15$biasAll, BTrial16$biasAll, BTrial17$biasAll, BTrial18$biasAll, 
          DTrial16$biasAll, DTrial17$biasAll, DTrial18$biasAll),
        c(ATrial19$biasAll, ATrial20$biasAll, ATrial21$biasAll, ATrial22$biasAll, ATrial23$biasAll, ATrial24$biasAll, 
          CTrial22$biasAll, CTrial23$biasAll, CTrial24$biasAll,  
          BTrial19$biasAll, BTrial20$biasAll, BTrial21$biasAll, BTrial22$biasAll, BTrial23$biasAll, BTrial24$biasAll, 
          DTrial22$biasAll, DTrial23$biasAll, DTrial24$biasAll),
        ylim=c(-4.5,4.5), col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
        names=c("","","",""), ylab="")
abline(h=0)
mtext('Relative Error',side=2,line=1.1,cex=1.1)
title('Performance by F Pattern', line=-1.2, font.main=1, cex=1)
# abline(h=0)

par(mar=c(2.1, 2.1, 0.3, 0.1))
vioplot(c(ATrial1$RMSE, ATrial2$RMSE, ATrial3$RMSE, ATrial4$RMSE, ATrial5$RMSE, ATrial6$RMSE, 
          CTrial4$RMSE, CTrial5$RMSE, CTrial6$RMSE, 
          BTrial1$RMSE, BTrial2$RMSE, BTrial3$RMSE, BTrial4$RMSE, BTrial5$RMSE, BTrial6$RMSE,
          DTrial4$RMSE, DTrial5$RMSE, DTrial6$RMSE), 
        c(ATrial7$RMSE, ATrial8$RMSE, ATrial9$RMSE, ATrial10$RMSE, ATrial11$RMSE, ATrial12$RMSE, 
          CTrial10$RMSE, CTrial11$RMSE, CTrial12$RMSE,  
          BTrial7$RMSE, BTrial8$RMSE, BTrial9$RMSE, BTrial10$RMSE, BTrial11$RMSE, BTrial12$RMSE, 
          DTrial10$RMSE, DTrial11$RMSE, DTrial12$RMSE),
        c(ATrial13$RMSE, ATrial14$RMSE, ATrial15$RMSE, ATrial16$RMSE, ATrial17$RMSE, ATrial18$RMSE, 
          CTrial16$RMSE, CTrial17$RMSE, CTrial18$RMSE, 
          BTrial13$RMSE, BTrial14$RMSE, BTrial15$RMSE, BTrial16$RMSE, BTrial17$RMSE, BTrial18$RMSE, 
          DTrial16$RMSE, DTrial17$RMSE, DTrial18$RMSE),
        c(ATrial19$RMSE, ATrial20$RMSE, ATrial21$RMSE, ATrial22$RMSE, ATrial23$RMSE, ATrial24$RMSE, 
          CTrial22$RMSE, CTrial23$RMSE, CTrial24$RMSE,  
          BTrial19$RMSE, BTrial20$RMSE, BTrial21$RMSE, BTrial22$RMSE, BTrial23$RMSE, BTrial24$RMSE, 
          DTrial22$RMSE, DTrial23$RMSE, DTrial24$RMSE),
        ylim=c(0,2), ylab="", col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
mtext('RMSE',side=2,line=1.1,cex=1.1)
# title('RMSE', line=-1.4)


# 
# vioplot(c(apply(ATrial1$FitRatios, 1, FUN=mean), apply(ATrial2$FitRatios, 1, FUN=mean), 
#           apply(ATrial3$FitRatios, 1, FUN=mean), apply(ATrial4$FitRatios, 1, FUN=mean),
#           apply(ATrial5$FitRatios, 1, FUN=mean), apply(ATrial6$FitRatios, 1, FUN=mean), 
#           apply(CTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean),
#           apply(CTrial6$FitRatios, 1, FUN=mean),
#           apply(BTrial1$FitRatios, 1, FUN=mean), apply(BTrial2$FitRatios, 1, FUN=mean), 
#           apply(BTrial3$FitRatios, 1, FUN=mean), apply(BTrial4$FitRatios, 1, FUN=mean),
#           apply(BTrial5$FitRatios, 1, FUN=mean), apply(BTrial6$FitRatios, 1, FUN=mean), 
#           apply(DTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean),
#           apply(DTrial6$FitRatios, 1, FUN=mean)),
#         
#         c(apply(ATrial7$FitRatios, 1, FUN=mean), apply(ATrial8$FitRatios, 1, FUN=mean), 
#           apply(ATrial9$FitRatios, 1, FUN=mean), apply(ATrial10$FitRatios, 1, FUN=mean),
#           apply(ATrial11$FitRatios, 1, FUN=mean), apply(ATrial12$FitRatios, 1, FUN=mean), 
#           apply(CTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean),
#           apply(CTrial12$FitRatios, 1, FUN=mean),
#           apply(BTrial7$FitRatios, 1, FUN=mean), apply(BTrial8$FitRatios, 1, FUN=mean), 
#           apply(BTrial9$FitRatios, 1, FUN=mean), apply(BTrial10$FitRatios, 1, FUN=mean),
#           apply(BTrial11$FitRatios, 1, FUN=mean), apply(BTrial12$FitRatios, 1, FUN=mean), 
#           apply(DTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean),
#           apply(DTrial12$FitRatios, 1, FUN=mean)),
#         
#         c(apply(ATrial13$FitRatios, 1, FUN=mean), apply(ATrial14$FitRatios, 1, FUN=mean), 
#           apply(ATrial15$FitRatios, 1, FUN=mean), apply(ATrial16$FitRatios, 1, FUN=mean),
#           apply(ATrial17$FitRatios, 1, FUN=mean), apply(ATrial18$FitRatios, 1, FUN=mean), 
#           apply(CTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean),
#           apply(CTrial18$FitRatios, 1, FUN=mean),
#           apply(BTrial13$FitRatios, 1, FUN=mean), apply(BTrial14$FitRatios, 1, FUN=mean), 
#           apply(BTrial15$FitRatios, 1, FUN=mean), apply(BTrial16$FitRatios, 1, FUN=mean),
#           apply(BTrial17$FitRatios, 1, FUN=mean), apply(BTrial18$FitRatios, 1, FUN=mean), 
#           apply(DTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean),
#           apply(DTrial18$FitRatios, 1, FUN=mean)),
#         
#         c(apply(ATrial19$FitRatios, 1, FUN=mean), apply(ATrial20$FitRatios, 1, FUN=mean), 
#           apply(ATrial21$FitRatios, 1, FUN=mean), apply(ATrial22$FitRatios, 1, FUN=mean),
#           apply(ATrial23$FitRatios, 1, FUN=mean), apply(ATrial24$FitRatios, 1, FUN=mean), 
#           apply(CTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean),
#           apply(CTrial24$FitRatios, 1, FUN=mean),
#           apply(BTrial19$FitRatios, 1, FUN=mean), apply(BTrial20$FitRatios, 1, FUN=mean), 
#           apply(BTrial21$FitRatios, 1, FUN=mean), apply(BTrial22$FitRatios, 1, FUN=mean),
#           apply(BTrial23$FitRatios, 1, FUN=mean), apply(BTrial24$FitRatios, 1, FUN=mean), 
#           apply(DTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean),
#           apply(DTrial24$FitRatios, 1, FUN=mean)),
#         
#         ylim=c(0,1), col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
#         names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
# title('Fit Ratio', line=-1.4)
###############
dev.off()


png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\BackTransformed_ByFPattern.png", 
    type="cairo",
    units="mm", 
    width=250, 
    height=175, 
    pointsize=16, 
    res=600)
###############
par(mfrow=c(2,1), mar=c(0.6, 2.1, 0.3, 0.1),tcl = -0.1, mgp = c(0.8, 0.2, 0), cex=1)

vioplot(c(ATrial1$biasAll.BT, ATrial2$biasAll.BT, ATrial3$biasAll.BT, ATrial4$biasAll.BT, ATrial5$biasAll.BT, ATrial6$biasAll.BT, 
          CTrial4$biasAll.BT, CTrial5$biasAll.BT, CTrial6$biasAll.BT, 
          BTrial1$biasAll.BT, BTrial2$biasAll.BT, BTrial3$biasAll.BT, BTrial4$biasAll.BT, BTrial5$biasAll.BT, BTrial6$biasAll.BT,
          DTrial4$biasAll.BT, DTrial5$biasAll.BT, DTrial6$biasAll.BT),
        c(ATrial7$biasAll.BT, ATrial8$biasAll.BT, ATrial9$biasAll.BT, ATrial10$biasAll.BT, ATrial11$biasAll.BT, ATrial12$biasAll.BT, 
          CTrial10$biasAll.BT, CTrial11$biasAll.BT, CTrial12$biasAll.BT,  
          BTrial7$biasAll.BT, BTrial8$biasAll.BT, BTrial9$biasAll.BT, BTrial10$biasAll.BT, BTrial11$biasAll.BT, BTrial12$biasAll.BT, 
          DTrial10$biasAll.BT, DTrial11$biasAll.BT, DTrial12$biasAll.BT),
        c(ATrial13$biasAll.BT, ATrial14$biasAll.BT, ATrial15$biasAll.BT, ATrial16$biasAll.BT, ATrial17$biasAll.BT, ATrial18$biasAll.BT, 
          CTrial16$biasAll.BT, CTrial17$biasAll.BT, CTrial18$biasAll.BT, 
          BTrial13$biasAll.BT, BTrial14$biasAll.BT, BTrial15$biasAll.BT, BTrial16$biasAll.BT, BTrial17$biasAll.BT, BTrial18$biasAll.BT, 
          DTrial16$biasAll.BT, DTrial17$biasAll.BT, DTrial18$biasAll.BT),
        c(ATrial19$biasAll.BT, ATrial20$biasAll.BT, ATrial21$biasAll.BT, ATrial22$biasAll.BT, ATrial23$biasAll.BT, ATrial24$biasAll.BT, 
          CTrial22$biasAll.BT, CTrial23$biasAll.BT, CTrial24$biasAll.BT,  
          BTrial19$biasAll.BT, BTrial20$biasAll.BT, BTrial21$biasAll.BT, BTrial22$biasAll.BT, BTrial23$biasAll.BT, BTrial24$biasAll.BT, 
          DTrial22$biasAll.BT, DTrial23$biasAll.BT, DTrial24$biasAll.BT),
        ylim=c(-4.5,4.5), col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
        names=c("","","",""), ylab="")
abline(h=0)
mtext('Relative Error',side=2,line=1.1,cex=1.1)
title('Back-transformed Performance by F Pattern', line=-1.2, font.main=1, cex=1)
# abline(h=0)

par(mar=c(2.1, 2.1, 0.3, 0.1))
vioplot(c(ATrial1$RMSE.BT, ATrial2$RMSE.BT, ATrial3$RMSE.BT, ATrial4$RMSE.BT, ATrial5$RMSE.BT, ATrial6$RMSE.BT, 
          CTrial4$RMSE.BT, CTrial5$RMSE.BT, CTrial6$RMSE.BT, 
          BTrial1$RMSE.BT, BTrial2$RMSE.BT, BTrial3$RMSE.BT, BTrial4$RMSE.BT, BTrial5$RMSE.BT, BTrial6$RMSE.BT,
          DTrial4$RMSE.BT, DTrial5$RMSE.BT, DTrial6$RMSE.BT), 
        c(ATrial7$RMSE.BT, ATrial8$RMSE.BT, ATrial9$RMSE.BT, ATrial10$RMSE.BT, ATrial11$RMSE.BT, ATrial12$RMSE.BT, 
          CTrial10$RMSE.BT, CTrial11$RMSE.BT, CTrial12$RMSE.BT,  
          BTrial7$RMSE.BT, BTrial8$RMSE.BT, BTrial9$RMSE.BT, BTrial10$RMSE.BT, BTrial11$RMSE.BT, BTrial12$RMSE.BT, 
          DTrial10$RMSE.BT, DTrial11$RMSE.BT, DTrial12$RMSE.BT),
        c(ATrial13$RMSE.BT, ATrial14$RMSE.BT, ATrial15$RMSE.BT, ATrial16$RMSE.BT, ATrial17$RMSE.BT, ATrial18$RMSE.BT, 
          CTrial16$RMSE.BT, CTrial17$RMSE.BT, CTrial18$RMSE.BT, 
          BTrial13$RMSE.BT, BTrial14$RMSE.BT, BTrial15$RMSE.BT, BTrial16$RMSE.BT, BTrial17$RMSE.BT, BTrial18$RMSE.BT, 
          DTrial16$RMSE.BT, DTrial17$RMSE.BT, DTrial18$RMSE.BT),
        c(ATrial19$RMSE.BT, ATrial20$RMSE.BT, ATrial21$RMSE.BT, ATrial22$RMSE.BT, ATrial23$RMSE.BT, ATrial24$RMSE.BT, 
          CTrial22$RMSE.BT, CTrial23$RMSE.BT, CTrial24$RMSE.BT,  
          BTrial19$RMSE.BT, BTrial20$RMSE.BT, BTrial21$RMSE.BT, BTrial22$RMSE.BT, BTrial23$RMSE.BT, BTrial24$RMSE.BT, 
          DTrial22$RMSE.BT, DTrial23$RMSE.BT, DTrial24$RMSE.BT),
        ylim=c(0,2), ylab="", col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
mtext('RMSE',side=2,line=1.1,cex=1.1)
# title('RMSE', line=-1.4)

##################
dev.off()





png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\ByFPattern_2.png", 
    type="cairo",
    units="mm", 
    width=150, 
    height=200, 
    pointsize=16, 
    res=600)
############
par(mfrow=c(3,1), mar=c(1.6, 2.6, 0.6, 0.1),tcl = -0.1, mgp = c(1.1, 0.3, 0))

vioplot(c(ATrial1$biasAll, ATrial2$biasAll, ATrial3$biasAll, ATrial4$biasAll, ATrial5$biasAll, ATrial6$biasAll, 
          CTrial4$biasAll, CTrial5$biasAll, CTrial6$biasAll, 
          BTrial1$biasAll, BTrial2$biasAll, BTrial3$biasAll, BTrial4$biasAll, BTrial5$biasAll, BTrial6$biasAll,
          DTrial4$biasAll, DTrial5$biasAll, DTrial6$biasAll),
        c(ATrial7$biasAll, ATrial8$biasAll, ATrial9$biasAll, ATrial10$biasAll, ATrial11$biasAll, ATrial12$biasAll, 
          CTrial10$biasAll, CTrial11$biasAll, CTrial12$biasAll,  
          BTrial7$biasAll, BTrial8$biasAll, BTrial9$biasAll, BTrial10$biasAll, BTrial11$biasAll, BTrial12$biasAll, 
          DTrial10$biasAll, DTrial11$biasAll, DTrial12$biasAll),
        c(ATrial13$biasAll, ATrial14$biasAll, ATrial15$biasAll, ATrial16$biasAll, ATrial17$biasAll, ATrial18$biasAll, 
          CTrial16$biasAll, CTrial17$biasAll, CTrial18$biasAll, 
          BTrial13$biasAll, BTrial14$biasAll, BTrial15$biasAll, BTrial16$biasAll, BTrial17$biasAll, BTrial18$biasAll, 
          DTrial16$biasAll, DTrial17$biasAll, DTrial18$biasAll),
        c(ATrial19$biasAll, ATrial20$biasAll, ATrial21$biasAll, ATrial22$biasAll, ATrial23$biasAll, ATrial24$biasAll, 
          CTrial22$biasAll, CTrial23$biasAll, CTrial24$biasAll,  
          BTrial19$biasAll, BTrial20$biasAll, BTrial21$biasAll, BTrial22$biasAll, BTrial23$biasAll, BTrial24$biasAll, 
          DTrial22$biasAll, DTrial23$biasAll, DTrial24$biasAll),
        ylim=c(-4.5,4.5), col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'),
        ylab="")
mtext('Relative Error', side=2,line=1.2)
# abline(h=0)

vioplot(c(ATrial1$RMSE, ATrial2$RMSE, ATrial3$RMSE, ATrial4$RMSE, ATrial5$RMSE, ATrial6$RMSE, 
          CTrial4$RMSE, CTrial5$RMSE, CTrial6$RMSE, 
          BTrial1$RMSE, BTrial2$RMSE, BTrial3$RMSE, BTrial4$RMSE, BTrial5$RMSE, BTrial6$RMSE,
          DTrial4$RMSE, DTrial5$RMSE, DTrial6$RMSE), 
        c(ATrial7$RMSE, ATrial8$RMSE, ATrial9$RMSE, ATrial10$RMSE, ATrial11$RMSE, ATrial12$RMSE, 
          CTrial10$RMSE, CTrial11$RMSE, CTrial12$RMSE,  
          BTrial7$RMSE, BTrial8$RMSE, BTrial9$RMSE, BTrial10$RMSE, BTrial11$RMSE, BTrial12$RMSE, 
          DTrial10$RMSE, DTrial11$RMSE, DTrial12$RMSE),
        c(ATrial13$RMSE, ATrial14$RMSE, ATrial15$RMSE, ATrial16$RMSE, ATrial17$RMSE, ATrial18$RMSE, 
          CTrial16$RMSE, CTrial17$RMSE, CTrial18$RMSE, 
          BTrial13$RMSE, BTrial14$RMSE, BTrial15$RMSE, BTrial16$RMSE, BTrial17$RMSE, BTrial18$RMSE, 
          DTrial16$RMSE, DTrial17$RMSE, DTrial18$RMSE),
        c(ATrial19$RMSE, ATrial20$RMSE, ATrial21$RMSE, ATrial22$RMSE, ATrial23$RMSE, ATrial24$RMSE, 
          CTrial22$RMSE, CTrial23$RMSE, CTrial24$RMSE,  
          BTrial19$RMSE, BTrial20$RMSE, BTrial21$RMSE, BTrial22$RMSE, BTrial23$RMSE, BTrial24$RMSE, 
          DTrial22$RMSE, DTrial23$RMSE, DTrial24$RMSE),
        ylim=c(0,2), col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
mtext('RMSE', side=2,line=1.2)




vioplot(c(apply(ATrial1$FitRatios, 1, FUN=mean), apply(ATrial2$FitRatios, 1, FUN=mean), 
          apply(ATrial3$FitRatios, 1, FUN=mean), apply(ATrial4$FitRatios, 1, FUN=mean),
          apply(ATrial5$FitRatios, 1, FUN=mean), apply(ATrial6$FitRatios, 1, FUN=mean), 
          apply(CTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean),
          apply(CTrial6$FitRatios, 1, FUN=mean),
          apply(BTrial1$FitRatios, 1, FUN=mean), apply(BTrial2$FitRatios, 1, FUN=mean), 
          apply(BTrial3$FitRatios, 1, FUN=mean), apply(BTrial4$FitRatios, 1, FUN=mean),
          apply(BTrial5$FitRatios, 1, FUN=mean), apply(BTrial6$FitRatios, 1, FUN=mean), 
          apply(DTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean),
          apply(DTrial6$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial7$FitRatios, 1, FUN=mean), apply(ATrial8$FitRatios, 1, FUN=mean), 
          apply(ATrial9$FitRatios, 1, FUN=mean), apply(ATrial10$FitRatios, 1, FUN=mean),
          apply(ATrial11$FitRatios, 1, FUN=mean), apply(ATrial12$FitRatios, 1, FUN=mean), 
          apply(CTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean),
          apply(CTrial12$FitRatios, 1, FUN=mean),
          apply(BTrial7$FitRatios, 1, FUN=mean), apply(BTrial8$FitRatios, 1, FUN=mean), 
          apply(BTrial9$FitRatios, 1, FUN=mean), apply(BTrial10$FitRatios, 1, FUN=mean),
          apply(BTrial11$FitRatios, 1, FUN=mean), apply(BTrial12$FitRatios, 1, FUN=mean), 
          apply(DTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean),
          apply(DTrial12$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial13$FitRatios, 1, FUN=mean), apply(ATrial14$FitRatios, 1, FUN=mean), 
          apply(ATrial15$FitRatios, 1, FUN=mean), apply(ATrial16$FitRatios, 1, FUN=mean),
          apply(ATrial17$FitRatios, 1, FUN=mean), apply(ATrial18$FitRatios, 1, FUN=mean), 
          apply(CTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean),
          apply(CTrial18$FitRatios, 1, FUN=mean),
          apply(BTrial13$FitRatios, 1, FUN=mean), apply(BTrial14$FitRatios, 1, FUN=mean), 
          apply(BTrial15$FitRatios, 1, FUN=mean), apply(BTrial16$FitRatios, 1, FUN=mean),
          apply(BTrial17$FitRatios, 1, FUN=mean), apply(BTrial18$FitRatios, 1, FUN=mean), 
          apply(DTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean),
          apply(DTrial18$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial19$FitRatios, 1, FUN=mean), apply(ATrial20$FitRatios, 1, FUN=mean), 
          apply(ATrial21$FitRatios, 1, FUN=mean), apply(ATrial22$FitRatios, 1, FUN=mean),
          apply(ATrial23$FitRatios, 1, FUN=mean), apply(ATrial24$FitRatios, 1, FUN=mean), 
          apply(CTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean),
          apply(CTrial24$FitRatios, 1, FUN=mean),
          apply(BTrial19$FitRatios, 1, FUN=mean), apply(BTrial20$FitRatios, 1, FUN=mean), 
          apply(BTrial21$FitRatios, 1, FUN=mean), apply(BTrial22$FitRatios, 1, FUN=mean),
          apply(BTrial23$FitRatios, 1, FUN=mean), apply(BTrial24$FitRatios, 1, FUN=mean), 
          apply(DTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean),
          apply(DTrial24$FitRatios, 1, FUN=mean)),
        
        ylim=c(0,0.8), col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
mtext('Mean Fit Ratio', side=2,line=1.2)
###################

dev.off()



#### BY F PATTERN ##########
png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\ByF_#S_q.png", 
    type="cairo",
    units="mm", 
    width=250, 
    height=300, 
    pointsize=16, 
    res=600)
###################
par(mfrow=c(3,1), mar=c(0.6, 1.6, 0.6, 0.1),tcl = -0.1, mgp = c(2, 0.3, 0))

vioplot(c(ATrial1$biasAll, ATrial2$biasAll, ATrial3$biasAll), 
        c(ATrial4$biasAll, ATrial5$biasAll, ATrial6$biasAll, CTrial4$biasAll, CTrial5$biasAll, CTrial6$biasAll), 
        c(BTrial1$biasAll, BTrial2$biasAll, BTrial3$biasAll), 
        c(BTrial4$biasAll, BTrial5$biasAll, BTrial6$biasAll, DTrial4$biasAll, DTrial5$biasAll, DTrial6$biasAll),
        c(ATrial7$biasAll, ATrial8$biasAll, ATrial9$biasAll), 
        c(ATrial10$biasAll, ATrial11$biasAll, ATrial12$biasAll, CTrial10$biasAll, CTrial11$biasAll, CTrial12$biasAll),  
        c(BTrial7$biasAll, BTrial8$biasAll, BTrial9$biasAll), 
        c(BTrial10$biasAll, BTrial11$biasAll, BTrial12$biasAll, DTrial10$biasAll, DTrial11$biasAll, DTrial12$biasAll),
        c(ATrial13$biasAll, ATrial14$biasAll, ATrial15$biasAll), 
        c(ATrial16$biasAll, ATrial17$biasAll, ATrial18$biasAll, CTrial16$biasAll, CTrial17$biasAll, CTrial18$biasAll), 
        c(BTrial13$biasAll, BTrial14$biasAll, BTrial15$biasAll), 
        c(BTrial16$biasAll, BTrial17$biasAll, BTrial18$biasAll, DTrial16$biasAll, DTrial17$biasAll, DTrial18$biasAll),
        c(ATrial19$biasAll, ATrial20$biasAll, ATrial21$biasAll), 
        c(ATrial22$biasAll, ATrial23$biasAll, ATrial24$biasAll, CTrial22$biasAll, CTrial23$biasAll, CTrial24$biasAll),  
        c(BTrial19$biasAll, BTrial20$biasAll, BTrial21$biasAll), 
        c(BTrial22$biasAll, BTrial23$biasAll, BTrial24$biasAll, DTrial22$biasAll, DTrial23$biasAll, DTrial24$biasAll),
        ylim=c(-4.5,4.5), col=c('darkolivegreen1', 'darkolivegreen3','darkolivegreen4','darkolivegreen',
                            'skyblue','deepskyblue2','deepskyblue3','deepskyblue4',
                            'orange','darkorange','darkorange2','darkorange3',
                            'mediumorchid','darkorchid1','darkorchid','darkorchid4'), names=FALSE)
title('Relative Error', line=-1.4)
# abline(h=0)

vioplot(c(ATrial1$RMSE, ATrial2$RMSE, ATrial3$RMSE), 
        c(ATrial4$RMSE, ATrial5$RMSE, ATrial6$RMSE, CTrial4$RMSE, CTrial5$RMSE, CTrial6$RMSE), 
        c(BTrial1$RMSE, BTrial2$RMSE, BTrial3$RMSE), 
        c(BTrial4$RMSE, BTrial5$RMSE, BTrial6$RMSE, DTrial4$RMSE, DTrial5$RMSE, DTrial6$RMSE), 
        c(ATrial7$RMSE, ATrial8$RMSE, ATrial9$RMSE), 
        c(ATrial10$RMSE, ATrial11$RMSE, ATrial12$RMSE, CTrial10$RMSE, CTrial11$RMSE, CTrial12$RMSE),  
        c(BTrial7$RMSE, BTrial8$RMSE, BTrial9$RMSE), 
        c(BTrial10$RMSE, BTrial11$RMSE, BTrial12$RMSE, DTrial10$RMSE, DTrial11$RMSE, DTrial12$RMSE),
        c(ATrial13$RMSE, ATrial14$RMSE, ATrial15$RMSE), 
        c(ATrial16$RMSE, ATrial17$RMSE, ATrial18$RMSE, CTrial16$RMSE, CTrial17$RMSE, CTrial18$RMSE), 
        c(BTrial13$RMSE, BTrial14$RMSE, BTrial15$RMSE), 
        c(BTrial16$RMSE, BTrial17$RMSE, BTrial18$RMSE, DTrial16$RMSE, DTrial17$RMSE, DTrial18$RMSE),
        c(ATrial19$RMSE, ATrial20$RMSE, ATrial21$RMSE), 
        c(ATrial22$RMSE, ATrial23$RMSE, ATrial24$RMSE, CTrial22$RMSE, CTrial23$RMSE, CTrial24$RMSE),  
        c(BTrial19$RMSE, BTrial20$RMSE, BTrial21$RMSE), 
        c(BTrial22$RMSE, BTrial23$RMSE, BTrial24$RMSE, DTrial22$RMSE, DTrial23$RMSE, DTrial24$RMSE),
        ylim=c(0,2), col=c('darkolivegreen1', 'darkolivegreen3','darkolivegreen4','darkolivegreen',
                           'skyblue','deepskyblue2','deepskyblue3','deepskyblue4',
                           'orange','darkorange','darkorange2','darkorange3',
                           'mediumorchid','darkorchid1','darkorchid','darkorchid4'), names=FALSE)
title('RMSE', line=-1.4)


par( mar=c(1.6, 1.6, 0.6, 0.1))
vioplot(c(apply(ATrial1$FitRatios, 1, FUN=mean), apply(ATrial2$FitRatios, 1, FUN=mean), apply(ATrial3$FitRatios, 1, FUN=mean)), 
        c(apply(ATrial4$FitRatios, 1, FUN=mean), apply(ATrial5$FitRatios, 1, FUN=mean), apply(ATrial6$FitRatios, 1, FUN=mean), 
          apply(CTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean), apply(CTrial6$FitRatios, 1, FUN=mean)),
        c(apply(BTrial1$FitRatios, 1, FUN=mean), apply(BTrial2$FitRatios, 1, FUN=mean), apply(BTrial3$FitRatios, 1, FUN=mean)), 
        c(apply(BTrial4$FitRatios, 1, FUN=mean), apply(BTrial5$FitRatios, 1, FUN=mean), apply(BTrial6$FitRatios, 1, FUN=mean), 
          apply(DTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean), apply(DTrial6$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial7$FitRatios, 1, FUN=mean), apply(ATrial8$FitRatios, 1, FUN=mean), apply(ATrial9$FitRatios, 1, FUN=mean)), 
        c(apply(ATrial10$FitRatios, 1, FUN=mean), apply(ATrial11$FitRatios, 1, FUN=mean), apply(ATrial12$FitRatios, 1, FUN=mean), 
          apply(CTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean), apply(CTrial12$FitRatios, 1, FUN=mean)),
        c(apply(BTrial7$FitRatios, 1, FUN=mean), apply(BTrial8$FitRatios, 1, FUN=mean), apply(BTrial9$FitRatios, 1, FUN=mean)), 
        c(apply(BTrial10$FitRatios, 1, FUN=mean), apply(BTrial11$FitRatios, 1, FUN=mean), apply(BTrial12$FitRatios, 1, FUN=mean), 
          apply(DTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean), apply(DTrial12$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial13$FitRatios, 1, FUN=mean), apply(ATrial14$FitRatios, 1, FUN=mean), apply(ATrial15$FitRatios, 1, FUN=mean)), 
        c(apply(ATrial16$FitRatios, 1, FUN=mean), apply(ATrial17$FitRatios, 1, FUN=mean), apply(ATrial18$FitRatios, 1, FUN=mean), 
          apply(CTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean), apply(CTrial18$FitRatios, 1, FUN=mean)),
        c(apply(BTrial13$FitRatios, 1, FUN=mean), apply(BTrial14$FitRatios, 1, FUN=mean), apply(BTrial15$FitRatios, 1, FUN=mean)), 
        c(apply(BTrial16$FitRatios, 1, FUN=mean), apply(BTrial17$FitRatios, 1, FUN=mean), apply(BTrial18$FitRatios, 1, FUN=mean), 
          apply(DTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean), apply(DTrial18$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial19$FitRatios, 1, FUN=mean), apply(ATrial20$FitRatios, 1, FUN=mean), apply(ATrial21$FitRatios, 1, FUN=mean)), 
        c(apply(ATrial22$FitRatios, 1, FUN=mean), apply(ATrial23$FitRatios, 1, FUN=mean), apply(ATrial24$FitRatios, 1, FUN=mean), 
          apply(CTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean), apply(CTrial24$FitRatios, 1, FUN=mean)),
        c(apply(BTrial19$FitRatios, 1, FUN=mean), apply(BTrial20$FitRatios, 1, FUN=mean), apply(BTrial21$FitRatios, 1, FUN=mean)), 
        c(apply(BTrial22$FitRatios, 1, FUN=mean), apply(BTrial23$FitRatios, 1, FUN=mean), apply(BTrial24$FitRatios, 1, FUN=mean), 
          apply(DTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean), apply(DTrial24$FitRatios, 1, FUN=mean)),
        
        ylim=c(0,0.8), col=c('darkolivegreen1', 'darkolivegreen3','darkolivegreen4','darkolivegreen',
                            'skyblue','deepskyblue2','deepskyblue3','deepskyblue4',
                            'orange','darkorange','darkorange2','darkorange3',
                            'mediumorchid','darkorchid1','darkorchid','darkorchid4'),
        names=c('3S', expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')),
                '3S',expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')),
                '3S',expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')),
                '3S',expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')) ))
title('Fit Ratio', line=-1.4)
#################
dev.off()




#### BY F PATTERN ##########
png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\ByF_#S_q2.png", 
    type="cairo",
    units="mm", 
    width=250, 
    height=200, 
    pointsize=16, 
    res=600)
##########
par(mfrow=c(2,1), mar=c(0.6, 1.6, 0.6, 0.1),tcl = -0.1, mgp = c(2, 0.3, 0))

vioplot(c(ATrial1$biasAll, ATrial2$biasAll, ATrial3$biasAll), 
        c(ATrial4$biasAll, ATrial5$biasAll, ATrial6$biasAll, CTrial4$biasAll, CTrial5$biasAll, CTrial6$biasAll), 
        c(BTrial1$biasAll, BTrial2$biasAll, BTrial3$biasAll), 
        c(BTrial4$biasAll, BTrial5$biasAll, BTrial6$biasAll, DTrial4$biasAll, DTrial5$biasAll, DTrial6$biasAll),
        c(ATrial7$biasAll, ATrial8$biasAll, ATrial9$biasAll), 
        c(ATrial10$biasAll, ATrial11$biasAll, ATrial12$biasAll, CTrial10$biasAll, CTrial11$biasAll, CTrial12$biasAll),  
        c(BTrial7$biasAll, BTrial8$biasAll, BTrial9$biasAll), 
        c(BTrial10$biasAll, BTrial11$biasAll, BTrial12$biasAll, DTrial10$biasAll, DTrial11$biasAll, DTrial12$biasAll),
        c(ATrial13$biasAll, ATrial14$biasAll, ATrial15$biasAll), 
        c(ATrial16$biasAll, ATrial17$biasAll, ATrial18$biasAll, CTrial16$biasAll, CTrial17$biasAll, CTrial18$biasAll), 
        c(BTrial13$biasAll, BTrial14$biasAll, BTrial15$biasAll), 
        c(BTrial16$biasAll, BTrial17$biasAll, BTrial18$biasAll, DTrial16$biasAll, DTrial17$biasAll, DTrial18$biasAll),
        c(ATrial19$biasAll, ATrial20$biasAll, ATrial21$biasAll), 
        c(ATrial22$biasAll, ATrial23$biasAll, ATrial24$biasAll, CTrial22$biasAll, CTrial23$biasAll, CTrial24$biasAll),  
        c(BTrial19$biasAll, BTrial20$biasAll, BTrial21$biasAll), 
        c(BTrial22$biasAll, BTrial23$biasAll, BTrial24$biasAll, DTrial22$biasAll, DTrial23$biasAll, DTrial24$biasAll),
        ylim=c(-4.5,4.5), col=c('darkolivegreen1', 'darkolivegreen3','darkolivegreen4','darkolivegreen',
                            'skyblue','deepskyblue2','deepskyblue3','deepskyblue4',
                            'orange','darkorange','darkorange2','darkorange3',
                            'mediumorchid','darkorchid1','darkorchid','darkorchid4'), names=FALSE)
mtext('Relative Error', 2, line=1.1)
mtext('Atlantic sharpnose shark', 3, line=-1.1)
# abline(h=0)

par( mar=c(1.6, 1.6, 0.6, 0.1))

vioplot(c(ATrial1$RMSE, ATrial2$RMSE, ATrial3$RMSE), 
        c(ATrial4$RMSE, ATrial5$RMSE, ATrial6$RMSE, CTrial4$RMSE, CTrial5$RMSE, CTrial6$RMSE), 
        c(BTrial1$RMSE, BTrial2$RMSE, BTrial3$RMSE), 
        c(BTrial4$RMSE, BTrial5$RMSE, BTrial6$RMSE, DTrial4$RMSE, DTrial5$RMSE, DTrial6$RMSE), 
        c(ATrial7$RMSE, ATrial8$RMSE, ATrial9$RMSE), 
        c(ATrial10$RMSE, ATrial11$RMSE, ATrial12$RMSE, CTrial10$RMSE, CTrial11$RMSE, CTrial12$RMSE),  
        c(BTrial7$RMSE, BTrial8$RMSE, BTrial9$RMSE), 
        c(BTrial10$RMSE, BTrial11$RMSE, BTrial12$RMSE, DTrial10$RMSE, DTrial11$RMSE, DTrial12$RMSE),
        c(ATrial13$RMSE, ATrial14$RMSE, ATrial15$RMSE), 
        c(ATrial16$RMSE, ATrial17$RMSE, ATrial18$RMSE, CTrial16$RMSE, CTrial17$RMSE, CTrial18$RMSE), 
        c(BTrial13$RMSE, BTrial14$RMSE, BTrial15$RMSE), 
        c(BTrial16$RMSE, BTrial17$RMSE, BTrial18$RMSE, DTrial16$RMSE, DTrial17$RMSE, DTrial18$RMSE),
        c(ATrial19$RMSE, ATrial20$RMSE, ATrial21$RMSE), 
        c(ATrial22$RMSE, ATrial23$RMSE, ATrial24$RMSE, CTrial22$RMSE, CTrial23$RMSE, CTrial24$RMSE),  
        c(BTrial19$RMSE, BTrial20$RMSE, BTrial21$RMSE), 
        c(BTrial22$RMSE, BTrial23$RMSE, BTrial24$RMSE, DTrial22$RMSE, DTrial23$RMSE, DTrial24$RMSE),
        ylim=c(0,2), col=c('darkolivegreen1', 'darkolivegreen3','darkolivegreen4','darkolivegreen',
                           'skyblue','deepskyblue2','deepskyblue3','deepskyblue4',
                           'orange','darkorange','darkorange2','darkorange3',
                           'mediumorchid','darkorchid1','darkorchid','darkorchid4'), 
        names=c('3S', expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')),
                '3S',expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')),
                '3S',expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')),
                '3S',expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')) ))
mtext('RMSE', 2, line=1.1)

##########

dev.off()



#### BY F PATTERN ##########
png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\ByFPattern.png", 
    type="cairo",
    units="mm", 
    width=500, 
    height=167, 
    pointsize=16, 
    res=600)
#############
par(mfrow=c(1,3), mar=c(1.6, 1.6, 1.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))

vioplot(c(ATrial1$biasAll, ATrial2$biasAll, ATrial3$biasAll, ATrial4$biasAll, ATrial5$biasAll, ATrial6$biasAll, 
          CTrial4$biasAll, CTrial5$biasAll, CTrial6$biasAll, 
          BTrial1$biasAll, BTrial2$biasAll, BTrial3$biasAll, BTrial4$biasAll, BTrial5$biasAll, BTrial6$biasAll,
          DTrial4$biasAll, DTrial5$biasAll, DTrial6$biasAll),
        c(ATrial7$biasAll, ATrial8$biasAll, ATrial9$biasAll, ATrial10$biasAll, ATrial11$biasAll, ATrial12$biasAll, 
          CTrial10$biasAll, CTrial11$biasAll, CTrial12$biasAll,  
          BTrial7$biasAll, BTrial8$biasAll, BTrial9$biasAll, BTrial10$biasAll, BTrial11$biasAll, BTrial12$biasAll, 
          DTrial10$biasAll, DTrial11$biasAll, DTrial12$biasAll),
        c(ATrial13$biasAll, ATrial14$biasAll, ATrial15$biasAll, ATrial16$biasAll, ATrial17$biasAll, ATrial18$biasAll, 
          CTrial16$biasAll, CTrial17$biasAll, CTrial18$biasAll, 
          BTrial13$biasAll, BTrial14$biasAll, BTrial15$biasAll, BTrial16$biasAll, BTrial17$biasAll, BTrial18$biasAll, 
          DTrial16$biasAll, DTrial17$biasAll, DTrial18$biasAll),
        c(ATrial19$biasAll, ATrial20$biasAll, ATrial21$biasAll, ATrial22$biasAll, ATrial23$biasAll, ATrial24$biasAll, 
          CTrial22$biasAll, CTrial23$biasAll, CTrial24$biasAll,  
          BTrial19$biasAll, BTrial20$biasAll, BTrial21$biasAll, BTrial22$biasAll, BTrial23$biasAll, BTrial24$biasAll, 
          DTrial22$biasAll, DTrial23$biasAll, DTrial24$biasAll),
        ylim=c(-4.5,4.5), col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
title('Relative Error', line=-1.4)
# abline(h=0)

vioplot(c(ATrial1$RMSE, ATrial2$RMSE, ATrial3$RMSE, ATrial4$RMSE, ATrial5$RMSE, ATrial6$RMSE, 
          CTrial4$RMSE, CTrial5$RMSE, CTrial6$RMSE, 
          BTrial1$RMSE, BTrial2$RMSE, BTrial3$RMSE, BTrial4$RMSE, BTrial5$RMSE, BTrial6$RMSE,
          DTrial4$RMSE, DTrial5$RMSE, DTrial6$RMSE), 
        c(ATrial7$RMSE, ATrial8$RMSE, ATrial9$RMSE, ATrial10$RMSE, ATrial11$RMSE, ATrial12$RMSE, 
          CTrial10$RMSE, CTrial11$RMSE, CTrial12$RMSE,  
          BTrial7$RMSE, BTrial8$RMSE, BTrial9$RMSE, BTrial10$RMSE, BTrial11$RMSE, BTrial12$RMSE, 
          DTrial10$RMSE, DTrial11$RMSE, DTrial12$RMSE),
        c(ATrial13$RMSE, ATrial14$RMSE, ATrial15$RMSE, ATrial16$RMSE, ATrial17$RMSE, ATrial18$RMSE, 
          CTrial16$RMSE, CTrial17$RMSE, CTrial18$RMSE, 
          BTrial13$RMSE, BTrial14$RMSE, BTrial15$RMSE, BTrial16$RMSE, BTrial17$RMSE, BTrial18$RMSE, 
          DTrial16$RMSE, DTrial17$RMSE, DTrial18$RMSE),
        c(ATrial19$RMSE, ATrial20$RMSE, ATrial21$RMSE, ATrial22$RMSE, ATrial23$RMSE, ATrial24$RMSE, 
          CTrial22$RMSE, CTrial23$RMSE, CTrial24$RMSE,  
          BTrial19$RMSE, BTrial20$RMSE, BTrial21$RMSE, BTrial22$RMSE, BTrial23$RMSE, BTrial24$RMSE, 
          DTrial22$RMSE, DTrial23$RMSE, DTrial24$RMSE),
        ylim=c(0,2), col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
title('RMSE', line=-1.4)



vioplot(c(apply(ATrial1$FitRatios, 1, FUN=mean), apply(ATrial2$FitRatios, 1, FUN=mean), 
          apply(ATrial3$FitRatios, 1, FUN=mean), apply(ATrial4$FitRatios, 1, FUN=mean),
          apply(ATrial5$FitRatios, 1, FUN=mean), apply(ATrial6$FitRatios, 1, FUN=mean), 
          apply(CTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean),
          apply(CTrial6$FitRatios, 1, FUN=mean),
          apply(BTrial1$FitRatios, 1, FUN=mean), apply(BTrial2$FitRatios, 1, FUN=mean), 
          apply(BTrial3$FitRatios, 1, FUN=mean), apply(BTrial4$FitRatios, 1, FUN=mean),
          apply(BTrial5$FitRatios, 1, FUN=mean), apply(BTrial6$FitRatios, 1, FUN=mean), 
          apply(DTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean),
          apply(DTrial6$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial7$FitRatios, 1, FUN=mean), apply(ATrial8$FitRatios, 1, FUN=mean), 
          apply(ATrial9$FitRatios, 1, FUN=mean), apply(ATrial10$FitRatios, 1, FUN=mean),
          apply(ATrial11$FitRatios, 1, FUN=mean), apply(ATrial12$FitRatios, 1, FUN=mean), 
          apply(CTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean),
          apply(CTrial12$FitRatios, 1, FUN=mean),
          apply(BTrial7$FitRatios, 1, FUN=mean), apply(BTrial8$FitRatios, 1, FUN=mean), 
          apply(BTrial9$FitRatios, 1, FUN=mean), apply(BTrial10$FitRatios, 1, FUN=mean),
          apply(BTrial11$FitRatios, 1, FUN=mean), apply(BTrial12$FitRatios, 1, FUN=mean), 
          apply(DTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean),
          apply(DTrial12$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial13$FitRatios, 1, FUN=mean), apply(ATrial14$FitRatios, 1, FUN=mean), 
          apply(ATrial15$FitRatios, 1, FUN=mean), apply(ATrial16$FitRatios, 1, FUN=mean),
          apply(ATrial17$FitRatios, 1, FUN=mean), apply(ATrial18$FitRatios, 1, FUN=mean), 
          apply(CTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean),
          apply(CTrial18$FitRatios, 1, FUN=mean),
          apply(BTrial13$FitRatios, 1, FUN=mean), apply(BTrial14$FitRatios, 1, FUN=mean), 
          apply(BTrial15$FitRatios, 1, FUN=mean), apply(BTrial16$FitRatios, 1, FUN=mean),
          apply(BTrial17$FitRatios, 1, FUN=mean), apply(BTrial18$FitRatios, 1, FUN=mean), 
          apply(DTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean),
          apply(DTrial18$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial19$FitRatios, 1, FUN=mean), apply(ATrial20$FitRatios, 1, FUN=mean), 
          apply(ATrial21$FitRatios, 1, FUN=mean), apply(ATrial22$FitRatios, 1, FUN=mean),
          apply(ATrial23$FitRatios, 1, FUN=mean), apply(ATrial24$FitRatios, 1, FUN=mean), 
          apply(CTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean),
          apply(CTrial24$FitRatios, 1, FUN=mean),
          apply(BTrial19$FitRatios, 1, FUN=mean), apply(BTrial20$FitRatios, 1, FUN=mean), 
          apply(BTrial21$FitRatios, 1, FUN=mean), apply(BTrial22$FitRatios, 1, FUN=mean),
          apply(BTrial23$FitRatios, 1, FUN=mean), apply(BTrial24$FitRatios, 1, FUN=mean), 
          apply(DTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean),
          apply(DTrial24$FitRatios, 1, FUN=mean)),
        
        col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F')) #ylim=c(0,2), 
title('Fit Ratio', line=-1.4)
#################
dev.off()


#### BY F PATTERN ##########
png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\ByFPattern_LgFont.png", 
    type="cairo",
    units="mm", 
    width=250, 
    height=175, 
    pointsize=16, 
    res=600)
################
par(mfrow=c(2,1), mar=c(0.6, 2.1, 0.3, 0.1),tcl = -0.1, mgp = c(0.8, 0.2, 0), cex=1)

vioplot(c(ATrial1$biasAll, ATrial2$biasAll, ATrial3$biasAll, ATrial4$biasAll, ATrial5$biasAll, ATrial6$biasAll, 
          CTrial4$biasAll, CTrial5$biasAll, CTrial6$biasAll, 
          BTrial1$biasAll, BTrial2$biasAll, BTrial3$biasAll, BTrial4$biasAll, BTrial5$biasAll, BTrial6$biasAll,
          DTrial4$biasAll, DTrial5$biasAll, DTrial6$biasAll),
        c(ATrial7$biasAll, ATrial8$biasAll, ATrial9$biasAll, ATrial10$biasAll, ATrial11$biasAll, ATrial12$biasAll, 
          CTrial10$biasAll, CTrial11$biasAll, CTrial12$biasAll,  
          BTrial7$biasAll, BTrial8$biasAll, BTrial9$biasAll, BTrial10$biasAll, BTrial11$biasAll, BTrial12$biasAll, 
          DTrial10$biasAll, DTrial11$biasAll, DTrial12$biasAll),
        c(ATrial13$biasAll, ATrial14$biasAll, ATrial15$biasAll, ATrial16$biasAll, ATrial17$biasAll, ATrial18$biasAll, 
          CTrial16$biasAll, CTrial17$biasAll, CTrial18$biasAll, 
          BTrial13$biasAll, BTrial14$biasAll, BTrial15$biasAll, BTrial16$biasAll, BTrial17$biasAll, BTrial18$biasAll, 
          DTrial16$biasAll, DTrial17$biasAll, DTrial18$biasAll),
        c(ATrial19$biasAll, ATrial20$biasAll, ATrial21$biasAll, ATrial22$biasAll, ATrial23$biasAll, ATrial24$biasAll, 
          CTrial22$biasAll, CTrial23$biasAll, CTrial24$biasAll,  
          BTrial19$biasAll, BTrial20$biasAll, BTrial21$biasAll, BTrial22$biasAll, BTrial23$biasAll, BTrial24$biasAll, 
          DTrial22$biasAll, DTrial23$biasAll, DTrial24$biasAll),
        ylim=c(-4.5,4.5), col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
        names=c("","","",""), ylab=" ")
abline(h=0)
mtext('Relative Error',side=2,line=1.1,cex=1.1)
title('Performance by F Pattern', line=-1.2, font.main=1, cex=1)
# abline(h=0)

par(mar=c(2.1, 2.1, 0.3, 0.1))
vioplot(c(ATrial1$RMSE, ATrial2$RMSE, ATrial3$RMSE, ATrial4$RMSE, ATrial5$RMSE, ATrial6$RMSE, 
          CTrial4$RMSE, CTrial5$RMSE, CTrial6$RMSE, 
          BTrial1$RMSE, BTrial2$RMSE, BTrial3$RMSE, BTrial4$RMSE, BTrial5$RMSE, BTrial6$RMSE,
          DTrial4$RMSE, DTrial5$RMSE, DTrial6$RMSE), 
        c(ATrial7$RMSE, ATrial8$RMSE, ATrial9$RMSE, ATrial10$RMSE, ATrial11$RMSE, ATrial12$RMSE, 
          CTrial10$RMSE, CTrial11$RMSE, CTrial12$RMSE,  
          BTrial7$RMSE, BTrial8$RMSE, BTrial9$RMSE, BTrial10$RMSE, BTrial11$RMSE, BTrial12$RMSE, 
          DTrial10$RMSE, DTrial11$RMSE, DTrial12$RMSE),
        c(ATrial13$RMSE, ATrial14$RMSE, ATrial15$RMSE, ATrial16$RMSE, ATrial17$RMSE, ATrial18$RMSE, 
          CTrial16$RMSE, CTrial17$RMSE, CTrial18$RMSE, 
          BTrial13$RMSE, BTrial14$RMSE, BTrial15$RMSE, BTrial16$RMSE, BTrial17$RMSE, BTrial18$RMSE, 
          DTrial16$RMSE, DTrial17$RMSE, DTrial18$RMSE),
        c(ATrial19$RMSE, ATrial20$RMSE, ATrial21$RMSE, ATrial22$RMSE, ATrial23$RMSE, ATrial24$RMSE, 
          CTrial22$RMSE, CTrial23$RMSE, CTrial24$RMSE,  
          BTrial19$RMSE, BTrial20$RMSE, BTrial21$RMSE, BTrial22$RMSE, BTrial23$RMSE, BTrial24$RMSE, 
          DTrial22$RMSE, DTrial23$RMSE, DTrial24$RMSE),
        ylim=c(0,2), ylab="", col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
mtext('RMSE',side=2,line=1.1,cex=1.1)
# title('RMSE', line=-1.4)


# 
# vioplot(c(apply(ATrial1$FitRatios, 1, FUN=mean), apply(ATrial2$FitRatios, 1, FUN=mean), 
#           apply(ATrial3$FitRatios, 1, FUN=mean), apply(ATrial4$FitRatios, 1, FUN=mean),
#           apply(ATrial5$FitRatios, 1, FUN=mean), apply(ATrial6$FitRatios, 1, FUN=mean), 
#           apply(CTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean),
#           apply(CTrial6$FitRatios, 1, FUN=mean),
#           apply(BTrial1$FitRatios, 1, FUN=mean), apply(BTrial2$FitRatios, 1, FUN=mean), 
#           apply(BTrial3$FitRatios, 1, FUN=mean), apply(BTrial4$FitRatios, 1, FUN=mean),
#           apply(BTrial5$FitRatios, 1, FUN=mean), apply(BTrial6$FitRatios, 1, FUN=mean), 
#           apply(DTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean),
#           apply(DTrial6$FitRatios, 1, FUN=mean)),
#         
#         c(apply(ATrial7$FitRatios, 1, FUN=mean), apply(ATrial8$FitRatios, 1, FUN=mean), 
#           apply(ATrial9$FitRatios, 1, FUN=mean), apply(ATrial10$FitRatios, 1, FUN=mean),
#           apply(ATrial11$FitRatios, 1, FUN=mean), apply(ATrial12$FitRatios, 1, FUN=mean), 
#           apply(CTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean),
#           apply(CTrial12$FitRatios, 1, FUN=mean),
#           apply(BTrial7$FitRatios, 1, FUN=mean), apply(BTrial8$FitRatios, 1, FUN=mean), 
#           apply(BTrial9$FitRatios, 1, FUN=mean), apply(BTrial10$FitRatios, 1, FUN=mean),
#           apply(BTrial11$FitRatios, 1, FUN=mean), apply(BTrial12$FitRatios, 1, FUN=mean), 
#           apply(DTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean),
#           apply(DTrial12$FitRatios, 1, FUN=mean)),
#         
#         c(apply(ATrial13$FitRatios, 1, FUN=mean), apply(ATrial14$FitRatios, 1, FUN=mean), 
#           apply(ATrial15$FitRatios, 1, FUN=mean), apply(ATrial16$FitRatios, 1, FUN=mean),
#           apply(ATrial17$FitRatios, 1, FUN=mean), apply(ATrial18$FitRatios, 1, FUN=mean), 
#           apply(CTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean),
#           apply(CTrial18$FitRatios, 1, FUN=mean),
#           apply(BTrial13$FitRatios, 1, FUN=mean), apply(BTrial14$FitRatios, 1, FUN=mean), 
#           apply(BTrial15$FitRatios, 1, FUN=mean), apply(BTrial16$FitRatios, 1, FUN=mean),
#           apply(BTrial17$FitRatios, 1, FUN=mean), apply(BTrial18$FitRatios, 1, FUN=mean), 
#           apply(DTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean),
#           apply(DTrial18$FitRatios, 1, FUN=mean)),
#         
#         c(apply(ATrial19$FitRatios, 1, FUN=mean), apply(ATrial20$FitRatios, 1, FUN=mean), 
#           apply(ATrial21$FitRatios, 1, FUN=mean), apply(ATrial22$FitRatios, 1, FUN=mean),
#           apply(ATrial23$FitRatios, 1, FUN=mean), apply(ATrial24$FitRatios, 1, FUN=mean), 
#           apply(CTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean),
#           apply(CTrial24$FitRatios, 1, FUN=mean),
#           apply(BTrial19$FitRatios, 1, FUN=mean), apply(BTrial20$FitRatios, 1, FUN=mean), 
#           apply(BTrial21$FitRatios, 1, FUN=mean), apply(BTrial22$FitRatios, 1, FUN=mean),
#           apply(BTrial23$FitRatios, 1, FUN=mean), apply(BTrial24$FitRatios, 1, FUN=mean), 
#           apply(DTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean),
#           apply(DTrial24$FitRatios, 1, FUN=mean)),
#         
#         ylim=c(0,1), col=c('darkolivegreen3','deepskyblue','orange','darkorchid'),
#         names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
# title('Fit Ratio', line=-1.4)
###############
dev.off()





#### BY F PATTERN ##########
png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\Backtransformed_ByF_#S_q.png", 
    type="cairo",
    units="mm", 
    width=170, 
    height=120, 
    pointsize=8, 
    res=600)
############
par(mfrow=c(2,1),mar=c(0.5, 2.3, 0.1, 0.1),tcl = -0.1, mgp = c(0.8, 0.2, 0))


vioplot(c(ATrial1$biasAll.BT, ATrial2$biasAll.BT, ATrial3$biasAll.BT), 
        c(ATrial4$biasAll.BT, ATrial5$biasAll.BT, ATrial6$biasAll.BT, CTrial4$biasAll.BT, CTrial5$biasAll.BT, CTrial6$biasAll.BT), 
        c(BTrial1$biasAll.BT, BTrial2$biasAll.BT, BTrial3$biasAll.BT), 
        c(BTrial4$biasAll.BT, BTrial5$biasAll.BT, BTrial6$biasAll.BT, DTrial4$biasAll.BT, DTrial5$biasAll.BT, DTrial6$biasAll.BT),
        c(ATrial7$biasAll.BT, ATrial8$biasAll.BT, ATrial9$biasAll.BT), 
        c(ATrial10$biasAll.BT, ATrial11$biasAll.BT, ATrial12$biasAll.BT, CTrial10$biasAll.BT, CTrial11$biasAll.BT, CTrial12$biasAll.BT),  
        c(BTrial7$biasAll.BT, BTrial8$biasAll.BT, BTrial9$biasAll.BT), 
        c(BTrial10$biasAll.BT, BTrial11$biasAll.BT, BTrial12$biasAll.BT, DTrial10$biasAll.BT, DTrial11$biasAll.BT, DTrial12$biasAll.BT),
        c(ATrial13$biasAll.BT, ATrial14$biasAll.BT, ATrial15$biasAll.BT), 
        c(ATrial16$biasAll.BT, ATrial17$biasAll.BT, ATrial18$biasAll.BT, CTrial16$biasAll.BT, CTrial17$biasAll.BT, CTrial18$biasAll.BT), 
        c(BTrial13$biasAll.BT, BTrial14$biasAll.BT, BTrial15$biasAll.BT), 
        c(BTrial16$biasAll.BT, BTrial17$biasAll.BT, BTrial18$biasAll.BT, DTrial16$biasAll.BT, DTrial17$biasAll.BT, DTrial18$biasAll.BT),
        c(ATrial19$biasAll.BT, ATrial20$biasAll.BT, ATrial21$biasAll.BT), 
        c(ATrial22$biasAll.BT, ATrial23$biasAll.BT, ATrial24$biasAll.BT, CTrial22$biasAll.BT, CTrial23$biasAll.BT, CTrial24$biasAll.BT),  
        c(BTrial19$biasAll.BT, BTrial20$biasAll.BT, BTrial21$biasAll.BT), 
        c(BTrial22$biasAll.BT, BTrial23$biasAll.BT, BTrial24$biasAll.BT, DTrial22$biasAll.BT, DTrial23$biasAll.BT, DTrial24$biasAll.BT),
        ylim=c(-4.5,4.5), col=c('darkolivegreen1', 'darkolivegreen3','darkolivegreen4','darkolivegreen',
                            'skyblue','deepskyblue2','deepskyblue3','deepskyblue4',
                            'orange','darkorange','darkorange2','darkorange3',
                            'mediumorchid','darkorchid1','darkorchid','darkorchid4'), names=FALSE)
mtext('Back-transformed Relative Error', 2, line=1.1, cex=1.2)
mtext('Atlantic sharpnose shark', 3, line=-1.5, cex=1.5)
abline(h=0)

par( mar=c(1.6, 2.3, 0.1, 0.1))

vioplot(c(ATrial1$RMSE.BT, ATrial2$RMSE.BT, ATrial3$RMSE.BT), 
        c(ATrial4$RMSE.BT, ATrial5$RMSE.BT, ATrial6$RMSE.BT, CTrial4$RMSE.BT, CTrial5$RMSE.BT, CTrial6$RMSE.BT), 
        c(BTrial1$RMSE.BT, BTrial2$RMSE.BT, BTrial3$RMSE.BT), 
        c(BTrial4$RMSE.BT, BTrial5$RMSE.BT, BTrial6$RMSE.BT, DTrial4$RMSE.BT, DTrial5$RMSE.BT, DTrial6$RMSE.BT), 
        c(ATrial7$RMSE.BT, ATrial8$RMSE.BT, ATrial9$RMSE.BT), 
        c(ATrial10$RMSE.BT, ATrial11$RMSE.BT, ATrial12$RMSE.BT, CTrial10$RMSE.BT, CTrial11$RMSE.BT, CTrial12$RMSE.BT),  
        c(BTrial7$RMSE.BT, BTrial8$RMSE.BT, BTrial9$RMSE.BT), 
        c(BTrial10$RMSE.BT, BTrial11$RMSE.BT, BTrial12$RMSE.BT, DTrial10$RMSE.BT, DTrial11$RMSE.BT, DTrial12$RMSE.BT),
        c(ATrial13$RMSE.BT, ATrial14$RMSE.BT, ATrial15$RMSE.BT), 
        c(ATrial16$RMSE.BT, ATrial17$RMSE.BT, ATrial18$RMSE.BT, CTrial16$RMSE.BT, CTrial17$RMSE.BT, CTrial18$RMSE.BT), 
        c(BTrial13$RMSE.BT, BTrial14$RMSE.BT, BTrial15$RMSE.BT), 
        c(BTrial16$RMSE.BT, BTrial17$RMSE.BT, BTrial18$RMSE.BT, DTrial16$RMSE.BT, DTrial17$RMSE.BT, DTrial18$RMSE.BT),
        c(ATrial19$RMSE.BT, ATrial20$RMSE.BT, ATrial21$RMSE.BT), 
        c(ATrial22$RMSE.BT, ATrial23$RMSE.BT, ATrial24$RMSE.BT, CTrial22$RMSE.BT, CTrial23$RMSE.BT, CTrial24$RMSE.BT),  
        c(BTrial19$RMSE.BT, BTrial20$RMSE.BT, BTrial21$RMSE.BT), 
        c(BTrial22$RMSE.BT, BTrial23$RMSE.BT, BTrial24$RMSE.BT, DTrial22$RMSE.BT, DTrial23$RMSE.BT, DTrial24$RMSE.BT),
        ylim=c(0,2), col=c('darkolivegreen1', 'darkolivegreen3','darkolivegreen4','darkolivegreen',
                           'skyblue','deepskyblue2','deepskyblue3','deepskyblue4',
                           'orange','darkorange','darkorange2','darkorange3',
                           'mediumorchid','darkorchid1','darkorchid','darkorchid4'), 
        names=FALSE)
        # names=c('3S', expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')),
        #         '3S',expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')),
        #         '3S',expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')),
        #         '3S',expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')) ))
axis(1, at=seq(from=1.5,to=15.5,length=8), 
     c(rep(c('3S','4S'), 4)), tck=0, cex=0.9, line=NA, pos=)
mtext('Back-transformed RMSE', 2, line=1.1, cex=1.2)

legend(12.75,2.1,c('No Change in q','Change in q','','','ConstF','IncF','DecF','UF'),pch=15,col=c('grey','dimgrey',NA,NA,'darkolivegreen3','deepskyblue2','darkorange2','darkorchid2'), bty='n', cex=1.0, pt.cex=2, ncol=2)

# legend('right',c(), pch=15, col=c('darkolivegreen3','deepskyblue2','darkorange2','darkorchid2'), pt.cex=2, bty='n', cex=1.0)

######

dev.off()


#### BY F PATTERN ##########
tiff(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\Figure_3.tiff", 
    type="cairo",
    units="mm", 
    width=170, 
    height=120, 
    pointsize=10, 
    res=600)
############
par(mfrow=c(2,1),mar=c(0.5, 2.3, 0.1, 0.1),tcl = -0.1, mgp = c(0.8, 0.2, 0))


vioplot(c(ATrial1$biasAll.BT, ATrial2$biasAll.BT, ATrial3$biasAll.BT), 
        c(ATrial4$biasAll.BT, ATrial5$biasAll.BT, ATrial6$biasAll.BT, CTrial4$biasAll.BT, CTrial5$biasAll.BT, CTrial6$biasAll.BT), 
        c(BTrial1$biasAll.BT, BTrial2$biasAll.BT, BTrial3$biasAll.BT), 
        c(BTrial4$biasAll.BT, BTrial5$biasAll.BT, BTrial6$biasAll.BT, DTrial4$biasAll.BT, DTrial5$biasAll.BT, DTrial6$biasAll.BT),
        c(ATrial7$biasAll.BT, ATrial8$biasAll.BT, ATrial9$biasAll.BT), 
        c(ATrial10$biasAll.BT, ATrial11$biasAll.BT, ATrial12$biasAll.BT, CTrial10$biasAll.BT, CTrial11$biasAll.BT, CTrial12$biasAll.BT),  
        c(BTrial7$biasAll.BT, BTrial8$biasAll.BT, BTrial9$biasAll.BT), 
        c(BTrial10$biasAll.BT, BTrial11$biasAll.BT, BTrial12$biasAll.BT, DTrial10$biasAll.BT, DTrial11$biasAll.BT, DTrial12$biasAll.BT),
        c(ATrial13$biasAll.BT, ATrial14$biasAll.BT, ATrial15$biasAll.BT), 
        c(ATrial16$biasAll.BT, ATrial17$biasAll.BT, ATrial18$biasAll.BT, CTrial16$biasAll.BT, CTrial17$biasAll.BT, CTrial18$biasAll.BT), 
        c(BTrial13$biasAll.BT, BTrial14$biasAll.BT, BTrial15$biasAll.BT), 
        c(BTrial16$biasAll.BT, BTrial17$biasAll.BT, BTrial18$biasAll.BT, DTrial16$biasAll.BT, DTrial17$biasAll.BT, DTrial18$biasAll.BT),
        c(ATrial19$biasAll.BT, ATrial20$biasAll.BT, ATrial21$biasAll.BT), 
        c(ATrial22$biasAll.BT, ATrial23$biasAll.BT, ATrial24$biasAll.BT, CTrial22$biasAll.BT, CTrial23$biasAll.BT, CTrial24$biasAll.BT),  
        c(BTrial19$biasAll.BT, BTrial20$biasAll.BT, BTrial21$biasAll.BT), 
        c(BTrial22$biasAll.BT, BTrial23$biasAll.BT, BTrial24$biasAll.BT, DTrial22$biasAll.BT, DTrial23$biasAll.BT, DTrial24$biasAll.BT),
        ylim=c(-4.5,4.5), col=c('darkolivegreen1', 'darkolivegreen3','darkolivegreen4','darkolivegreen',
                                'skyblue','deepskyblue2','deepskyblue3','deepskyblue4',
                                'orange','darkorange','darkorange2','darkorange3',
                                'mediumorchid','darkorchid1','darkorchid','darkorchid4'), names=FALSE)
mtext('Back-transformed Relative Error', 2, line=1.1, cex=1)
mtext('Atlantic sharpnose shark', 3, line=-1.5, cex=1.1)
abline(h=0)

par( mar=c(1.6, 2.3, 0.1, 0.1))

vioplot(c(ATrial1$RMSE.BT, ATrial2$RMSE.BT, ATrial3$RMSE.BT), 
        c(ATrial4$RMSE.BT, ATrial5$RMSE.BT, ATrial6$RMSE.BT, CTrial4$RMSE.BT, CTrial5$RMSE.BT, CTrial6$RMSE.BT), 
        c(BTrial1$RMSE.BT, BTrial2$RMSE.BT, BTrial3$RMSE.BT), 
        c(BTrial4$RMSE.BT, BTrial5$RMSE.BT, BTrial6$RMSE.BT, DTrial4$RMSE.BT, DTrial5$RMSE.BT, DTrial6$RMSE.BT), 
        c(ATrial7$RMSE.BT, ATrial8$RMSE.BT, ATrial9$RMSE.BT), 
        c(ATrial10$RMSE.BT, ATrial11$RMSE.BT, ATrial12$RMSE.BT, CTrial10$RMSE.BT, CTrial11$RMSE.BT, CTrial12$RMSE.BT),  
        c(BTrial7$RMSE.BT, BTrial8$RMSE.BT, BTrial9$RMSE.BT), 
        c(BTrial10$RMSE.BT, BTrial11$RMSE.BT, BTrial12$RMSE.BT, DTrial10$RMSE.BT, DTrial11$RMSE.BT, DTrial12$RMSE.BT),
        c(ATrial13$RMSE.BT, ATrial14$RMSE.BT, ATrial15$RMSE.BT), 
        c(ATrial16$RMSE.BT, ATrial17$RMSE.BT, ATrial18$RMSE.BT, CTrial16$RMSE.BT, CTrial17$RMSE.BT, CTrial18$RMSE.BT), 
        c(BTrial13$RMSE.BT, BTrial14$RMSE.BT, BTrial15$RMSE.BT), 
        c(BTrial16$RMSE.BT, BTrial17$RMSE.BT, BTrial18$RMSE.BT, DTrial16$RMSE.BT, DTrial17$RMSE.BT, DTrial18$RMSE.BT),
        c(ATrial19$RMSE.BT, ATrial20$RMSE.BT, ATrial21$RMSE.BT), 
        c(ATrial22$RMSE.BT, ATrial23$RMSE.BT, ATrial24$RMSE.BT, CTrial22$RMSE.BT, CTrial23$RMSE.BT, CTrial24$RMSE.BT),  
        c(BTrial19$RMSE.BT, BTrial20$RMSE.BT, BTrial21$RMSE.BT), 
        c(BTrial22$RMSE.BT, BTrial23$RMSE.BT, BTrial24$RMSE.BT, DTrial22$RMSE.BT, DTrial23$RMSE.BT, DTrial24$RMSE.BT),
        ylim=c(0,2), col=c('darkolivegreen1', 'darkolivegreen3','darkolivegreen4','darkolivegreen',
                           'skyblue','deepskyblue2','deepskyblue3','deepskyblue4',
                           'orange','darkorange','darkorange2','darkorange3',
                           'mediumorchid','darkorchid1','darkorchid','darkorchid4'), 
        names=c('3S', expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')),
                '3S',expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')),
                '3S',expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')),
                '3S',expression(paste('3S-',Delta,'q')),'4S', expression(paste('4S-',Delta,'q')) ))
mtext('Back-transformed RMSE', 2, line=1.1, cex=1.0)

legend('topright',c('ConstF','IncF','DecF','UF'), pch=15, col=c('darkolivegreen3','deepskyblue2','darkorange2','darkorchid2'), pt.cex=2, bty='n', cex=1.0)
######

dev.off()





#### BY F PATTERN ##########
png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\Backtransformed_ByF_#S_q_4panel.png", 
    type="cairo",
    units="mm", 
    width=170, 
    height=120, 
    pointsize=9, 
    res=600)
############
par(mfrow=c(2,2),mar=c(0.5, 2.4, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0), oma=c(0,0,1.5,0))


vioplot(c(ATrial1$biasAll.BT, ATrial2$biasAll.BT, ATrial3$biasAll.BT), 
        c(ATrial4$biasAll.BT, ATrial5$biasAll.BT, ATrial6$biasAll.BT, CTrial4$biasAll.BT, CTrial5$biasAll.BT, CTrial6$biasAll.BT), 
        
        c(ATrial7$biasAll.BT, ATrial8$biasAll.BT, ATrial9$biasAll.BT), 
        c(ATrial10$biasAll.BT, ATrial11$biasAll.BT, ATrial12$biasAll.BT, CTrial10$biasAll.BT, CTrial11$biasAll.BT, CTrial12$biasAll.BT),  
        
        c(ATrial13$biasAll.BT, ATrial14$biasAll.BT, ATrial15$biasAll.BT), 
        c(ATrial16$biasAll.BT, ATrial17$biasAll.BT, ATrial18$biasAll.BT, CTrial16$biasAll.BT, CTrial17$biasAll.BT, CTrial18$biasAll.BT), 
        
        c(ATrial19$biasAll.BT, ATrial20$biasAll.BT, ATrial21$biasAll.BT), 
        c(ATrial22$biasAll.BT, ATrial23$biasAll.BT, ATrial24$biasAll.BT, CTrial22$biasAll.BT, CTrial23$biasAll.BT, CTrial24$biasAll.BT),  
        
        
        ylim=c(-4.5,4.5), col=c('darkolivegreen1','darkolivegreen4',
                                'deepskyblue','deepskyblue4',
                                'orange','darkorange3',
                                'darkorchid1','darkorchid4'), names=FALSE, lwd=0.8)
mtext('Back-transformed Relative Error', 2, line=1.1, cex=1.2)
mtext('3 Surveys', 3, line=-1.5, cex=1.2) 

mtext('Atlantic sharpnose shark', 3, line=0, cex=1.5, outer=T)
# col=c('darkolivegreen1', 'darkolivegreen3','darkolivegreen4','darkolivegreen',
#       'skyblue','deepskyblue2','deepskyblue3','deepskyblue4',
#       'orange','darkorange','darkorange2','darkorange3',
#       'mediumorchid','darkorchid1','darkorchid','darkorchid4')
abline(h=0)

par(mar=c(0.5, 1.1, 0.1, 0.3))

vioplot(c(BTrial1$biasAll.BT, BTrial2$biasAll.BT, BTrial3$biasAll.BT), 
        c(BTrial4$biasAll.BT, BTrial5$biasAll.BT, BTrial6$biasAll.BT, DTrial4$biasAll.BT, DTrial5$biasAll.BT, DTrial6$biasAll.BT),
        
        c(BTrial7$biasAll.BT, BTrial8$biasAll.BT, BTrial9$biasAll.BT), 
        c(BTrial10$biasAll.BT, BTrial11$biasAll.BT, BTrial12$biasAll.BT, DTrial10$biasAll.BT, DTrial11$biasAll.BT, DTrial12$biasAll.BT),
        
        c(BTrial13$biasAll.BT, BTrial14$biasAll.BT, BTrial15$biasAll.BT), 
        c(BTrial16$biasAll.BT, BTrial17$biasAll.BT, BTrial18$biasAll.BT, DTrial16$biasAll.BT, DTrial17$biasAll.BT, DTrial18$biasAll.BT),
        
        c(BTrial19$biasAll.BT, BTrial20$biasAll.BT, BTrial21$biasAll.BT), 
        c(BTrial22$biasAll.BT, BTrial23$biasAll.BT, BTrial24$biasAll.BT, DTrial22$biasAll.BT, DTrial23$biasAll.BT, DTrial24$biasAll.BT),
        ylim=c(-4.5,4.5),  col=c('darkolivegreen1','darkolivegreen4',
                                 'deepskyblue','deepskyblue4',
                                 'orange','darkorange3',
                                 'darkorchid1','darkorchid4'), names=FALSE, lwd=0.8)
mtext('4 Surveys', 3, line=-1.5, cex=1.2)
abline(h=0)

par( mar=c(1.6, 2.4, 0.1, 0.1))

vioplot(c(ATrial1$RMSE.BT, ATrial2$RMSE.BT, ATrial3$RMSE.BT), 
        c(ATrial4$RMSE.BT, ATrial5$RMSE.BT, ATrial6$RMSE.BT, CTrial4$RMSE.BT, CTrial5$RMSE.BT, CTrial6$RMSE.BT), 
        
        c(ATrial7$RMSE.BT, ATrial8$RMSE.BT, ATrial9$RMSE.BT), 
        c(ATrial10$RMSE.BT, ATrial11$RMSE.BT, ATrial12$RMSE.BT, CTrial10$RMSE.BT, CTrial11$RMSE.BT, CTrial12$RMSE.BT),  
        
        c(ATrial13$RMSE.BT, ATrial14$RMSE.BT, ATrial15$RMSE.BT), 
        c(ATrial16$RMSE.BT, ATrial17$RMSE.BT, ATrial18$RMSE.BT, CTrial16$RMSE.BT, CTrial17$RMSE.BT, CTrial18$RMSE.BT), 
        
        c(ATrial19$RMSE.BT, ATrial20$RMSE.BT, ATrial21$RMSE.BT), 
        c(ATrial22$RMSE.BT, ATrial23$RMSE.BT, ATrial24$RMSE.BT, CTrial22$RMSE.BT, CTrial23$RMSE.BT, CTrial24$RMSE.BT),  
        
        ylim=c(0,2), col=c('darkolivegreen1','darkolivegreen4',
                           'deepskyblue','deepskyblue4',
                           'orange','darkorange3',
                           'darkorchid1','darkorchid4'), 
        names=FALSE , lwd=0.8)
axis(1, at=seq(1.5,7.5, length=4), 
     c('Const F','Inc F','Dec F','U F'), tck=0, cex=0.9)
mtext('Back-transformed RMSE', 2, line=1.1, cex=1.2)


par( mar=c(1.6, 1.1, 0.1, 0.1))

vioplot(c(BTrial1$RMSE.BT, BTrial2$RMSE.BT, BTrial3$RMSE.BT), 
        c(BTrial4$RMSE.BT, BTrial5$RMSE.BT, BTrial6$RMSE.BT, DTrial4$RMSE.BT, DTrial5$RMSE.BT, DTrial6$RMSE.BT), 
        
        c(BTrial7$RMSE.BT, BTrial8$RMSE.BT, BTrial9$RMSE.BT), 
        c(BTrial10$RMSE.BT, BTrial11$RMSE.BT, BTrial12$RMSE.BT, DTrial10$RMSE.BT, DTrial11$RMSE.BT, DTrial12$RMSE.BT),
        
        c(BTrial13$RMSE.BT, BTrial14$RMSE.BT, BTrial15$RMSE.BT), 
        c(BTrial16$RMSE.BT, BTrial17$RMSE.BT, BTrial18$RMSE.BT, DTrial16$RMSE.BT, DTrial17$RMSE.BT, DTrial18$RMSE.BT),
        
        c(BTrial19$RMSE.BT, BTrial20$RMSE.BT, BTrial21$RMSE.BT), 
        c(BTrial22$RMSE.BT, BTrial23$RMSE.BT, BTrial24$RMSE.BT, DTrial22$RMSE.BT, DTrial23$RMSE.BT, DTrial24$RMSE.BT),
        ylim=c(0,2), col=c('darkolivegreen1','darkolivegreen4',
                           'deepskyblue','deepskyblue4',
                           'orange','darkorange3',
                           'darkorchid1','darkorchid4'), 
        names=FALSE, lwd=0.8)

axis(1, at=seq(1.5,7.5, length=4), 
     c('Const F','Inc F','Dec F','U F'), tck=0, cex=0.9)

legend('topright',c('No Change in q','Change in q'),pch=15,col=c('grey','dimgrey'), bty='n', cex=1.0, pt.cex=2)
######

dev.off()




#### BY F PATTERN ##########
png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\Backtransformed_Just3UF.png", 
    type="cairo",
    units="mm", 
    width=85, 
    height=80, 
    pointsize=6, 
    res=600)
############
par(mfrow=c(2,1),mar=c(0.5, 2.4, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0), oma=c(1,0,0,0))


vioplot(
        c(ATrial19$biasAll.BT), c(ATrial20$biasAll.BT), c( ATrial21$biasAll.BT), 
        c(ATrial22$biasAll.BT), c(ATrial23$biasAll.BT), c(ATrial24$biasAll.BT), 
        c(CTrial22$biasAll.BT), c(CTrial23$biasAll.BT), c(CTrial24$biasAll.BT),  
        
        
        ylim=c(-4,3), col=c('orchid','orchid','orchid',
                                'darkorchid2','darkorchid2','darkorchid2',
                                'darkorchid4','darkorchid4','darkorchid4'), names=FALSE, lwd=0.8)

mtext('Back-transformed Relative Error', 2, line=1.1, cex=1.2)
# 'darkolivegreen1','darkolivegreen4',
# 'deepskyblue','deepskyblue4',
# 'orange','darkorange3',
# 'darkorchid1','darkorchid4'
mtext('Atlantic sharpnose shark', 3, line=-1.5, cex=1.5)
mtext('3 Surveys / UF', 3, line=-2.5, cex=1.2) 

abline(h=0)

par( mar=c(1.6, 2.4, 0.1, 0.1))


vioplot(
        c(ATrial19$RMSE.BT), c(ATrial20$RMSE.BT), c(ATrial21$RMSE.BT), 
        
        c(ATrial22$RMSE.BT), c(ATrial23$RMSE.BT), c(ATrial24$RMSE.BT), 
        c(CTrial22$RMSE.BT), c(CTrial23$RMSE.BT), c(CTrial24$RMSE.BT),  
        
        ylim=c(0,1.6), col=c('orchid','orchid','orchid',
                           'darkorchid2','darkorchid2','darkorchid2',
                           'darkorchid4','darkorchid4','darkorchid4'), 
        names=FALSE , lwd=0.8)
axis(1, labels=FALSE, tck=0)
text(seq(1.4, 9.4, by=1), par("usr")[3] - 0.03, 
     labels = c('UF_const_1','UF_const_2','UF_const_3',
                'UF_knife_1','UF_knife_2','UF_knife_3',
                'UF_grad_1','UF_grad_2','UF_grad_3'), 
     srt = 25, xpd = NA, adj=1, cex=1.0)

mtext('Back-transformed RMSE', 2, line=1.1, cex=1.2)

# legend('topright',c('No Change in q','Change in q'),pch=15,col=c('grey','dimgrey'), bty='n', cex=1.0, pt.cex=2)
######

dev.off()




#### BY F PATTERN ##########
tiff(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\Backtransformed_ByF_#S_q_4panel.tiff", 
    type="cairo",
    units="mm", 
    width=170, 
    height=120, 
    pointsize=9, 
    res=600)
############
par(mfrow=c(2,2),mar=c(0.5, 2.4, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0), oma=c(0,0,1.5,0))


vioplot(c(ATrial1$biasAll.BT, ATrial2$biasAll.BT, ATrial3$biasAll.BT), 
        c(ATrial4$biasAll.BT, ATrial5$biasAll.BT, ATrial6$biasAll.BT, CTrial4$biasAll.BT, CTrial5$biasAll.BT, CTrial6$biasAll.BT), 
        
        c(ATrial7$biasAll.BT, ATrial8$biasAll.BT, ATrial9$biasAll.BT), 
        c(ATrial10$biasAll.BT, ATrial11$biasAll.BT, ATrial12$biasAll.BT, CTrial10$biasAll.BT, CTrial11$biasAll.BT, CTrial12$biasAll.BT),  
        
        c(ATrial13$biasAll.BT, ATrial14$biasAll.BT, ATrial15$biasAll.BT), 
        c(ATrial16$biasAll.BT, ATrial17$biasAll.BT, ATrial18$biasAll.BT, CTrial16$biasAll.BT, CTrial17$biasAll.BT, CTrial18$biasAll.BT), 
        
        c(ATrial19$biasAll.BT, ATrial20$biasAll.BT, ATrial21$biasAll.BT), 
        c(ATrial22$biasAll.BT, ATrial23$biasAll.BT, ATrial24$biasAll.BT, CTrial22$biasAll.BT, CTrial23$biasAll.BT, CTrial24$biasAll.BT),  
        
        
        ylim=c(-4.5,4.5), col=c('darkolivegreen1','darkolivegreen4',
                                'deepskyblue','deepskyblue4',
                                'orange','darkorange3',
                                'darkorchid1','darkorchid4'), names=FALSE, lwd=0.8)
mtext('Back-transformed Relative Error', 2, line=1.1, cex=1.2)
mtext('3 Surveys', 3, line=-1.5, cex=1.2) 

mtext('Atlantic sharpnose shark', 3, line=0, cex=1.5, outer=T)
# col=c('darkolivegreen1', 'darkolivegreen3','darkolivegreen4','darkolivegreen',
#       'skyblue','deepskyblue2','deepskyblue3','deepskyblue4',
#       'orange','darkorange','darkorange2','darkorange3',
#       'mediumorchid','darkorchid1','darkorchid','darkorchid4')
abline(h=0)

par(mar=c(0.5, 1.1, 0.1, 0.3))

vioplot(c(BTrial1$biasAll.BT, BTrial2$biasAll.BT, BTrial3$biasAll.BT), 
        c(BTrial4$biasAll.BT, BTrial5$biasAll.BT, BTrial6$biasAll.BT, DTrial4$biasAll.BT, DTrial5$biasAll.BT, DTrial6$biasAll.BT),
        
        c(BTrial7$biasAll.BT, BTrial8$biasAll.BT, BTrial9$biasAll.BT), 
        c(BTrial10$biasAll.BT, BTrial11$biasAll.BT, BTrial12$biasAll.BT, DTrial10$biasAll.BT, DTrial11$biasAll.BT, DTrial12$biasAll.BT),
        
        c(BTrial13$biasAll.BT, BTrial14$biasAll.BT, BTrial15$biasAll.BT), 
        c(BTrial16$biasAll.BT, BTrial17$biasAll.BT, BTrial18$biasAll.BT, DTrial16$biasAll.BT, DTrial17$biasAll.BT, DTrial18$biasAll.BT),
        
        c(BTrial19$biasAll.BT, BTrial20$biasAll.BT, BTrial21$biasAll.BT), 
        c(BTrial22$biasAll.BT, BTrial23$biasAll.BT, BTrial24$biasAll.BT, DTrial22$biasAll.BT, DTrial23$biasAll.BT, DTrial24$biasAll.BT),
        ylim=c(-4.5,4.5),  col=c('darkolivegreen1','darkolivegreen4',
                                 'deepskyblue','deepskyblue4',
                                 'orange','darkorange3',
                                 'darkorchid1','darkorchid4'), names=FALSE, lwd=0.8)
mtext('4 Surveys', 3, line=-1.5, cex=1.2)
abline(h=0)

par( mar=c(1.6, 2.4, 0.1, 0.1))

vioplot(c(ATrial1$RMSE.BT, ATrial2$RMSE.BT, ATrial3$RMSE.BT), 
        c(ATrial4$RMSE.BT, ATrial5$RMSE.BT, ATrial6$RMSE.BT, CTrial4$RMSE.BT, CTrial5$RMSE.BT, CTrial6$RMSE.BT), 
        
        c(ATrial7$RMSE.BT, ATrial8$RMSE.BT, ATrial9$RMSE.BT), 
        c(ATrial10$RMSE.BT, ATrial11$RMSE.BT, ATrial12$RMSE.BT, CTrial10$RMSE.BT, CTrial11$RMSE.BT, CTrial12$RMSE.BT),  
        
        c(ATrial13$RMSE.BT, ATrial14$RMSE.BT, ATrial15$RMSE.BT), 
        c(ATrial16$RMSE.BT, ATrial17$RMSE.BT, ATrial18$RMSE.BT, CTrial16$RMSE.BT, CTrial17$RMSE.BT, CTrial18$RMSE.BT), 
        
        c(ATrial19$RMSE.BT, ATrial20$RMSE.BT, ATrial21$RMSE.BT), 
        c(ATrial22$RMSE.BT, ATrial23$RMSE.BT, ATrial24$RMSE.BT, CTrial22$RMSE.BT, CTrial23$RMSE.BT, CTrial24$RMSE.BT),  
        
        ylim=c(0,2), col=c('darkolivegreen1','darkolivegreen4',
                           'deepskyblue','deepskyblue4',
                           'orange','darkorange3',
                           'darkorchid1','darkorchid4'), 
        names=FALSE , lwd=0.8)
axis(1, at=seq(1.5,7.5, length=4), 
     c('Const F','Inc F','Dec F','U F'), tck=0, cex=0.9)
mtext('Back-transformed RMSE', 2, line=1.1, cex=1.2)


par( mar=c(1.6, 1.1, 0.1, 0.1))

vioplot(c(BTrial1$RMSE.BT, BTrial2$RMSE.BT, BTrial3$RMSE.BT), 
        c(BTrial4$RMSE.BT, BTrial5$RMSE.BT, BTrial6$RMSE.BT, DTrial4$RMSE.BT, DTrial5$RMSE.BT, DTrial6$RMSE.BT), 
        
        c(BTrial7$RMSE.BT, BTrial8$RMSE.BT, BTrial9$RMSE.BT), 
        c(BTrial10$RMSE.BT, BTrial11$RMSE.BT, BTrial12$RMSE.BT, DTrial10$RMSE.BT, DTrial11$RMSE.BT, DTrial12$RMSE.BT),
        
        c(BTrial13$RMSE.BT, BTrial14$RMSE.BT, BTrial15$RMSE.BT), 
        c(BTrial16$RMSE.BT, BTrial17$RMSE.BT, BTrial18$RMSE.BT, DTrial16$RMSE.BT, DTrial17$RMSE.BT, DTrial18$RMSE.BT),
        
        c(BTrial19$RMSE.BT, BTrial20$RMSE.BT, BTrial21$RMSE.BT), 
        c(BTrial22$RMSE.BT, BTrial23$RMSE.BT, BTrial24$RMSE.BT, DTrial22$RMSE.BT, DTrial23$RMSE.BT, DTrial24$RMSE.BT),
        ylim=c(0,2), col=c('darkolivegreen1','darkolivegreen4',
                           'deepskyblue','deepskyblue4',
                           'orange','darkorange3',
                           'darkorchid1','darkorchid4'), 
        names=FALSE, lwd=0.8)

axis(1, at=seq(1.5,7.5, length=4), 
     c('Const F','Inc F','Dec F','U F'), tck=0, cex=0.9)

legend('topright',c('No Change in q','Change in q'),pch=15,col=c('grey','dimgrey'), bty='n', cex=1.0, pt.cex=2)
######

dev.off()




#### BY F PATTERN ##########
tiff(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\Backtransformed_Just3UF.tiff", 
    type="cairo",
    units="mm", 
    width=85, 
    height=80, 
    pointsize=6, 
    res=600)
############
par(mfrow=c(2,1),mar=c(0.5, 2.4, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0), oma=c(1,0,0,0))


vioplot(
  c(ATrial19$biasAll.BT), c(ATrial20$biasAll.BT), c( ATrial21$biasAll.BT), 
  c(ATrial22$biasAll.BT), c(ATrial23$biasAll.BT), c(ATrial24$biasAll.BT), 
  c(CTrial22$biasAll.BT), c(CTrial23$biasAll.BT), c(CTrial24$biasAll.BT),  
  
  
  ylim=c(-4,3), col=c('orchid','orchid','orchid',
                      'darkorchid2','darkorchid2','darkorchid2',
                      'darkorchid4','darkorchid4','darkorchid4'), names=FALSE, lwd=0.8)

mtext('Back-transformed Relative Error', 2, line=1.1, cex=1.2)
# 'darkolivegreen1','darkolivegreen4',
# 'deepskyblue','deepskyblue4',
# 'orange','darkorange3',
# 'darkorchid1','darkorchid4'
mtext('Atlantic sharpnose shark', 3, line=-1.5, cex=1.5)
mtext('3 Surveys / UF', 3, line=-2.5, cex=1.2) 

abline(h=0)

par( mar=c(1.6, 2.4, 0.1, 0.1))


vioplot(
  c(ATrial19$RMSE.BT), c(ATrial20$RMSE.BT), c(ATrial21$RMSE.BT), 
  
  c(ATrial22$RMSE.BT), c(ATrial23$RMSE.BT), c(ATrial24$RMSE.BT), 
  c(CTrial22$RMSE.BT), c(CTrial23$RMSE.BT), c(CTrial24$RMSE.BT),  
  
  ylim=c(0,1.6), col=c('orchid','orchid','orchid',
                       'darkorchid2','darkorchid2','darkorchid2',
                       'darkorchid4','darkorchid4','darkorchid4'), 
  names=FALSE , lwd=0.8)
axis(1, labels=FALSE, tck=0)
text(seq(1.4, 9.4, by=1), par("usr")[3] - 0.03, 
     labels = c('UF_const_1','UF_const_2','UF_const_3',
                'UF_knife_1','UF_knife_2','UF_knife_3',
                'UF_grad_1','UF_grad_2','UF_grad_3'), 
     srt = 25, xpd = NA, adj=1, cex=1.0)

mtext('Back-transformed RMSE', 2, line=1.1, cex=1.2)

# legend('topright',c('No Change in q','Change in q'),pch=15,col=c('grey','dimgrey'), bty='n', cex=1.0, pt.cex=2)
######

dev.off()





png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\AllBiasRMSE_FitRatios.png", 
     type="cairo",
     units="mm", 
     width=500, 
     height=167, 
     pointsize=16, 
     res=600)
######
par(mfrow=c(1,3), mar=c(1.6, 1.6, 1.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
# 
# vioplot(c(ATrial1$AbsBias, ATrial2$AbsBias, ATrial3$AbsBias, ATrial4$AbsBias, ATrial5$AbsBias, ATrial6$AbsBias, 
#           CTrial4$AbsBias, CTrial5$AbsBias, CTrial6$AbsBias, 
#           BTrial1$AbsBias, BTrial2$AbsBias, BTrial3$AbsBias, BTrial4$AbsBias, BTrial5$AbsBias, BTrial6$AbsBias,
#           DTrial4$AbsBias, DTrial5$AbsBias, DTrial6$AbsBias),
#         c(ATrial7$AbsBias, ATrial8$AbsBias, ATrial9$AbsBias, ATrial10$AbsBias, ATrial11$AbsBias, ATrial12$AbsBias, 
#           CTrial10$AbsBias, CTrial11$AbsBias, CTrial12$AbsBias,  
#           BTrial7$AbsBias, BTrial8$AbsBias, BTrial9$AbsBias, BTrial10$AbsBias, BTrial11$AbsBias, BTrial12$AbsBias, 
#           DTrial10$AbsBias, DTrial11$AbsBias, DTrial12$AbsBias),
#         c(ATrial13$AbsBias, ATrial14$AbsBias, ATrial15$AbsBias, ATrial16$AbsBias, ATrial17$AbsBias, ATrial18$AbsBias, 
#           CTrial16$AbsBias, CTrial17$AbsBias, CTrial18$AbsBias, 
#           BTrial13$AbsBias, BTrial14$AbsBias, BTrial15$AbsBias, BTrial16$AbsBias, BTrial17$AbsBias, BTrial18$AbsBias, 
#           DTrial16$AbsBias, DTrial17$AbsBias, DTrial18$AbsBias),
#         c(ATrial19$AbsBias, ATrial20$AbsBias, ATrial21$AbsBias, ATrial22$AbsBias, ATrial23$AbsBias, ATrial24$AbsBias, 
#           CTrial22$AbsBias, CTrial23$AbsBias, CTrial24$AbsBias,  
#           BTrial19$AbsBias, BTrial20$AbsBias, BTrial21$AbsBias, BTrial22$AbsBias, BTrial23$AbsBias, BTrial24$AbsBias, 
#           DTrial22$AbsBias, DTrial23$AbsBias, DTrial24$AbsBias),
#         ylim=c(0,50), col=c(rep('deepskyblue',2),'lightskyblue','darkolivegreen3','darkolivegreen1'),
#         names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
# title('Relative Error', line=-1.4)

vioplot(c(ATrial1$biasAll, ATrial2$biasAll, ATrial3$biasAll, ATrial4$biasAll, ATrial5$biasAll, ATrial6$biasAll, 
          CTrial4$biasAll, CTrial5$biasAll, CTrial6$biasAll, 
          BTrial1$biasAll, BTrial2$biasAll, BTrial3$biasAll, BTrial4$biasAll, BTrial5$biasAll, BTrial6$biasAll,
          DTrial4$biasAll, DTrial5$biasAll, DTrial6$biasAll),
        c(ATrial7$biasAll, ATrial8$biasAll, ATrial9$biasAll, ATrial10$biasAll, ATrial11$biasAll, ATrial12$biasAll, 
          CTrial10$biasAll, CTrial11$biasAll, CTrial12$biasAll,  
          BTrial7$biasAll, BTrial8$biasAll, BTrial9$biasAll, BTrial10$biasAll, BTrial11$biasAll, BTrial12$biasAll, 
          DTrial10$biasAll, DTrial11$biasAll, DTrial12$biasAll),
        c(ATrial13$biasAll, ATrial14$biasAll, ATrial15$biasAll, ATrial16$biasAll, ATrial17$biasAll, ATrial18$biasAll, 
          CTrial16$biasAll, CTrial17$biasAll, CTrial18$biasAll, 
          BTrial13$biasAll, BTrial14$biasAll, BTrial15$biasAll, BTrial16$biasAll, BTrial17$biasAll, BTrial18$biasAll, 
          DTrial16$biasAll, DTrial17$biasAll, DTrial18$biasAll),
        c(ATrial19$biasAll, ATrial20$biasAll, ATrial21$biasAll, ATrial22$biasAll, ATrial23$biasAll, ATrial24$biasAll, 
          CTrial22$biasAll, CTrial23$biasAll, CTrial24$biasAll,  
          BTrial19$biasAll, BTrial20$biasAll, BTrial21$biasAll, BTrial22$biasAll, BTrial23$biasAll, BTrial24$biasAll, 
          DTrial22$biasAll, DTrial23$biasAll, DTrial24$biasAll),
        ylim=c(-4.5,4.5), col=c('deepskyblue','lightskyblue','darkolivegreen3','darkolivegreen1'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
title('Relative Error', line=-1.4)
# abline(h=0)

vioplot(c(ATrial1$RMSE, ATrial2$RMSE, ATrial3$RMSE, ATrial4$RMSE, ATrial5$RMSE, ATrial6$RMSE, 
          CTrial4$RMSE, CTrial5$RMSE, CTrial6$RMSE, 
          BTrial1$RMSE, BTrial2$RMSE, BTrial3$RMSE, BTrial4$RMSE, BTrial5$RMSE, BTrial6$RMSE,
          DTrial4$RMSE, DTrial5$RMSE, DTrial6$RMSE), 
        c(ATrial7$RMSE, ATrial8$RMSE, ATrial9$RMSE, ATrial10$RMSE, ATrial11$RMSE, ATrial12$RMSE, 
          CTrial10$RMSE, CTrial11$RMSE, CTrial12$RMSE,  
          BTrial7$RMSE, BTrial8$RMSE, BTrial9$RMSE, BTrial10$RMSE, BTrial11$RMSE, BTrial12$RMSE, 
          DTrial10$RMSE, DTrial11$RMSE, DTrial12$RMSE),
        c(ATrial13$RMSE, ATrial14$RMSE, ATrial15$RMSE, ATrial16$RMSE, ATrial17$RMSE, ATrial18$RMSE, 
          CTrial16$RMSE, CTrial17$RMSE, CTrial18$RMSE, 
          BTrial13$RMSE, BTrial14$RMSE, BTrial15$RMSE, BTrial16$RMSE, BTrial17$RMSE, BTrial18$RMSE, 
          DTrial16$RMSE, DTrial17$RMSE, DTrial18$RMSE),
        c(ATrial19$RMSE, ATrial20$RMSE, ATrial21$RMSE, ATrial22$RMSE, ATrial23$RMSE, ATrial24$RMSE, 
          CTrial22$RMSE, CTrial23$RMSE, CTrial24$RMSE,  
          BTrial19$RMSE, BTrial20$RMSE, BTrial21$RMSE, BTrial22$RMSE, BTrial23$RMSE, BTrial24$RMSE, 
          DTrial22$RMSE, DTrial23$RMSE, DTrial24$RMSE),
        ylim=c(0,2), col=c('deepskyblue','lightskyblue','darkolivegreen3','darkolivegreen1'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
title('RMSE', line=-1.4)



vioplot(c(apply(ATrial1$FitRatios, 1, FUN=mean), apply(ATrial2$FitRatios, 1, FUN=mean), 
          apply(ATrial3$FitRatios, 1, FUN=mean), apply(ATrial4$FitRatios, 1, FUN=mean),
          apply(ATrial5$FitRatios, 1, FUN=mean), apply(ATrial6$FitRatios, 1, FUN=mean), 
          apply(CTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean),
          apply(CTrial6$FitRatios, 1, FUN=mean),
          apply(BTrial1$FitRatios, 1, FUN=mean), apply(BTrial2$FitRatios, 1, FUN=mean), 
          apply(BTrial3$FitRatios, 1, FUN=mean), apply(BTrial4$FitRatios, 1, FUN=mean),
          apply(BTrial5$FitRatios, 1, FUN=mean), apply(BTrial6$FitRatios, 1, FUN=mean), 
          apply(DTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean),
          apply(DTrial6$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial7$FitRatios, 1, FUN=mean), apply(ATrial8$FitRatios, 1, FUN=mean), 
          apply(ATrial9$FitRatios, 1, FUN=mean), apply(ATrial10$FitRatios, 1, FUN=mean),
          apply(ATrial11$FitRatios, 1, FUN=mean), apply(ATrial12$FitRatios, 1, FUN=mean), 
          apply(CTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean),
          apply(CTrial12$FitRatios, 1, FUN=mean),
          apply(BTrial7$FitRatios, 1, FUN=mean), apply(BTrial8$FitRatios, 1, FUN=mean), 
          apply(BTrial9$FitRatios, 1, FUN=mean), apply(BTrial10$FitRatios, 1, FUN=mean),
          apply(BTrial11$FitRatios, 1, FUN=mean), apply(BTrial12$FitRatios, 1, FUN=mean), 
          apply(DTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean),
          apply(DTrial12$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial13$FitRatios, 1, FUN=mean), apply(ATrial14$FitRatios, 1, FUN=mean), 
          apply(ATrial15$FitRatios, 1, FUN=mean), apply(ATrial16$FitRatios, 1, FUN=mean),
          apply(ATrial17$FitRatios, 1, FUN=mean), apply(ATrial18$FitRatios, 1, FUN=mean), 
          apply(CTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean),
          apply(CTrial18$FitRatios, 1, FUN=mean),
          apply(BTrial13$FitRatios, 1, FUN=mean), apply(BTrial14$FitRatios, 1, FUN=mean), 
          apply(BTrial15$FitRatios, 1, FUN=mean), apply(BTrial16$FitRatios, 1, FUN=mean),
          apply(BTrial17$FitRatios, 1, FUN=mean), apply(BTrial18$FitRatios, 1, FUN=mean), 
          apply(DTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean),
          apply(DTrial18$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial19$FitRatios, 1, FUN=mean), apply(ATrial20$FitRatios, 1, FUN=mean), 
          apply(ATrial21$FitRatios, 1, FUN=mean), apply(ATrial22$FitRatios, 1, FUN=mean),
          apply(ATrial23$FitRatios, 1, FUN=mean), apply(ATrial24$FitRatios, 1, FUN=mean), 
          apply(CTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean),
          apply(CTrial24$FitRatios, 1, FUN=mean),
          apply(BTrial19$FitRatios, 1, FUN=mean), apply(BTrial20$FitRatios, 1, FUN=mean), 
          apply(BTrial21$FitRatios, 1, FUN=mean), apply(BTrial22$FitRatios, 1, FUN=mean),
          apply(BTrial23$FitRatios, 1, FUN=mean), apply(BTrial24$FitRatios, 1, FUN=mean), 
          apply(DTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean),
          apply(DTrial24$FitRatios, 1, FUN=mean)),
        
        ylim=c(0,0.8), col=c('deepskyblue','lightskyblue','darkolivegreen3','darkolivegreen1'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
title('Fit Ratio', line=-1.4)
############
dev.off()

par(mfrow=c(1,1))
vioplot(apply(ATrial1$FitRatios, 1, FUN=mean), ylim=c(0,1), col='deepskyblue', names=c('Trial Base 1'))



png(filename="D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\FitRatios.png", 
     type="cairo",
     units="mm", 
     width=250, 
     height=250, 
     pointsize=16, 
     res=600)
#######
par(mfrow=c(1,1), mar=c(1.6, 1.6, 1.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))


vioplot(c(apply(ATrial1$FitRatios, 1, FUN=mean), apply(ATrial2$FitRatios, 1, FUN=mean), 
          apply(ATrial3$FitRatios, 1, FUN=mean), apply(ATrial4$FitRatios, 1, FUN=mean),
          apply(ATrial5$FitRatios, 1, FUN=mean), apply(ATrial6$FitRatios, 1, FUN=mean), 
          apply(CTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean),
          apply(CTrial6$FitRatios, 1, FUN=mean),
          apply(BTrial1$FitRatios, 1, FUN=mean), apply(BTrial2$FitRatios, 1, FUN=mean), 
          apply(BTrial3$FitRatios, 1, FUN=mean), apply(BTrial4$FitRatios, 1, FUN=mean),
          apply(BTrial5$FitRatios, 1, FUN=mean), apply(BTrial6$FitRatios, 1, FUN=mean), 
          apply(DTrial4$FitRatios, 1, FUN=mean), apply(CTrial5$FitRatios, 1, FUN=mean),
          apply(DTrial6$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial7$FitRatios, 1, FUN=mean), apply(ATrial8$FitRatios, 1, FUN=mean), 
          apply(ATrial9$FitRatios, 1, FUN=mean), apply(ATrial10$FitRatios, 1, FUN=mean),
          apply(ATrial11$FitRatios, 1, FUN=mean), apply(ATrial12$FitRatios, 1, FUN=mean), 
          apply(CTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean),
          apply(CTrial12$FitRatios, 1, FUN=mean),
          apply(BTrial7$FitRatios, 1, FUN=mean), apply(BTrial8$FitRatios, 1, FUN=mean), 
          apply(BTrial9$FitRatios, 1, FUN=mean), apply(BTrial10$FitRatios, 1, FUN=mean),
          apply(BTrial11$FitRatios, 1, FUN=mean), apply(BTrial12$FitRatios, 1, FUN=mean), 
          apply(DTrial10$FitRatios, 1, FUN=mean), apply(CTrial11$FitRatios, 1, FUN=mean),
          apply(DTrial12$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial13$FitRatios, 1, FUN=mean), apply(ATrial14$FitRatios, 1, FUN=mean), 
          apply(ATrial15$FitRatios, 1, FUN=mean), apply(ATrial16$FitRatios, 1, FUN=mean),
          apply(ATrial17$FitRatios, 1, FUN=mean), apply(ATrial18$FitRatios, 1, FUN=mean), 
          apply(CTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean),
          apply(CTrial18$FitRatios, 1, FUN=mean),
          apply(BTrial13$FitRatios, 1, FUN=mean), apply(BTrial14$FitRatios, 1, FUN=mean), 
          apply(BTrial15$FitRatios, 1, FUN=mean), apply(BTrial16$FitRatios, 1, FUN=mean),
          apply(BTrial17$FitRatios, 1, FUN=mean), apply(BTrial18$FitRatios, 1, FUN=mean), 
          apply(DTrial16$FitRatios, 1, FUN=mean), apply(CTrial17$FitRatios, 1, FUN=mean),
          apply(DTrial18$FitRatios, 1, FUN=mean)),
        
        c(apply(ATrial19$FitRatios, 1, FUN=mean), apply(ATrial20$FitRatios, 1, FUN=mean), 
          apply(ATrial21$FitRatios, 1, FUN=mean), apply(ATrial22$FitRatios, 1, FUN=mean),
          apply(ATrial23$FitRatios, 1, FUN=mean), apply(ATrial24$FitRatios, 1, FUN=mean), 
          apply(CTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean),
          apply(CTrial24$FitRatios, 1, FUN=mean),
          apply(BTrial19$FitRatios, 1, FUN=mean), apply(BTrial20$FitRatios, 1, FUN=mean), 
          apply(BTrial21$FitRatios, 1, FUN=mean), apply(BTrial22$FitRatios, 1, FUN=mean),
          apply(BTrial23$FitRatios, 1, FUN=mean), apply(BTrial24$FitRatios, 1, FUN=mean), 
          apply(DTrial22$FitRatios, 1, FUN=mean), apply(CTrial23$FitRatios, 1, FUN=mean),
          apply(DTrial24$FitRatios, 1, FUN=mean)),
        
        ylim=c(0,0.8), col=c('deepskyblue','lightskyblue','darkolivegreen3','darkolivegreen1'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
title('Fit Ratio', line=-1.4)
#######
dev.off()


##### SAME SCALE ############
setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\Plots\\MS1\\")

png(filename="SN_BiasRMSE.png", 
    type="cairo",
    units="mm", 
    width=300, 
    height=375, 
    pointsize=16, 
    res=600)
###########
par(mfrow=c(4,2), mar=c(0.6,1.6,1.1,0), tcl=-0.25, mgp=c(2, 0.4, 0))
vioplot(c(ATrial1$biasAll), c(ATrial2$biasAll), c( ATrial3$biasAll), c( ATrial4$biasAll), c( ATrial5$biasAll), c( ATrial6$biasAll), 
        c(CTrial4$biasAll), c( CTrial5$biasAll), c( CTrial6$biasAll), 
        c(BTrial1$biasAll), c( BTrial2$biasAll), c( BTrial3$biasAll), c( BTrial4$biasAll), c( BTrial5$biasAll), c( BTrial6$biasAll),
        c(DTrial4$biasAll), c( DTrial5$biasAll), c( DTrial6$biasAll),
        ylim=c(-4.5,4.5),
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('No Trend', line=-1)
title('Relative Error', line=0.4, cex=1.5)


vioplot(ATrial1$RMSE, ATrial2$RMSE, ATrial3$RMSE, ATrial4$RMSE, ATrial5$RMSE, ATrial6$RMSE, 
        CTrial4$RMSE, CTrial5$RMSE, CTrial6$RMSE, 
        BTrial1$RMSE, BTrial2$RMSE, BTrial3$RMSE, BTrial4$RMSE, BTrial5$RMSE, BTrial6$RMSE,
        DTrial4$RMSE, DTrial5$RMSE, DTrial6$RMSE,
        ylim=c(0,2), 
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('No Trend', line=-1)
title('RMSE', line=0.4, cex=1.5)



par(mar=c(0.6,1.6,0,0), tcl=-0.25, mgp=c(2, 0.4, 0))
vioplot(c(ATrial7$biasAll), c( ATrial8$biasAll), c( ATrial9$biasAll), c( ATrial10$biasAll), c( ATrial11$biasAll), c( ATrial12$biasAll), 
        c(CTrial10$biasAll), c( CTrial11$biasAll), c( CTrial12$biasAll),  
        c(BTrial7$biasAll), c( BTrial8$biasAll), c( BTrial9$biasAll), c( BTrial10$biasAll), c( BTrial11$biasAll), c( BTrial12$biasAll), 
        c(DTrial10$biasAll), c( DTrial11$biasAll), c( DTrial12$biasAll),
        ylim=c(-4.5,4.5),
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Decreasing Trend', line=-1)

vioplot(ATrial7$RMSE, ATrial8$RMSE, ATrial9$RMSE, ATrial10$RMSE, ATrial11$RMSE, ATrial12$RMSE, 
        CTrial10$RMSE, CTrial11$RMSE, CTrial12$RMSE,  
        BTrial7$RMSE, BTrial8$RMSE, BTrial9$RMSE, BTrial10$RMSE, BTrial11$RMSE, BTrial12$RMSE, 
        DTrial10$RMSE, DTrial11$RMSE, DTrial12$RMSE,
        ylim=c(0,2), 
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Decreasing Trend', line=-1)




vioplot(c(ATrial13$biasAll),c(ATrial14$biasAll),c( ATrial15$biasAll),c( ATrial16$biasAll),c( ATrial17$biasAll),c( ATrial18$biasAll), 
        c(CTrial16$biasAll),c( CTrial17$biasAll),c( CTrial18$biasAll), 
        c(BTrial13$biasAll),c( BTrial14$biasAll),c( BTrial15$biasAll),c( BTrial16$biasAll),c( BTrial17$biasAll),c( BTrial18$biasAll), 
        c(DTrial16$biasAll),c( DTrial17$biasAll),c( DTrial18$biasAll),
        ylim=c(-5,5), 
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Increasing Trend', line=-1)


vioplot(ATrial13$RMSE, ATrial14$RMSE, ATrial15$RMSE, ATrial16$RMSE, ATrial17$RMSE, ATrial18$RMSE, 
        CTrial16$RMSE, CTrial17$RMSE, CTrial18$RMSE, 
        BTrial13$RMSE, BTrial14$RMSE, BTrial15$RMSE, BTrial16$RMSE, BTrial17$RMSE, BTrial18$RMSE, 
        DTrial16$RMSE, DTrial17$RMSE, DTrial18$RMSE,
        ylim=c(0,2), 
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Increasing Trend', line=-1)

vioplot(c(ATrial19$biasAll),c( ATrial20$biasAll),c( ATrial21$biasAll),c( ATrial22$biasAll),c( ATrial23$biasAll),c( ATrial24$biasAll), 
        c(CTrial22$biasAll),c( CTrial23$biasAll),c( CTrial24$biasAll),  
        c(BTrial19$biasAll),c( BTrial20$biasAll),c( BTrial21$biasAll),c( BTrial22$biasAll),c( BTrial23$biasAll),c(BTrial24$biasAll), 
        c(DTrial22$biasAll),c( DTrial23$biasAll),c( DTrial24$biasAll),
        ylim=c(-5,5),
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Decrease - Increasing Trend', line=-1)


vioplot(ATrial19$RMSE, ATrial20$RMSE, ATrial21$RMSE, ATrial22$RMSE, ATrial23$RMSE, ATrial24$RMSE, 
        CTrial22$RMSE, CTrial23$RMSE, CTrial24$RMSE,  
        BTrial19$RMSE, BTrial20$RMSE, BTrial21$RMSE, BTrial22$RMSE, BTrial23$RMSE, BTrial24$RMSE, 
        DTrial22$RMSE, DTrial23$RMSE, DTrial24$RMSE,
        ylim=c(0,2), 
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Decrease - Increasing Trend', line=-1)
###################
dev.off()









##### CHANGE SCALE ############
# setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Plots")

png(filename="SN_BiasRMSE_Scale.png", 
    type="cairo",
    units="mm", 
    width=300, 
    height=375, 
    pointsize=16, 
    res=600)
###########
par(mfrow=c(4,2), mar=c(0.6,1.6,1.1,0), tcl=-0.25, mgp=c(2, 0.4, 0))
vioplot(c(ATrial1$biasAll), c(ATrial2$biasAll), c( ATrial3$biasAll), c( ATrial4$biasAll), c( ATrial5$biasAll), c( ATrial6$biasAll), 
        c(CTrial4$biasAll), c( CTrial5$biasAll), c( CTrial6$biasAll), 
        c(BTrial1$biasAll), c( BTrial2$biasAll), c( BTrial3$biasAll), c( BTrial4$biasAll), c( BTrial5$biasAll), c( BTrial6$biasAll),
        c(DTrial4$biasAll), c( DTrial5$biasAll), c( DTrial6$biasAll),
        ylim=c(-4.5,4.5),
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
abline(h=0)
title('No Trend', line=-1)
title('Relative Error', line=0.4, cex=1.5)


vioplot(ATrial1$RMSE, ATrial2$RMSE, ATrial3$RMSE, ATrial4$RMSE, ATrial5$RMSE, ATrial6$RMSE, 
        CTrial4$RMSE, CTrial5$RMSE, CTrial6$RMSE, 
        BTrial1$RMSE, BTrial2$RMSE, BTrial3$RMSE, BTrial4$RMSE, BTrial5$RMSE, BTrial6$RMSE,
        DTrial4$RMSE, DTrial5$RMSE, DTrial6$RMSE,
        ylim=c(0,2), 
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('No Trend', line=-1)
title('RMSE', line=0.4, cex=1.5)



par(mar=c(0.6,1.6,0,0), tcl=-0.25, mgp=c(2, 0.4, 0))
vioplot(c(ATrial7$biasAll), c( ATrial8$biasAll), c( ATrial9$biasAll), c( ATrial10$biasAll), c( ATrial11$biasAll), c( ATrial12$biasAll), 
        c(CTrial10$biasAll), c( CTrial11$biasAll), c( CTrial12$biasAll),  
        c(BTrial7$biasAll), c( BTrial8$biasAll), c( BTrial9$biasAll), c( BTrial10$biasAll), c( BTrial11$biasAll), c( BTrial12$biasAll), 
        c(DTrial10$biasAll), c( DTrial11$biasAll), c( DTrial12$biasAll),
        ylim=c(-1,1),
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
abline(h=0)
title('Decreasing Trend', line=-1)

vioplot(ATrial7$RMSE, ATrial8$RMSE, ATrial9$RMSE, ATrial10$RMSE, ATrial11$RMSE, ATrial12$RMSE, 
        CTrial10$RMSE, CTrial11$RMSE, CTrial12$RMSE,  
        BTrial7$RMSE, BTrial8$RMSE, BTrial9$RMSE, BTrial10$RMSE, BTrial11$RMSE, BTrial12$RMSE, 
        DTrial10$RMSE, DTrial11$RMSE, DTrial12$RMSE,
        ylim=c(0,0.5), 
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Decreasing Trend', line=-1)




vioplot(c(ATrial13$biasAll),c(ATrial14$biasAll),c( ATrial15$biasAll),c( ATrial16$biasAll),c( ATrial17$biasAll),c( ATrial18$biasAll), 
        c(CTrial16$biasAll),c( CTrial17$biasAll),c( CTrial18$biasAll), 
        c(BTrial13$biasAll),c( BTrial14$biasAll),c( BTrial15$biasAll),c( BTrial16$biasAll),c( BTrial17$biasAll),c( BTrial18$biasAll), 
        c(DTrial16$biasAll),c( DTrial17$biasAll),c( DTrial18$biasAll),
        ylim=c(-1,1), 
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
abline(h=0)
title('Increasing Trend', line=-1)


vioplot(ATrial13$RMSE, ATrial14$RMSE, ATrial15$RMSE, ATrial16$RMSE, ATrial17$RMSE, ATrial18$RMSE, 
        CTrial16$RMSE, CTrial17$RMSE, CTrial18$RMSE, 
        BTrial13$RMSE, BTrial14$RMSE, BTrial15$RMSE, BTrial16$RMSE, BTrial17$RMSE, BTrial18$RMSE, 
        DTrial16$RMSE, DTrial17$RMSE, DTrial18$RMSE,
        ylim=c(0,0.5), 
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Increasing Trend', line=-1)

vioplot(c(ATrial19$biasAll),c( ATrial20$biasAll),c( ATrial21$biasAll),c( ATrial22$biasAll),c( ATrial23$biasAll),c( ATrial24$biasAll), 
        c(CTrial22$biasAll),c( CTrial23$biasAll),c( CTrial24$biasAll),  
        c(BTrial19$biasAll),c( BTrial20$biasAll),c( BTrial21$biasAll),c( BTrial22$biasAll),c( BTrial23$biasAll),c(BTrial24$biasAll), 
        c(DTrial22$biasAll),c( DTrial23$biasAll),c( DTrial24$biasAll),
        ylim=c(-2,2),
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
abline(h=0)
title('Decrease - Increasing Trend', line=-1)


vioplot(ATrial19$RMSE, ATrial20$RMSE, ATrial21$RMSE, ATrial22$RMSE, ATrial23$RMSE, ATrial24$RMSE, 
        CTrial22$RMSE, CTrial23$RMSE, CTrial24$RMSE,  
        BTrial19$RMSE, BTrial20$RMSE, BTrial21$RMSE, BTrial22$RMSE, BTrial23$RMSE, BTrial24$RMSE, 
        DTrial22$RMSE, DTrial23$RMSE, DTrial24$RMSE,
        ylim=c(0,1), 
        col=c(rep('deepskyblue3',3),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Decrease - Increasing Trend', line=-1)
################
dev.off()

load("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial15\\ResultsList.RData")
load("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Base\\Trial15\\DFA_Results.RData")


##### Legend ############
setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\Plots")

png(filename="SN_BiasRMSE_Legend.png", 
    type="cairo",
    units="mm", 
    width=225, 
    height=50, 
    pointsize=16, 
    res=600)

par(mfrow=c(1,1), mar=c(0.6,1.6,1.1,0), tcl=-0.25, mgp=c(2, 0.4, 0))
plot(NULL)

legend('center', c('Constant q, 3 surveys', 'Knife-edge change in q, 3 surveys', 'Gradual change in q, 3 surveys',
                   'Constant q, 4 surveys', 'Knife-edge change in q, 4 surveys', 'Gradual change in q, 4 surveys'), 
       fill=c("deepskyblue3",'deepskyblue','lightskyblue1','darkolivegreen4','darkolivegreen3','darkolivegreen1'),ncol=2)

dev.off()







########################## EXTRA PLOTTING FOR POSTER #################

setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\FourIndicesGradual\\Trial22")

DFATrends = read.csv('DFATrends.csv')
N_All = read.csv('N77.csv')
FitRatios = read.csv('FitRatios.csv')
Iy1 = read.csv('Iy1.csv')
Iy2 = read.csv('Iy2.csv')
Iy3 = read.csv('Iy3.csv')
Iy4 = read.csv('Iy4.csv')
upCI_DFA = read.csv('upCI_DFATrends.csv')
lowCI_DFA = read.csv('lowCI_DFATrends.csv')

FitRatios[77,2:5]

DFAT = DFATrends[77,2:27]
Nall = N_All[2:27, 2]
I1 = Iy1[77,41:66]
I2 = Iy2[77,41:66]
I3 = Iy3[77,41:66]
I4 = c(rep(NA, 15), as.matrix(Iy4[77,56:66]))
upCI = upCI_DFA[77, 2:27]
lowCI = lowCI_DFA[77, 2:27]


par(mfrow=c(1,1), mar=c(3,3,1,1))
plot(x=1:26, y=I1, type='l', ylim=c(0,10), lwd=2, ylab='Index of abundance', xlab='Year')
lines(x=1:26, y=I2, type='l', lwd=2)
lines(x=1:26, y=I3, type='l', lwd=2)
lines(x=1:26, y=I4, type='l', lwd=2)




setwd("D:\\vspace1\\DFA_Simulation\\Atl_SN\\NEW_ANALYSES_2019\\\\Plots")

png(filename="SN_TrialD22_77_Indices.png", 
    type="cairo",
    units="mm", 
    width=225, 
    height=150, 
    pointsize=16, 
    res=600)

par(mfrow=c(1,1), mar=c(3,3,1,1), tcl=-0.2, mgp=c(1.5, 0.3, 0))
plot(x=1:26, y=I3-mean(as.matrix(I3)), type='l', lwd=2, col='steelblue4', ylab='Mean-scaled index of abundance', xlab='Year', ylim=c(-4, 4))
lines(x=1:26, y=I1-mean(as.matrix(I1)), type='l', lwd=2, col='darkolivegreen3')
lines(x=1:26, y=I2-mean(as.matrix(I2)), type='l', lwd=2, col='darkorchid')
lines(x=1:26, y=I4-mean(as.matrix(I4),na.rm=T), type='l', lwd=2, col='deepskyblue3')

dev.off()


png(filename="SN_TrialD22_77_DFAResults.png", 
    type="cairo",
    units="mm", 
    width=225, 
    height=150, 
    pointsize=16, 
    res=600)

par(mfrow=c(1,1), mar=c(3,3,1,1), tcl=-0.2, mgp=c(1.5, 0.3, 0))
years=1:26
df <- data.frame(years, index = c(as.matrix(DFAT)), lowerCI = c(as.matrix(lowCI)), upperCI = c(as.matrix(upCI)))
with(df, plot(index ~ years, type='n', bty='L',xlab='Year', ylab='Relative Abundance', ylim=c(min(lowerCI), max(upperCI)), axes=FALSE))
axis(1,labels=TRUE)
axis(2,labels=FALSE)
x <- as.numeric(df$years)
polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
with(df, lines(x, index, type='l', lwd=2, pch=16))
# title(paste("Iteration",l,sep=" "))
# abline(h=0)
lines(years, rescale(Nall, to=c(min(df$index), max(df$index))), type='l',col="darkorchid", lwd=2)
box()
legend('topright', c("DFA predicted trend +/- 95% CIs",'Rescaled simulated abundance'), lwd=2, col=c('black','darkorchid'))

dev.off()