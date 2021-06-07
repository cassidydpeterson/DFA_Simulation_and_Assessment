#####################################
# Sandbar shark simulation
# DFA Simulation/Management Test
# WORKING VERSION: Updated Feb 2019
#####################################

save.image("N:/Documents/DFA_Simulation/SB/SB_Sim.RData")
load("N:/Documents/DFA_Simulation/SB/SB_Sim.RData")

########################################### READ IN ########################################################################
library(scales)
library(MARSS)
library(vioplot)
# from: http://blog.revolutionanalytics.com/2016/08/simulating-form-the-bivariate-normal-distribution-in-r-1.html
rbvn<-function (n, m1, s1, m2, s2, rho) {
  X1 <- rnorm(n, m1, s1)
  X2 <- rnorm(n, m2 + (s2/s1) * rho * 
                (X1 - m1), sqrt((1 - rho^2)*s2^2))
  data.frame(Linf = X1, K = X2)
}
Linf_M=175.5
K_M=0.143
t0_M=-2.388
Linf_F=183.3
K_F=0.124
t0_F=-3.098
# vonBert = rbvn(nrow(df.expanded),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
ages = 1:31
# vonBert$L0 = 

# VON BERT:  Lt = Linf *(1-exp(-K*(t-t0)))
# SCHNUTE: L1 + (L2 - L1) * ((1-exp(-K*(t-t1)))/(1-exp(-K*(t2-t1))))
(L0_F = Linf_F * (1-exp(-K_F*(0-t0_F))) )# L0 = 58.46759
(L1_F = Linf_F * (1-exp(-K_F*(1-t0_F))) )# L1 = 73.02556
(L0_M = Linf_M * (1-exp(-K_M*(0-t0_M))) )# L0 = 50.76955
(L1_M = Linf_M * (1-exp(-K_M*(1-t0_M))) )# L1 = 67.38937

# Sel = read.csv("D:\\vspace1\\DFA_Simulation\\SB\\Selectivity_new.csv")
Sel = read.csv("D:\\vspace1\\DFA_Simulation\\SB\\Selectivity_new.csv")
FullSel = Sel
Sel = Sel[-1,]

Equilibrium_Nay_F = c(500.00065, 425.90205, 362.78465, 309.02105, 263.22506, 224.21589, 191.48498, 170.37628, 151.59453, 
                      134.88322, 120.01412, 106.78413,  95.01259,  84.53869,  75.21941,  66.92746, 59.54959,  52.98503,  
                      47.14412,  41.94710,  37.32298,  33.20861,  29.54780,  26.29054,  23.39236,  20.81366,  18.51922,  
                      16.47772,  14.66127,  13.04506,  11.60701)

Equilibrium_Nay_M = c(500.00065, 425.90205, 362.78465, 309.02105, 263.22506, 224.21589, 191.48498, 170.37628, 151.59453, 
                      134.88322, 120.01412, 106.78413,  95.01259,  84.53869,  75.21941,  66.92746, 59.54959,  52.98503,  
                      47.14412,  41.94710,  37.32298,  33.20861,  29.54780,  26.29054,  23.39236,  20.81366,  18.51922,  
                      16.47772,  14.66127,  13.04506,  11.60701)

ages = 1:31
A = max(ages)
years=100

# Define maturity and fecundity  
mat_F = as.vector(c(0, 0, 0, 0, 0, 0.01, 0.02, 0.03, 0.06, 0.12, 0.21, 0.33, 0.49, 0.65,  0.78, 0.88,  0.93, 0.96, 0.98, 0.99, 0.99, 1, 1, 1, 1,  1, 1, 1, 1, 1, 1 ))
fec = c(6.61072819,	7.027399402,	7.39547835,	7.720631873,	8.00786594,	8.261602725,	8.485748685,	8.683754708,	8.858669237,	
        9.013185206,	9.149681498,	9.27025957,	9.376775809,	9.470870107,	9.553991113,	9.627418534,	9.692282837,	9.749582655,	9.80020016,	
        9.844914642,	9.884414515,	9.919307906,	9.950132025,	9.97736143,	10.00141534,	10.02266407,	10.04143478,	10.05801644,	10.07266435,	
        10.08560401,	10.09703465)

M_const = c(0.1604,	0.1604,	0.1604,	0.1604,	0.1604,	0.1578,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	
            0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	
            0.1168,	0.1168,	0.1168,	0.1168)

FReal = as.vector(c(rep(0,44),0.001,	0.001,	0.001,	0.001,	0.001,	0.001,	0.001,	0.002,	
                    0.004,	0.014,	0.055,	0.248,	0.209,	0.386,	0.271,	0.226,	0.318,	0.234,	0.369,	0.347,	
                    0.418,	0.38,	0.395,	0.353,	0.451,	0.441,	0.498,	0.387,	0.374,	0.374,	0.195, 0.235,	
                    0.215,	0.182,	0.157,	0.144,	0.17,	0.108,	0.06,	0.068,	0.073,	0.052,	0.059,	0.056,	
                    0.05,	0.054, rep(0.054,10)))  # NOTE THIS IS REAL!!! ACCORDING TO SEDAR 54
F0 = as.vector(c(rep(0.0, 100)))
F1 = as.vector(c(rep(0,45), rep(0.1, 10), rep(0.3, 10), rep(0.2, 10), rep(0.05, 25))) # F1
F2 = as.vector(c(rep(0,49), rep(0.4, 7), rep(0.025, 25), rep(0.2, 10), rep(0.05, 9))) # F2

load(file = "D:\\vspace1\\DFA_Simulation\\SB\\Iay_eq_F0.RData")
load(file = "D:\\vspace1\\DFA_Simulation\\SB\\Iay_eq_F1.RData")
load(file = "D:\\vspace1\\DFA_Simulation\\SB\\Iay_eq_F2.RData")
load(file = "D:\\vspace1\\DFA_Simulation\\SB\\Iay_eq_Real.RData")

q0 = rep(0.01, 100)                                                                     # q0 = constant
q1 = rep(c(0.01, 0.02), 50)                                                             # q1 = period 2
q2 = rep(c(0.01, 0.02, 0.015), length=100)                                              # q2 = period 3
q3 = c(rep(0.01, 50), rep(seq(0.01, 0.045, by=0.005), each = 5), rep(0.045, 10))        # q3 = step increase (step=5)
q4 = c(rep(0.05, 50), rep(seq(0.05, 0.0325, by=-0.0025), each = 5), rep(0.0325, 10))    # q4 = step decrease (step=5)
q5 = c(rep(seq(0.03, 0.05, by=0.005), seq(0.05, 0.03, by=-0.005), length=100))          #q5 = increase then decrease
q6 = c(rep(0.01, 50), seq(0.01, 0.049, by=0.001), rep(0.049,10))                        # q6 = increase
q7 = c(rep(0.049, 50), seq(0.049, 0.01, by=-0.001), rep(0.01,10))                       # q7 = decrease
set.seed(430)
q8 = runif(100, 0.03, 0.05)                                                             # q8 = random between 0.01 and 0.05
set.seed(430)
q9 = jitter(q6, amount=0.005)                                                           # q9 = jitter increase
set.seed(431)
q10 = jitter(q7, amount=0.005)                                                          # q10 = jitter decrease
set.seed(432)
q11 = jitter(q1, amount=0.005)                                                         # q11 = jitter period 2
set.seed(433)
q12 = jitter(q2, amount=0.005)                                                         # q12 = jitter period 3                                                         # jitter decrease
set.seed(434)
q13 = jitter(q0, amount=0.005)                                                         # q13 = jitter constant
set.seed(435)
q14 = jitter(q0, amount=0.003)                                                         # q14 = jitter constant 2
set.seed(436) 
q15 = jitter(q0, amount=0.002)                                                         # q15 = jitter constant 3
set.seed(437) 
q16 = jitter(q0, amount=0.004)                                                         # q16 = jitter constant 3
set.seed(438) 
q17 = jitter(q0, amount=0.003)                                                         # q17 = jitter constant 4
set.seed(439) 
q18 = jitter(q0, amount=0.002)                                                         # q18 = jitter constant 5
set.seed(440) 
q19 = jitter(q0, amount=0.004)                                                         # q19 = jitter constant 6
set.seed(441) 
q20 = jitter(q0, amount=0.003)                                                         # q20 = jitter constant 7
set.seed(441) 
q21 = jitter(q0, amount=0.002)                                                         # q21 = jitter constant 8



# F1 = COM GOM
# F2 = COM SA
# F3 = REC MEX
# F4 = MEN DISC

# Survey 1 = VIMS
# Survey 2 = LPS
# Survey 3 = PLLOP
# Survey 4 = NMFS LL SE
# Survey 5 = CST NE LL
# Survey 6 = BLLOP 1
# Survey 7 = NMFS NE

CV1=0.38
CV2=0.48
CV3=0.65
CV4=0.24
CV5=0.30
CV6=0.36
CV7=0.40

############ LOAD SAVED TRIALS #################
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial1\\SB1.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial2\\SB2.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial3\\SB3.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial4\\SB4.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial5\\SB5.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial6\\SB6.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial7\\SB7.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial8\\SB8.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial9\\SB9.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial10\\SB10.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial11\\SB11.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial12\\SB12.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial13\\SB13.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial14\\SB14.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial15\\SB15.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial16\\SB16.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial17\\SB17.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial18\\SB18.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial19\\SB19.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial20\\SB20.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial21\\SB21.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial22\\SB22.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial23\\SB23.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial24\\SB24.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial25\\SB25.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial26\\SB26.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial27\\SB27.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial28\\SB28.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial29\\SB29.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial30\\SB30.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial31\\SB31.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial32\\SB32.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial33\\SB33.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial34\\SB34.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial35\\SB35.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial36\\SB36.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial37\\SB37.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial38\\SB38.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial39\\SB39.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial101\\SB101.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial102\\SB102.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial103\\SB103.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial104\\SB104.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial105\\SB105.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial106\\SB106.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial107\\SB107.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial108\\SB108.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial109\\SB109.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial110\\SB110.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial111\\SB111.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial112\\SB112.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial113\\SB113.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial114\\SB114.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial115\\SB115.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial116\\SB116.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial117\\SB117.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial118\\SB118.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial119\\SB119.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial120\\SB120.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial121\\SB121.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial122\\SB122.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial123\\SB123.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial124\\SB124.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial125\\SB125.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial126\\SB126.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial127\\SB127.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial128\\SB128.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial129\\SB129.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial130\\SB130.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial131\\SB131.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial132\\SB132.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial133\\SB133.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial134\\SB134.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial135\\SB135.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial136\\SB136.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial137\\SB137.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial138\\SB138.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial139\\SB139.RData")


load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial40\\SB40.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial41\\SB41.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial140\\SB140.RData")
load( file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial141\\SB141.RData")

############# AGE-BASED SELECTIVITY: #################
Sel = read.csv("D:\\vspace1\\DFA_Simulation\\SB\\Selectivity_new.csv")
FullSel = Sel
Sel = Sel[-1,]
# Sel[1,]
# Age classes: 31
# Longevity: 31
# Age at 50% maturity; A50 = 13
# GOM: A50 = 1.6 -- used in 2013 assessment

# Virgin conditions: 1960
# Consider stock-recruitment relationship following a functional form...

# Reference Porch 2003: Preliminary assessment of atlantic white marlin (Tetrapturus albidus) using a state-space 
#   implementation of an age-structured production model

# SEDAR 34 estimates steepness of 0.3. 



# catch streams included: 
#   4 commercial: 
#   10 surveys


#age structure
ages=1:31
A=max(ages)

#initial population
# No=1000

#length of simulation
# yrs=1000


#F= fishing mortality (estimated Fmsy=0.06)
# FM=0.00

# iterations=1000




#ages
ages = 1:31
A = max(ages)

# Define maturity function 
mat_F = as.vector(c(0, 0, 0, 0, 0, 0.01, 0.02, 0.03, 0.06, 0.12, 0.21, 0.33, 0.49, 0.65,  0.78, 0.88,  0.93, 0.96, 0.98, 0.99, 0.99, 1, 1, 1, 1,  1, 1, 1, 1, 1, 1 ))
# mat_F = as.vector(c(0, 0, 0, 0, 0, 0, 0.01, 0.02, 0.03, 0.06, 0.12, 0.21, 0.33, 0.49, 0.65,  0.78, 0.88,  0.93, 0.96, 0.98, 0.99, 0.99, 1, 1, 1, 1,  1, 1, 1, 1, 1, 1 ))

# Define fecundity (ALL pups_ sex ratio = 1:1)
# fec = c(6.139049766, 6.61072819,	7.027399402,	7.39547835,	7.720631873,	8.00786594,	8.261602725,	8.485748685,	8.683754708,	8.858669237,	
#         9.013185206,	9.149681498,	9.27025957,	9.376775809,	9.470870107,	9.553991113,	9.627418534,	9.692282837,	9.749582655,	9.80020016,	
#         9.844914642,	9.884414515,	9.919307906,	9.950132025,	9.97736143,	10.00141534,	10.02266407,	10.04143478,	10.05801644,	10.07266435,	
#         10.08560401,	10.09703465)
fec = c(6.61072819,	7.027399402,	7.39547835,	7.720631873,	8.00786594,	8.261602725,	8.485748685,	8.683754708,	8.858669237,	
        9.013185206,	9.149681498,	9.27025957,	9.376775809,	9.470870107,	9.553991113,	9.627418534,	9.692282837,	9.749582655,	9.80020016,	
        9.844914642,	9.884414515,	9.919307906,	9.950132025,	9.97736143,	10.00141534,	10.02266407,	10.04143478,	10.05801644,	10.07266435,	
        10.08560401,	10.09703465)
#

# length of simulation
years = 1000


N0 = 1000
# N0 = 1220.88025

################## CHANGE INPUTS 
#__________________________________________________


################# define natural mortality
# WITHOUT M0
M_const = c(0.1604,	0.1604,	0.1604,	0.1604,	0.1604,	0.1578,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	
            0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	
            0.1168,	0.1168,	0.1168,	0.1168)
# WITH M0
# M_const = c(0.1604,	0.1604,	0.1604,	0.1604,	0.1604,	0.1604,	0.1578,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	
#             0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	
#             0.1168,	0.1168,	0.1168,	0.1168)


################# Assume constant fishing mortality
F_const1 = as.vector(rep(0.0,years))
F_const2 = as.vector(rep(0.0,years))
F_const3 = as.vector(rep(0.0,years))
F_const4 = as.vector(rep(0.0,years))
F_const = rbind(F_const1, F_const2, F_const3, F_const4)
F_Sum = apply(F_const, 2, FUN=sum)
# Store results



#################### Run SImulation 
# Create matrices to store results
Nay_M = matrix(nrow=A, ncol=years)
Nay_F = matrix(nrow=A, ncol=years)
Cay1_M = matrix(nrow=A, ncol=years)
Cay2_M = matrix(nrow=A, ncol=years)
Cay3_M = matrix(nrow=A, ncol=years)
Cay4_M = matrix(nrow=A, ncol=years)
Cay1_F = matrix(nrow=A, ncol=years)
Cay2_F = matrix(nrow=A, ncol=years)
Cay3_F = matrix(nrow=A, ncol=years)
Cay4_F = matrix(nrow=A, ncol=years)
M0 = vector(length=years)
Npups = vector(length=years)




####### LFSR - FOR STOCK SYNTHESIS ####################


# Inputs: 
M = M_const
# M = M_const/2
FM = F_const


########## without plus group 
# 
# Nay_M[1,1] = N0/2
# Nay_F[1,1] = N0/2
# 
# # populate first row
# for (i in 1:(A-1)){
#   Nay_M[i+1,1] = Nay_M[i,1]*exp(-M[i])
#   Nay_F[i+1,1] = Nay_F[i,1]*exp(-M[i])
# }
# #                                           WITH OR WITHOUT PLUS GROUP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Nay_M[A,1] = ( Nay_M[A-1,1]*exp(-M[A-1]) ) #+ ( Nay_M[A,1]*exp(-M[A]) )#####################################################################################
# Nay_F[A,1] = ( Nay_F[A-1,1]*exp(-M[A-1]) ) #+ ( Nay_F[A,1]*exp(-M[A]) )#####################################################################################
# 
# 
# #
# # TO ESTIMATE Survival-based Recruitmetn parameters from (Taylor et al. 2012)
# # 1. generate initial N at age vector (Na0), given M at age vector and arbitrary N0
# # 2. calculate Npups0, given maturity & fecundity at age vectors and Na0 as estimated above
# # 3. Z0 = -ln(N0/Npups0)
# # 4. set Zmin as fixed M for age 0
# # 5. zfrac = 1- (zmin/z0)
# # 6. beta =  ( ln( 1 - (ln(h/0.2) / (z0 * zfrac) ) ) ) / (ln(0.2))
# 
# 
# Npups0 = sum((Nay_F[,1]/2.5)*mat_F*fec); Npups0 #2919.848
# # Npups0 = sum(Nay_F[,1]*mat_F*fec)
# z0 = -log(N0/Npups0 ); z0
# #z0 = 1.071531
# (s0 = exp(-z0)) #0.3424836
# zmin = 0.1604
# (smax = exp(-zmin)) #0.851803
# ( zfrac = (log(smax) - log(s0))/(-log(s0)) ) # check. zfrac = 0.8503077
# ( zfrac = 1 - (zmin/z0) )
# h=0.3
# # h=0.205
# # h=0.25# USE THIS ONE FOR LOH
# beta = (log(1-(log(h/0.2)/(z0*zfrac))))/log(0.2)
# beta #0.3658483
# 
# ( S_frac = (zmin - z0)/(0-z0) ) #0.8503077
# 
# Z_max = z0 + S_frac*(0.0-z0) #0.1604
# 
# 
# 
# 
# Z0 = 1.071531
# Zmin = 0.1604
# zfrac = 0.8503077
# beta = 0.3658483 
# Npups_eq = 2919.848 #Npups0





############### WITH PLUS GROUP
Nay_M[1,1] = N0/2
Nay_F[1,1] = N0/2

# populate first row
for (i in 1:(A-1)){
  Nay_M[i+1,1] = Nay_M[i,1]*exp(-M[i])
  Nay_F[i+1,1] = Nay_F[i,1]*exp(-M[i])
}
#                                           WITH OR WITHOUT PLUS GROUP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Nay_M[A,1] = ( Nay_M[A-1,1]*exp(-M[A-1]) ) + ( Nay_M[A,1]*exp(-M[A]) )#####################################################################################
Nay_F[A,1] = ( Nay_F[A-1,1]*exp(-M[A-1]) ) + ( Nay_F[A,1]*exp(-M[A]) )#####################################################################################

####### FOR STOCK SYNTHESIS 
#
# TO ESTIMATE Survival-based Recruitmetn parameters from (Taylor et al. 2012)
# 1. generate initial N at age vector (Na0), given M at age vector and arbitrary N0
# 2. calculate Npups0, given maturity & fecundity at age vectors and Na0 as estimated above
# 3. Z0 = -ln(N0/Npups0)
# 4. set Zmin as fixed M for age 0
# 5. zfrac = 1- (zmin/z0)
# 6. beta =  ( ln( 1 - (ln(h/0.2) / (z0 * zfrac) ) ) ) / (ln(0.2))


Npups0 = sum((Nay_F[,1]/2.5)*mat_F*fec); Npups0 #2961.558
# Npups0 = sum(Nay_F[,1]*mat_F*fec)
z0 = -log(N0/Npups0 ); z0
#z0 = 1.085716
(s0 = exp(-z0)) #0.3424836
zmin = 0.1604
(smax = exp(-zmin)) #0.851803
( zfrac = (log(smax) - log(s0))/(-log(s0)) ) # check. zfrac = 0.8522633
( zfrac = 1 - (zmin/z0) )
h=0.3
# h=0.205
# h=0.25
beta = (log(1-(log(h/0.2)/(z0*zfrac))))/log(0.2)
beta #0.3582577

( S_frac = (zmin - z0)/(0-z0) ) #0.8522633

Z_max = z0 + S_frac*(0.0-z0) #0.1604




Z0 = 1.085716
Zmin = 0.1604
zfrac = 0.8522633
beta = 0.3582577
Npups_eq = 2961.558 #Npups0
# REAL Npups_eq = 3617.708
# REAL R0 = 1220.88025
# REAL lnR0 = 7.1073

#Define F vector: F FOR EACH FISHERY
#
# years = 10000

for (i in 1:(years-1)){
  
  Npups[i] = sum((Nay_F[,i]/2.5)*mat_F*fec) # equilibrium Npups = 2919.848 
  
  M0[i] = ( (1- (Npups[i]/Npups_eq)^beta ) * (Zmin - Z0) ) + Z0
  
  
  Nay_F[1,i+1] = 0.5 * Npups[i]*exp(-M0[i])  # pups produced that survive to age 1; Nay[1,] equilibrium = xxxx
  Nay_M[1,i+1] = 0.5 * Npups[i]*exp(-M0[i])  # pups produced that survive to age 1; Nay[1,] equilibrium = xxxx
  
 
  Zay_F = vector(length=length(M))
  Zay_M = vector(length=length(M))
  for (j in 1:A){
    Zay_F[j] = M[j] + (Sel[j,2]*FM[1,i]) + (Sel[j,3]*FM[2,i]) + (Sel[j,4]*FM[3,i]) + (Sel[j,5]*FM[4,i])
    Zay_M[j] = M[j] + (Sel[j,2]*FM[1,i]) + (Sel[j,3]*FM[2,i]) + (Sel[j,4]*FM[3,i]) + (Sel[j,5]*FM[4,i])
  } # END J LOOP
  
  for (j in 1:(A-2)){
    Nay_F[j+1,i+1] = Nay_F[j,i]*exp(-Zay_F[j]) 
    Nay_M[j+1,i+1] = Nay_M[j,i]*exp(-Zay_M[j]) 
    
    
    Cay1_F[j,i+1] = Nay_F[j,i]*((Sel[j,2]*FM[1,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
    Cay2_F[j,i+1] = Nay_F[j,i]*((Sel[j,4]*FM[2,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
    Cay3_F[j,i+1] = Nay_F[j,i]*((Sel[j,6]*FM[3,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
    Cay4_F[j,i+1] = Nay_F[j,i]*((Sel[j,8]*FM[4,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
    
    Cay1_M[j,i+1] = Nay_M[j,i]*((Sel[j,3]*FM[1,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
    Cay2_M[j,i+1] = Nay_M[j,i]*((Sel[j,5]*FM[2,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
    Cay3_M[j,i+1] = Nay_M[j,i]*((Sel[j,7]*FM[3,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
    Cay4_M[j,i+1] = Nay_M[j,i]*((Sel[j,8]*FM[4,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
    
    
  } # END J LOOP
  
  Nay_F[A,i+1] = ( Nay_F[A-1,i]*exp(-Zay_F[A-1]) ) #+ ( Nay_F[A,i]*exp(-Zay_F[A]) )
  Nay_M[A,i+1] = ( Nay_M[A-1,i]*exp(-Zay_M[A-1]) ) #+ ( Nay_M[A,i]*exp(-Zay_M[A]) )
  
  Cay1_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,2]*FM[1,i])/Nay_F[j])*(1-exp(-Nay_F[j])) 
  Cay2_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,4]*FM[2,i])/Nay_F[j])*(1-exp(-Nay_F[j])) 
  Cay3_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,6]*FM[3,i])/Nay_F[j])*(1-exp(-Nay_F[j])) 
  Cay4_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,8]*FM[4,i])/Nay_F[j])*(1-exp(-Nay_F[j])) 
  
  Cay1_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,3]*FM[1,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
  Cay2_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,5]*FM[2,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
  Cay3_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,7]*FM[3,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
  Cay4_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,8]*FM[4,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
  
 
   
  Cay1_F[A,i]=Nay_F[A,i]*((Sel[A,2]*FM[1,i])/Nay_F[j])*(1-exp(-Nay_F[j])) 
  Cay2_F[A,i]=Nay_F[A,i]*((Sel[A,4]*FM[2,i])/Nay_F[j])*(1-exp(-Nay_F[j])) 
  Cay3_F[A,i]=Nay_F[A,i]*((Sel[A,6]*FM[3,i])/Nay_F[j])*(1-exp(-Nay_F[j])) 
  Cay4_F[A,i]=Nay_F[A,i]*((Sel[A,8]*FM[4,i])/Nay_F[j])*(1-exp(-Nay_F[j])) 
  
  Cay1_M[A,i]=Nay_M[A,i]*((Sel[A,3]*FM[1,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
  Cay2_M[A,i]=Nay_M[A,i]*((Sel[A,5]*FM[2,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
  Cay3_M[A,i]=Nay_M[A,i]*((Sel[A,7]*FM[3,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
  Cay4_M[A,i]=Nay_M[A,i]*((Sel[A,8]*FM[4,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
  
} # END I LOOP



plot(1:years, colSums(Nay_M) + colSums(Nay_F),type='l', ylim=c(0, 15000))#, ylim=c(0,max(colSums(Nay))+100))
lines(1:years, Npups,type='l')#, ylim=c(0,max(colSums(Nay))+100))
lines(1:years, colSums(Nay_M) + colSums(Nay_F),type='l',col='red')
tail(colSums(Nay_M))
tail(colSums(Nay_F))

# View(tail(Nay_M))

# Equilibrium conditions:
Nay_F[,990:1000]
Nay_M[,990:1000]
Equilibrium_Nay_F = Nay_F[,1000] 
Equilibrium_Nay_M = Nay_M[,1000] 


### ALTER INPUTS 
# length of simulation
years = 1000


# Fishing mortality will vary over time. For now, set it to a constant 0.1

# define natural mortality
M_const = c(0.1604,	0.1604,	0.1604,	0.1604,	0.1604,	0.1578,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	
            0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	0.1168,	
            0.1168,	0.1168,	0.1168,	0.1168)


# FOR SN !!!
# Assume constant fishing mortality
# F_var = as.vector(rep(0.0,years))
# F_var = as.vector(rep(0.2,years)) # Trials 1-6
# F_var = as.vector(rep(0.1,years)) # 
# # F_var = as.vector(c(rep(0.1, 20), rep(0.2, 20), rep(0.3, 20), rep(0.4, 20), rep(0.2, 20)))
# # F_var = as.vector(c(rep(0.1, 50), rep(0.2, 50)))
# # F_var = as.vector(c(rep(0.4, 50), rep(0.2, 50)))
# F_var = as.vector(c(rep(0.4, 50), rep(0, 50))) # I like this one. Trials 13-18
# F_var = as.vector(c(rep(0.0, 50), rep(0.4, 50))) # I like this one. Trials 7-12
# F_var = as.vector(c(rep(0,40), rep(0.4, 10), rep(0.2, 5), rep(0.05, 45))) # Trials 19-24
# # F_var = as.vector(c(rep(0,30), rep(0.4, 20), rep(0.2, 5), rep(0.05, 45))) # 
# # F_var = as.vector(c(rep(0.0, 50), rep(0.3, 50)))
# # F_var = as.vector(c(rep(0.0, 50), rep(0.2, 50)))
# # F_var = as.vector(c(rep(0.6, 25), rep(0.2, 50), rep(0.6, 25)))
# # F_var = as.vector(c(rep(0.1, 50), rep(0.2, 50)))
# # F_var = as.vector(c(rep(0.05, 50), rep(0.1, 50))) # BASE CASE
# # F_var = as.vector(c(rep(0.05, 50), rep(0.0, 50)))

# F_var = as.vector(rep(1,years)) # Test

F_var = as.vector(c(rep(0.2, 100)))
F_var = as.vector(c(rep(0.1, years)))
F_var = as.vector(c(rep(0,45), rep(0.05, 10), seq(0.05,0.45, by=0.05), seq(0.4, 0.1, by=-0.1), rep(0.05, 32))) # F?
F_var = as.vector(c(rep(0.3,55), 0.35, 0.35, seq(0.4, 0.1, by=-0.1), seq(0.05, 0.01, by=-0.005), rep(0.01, 30) ) )  #, rep(0.05, 32)))
F_var = as.vector(c(rep(0,45), rep(0.4, 10), rep(0.025, 20), rep(0.2, 10), rep(0.05, 15))) #
F_var = as.vector(c(rep(4,5), rep(0.025, 95))) # F?

FReal = as.vector(c(rep(0,44),0.001,	0.001,	0.001,	0.001,	0.001,	0.001,	0.001,	0.002,	
                    0.004,	0.014,	0.055,	0.248,	0.209,	0.386,	0.271,	0.226,	0.318,	0.234,	0.369,	0.347,	
                    0.418,	0.38,	0.395,	0.353,	0.451,	0.441,	0.498,	0.387,	0.374,	0.374,	0.195, 0.235,	
                    0.215,	0.182,	0.157,	0.144,	0.17,	0.108,	0.06,	0.068,	0.073,	0.052,	0.059,	0.056,	
                    0.05,	0.054, rep(0.054,10)))  # NOTE THIS IS REAL!!! ACCORDING TO SEDAR 54
F0 = as.vector(c(rep(0.0, 100)))
F1 = as.vector(c(rep(0,45), rep(0.1, 10), rep(0.3, 10), rep(0.2, 10), rep(0.05, 25))) # F1
F2 = as.vector(c(rep(0,49), rep(0.4, 7), rep(0.025, 25), rep(0.2, 10), rep(0.05, 9))) # F2

F_var = FReal
F_var = F1

# length(F_var)
# plot(F_var)
# abline(v=51, lwd=2, col="red")
# abline(v=90, lwd=2, col="red")

Prop_F = c(0.2, 0.044, 0.754, 0.002)
sum(Prop_F)
Tot_F = F_var
F_const1 = Tot_F * Prop_F[1]
F_const2 = Tot_F * Prop_F[2]
F_const3 = Tot_F * Prop_F[3]
F_const4 = Tot_F * Prop_F[4]
F_const = rbind(F_const1, F_const2, F_const3, F_const4)
apply(F_const, 2, FUN=sum)
FM = F_const

#### INITIALIZE EQUILIBRIUM CONDITIONS 
years=100
# Create matrices to store results
Nay_M = matrix(nrow=A, ncol=years)
Nay_F = matrix(nrow=A, ncol=years)
Cay1_M = matrix(nrow=A, ncol=years)
Cay2_M = matrix(nrow=A, ncol=years)
Cay3_M = matrix(nrow=A, ncol=years)
Cay4_M = matrix(nrow=A, ncol=years)
Cay1_F = matrix(nrow=A, ncol=years)
Cay2_F = matrix(nrow=A, ncol=years)
Cay3_F = matrix(nrow=A, ncol=years)
Cay4_F = matrix(nrow=A, ncol=years)
Iay_M = matrix(nrow=A, ncol=years)
Iay_F = matrix(nrow=A, ncol=years)


q_survey=0.01
CV=0.5


M0 = vector(length=years)
Npups = vector(length=years)
# REMEMBER: S (annual survival) is related to exp(-Z), where Z is total mortality. So, M0 ~= Z0 : S0 ~= exp(-M0). 

# Inputs: 
M = M_const
# FM = F_const
# FM = F_var

# Set up equilibrium conditions
Nay_F[,1] = Equilibrium_Nay_F
Nay_M[,1] = Equilibrium_Nay_M


Z0 = 1.085716
Zmin = 0.1604
zfrac = 0.8522633
beta = 0.3582577
Npups_eq = 2961.558 #Npups0


for (i in 1:(years-1)){
  
  Npups[i] = sum((Nay_F[,i]/2.5)*mat_F*fec) # equilibrium Npups = 6990.286
  
  M0[i] = ( (1- (Npups[i]/Npups_eq)^beta ) * (Zmin - Z0) ) + Z0
  
  
  Nay_F[1,i+1] = 0.5 * Npups[i]*exp(-M0[i])  # pups produced that survive to age 1; Nay[1,] equilibrium = xxxx
  Nay_M[1,i+1] = 0.5 * Npups[i]*exp(-M0[i])  # pups produced that survive to age 1; Nay[1,] equilibrium = xxxx
  
  
  Zay_F = vector(length=length(M))
  Zay_M = vector(length=length(M))
  for (j in 1:A){
    Zay_F[j] = M[j] + (Sel[j,2]*FM[1,i]) + (Sel[j,4]*FM[2,i]) + (Sel[j,6]*FM[3,i]) + (Sel[j,8]*FM[4,i])
    Zay_M[j] = M[j] + (Sel[j,3]*FM[1,i]) + (Sel[j,5]*FM[2,i]) + (Sel[j,7]*FM[3,i]) + (Sel[j,8]*FM[4,i])
  } # END J LOOP
  
  for (j in 1:(A-2)){
    Nay_F[j+1,i+1] = Nay_F[j,i]*exp(-Zay_F[j]) 
    Nay_M[j+1,i+1] = Nay_M[j,i]*exp(-Zay_M[j]) 
    
    
    Cay1_F[j,i] = Nay_F[j,i]*((Sel[j,2]*FM[1,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
    Cay2_F[j,i] = Nay_F[j,i]*((Sel[j,4]*FM[2,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
    Cay3_F[j,i] = Nay_F[j,i]*((Sel[j,6]*FM[3,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
    Cay4_F[j,i] = Nay_F[j,i]*((Sel[j,8]*FM[4,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
    
    Cay1_M[j,i] = Nay_M[j,i]*((Sel[j,3]*FM[1,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
    Cay2_M[j,i] = Nay_M[j,i]*((Sel[j,5]*FM[2,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
    Cay3_M[j,i] = Nay_M[j,i]*((Sel[j,7]*FM[3,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
    Cay4_M[j,i] = Nay_M[j,i]*((Sel[j,8]*FM[4,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
    
    
  } # END J LOOP
  
  Nay_F[A,i+1] = Nay_F[A-1,i]*exp(-Zay_F[A-1]) #+ Nay_F[A,i]*exp(-Zay[A])
  Nay_M[A,i+1] = Nay_M[A-1,i]*exp(-Zay_M[A-1]) #+ Nay_M[A,i]*exp(-Zay[A])
  
  Cay1_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,2]*FM[1,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
  Cay2_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,4]*FM[2,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
  Cay3_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,6]*FM[3,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
  Cay4_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,8]*FM[4,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
  
  Cay1_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,3]*FM[1,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
  Cay2_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,5]*FM[2,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
  Cay3_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,7]*FM[3,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
  Cay4_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,8]*FM[4,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
  
  
  
  Cay1_F[A,i]=Nay_F[A,i]*((Sel[A,2]*FM[1,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
  Cay2_F[A,i]=Nay_F[A,i]*((Sel[A,4]*FM[2,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
  Cay3_F[A,i]=Nay_F[A,i]*((Sel[A,6]*FM[3,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
  Cay4_F[A,i]=Nay_F[A,i]*((Sel[A,8]*FM[4,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
  
  Cay1_M[A,i]=Nay_M[A,i]*((Sel[A,3]*FM[1,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
  Cay2_M[A,i]=Nay_M[A,i]*((Sel[A,5]*FM[2,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
  Cay3_M[A,i]=Nay_M[A,i]*((Sel[A,7]*FM[3,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
  Cay4_M[A,i]=Nay_M[A,i]*((Sel[A,8]*FM[4,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
  
  for(j in 1:A){
    #------------------------------------------------------------ LNORM --------------------------------------------------- .
    # NOTE: ASSUME FOR THIS QUICK EXAMPLE THAT SELECTIVITY = 1 FOR ALL AGES. 
    Iay_F[j,i] = max(q_survey * Nay_F[j,i] +  rnorm(1, 0, (q_survey *  Nay_F[j,i])*CV),0) #LOTS OF UNCERTAINTY!based off CV;
    Iay_M[j,i] = max(q_survey * Nay_M[j,i] +  rnorm(1, 0, (q_survey *  Nay_M[j,i])*CV),0) #LOTS OF UNCERTAINTY!based off CV;
  }
  Iay = Iay_F + Iay_M
  
}

Iay_eq = Iay

Nay = Nay_F + Nay_M
Ny = apply(Nay, 2, sum)

par(mfrow=c(1,1))
EquilibriumPop = colSums(Nay_F)+colSums(Nay_M)
plot(1:years, EquilibriumPop,type='l', ylim=c(0,max(EquilibriumPop)+1500))
lines(1:years, EquilibriumPop, type='l', col='blue')
abline(v=51, lwd=2, col="red")
abline(v=90, lwd=2, col="red")
lines(1:years, colSums(Nay),type='l',col='blue')

# Iay
par(mfrow=c(1,1))
plot(1:years, apply(Iay_eq, 2, sum), type='l',ylim=c(0,100))
abline(v=51, lwd=2, col="red")
abline(v=90, lwd=2, col="red")
lines(1:years, EquilibriumPop/100, type='l', col='blue')

Cay1_F + Cay2_F + Cay3_F + Cay4_F + Cay1_M + Cay2_M + Cay3_M + Cay4_M
Nay_F+ Nay_M

# sum(Cay1_F[,999]+Cay1_M[,999])
# sum(Cay2_F[,999]+Cay2_M[,999])
# sum(Cay3_F[,999]+Cay3_M[,999])
# sum(Cay4_F[,999]+Cay4_M[,999])

# Iay_eq_F0 = Iay_eq
# Iay_eq_F1 = Iay_eq
# Iay_eq_F2 = Iay_eq
# Iay_eq_Real = Iay_eq
# save(Iay_eq_F0, file = "D:\\vspace1\\DFA_Simulation\\SB\\Iay_eq_F0.RData")
# save(Iay_eq_F1, file = "D:\\vspace1\\DFA_Simulation\\SB\\Iay_eq_F1.RData")
# save(Iay_eq_F2, file = "D:\\vspace1\\DFA_Simulation\\SB\\Iay_eq_F2.RData")
# save(Iay_eq_Real, file = "D:\\vspace1\\DFA_Simulation\\SB\\Iay_eq_Real.RData")


##### TRIALS ################
CV1=0.5
CV2=0.5
CV3=0.5
CV4=0.5
CV5=0.5
CV6=0.5
CV7=0.5

q_survey1 = rep(0.01, length=years)
q_survey2 = rep(0.01, length=years)
q_survey3 = rep(0.01, length=years)
q_survey4 = rep(0.01, length=years)
q_survey5 = rep(0.01, length=years)
q_survey6 = rep(0.01, length=years)
q_survey7 = rep(0.01, length=years)


q0 = rep(0.01, 100)                                                                     # q0 = constant
q1 = rep(c(0.01, 0.02), 50)                                                             # q1 = period 2
q2 = rep(c(0.01, 0.02, 0.015), length=100)                                              # q2 = period 3
q3 = c(rep(0.01, 50), rep(seq(0.01, 0.045, by=0.005), each = 5), rep(0.045, 10))        # q3 = step increase (step=5)
q4 = c(rep(0.05, 50), rep(seq(0.05, 0.0325, by=-0.0025), each = 5), rep(0.0325, 10))    # q4 = step decrease (step=5)
q5 = c(rep(seq(0.03, 0.05, by=0.005), seq(0.05, 0.03, by=-0.005), length=100))          #q5 = increase then decrease
q6 = c(rep(0.01, 50), seq(0.01, 0.049, by=0.001), rep(0.049,10))                        # q6 = increase
q7 = c(rep(0.049, 50), seq(0.049, 0.01, by=-0.001), rep(0.01,10))                       # q7 = decrease
set.seed(430)
q8 = runif(100, 0.03, 0.05)                                                             # q8 = random between 0.01 and 0.05
set.seed(430)
q9 = jitter(q6, amount=0.005)                                                           # q9 = jitter increase
set.seed(431)
q10 = jitter(q7, amount=0.005)                                                          # q10 = jitter decrease
set.seed(432)
q11 = jitter(q1, amount=0.005)                                                         # q11 = jitter period 2
set.seed(433)
q12 = jitter(q2, amount=0.005)                                                         # q12 = jitter period 3                                                         # jitter decrease
set.seed(434)
q13 = jitter(q0, amount=0.005)                                                         # q13 = jitter constant
set.seed(435)
q14 = jitter(q0, amount=0.003)                                                         # q14 = jitter constant 2
set.seed(436) 
q15 = jitter(q0, amount=0.002)                                                         # q15 = jitter constant 3
set.seed(437) 
q16 = jitter(q0, amount=0.004)                                                         # q16 = jitter constant 3
set.seed(438) 
q17 = jitter(q0, amount=0.003)                                                         # q17 = jitter constant 4
set.seed(439) 
q18 = jitter(q0, amount=0.002)                                                         # q18 = jitter constant 5
set.seed(440) 
q19 = jitter(q0, amount=0.004)                                                         # q19 = jitter constant 6
set.seed(441) 
q20 = jitter(q0, amount=0.003)                                                         # q20 = jitter constant 7
set.seed(441) 
q21 = jitter(q0, amount=0.002)                                                         # q21 = jitter constant 8

# q_survey = c(rep(0.00025, 49), seq(0.0004, 0.0014, by=0.0001), rep(0.0015, 40))
# q_survey = c(rep(0.002, 49), seq(0.0018, 0.0006, by=-0.0002), rep(0.0005, 44))
# q_survey1 = c(rep(0.001, length=years/5), rep(0.0015, length=years/5),rep(0.002, length=years/5),rep(0.0025, length=years/5),rep(0.003, length=years/5))
# q_survey2 = c(rep(0.00005, length=years/2), rep(0.003, length=years/2))
# q_survey3 = c(rep(0.003, length=years/2), rep(0.00001, length=years/2))
# q_survey3 = c(rep(0.003, length=years))


Prop_F = c(0.2, 0.044, 0.754, 0.002)
sum(Prop_F)
Tot_F = as.vector(c(rep(0,45), rep(0.1, 10), rep(0.3, 10), rep(0.2, 10), rep(0.05, 25))) # F1
F_const1 = Tot_F * Prop_F[1]
F_const2 = Tot_F * Prop_F[2]
F_const3 = Tot_F * Prop_F[3]
F_const4 = Tot_F * Prop_F[4]
F_const = rbind(F_const1, F_const2, F_const3, F_const4)
apply(F_const, 2, FUN=sum)
F_var = F_const
FM = F_const

years = 100

# Survey 1 = VIMS
# Survey 2 = LPS
# Survey 3 = PLLOP
# Survey 4 = NMFS LL SE
# Survey 5 = CST NE LL
# Survey 6 = BLLOP 1
# Survey 7 = NMFS NE

CV1=0.38
CV2=0.48
CV3=0.65
CV4=0.24
CV5=0.30
CV6=0.36
CV7=0.40


Tot_F = F1

q_survey1 = q0
q_survey2 = q10
q_survey3 = q0
q_survey4 = q0
q_survey5 = q0
q_survey6 = q0
q_survey7 = q0
Iay_eq = Iay_eq_F1
Niters=100


########################## SIMULATION FUNCTION ##########################################################################################################################################
DFA_Sim = function(Tot_F, CV1, CV2, CV3, CV4, CV5, CV6, CV7, 
                   q_survey1, q_survey2, q_survey3, q_survey4, q_survey5, q_survey6,  q_survey7, Iay_eq = Iay_eq, Niters=100, c=5){
  
  # Inputs: 
  M = M_const
  # GET F
  Prop_F = c(0.2, 0.044, 0.754, 0.002)
  F_const1 = Tot_F * Prop_F[1]
  F_const2 = Tot_F * Prop_F[2]
  F_const3 = Tot_F * Prop_F[3]
  F_const4 = Tot_F * Prop_F[4]
  FM = rbind(F_const1, F_const2, F_const3, F_const4)
  My_LFreqs = function(data, multC = 100, Lbins = seq(30, 260, by=2), Lyrs=c(51:90), Linf=Linf_F, K=K_F, t0=t0_F){
    data_N = cbind(ages, round(data*multC))
    LF = data.frame(Lbins)
    for( l in Lyrs){
      df = as.data.frame(data_N[,c(1,l+1)])
      colnames(df) = c("ages","count") 
      df.expanded = df[rep(row.names(df), df$count),1:2]
      vonBert = rbvn(nrow(df.expanded),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expanded$L = vonBert$Linf * ( 1 - exp(-(vonBert$K) * (jitter(df.expanded$ages, amount=0.5) - t0) ))
      hst=hist(df.expanded$L, breaks=c( seq(30, 262, by=2)), plot=F)
      LF = cbind(LF, hst$counts)
    }
    colnames(LF) = c("LBins", as.character(Lyrs)) 
    return(LF)
  }
  
  
  # Iay
  par(mfrow=c(1,1))
  plot(1:years, apply(Iay_eq, 2, sum), type='l',ylim=c(0,max(apply(Iay_eq, 2, sum),na.rm=T)*3))
  abline(v=51, lwd=2, col="red")
  abline(v=90, lwd=2, col="red")
  # plot(1:years, apply(Iay, 2, sum), type='l')
  # Iy = apply(Iay, 2, sum)
  
  
  ####### ADD STOCHASTICITY ###########
  Iay1_F = matrix(nrow=A, ncol=years)
  Iay2_F = matrix(nrow=A, ncol=years)
  Iay3_F = matrix(nrow=A, ncol=years)
  Iay4_F = matrix(nrow=A, ncol=years)
  Iay5_F = matrix(nrow=A, ncol=years)
  Iay6_F = matrix(nrow=A, ncol=years)
  Iay7_F = matrix(nrow=A, ncol=years)
  Iay1_M = matrix(nrow=A, ncol=years)
  Iay2_M = matrix(nrow=A, ncol=years)
  Iay3_M = matrix(nrow=A, ncol=years)
  Iay4_M = matrix(nrow=A, ncol=years)
  Iay5_M = matrix(nrow=A, ncol=years)
  Iay6_M = matrix(nrow=A, ncol=years)
  Iay7_M = matrix(nrow=A, ncol=years)
  Iay1_CV = matrix(nrow=A, ncol=years)
  Iay2_CV = matrix(nrow=A, ncol=years)
  Iay3_CV = matrix(nrow=A, ncol=years)
  Iay4_CV = matrix(nrow=A, ncol=years)
  Iay5_CV = matrix(nrow=A, ncol=years)
  Iay6_CV = matrix(nrow=A, ncol=years)
  Iay7_CV = matrix(nrow=A, ncol=years)
  Iy1 = matrix(nrow=Niters, ncol=years)
  Iy2 = matrix(nrow=Niters, ncol=years)
  Iy3 = matrix(nrow=Niters, ncol=years)
  Iy4 = matrix(nrow=Niters, ncol=years)
  Iy5 = matrix(nrow=Niters, ncol=years)
  Iy6 = matrix(nrow=Niters, ncol=years)
  Iy7 = matrix(nrow=Niters, ncol=years)
  Iy1_CV = matrix(nrow=Niters, ncol=years)
  Iy2_CV = matrix(nrow=Niters, ncol=years)
  Iy3_CV = matrix(nrow=Niters, ncol=years)
  Iy4_CV = matrix(nrow=Niters, ncol=years)
  Iy5_CV = matrix(nrow=Niters, ncol=years)
  Iy6_CV = matrix(nrow=Niters, ncol=years)
  Iy7_CV = matrix(nrow=Niters, ncol=years)
  Cy_1 = matrix(nrow=Niters, ncol=years)
  Cy_2 = matrix(nrow=Niters, ncol=years)
  Cy_3 = matrix(nrow=Niters, ncol=years)
  Cy_4 = matrix(nrow=Niters, ncol=years)
  Ny_F = matrix(nrow=Niters, ncol=years)
  Ny_M = matrix(nrow=Niters, ncol=years)
  Cy = matrix(nrow=Niters, ncol=years)
  Ny = matrix(nrow=Niters, ncol=years)
  
  
  M0 = matrix(nrow=Niters, ncol=years)
  Npups = matrix(nrow=Niters, ncol=years)
  FSS = matrix(nrow=Niters, ncol=years) # Female spawning stock 
  
  ResultsList = list()

  
  Z0 = 1.085716
  Zmin = 0.1604
  zfrac = 0.8522633
  beta = 0.3582577
  Npups_eq = 2961.558 #Npups0

  
  
  set.seed(430)
  
  # F_const = rep(0, yrs)
  
  for(k in 1:Niters){
    # paste("Iteration",k,sep=" ")
    # Create matrices to store results
    Nay_F = matrix(nrow=A, ncol=years)
    Nay_M = matrix(nrow=A, ncol=years)
    Cay1_F = matrix(nrow=A, ncol=years)
    Cay2_F = matrix(nrow=A, ncol=years)
    Cay3_F = matrix(nrow=A, ncol=years)
    Cay4_F = matrix(nrow=A, ncol=years)
    Cay1_M = matrix(nrow=A, ncol=years)
    Cay2_M = matrix(nrow=A, ncol=years)
    Cay3_M = matrix(nrow=A, ncol=years)
    Cay4_M = matrix(nrow=A, ncol=years)


    
    # Set up equilibrium conditions
    Nay_F[,1] = Equilibrium_Nay_F
    Nay_M[,1] = Equilibrium_Nay_M
    
    for (i in 1:(years-1)){
      # i=i+1
      
      matjit = rnorm(length(mat_F),mat_F,mat_F*0.01)
      matjit = ifelse(matjit>1, 1, matjit)
      FSS[k,i] = sum(Nay_F[,i]*matjit) #Fecund Stock size
      Npups[k,i] = sum((Nay_F[,i]/2.5)*matjit*rnorm(length(fec), fec, fec*0.1)) # 
      M0[k,i] = rnorm(1, ( ( (1- (Npups[k,i]/Npups_eq)^beta ) * (Zmin - Z0) ) + Z0 ) , 0.1)
      # M0[i] = ( (1- (Npups[i]/Npups_eq)^beta ) * (Zmin - Z0) ) + Z0
      # M0[i] = 1.9044160 # THIS IS EQUILIBRIUM!!!
      
      
      Nay_F[1,i+1] = 0.5 * Npups[k,i]*exp(-M0[k,i])  # pups produced that survive to age 1; Recruits
      Nay_M[1,i+1] = 0.5 * Npups[k,i]*exp(-M0[k,i])  # pups produced that survive to age 1; Recruits
      
      
      Zay_F = vector(length=length(M))
      Zay_M = vector(length=length(M))
      for (j in 1:A){
        Zay_F[j] = M[j] + (Sel[j,2]*FM[1,i]) + (Sel[j,4]*FM[2,i]) + (Sel[j,6]*FM[3,i]) + (Sel[j,8]*FM[4,i])
        Zay_M[j] = M[j] + (Sel[j,3]*FM[1,i]) + (Sel[j,5]*FM[2,i]) + (Sel[j,7]*FM[3,i]) + (Sel[j,8]*FM[4,i])
      } # END J LOOP
      
      for (j in 1:(A-2)){
        Mjit = rnorm(length(M), M, M*0.1) 
        Nay_F[j+1,i+1] = Nay_F[j,i]*exp(-Zay_F[j]) 
        Nay_M[j+1,i+1] = Nay_M[j,i]*exp(-Zay_M[j]) 
        
        
        Cay1_F[j,i] = Nay_F[j,i]*((Sel[j,2]*FM[1,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
        Cay2_F[j,i] = Nay_F[j,i]*((Sel[j,4]*FM[2,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
        Cay3_F[j,i] = Nay_F[j,i]*((Sel[j,6]*FM[3,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
        Cay4_F[j,i] = Nay_F[j,i]*((Sel[j,8]*FM[4,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
        
        Cay1_M[j,i] = Nay_M[j,i]*((Sel[j,3]*FM[1,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
        Cay2_M[j,i] = Nay_M[j,i]*((Sel[j,5]*FM[2,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
        Cay3_M[j,i] = Nay_M[j,i]*((Sel[j,7]*FM[3,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
        Cay4_M[j,i] = Nay_M[j,i]*((Sel[j,8]*FM[4,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
        
        
      } # END J LOOP
      
      Nay_F[A,i+1] = Nay_F[A-1,i]*exp(-Zay_F[A-1]) #+ Nay_F[A,i]*exp(-Zay[A])
      Nay_M[A,i+1] = Nay_M[A-1,i]*exp(-Zay_M[A-1]) #+ Nay_M[A,i]*exp(-Zay[A])
      
      Cay1_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,2]*FM[1,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      Cay2_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,4]*FM[2,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      Cay3_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,6]*FM[3,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      Cay4_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,8]*FM[4,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      
      Cay1_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,3]*FM[1,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
      Cay2_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,5]*FM[2,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
      Cay3_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,7]*FM[3,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
      Cay4_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,8]*FM[4,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
      
      
      
      Cay1_F[A,i]=Nay_F[A,i]*((Sel[A,2]*FM[1,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      Cay2_F[A,i]=Nay_F[A,i]*((Sel[A,4]*FM[2,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      Cay3_F[A,i]=Nay_F[A,i]*((Sel[A,6]*FM[3,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      Cay4_F[A,i]=Nay_F[A,i]*((Sel[A,8]*FM[4,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      
      Cay1_M[A,i]=Nay_M[A,i]*((Sel[A,3]*FM[1,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
      Cay2_M[A,i]=Nay_M[A,i]*((Sel[A,5]*FM[2,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
      Cay3_M[A,i]=Nay_M[A,i]*((Sel[A,7]*FM[3,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
      Cay4_M[A,i]=Nay_M[A,i]*((Sel[A,8]*FM[4,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
      
      
      for(h in 1:A){
        # h=h+1
        Iay1_CV[h,i] = runif(1, min=CV1-0.1, max=CV1+0.1)
        Iay1_F[h,i] = max(q_survey1[i]*Sel[h,17] * Nay_F[h,i] * 
                            exp(rnorm(1, 0, (q_survey1[i]*Sel[h,17]*Nay_F[h,i])*Iay1_CV[h,i]) - 
                                  ( ((q_survey1[i]*Sel[h,17]*Nay_F[h,i]*Iay1_CV[h,i])^2)/2) ), 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay1_M[h,i] = max(q_survey1[i]*Sel[h,18] * Nay_M[h,i] * 
                            exp(rnorm(1, 0, (q_survey1[i]*Sel[h,18]*Nay_M[h,i])*Iay1_CV[h,i]) - 
                                  ( ((q_survey1[i]*Sel[h,18]*Nay_M[h,i]*Iay1_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        
        Iay2_CV[h,i] = runif(1, min=CV2-0.1, max=CV2+0.1)
        Iay2_F[h,i] = max(q_survey2[i]*Sel[h,9] * Nay_F[h,i] *
                            exp(rnorm(1, 0, (q_survey2[i]*Sel[h,9]*Nay_F[h,i])*Iay2_CV[h,i]) -
                                  ( ((q_survey2[i]*Sel[h,9]*Nay_F[h,i]*Iay2_CV[h,i])^2)/2) ), 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay2_M[h,i] = max(q_survey2[i]*Sel[h,10] * Nay_M[h,i] *
                            exp(rnorm(1, 0, (q_survey2[i]*Sel[h,10]*Nay_M[h,i])*Iay2_CV[h,i]) -
                                  ( ((q_survey2[i]*Sel[h,10]*Nay_M[h,i]*Iay2_CV[h,i])^2)/2) ), 0)
        
        Iay3_CV[h,i] = runif(1, min=CV3-0.1, max=CV3+0.1)
        Iay3_F[h,i] = max(q_survey3[i]*Sel[h,28] * Nay_F[h,i] *
                            exp(rnorm(1, 0, (q_survey3[i]*Sel[h,28]*Nay_F[h,i])*Iay3_CV[h,i]) -
                                  ( ((q_survey3[i]*Sel[h,28]*Nay_F[h,i]*Iay3_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay3_M[h,i] = max(q_survey3[i]*Sel[h,29] * Nay_M[h,i] *
                            exp(rnorm(1, 0, (q_survey3[i]*Sel[h,29]*Nay_M[h,i])*Iay3_CV[h,i]) -
                                  ( ((q_survey3[i]*Sel[h,29]*Nay_M[h,i]*Iay3_CV[h,i]) ^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        
        Iay4_CV[h,i] = runif(1, min=CV4-0.1, max=CV4+0.1)
        Iay4_F[h,i] = max(q_survey4[i]*Sel[h,19] * Nay_F[h,i] *
                            exp(rnorm(1, 0, (q_survey4[i]*Sel[h,19]*Nay_F[h,i])*Iay4_CV[h,i]) -
                                  ( ((q_survey4[i]*Sel[h,19]*Nay_F[h,i]*Iay4_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay4_M[h,i] = max(q_survey4[i]*Sel[h,21] * Nay_M[h,i] * 
                            exp(rnorm(1, 0, (q_survey4[i]*Sel[h,21]*Nay_M[h,i])*Iay4_CV[h,i]) -
                                  ( ((q_survey4[i]*Sel[h,21]*Nay_M[h,i]*Iay4_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        
        Iay5_CV[h,i] = runif(1, min=CV5-0.1, max=CV5+0.1)
        Iay5_F[h,i] = max(q_survey5[i]*Sel[h,23] * Nay_F[h,i] *
                            exp(rnorm(1, 0, (q_survey5[i]*Sel[h,23]*Nay_F[h,i])*Iay5_CV[h,i]) -
                                  ( ((q_survey5[i]*Sel[h,23]*Nay_F[h,i]*Iay5_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay5_M[h,i] = max(q_survey5[i]*Sel[h,24] * Nay_M[h,i] *
                            exp(rnorm(1, 0, (q_survey5[i]*Sel[h,24]*Nay_M[h,i])*Iay5_CV[h,i]) -
                                  ( ((q_survey5[i]*Sel[h,24]*Nay_M[h,i]*Iay5_CV[h,i])^2 )/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        
        Iay6_CV[h,i] = runif(1, min=CV6-0.1, max=CV6+0.1)
        Iay6_F[h,i] = max(q_survey6[i]*Sel[h,13] * Nay_F[h,i] *
                            exp(rnorm(1, 0, (q_survey6[i]*Sel[h,13]*Nay_F[h,i])*Iay6_CV[h,i]) -
                                  ( ((q_survey6[i]*Sel[h,13]*Nay_F[h,i]*Iay6_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay6_M[h,i] = max(q_survey6[i]*Sel[h,12] * Nay_M[h,i] *
                            exp(rnorm(1, 0, (q_survey6[i]*Sel[h,12]*Nay_M[h,i])*Iay6_CV[h,i]) -
                                  ( ((q_survey6[i]*Sel[h,12]*Nay_M[h,i]*Iay6_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        
        Iay7_CV[h,i] = runif(1, min=CV7-0.1, max=CV7+0.1)
        Iay7_F[h,i] = max(q_survey7[i]*Sel[h,27] * Nay_F[h,i] *
                            exp(rnorm(1, 0, (q_survey7[i]*Sel[h,27]*Nay_F[h,i])*Iay7_CV[h,i]) -
                                  ( ((q_survey7[i]*Sel[h,27]*Nay_F[h,i]*Iay7_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay7_M[h,i] = max(q_survey7[i]*Sel[h,26] * Nay_M[h,i] *
                            exp(rnorm(1, 0, (q_survey7[i]*Sel[h,26]*Nay_M[h,i])*Iay7_CV[h,i]) -
                                  ( ((q_survey7[i]*Sel[h,26]*Nay_M[h,i]*Iay7_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;} # end h
      } # end h
      
      
    }# ends i loop
    
    
    # USE VON BERT CURVE TO ESTIMATE LENGTH COMP. Let von Bert parameters vary by Niter...???? 
    #   year (implement in i loop)... 
    #   Or should they vary by year class? Could create a matrix of von bert parameters based on age class 
    #       and pull out appropriate values as needed
    # Lt = Linf * (1 - exp(-K*(t - t0)))
    # For sexes combined.

    Linf_M=175.5
    K_M=0.143
    t0_M=-2.388
    Linf_F=183.3
    K_F=0.124
    t0_F=-3.098
    # L0 = Linf * (1 - exp(-K*(0-t0)))
    
    
    #####  GET LENGTH FREQUENCIES #####
    
    LFreqList = list()

    ## Female Fishery LFreqs ##
    LFreqList[[paste("LF_Cay1F",k,sep="_")]] = My_LFreqs(data=Cay1_F, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90))
    LFreqList[[paste("LF_Cay2F",k,sep="_")]] = My_LFreqs(data=Cay2_F, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90))
    LFreqList[[paste("LF_Cay3F",k,sep="_")]] = My_LFreqs(data=Cay3_F, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90))
    LFreqList[[paste("LF_Cay4F",k,sep="_")]] = My_LFreqs(data=Cay4_F, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90))
    
    ## Male Fishery LFreqs ##function(data, multC = 100, Lbins = seq(30, 260, by=2), Lyrs=c(51:90), Linf=Linf_F, K=K_F, t0=t0_F)
    LFreqList[[paste("LF_Cay1M",k,sep="_")]] = My_LFreqs(data=Cay1_M, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Cay2M",k,sep="_")]] = My_LFreqs(data=Cay2_M, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Cay3M",k,sep="_")]] = My_LFreqs(data=Cay3_M, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Cay4M",k,sep="_")]] = My_LFreqs(data=Cay4_M, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    
    ## Female Survey LFreqs ##
    LFreqList[[paste("LF_Iay1F",k,sep="_")]] = My_LFreqs(data=Iay1_F, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:53,59:64,67:90))
    LFreqList[[paste("LF_Iay2F",k,sep="_")]] = My_LFreqs(data=Iay2_F, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(58:90))
    LFreqList[[paste("LF_Iay3F",k,sep="_")]] = My_LFreqs(data=Iay3_F, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(64:90))
    LFreqList[[paste("LF_Iay4F",k,sep="_")]] = My_LFreqs(data=Iay4_F, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(67:90))
    LFreqList[[paste("LF_Iay5F",k,sep="_")]] = My_LFreqs(data=Iay5_F, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(73:90))
    LFreqList[[paste("LF_Iay6F",k,sep="_")]] = My_LFreqs(data=Iay6_F, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(66:79))
    LFreqList[[paste("LF_Iay7F",k,sep="_")]] = My_LFreqs(data=Iay7_F, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(seq(68, 89, by=3)))
    
    ## Male Survey LFreqs ##
    LFreqList[[paste("LF_Iay1M",k,sep="_")]] = My_LFreqs(data=Iay1_M, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:53,59:64,67:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Iay2M",k,sep="_")]] = My_LFreqs(data=Iay2_M, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(58:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Iay3M",k,sep="_")]] = My_LFreqs(data=Iay3_M, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(64:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Iay4M",k,sep="_")]] = My_LFreqs(data=Iay4_M, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(67:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Iay5M",k,sep="_")]] = My_LFreqs(data=Iay5_M, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(73:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Iay6M",k,sep="_")]] = My_LFreqs(data=Iay6_M, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(66:79),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Iay7M",k,sep="_")]] = My_LFreqs(data=Iay7_M, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(seq(68, 89, by=3)),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    
    
    save(LFreqList, file=paste("LFreqList_",k,".RData",sep=""))
    
    #### END LFREQS #####
    
    
    
    Iy1[k,] = apply(Iay1_F, 2, sum) + apply(Iay1_M, 2, sum)
    Iy2[k,] = apply(Iay2_F, 2, sum) + apply(Iay2_M, 2, sum)
    Iy3[k,] = apply(Iay3_F, 2, sum) + apply(Iay3_M, 2, sum)
    Iy4[k,] = apply(Iay4_F, 2, sum) + apply(Iay4_M, 2, sum)
    Iy5[k,] = apply(Iay5_F, 2, sum) + apply(Iay5_M, 2, sum)
    Iy6[k,] = apply(Iay6_F, 2, sum) + apply(Iay6_M, 2, sum)
    Iy7[k,] = apply(Iay7_F, 2, sum) + apply(Iay7_M, 2, sum)
    
    lines(1:years, Iy1[k,], type='l', col='grey')
    lines(1:years, Iy2[k,], type='l', col='blue')
    lines(1:years, Iy3[k,], type='l', col='red')
    lines(1:years, Iy4[k,], type='l', col='green')
    lines(1:years, Iy5[k,], type='l', col='orange')
    lines(1:years, Iy6[k,], type='l', col='mediumorchid2')
    lines(1:years, Iy7[k,], type='l', col='cyan')
    Iy1_CV[k,] = apply(Iay1_CV, 2, mean, na.rm=T)
    Iy2_CV[k,] = apply(Iay2_CV, 2, mean, na.rm=T)
    Iy3_CV[k,] = apply(Iay3_CV, 2, mean, na.rm=T)
    Iy4_CV[k,] = apply(Iay4_CV, 2, mean, na.rm=T)
    Iy5_CV[k,] = apply(Iay5_CV, 2, mean, na.rm=T)
    Iy6_CV[k,] = apply(Iay6_CV, 2, mean, na.rm=T)
    Iy7_CV[k,] = apply(Iay7_CV, 2, mean, na.rm=T)
    
    Cy_1[k,] = apply(Cay1_F, 2, sum, na.rm=T) + apply(Cay1_M, 2, sum, na.rm=T) 
    Cy_2[k,] = apply(Cay2_F, 2, sum, na.rm=T) + apply(Cay2_M, 2, sum, na.rm=T) 
    Cy_3[k,] = apply(Cay3_F, 2, sum, na.rm=T) + apply(Cay3_M, 2, sum, na.rm=T) 
    Cy_4[k,] = apply(Cay4_F, 2, sum, na.rm=T) + apply(Cay4_M, 2, sum, na.rm=T) 
    Cy[k,] = Cy_1[k,] + Cy_2[k,] + Cy_3[k,] + Cy_4[k,]
    
    Ny_F[k,] = apply(Nay_F, 2, sum, na.rm=T)
    Ny_M[k,] = apply(Nay_M, 2, sum, na.rm=T)
    Ny[k,] = Ny_F[k,] + Ny_M[k,]
    # lines(1:years, Ny[k,]/100, type='l', lwd=2)
    
    
    #write 
    AnnualNC = list()
    # write.csv(Nay_F, paste(paste("NayF", k, sep="_"),".csv",sep=""))
    # write.csv(Nay_M, paste(paste("NayM", k, sep="_"),".csv",sep=""))
    assign(paste("NayF", k, sep="_"), Nay_F)
    assign(paste("NayM", k, sep="_"), Nay_M)
    AnnualNC[[paste("NayF",k,sep="_")]] = Nay_F
    AnnualNC[[paste("NayM",k,sep="_")]] = Nay_M
    
    # write.csv(Cay_F, paste(paste("CayF", k, sep="_"),".csv",sep=""))
    # write.csv(Cay_M, paste(paste("CayM", k, sep="_"),".csv",sep=""))
    assign(paste("Cay1F", k, sep="_"), Cay1_F)
    assign(paste("Cay2F", k, sep="_"), Cay2_F)
    assign(paste("Cay3F", k, sep="_"), Cay3_F)
    assign(paste("Cay4F", k, sep="_"), Cay4_F)
    assign(paste("Cay1M", k, sep="_"), Cay1_M)
    assign(paste("Cay2M", k, sep="_"), Cay2_M)
    assign(paste("Cay3M", k, sep="_"), Cay3_M)
    assign(paste("Cay4M", k, sep="_"), Cay4_M)
    AnnualNC[[paste("Cay1F",k,sep="_")]] = Cay1_F
    AnnualNC[[paste("Cay2F",k,sep="_")]] = Cay2_F
    AnnualNC[[paste("Cay3F",k,sep="_")]] = Cay3_F
    AnnualNC[[paste("Cay4F",k,sep="_")]] = Cay4_F
    AnnualNC[[paste("Cay1M",k,sep="_")]] = Cay1_M
    AnnualNC[[paste("Cay2M",k,sep="_")]] = Cay2_M
    AnnualNC[[paste("Cay3M",k,sep="_")]] = Cay3_M
    AnnualNC[[paste("Cay4M",k,sep="_")]] = Cay4_M

    save(AnnualNC, file=paste("AnnualNC_",k,".RData",sep=""))
    
    # write.csv() ########### WRITE CSVs to SAVE Nay AND Cay!!!
    
  } # end k loop
  
  
  # 
  
  ResultsList[["Iy1"]] = Iy1
  ResultsList[["Iy2"]] = Iy2
  ResultsList[["Iy3"]] = Iy3
  ResultsList[["Iy4"]] = Iy4
  ResultsList[["Iy5"]] = Iy5
  ResultsList[["Iy6"]] = Iy6
  ResultsList[["Iy7"]] = Iy7
  ResultsList[["Iy1_CV"]] = Iy1_CV
  ResultsList[["Iy2_CV"]] = Iy2_CV
  ResultsList[["Iy3_CV"]] = Iy3_CV
  ResultsList[["Iy4_CV"]] = Iy4_CV
  ResultsList[["Iy5_CV"]] = Iy5_CV
  ResultsList[["Iy6_CV"]] = Iy6_CV
  ResultsList[["Iy7_CV"]] = Iy7_CV
  
  ResultsList[["M0"]] = M0
  ResultsList[["FSS"]] = FSS
  ResultsList[["Npups"]] = Npups
  
  ResultsList[["Cy"]] = Cy
  ResultsList[["Cy_1"]] = Cy_1
  ResultsList[["Cy_2"]] = Cy_2
  ResultsList[["Cy_3"]] = Cy_3
  ResultsList[["Cy_4"]] = Cy_4
  
  ResultsList[["Ny"]] = Ny
  
  save(ResultsList, file="ResultsList.RData")
  
  
  
  
  ########## TEST DFA #############
  # 51 - 90; 
  # library(MARSS)
  yrs = 51:90
  biasAll = matrix(nrow=Niters, ncol=length(yrs))
  bias = vector(length=Niters)
  RMSE = vector(length=Niters)
  biasAll.BT = matrix(nrow=Niters, ncol=length(yrs))
  RMSE.BT = vector(length=Niters)
 
  FitRatios = matrix(nrow=Niters, ncol=7)
  FactorLoadings=vector()
  DFATrends = vector()
  DFATrendsBT = vector()
  DFATrendsSE = vector()
  DFATrendsSEBT = vector()
  upCI_DFATrends = vector()
  lowCI_DFATrends = vector()
  
  datz_SD = matrix(nrow=Niters, ncol=7)
  
  # i=1
  
  for(i in 1:Niters){
    I1 = c(Iy1[i,51:53], rep(NA, 5), Iy1[i,59:64], NA, NA, Iy1[i,67:90])
    I2 = c(rep(NA, 7), Iy2[i,58:90])
    I3 = c(rep(NA, 13), Iy3[i,64:90])
    I4 = c(rep(NA, 16), Iy4[i,67:90])
    I5 = c(rep(NA, 22), Iy5[i,73:90])
    I6 = c(rep(NA, 15), Iy6[i,66:79], rep(NA, 11))
    I7 = c(rep(NA, 17), Iy7[i,68], NA, NA, Iy7[i,71], NA, NA, Iy7[i,74], NA, NA, Iy7[i,77], NA, NA, Iy7[i,80], 
           NA, NA, Iy7[i,83], NA, NA, Iy7[i,86], NA, NA, Iy7[i,89], NA)
    
    assign(paste("dat",i,sep=""), rbind(I1, I2, I3, I4, I5, I6, I7))
    # assign(paste("dat",i,sep=""), rbind(Iy1[i,51:90], Iy2[i,51:90], Iy3[i,51:90], Iy4[i,51:90], Iy5[i,51:90], Iy6[i,51:90], Iy7[i,51:90]))
    
    dat.a = get(paste("dat",i,sep=""))
    # c = c(30, 10, 2, 1, 25, 2, 1)
    # c = c(100, 15, 2, 0.7, 40, 1.5, 0.6)
    # c = c(1, 1, 1, 1, 1, 1, 1)
    # c = c(12, 3, 1, 0.4, 11, 0.8, 0.35)
    # c = c(12, 3, 200, 0.4, 200, 0.8, 250)
    
    # c = c(12, 3, 1, 0.4, 11, 0.8, 0.35)
    # c = c(8, 20, 0.7, 0.45, 11, 0.8, 0.6)
    # c=c(29, 2.9, 21, 0.4, 2100, 1.6, 2) #40
    # c=c(29, 2.9, 21, 0.4, 2100, 1.6, 2) #41
    # c = c(45, 4.0, 5, 4, 2000, 0.3, 1.2) #40
    
    dat = dat.a * c
    
    TT=ncol(dat)
    N.ts = nrow(dat)
    # Standardize data
    datL=log(dat)
    y.bar = apply(datL, 1, mean, na.rm=TRUE)
    dat.dm = (datL-y.bar) / y.bar
    gsd = sd(c(dat.dm), na.rm=T)
    dat.z = dat.dm/gsd
    y.bar
    gsd
    datz_SD[i,] = apply(dat.z, 1, sd, na.rm=T)
    datz_SD[i,]
    min(datL, na.rm=T)
    
    dat.z = as.matrix(dat.z)
    rownames(dat.z) = c("Survey1","Survey2","Survey3","Survey4","Survey5","Survey6","Survey7")
    
    
    
    
    ##### SB DELTA LOGNORMAL ####
    cntl.list = list(minit=200, maxit=500000, abstol=0.02, allow.degen=FALSE, conv.test.slope.tol = 0.05)
    R = diag(c(CV1, CV2, CV3, CV4, CV5, CV6, CV7), nrow=7, ncol=7)
    # dfa = MARSS(dat.z, model = list(m=1,R="diagonal and equal"), control=cntl.list, form="dfa", z.score=TRUE)
    dfa = MARSS(dat.z, model = list(m=1,R=R), control=cntl.list, form="dfa", z.score=FALSE)
    
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
    
    
    # plot factor loadings
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
    
    par(mfrow=c(4,2))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='Year', ylab='Abundance', ylim=c(min(lowerCI), max(upperCI))))
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, index, type='l', lwd=2, pch=16))
    title(paste("Iteration",i,sep=" "))
    abline(h=0)
    
    # Nay = get(paste("NayF", i, sep="_")) + get(paste("NayM", i, sep="_"))
    # N = apply(Nay[,yrs],2,sum)
    N = Ny[i,51:90]
    Nscale = (N-mean(N))/sd(N)
    NL = log(N)
    NLscale = (NL - mean(NL))/sd(NL)
    index.z = (index - mean(index))/sd(index)
    index.z = ifelse(is.na(index.z), 0, index.z)
    index.z = ifelse(abs(index.z)==Inf, 0, index.z)
    
    
    indexSEBT = indexSE * gsd
    assign(paste("indexSEBT",i,sep=""), indexSEBT)
    
    indexBT = exp(index*gsd + ((indexSEBT^2)/2) )
    assign(paste("indexBT",i,sep=""), indexBT)
    

    
    lines(yrs,rescale(NL,to=c(min(index),max(index))), type='l', col='red', lwd=2)
    assign(paste("N",i,sep=""), N)
    assign(paste("Nscale",i,sep=""),Nscale)
    assign(paste("NL",i,sep=""), NL)
    assign(paste("NLscale",i,sep=""),NLscale)
    # write.csv(N,paste("N",i,".csv",sep=""))
    
    
    biasAll[i,] = index.z - NLscale
    bias[i] = mean(index.z - NLscale)
    RMSE[i] = sqrt(sum((index.z-NLscale)^2)/length(index.z))
    # RMSE_yr[i] = sqrt(sum((Nscale - indexscale)^2)/length(indexscale))
    
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
      title(paste("Sandbar",survey[n], sep=" "))
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
  
  colnames(FactorLoadings) = c("FL1","FL2","FL3","FL4","FL5","FL6","FL7",
                               "lowCI_FL1","lowCI_FL2","lowCI_FL3","lowCI_FL4","lowCI_FL5","lowCI_FL6","lowCI_FL7",
                               "upCI_FL1","upCI_FL2","upCI_FL3","upCI_FL4","upCI_FL5","upCI_FL6","upCI_FL7")
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
  # write.csv(DFATrends, "DFATrends.csv")
  # write.csv(DFATrendsSE, "DFATrendsSE.csv")
  # write.csv(upCI_DFATrends, "upCI_DFATrends.csv")
  # write.csv(lowCI_DFATrends, "lowCI_DFATrends.csv")
  
  # write.csv(FLoadings, "FactorLoadings_sum.csv")
  
  # write.csv(bias, "bias.csv")
  # write.csv(biasAll, "biasAll.csv")
  # write.csv(RMSE, "RMSE.csv")
  
  
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
    text(x = min(yrs)+4, y = 0.8*max(min(lowerCI)), labels = format( FR[l], digits=3 ) , cex=1)
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
    plot(51:90, c(Iy1[l,51:53], rep(NA, 5), Iy1[l,59:64], NA, NA, Iy1[l,67:90]), type='l', col='grey', xlim=c(51, 90), axes=F,
         ylim=c(0, max(c(Iy1[l,], Iy2[l,], Iy3[l,],Iy4[l,], Iy5[l,], Iy6[l,], Iy7[l,]), na.rm=T)+0.1))
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    lines(51:90, c(rep(NA, 7), Iy2[l,58:90]), type='l', col='blue')
    lines(51:90, c(rep(NA, 13), Iy3[l,64:90]), type='l', col='red')
    lines(51:90, c(rep(NA, 16), Iy4[l,67:90]), type='l', col='green')
    lines(51:90, c(rep(NA, 22), Iy5[l,73:90]), type='l', col='orange')
    lines(51:90, c(rep(NA, 15), Iy6[l,66:79], rep(NA, 11)), type='l', col='mediumorchid2')
    points(51:90, c(rep(NA, 17), Iy7[l,68], NA, NA, Iy7[l,71], NA, NA, Iy7[l,74], NA, NA, Iy7[l,77], NA, NA, Iy7[l,80], 
                    NA, NA, Iy7[l,83], NA, NA, Iy7[l,86], NA, NA, Iy7[l,89], NA), type='l', col='cyan')
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
    text(x = min(yrs)+4, y = 0.8*min(lowerCI), labels = format( FR[l], digits=3 ) , cex=1)
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
    plot(yrs, N, type='l',col="red", axes=FALSE)
    axis(1, labels=FALSE)
    axis(2, labels=FALSE)
    lines(yrs, rescale(index, to=c(min(N),max(N))))
    text(x = min(yrs)+4, y = 1.05*min(N), labels = format( FR[l], digits=3 ) , cex=1)
    
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
    text(x = min(yrs)+4, y = 0.8*min((index-mean(index))/sd(index)), labels = format( FR[l], digits=3 ) , cex=1)
    
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
  #   text(x = min(yrs)+4, y = 0.8*min(lowerCI), labels = format( FR[l], digits=3 ) , cex=1)
  #   par(new=T)
  #   NL = get(paste("NL",l,sep=""))
  #   plot(yrs,NL, type='l', col='red', axes=FALSE, lwd=1.5)
  #   axis(4, labels=FALSE)
  # }
  # 
  # dev.off()
  
  
  
  
  
  
  
  
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
    text(x = min(yrs)+4, y = 1.2*min(NL), labels = format( FR[l], digits=3 ) , cex=1)
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
  
  return(list(RMSE=RMSE, biasAll=biasAll, FitRatios = FitRatios, RMSE.BT=RMSE.BT, biasAll.BT=biasAll.BT))} # END FUNCTION




#################
# Tot_F = as.vector(rep(0,100))
setwd("C:\\~")


Tot_F = F1
CV1=0.38; CV2=0.48; CV3=0.65; CV4=0.24; CV5=0.30; CV6=0.36; CV7=0.40
q_survey1 = q0
q_survey2 = q0
q_survey3 = q0
q_survey4 = q0
q_survey5 = q0
q_survey6 = q0
q_survey7 = q0
Iay_eq = Iay_eq_F1

#q10
#q9
# 9 & 3; then 10 & 4
# 9 & 4; then 10 & 3
# 3 9 9 ; then 4 10 10
# 4 9 9 ; then 3 10 10 
q_survey1 = q4
q_survey3 = q9
q_survey5 = q9

Niters=15

Tot_F = F1
CV1=0.38; CV2=0.48; CV3=0.65; CV4=0.24; CV5=0.30; CV6=0.36; CV7=0.40;
q_survey1 = q0
q_survey2 = q3
q_survey3 = q0
q_survey4 = q9
q_survey5 = q0
q_survey6 = q9
q_survey7 = q0
Iay_eq = Iay_eq_F1
c=c(3.5, 0.7, 4.5, 0.4, 6, 0.6, 3.5)

Tot_F = F1
CV1=0.38; CV2=0.48; CV3=0.65; CV4=0.24; CV5=0.30; CV6=0.36; CV7=0.40;
q_survey1 = q0
q_survey2 = q0;
q_survey3 = q0;
q_survey4 = q0;
q_survey5 = q10;
q_survey6 = q0;
q_survey7 = q4;
Iay_eq = Iay_eq_F1;
Niters=100;
c=c(3.5, 1.3, 4.5, 6, 10, 5.5, 2.8)


Tot_F = F1;
CV1=0.38; CV2=0.48; CV3=0.65; CV4=0.24; CV5=0.30; CV6=0.36; CV7=0.40;
q_survey1 = q0;
q_survey2 = q0;
q_survey3 = q0;
q_survey4 = q0;
q_survey5 = q0;
q_survey6 = q0;
q_survey7 = q0;
Iay_eq = Iay_eq_F1;
Niters=100;
c=c(3.5, 1.3, 4.5, 6, 6, 5.5, 3.5)

# GET SB101 Age-LENGTH information. 
Linf_M=175.5
K_M=0.143
t0_M=-2.388
Linf_F=183.3
K_F=0.124
t0_F=-3.098
data=Nay_M; multC = 1; Lbins = seq(30, 260, by=2); Lyrs=c(51:90)
Linf=Linf_M; K=K_M; t0=t0_M
# df.expandedF = df.expanded
# df.expandedM = df.expanded
ALen = list()
ALen[["df.expandedF"]] <- df.expandedF
ALen[["df.expandedM"]] <- df.expandedM
save(ALen, file="D:\\vspace1\\DFA_Simulation\\SB\\ALen_1.RData")



                         
setwd("C:/~")                      
################# FULL DATA SIMULATION FUNCTION ##############################################################################

DFA_Sim_Full = function(Tot_F, CV1, CV2, CV3, CV4, CV5, CV6, CV7, 
                   q_survey1, q_survey2, q_survey3, q_survey4, q_survey5, q_survey6,  q_survey7, Iay_eq = Iay_eq, Niters=100, c=5){
  
  # Inputs: 
  M = M_const
  # GET F
  Prop_F = c(0.2, 0.044, 0.754, 0.002)
  F_const1 = Tot_F * Prop_F[1]
  F_const2 = Tot_F * Prop_F[2]
  F_const3 = Tot_F * Prop_F[3]
  F_const4 = Tot_F * Prop_F[4]
  FM = rbind(F_const1, F_const2, F_const3, F_const4)
  My_LFreqs = function(data, multC = 100, Lbins = seq(30, 260, by=2), Lyrs=c(51:90), Linf=Linf_F, K=K_F, t0=t0_F){
    data_N = cbind(ages, round(data*multC))
    LF = data.frame(Lbins)
    for( l in Lyrs){
      df = as.data.frame(data_N[,c(1,l+1)])
      colnames(df) = c("ages","count") 
      df.expanded = df[rep(row.names(df), df$count),1:2]
      vonBert = rbvn(nrow(df.expanded),m1 = Linf,s1 = Linf*0.1, m2 = K, s2 = K*0.1, rho = -0.9)
      df.expanded$L = vonBert$Linf * ( 1 - exp(-(vonBert$K) * (jitter(df.expanded$ages, amount=0.5) - t0) ))
      hst=hist(df.expanded$L, breaks=c( seq(30, 262, by=2)), plot=F)
      LF = cbind(LF, hst$counts)
    }
    colnames(LF) = c("LBins", as.character(Lyrs)) 
    return(LF)
  }

  
  # Iay
  par(mfrow=c(1,1))
  plot(1:years, apply(Iay_eq, 2, sum), type='l',ylim=c(0,max(apply(Iay_eq, 2, sum),na.rm=T)*3))
  abline(v=51, lwd=2, col="red")
  abline(v=90, lwd=2, col="red")
  # plot(1:years, apply(Iay, 2, sum), type='l')
  # Iy = apply(Iay, 2, sum)
  
  
  ####### ADD STOCHASTICITY ###########
  Iay1_F = matrix(nrow=A, ncol=years)
  Iay2_F = matrix(nrow=A, ncol=years)
  Iay3_F = matrix(nrow=A, ncol=years)
  Iay4_F = matrix(nrow=A, ncol=years)
  Iay5_F = matrix(nrow=A, ncol=years)
  Iay6_F = matrix(nrow=A, ncol=years)
  Iay7_F = matrix(nrow=A, ncol=years)
  Iay1_M = matrix(nrow=A, ncol=years)
  Iay2_M = matrix(nrow=A, ncol=years)
  Iay3_M = matrix(nrow=A, ncol=years)
  Iay4_M = matrix(nrow=A, ncol=years)
  Iay5_M = matrix(nrow=A, ncol=years)
  Iay6_M = matrix(nrow=A, ncol=years)
  Iay7_M = matrix(nrow=A, ncol=years)
  Iay1_CV = matrix(nrow=A, ncol=years)
  Iay2_CV = matrix(nrow=A, ncol=years)
  Iay3_CV = matrix(nrow=A, ncol=years)
  Iay4_CV = matrix(nrow=A, ncol=years)
  Iay5_CV = matrix(nrow=A, ncol=years)
  Iay6_CV = matrix(nrow=A, ncol=years)
  Iay7_CV = matrix(nrow=A, ncol=years)
  Iy1 = matrix(nrow=Niters, ncol=years)
  Iy2 = matrix(nrow=Niters, ncol=years)
  Iy3 = matrix(nrow=Niters, ncol=years)
  Iy4 = matrix(nrow=Niters, ncol=years)
  Iy5 = matrix(nrow=Niters, ncol=years)
  Iy6 = matrix(nrow=Niters, ncol=years)
  Iy7 = matrix(nrow=Niters, ncol=years)
  Iy1_CV = matrix(nrow=Niters, ncol=years)
  Iy2_CV = matrix(nrow=Niters, ncol=years)
  Iy3_CV = matrix(nrow=Niters, ncol=years)
  Iy4_CV = matrix(nrow=Niters, ncol=years)
  Iy5_CV = matrix(nrow=Niters, ncol=years)
  Iy6_CV = matrix(nrow=Niters, ncol=years)
  Iy7_CV = matrix(nrow=Niters, ncol=years)
  Cy_1 = matrix(nrow=Niters, ncol=years)
  Cy_2 = matrix(nrow=Niters, ncol=years)
  Cy_3 = matrix(nrow=Niters, ncol=years)
  Cy_4 = matrix(nrow=Niters, ncol=years)
  Ny_F = matrix(nrow=Niters, ncol=years)
  Ny_M = matrix(nrow=Niters, ncol=years)
  Cy = matrix(nrow=Niters, ncol=years)
  Ny = matrix(nrow=Niters, ncol=years)
  
  
  M0 = matrix(nrow=Niters, ncol=years)
  Npups = matrix(nrow=Niters, ncol=years)
  FSS = matrix(nrow=Niters, ncol=years) # Female spawning stock 
  
  ResultsList = list()
  
  
  Z0 = 1.085716
  Zmin = 0.1604
  zfrac = 0.8522633
  beta = 0.3582577
  Npups_eq = 2961.558 #Npups0
  
  
  set.seed(430)
  
  # F_const = rep(0, yrs)
  
  for(k in 1:Niters){
    # paste("Iteration",k,sep=" ")
    # Create matrices to store results
    Nay_F = matrix(nrow=A, ncol=years)
    Nay_M = matrix(nrow=A, ncol=years)
    Cay1_F = matrix(nrow=A, ncol=years)
    Cay2_F = matrix(nrow=A, ncol=years)
    Cay3_F = matrix(nrow=A, ncol=years)
    Cay4_F = matrix(nrow=A, ncol=years)
    Cay1_M = matrix(nrow=A, ncol=years)
    Cay2_M = matrix(nrow=A, ncol=years)
    Cay3_M = matrix(nrow=A, ncol=years)
    Cay4_M = matrix(nrow=A, ncol=years)
    
    
    
    # Set up equilibrium conditions
    Nay_F[,1] = Equilibrium_Nay_F
    Nay_M[,1] = Equilibrium_Nay_M
    
    for (i in 1:(years-1)){
      # i=i+1
      
      matjit = rnorm(length(mat_F),mat_F,mat_F*0.01)
      matjit = ifelse(matjit>1, 1, matjit)
      FSS[k,i] = sum(Nay_F[,i]*matjit) #Fecund Stock size
      Npups[k,i] = sum((Nay_F[,i]/2.5)*matjit*rnorm(length(fec), fec, fec*0.1)) # 
      M0[k,i] = rnorm(1, ( ( (1- (Npups[k,i]/Npups_eq)^beta ) * (Zmin - Z0) ) + Z0 ) , 0.1)
      # M0[i] = ( (1- (Npups[i]/Npups_eq)^beta ) * (Zmin - Z0) ) + Z0
      # M0[i] = 1.9044160 # THIS IS EQUILIBRIUM!!!
      
      
      Nay_F[1,i+1] = 0.5 * Npups[k,i]*exp(-M0[k,i])  # pups produced that survive to age 1; Recruits
      Nay_M[1,i+1] = 0.5 * Npups[k,i]*exp(-M0[k,i])  # pups produced that survive to age 1; Recruits
      
      
      Zay_F = vector(length=length(M))
      Zay_M = vector(length=length(M))
      for (j in 1:A){
        Zay_F[j] = M[j] + (Sel[j,2]*FM[1,i]) + (Sel[j,4]*FM[2,i]) + (Sel[j,6]*FM[3,i]) + (Sel[j,8]*FM[4,i])
        Zay_M[j] = M[j] + (Sel[j,3]*FM[1,i]) + (Sel[j,5]*FM[2,i]) + (Sel[j,7]*FM[3,i]) + (Sel[j,8]*FM[4,i])
      } # END J LOOP
      
      for (j in 1:(A-2)){
        Mjit = rnorm(length(M), M, M*0.1) 
        Nay_F[j+1,i+1] = Nay_F[j,i]*exp(-Zay_F[j]) 
        Nay_M[j+1,i+1] = Nay_M[j,i]*exp(-Zay_M[j]) 
        
        
        Cay1_F[j,i] = Nay_F[j,i]*((Sel[j,2]*FM[1,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
        Cay2_F[j,i] = Nay_F[j,i]*((Sel[j,4]*FM[2,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
        Cay3_F[j,i] = Nay_F[j,i]*((Sel[j,6]*FM[3,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
        Cay4_F[j,i] = Nay_F[j,i]*((Sel[j,8]*FM[4,i])/Zay_F[j])*(1-exp(-Zay_F[j]))
        
        Cay1_M[j,i] = Nay_M[j,i]*((Sel[j,3]*FM[1,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
        Cay2_M[j,i] = Nay_M[j,i]*((Sel[j,5]*FM[2,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
        Cay3_M[j,i] = Nay_M[j,i]*((Sel[j,7]*FM[3,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
        Cay4_M[j,i] = Nay_M[j,i]*((Sel[j,8]*FM[4,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
        
        
      } # END J LOOP
      
      Nay_F[A,i+1] = Nay_F[A-1,i]*exp(-Zay_F[A-1]) #+ Nay_F[A,i]*exp(-Zay[A])
      Nay_M[A,i+1] = Nay_M[A-1,i]*exp(-Zay_M[A-1]) #+ Nay_M[A,i]*exp(-Zay[A])
      
      Cay1_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,2]*FM[1,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      Cay2_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,4]*FM[2,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      Cay3_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,6]*FM[3,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      Cay4_F[A-1,i]=Nay_F[A-1,i]*((Sel[A-1,8]*FM[4,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      
      Cay1_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,3]*FM[1,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
      Cay2_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,5]*FM[2,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
      Cay3_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,7]*FM[3,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
      Cay4_M[A-1,i]=Nay_M[A-1,i]*((Sel[A-1,8]*FM[4,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
      
      
      
      Cay1_F[A,i]=Nay_F[A,i]*((Sel[A,2]*FM[1,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      Cay2_F[A,i]=Nay_F[A,i]*((Sel[A,4]*FM[2,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      Cay3_F[A,i]=Nay_F[A,i]*((Sel[A,6]*FM[3,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      Cay4_F[A,i]=Nay_F[A,i]*((Sel[A,8]*FM[4,i])/Zay_F[j])*(1-exp(-Zay_F[j])) 
      
      Cay1_M[A,i]=Nay_M[A,i]*((Sel[A,3]*FM[1,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
      Cay2_M[A,i]=Nay_M[A,i]*((Sel[A,5]*FM[2,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
      Cay3_M[A,i]=Nay_M[A,i]*((Sel[A,7]*FM[3,i])/Zay_M[j])*(1-exp(-Zay_M[j])) 
      Cay4_M[A,i]=Nay_M[A,i]*((Sel[A,8]*FM[4,i])/Zay_M[j])*(1-exp(-Zay_M[j]))
      
      
      for(h in 1:A){
        # h=h+1
        Iay1_CV[h,i] = runif(1, min=CV1-0.1, max=CV1+0.1)
        Iay1_F[h,i] = max(q_survey1[i]*Sel[h,17] * Nay_F[h,i] * 
                            exp(rnorm(1, 0, (q_survey1[i]*Sel[h,17]*Nay_F[h,i])*Iay1_CV[h,i]) - 
                                  ( ((q_survey1[i]*Sel[h,17]*Nay_F[h,i]*Iay1_CV[h,i])^2)/2) ), 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay1_M[h,i] = max(q_survey1[i]*Sel[h,18] * Nay_M[h,i] * 
                            exp(rnorm(1, 0, (q_survey1[i]*Sel[h,18]*Nay_M[h,i])*Iay1_CV[h,i]) - 
                                  ( ((q_survey1[i]*Sel[h,18]*Nay_M[h,i]*Iay1_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        
        Iay2_CV[h,i] = runif(1, min=CV2-0.1, max=CV2+0.1)
        Iay2_F[h,i] = max(q_survey2[i]*Sel[h,9] * Nay_F[h,i] *
                            exp(rnorm(1, 0, (q_survey2[i]*Sel[h,9]*Nay_F[h,i])*Iay2_CV[h,i]) -
                                  ( ((q_survey2[i]*Sel[h,9]*Nay_F[h,i]*Iay2_CV[h,i])^2)/2) ), 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay2_M[h,i] = max(q_survey2[i]*Sel[h,10] * Nay_M[h,i] *
                            exp(rnorm(1, 0, (q_survey2[i]*Sel[h,10]*Nay_M[h,i])*Iay2_CV[h,i]) -
                                  ( ((q_survey2[i]*Sel[h,10]*Nay_M[h,i]*Iay2_CV[h,i])^2)/2) ), 0)
        
        Iay3_CV[h,i] = runif(1, min=CV3-0.1, max=CV3+0.1)
        Iay3_F[h,i] = max(q_survey3[i]*Sel[h,28] * Nay_F[h,i] *
                            exp(rnorm(1, 0, (q_survey3[i]*Sel[h,28]*Nay_F[h,i])*Iay3_CV[h,i]) -
                                  ( ((q_survey3[i]*Sel[h,28]*Nay_F[h,i]*Iay3_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay3_M[h,i] = max(q_survey3[i]*Sel[h,29] * Nay_M[h,i] *
                            exp(rnorm(1, 0, (q_survey3[i]*Sel[h,29]*Nay_M[h,i])*Iay3_CV[h,i]) -
                                  ( ((q_survey3[i]*Sel[h,29]*Nay_M[h,i]*Iay3_CV[h,i]) ^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        
        Iay4_CV[h,i] = runif(1, min=CV4-0.1, max=CV4+0.1)
        Iay4_F[h,i] = max(q_survey4[i]*Sel[h,19] * Nay_F[h,i] *
                            exp(rnorm(1, 0, (q_survey4[i]*Sel[h,19]*Nay_F[h,i])*Iay4_CV[h,i]) -
                                  ( ((q_survey4[i]*Sel[h,19]*Nay_F[h,i]*Iay4_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay4_M[h,i] = max(q_survey4[i]*Sel[h,21] * Nay_M[h,i] * 
                            exp(rnorm(1, 0, (q_survey4[i]*Sel[h,21]*Nay_M[h,i])*Iay4_CV[h,i]) -
                                  ( ((q_survey4[i]*Sel[h,21]*Nay_M[h,i]*Iay4_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        
        Iay5_CV[h,i] = runif(1, min=CV5-0.1, max=CV5+0.1)
        Iay5_F[h,i] = max(q_survey5[i]*Sel[h,23] * Nay_F[h,i] *
                            exp(rnorm(1, 0, (q_survey5[i]*Sel[h,23]*Nay_F[h,i])*Iay5_CV[h,i]) -
                                  ( ((q_survey5[i]*Sel[h,23]*Nay_F[h,i]*Iay5_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay5_M[h,i] = max(q_survey5[i]*Sel[h,24] * Nay_M[h,i] *
                            exp(rnorm(1, 0, (q_survey5[i]*Sel[h,24]*Nay_M[h,i])*Iay5_CV[h,i]) -
                                  ( ((q_survey5[i]*Sel[h,24]*Nay_M[h,i]*Iay5_CV[h,i])^2 )/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        
        Iay6_CV[h,i] = runif(1, min=CV6-0.1, max=CV6+0.1)
        Iay6_F[h,i] = max(q_survey6[i]*Sel[h,13] * Nay_F[h,i] *
                            exp(rnorm(1, 0, (q_survey6[i]*Sel[h,13]*Nay_F[h,i])*Iay6_CV[h,i]) -
                                  ( ((q_survey6[i]*Sel[h,13]*Nay_F[h,i]*Iay6_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay6_M[h,i] = max(q_survey6[i]*Sel[h,12] * Nay_M[h,i] *
                            exp(rnorm(1, 0, (q_survey6[i]*Sel[h,12]*Nay_M[h,i])*Iay6_CV[h,i]) -
                                  ( ((q_survey6[i]*Sel[h,12]*Nay_M[h,i]*Iay6_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        
        Iay7_CV[h,i] = runif(1, min=CV7-0.1, max=CV7+0.1)
        Iay7_F[h,i] = max(q_survey7[i]*Sel[h,27] * Nay_F[h,i] *
                            exp(rnorm(1, 0, (q_survey7[i]*Sel[h,27]*Nay_F[h,i])*Iay7_CV[h,i]) -
                                  ( ((q_survey7[i]*Sel[h,27]*Nay_F[h,i]*Iay7_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
        Iay7_M[h,i] = max(q_survey7[i]*Sel[h,26] * Nay_M[h,i] *
                            exp(rnorm(1, 0, (q_survey7[i]*Sel[h,26]*Nay_M[h,i])*Iay7_CV[h,i]) -
                                  ( ((q_survey7[i]*Sel[h,26]*Nay_M[h,i]*Iay7_CV[h,i])^2)/2) ) , 0) #LOTS OF UNCERTAINTY!based off CV;
      } # end h
      
      
    }# ends i loop
    
    
    
    # USE VON BERT CURVE TO ESTIMATE LENGTH COMP. Let von Bert parameters vary by Niter...???? 
    #   year (implement in i loop)... 
    #   Or should they vary by year class? Could create a matrix of von bert parameters based on age class 
    #       and pull out appropriate values as needed
    # Lt = Linf * (1 - exp(-K*(t - t0)))
    # For sexes combined.
    
    Linf_M=175.5
    K_M=0.143
    t0_M=-2.388
    Linf_F=183.3
    K_F=0.124
    t0_F=-3.098
    # L0 = Linf * (1 - exp(-K*(0-t0)))
    
    
    #####  GET LENGTH FREQUENCIES #####
    
    LFreqList = list()
    
    ## Female Fishery LFreqs ##
    LFreqList[[paste("LF_Cay1F",k,sep="_")]] = My_LFreqs(data=Cay1_F, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90))
    LFreqList[[paste("LF_Cay2F",k,sep="_")]] = My_LFreqs(data=Cay2_F, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90))
    LFreqList[[paste("LF_Cay3F",k,sep="_")]] = My_LFreqs(data=Cay3_F, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90))
    LFreqList[[paste("LF_Cay4F",k,sep="_")]] = My_LFreqs(data=Cay4_F, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90))
    
    ## Male Fishery LFreqs ##
    LFreqList[[paste("LF_Cay1M",k,sep="_")]] = My_LFreqs(data=Cay1_M, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Cay2M",k,sep="_")]] = My_LFreqs(data=Cay2_M, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Cay3M",k,sep="_")]] = My_LFreqs(data=Cay3_M, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Cay4M",k,sep="_")]] = My_LFreqs(data=Cay4_M, multC = 10, Lbins = seq(30, 260, by=2), Lyrs=c(51:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    
    ## Female Survey LFreqs ##
    LFreqList[[paste("LF_Iay1F",k,sep="_")]] = My_LFreqs(data=Iay1_F, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:90))
    LFreqList[[paste("LF_Iay2F",k,sep="_")]] = My_LFreqs(data=Iay2_F, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:90))
    LFreqList[[paste("LF_Iay3F",k,sep="_")]] = My_LFreqs(data=Iay3_F, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:90))
    LFreqList[[paste("LF_Iay4F",k,sep="_")]] = My_LFreqs(data=Iay4_F, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:90))
    LFreqList[[paste("LF_Iay5F",k,sep="_")]] = My_LFreqs(data=Iay5_F, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:90))
    LFreqList[[paste("LF_Iay6F",k,sep="_")]] = My_LFreqs(data=Iay6_F, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:90))
    LFreqList[[paste("LF_Iay7F",k,sep="_")]] = My_LFreqs(data=Iay7_F, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:90))
    
    ## Male Survey LFreqs ##
    LFreqList[[paste("LF_Iay1M",k,sep="_")]] = My_LFreqs(data=Iay1_M, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Iay2M",k,sep="_")]] = My_LFreqs(data=Iay2_M, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Iay3M",k,sep="_")]] = My_LFreqs(data=Iay3_M, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Iay4M",k,sep="_")]] = My_LFreqs(data=Iay4_M, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Iay5M",k,sep="_")]] = My_LFreqs(data=Iay5_M, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Iay6M",k,sep="_")]] = My_LFreqs(data=Iay6_M, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    LFreqList[[paste("LF_Iay7M",k,sep="_")]] = My_LFreqs(data=Iay7_M, multC = 1, Lbins = seq(30, 260, by=2), Lyrs=c(51:90),
                                                         Linf=Linf_M, K=K_M, t0=t0_F)
    
    
    save(LFreqList, file=paste("LFreqList_",k,".RData",sep=""))
    
    #### END LFREQS #####
    
    
    
    Iy1[k,] = apply(Iay1_F, 2, sum) + apply(Iay1_M, 2, sum)
    Iy2[k,] = apply(Iay2_F, 2, sum) + apply(Iay2_M, 2, sum)
    Iy3[k,] = apply(Iay3_F, 2, sum) + apply(Iay3_M, 2, sum)
    Iy4[k,] = apply(Iay4_F, 2, sum) + apply(Iay4_M, 2, sum)
    Iy5[k,] = apply(Iay5_F, 2, sum) + apply(Iay5_M, 2, sum)
    Iy6[k,] = apply(Iay6_F, 2, sum) + apply(Iay6_M, 2, sum)
    Iy7[k,] = apply(Iay7_F, 2, sum) + apply(Iay7_M, 2, sum)
    
    lines(1:years, Iy1[k,], type='l', col='grey')
    lines(1:years, Iy2[k,], type='l', col='blue')
    lines(1:years, Iy3[k,], type='l', col='red')
    lines(1:years, Iy4[k,], type='l', col='green')
    lines(1:years, Iy5[k,], type='l', col='orange')
    lines(1:years, Iy6[k,], type='l', col='mediumorchid2')
    lines(1:years, Iy7[k,], type='l', col='cyan')
    Iy1_CV[k,] = apply(Iay1_CV, 2, mean, na.rm=T)
    Iy2_CV[k,] = apply(Iay2_CV, 2, mean, na.rm=T)
    Iy3_CV[k,] = apply(Iay3_CV, 2, mean, na.rm=T)
    Iy4_CV[k,] = apply(Iay4_CV, 2, mean, na.rm=T)
    Iy5_CV[k,] = apply(Iay5_CV, 2, mean, na.rm=T)
    Iy6_CV[k,] = apply(Iay6_CV, 2, mean, na.rm=T)
    Iy7_CV[k,] = apply(Iay7_CV, 2, mean, na.rm=T)
    
    Cy_1[k,] = apply(Cay1_F, 2, sum, na.rm=T) + apply(Cay1_M, 2, sum, na.rm=T) 
    Cy_2[k,] = apply(Cay2_F, 2, sum, na.rm=T) + apply(Cay2_M, 2, sum, na.rm=T) 
    Cy_3[k,] = apply(Cay3_F, 2, sum, na.rm=T) + apply(Cay3_M, 2, sum, na.rm=T) 
    Cy_4[k,] = apply(Cay4_F, 2, sum, na.rm=T) + apply(Cay4_M, 2, sum, na.rm=T) 
    Cy[k,] = Cy_1[k,] + Cy_2[k,] + Cy_3[k,] + Cy_4[k,]
    
    Ny_F[k,] = apply(Nay_F, 2, sum, na.rm=T)
    Ny_M[k,] = apply(Nay_M, 2, sum, na.rm=T)
    Ny[k,] = Ny_F[k,] + Ny_M[k,]
    # lines(1:years, Ny[k,]/100, type='l', lwd=2)
    
    
    #write 
    AnnualNC = list()
    # write.csv(Nay_F, paste(paste("NayF", k, sep="_"),".csv",sep=""))
    # write.csv(Nay_M, paste(paste("NayM", k, sep="_"),".csv",sep=""))
    assign(paste("NayF", k, sep="_"), Nay_F)
    assign(paste("NayM", k, sep="_"), Nay_M)
    AnnualNC[[paste("NayF",k,sep="_")]] = Nay_F
    AnnualNC[[paste("NayM",k,sep="_")]] = Nay_M
    
    # write.csv(Cay_F, paste(paste("CayF", k, sep="_"),".csv",sep=""))
    # write.csv(Cay_M, paste(paste("CayM", k, sep="_"),".csv",sep=""))
    assign(paste("Cay1F", k, sep="_"), Cay1_F)
    assign(paste("Cay2F", k, sep="_"), Cay2_F)
    assign(paste("Cay3F", k, sep="_"), Cay3_F)
    assign(paste("Cay4F", k, sep="_"), Cay4_F)
    assign(paste("Cay1M", k, sep="_"), Cay1_M)
    assign(paste("Cay2M", k, sep="_"), Cay2_M)
    assign(paste("Cay3M", k, sep="_"), Cay3_M)
    assign(paste("Cay4M", k, sep="_"), Cay4_M)
    AnnualNC[[paste("Cay1F",k,sep="_")]] = Cay1_F
    AnnualNC[[paste("Cay2F",k,sep="_")]] = Cay2_F
    AnnualNC[[paste("Cay3F",k,sep="_")]] = Cay3_F
    AnnualNC[[paste("Cay4F",k,sep="_")]] = Cay4_F
    AnnualNC[[paste("Cay1M",k,sep="_")]] = Cay1_M
    AnnualNC[[paste("Cay2M",k,sep="_")]] = Cay2_M
    AnnualNC[[paste("Cay3M",k,sep="_")]] = Cay3_M
    AnnualNC[[paste("Cay4M",k,sep="_")]] = Cay4_M
    
    save(AnnualNC, file=paste("AnnualNC_",k,".RData",sep=""))
    
    # write.csv() ########### WRITE CSVs to SAVE Nay AND Cay!!!
    
  } # end k loop
  
  
  # 
  
  ResultsList[["Iy1"]] = Iy1
  ResultsList[["Iy2"]] = Iy2
  ResultsList[["Iy3"]] = Iy3
  ResultsList[["Iy4"]] = Iy4
  ResultsList[["Iy5"]] = Iy5
  ResultsList[["Iy6"]] = Iy6
  ResultsList[["Iy7"]] = Iy7
  ResultsList[["Iy1_CV"]] = Iy1_CV
  ResultsList[["Iy2_CV"]] = Iy2_CV
  ResultsList[["Iy3_CV"]] = Iy3_CV
  ResultsList[["Iy4_CV"]] = Iy4_CV
  ResultsList[["Iy5_CV"]] = Iy5_CV
  ResultsList[["Iy6_CV"]] = Iy6_CV
  ResultsList[["Iy7_CV"]] = Iy7_CV
  
  ResultsList[["M0"]] = M0
  ResultsList[["FSS"]] = FSS
  ResultsList[["Npups"]] = Npups
  
  ResultsList[["Cy"]] = Cy
  ResultsList[["Cy_1"]] = Cy_1
  ResultsList[["Cy_2"]] = Cy_2
  ResultsList[["Cy_3"]] = Cy_3
  ResultsList[["Cy_4"]] = Cy_4
  
  ResultsList[["Ny"]] = Ny
  
  save(ResultsList, file="ResultsList.RData")
  
  # 
  # write.csv(Iy1, "Iy1.csv")
  # write.csv(Iy2, "Iy2.csv")
  # write.csv(Iy3, "Iy3.csv")  
  # write.csv(Iy4, "Iy4.csv")
  # write.csv(Iy5, "Iy5.csv")
  # write.csv(Iy6, "Iy6.csv")  
  # write.csv(Iy7, "Iy7.csv")  
  # 
  # write.csv(M0, "M0.csv") # populates where each row is an iteration and each column is one year
  # write.csv(FSS, "FSS.csv") # populates where each row is an iteration and each column is one year
  # write.csv(Npups, "Npups.csv") # populates where each row is an iteration and each column is one year
  # 
  # write.csv(Iy1_CV, "Iy1_CV.csv")
  # write.csv(Iy2_CV, "Iy2_CV.csv")
  # write.csv(Iy3_CV, "Iy3_CV.csv")
  # write.csv(Iy4_CV, "Iy4_CV.csv")
  # write.csv(Iy5_CV, "Iy5_CV.csv")
  # write.csv(Iy6_CV, "Iy6_CV.csv")
  # write.csv(Iy7_CV, "Iy7_CV.csv")
  # 
  # write.csv(Cy, "Cy.csv")
  # write.csv(Ny, "Ny.csv")
  
  
  
  ########## TEST DFA #############
  # 51 - 90; 
  # library(MARSS)
  yrs = 51:90
  biasAll = matrix(nrow=Niters, ncol=length(yrs))
  bias = vector(length=Niters)
  RMSE = vector(length=Niters)
  biasAll.BT = matrix(nrow=Niters, ncol=length(yrs))
  RMSE.BT = vector(length=Niters)
  
  FitRatios = matrix(nrow=Niters, ncol=7)
  FactorLoadings=vector()
  DFATrends = vector()
  DFATrendsBT = vector()
  DFATrendsSE = vector()
  DFATrendsSEBT = vector()
  upCI_DFATrends = vector()
  lowCI_DFATrends = vector()
  
  datz_SD = matrix(nrow=Niters, ncol=7)
  
  # i=1
  
  
  for(i in 1:Niters){
    I1 = Iy1[i,51:90]
    I2 = Iy2[i,51:90]
    I3 = Iy3[i,51:90]
    I4 = Iy4[i,51:90]
    I5 = Iy5[i,51:90]
    I6 = Iy6[i,51:90]
    I7 = Iy7[i,51:90]
    # I1 = c(Iy1[i,51:53], rep(NA, 5), Iy1[i,59:64], NA, NA, Iy1[i,67:90])
    # I2 = c(rep(NA, 7), Iy2[i,58:90])
    # I3 = c(rep(NA, 13), Iy3[i,64:90])
    # I4 = c(rep(NA, 16), Iy4[i,67:90])
    # I5 = c(rep(NA, 22), Iy5[i,73:90])
    # I6 = c(rep(NA, 15), Iy6[i,66:79], rep(NA, 11))
    # I7 = c(rep(NA, 17), Iy7[i,68], NA, NA, Iy7[i,71], NA, NA, Iy7[i,74], NA, NA, Iy7[i,77], NA, NA, Iy7[i,80], 
    #        NA, NA, Iy7[i,83], NA, NA, Iy7[i,86], NA, NA, Iy7[i,89], NA)
    # 
    assign(paste("dat",i,sep=""), rbind(I1, I2, I3, I4, I5, I6, I7))
    # assign(paste("dat",i,sep=""), rbind(Iy1[i,51:90], Iy2[i,51:90], Iy3[i,51:90], Iy4[i,51:90], Iy5[i,51:90], Iy6[i,51:90], Iy7[i,51:90]))
    
    dat.a = get(paste("dat",i,sep=""))

    
    
    # c=c(3.5, 1.3, 4.5, 6, 10, 5.5, 2.8)
    # c=c(3.0, 7, 4.5, 6, 5000, 5.5, 2.8) # SB127; i=23
    # c=c(3.5, 1.3, 4.5, 6, 5000, 5.5, 2.8) # SB127; i=24
        # c=c(3.5, 0.8, 4, 0.4, 7, 0.6, 2.5)
    # c=c(14, 1.7, 1.3, 50, 55, 0.5, 0.9) # SB140
    # # c=c(200, 4, 400, 0.3, 500, 2000, 200) #SB141
    # c=c(3.5, 0.7, 4.5, 200, 6, 300, 3) #SB138
    
    dat=dat.a * c
    
    TT=ncol(dat)
    N.ts = nrow(dat)
    # Standardize data
    datL=log(dat)
    y.bar = apply(datL, 1, mean, na.rm=TRUE)
    dat.dm = (datL-y.bar) / y.bar
    gsd = sd(c(dat.dm), na.rm=T)
    dat.z = dat.dm/gsd
    y.bar
    gsd
    datz_SD[i,] = apply(dat.z, 1, sd, na.rm=T) ; datz_SD[i,]
    
    
    
    dat.z = as.matrix(dat.z)
    rownames(dat.z) = c("Survey1","Survey2","Survey3","Survey4","Survey5","Survey6","Survey7")
    
    
    
    
    ##### SB DELTA LOGNORMAL ####
    cntl.list = list(minit=200, maxit=500000, abstol=0.02, allow.degen=FALSE, conv.test.slope.tol = 0.05)
    R = diag(c(CV1, CV2, CV3, CV4, CV5, CV6, CV7), nrow=7, ncol=7)
    # dfa = MARSS(dat.z, model = list(m=1,R="diagonal and equal"), control=cntl.list, form="dfa", z.score=TRUE)
    dfa = MARSS(dat.z, model = list(m=1,R=R), control=cntl.list, form="dfa", z.score=FALSE)
    
   
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
    
    
    # plot factor loadings
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
    
    par(mfrow=c(4,2))
    df <- data.frame(yrs, index, lowerCI, upperCI)
    with(df, plot(index ~ yrs, type='n', bty='L',xlab='Year', ylab='Abundance', ylim=c(min(lowerCI), max(upperCI))))
    x <- as.numeric(df$yrs)
    polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
    with(df, lines(x, index, type='l', lwd=2, pch=16))
    title(paste("Iteration",i,sep=" "))
    abline(h=0)
    
    # Nay = get(paste("NayF", i, sep="_")) + get(paste("NayM", i, sep="_"))
    # N = apply(Nay[,yrs],2,sum)
    N = Ny[i,51:90]
    Nscale = (N-mean(N))/sd(N)
    NL = log(N)
    NLscale = (NL - mean(NL))/sd(NL)
    index.z = (index - mean(index))/sd(index)
    index.z = ifelse(is.na(index.z), 0, index.z)
    index.z = ifelse(abs(index.z)==Inf, 0, index.z)
    
    # indexBT = exp(index*gsd)
    # assign(paste("indexBT",i,sep=""), indexBT)
    # 
    # indexSEBT = indexSE * gsd
    # assign(paste("indexSEBT",i,sep=""), indexSEBT)
    
    indexSEBT = indexSE * gsd
    assign(paste("indexSEBT",i,sep=""), indexSEBT)
    
    indexBT = exp(index*gsd + ((indexSEBT^2)/2) )
    assign(paste("indexBT",i,sep=""), indexBT)
    
    
    lines(yrs,rescale(NL,to=c(min(index),max(index))), type='l', col='red', lwd=2)
    assign(paste("N",i,sep=""), N)
    assign(paste("Nscale",i,sep=""),Nscale)
    assign(paste("NL",i,sep=""), NL)
    assign(paste("NLscale",i,sep=""),NLscale)
    # write.csv(N,paste("N",i,".csv",sep=""))
    
    
    biasAll[i,] = index.z - NLscale
    bias[i] = mean(index.z - NLscale)
    RMSE[i] = sqrt(sum((index.z-NLscale)^2)/length(index.z))
    # RMSE_yr[i] = sqrt(sum((Nscale - indexscale)^2)/length(indexscale))
    
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
      title(paste("Sandbar",survey[n], sep=" "))
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
  
  colnames(FactorLoadings) = c("FL1","FL2","FL3","FL4","FL5","FL6","FL7",
                               "lowCI_FL1","lowCI_FL2","lowCI_FL3","lowCI_FL4","lowCI_FL5","lowCI_FL6","lowCI_FL7",
                               "upCI_FL1","upCI_FL2","upCI_FL3","upCI_FL4","upCI_FL5","upCI_FL6","upCI_FL7")
  FR = apply(FitRatios,1,mean)
  FLoadings = cbind(FactorLoadings, meanFactorLoading = FR)
  
  DFA_Results = list()
  DFA_Results[["FitRatios"]] = FitRatios
  DFA_Results[["DFATrends"]] = DFATrends
  DFA_Results[["DFATrendsSE"]] = DFATrendsSE
  DFA_Results[["DFATrendsBT"]] = DFATrendsBT
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
  # write.csv(DFATrends, "DFATrends.csv")
  # write.csv(DFATrendsSE, "DFATrendsSE.csv")
  # write.csv(upCI_DFATrends, "upCI_DFATrends.csv")
  # write.csv(lowCI_DFATrends, "lowCI_DFATrends.csv")
  
  # write.csv(FLoadings, "FactorLoadings_sum.csv")
  
  # write.csv(bias, "bias.csv")
  # write.csv(biasAll, "biasAll.csv")
  # write.csv(RMSE, "RMSE.csv")
  
  
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
    text(x = min(yrs)+4, y = 0.8*max(min(lowerCI)), labels = format( FR[l], digits=3 ) , cex=1)
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
    plot(51:90, Iy1[l,51:90], type='l', col='grey', xlim=c(51, 90), axes=F,
         ylim=c(0, max(c(Iy1[l,], Iy2[l,], Iy3[l,],Iy4[l,], Iy5[l,], Iy6[l,], Iy7[l,]), na.rm=T)+0.1) )
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    lines(51:90, Iy2[l,51:90], type='l', col='blue')
    lines(51:90, Iy3[l,51:90], type='l', col='red')
    lines(51:90, Iy4[l,51:90], type='l', col='green')
    lines(51:90, Iy5[l,51:90], type='l', col='orange')
    lines(51:90, Iy6[l,51:90], type='l', col='mediumorchid2')
    lines(51:90, Iy7[l,51:90], type='l', col='cyan')
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
    text(x = min(yrs)+4, y = 0.8*min(lowerCI), labels = format( FR[l], digits=3 ) , cex=1)
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
    plot(yrs, N, type='l',col="red", axes=FALSE)
    axis(1, labels=FALSE)
    axis(2, labels=FALSE)
    lines(yrs, rescale(index, to=c(min(N),max(N))))
    text(x = min(yrs)+4, y = 1.05*min(N), labels = format( FR[l], digits=3 ) , cex=1)
    
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
    text(x = min(yrs)+4, y = 0.8*min((index-mean(index))/sd(index)), labels = format( FR[l], digits=3 ) , cex=1)
    
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
  #   text(x = min(yrs)+4, y = 0.8*min(lowerCI), labels = format( FR[l], digits=3 ) , cex=1)
  #   par(new=T)
  #   NL = get(paste("NL",l,sep=""))
  #   plot(yrs,NL, type='l', col='red', axes=FALSE, lwd=1.5)
  #   axis(4, labels=FALSE)
  # }
  # 
  # dev.off()
  
  
  
  
  
  
  
  
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
    text(x = min(yrs)+4, y = 1.2*min(NL), labels = format( FR[l], digits=3 ) , cex=1)
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
  
  return(list(RMSE=RMSE, biasAll=biasAll, FitRatios = FitRatios, RMSE.BT=RMSE.BT, biasAll.BT=biasAll.BT))} # END FUNCTION





########################
# Survey 1 = VIMS
# Survey 2 = LPS
# Survey 3 = PLLOP
# Survey 4 = NMFS LL SE
# Survey 5 = CST NE LL
# Survey 6 = BLLOP 1
# Survey 7 = NMFS NE
#---
# q0 = constant
# q1 = period 2
# q2 = period 3
# q3 = step increase (step=5)
# q4 = step decrease (step=5)
# q5 = increase then decrease
# q6 = increase
# q7 = decrease
# q8 = random between 0.025 and 0.05
# q9 = jitter increase
# q10 = jitter decrease
# q11 = jitter period 2
# q12 = jitter period 3
# q13 = jitter constant
# q14 = jitter constant2
# q15 = jitter constant3
# q16 = jitter constant4
# q17 = jitter constant5
# q18 = jitter constant6
# q19 = jitter constant7
# q20 = jitter constant8
# q21 = jitter constant9

# save.image("N:/Documents/DFA_Simulation/SB/SB_Sim.RData")


#### MISSING DATA #############

set.seed(111)
q9a = jitter(q6, amount=0.005)                                                           # q9 = jitter increase
set.seed(222)
q9b = jitter(q6, amount=0.005)                                                           # q9 = jitter increase

# SB40
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial40")
# c = c(45, 4.0, 5, 4, 2000, 0.3, 1.2) #40
SB40 = DFA_Sim(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q9, # jitter increase
                     q_survey2 = q10, # jitter decrease
                     q_survey3 = q6, # increase
                     q_survey4 = q7, # decrease
                     q_survey5 = q3, # step increase
                     q_survey6 = q9a, # jitter increase
                     q_survey7 = q9b, # jitter increase
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
               c=c(45, 4.0, 5, 4, 2000, 0.3, 1.2)  )
save(SB40, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial40\\SB40.RData")


# SB41
set.seed(111)
q10a = jitter(q7, amount=0.005)                                                           # q9 = jitter increase
set.seed(222)
q10b = jitter(q7, amount=0.005)                                                           # q9 = jitter increase

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial41")
# c=c(29, 2.9, 21, 0.4, 2100, 1.6, 2) #41
SB41 = DFA_Sim(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q10, # jitter decrease
                     q_survey2 = q9, # jitter increase
                     q_survey3 = q7, # decrease
                     q_survey4 = q6, # increase
                     q_survey5 = q4, # step decrease
                     q_survey6 = q10a, # jitter decrease
                     q_survey7 = q10b, # jitter decrease
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
               c=c(29, 2.9, 21, 0.4, 2100, 1.6, 2) )
save(SB41, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial41\\SB41.RData")


# NO CHANGE IN Q ---------------------------------------------------------------------------------------------------------
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial1")
SB1 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 3, 1, 0.4, 11, 0.8, 0.35))
save(SB1, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial1\\SB1.RData")

summary(SB1$RMSE)
summary(apply(SB1$FitRatios,1,FUN=mean))


load(file="DFA_Results.RData")

# CHANGE increase 1 Q ---------------------------------------------------------------------------------------------------------
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial2")
SB2 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q10,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(500, 3, 1, 0.4, 11, 0.8, 0.35))
save(SB2, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial2\\SB2.RData")
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial3")
SB3 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q10,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 100, 1, 0.4, 11, 0.8, 0.35))
save(SB3, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial3\\SB3.RData")


setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial4")
SB4 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q10,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 3, 200, 0.4, 11, 0.8, 0.35))
save(SB4, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial4\\SB4.RData")


setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial5")
SB5 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q10,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 3, 1, 100, 11, 0.8, 0.35))
save(SB5, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial5\\SB5.RData")


setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial6")
SB6 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q10,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 3, 1, 0.4, 150, 0.8, 0.35))
save(SB6, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial6\\SB6.RData")


setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial7")
SB7 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q10,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 3, 1, 0.4, 11, 4, 0.35))
save(SB7, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial7\\SB7.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial8")
SB8 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q10,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 3, 1, 0.4, 11, 0.8, 400))

save(SB8, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial8\\SB8.RData")


# PLOT -------------------------------------------------------------------------------
par(mfrow=c(3,1))
vioplot(c(SB1$biasAll), c(SB2$biasAll), c(SB3$biasAll), c(SB4$biasAll), c(SB5$biasAll), c(SB6$biasAll), c(SB7$biasAll), c(SB8$biasAll), 
        col='deepskyblue')
abline(h=mean(c(SB1$biasAll)))
abline(h=summary(c(SB1$biasAll))[2], lty=2)
abline(h=summary(c(SB1$biasAll))[5], lty=2)
vioplot(SB1$RMSE, SB2$RMSE, SB3$RMSE, SB4$RMSE, SB5$RMSE, SB6$RMSE, SB7$RMSE, SB8$RMSE, col='deepskyblue')
abline(h=mean(SB1$RMSE))
abline(h=summary(SB1$RMSE)[2], lty=2)
abline(h=summary(SB1$RMSE)[5], lty=2)
vioplot(apply(SB1$FitRatios, 1, mean), apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), 
        apply(SB5$FitRatios, 1, mean), apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean), 
        col='deepskyblue')
abline(h=mean(apply(SB1$FitRatios, 1, mean)))
abline(h=summary(apply(SB1$FitRatios, 1, mean))[2], lty=2)
abline(h=summary(apply(SB1$FitRatios, 1, mean))[5], lty=2)





# CHANGE decrease 1 Q ---------------------------------------------------------------------------------------------------------
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial9")
SB9 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q9,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(7000, 3, 1, 0.4, 11, 0.8, 0.35))
save(SB9, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial9\\SB9.RData")
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial10")
SB10 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q9,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 12, 1, 0.4, 11, 0.8, 0.35))
save(SB10, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial10\\SB10.RData")
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial11")
SB11 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q9,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 3, 100, 0.4, 11, 0.8, 0.35))
save(SB11, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial11\\SB11.RData")
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial12")
SB12 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q9,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 3, 1, 0.9, 11, 0.8, 0.35))
save(SB12, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial12\\SB12.RData")
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial13")
SB13 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q9,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 3, 1, 0.4, 11000000, 0.8, 0.35))
save(SB13, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial13\\SB13.RData")
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial14")
SB14 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q9,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 3, 1, 0.4, 11, 0.8, 0.35))

save(SB14, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial14\\SB14.RData")
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial15")
SB15 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q9,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 3, 1, 0.4, 11, 0.8, 50))

save(SB15, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial15\\SB15.RData")



# Plot ---------------------------------------------------------------------------------------------------------
par(mfrow=c(3,2), mar=c(0.9, 1.4, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB1$biasAll), c(SB2$biasAll), c(SB3$biasAll), c(SB4$biasAll), c(SB5$biasAll), c(SB6$biasAll), c(SB7$biasAll), c(SB8$biasAll), 
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue', 'deepskyblue1','deepskyblue2','deepskyblue3','deepskyblue4'),
        names=FALSE, ylim=c(-1.25,1.25))
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
abline(h=0, col='black')
title('Decreasing q', line=-0.9)
vioplot(c(SB1$biasAll), c(SB9$biasAll), c(SB10$biasAll), c(SB11$biasAll), c(SB12$biasAll), c(SB13$biasAll), c(SB14$biasAll), c(SB15$biasAll), 
        col=c('grey','olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen'),
        names=FALSE, ylim=c(-5, 2))
abline(h=0, col='black')
title('Increasing q', line=-0.9)

vioplot(SB1$RMSE, SB2$RMSE, SB3$RMSE, SB4$RMSE, SB5$RMSE, SB6$RMSE, SB7$RMSE, SB8$RMSE,  
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue', 'deepskyblue1','deepskyblue2','deepskyblue3','deepskyblue4'),
        names=FALSE, ylim=c(0,0.6))
mtext('RMSE',side=2,line=1.1,cex=0.8)
abline(h=mean(SB1$RMSE), col='white')
vioplot(SB1$RMSE, SB9$RMSE, SB10$RMSE, SB11$RMSE, SB12$RMSE, SB13$RMSE, SB14$RMSE, SB15$RMSE, 
        col=c('grey','olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen'),
        names=FALSE, ylim=c(0,1.7))
abline(h=mean(SB1$RMSE), col='white')

vioplot(apply(SB1$FitRatios, 1, mean), apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), 
        apply(SB5$FitRatios, 1, mean), apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean), 
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue', 'deepskyblue1','deepskyblue2','deepskyblue3','deepskyblue4'),
        names=c('None','S1','S2','S3','S4','S5','S6','S7'), ylim=c(0,1))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
abline(h=mean(apply(SB1$FitRatios, 1, mean)), col='white')
vioplot(apply(SB1$FitRatios, 1, mean), apply(SB9$FitRatios, 1, mean), apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), 
        apply(SB12$FitRatios, 1, mean), apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean), 
        col=c('grey','olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen') ,
        names=c('None','S1','S2','S3','S4','S5','S6','S7'), ylim=c(0,1))
abline(h=mean(apply(SB1$FitRatios, 1, mean)), col='white')







png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Changing_q.png",
     type="cairo",
     units="mm",
     width=300,
     height=250,
     pointsize=16,
     res=600)
# STANDARDIZE AXES
par(mfrow=c(3,2), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB1$biasAll), c(SB2$biasAll), c(SB3$biasAll), c(SB4$biasAll), c(SB5$biasAll), c(SB6$biasAll), c(SB7$biasAll), c(SB8$biasAll), 
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue', 'deepskyblue1','deepskyblue2','deepskyblue3','deepskyblue4'),
        names=FALSE, ylim=c(-5,2))
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
abline(h=0, col='black')
title('Decreasing q', line=-0.9)
vioplot(c(SB1$biasAll), c(SB9$biasAll), c(SB10$biasAll), c(SB11$biasAll), c(SB12$biasAll), c(SB13$biasAll), c(SB14$biasAll), c(SB15$biasAll), 
        col=c('grey','olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen'),
        names=FALSE, ylim=c(-5, 2))
abline(h=0, col='black')
title('Increasing q', line=-0.9)

vioplot(SB1$RMSE, SB2$RMSE, SB3$RMSE, SB4$RMSE, SB5$RMSE, SB6$RMSE, SB7$RMSE, SB8$RMSE,  
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue', 'deepskyblue1','deepskyblue2','deepskyblue3','deepskyblue4'),
        names=FALSE, ylim=c(0,1.7))
mtext('RMSE',side=2,line=1.1,cex=0.8)
abline(h=mean(SB1$RMSE), col='white')
vioplot(SB1$RMSE, SB9$RMSE, SB10$RMSE, SB11$RMSE, SB12$RMSE, SB13$RMSE, SB14$RMSE, SB15$RMSE, 
        col=c('grey','olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen'),
        names=FALSE, ylim=c(0,1.7))
abline(h=mean(SB1$RMSE), col='white')

vioplot(apply(SB1$FitRatios, 1, mean), apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), 
        apply(SB5$FitRatios, 1, mean), apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean), 
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue', 'deepskyblue1','deepskyblue2','deepskyblue3','deepskyblue4'),
        names=c('None','S1','S2','S3','S4','S5','S6','S7'), ylim=c(0,0.8))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
abline(h=mean(apply(SB1$FitRatios, 1, mean)), col='white')
vioplot(apply(SB1$FitRatios, 1, mean), apply(SB9$FitRatios, 1, mean), apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), 
        apply(SB12$FitRatios, 1, mean), apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean), 
        col=c('grey','olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen') ,
        names=c('None','S1','S2','S3','S4','S5','S6','S7'), ylim=c(0,0.8))
abline(h=mean(apply(SB1$FitRatios, 1, mean)), col='white')

dev.off()

##### NO MISSING DATA ########
# q3 = step increase (step=5)
# q4 = step decrease (step=5)
# q5 = increase then decrease
# q6 = increase
# q7 = decrease
# q8 = random between 0.025 and 0.05
# q9 = jitter increase
# q10 = jitter decrease

# BASE All CHANGE IN Q ---------------------------------------------------------------------------------------------------------

### EXTRA RUNS ####
### SB140
set.seed(111)
q9a = jitter(q6, amount=0.005)                                                           # q9 = jitter increase
set.seed(222)
q9b = jitter(q6, amount=0.005)                                                           # q9 = jitter increase

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial140")
# c=c(14, 1.7, 1.3, 50, 55, 0.5, 0.9) # SB140
SB140 = DFA_Sim_Full(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q9, # jitter increase
                     q_survey2 = q10, # jitter decrease
                     q_survey3 = q6, # increase
                     q_survey4 = q7, # decrease
                     q_survey5 = q3, # step increase
                     q_survey6 = q9a, # jitter increase
                     q_survey7 = q9b, # jitter increase
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
                     c=c(14, 1.7, 1.3, 50, 55, 0.5, 0.9))
save(SB140, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial140\\SB140.RData")


# SB141
set.seed(111)
q10a = jitter(q7, amount=0.005)                                                           # q9 = jitter increase
set.seed(222)
q10b = jitter(q7, amount=0.005)                                                           # q9 = jitter increase

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial141")
# c=c(200, 4, 400, 0.3, 500, 2000, 200) #SB141
SB141 = DFA_Sim_Full(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q10, # jitter decrease
                     q_survey2 = q9, # jitter increase
                     q_survey3 = q7, # decrease
                     q_survey4 = q6, # increase
                     q_survey5 = q4, # step decrease
                     q_survey6 = q10a, # jitter decrease
                     q_survey7 = q10b, # jitter decrease
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
                     c=c(200, 4, 400, 0.3, 500, 2000, 200) )
save(SB141, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial141\\SB141.RData")

# BASE NO CHANGE IN Q ---------------------------------------------------------------------------------------------------------
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial101")
SB101 = DFA_Sim_Full(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c=c(3.5, 1.3, 4.5, 6, 6, 5.5, 3.5))
save(SB101, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial101\\SB101.RData")

#####
summary(SB1$RMSE)
summary(apply(SB1$FitRatios,1,FUN=mean))

summary(SB101$RMSE)
summary(apply(SB101$FitRatios,1,FUN=mean))

png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Missing_Data.png", 
    type="cairo",
    units="mm", 
    width=100, 
    height=200, 
    pointsize=16, 
    res=600)
par(mfrow=c(3,1), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll), c(SB1$biasAll), 
        col=c('cadetblue1','red'),
        names=FALSE, ylim=c(-2, 1.5))
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
abline(h=0, col='black')
title('Effect of Missing Data', line=-0.9)

vioplot(SB101$RMSE, SB1$RMSE,   
        col=c('cadetblue1','red'),
        names=FALSE, ylim=c(0,0.8))
mtext('RMSE',side=2,line=1.1,cex=0.8)

vioplot(apply(SB101$FitRatios, 1, mean), apply(SB1$FitRatios, 1, mean), 
        col=c('cadetblue1','red'),
        names=c('Complete','Missing Data'), ylim=c(0,3))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
# abline(h=mean(apply(SB1$FitRatios, 1, mean)), col='white')
dev.off()

# CHANGE increase 1 Q ---------------------------------------------------------------------------------------------------------
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial102")
SB102 = DFA_Sim_Full(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q10,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c=c(30, 1.3, 4.5, 6, 6, 5.5, 3.5))
save(SB102, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial102\\SB102.RData")
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial103")
SB103 = DFA_Sim_Full(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q10,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c=c(3.5, 3, 4.5, 6, 6, 5.5, 3.5))
save(SB103, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial103\\SB103.RData")
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial104")
SB104 = DFA_Sim_Full(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q10,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c=c(3.5, 1.3, 50, 6, 6, 5.5, 3.5))
save(SB104, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial104\\SB104.RData")
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial105")
SB105 = DFA_Sim_Full(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q10,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c=c(3.5, 1.3, 4.5, 200, 6, 5.5, 3.5))
save(SB105, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial105\\SB105.RData")
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial106")
SB106 = DFA_Sim_Full(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q10,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c=c(3.5, 1.3, 4.5, 6, 15, 5.5, 3.5))
save(SB106, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial106\\SB106.RData")
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial107")
SB107 = DFA_Sim_Full(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q10,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c=c(3.5, 1.3, 4.5, 6, 6, 300, 3.5))
save(SB107, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial107\\SB107.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial108")
SB108 = DFA_Sim_Full(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q10,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c=c(3.5, 1.3, 4.5, 6, 6, 5.5, 60))
save(SB108, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial108\\SB108.RData")


par(mfrow=c(3,1))
vioplot(c(SB101$biasAll), c(SB102$biasAll), c(SB103$biasAll), c(SB104$biasAll), c(SB105$biasAll), c(SB106$biasAll), c(SB107$biasAll), c(SB108$biasAll), 
        col='deepskyblue')
abline(h=mean(c(SB101$biasAll)))
vioplot(SB101$RMSE, SB102$RMSE, SB103$RMSE, SB104$RMSE, SB105$RMSE, SB106$RMSE, SB107$RMSE, SB108$RMSE, col='deepskyblue')
abline(h=mean(SB101$RMSE))
vioplot(apply(SB101$FitRatios, 1, mean), apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), 
        apply(SB105$FitRatios, 1, mean), apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean), 
        col='deepskyblue')
abline(h=mean(apply(SB101$FitRatios, 1, mean)))





# CHANGE decrease 1 Q ---------------------------------------------------------------------------------------------------------
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial109")
SB109 = DFA_Sim_Full(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q9,
              q_survey2 = q0,
              q_survey3 = q0,
              q_survey4 = q0,
              q_survey5 = q0,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c=c(40, 1.3, 4.5, 6, 6, 5.5, 3.5))
save(SB109, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial109\\SB109.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial110")
SB110 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q9,
               q_survey3 = q0,
               q_survey4 = q0,
               q_survey5 = q0,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(3.5, 1.3, 4.5, 6, 6, 5.5, 3.5))
save(SB110, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial110\\SB110.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial111")
SB111 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q9,
               q_survey4 = q0,
               q_survey5 = q0,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(3.5, 1.3, 3.5, 6, 6, 5.5, 3.5))
save(SB111, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial111\\SB111.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial112")
SB112 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q0,
               q_survey4 = q9,
               q_survey5 = q0,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(3.5, 1.3, 4.5, 0.4, 6, 5.5, 3.5))
save(SB112, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial112\\SB112.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial113")
SB113 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q0,
               q_survey4 = q0,
               q_survey5 = q9,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(3.5, 1.3, 4.5, 6, 350, 5.5, 3.5))
save(SB113, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial113\\SB113.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial114")
SB114 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q0,
               q_survey4 = q0,
               q_survey5 = q0,
               q_survey6 = q9,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(3.5, 1.3, 4.5, 6, 6, 0.6, 3.5))

save(SB114, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial114\\SB114.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial115")
SB115 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q0,
               q_survey4 = q0,
               q_survey5 = q0,
               q_survey6 = q0,
               q_survey7 = q9,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(3.5, 1.3, 4.5, 6, 6, 5.5, 1.5))

save(SB115, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial115\\SB115.RData")




# PLOT  ---------------------------------------------------------------------------------------------------------
par(mfrow=c(3,2), mar=c(0.9, 1.4, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll), c(SB102$biasAll), c(SB103$biasAll), c(SB104$biasAll), c(SB105$biasAll), c(SB106$biasAll), c(SB107$biasAll), c(SB108$biasAll), 
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue', 'deepskyblue1','deepskyblue2','deepskyblue3','deepskyblue4'),
        names=FALSE, ylim=c(-0.6,0.6))
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
abline(h=0, col='black')
title('Decreasing q', line=-0.9)
vioplot(c(SB101$biasAll), c(SB109$biasAll), c(SB110$biasAll), c(SB111$biasAll), c(SB112$biasAll), c(SB113$biasAll), c(SB114$biasAll), c(SB115$biasAll), 
        col=c('grey','olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen'),
        names=FALSE, ylim=c(-0.6,0.6))
abline(h=0, col='black')
title('Increasing q', line=-0.9)

vioplot(SB101$RMSE, SB102$RMSE, SB103$RMSE, SB104$RMSE, SB105$RMSE, SB106$RMSE, SB107$RMSE, SB108$RMSE,  
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue', 'deepskyblue1','deepskyblue2','deepskyblue3','deepskyblue4'),
        names=FALSE, ylim=c(0,0.25))
mtext('RMSE',side=2,line=1.1,cex=0.8)
abline(h=mean(SB101$RMSE), col='white')
vioplot(SB101$RMSE, SB109$RMSE, SB110$RMSE, SB111$RMSE, SB112$RMSE, SB113$RMSE, SB114$RMSE, SB115$RMSE, 
        col=c('grey','olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen'),
        names=FALSE, ylim=c(0,0.25))
abline(h=mean(SB101$RMSE), col='white')

vioplot(apply(SB101$FitRatios, 1, mean), apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), 
        apply(SB105$FitRatios, 1, mean), apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean), 
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue', 'deepskyblue1','deepskyblue2','deepskyblue3','deepskyblue4'),
        names=c('None','S1','S2','S3','S4','S5','S6','S7'), ylim=c(0,0.4))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
abline(h=mean(apply(SB101$FitRatios, 1, mean)), col='white')
vioplot(apply(SB101$FitRatios, 1, mean), apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), 
        apply(SB112$FitRatios, 1, mean), apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean), 
        col=c('grey','olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen') ,
        names=c('None','S1','S2','S3','S4','S5','S6','S7'), ylim=c(0,0.4))
abline(h=mean(apply(SB101$FitRatios, 1, mean)), col='white')







png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Changing_q_CompleteData.png", 
    type="cairo",
    units="mm", 
    width=300, 
    height=250, 
    pointsize=16, 
    res=600)
# STANDARDIZE AXES
par(mfrow=c(3,2), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll), c(SB102$biasAll), c(SB103$biasAll), c(SB104$biasAll), c(SB105$biasAll), c(SB106$biasAll), c(SB107$biasAll), c(SB108$biasAll), 
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue', 'deepskyblue1','deepskyblue2','deepskyblue3','deepskyblue4'),
        names=FALSE, ylim=c(-0.6,0.6))
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
abline(h=0, col='black')
title('Decreasing q', line=-0.9)
vioplot(c(SB101$biasAll), c(SB109$biasAll), c(SB110$biasAll), c(SB111$biasAll), c(SB112$biasAll), c(SB113$biasAll), c(SB114$biasAll), c(SB115$biasAll), 
        col=c('grey','olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen'),
        names=FALSE, ylim=c(-0.6,0.6))
abline(h=0, col='black')
title('Increasing q', line=-0.9)

vioplot(SB101$RMSE, SB102$RMSE, SB103$RMSE, SB104$RMSE, SB105$RMSE, SB106$RMSE, SB107$RMSE, SB108$RMSE,  
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue', 'deepskyblue1','deepskyblue2','deepskyblue3','deepskyblue4'),
        names=FALSE, ylim=c(0,0.25))
mtext('RMSE',side=2,line=1.1,cex=0.8)
abline(h=mean(SB101$RMSE), col='white')
vioplot(SB101$RMSE, SB109$RMSE, SB110$RMSE, SB111$RMSE, SB112$RMSE, SB113$RMSE, SB114$RMSE, SB115$RMSE, 
        col=c('grey','olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen'),
        names=FALSE, ylim=c(0,0.25))
abline(h=mean(SB101$RMSE), col='white')

vioplot(apply(SB101$FitRatios, 1, mean), apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), 
        apply(SB105$FitRatios, 1, mean), apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean), 
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue', 'deepskyblue1','deepskyblue2','deepskyblue3','deepskyblue4'),
        names=c('None','S1','S2','S3','S4','S5','S6','S7'), ylim=c(0,0.4))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
abline(h=mean(apply(SB101$FitRatios, 1, mean)), col='white')
vioplot(apply(SB101$FitRatios, 1, mean), apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), 
        apply(SB112$FitRatios, 1, mean), apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean), 
        col=c('grey','olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen') ,
        names=c('None','S1','S2','S3','S4','S5','S6','S7'), ylim=c(0,0.4))
abline(h=mean(apply(SB101$FitRatios, 1, mean)), col='white')

dev.off()





png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\FullVSMissingData.png", 
    type="cairo",
    units="mm", 
    width=200, 
    height=250, 
    pointsize=16, 
    res=600)
# STANDARDIZE AXES
par(mfrow=c(3,2), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll), c(SB102$biasAll, SB103$biasAll, SB104$biasAll, SB105$biasAll, SB106$biasAll, SB107$biasAll, SB108$biasAll),
        c(SB1$biasAll), c(SB2$biasAll, SB3$biasAll, SB4$biasAll, SB5$biasAll, SB6$biasAll, SB7$biasAll, SB8$biasAll),
        col=c('cadetblue1','cadetblue3','red','red3'),
        names=FALSE, ylim=c(-5,2))
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
abline(h=0, col='black')
title('Decreasing q', line=-0.9)
vioplot(c(SB101$biasAll), c(SB109$biasAll, SB110$biasAll, SB111$biasAll, SB112$biasAll, SB113$biasAll, SB114$biasAll, SB115$biasAll),
        c(SB1$biasAll), c(SB9$biasAll, SB10$biasAll, SB11$biasAll, SB12$biasAll, SB13$biasAll, SB14$biasAll, SB15$biasAll),
        col=c('cadetblue1','cadetblue3','red','red3'),
        names=FALSE, ylim=c(-5, 2))
abline(h=0, col='black')
title('Increasing q', line=-0.9)

vioplot(SB101$RMSE, c(SB102$RMSE, SB103$RMSE, SB104$RMSE, SB105$RMSE, SB106$RMSE, SB107$RMSE, SB108$RMSE),
        SB1$RMSE, c(SB2$RMSE, SB3$RMSE, SB4$RMSE, SB5$RMSE, SB6$RMSE, SB7$RMSE, SB8$RMSE),
        col=c('cadetblue1','cadetblue3','red','red3'),
        names=FALSE, ylim=c(0,1.7))
mtext('RMSE',side=2,line=1.1,cex=0.8)
vioplot(SB101$RMSE, c(SB109$RMSE, SB110$RMSE, SB111$RMSE, SB112$RMSE, SB113$RMSE, SB114$RMSE, SB115$RMSE), 
        SB1$RMSE, c(SB9$RMSE, SB10$RMSE, SB11$RMSE, SB12$RMSE, SB13$RMSE, SB14$RMSE, SB15$RMSE), 
        col=c('cadetblue1','cadetblue3','red','red3'),
        names=FALSE, ylim=c(0,1.7))


vioplot(apply(SB101$FitRatios, 1, mean),
        c(apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), apply(SB105$FitRatios, 1, mean), 
          apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean)), 
        apply(SB1$FitRatios, 1, mean),
        c(apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), apply(SB5$FitRatios, 1, mean), 
          apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean)), 
        col=c('cadetblue1','cadetblue3','red','red3'),
        names=c(
'Full
const q',
'Full
chang q',
'Missing
const q',
'Missing
chang q'), ylim=c(0,0.8))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)

vioplot(apply(SB101$FitRatios, 1, mean), 
        c(apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), apply(SB112$FitRatios, 1, mean), 
          apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean)), 
        apply(SB1$FitRatios, 1, mean), 
        c(apply(SB9$FitRatios, 1, mean), apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), apply(SB12$FitRatios, 1, mean), 
          apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean)), 
        col=c('cadetblue1','cadetblue3','red','red3'),
        names=c(
'Full
const q',
'Full
chang q',
'Missing
const q',
'Missing
chang q'), ylim=c(0,0.8))


dev.off()










######## ROUND 2 ##################

# Change 2 qs missing data ---------------------------------------------------------------------------------------------------------------------
# both qs increase
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial16")
SB16 = DFA_Sim(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q9,
               q_survey2 = q0,
               q_survey3 = q3,
               q_survey4 = q0,
               q_survey5 = q0,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c = c(7000, 3, 100, 0.4, 11, 0.8, 0.35))

save(SB16, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial16\\SB16.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial17")
SB17 = DFA_Sim(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q9,
               q_survey4 = q0,
               q_survey5 = q3,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c = c(12, 3, 100, 0.4, 1100000, 0.8, 0.35))
save(SB17, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial17\\SB17.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial18")
SB18 = DFA_Sim(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q0,
               q_survey4 = q0,
               q_survey5 = q9,
               q_survey6 = q0,
               q_survey7 = q3,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c = c(12, 3, 1, 0.4, 11000000, 0.8, 5))
save(SB18, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial18\\SB18.RData")


# both qs decrease
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial19")
SB19 = DFA_Sim(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q10,
               q_survey2 = q0,
               q_survey3 = q4,
               q_survey4 = q0,
               q_survey5 = q0,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c = c(500, 3, 25, 0.4, 11, 0.8, 0.35))
save(SB19, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial19\\SB19.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial20")
SB20 = DFA_Sim(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q10,
               q_survey4 = q0,
               q_survey5 = q4,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c = c(12, 3, 200, 0.4, 50000, 0.8, 0.35))
save(SB20, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial20\\SB20.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial21")
SB21 = DFA_Sim(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q0,
               q_survey4 = q0,
               q_survey5 = q10,
               q_survey6 = q0,
               q_survey7 = q4,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c = c(12, 3, 1, 0.4, 200, 0.8, 20))
save(SB21, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial21\\SB21.RData")

# q incr - decr
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial22")
SB22 = DFA_Sim(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q9,
               q_survey2 = q0,
               q_survey3 = q4,
               q_survey4 = q0,
               q_survey5 = q0,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c = c(7000, 3, 22, 0.4, 11, 0.8, 0.35))
save(SB22, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial22\\SB22.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial23")
SB23 = DFA_Sim(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q9,
               q_survey4 = q0,
               q_survey5 = q4,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c = c(12, 3, 150, 0.4, 50000, 0.8, 0.35))
save(SB23, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial23\\SB23.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial24")
SB24 = DFA_Sim(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q0,
               q_survey4 = q0,
               q_survey5 = q9,
               q_survey6 = q0,
               q_survey7 = q4,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c = c(12, 3, 1, 0.4, 11000000, 0.8, 20))
save(SB24, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial24\\SB24.RData")

# q decr - incr
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial25")
SB25 = DFA_Sim(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q10,
               q_survey2 = q0,
               q_survey3 = q3,
               q_survey4 = q0,
               q_survey5 = q0,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c = c(400, 3, 150, 0.4, 11, 0.8, 0.35))
save(SB25, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial25\\SB25.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial26")
SB26 = DFA_Sim(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q10,
               q_survey4 = q0,
               q_survey5 = q3,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c = c(12, 3, 150, 0.4, 1100000, 0.8, 0.35))
save(SB26, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial26\\SB26.RData")

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



# Change 2 qs full data ---------------------------------------------------------------------------------------------------------------------
# both qs increase
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial116")
SB116 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q9,
               q_survey2 = q0,
               q_survey3 = q3,
               q_survey4 = q0,
               q_survey5 = q0,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(40, 1.3, 3, 6, 6, 5.5, 3.5))
save(SB116, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial116\\SB116.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial117")
SB117 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q9,
               q_survey4 = q0,
               q_survey5 = q3,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(3.5, 1.3, 3, 6, 350, 5.5, 3.5))
save(SB117, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial117\\SB117.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial118")
SB118 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q0,
               q_survey4 = q0,
               q_survey5 = q9,
               q_survey6 = q0,
               q_survey7 = q3,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(3.5, 1.3, 4.5, 6, 350, 5.5, 0.7) )
save(SB118, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial118\\SB118.RData")


# both qs decrease
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial119")
SB119 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q10,
               q_survey2 = q0,
               q_survey3 = q4,
               q_survey4 = q0,
               q_survey5 = q0,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(35, 1.3, 4.5, 6, 6, 5.5, 3.5))
save(SB119, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial119\\SB119.RData")


setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial120")
SB120 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q10,
               q_survey4 = q0,
               q_survey5 = q4,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(3.5, 1.3, 45, 6, 150, 5.5, 3.5))
save(SB120, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial120\\SB120.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial121")
SB121 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q0,
               q_survey4 = q0,
               q_survey5 = q10,
               q_survey6 = q0,
               q_survey7 = q4,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(3.5, 1.3, 4.5, 6, 13, 5.5, 3))
save(SB121, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial121\\SB121.RData")

# q incr - decr
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial122")
SB122 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q9,
               q_survey2 = q0,
               q_survey3 = q4,
               q_survey4 = q0,
               q_survey5 = q0,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(40, 1.3, 4, 6, 6, 5.5, 3.5))
save(SB122, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial122\\SB122.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial123")
SB123 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q9,
               q_survey4 = q0,
               q_survey5 = q4,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(3.5, 1.3, 3.5, 6, 150, 5.5, 3.5))
save(SB123, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial123\\SB123.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial124")
SB124 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q0,
               q_survey4 = q0,
               q_survey5 = q9,
               q_survey6 = q0,
               q_survey7 = q4,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(3.5, 1.3, 4.5, 6, 350, 5.5, 2.5))
save(SB124, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial124\\SB124.RData")

# q decr - incr
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial125")
SB125 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q10,
               q_survey2 = q0,
               q_survey3 = q3,
               q_survey4 = q0,
               q_survey5 = q0,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(35, 1.3, 4.5, 6, 6, 5.5, 3.5))
save(SB125, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial125\\SB125.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial126")
SB126 = DFA_Sim_Full(Tot_F = F1,
               CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
               q_survey1 = q0,
               q_survey2 = q0,
               q_survey3 = q10,
               q_survey4 = q0,
               q_survey5 = q3,
               q_survey6 = q0,
               q_survey7 = q0,
               Iay_eq = Iay_eq_F1,
               Niters=100,
               c=c(3.5, 1.3, 50, 6, 100, 5.5, 3.5))
save(SB126, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial126\\SB126.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial127")
SB127 = DFA_Sim_Full(Tot_F = F1,
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
               c=c(3.5, 1.3, 4.5, 6, 10, 5.5, 2.8))

save(SB127, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial127\\SB127.RData")




# PLOT  ---------------------------------------------------------------------------------------------------------

png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\CompleteData_ChangeQ2.png", 
    type="cairo",
    units="mm", 
    width=600, 
    height=250, 
    pointsize=16, 
    res=600)
# STANDARDIZE AXES
par(mfrow=c(3,1), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll), c(SB102$biasAll), c(SB103$biasAll), c(SB104$biasAll), c(SB105$biasAll), c(SB106$biasAll), c(SB107$biasAll), c(SB108$biasAll),
        c(SB109$biasAll), c(SB110$biasAll), c(SB111$biasAll), c(SB112$biasAll), c(SB113$biasAll), c(SB114$biasAll), c(SB115$biasAll),
        c(SB119$biasAll), c(SB120$biasAll), c(SB121$biasAll), c(SB116$biasAll), c(SB117$biasAll), c(SB118$biasAll), 
        c(SB122$biasAll), c(SB123$biasAll), c(SB124$biasAll), c(SB125$biasAll), c(SB126$biasAll), c(SB127$biasAll),
        
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'darkslategray1','darkslategray3','darkslategray4','darkseagreen1','darkseagreen2','darkseagreen3',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4'),
        names=FALSE, ylim=c(-0.7,0.7) )
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
abline(h=0, col='black')
title('Complete Data', line=-0.9)

vioplot(SB101$RMSE, SB102$RMSE, SB103$RMSE, SB104$RMSE, SB105$RMSE, SB106$RMSE, SB107$RMSE, SB108$RMSE,
        SB109$RMSE, SB110$RMSE, SB111$RMSE, SB112$RMSE, SB113$RMSE, SB114$RMSE, SB115$RMSE,
        SB119$RMSE, SB120$RMSE, SB121$RMSE, SB116$RMSE, SB117$RMSE, SB118$RMSE, 
        SB122$RMSE, SB123$RMSE, SB124$RMSE, SB125$RMSE, SB126$RMSE, SB127$RMSE,
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'darkslategray1','darkslategray3','darkslategray4','darkseagreen1','darkseagreen2','darkseagreen3',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4'),
        names=FALSE, ylim=c(0,0.30))
mtext('RMSE',side=2,line=1.1,cex=0.8)
# abline(h=mean(SB101$RMSE), col='white')


vioplot(apply(SB101$FitRatios, 1, mean), apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), 
        apply(SB105$FitRatios, 1, mean), apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean), 
        apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), 
        apply(SB112$FitRatios, 1, mean), apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean),
        
        apply(SB119$FitRatios, 1, mean), apply(SB120$FitRatios, 1, mean), apply(SB121$FitRatios, 1, mean), 
        apply(SB116$FitRatios, 1, mean), apply(SB117$FitRatios, 1, mean), apply(SB118$FitRatios, 1, mean), 
        apply(SB122$FitRatios, 1, mean), apply(SB123$FitRatios, 1, mean), apply(SB124$FitRatios, 1, mean), 
        apply(SB125$FitRatios, 1, mean), apply(SB126$FitRatios, 1, mean), apply(SB127$FitRatios, 1, mean), 
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'darkslategray1','darkslategray3','darkslategray4','darkseagreen1','darkseagreen2','darkseagreen3',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4'),
        names=c('None','D-S1','D-S2','D-S3','D-S4','D-S5','D-S6','D-S7','I-S1','I-S2','I-S3','I-S4','I-S5','I-S6','I-S7',
                'D2','D2','D2','I2','I2','I2','I1-D1','I1-D1','I1-D1','D1-I1','D1-I1','D1-I1'))#, ylim=c(0,0.45)
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
# abline(h=mean(apply(SB101$FitRatios, 1, mean)), col='white')


dev.off()




png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\MissingData_ChangeQ2.png", 
    type="cairo",
    units="mm", 
    width=600, 
    height=250, 
    pointsize=16, 
    res=600)
# STANDARDIZE AXES
par(mfrow=c(3,1), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB1$biasAll), c(SB2$biasAll), c(SB3$biasAll), c(SB4$biasAll), c(SB5$biasAll), c(SB6$biasAll), c(SB7$biasAll), c(SB8$biasAll),
        c(SB9$biasAll), c(SB10$biasAll), c(SB11$biasAll), c(SB12$biasAll), c(SB13$biasAll), c(SB14$biasAll), c(SB15$biasAll),
        c(SB19$biasAll), c(SB20$biasAll), c(SB21$biasAll), c(SB16$biasAll), c(SB17$biasAll), c(SB18$biasAll), 
        c(SB22$biasAll), c(SB23$biasAll), c(SB24$biasAll), c(SB25$biasAll), c(SB26$biasAll), c(SB27$biasAll),
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'darkslategray1','darkslategray3','darkslategray4','darkseagreen1','darkseagreen2','darkseagreen3',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4'),
        names=FALSE, ylim=c(-5.25,2) )
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
abline(h=0, col='black')
title('Missing Data', line=-0.9)

vioplot(SB1$RMSE, SB2$RMSE, SB3$RMSE, SB4$RMSE, SB5$RMSE, SB6$RMSE, SB7$RMSE, SB8$RMSE,
        SB9$RMSE,  SB10$RMSE, SB11$RMSE, SB12$RMSE, SB13$RMSE, SB14$RMSE, SB15$RMSE,
        SB19$RMSE, SB20$RMSE, SB21$RMSE, SB16$RMSE, SB17$RMSE, SB18$RMSE, 
        SB22$RMSE, SB23$RMSE, SB24$RMSE, SB25$RMSE, SB26$RMSE, SB27$RMSE,
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'darkslategray1','darkslategray3','darkslategray4','darkseagreen1','darkseagreen2','darkseagreen3',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4'),
        names=FALSE, ylim=c(0,1.8))
mtext('RMSE',side=2,line=1.1,cex=0.8)
# abline(h=mean(SB101$RMSE), col='white')


vioplot(apply(SB1$FitRatios, 1, mean), apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), 
        apply(SB5$FitRatios, 1, mean), apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean), 
        apply(SB9$FitRatios, 1, mean),  apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), 
        apply(SB12$FitRatios, 1, mean), apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean),
        
        apply(SB19$FitRatios, 1, mean), apply(SB20$FitRatios, 1, mean), apply(SB21$FitRatios, 1, mean), 
        apply(SB16$FitRatios, 1, mean), apply(SB17$FitRatios, 1, mean), apply(SB18$FitRatios, 1, mean), 
        apply(SB22$FitRatios, 1, mean), apply(SB23$FitRatios, 1, mean), apply(SB24$FitRatios, 1, mean), 
        apply(SB25$FitRatios, 1, mean), apply(SB26$FitRatios, 1, mean), apply(SB27$FitRatios, 1, mean), 
        col=c('grey','skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'darkslategray1','darkslategray3','darkslategray4','darkseagreen1','darkseagreen2','darkseagreen3',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4'),
        names=c('None','D-S1','D-S2','D-S3','D-S4','D-S5','D-S6','D-S7','I-S1','I-S2','I-S3','I-S4','I-S5','I-S6','I-S7',
                'D2','D2','D2','I2','I2','I2','I1-D1','I1-D1','I1-D1','D1-I1','D1-I1','D1-I1'))#, ylim=c(0,0.8)
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
# abline(h=mean(apply(SB101$FitRatios, 1, mean)), col='white')


dev.off()







png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\FullVSMissingData_2.png", 
    type="cairo",
    units="mm", 
    width=200, 
    height=250, 
    pointsize=16, 
    res=600)
# STANDARDIZE AXES
par(mfrow=c(3,2), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll), 
        c(SB102$biasAll, SB103$biasAll, SB104$biasAll, SB105$biasAll, SB106$biasAll, SB107$biasAll, SB108$biasAll),
        c(SB109$biasAll, SB110$biasAll, SB111$biasAll, SB112$biasAll, SB113$biasAll, SB114$biasAll, SB115$biasAll),
        c(SB119$biasAll, SB120$biasAll, SB121$biasAll), c(SB116$biasAll, SB117$biasAll, SB118$biasAll),
        c(SB122$biasAll, SB123$biasAll, SB124$biasAll), c(SB125$biasAll, SB126$biasAll, SB127$biasAll),
        col=c('grey','deepskyblue', 'darkolivegreen1', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3'),
        names=FALSE, ylim=c(-5,2))
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
abline(h=0, col='black')
title('Complete Data', line=-0.9)
vioplot(c(SB1$biasAll), 
        c(SB2$biasAll, SB3$biasAll, SB4$biasAll, SB5$biasAll, SB6$biasAll, SB7$biasAll, SB8$biasAll),
        c(SB9$biasAll,  SB10$biasAll, SB11$biasAll, SB12$biasAll, SB13$biasAll, SB14$biasAll, SB15$biasAll),
        c(SB19$biasAll, SB20$biasAll, SB21$biasAll), c(SB16$biasAll, SB17$biasAll, SB18$biasAll),
        c(SB22$biasAll, SB23$biasAll, SB24$biasAll), c(SB25$biasAll, SB26$biasAll, SB27$biasAll),
        col=c('grey','deepskyblue', 'darkolivegreen1', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3'),
        names=FALSE, ylim=c(-5, 2))
abline(h=0, col='black')
title('Missing Data', line=-0.9)

vioplot(SB101$RMSE, 
        c(SB102$RMSE, SB103$RMSE, SB104$RMSE, SB105$RMSE, SB106$RMSE, SB107$RMSE, SB108$RMSE),
        c(SB109$RMSE, SB110$RMSE, SB111$RMSE, SB112$RMSE, SB113$RMSE, SB114$RMSE, SB115$RMSE),
        c(SB119$RMSE, SB120$RMSE, SB121$RMSE), c(SB116$RMSE, SB117$RMSE, SB118$RMSE),
        c(SB122$RMSE, SB123$RMSE, SB124$RMSE), c(SB125$RMSE, SB126$RMSE, SB127$RMSE),
        col=c('grey','deepskyblue', 'darkolivegreen1', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3'),
        names=FALSE, ylim=c(0,1.7))
mtext('RMSE',side=2,line=1.1,cex=0.8)
vioplot(SB1$RMSE, c(SB2$RMSE, SB3$RMSE, SB4$RMSE, SB5$RMSE, SB6$RMSE, SB7$RMSE, SB8$RMSE), 
        c(SB9$RMSE, SB10$RMSE, SB11$RMSE, SB12$RMSE, SB13$RMSE, SB14$RMSE, SB15$RMSE),
        c(SB19$RMSE, SB20$RMSE, SB21$RMSE), c(SB16$RMSE, SB17$RMSE, SB18$RMSE),
        c(SB22$RMSE, SB23$RMSE, SB24$RMSE), c(SB25$RMSE, SB26$RMSE, SB27$RMSE),
        col=c('grey','deepskyblue', 'darkolivegreen1', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3'),
        names=FALSE, ylim=c(0,1.7))


vioplot(apply(SB101$FitRatios, 1, mean),
        c(apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), apply(SB105$FitRatios, 1, mean), 
          apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean)), 
        c(apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), apply(SB112$FitRatios, 1, mean), 
          apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean)), 
        c(apply(SB119$FitRatios, 1, mean), apply(SB120$FitRatios, 1, mean), apply(SB121$FitRatios, 1, mean)), 
        c(apply(SB116$FitRatios, 1, mean), apply(SB117$FitRatios, 1, mean), apply(SB118$FitRatios, 1, mean)), 
        c(apply(SB122$FitRatios, 1, mean), apply(SB123$FitRatios, 1, mean), apply(SB124$FitRatios, 1, mean)), 
        c(apply(SB125$FitRatios, 1, mean), apply(SB126$FitRatios, 1, mean), apply(SB127$FitRatios, 1, mean)), 
        col=c('grey','deepskyblue', 'darkolivegreen1', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3'),
        names=c(
'No
change','D1','I1','D2','I2','I1-D1','D1-I1'), ylim=c(0,0.8))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)

vioplot(apply(SB1$FitRatios, 1, mean),
        c(apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), apply(SB5$FitRatios, 1, mean), 
          apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean)), 
        c(apply(SB9$FitRatios, 1, mean), apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), apply(SB12$FitRatios, 1, mean), 
          apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean)), 
        c(apply(SB19$FitRatios, 1, mean), apply(SB20$FitRatios, 1, mean), apply(SB21$FitRatios, 1, mean)), 
        c(apply(SB16$FitRatios, 1, mean), apply(SB17$FitRatios, 1, mean), apply(SB18$FitRatios, 1, mean)), 
        c(apply(SB22$FitRatios, 1, mean), apply(SB23$FitRatios, 1, mean), apply(SB24$FitRatios, 1, mean)), 
        c(apply(SB25$FitRatios, 1, mean), apply(SB26$FitRatios, 1, mean), apply(SB27$FitRatios, 1, mean)), 
        col=c('grey','deepskyblue', 'darkolivegreen1', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3'),
        names=c(
'No
change','D1','I1','D2','I2','I1-D1','D1-I1'),  ylim=c(0,0.8))


dev.off()







png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\FullVSMissingData_3.png", 
    type="cairo",
    units="mm", 
    width=200, 
    height=250, 
    pointsize=16, 
    res=600)
# STANDARDIZE AXES
par(mfrow=c(3,1), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll), 
        c(SB1$biasAll), 
        c(SB102$biasAll, SB103$biasAll, SB104$biasAll, SB105$biasAll, SB106$biasAll, SB107$biasAll, SB108$biasAll),
        c(SB2$biasAll, SB3$biasAll, SB4$biasAll, SB5$biasAll, SB6$biasAll, SB7$biasAll, SB8$biasAll),
        c(SB109$biasAll, SB110$biasAll, SB111$biasAll, SB112$biasAll, SB113$biasAll, SB114$biasAll, SB115$biasAll),
        c(SB9$biasAll,  SB10$biasAll, SB11$biasAll, SB12$biasAll, SB13$biasAll, SB14$biasAll, SB15$biasAll),
        c(SB119$biasAll, SB120$biasAll, SB121$biasAll), 
        c(SB19$biasAll, SB20$biasAll, SB21$biasAll), 
        c(SB116$biasAll, SB117$biasAll, SB118$biasAll),
        c(SB16$biasAll, SB17$biasAll, SB18$biasAll),
        c(SB122$biasAll, SB123$biasAll, SB124$biasAll), 
        c(SB22$biasAll, SB23$biasAll, SB24$biasAll), 
        c(SB125$biasAll, SB126$biasAll, SB127$biasAll),
        c(SB25$biasAll, SB26$biasAll, SB27$biasAll),
        col=c('darkgrey','grey','skyblue1','deepskyblue','olivedrab1','olivedrab3','darkslategray1','darkslategray4','darkseagreen1',
              'darkseagreen4','orchid','mediumorchid1','darkorchid1','darkorchid4'),
        names=FALSE, ylim=c(-5.5,2.5))
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
title('Complete vs. Missing Data', line=-0.9)
abline(h=0, col='black')



vioplot(SB101$RMSE, 
        SB1$RMSE,
        c(SB102$RMSE, SB103$RMSE, SB104$RMSE, SB105$RMSE, SB106$RMSE, SB107$RMSE, SB108$RMSE),
        c(SB2$RMSE, SB3$RMSE, SB4$RMSE, SB5$RMSE, SB6$RMSE, SB7$RMSE, SB8$RMSE), 
        c(SB109$RMSE, SB110$RMSE, SB111$RMSE, SB112$RMSE, SB113$RMSE, SB114$RMSE, SB115$RMSE),
        c(SB9$RMSE, SB10$RMSE, SB11$RMSE, SB12$RMSE, SB13$RMSE, SB14$RMSE, SB15$RMSE),
        c(SB119$RMSE, SB120$RMSE, SB121$RMSE), 
        c(SB19$RMSE, SB20$RMSE, SB21$RMSE),
        c(SB116$RMSE, SB117$RMSE, SB118$RMSE),
        c(SB16$RMSE, SB17$RMSE, SB18$RMSE),
        c(SB122$RMSE, SB123$RMSE, SB124$RMSE), 
        c(SB22$RMSE, SB23$RMSE, SB24$RMSE),
        c(SB125$RMSE, SB126$RMSE, SB127$RMSE),
        c(SB25$RMSE, SB26$RMSE, SB27$RMSE),
        col=c('darkgrey','grey','skyblue1','deepskyblue','olivedrab1','olivedrab3','darkslategray1','darkslategray4','darkseagreen1',
              'darkseagreen4','orchid','mediumorchid1','darkorchid1','darkorchid4'),
        names=FALSE, ylim=c(0,1.7))
mtext('RMSE',side=2,line=1.1,cex=0.8)


vioplot(apply(SB101$FitRatios, 1, mean),
        apply(SB1$FitRatios, 1, mean),
        c(apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), apply(SB105$FitRatios, 1, mean), 
          apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean)), 
        c(apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), apply(SB5$FitRatios, 1, mean), 
          apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean)), 
        c(apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), apply(SB112$FitRatios, 1, mean), 
          apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean)), 
        c(apply(SB9$FitRatios, 1, mean), apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), apply(SB12$FitRatios, 1, mean), 
          apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean)), 
        c(apply(SB119$FitRatios, 1, mean), apply(SB120$FitRatios, 1, mean), apply(SB121$FitRatios, 1, mean)), 
        c(apply(SB19$FitRatios, 1, mean), apply(SB20$FitRatios, 1, mean), apply(SB21$FitRatios, 1, mean)), 
        c(apply(SB116$FitRatios, 1, mean), apply(SB117$FitRatios, 1, mean), apply(SB118$FitRatios, 1, mean)), 
        c(apply(SB16$FitRatios, 1, mean), apply(SB17$FitRatios, 1, mean), apply(SB18$FitRatios, 1, mean)), 
        c(apply(SB122$FitRatios, 1, mean), apply(SB123$FitRatios, 1, mean), apply(SB124$FitRatios, 1, mean)), 
        c(apply(SB22$FitRatios, 1, mean), apply(SB23$FitRatios, 1, mean), apply(SB24$FitRatios, 1, mean)), 
        c(apply(SB125$FitRatios, 1, mean), apply(SB126$FitRatios, 1, mean), apply(SB127$FitRatios, 1, mean)), 
        c(apply(SB25$FitRatios, 1, mean), apply(SB26$FitRatios, 1, mean), apply(SB27$FitRatios, 1, mean)), 
        col=c('darkgrey','grey','skyblue1','deepskyblue','olivedrab1','olivedrab3','darkslategray1','darkslategray4','darkseagreen1',
              'darkseagreen4','orchid','mediumorchid1','darkorchid1','darkorchid4'),
        names=c(
'No
change','','D1','','I1','','D2','','I2','','I1-D1','','D1-I1',''), ylim=c(0,1.4))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)


dev.off()








######################################################################################################






################### ROUND 3 ##############################################################
# Survey 1 = VIMS
# Survey 2 = LPS
# Survey 3 = PLLOP
# Survey 4 = NMFS LL SE
# Survey 5 = CST NE LL
# Survey 6 = BLLOP 1
# Survey 7 = NMFS NE
#---
# q0 = constant
# q1 = period 2
# q2 = period 3
# q3 = step increase (step=5)
# q4 = step decrease (step=5)
# q5 = increase then decrease
# q6 = increase
# q7 = decrease
# q8 = random between 0.025 and 0.05
# q9 = jitter increase
# q10 = jitter decrease
# q11 = jitter period 2
# q12 = jitter period 3


# Change 3 qs missing data -------------------------------------------------------------------------------------------
# Increase
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial28")
SB28 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q3,
              q_survey2 = q0,
              q_survey3 = q9,
              q_survey4 = q0,
              q_survey5 = q9,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(7000, 3, 100, 0.4, 11000000, 0.8, 0.35))
save(SB28, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial28\\SB28.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial29")
SB29 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q3,
              q_survey3 = q0,
              q_survey4 = q9,
              q_survey5 = q0,
              q_survey6 = q9,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 4, 1, 0.9, 11, 0.8, 0.35))
save(SB29, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial29\\SB29.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial30")
SB30 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q3,
              q_survey4 = q0,
              q_survey5 = q9,
              q_survey6 = q0,
              q_survey7 = q9,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 3, 100, 0.4, 11000000, 0.8, 80))
save(SB30, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial30\\SB30.RData")

# decrease
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial31")
SB31 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q4,
              q_survey2 = q0,
              q_survey3 = q10,
              q_survey4 = q0,
              q_survey5 = q10,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(700, 3, 200, 0.4, 150, 0.8, 0.35))
save(SB31, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial31\\SB31.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial32")
SB32 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q4,
              q_survey3 = q0,
              q_survey4 = q10,
              q_survey5 = q0,
              q_survey6 = q10,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 25, 1, 100, 11, 4, 0.35))
save(SB32, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial32\\SB32.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial33")
SB33 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q4,
              q_survey4 = q0,
              q_survey5 = q10,
              q_survey6 = q0,
              q_survey7 = q10,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 3, 20, 0.4, 150, 0.8, 400))
save(SB33, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial33\\SB33.RData")

# D - I2
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial34")
SB34 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q4,
              q_survey2 = q0,
              q_survey3 = q9,
              q_survey4 = q0,
              q_survey5 = q9,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(1500, 3, 200, 0.4, 11000000, 0.8, 0.35))
save(SB34, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial34\\SB34.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial35")
SB35 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q4,
              q_survey3 = q0,
              q_survey4 = q9,
              q_survey5 = q0,
              q_survey6 = q9,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 25, 1, 1, 11, 0.8, 0.35))
save(SB35, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial35\\SB35.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial36")
SB36 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q4,
              q_survey4 = q0,
              q_survey5 = q9,
              q_survey6 = q0,
              q_survey7 = q9,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 3, 20, 0.4, 11000000, 0.8, 80))
save(SB36, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial36\\SB36.RData")

# I - D2
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial37")
SB37 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q3,
              q_survey2 = q0,
              q_survey3 = q10,
              q_survey4 = q0,
              q_survey5 = q10,
              q_survey6 = q0,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(10000, 3, 200, 0.4, 200, 0.8, 0.35))
save(SB37, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial37\\SB37.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial38")
SB38 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q3,
              q_survey3 = q0,
              q_survey4 = q10,
              q_survey5 = q0,
              q_survey6 = q10,
              q_survey7 = q0,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 5, 1, 100, 11, 4, 0.35))
save(SB38, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial38\\SB38.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial39")
SB39 = DFA_Sim(Tot_F = F1,
              CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
              q_survey1 = q0,
              q_survey2 = q0,
              q_survey3 = q3,
              q_survey4 = q0,
              q_survey5 = q10,
              q_survey6 = q0,
              q_survey7 = q10,
              Iay_eq = Iay_eq_F1,
              Niters=100,
              c = c(12, 3, 200, 0.4, 200, 0.8, 250))

save(SB39, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial39\\SB39.RData")




# Change 3 qs complete data -------------------------------------------------------------------------------------------
# Increase
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial128")
SB128 = DFA_Sim_Full(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q3,
                     q_survey2 = q0,
                     q_survey3 = q9,
                     q_survey4 = q0,
                     q_survey5 = q9,
                     q_survey6 = q0,
                     q_survey7 = q0,
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
                     c=c(40, 1.3, 3, 6, 350, 5.5, 3.5))
save(SB128, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial128\\SB128.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial129")
SB129 = DFA_Sim_Full(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q0,
                     q_survey2 = q3,
                     q_survey3 = q0,
                     q_survey4 = q9,
                     q_survey5 = q0,
                     q_survey6 = q9,
                     q_survey7 = q0,
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
                     c=c(2.5, 3, 2.5, 0.4, 15, 0.6, 2) )#c=c(3.5, 0.7, 4.5, 0.4, 6, 0.6, 3.5))
save(SB129, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial129\\SB129.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial130")
SB130 = DFA_Sim_Full(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q0,
                     q_survey2 = q0,
                     q_survey3 = q3,
                     q_survey4 = q0,
                     q_survey5 = q9,
                     q_survey6 = q0,
                     q_survey7 = q9,
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
                     c=c(3.5, 1.3, 3, 6, 350, 5.5, 1.5))
save(SB130, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial130\\SB130.RData")

# decrease
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial131")
SB131 = DFA_Sim_Full(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q4,
                     q_survey2 = q0,
                     q_survey3 = q10,
                     q_survey4 = q0,
                     q_survey5 = q10,
                     q_survey6 = q0,
                     q_survey7 = q0,
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
                     c=c(3.5, 1.3, 3, 6, 350, 5.5, 1.5))
save(SB131, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial131\\SB131.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial132")
SB132 = DFA_Sim_Full(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q0,
                     q_survey2 = q4,
                     q_survey3 = q0,
                     q_survey4 = q10,
                     q_survey5 = q0,
                     q_survey6 = q10,
                     q_survey7 = q0,
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
                     c=c(3.5, 1, 4.5, 200, 6, 300, 3))
save(SB132, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial132\\SB132.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial133")
SB133 = DFA_Sim_Full(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q0,
                     q_survey2 = q0,
                     q_survey3 = q4,
                     q_survey4 = q0,
                     q_survey5 = q10,
                     q_survey6 = q0,
                     q_survey7 = q10,
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
                     c=c(3.5, 1, 3, 5, 10, 4.5, 30))
save(SB133, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial133\\SB133.RData")

# D - I2
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial134")
SB134 = DFA_Sim_Full(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q4,
                     q_survey2 = q0,
                     q_survey3 = q9,
                     q_survey4 = q0,
                     q_survey5 = q9,
                     q_survey6 = q0,
                     q_survey7 = q0,
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
                     c=c(45, 1.3, 35, 6, 12, 5.5, 3))
save(SB134, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial134\\SB134.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial135")
SB135 = DFA_Sim_Full(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q0,
                     q_survey2 = q4,
                     q_survey3 = q0,
                     q_survey4 = q9,
                     q_survey5 = q0,
                     q_survey6 = q9,
                     q_survey7 = q0,
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
                     c=c(3.5, 0.6, 4.5, 200, 6, 300, 3))
save(SB135, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial135\\SB135.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial136")
SB136 = DFA_Sim_Full(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q0,
                     q_survey2 = q0,
                     q_survey3 = q4,
                     q_survey4 = q0,
                     q_survey5 = q9,
                     q_survey6 = q0,
                     q_survey7 = q9,
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
                     c=c(3.5, 1, 2.5, 5, 10, 4.5, 40) )
save(SB136, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial136\\SB136.RData")

# I - D2
setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial137")
SB137 = DFA_Sim_Full(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q3,
                     q_survey2 = q0,
                     q_survey3 = q10,
                     q_survey4 = q0,
                     q_survey5 = q10,
                     q_survey6 = q0,
                     q_survey7 = q0,
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
                     c=c(20, 1.3, 3.5, 6, 450, 5, 3))
save(SB137, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial137\\SB137.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial138")
SB138 = DFA_Sim_Full(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q0,
                     q_survey2 = q3,
                     q_survey3 = q0,
                     q_survey4 = q10,
                     q_survey5 = q0,
                     q_survey6 = q10,
                     q_survey7 = q0,
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
                     c=c(3.5, 0.7, 4.5, 200, 6, 300, 3) )
save(SB138, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial138\\SB138.RData")

setwd("D:\\vspace1\\DFA_Simulation\\SB\\Trial139")
SB139 = DFA_Sim_Full(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q0,
                     q_survey2 = q0,
                     q_survey3 = q3,
                     q_survey4 = q0,
                     q_survey5 = q10,
                     q_survey6 = q0,
                     q_survey7 = q10,
                     Iay_eq = Iay_eq_F1,
                     Niters=100,
                     c=c(3.5, 1.3, 3.5, 6, 400, 5, 1.5))


save(SB139, file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial139\\SB139.RData")



######## PLOT COMPLETE VS. MISSING DATA DFA TRENDS #########


load(file = "D:\\vspace1\\DFA_Simulation\\SB\\Trial1\\DFA_Results.RData")
names(DFA_Results)
DFA_Results$DFATrends


for(i in 1:39){
  # for(i in 16:39){
  load(file = paste("D:\\vspace1\\DFA_Simulation\\SB\\Trial",i,"\\DFA_Results.RData", sep=""))
  # assign(paste("DFA_",i,sep=""),DFA_Results$DFATrends)
  assign("x",DFA_Results$DFATrendsBT)
  
  j=i+100
  load(file = paste("D:\\vspace1\\DFA_Simulation\\SB\\Trial",j,"\\DFA_Results.RData", sep=""))
  assign("y",DFA_Results$DFATrendsBT)
  
  png(filename=paste("D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Compare_MissData\\Compare_",i,".png", sep=""),
      type="cairo",
      units="mm",
      width=300,
      height=300,
      pointsize=12,
      res=400)
  par(mfrow=c(10,10),  mar=c(0.1, 0.1, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
  for(l in 1:100){
    plot(1:40, (y[l,]-mean(y[l,]))/sd(y[l,]),lwd=2,type='l',col='black', ylim=c(-3,3), axes=F,ylab="", xlab="")
    axis(1,labels=FALSE)
    axis(2,labels=FALSE)
    lines(1:40, (x[l,]-mean(x[l,]))/sd(x[l,]), type='l', lwd=2,col='red')
    text(8,-2.5, labels = format( mean( get(paste("SB",j,sep=""))$FitRatios[l,]), digits=3 ) , col='black', cex=1.1)
    text(33,-2.5, labels = format( mean( get(paste("SB",i,sep=""))$FitRatios[l,]), digits=3 ) , col='red', cex=1.1)
    
  }
  dev.off()
  
}


cbind(apply(SB109$FitRatios, 1, mean), apply(SB9$FitRatios, 1, mean))

#BAD: 9, 16, 22
cbind(apply(SB109$FitRatios, 1, mean), apply(SB9$FitRatios, 1, mean))

#Potentially problematic: 5, 10, 12, (17) 






#######################################################################################################


###################
# EXAMPLE #
###################
#
 #### 128 ####
setwd("C:/~")
getwd() 
SB128 = DFA_Sim_Full(Tot_F = F1,
                     CV1=0.38, CV2=0.48, CV3=0.65, CV4=0.24, CV5=0.30, CV6=0.36, CV7=0.40,
                     q_survey1 = q3,
                     q_survey2 = q0,
                     q_survey3 = q9,
                     q_survey4 = q0,
                     q_survey5 = q9,
                     q_survey6 = q0,
                     q_survey7 = q0,
                     Iay_eq = Iay_eq_F1,
                     Niters=1,
                     c=c(40, 1.3, 3, 6, 350, 5.5, 3.5))
Tot_F = F1;
CV1=0.38; CV2=0.48; CV3=0.65; CV4=0.24; CV5=0.30; CV6=0.36; CV7=0.40;
q_survey1 = q3;
q_survey2 = q0;
q_survey3 = q9;
q_survey4 = q0;
q_survey5 = q9;
q_survey6 = q0;
q_survey7 = q0;
Iay_eq = Iay_eq_F1;
Niters=1;
c=c(40, 1.3, 3, 6, 350, 5.5, 3.5)
c=c(40, 1.3, 3, 6, 350, 5.5, 3.5)
# RUN MANUALLY...

yrs = 51:90

SB128
SB_128_List <- list()
SB_128_List[["Iy1"]]<- Iy1
SB_128_List[["Iy2"]]<- Iy2
SB_128_List[["Iy3"]]<- Iy3
SB_128_List[["Iy4"]]<- Iy4
SB_128_List[["Iy5"]]<- Iy5
SB_128_List[["Iy6"]]<- Iy6
SB_128_List[["Iy7"]]<- Iy7
SB_128_List[["Z.rot.1"]]<- Z.rot.1
SB_128_List[["trends.rot.1"]]<- trends.rot.1

SB_128_List[["index1"]]<- index1
SB_128_List[["lowerCI1"]]<- lowerCI1
SB_128_List[["upperCI1"]]<- upperCI1
SB_128_List[["indexBT1"]]<- indexBT1
SB_128_List[["indexSEBT1"]]<- indexSEBT1

save(SB_128_List, file="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\MS1\\SB_128_List.RData")



png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\MS1\\SB128_1_RawIndices.png", 
    type="cairo",
    units="mm", 
    width=85, 
    height=70, 
    pointsize=8, 
    res=600)
#####
par(mfrow=c(1,1), mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0))
l=1

plot(51:90, Iy1[l,51:90]/mean( Iy1[l,51:90]), type='l', col='grey', xlim=c(51, 90), 
     ylim=c(0, 6), lwd=1, ylab='',xlab='', lty=2)
lines(51:90, Iy2[l,51:90]/mean( Iy2[l,51:90]), type='l', col='blue', lwd=1)
lines(51:90, Iy4[l,51:90]/mean( Iy4[l,51:90]), type='l', col='green', lwd=1)
lines(51:90, Iy6[l,51:90]/mean( Iy6[l,51:90]), type='l', col='mediumorchid2', lwd=1)
lines(51:90, Iy7[l,51:90]/mean( Iy7[l,51:90]), type='l', col='cyan', lwd=1)
lines(51:90, Iy3[l,51:90]/mean( Iy3[l,51:90]), type='l', col='red', lwd=1, lty=2)
lines(51:90, Iy5[l,51:90]/mean( Iy5[l,51:90]), type='l', col='orange', lwd=1, lty=2)

mtext('Simulated Year', 1, cex=1, line=1.1)
mtext('Mean-standardized index', 2, cex=1, line=1.1)
mtext('Iteration 1 of Trial SB128', 3, cex=1, line=-1.1)


# title(paste("Iteration",l, sep=" "))
#####
dev.off() 




png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\MS1\\SB128_1_RawIndices_a.png", 
    type="cairo",
    units="mm", 
    width=85, 
    height=70, 
    pointsize=8, 
    res=600)
#####
par(mfrow=c(1,1), mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0))
l=1

plot(51:90, Iy1[l,51:90]/mean( Iy1[l,51:90]), type='l', col='grey', xlim=c(51, 90), 
     ylim=c(0, 15), lwd=1, ylab='',xlab='', lty=2)
lines(51:90, Iy2[l,51:90]/mean( Iy2[l,51:90]), type='l', col='blue', lwd=1)
lines(51:90, Iy4[l,51:90]/mean( Iy4[l,51:90]), type='l', col='green', lwd=1)
lines(51:90, Iy6[l,51:90]/mean( Iy6[l,51:90]), type='l', col='mediumorchid2', lwd=1)
lines(51:90, Iy7[l,51:90]/mean( Iy7[l,51:90]), type='l', col='cyan', lwd=1)
lines(51:90, Iy3[l,51:90]/mean( Iy3[l,51:90]), type='l', col='red', lwd=1, lty=2)
lines(51:90, Iy5[l,51:90]/mean( Iy5[l,51:90]), type='l', col='orange', lwd=1, lty=2)

mtext('Simulated Year', 1, cex=1, line=1.1)
mtext('Mean-standardized index', 2, cex=1, line=1.1)
mtext('Iteration 1 of Trial SB128', 3, cex=1, line=-1.1)

par(fig = c(0.07,0.7, 0.25, 0.925), new = T)  
plot(51:90, Iy1[l,51:90]/mean( Iy1[l,51:90]), type='l', col='grey', xlim=c(51, 90), 
     ylim=c(0, 3.5), lwd=1, ylab='',xlab='', lty=2)
lines(51:90, Iy2[l,51:90]/mean( Iy2[l,51:90]), type='l', col='blue', lwd=1)
lines(51:90, Iy4[l,51:90]/mean( Iy4[l,51:90]), type='l', col='green', lwd=1)
lines(51:90, Iy6[l,51:90]/mean( Iy6[l,51:90]), type='l', col='mediumorchid2', lwd=1)
lines(51:90, Iy7[l,51:90]/mean( Iy7[l,51:90]), type='l', col='cyan', lwd=1)
lines(51:90, Iy3[l,51:90]/mean( Iy3[l,51:90]), type='l', col='red', lwd=1, lty=2)
lines(51:90, Iy5[l,51:90]/mean( Iy5[l,51:90]), type='l', col='orange', lwd=1, lty=2)

# title(paste("Iteration",l, sep=" "))
#####
dev.off() 



png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\SB128_1_FactorLoadings.png", 
    type="cairo",
    units="mm", 
    width=85, 
    height=70, 
    pointsize=6, 
    res=600)
#####
par(mfrow=c(1,1), mar=c(0.3, 2.5, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)

survey = rownames(dat.z)
minZ = 0.00
for(l in 1:Niters){
  Z.rot = get(paste("Z.rot",l,sep="."))
  trends.rot = get(paste("trends.rot",l,sep="."))
  ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
  # par(mfrow=c(nrow(trends.rot),1), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
  for(h in 1:nrow(trends.rot)) {
    plot(c(1:N.ts)[abs(Z.rot[,h])>minZ], as.vector(Z.rot[abs(Z.rot[,h])>minZ,h]),
         type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1), axes=FALSE,
         col=c('grey','blue','red','green','orange','mediumorchid2','cyan'))
    axis(1, labels=FALSE)
    axis(2)
      for(j in 1:N.ts) {
        if(Z.rot[j,h] > minZ) {text(j, -0.02, paste(survey[j],": ", round(Z.rot[j],3),sep=""), srt=90, adj=1, cex=0.9)}
        if(Z.rot[j,h] < -minZ) {text(j, 0.02, paste(survey[j],": ", round(Z.rot[j],3),sep=""), srt=90, adj=0, cex=0.9)}
      } # end j loop
      mtext('Iteration 1 of Trial SB128',side=3,line=-1.3)
      mtext('Factor loadings on trend 1', side=2, line=1.2)
    abline(h=0, lwd=1, col="black")
    abline(h=0.2, col='gray', lty=2)
    abline(h=-0.2, col='gray', lty=2)
    box()
  } # end h loop
}
#####
dev.off()



png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\SB128_1_DFATrend.png", 
    type="cairo",
    units="mm", 
    width=85, 
    height=70, 
    pointsize=6, 
    res=600)
#####
par(mfrow=c(1,1), mar=c(2.1, 2.5, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)

for(l in 1:Niters){
  index = get(paste("index",l,sep=""))
  lowerCI = get(paste("lowerCI",l,sep=""))
  upperCI = get(paste("upperCI",l,sep=""))
  df <- data.frame(yrs, index, lowerCI, upperCI)
  with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)) ))
  abline(h=0)
  # axis(1,labels=FALSE)
  # axis(2,labels=FALSE)
  x <- as.numeric(df$yrs)
  polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='lightblue', border=NA)
  with(df, lines(x, index, type='l', lwd=1, pch=16, col="deepskyblue4"))
  # title(paste("Iteration",l,sep=" "))
  NL = get(paste("NL",l,sep=""))
  lines(yrs, rescale(NL, to=c(min(index), max(index))), type='l',col="red", lwd=2, lty=2)
  mtext('Iteration 1 of Trial SB128',side=3,line=-1.1)
  mtext('DFA trend', side=2, line=1.1)
  mtext('Simulated year', 1, line=1.1)
  # text(x = min(yrs)+4, y = 0.8*min(lowerCI), labels = format( FR[l], digits=3 ) , cex=1)
  legend("bottomleft",c("DFA trend","95% CI","Rescaled simulated abundance"),
         pch=c(0,15,0), pt.cex = c(0,3,0), lwd=c(1,0,1), lty=c(1,0,2), col=c('deepskyblue4','lightblue','red'), cex=0.8,
         bty='n')
}
#####

dev.off()






png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\SB128_1_Backtransformed_DFATrend.png", 
    type="cairo",
    units="mm", 
    width=85, 
    height=70, 
    pointsize=6, 
    res=600)
##########
par(mfrow=c(1,1), mar=c(2.1, 2.5, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)
for(l in 1:Niters){
  indexBT = get(paste("indexBT",l,sep=""))
  indexSEBT = get(paste("indexSEBT",l,sep=""))
  lowerCI = indexBT-(1.96*indexSEBT)
  upperCI = indexBT+(1.96*indexSEBT)
  df <- data.frame(yrs, indexBT, lowerCI, upperCI)
  with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI))))
  axis(1,labels=FALSE)
  axis(2,labels=FALSE)
  x <- as.numeric(df$yrs)
  polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
  with(df, lines(x, indexBT, type='l', lwd=1, pch=16))
  # title(paste("Iteration",l,sep=" "))
  # abline(h=0)
  N = get(paste("N",l,sep=""))
  lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT)+0.05)), type='l',col="red", lwd=2, lty=2)
  # lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT)+0.05)), type='l',col="red", lwd=2, lty=2)
  # lines(yrs, rescale(N, to=c(min(indexBT)-0.15, max(indexBT)+0.15)), type='l',col="red", lwd=2, lty=2)
  # lines(yrs, rescale(N, to=c(min(indexBT)-0.1, max(indexBT)+0.1)), type='l',col="blue", lwd=2, lty=2)
  # lines(yrs, rescale(N, to=c(min(indexBT)-0, max(indexBT)+0)), type='l',col="green", lwd=2, lty=2)
  # lines(yrs, rescale(N, to=c(min(indexBT)-0.1, max(indexBT)+0.3)), type='l',col="orange", lwd=2, lty=2)
  # lines(yrs, rescale(N, to=c(min(indexBT)-0.2, max(indexBT)+0.2)), type='l',col="orange", lwd=2, lty=2)
  mtext('Iteration 1 of Trial SB128',side=3,line=-1.1)
  mtext('Back-transformed DFA trend', side=2, line=1.1)
  mtext('Simulated year', 1, line=1.1)
  legend("bottomleft",c("DFA trend","95% CI","Rescaled simulated abundance"),
         pch=c(0,15,0), pt.cex = c(0,3,0), lwd=c(1,0,1), lty=c(1,0,2), col=c('black','grey','red'), cex=0.8,
         bty='n')
}

#########
dev.off()


png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\SB128_1_Backtransformed_DFATrend_b.png", 
    type="cairo",
    units="mm", 
    width=85, 
    height=70, 
    pointsize=6, 
    res=600)
########
par(mfrow=c(1,1), mar=c(2.1, 2.5, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)
for(l in 1:Niters){
  indexBT = get(paste("indexBT",l,sep=""))
  indexSEBT = get(paste("indexSEBT",l,sep=""))
  lowerCI = indexBT-(1.96*indexSEBT)
  upperCI = indexBT+(1.96*indexSEBT)
  df <- data.frame(yrs, indexBT, lowerCI, upperCI)
  with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI))))
  axis(1,labels=FALSE)
  axis(2,labels=FALSE)
  x <- as.numeric(df$yrs)
  polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
  with(df, lines(x, indexBT, type='l', lwd=1, pch=16))
  # title(paste("Iteration",l,sep=" "))
  # abline(h=0)
  N = get(paste("N",l,sep=""))
  # lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT)+0.05)), type='l',col="red", lwd=2, lty=2)
  # lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT)+0.05)), type='l',col="red", lwd=2, lty=2)
  # lines(yrs, rescale(N, to=c(min(indexBT)-0.15, max(indexBT)+0.15)), type='l',col="red", lwd=2, lty=2)
  # lines(yrs, rescale(N, to=c(min(indexBT)-0.1, max(indexBT)+0.1)), type='l',col="blue", lwd=2, lty=2)
  # lines(yrs, rescale(N, to=c(min(indexBT)-0, max(indexBT)+0)), type='l',col="green", lwd=2, lty=2)
  # lines(yrs, rescale(N, to=c(min(indexBT)-0.1, max(indexBT)+0.3)), type='l',col="orange", lwd=2, lty=2)
  # lines(yrs, rescale(N, to=c(min(indexBT)-0.2, max(indexBT)+0.2)), type='l',col="orange", lwd=2, lty=2)
  mtext('Iteration 1 of Trial SB128',side=3,line=-1.1)
  mtext('Back-transformed DFA trend', side=2, line=1.1)
  mtext('Simulated year', 1, line=1.1)
  legend("bottomleft",c("DFA trend","95% CI","Rescaled simulated abundance"),
         pch=c(0,15,0), pt.cex = c(0,3,0), lwd=c(1,0,1), lty=c(1,0,2), col=c('black','grey','red'), cex=0.8,
         bty='n')
}
###########

dev.off()



png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\MS1\\SB128_1_ALL.png", 
    type="cairo",
    units="mm", 
    width=170, 
    height=140, 
    pointsize=8, 
    res=600)
############
par(mfrow=c(2,2), mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0))
l=1

plot(51:90, Iy1[l,51:90]/mean( Iy1[l,51:90]), type='l', col='grey', xlim=c(51, 90), 
     ylim=c(0, 6), lwd=1, ylab='',xlab='', lty=2)
lines(51:90, Iy2[l,51:90]/mean( Iy2[l,51:90]), type='l', col='blue', lwd=1)
lines(51:90, Iy3[l,51:90]/mean( Iy3[l,51:90]), type='l', col='red', lwd=1, lty=2)
lines(51:90, Iy4[l,51:90]/mean( Iy4[l,51:90]), type='l', col='green', lwd=1)
lines(51:90, Iy5[l,51:90]/mean( Iy5[l,51:90]), type='l', col='orange', lwd=1, lty=2)
lines(51:90, Iy6[l,51:90]/mean( Iy6[l,51:90]), type='l', col='mediumorchid2', lwd=1)
lines(51:90, Iy7[l,51:90]/mean( Iy7[l,51:90]), type='l', col='cyan', lwd=1)
mtext('Simulated Year', 1, cex=1, line=1.1)
mtext('Mean-standardized index', 2, cex=1, line=1.1)
# mtext('Iteration 1 of Trial SB128', 3, cex=1, line=-1.1)
# title(paste("Iteration",l, sep=" "))

text(par("usr")[1]-1,par("usr")[4], "a.", srt = 0, xpd = TRUE, adj=1.75, font=2, cex=1.2)


par( mar=c(0.3, 2.5, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)

survey = rownames(dat.z)
minZ = 0.00
for(l in 1:Niters){
  Z.rot = get(paste("Z.rot",l,sep="."))
  trends.rot = get(paste("trends.rot",l,sep="."))
  ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
  # par(mfrow=c(nrow(trends.rot),1), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
  for(h in 1:nrow(trends.rot)) {
    plot(c(1:N.ts)[abs(Z.rot[,h])>minZ], as.vector(Z.rot[abs(Z.rot[,h])>minZ,h]),
         type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1),
         col=c('grey','blue','red','green','orange','mediumorchid2','cyan'))
    axis(1, labels=FALSE)
    axis(2)
    for(j in 1:N.ts) {
      if(Z.rot[j,h] > minZ) {text(j, -0.025, paste(survey[j],": ", round(Z.rot[j],3),sep=""), srt=90, adj=1, cex=0.9)}
      if(Z.rot[j,h] < -minZ) {text(j, 0.025, paste(survey[j],": ", round(Z.rot[j],3),sep=""), srt=90, adj=0, cex=0.9)}
    } # end j loop
    # mtext('Iteration 1 of Trial SB128',side=3,line=-1.3)
    mtext('Factor loadings on trend 1', side=2, line=1.2)
    abline(h=0, lwd=1, col="black")
    abline(h=0.2, col='gray', lty=2)
    abline(h=-0.2, col='gray', lty=2)
  } # end h loop
}

text(par("usr")[1]-0.216,par("usr")[4], "b.", srt = 0, xpd = TRUE, adj=1.75, font=2, cex=1.2)


par( mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)

for(l in 1:Niters){
  index = get(paste("index",l,sep=""))
  lowerCI = get(paste("lowerCI",l,sep=""))
  upperCI = get(paste("upperCI",l,sep=""))
  df <- data.frame(yrs, index, lowerCI, upperCI)
  with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)) ))
  abline(h=0)
  # axis(1,labels=FALSE)
  # axis(2,labels=FALSE)
  x <- as.numeric(df$yrs)
  polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='lightblue', border=NA)
  with(df, lines(x, index, type='l', lwd=1, pch=16, col="deepskyblue4"))
  # title(paste("Iteration",l,sep=" "))
  NL = get(paste("NL",l,sep=""))
  lines(yrs, rescale(NL, to=c(min(index), max(index))), type='l',col="red", lwd=2, lty=2)
  # mtext('Iteration 1 of Trial SB128',side=3,line=-1.1)
  mtext('DFA trend', side=2, line=1.1)
  mtext('Simulated year', 1, line=1.1)
  # text(x = min(yrs)+4, y = 0.8*min(lowerCI), labels = format( FR[l], digits=3 ) , cex=1)
  legend("bottomleft",c("DFA trend","95% CI","Rescaled simulated abundance"),
         pch=c(0,15,0), pt.cex = c(0,3,0), lwd=c(1,0,1), lty=c(1,0,2), col=c('deepskyblue4','lightblue','red'), cex=0.8,
         bty='n')
  box()
}

text(par("usr")[1]-1,par("usr")[4], "c.", srt = 0, xpd = TRUE, adj=1.65, font=2, cex=1.2)


par( mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)
for(l in 1:Niters){
  indexBT = get(paste("indexBT",l,sep=""))
  indexSEBT = get(paste("indexSEBT",l,sep=""))
  lowerCI = indexBT-(1.96*indexSEBT)
  upperCI = indexBT+(1.96*indexSEBT)
  df <- data.frame(yrs, indexBT, lowerCI, upperCI)
  with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI))))
  axis(1,labels=FALSE)
  axis(2,labels=FALSE)
  x <- as.numeric(df$yrs)
  polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
  with(df, lines(x, indexBT, type='l', lwd=1, pch=16))
  # title(paste("Iteration",l,sep=" "))
  # abline(h=0)
  N = get(paste("N",l,sep=""))
  lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT))), type='l',col="red", lwd=2, lty=2)
  # mtext('Iteration 1 of Trial SB128',side=3,line=-1.1)
  mtext('Back-transformed DFA trend', side=2, line=1.1)
  mtext('Simulated year', 1, line=1.1)
  legend("bottomleft",c("DFA trend","95% CI","Rescaled simulated abundance"),
         pch=c(0,15,0), pt.cex = c(0,3,0), lwd=c(1,0,1), lty=c(1,0,2), col=c('black','grey','red'), cex=0.8,
         bty='n')
  box()
}

text(par("usr")[1]-1,par("usr")[4], "d.", srt = 0, xpd = TRUE, adj=1.65, font=2, cex=1.2)
#####

dev.off()


png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\MS1\\SB128_1_ALL_a.png", 
    type="cairo",
    units="mm", 
    width=170, 
    height=140, 
    pointsize=8, 
    res=600)
#######
par(mfrow=c(2,2), mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0))
l=1

plot(51:90, Iy1[l,51:90]/mean( Iy1[l,51:90]), type='l', col='grey', xlim=c(51, 90), 
     ylim=c(0, 15), lwd=1, ylab='',xlab='', lty=2)
lines(51:90, Iy2[l,51:90]/mean( Iy2[l,51:90]), type='l', col='blue', lwd=1)
lines(51:90, Iy4[l,51:90]/mean( Iy4[l,51:90]), type='l', col='green', lwd=1)
lines(51:90, Iy6[l,51:90]/mean( Iy6[l,51:90]), type='l', col='mediumorchid2', lwd=1)
lines(51:90, Iy7[l,51:90]/mean( Iy7[l,51:90]), type='l', col='cyan', lwd=1)
lines(51:90, Iy3[l,51:90]/mean( Iy3[l,51:90]), type='l', col='red', lwd=1, lty=2)
lines(51:90, Iy5[l,51:90]/mean( Iy5[l,51:90]), type='l', col='orange', lwd=1, lty=2)

mtext('Simulated Year', 1, cex=1, line=1.1)
mtext('Mean-standardized index', 2, cex=1, line=1.1)
mtext('Iteration 1 of Trial SB128', 3, cex=1, line=-1.1)


text(par("usr")[1]-1,par("usr")[4], "a.", srt = 0, xpd = TRUE, adj=1.75, font=2, cex=1.2)

par( mar=c(0.3, 2.5, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)

survey = rownames(dat.z)
minZ = 0.00
for(l in 1:Niters){
  Z.rot = get(paste("Z.rot",l,sep="."))
  trends.rot = get(paste("trends.rot",l,sep="."))
  ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
  # par(mfrow=c(nrow(trends.rot),1), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
  for(h in 1:nrow(trends.rot)) {
    plot(c(1:N.ts)[abs(Z.rot[,h])>minZ], as.vector(Z.rot[abs(Z.rot[,h])>minZ,h]),
         type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1),
         col=c('grey','blue','red','green','orange','mediumorchid2','cyan'))
    axis(1, labels=FALSE)
    axis(2)
    for(j in 1:N.ts) {
      if(Z.rot[j,h] > minZ) {text(j, -0.025, paste(survey[j],": ", round(Z.rot[j],3),sep=""), srt=90, adj=1, cex=0.9)}
      if(Z.rot[j,h] < -minZ) {text(j, 0.025, paste(survey[j],": ", round(Z.rot[j],3),sep=""), srt=90, adj=0, cex=0.9)}
    } # end j loop
    # mtext('Iteration 1 of Trial SB128',side=3,line=-1.3)
    mtext('Factor loadings on trend 1', side=2, line=1.2)
    abline(h=0, lwd=1, col="black")
    abline(h=0.2, col='gray', lty=2)
    abline(h=-0.2, col='gray', lty=2)
  } # end h loop
}
mtext('Iteration 1 of Trial SB128', 3, cex=1, line=-1.1)

text(par("usr")[1]-0.216,par("usr")[4], "b.", srt = 0, xpd = TRUE, adj=1.75, font=2, cex=1.2)


par( mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)

for(l in 1:Niters){
  index = get(paste("index",l,sep=""))
  lowerCI = get(paste("lowerCI",l,sep=""))
  upperCI = get(paste("upperCI",l,sep=""))
  df <- data.frame(yrs, index, lowerCI, upperCI)
  with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)) ))
  abline(h=0)
  # axis(1,labels=FALSE)
  # axis(2,labels=FALSE)
  x <- as.numeric(df$yrs)
  polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='lightblue', border=NA)
  with(df, lines(x, index, type='l', lwd=1, pch=16, col="deepskyblue4"))
  # title(paste("Iteration",l,sep=" "))
  NL = get(paste("NL",l,sep=""))
  lines(yrs, rescale(NL, to=c(min(index), max(index))), type='l',col="red", lwd=2, lty=2)
  # mtext('Iteration 1 of Trial SB128',side=3,line=-1.1)
  mtext('DFA trend', side=2, line=1.1)
  mtext('Simulated year', 1, line=1.1)
  # text(x = min(yrs)+4, y = 0.8*min(lowerCI), labels = format( FR[l], digits=3 ) , cex=1)
  legend("bottomleft",c("DFA trend","95% CI","Rescaled simulated abundance"),
         pch=c(0,15,0), pt.cex = c(0,3,0), lwd=c(1,0,1), lty=c(1,0,2), col=c('deepskyblue4','lightblue','red'), cex=0.8,
         bty='n')
  box()
}
mtext('Iteration 1 of Trial SB128', 3, cex=1, line=-1.1)

text(par("usr")[1]-1,par("usr")[4], "c.", srt = 0, xpd = TRUE, adj=1.65, font=2, cex=1.2)


par( mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)
for(l in 1:Niters){
  indexBT = get(paste("indexBT",l,sep=""))
  indexSEBT = get(paste("indexSEBT",l,sep=""))
  lowerCI = indexBT-(1.96*indexSEBT)
  upperCI = indexBT+(1.96*indexSEBT)
  df <- data.frame(yrs, indexBT, lowerCI, upperCI)
  with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI))))
  axis(1,labels=FALSE)
  axis(2,labels=FALSE)
  x <- as.numeric(df$yrs)
  polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
  with(df, lines(x, indexBT, type='l', lwd=1, pch=16))
  # title(paste("Iteration",l,sep=" "))
  # abline(h=0)
  N = get(paste("N",l,sep=""))
  lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT))), type='l',col="red", lwd=2, lty=2)
  # mtext('Iteration 1 of Trial SB128',side=3,line=-1.1)
  mtext('Back-transformed DFA trend', side=2, line=1.1)
  mtext('Simulated year', 1, line=1.1)
  legend("bottomleft",c("DFA trend","95% CI","Rescaled simulated abundance"),
         pch=c(0,15,0), pt.cex = c(0,3,0), lwd=c(1,0,1), lty=c(1,0,2), col=c('black','grey','red'), cex=0.8,
         bty='n')
  box()
}
mtext('Iteration 1 of Trial SB128', 3, cex=1, line=-1.1)

text(par("usr")[1]-1,par("usr")[4], "d.", srt = 0, xpd = TRUE, adj=1.65, font=2, cex=1.2)


par(fig = c(0.03,0.35, 0.625, 0.96), new = T)  
plot(51:90, Iy1[l,51:90]/mean( Iy1[l,51:90]), type='l', col='grey', xlim=c(51, 90), 
     ylim=c(0, 3.5), lwd=1, ylab='',xlab='', lty=2)
lines(51:90, Iy2[l,51:90]/mean( Iy2[l,51:90]), type='l', col='blue', lwd=1)
lines(51:90, Iy4[l,51:90]/mean( Iy4[l,51:90]), type='l', col='green', lwd=1)
lines(51:90, Iy6[l,51:90]/mean( Iy6[l,51:90]), type='l', col='mediumorchid2', lwd=1)
lines(51:90, Iy7[l,51:90]/mean( Iy7[l,51:90]), type='l', col='cyan', lwd=1)
lines(51:90, Iy3[l,51:90]/mean( Iy3[l,51:90]), type='l', col='red', lwd=1, lty=2)
lines(51:90, Iy5[l,51:90]/mean( Iy5[l,51:90]), type='l', col='orange', lwd=1, lty=2)

#########
dev.off()


tiff(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\MS1\\Figure_2.tiff", 
    type="cairo",
    units="mm", 
    width=170, 
    height=140, 
    pointsize=8, 
    res=600)
#######
par(mfrow=c(2,2), mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0))
l=1

plot(51:90, Iy1[l,51:90]/mean( Iy1[l,51:90]), type='l', col='grey', xlim=c(51, 90), 
     ylim=c(0, 15), lwd=1, ylab='',xlab='', lty=2)
lines(51:90, Iy2[l,51:90]/mean( Iy2[l,51:90]), type='l', col='blue', lwd=1)
lines(51:90, Iy4[l,51:90]/mean( Iy4[l,51:90]), type='l', col='green', lwd=1)
lines(51:90, Iy6[l,51:90]/mean( Iy6[l,51:90]), type='l', col='mediumorchid2', lwd=1)
lines(51:90, Iy7[l,51:90]/mean( Iy7[l,51:90]), type='l', col='cyan', lwd=1)
lines(51:90, Iy3[l,51:90]/mean( Iy3[l,51:90]), type='l', col='red', lwd=1, lty=2)
lines(51:90, Iy5[l,51:90]/mean( Iy5[l,51:90]), type='l', col='orange', lwd=1, lty=2)

mtext('Simulated Year', 1, cex=1, line=1.1)
mtext('Mean-standardized index', 2, cex=1, line=1.1)
mtext('Iteration 1 of Trial SB128', 3, cex=1, line=-1.1)


text(par("usr")[1]-1,par("usr")[4], "a.", srt = 0, xpd = TRUE, adj=1.75, font=2, cex=1.2)

par( mar=c(0.3, 2.5, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)

survey = rownames(dat.z)
minZ = 0.00
for(l in 1:Niters){
  Z.rot = get(paste("Z.rot",l,sep="."))
  trends.rot = get(paste("trends.rot",l,sep="."))
  ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
  # par(mfrow=c(nrow(trends.rot),1), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
  for(h in 1:nrow(trends.rot)) {
    plot(c(1:N.ts)[abs(Z.rot[,h])>minZ], as.vector(Z.rot[abs(Z.rot[,h])>minZ,h]),
         type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1),
         col=c('grey','blue','red','green','orange','mediumorchid2','cyan'))
    axis(1, labels=FALSE)
    axis(2)
    for(j in 1:N.ts) {
      if(Z.rot[j,h] > minZ) {text(j, -0.025, paste(survey[j],": ", round(Z.rot[j],3),sep=""), srt=90, adj=1, cex=0.9)}
      if(Z.rot[j,h] < -minZ) {text(j, 0.025, paste(survey[j],": ", round(Z.rot[j],3),sep=""), srt=90, adj=0, cex=0.9)}
    } # end j loop
    # mtext('Iteration 1 of Trial SB128',side=3,line=-1.3)
    mtext('Factor loadings on trend 1', side=2, line=1.2)
    abline(h=0, lwd=1, col="black")
    abline(h=0.2, col='gray', lty=2)
    abline(h=-0.2, col='gray', lty=2)
  } # end h loop
}
mtext('Iteration 1 of Trial SB128', 3, cex=1, line=-1.1)

text(par("usr")[1]-0.216,par("usr")[4], "b.", srt = 0, xpd = TRUE, adj=1.75, font=2, cex=1.2)


par( mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)

for(l in 1:Niters){
  index = get(paste("index",l,sep=""))
  lowerCI = get(paste("lowerCI",l,sep=""))
  upperCI = get(paste("upperCI",l,sep=""))
  df <- data.frame(yrs, index, lowerCI, upperCI)
  with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)) ))
  abline(h=0)
  # axis(1,labels=FALSE)
  # axis(2,labels=FALSE)
  x <- as.numeric(df$yrs)
  polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='lightblue', border=NA)
  with(df, lines(x, index, type='l', lwd=1, pch=16, col="deepskyblue4"))
  # title(paste("Iteration",l,sep=" "))
  NL = get(paste("NL",l,sep=""))
  lines(yrs, rescale(NL, to=c(min(index), max(index))), type='l',col="red", lwd=2, lty=2)
  # mtext('Iteration 1 of Trial SB128',side=3,line=-1.1)
  mtext('DFA trend', side=2, line=1.1)
  mtext('Simulated year', 1, line=1.1)
  # text(x = min(yrs)+4, y = 0.8*min(lowerCI), labels = format( FR[l], digits=3 ) , cex=1)
  legend("bottomleft",c("DFA trend","95% CI","Rescaled simulated abundance"),
         pch=c(0,15,0), pt.cex = c(0,3,0), lwd=c(1,0,1), lty=c(1,0,2), col=c('deepskyblue4','lightblue','red'), cex=0.8,
         bty='n')
  box()
}
mtext('Iteration 1 of Trial SB128', 3, cex=1, line=-1.1)

text(par("usr")[1]-1,par("usr")[4], "c.", srt = 0, xpd = TRUE, adj=1.65, font=2, cex=1.2)


par( mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)
for(l in 1:Niters){
  indexBT = get(paste("indexBT",l,sep=""))
  indexSEBT = get(paste("indexSEBT",l,sep=""))
  lowerCI = indexBT-(1.96*indexSEBT)
  upperCI = indexBT+(1.96*indexSEBT)
  df <- data.frame(yrs, indexBT, lowerCI, upperCI)
  with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI))))
  axis(1,labels=FALSE)
  axis(2,labels=FALSE)
  x <- as.numeric(df$yrs)
  polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
  with(df, lines(x, indexBT, type='l', lwd=1, pch=16))
  # title(paste("Iteration",l,sep=" "))
  # abline(h=0)
  N = get(paste("N",l,sep=""))
  lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT))), type='l',col="red", lwd=2, lty=2)
  # mtext('Iteration 1 of Trial SB128',side=3,line=-1.1)
  mtext('Back-transformed DFA trend', side=2, line=1.1)
  mtext('Simulated year', 1, line=1.1)
  legend("bottomleft",c("DFA trend","95% CI","Rescaled simulated abundance"),
         pch=c(0,15,0), pt.cex = c(0,3,0), lwd=c(1,0,1), lty=c(1,0,2), col=c('black','grey','red'), cex=0.8,
         bty='n')
  box()
}
mtext('Iteration 1 of Trial SB128', 3, cex=1, line=-1.1)

text(par("usr")[1]-1,par("usr")[4], "d.", srt = 0, xpd = TRUE, adj=1.65, font=2, cex=1.2)


par(fig = c(0.03,0.35, 0.625, 0.96), new = T)  
plot(51:90, Iy1[l,51:90]/mean( Iy1[l,51:90]), type='l', col='grey', xlim=c(51, 90), 
     ylim=c(0, 3.5), lwd=1, ylab='',xlab='', lty=2)
lines(51:90, Iy2[l,51:90]/mean( Iy2[l,51:90]), type='l', col='blue', lwd=1)
lines(51:90, Iy4[l,51:90]/mean( Iy4[l,51:90]), type='l', col='green', lwd=1)
lines(51:90, Iy6[l,51:90]/mean( Iy6[l,51:90]), type='l', col='mediumorchid2', lwd=1)
lines(51:90, Iy7[l,51:90]/mean( Iy7[l,51:90]), type='l', col='cyan', lwd=1)
lines(51:90, Iy3[l,51:90]/mean( Iy3[l,51:90]), type='l', col='red', lwd=1, lty=2)
lines(51:90, Iy5[l,51:90]/mean( Iy5[l,51:90]), type='l', col='orange', lwd=1, lty=2)

#########
dev.off()


#### 130 ######

# try with 130
Tot_F = F1
CV1=0.38; CV2=0.48; CV3=0.65; CV4=0.24; CV5=0.30; CV6=0.36; CV7=0.40;
q_survey1 = q0
q_survey2 = q0
q_survey3 = q3
q_survey4 = q0
q_survey5 = q9
q_survey6 = q0
q_survey7 = q9
Iay_eq = Iay_eq_F1
Niters=1
c=c(3.5, 1.3, 3, 6, 350, 5.5, 1.5)


png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\MS1\\SB130_1_ALL.png", 
    type="cairo",
    units="mm", 
    width=170, 
    height=140, 
    pointsize=8, 
    res=600)
par(mfrow=c(2,2), mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0))
l=1

plot(51:90, Iy1[l,51:90]/mean( Iy1[l,51:90]), type='l', col='cyan', xlim=c(51, 90), 
     ylim=c(0, 6), lwd=1, ylab='',xlab='', lty=1)
lines(51:90, Iy2[l,51:90]/mean( Iy2[l,51:90]), type='l', col='blue', lwd=1)
lines(51:90, Iy3[l,51:90]/mean( Iy3[l,51:90]), type='l', col='red', lwd=1, lty=2)
lines(51:90, Iy4[l,51:90]/mean( Iy4[l,51:90]), type='l', col='green', lwd=1)
lines(51:90, Iy5[l,51:90]/mean( Iy5[l,51:90]), type='l', col='orange', lwd=1, lty=2)
lines(51:90, Iy6[l,51:90]/mean( Iy6[l,51:90]), type='l', col='mediumorchid2', lwd=1)
lines(51:90, Iy7[l,51:90]/mean( Iy7[l,51:90]), type='l', col='grey', lwd=1, lty=2)
mtext('Simulated Year', 1, cex=1, line=1.1)
mtext('Mean-standardized index', 2, cex=1, line=1.1)
# mtext('Iteration 1 of Trial SB128', 3, cex=1, line=-1.1)
# title(paste("Iteration",l, sep=" "))

text(par("usr")[1]-1,par("usr")[4], "a.", srt = 0, xpd = TRUE, adj=1.75, font=2, cex=1.2)


par( mar=c(0.3, 2.5, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)

survey = rownames(dat.z)
minZ = 0.00
for(l in 1:Niters){
  Z.rot = get(paste("Z.rot",l,sep="."))
  trends.rot = get(paste("trends.rot",l,sep="."))
  ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
  # par(mfrow=c(nrow(trends.rot),1), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
  for(h in 1:nrow(trends.rot)) {
    plot(c(1:N.ts)[abs(Z.rot[,h])>minZ], as.vector(Z.rot[abs(Z.rot[,h])>minZ,h]),
         type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1),
         col=c('cyan','blue','red','green','orange','mediumorchid2','grey'))
    axis(1, labels=FALSE)
    axis(2)
    for(j in 1:N.ts) {
      if(Z.rot[j,h] > minZ) {text(j, -0.025, paste(survey[j],": ", round(Z.rot[j],3),sep=""), srt=90, adj=1, cex=0.9)}
      if(Z.rot[j,h] < -minZ) {text(j, 0.025, paste(survey[j],": ", round(Z.rot[j],3),sep=""), srt=90, adj=0, cex=0.9)}
    } # end j loop
    # mtext('Iteration 1 of Trial SB128',side=3,line=-1.3)
    mtext('Factor loadings on trend 1', side=2, line=1.2)
    abline(h=0, lwd=1, col="black")
    abline(h=0.2, col='gray', lty=2)
    abline(h=-0.2, col='gray', lty=2)
  } # end h loop
}

text(par("usr")[1]-0.216,par("usr")[4], "b.", srt = 0, xpd = TRUE, adj=1.75, font=2, cex=1.2)


par( mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)

for(l in 1:Niters){
  index = get(paste("index",l,sep=""))
  lowerCI = get(paste("lowerCI",l,sep=""))
  upperCI = get(paste("upperCI",l,sep=""))
  df <- data.frame(yrs, index, lowerCI, upperCI)
  with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)) ))
  abline(h=0)
  # axis(1,labels=FALSE)
  # axis(2,labels=FALSE)
  x <- as.numeric(df$yrs)
  polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='lightblue', border=NA)
  with(df, lines(x, index, type='l', lwd=1, pch=16, col="deepskyblue4"))
  # title(paste("Iteration",l,sep=" "))
  NL = get(paste("NL",l,sep=""))
  lines(yrs, rescale(NL, to=c(min(index), max(index)+0.5)), type='l',col="red", lwd=2, lty=2)
  # mtext('Iteration 1 of Trial SB128',side=3,line=-1.1)
  mtext('DFA trend', side=2, line=1.1)
  mtext('Simulated year', 1, line=1.1)
  # text(x = min(yrs)+4, y = 0.8*min(lowerCI), labels = format( FR[l], digits=3 ) , cex=1)
  legend("bottomleft",c("DFA trend","95% CI","Rescaled simulated abundance"),
         pch=c(0,15,0), pt.cex = c(0,3,0), lwd=c(1,0,1), lty=c(1,0,2), col=c('deepskyblue4','lightblue','red'), cex=0.8,
         bty='n')
  box()
}

text(par("usr")[1]-1,par("usr")[4], "c.", srt = 0, xpd = TRUE, adj=1.65, font=2, cex=1.2)


par( mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)
for(l in 1:Niters){
  indexBT = get(paste("indexBT",l,sep=""))
  indexSEBT = get(paste("indexSEBT",l,sep=""))
  lowerCI = indexBT-(1.96*indexSEBT)
  upperCI = indexBT+(1.96*indexSEBT)
  df <- data.frame(yrs, indexBT, lowerCI, upperCI)
  with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI))))
  axis(1,labels=FALSE)
  axis(2,labels=FALSE)
  x <- as.numeric(df$yrs)
  polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
  with(df, lines(x, indexBT, type='l', lwd=1, pch=16))
  # title(paste("Iteration",l,sep=" "))
  # abline(h=0)
  N = get(paste("N",l,sep=""))
  lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT)+0.1)), type='l',col="red", lwd=2, lty=2)
  # mtext('Iteration 1 of Trial SB128',side=3,line=-1.1)
  mtext('Back-transformed DFA trend', side=2, line=1.1)
  mtext('Simulated year', 1, line=1.1)
  legend("bottomleft",c("DFA trend","95% CI","Rescaled simulated abundance"),
         pch=c(0,15,0), pt.cex = c(0,3,0), lwd=c(1,0,1), lty=c(1,0,2), col=c('black','grey','red'), cex=0.8,
         bty='n')
  box()
}

text(par("usr")[1]-1,par("usr")[4], "d.", srt = 0, xpd = TRUE, adj=1.65, font=2, cex=1.2)


dev.off()






# png(filename="Trend_Standardize.png", 
#     type="cairo",
#     units="mm", 
#     width=300, 
#     height=300, 
#     pointsize=12, 
#     res=600)

# par(mfrow=c(10, 10))
# par(mar=c(0.3, 0.3, 0.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
for(l in 1:Niters){
  
  index = get(paste("index",l,sep=""))
  NL = get(paste("NL",l,sep=""))
  plot(yrs, (NL-mean(NL))/sd(NL), type='l',col="red", axes=FALSE, lwd=1.5)
  axis(1, labels=FALSE)
  axis(2, labels=FALSE)
  lines(yrs, (index-mean(index))/sd(index), type='l', lwd=1.5)
  
  text(x = min(yrs)+4, y = 0.8*min((index-mean(index))/sd(index)), labels = format( FR[l], digits=3 ) , cex=1)
  
}

# dev.off()







# plot(51:90, Iy1[l,51:90], type='l', col='grey', xlim=c(51, 90), 
#      ylim=c(0, 100), lwd=2, ylab='',xlab='', lty=2)
# lines(51:90, Iy2[l,51:90], type='l', col='blue', lwd=2)
# lines(51:90, Iy3[l,51:90], type='l', col='red', lwd=2, lty=2)
# lines(51:90, Iy4[l,51:90], type='l', col='green', lwd=2)
# lines(51:90, Iy5[l,51:90], type='l', col='orange', lwd=2, lty=2)
# lines(51:90, Iy6[l,51:90], type='l', col='mediumorchid2', lwd=2)
# lines(51:90, Iy7[l,51:90], type='l', col='cyan', lwd=2)
# mtext('Year', 1, cex=1, line=1.1)
# mtext('Mean-standardized index', 2, cex=1, line=1.1)
# mtext('Iteration 1 of Trial SB128', 3, cex=1, line=-1.1)
# # title(paste("Iteration",l, sep=" "))

# par(mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0))
# l=1
# 
# plot(51:90, (Iy1[l,51:90] - mean(Iy1[51:90]) )/mean( Iy1[l,51:90]), type='l', col='grey', xlim=c(51, 90), 
#      ylim=c(0, 6), lwd=2, ylab='',xlab='')
# lines(51:90, (Iy2[l,51:90] - mean(Iy2[51:90]) )/mean( Iy2[l,51:90]), type='l', col='blue', lwd=2)
# lines(51:90, (Iy3[l,51:90] - mean(Iy3[51:90]) )/mean( Iy3[l,51:90]), type='l', col='red', lwd=2)
# lines(51:90, (Iy4[l,51:90] - mean(Iy4[51:90]) )/mean( Iy4[l,51:90]), type='l', col='green', lwd=2)
# lines(51:90, (Iy5[l,51:90] - mean(Iy5[51:90]) )/mean( Iy5[l,51:90]), type='l', col='orange', lwd=2)
# lines(51:90, (Iy6[l,51:90] - mean(Iy6[51:90]) )/mean( Iy6[l,51:90]), type='l', col='mediumorchid2', lwd=2)
# lines(51:90, (Iy7[l,51:90] - mean(Iy7[51:90]) )/mean( Iy7[l,51:90]), type='l', col='cyan', lwd=2)
# mtext('Year', 1, cex=1, line=1.1)
# mtext('Mean-standardized index', 2, cex=1, line=1.1)
# mtext('Iteration 1 of Trial SB128', 3, cex=1, line=-1.1)
# # title(paste("Iteration",l, sep=" "))


##### 129 ####

# try with 129
Tot_F = F1
CV1=0.38; CV2=0.48; CV3=0.65; CV4=0.24; CV5=0.30; CV6=0.36; CV7=0.40;
q_survey1 = q0
q_survey2 = q3
q_survey3 = q0
q_survey4 = q9
q_survey5 = q0
q_survey6 = q9
q_survey7 = q0
Iay_eq = Iay_eq_F1
Niters=1
c=c(2.5, 3, 2.5, 0.4, 15, 0.6, 2)


png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\MS1\\SB129_1_ALL.png", 
    type="cairo",
    units="mm", 
    width=170, 
    height=140, 
    pointsize=8, 
    res=600)
par(mfrow=c(2,2), mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0))
l=1

plot(51:90, Iy1[l,51:90]/mean( Iy1[l,51:90]), type='l', col='cyan', xlim=c(51, 90), 
     ylim=c(0, 4.25), lwd=1, ylab='',xlab='', lty=1)
lines(51:90, Iy2[l,51:90]/mean( Iy2[l,51:90]), type='l', col='red', lwd=1)
lines(51:90, Iy3[l,51:90]/mean( Iy3[l,51:90]), type='l', col='blue', lwd=1, lty=2)
lines(51:90, Iy4[l,51:90]/mean( Iy4[l,51:90]), type='l', col='orange', lwd=1)
lines(51:90, Iy5[l,51:90]/mean( Iy5[l,51:90]), type='l', col='green', lwd=1, lty=2)
lines(51:90, Iy6[l,51:90]/mean( Iy6[l,51:90]), type='l', col='grey', lwd=1)
lines(51:90, Iy7[l,51:90]/mean( Iy7[l,51:90]), type='l', col='mediumorchid2', lwd=1, lty=2)
mtext('Simulated Year', 1, cex=1, line=1.1)
mtext('Mean-standardized index', 2, cex=1, line=1.1)
# mtext('Iteration 1 of Trial SB128', 3, cex=1, line=-1.1)
# title(paste("Iteration",l, sep=" "))

text(par("usr")[1]-1,par("usr")[4], "a.", srt = 0, xpd = TRUE, adj=1.75, font=2, cex=1.2)


par( mar=c(0.3, 2.5, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)

survey = rownames(dat.z)
minZ = 0.00
for(l in 1:Niters){
  Z.rot = get(paste("Z.rot",l,sep="."))
  trends.rot = get(paste("trends.rot",l,sep="."))
  ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
  # par(mfrow=c(nrow(trends.rot),1), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
  for(h in 1:nrow(trends.rot)) {
    plot(c(1:N.ts)[abs(Z.rot[,h])>minZ], as.vector(Z.rot[abs(Z.rot[,h])>minZ,h]),
         type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1),
         col=c('cyan','red','blue','orange','green','grey','mediumorchid2'))
    axis(1, labels=FALSE)
    axis(2)
    for(j in 1:N.ts) {
      if(Z.rot[j,h] > minZ) {text(j, -0.025, paste(survey[j],": ", round(Z.rot[j],3),sep=""), srt=90, adj=1, cex=0.9)}
      if(Z.rot[j,h] < -minZ) {text(j, 0.025, paste(survey[j],": ", round(Z.rot[j],3),sep=""), srt=90, adj=0, cex=0.9)}
    } # end j loop
    # mtext('Iteration 1 of Trial SB128',side=3,line=-1.3)
    mtext('Factor loadings on trend 1', side=2, line=1.2)
    abline(h=0, lwd=1, col="black")
    abline(h=0.2, col='gray', lty=2)
    abline(h=-0.2, col='gray', lty=2)
  } # end h loop
}

text(par("usr")[1]-0.216,par("usr")[4], "b.", srt = 0, xpd = TRUE, adj=1.75, font=2, cex=1.2)


par( mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)

for(l in 1:Niters){
  index = get(paste("index",l,sep=""))
  lowerCI = get(paste("lowerCI",l,sep=""))
  upperCI = get(paste("upperCI",l,sep=""))
  df <- data.frame(yrs, index, lowerCI, upperCI)
  with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI)) ))
  abline(h=0)
  # axis(1,labels=FALSE)
  # axis(2,labels=FALSE)
  x <- as.numeric(df$yrs)
  polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='lightblue', border=NA)
  with(df, lines(x, index, type='l', lwd=1, pch=16, col="deepskyblue4"))
  # title(paste("Iteration",l,sep=" "))
  NL = get(paste("NL",l,sep=""))
  # lines(yrs, rescale(NL, to=c(min(index), max(index)+0.5)), type='l',col="red", lwd=2, lty=2)
  # lines(yrs, rescale(NL, to=c(min(index), max(index))), type='l',col="red", lwd=2, lty=2)
  lines(yrs, rescale(NL, to=c(min(index), max(index)-1)), type='l',col="red", lwd=2, lty=2)
  # mtext('Iteration 1 of Trial SB128',side=3,line=-1.1)
  mtext('DFA trend', side=2, line=1.1)
  mtext('Simulated year', 1, line=1.1)
  # text(x = min(yrs)+4, y = 0.8*min(lowerCI), labels = format( FR[l], digits=3 ) , cex=1)
  legend("bottomleft",c("DFA trend","95% CI","Rescaled simulated abundance"),
         pch=c(0,15,0), pt.cex = c(0,3,0), lwd=c(1,0,1), lty=c(1,0,2), col=c('deepskyblue4','lightblue','red'), cex=0.8,
         bty='n')
  box()
}

text(par("usr")[1]-1,par("usr")[4], "c.", srt = 0, xpd = TRUE, adj=1.65, font=2, cex=1.2)


par( mar=c(2.1, 2.1, 0.3, 0.3),tcl = -0.1, mgp = c(2, 0.2, 0), cex=0.9)
for(l in 1:Niters){
  indexBT = get(paste("indexBT",l,sep=""))
  indexSEBT = get(paste("indexSEBT",l,sep=""))
  lowerCI = indexBT-(1.96*indexSEBT)
  upperCI = indexBT+(1.96*indexSEBT)
  df <- data.frame(yrs, indexBT, lowerCI, upperCI)
  with(df, plot(index ~ yrs, type='n', bty='L',xlab='', ylab='', ylim=c(min(lowerCI), max(upperCI))))
  axis(1,labels=FALSE)
  axis(2,labels=FALSE)
  x <- as.numeric(df$yrs)
  polygon(c(x,rev(x)), c(df$lowerCI, rev(df$upperCI)), col='grey', border=NA)
  with(df, lines(x, indexBT, type='l', lwd=1, pch=16))
  # title(paste("Iteration",l,sep=" "))
  # abline(h=0)
  N = get(paste("N",l,sep=""))
  # lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT)+0.1)), type='l',col="red", lwd=2, lty=2)
  # lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT)-0.1)), type='l',col="red", lwd=2, lty=2)
  lines(yrs, rescale(N, to=c(min(indexBT), max(indexBT)-0.15)), type='l',col="red", lwd=2, lty=2)
  # mtext('Iteration 1 of Trial SB128',side=3,line=-1.1)
  mtext('Back-transformed DFA trend', side=2, line=1.1)
  mtext('Simulated year', 1, line=1.1)
  legend("bottomleft",c("DFA trend","95% CI","Rescaled simulated abundance"),
         pch=c(0,15,0), pt.cex = c(0,3,0), lwd=c(1,0,1), lty=c(1,0,2), col=c('black','grey','red'), cex=0.8,
         bty='n')
  box()
}

text(par("usr")[1]-1,par("usr")[4], "d.", srt = 0, xpd = TRUE, adj=1.65, font=2, cex=1.2)


dev.off()








# PLOT  ---------------------------------------------------------------------------------------------------------


png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Backtransformed_CompleteData_FULL.png", 
    type="cairo",
    units="mm", 
    width=400, 
    height=170, 
    pointsize=10, 
    res=600)
#########
# STANDARDIZE AXES
par(mfrow=c(2,1), mar=c(0.6, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll.BT), 
        c(SB109$biasAll.BT), c(SB110$biasAll.BT), c(SB111$biasAll.BT), c(SB112$biasAll.BT), c(SB113$biasAll.BT), c(SB114$biasAll.BT), c(SB115$biasAll.BT),
        c(SB102$biasAll.BT), c(SB103$biasAll.BT), c(SB104$biasAll.BT), c(SB105$biasAll.BT), c(SB106$biasAll.BT), c(SB107$biasAll.BT), c(SB108$biasAll.BT),
        
        c(SB116$biasAll.BT), c(SB117$biasAll.BT), c(SB118$biasAll.BT), c(SB119$biasAll.BT), c(SB120$biasAll.BT), c(SB121$biasAll.BT),  
        c(SB122$biasAll.BT), c(SB123$biasAll.BT), c(SB124$biasAll.BT), c(SB125$biasAll.BT), c(SB126$biasAll.BT), c(SB127$biasAll.BT),
        c(SB128$biasAll.BT), c(SB129$biasAll.BT), c(SB130$biasAll.BT), c(SB131$biasAll.BT), c(SB132$biasAll.BT), c(SB133$biasAll.BT), 
        
        c(SB137$biasAll.BT), c(SB138$biasAll.BT), c(SB139$biasAll.BT),
        c(SB134$biasAll.BT), c(SB135$biasAll.BT), c(SB136$biasAll.BT), 
        
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=FALSE, cex=0.7 , ylim=c(-1.25,2) )
mtext('Back-transformed Annual Bias',side=2,line=1,cex=1.2)
abline(h=0, col='black')
mtext('Sandbar shark: Complete Data', line=-1.4, cex=1.2)

par( mar=c(1.3, 2.1, 0.1, 0.3))
vioplot(SB101$RMSE.BT, 
        SB109$RMSE.BT, SB110$RMSE.BT, SB111$RMSE.BT, SB112$RMSE.BT, SB113$RMSE.BT, SB114$RMSE.BT, SB115$RMSE.BT,
        SB102$RMSE.BT, SB103$RMSE.BT, SB104$RMSE.BT, SB105$RMSE.BT, SB106$RMSE.BT, SB107$RMSE.BT, SB108$RMSE.BT,
        SB116$RMSE.BT, SB117$RMSE.BT, SB118$RMSE.BT, SB119$RMSE.BT, SB120$RMSE.BT, SB121$RMSE.BT,
        SB122$RMSE.BT, SB123$RMSE.BT, SB124$RMSE.BT, SB125$RMSE.BT, SB126$RMSE.BT, SB127$RMSE.BT,
        SB128$RMSE.BT, SB129$RMSE.BT, SB130$RMSE.BT, SB131$RMSE.BT, SB132$RMSE.BT, SB133$RMSE.BT,
         
        SB137$RMSE.BT, SB138$RMSE.BT, SB139$RMSE.BT,
        SB134$RMSE.BT, SB135$RMSE.BT, SB136$RMSE.BT,
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        cex=0.7, ylim=c(0,0.5),
        names=c('SB101','SB102','SB103','SB104','SB105','SB106','SB107','SB108','SB109','SB110','SB111','SB112','SB113','SB114','SB115',
                'SB116','SB117','SB118','SB119','SB120','SB121','SB122','SB123','SB124','SB125','SB126','SB127','SB128','SB129',
                'SB130','SB131','SB132','SB133','SB134','SB135','SB136','SB137','SB138','SB139') )
mtext('Back-transformed RMSE',side=2,line=1,cex=1.2)
# abline(h=mean(SB101$RMSE.BT), col='white')

##############

dev.off()




png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Backtransformed_MissingData_FULL.png", 
    type="cairo",
    units="mm", 
    width=400, 
    height=170, 
    pointsize=10, 
    res=600)
###########
# STANDARDIZE AXES
par(mfrow=c(2,1), mar=c(0.6, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB1$biasAll.BT), 
        c(SB9$biasAll.BT), c(SB10$biasAll.BT), c(SB11$biasAll.BT), c(SB12$biasAll.BT), c(SB13$biasAll.BT), c(SB14$biasAll.BT), c(SB15$biasAll.BT),
        c(SB2$biasAll.BT), c(SB3$biasAll.BT), c(SB4$biasAll.BT), c(SB5$biasAll.BT), c(SB6$biasAll.BT), c(SB7$biasAll.BT), c(SB8$biasAll.BT),
        c(SB16$biasAll.BT), c(SB17$biasAll.BT), c(SB18$biasAll.BT), c(SB19$biasAll.BT), c(SB20$biasAll.BT), c(SB21$biasAll.BT), 
        c(SB22$biasAll.BT), c(SB23$biasAll.BT), c(SB24$biasAll.BT), c(SB25$biasAll.BT), c(SB26$biasAll.BT), c(SB27$biasAll.BT),
        c(SB28$biasAll.BT), c(SB29$biasAll.BT), c(SB30$biasAll.BT), c(SB31$biasAll.BT), c(SB32$biasAll.BT), c(SB33$biasAll.BT), 
         c(SB37$biasAll.BT), c(SB38$biasAll.BT), c(SB39$biasAll.BT),
        c(SB34$biasAll.BT), c(SB35$biasAll.BT), c(SB36$biasAll.BT),
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=FALSE,  cex=0.7, ylim=c(-4,3.5) )
mtext('Back-transformed Annual Bias',side=2,line=1,cex=1.2)
abline(h=0, col='black')
mtext('Sandbar shark: Missing Data', line=-1.4, cex=1.2)

par( mar=c(1.3, 2.1, 0.1, 0.3))
vioplot(SB1$RMSE.BT, 
        SB9$RMSE.BT,  SB10$RMSE.BT, SB11$RMSE.BT, SB12$RMSE.BT, SB13$RMSE.BT, SB14$RMSE.BT, SB15$RMSE.BT,
        SB2$RMSE.BT, SB3$RMSE.BT, SB4$RMSE.BT, SB5$RMSE.BT, SB6$RMSE.BT, SB7$RMSE.BT, SB8$RMSE.BT,
        SB16$RMSE.BT, SB17$RMSE.BT, SB18$RMSE.BT, SB19$RMSE.BT, SB20$RMSE.BT, SB21$RMSE.BT,
        SB22$RMSE.BT, SB23$RMSE.BT, SB24$RMSE.BT, SB25$RMSE.BT, SB26$RMSE.BT, SB27$RMSE.BT,
        SB28$RMSE.BT, SB29$RMSE.BT, SB30$RMSE.BT, SB31$RMSE.BT, SB32$RMSE.BT, SB33$RMSE.BT,
        SB37$RMSE.BT, SB38$RMSE.BT, SB39$RMSE.BT,
        SB34$RMSE.BT, SB35$RMSE.BT, SB36$RMSE.BT, 
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
       ylim=c(0,2),cex=0.7,
       names=c('SB1','SB2','SB3','SB4','SB5','SB6','SB7','SB8','SB9','SB10','SB11','SB12','SB13','SB14','SB15',
               'SB16','SB17','SB18','SB19','SB20','SB21','SB22','SB23','SB24','SB25','SB26','SB27','SB28','SB29',
               'SB30','SB31','SB32','SB33','SB34','SB35','SB36','SB37','SB38','SB39') )
mtext('Back-transformed RMSE',side=2,line=1,cex=1.2)
# abline(h=mean(SB101$RMSE), col='white')

##############

dev.off()





png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Backtransformed_CompleteData_FULL_grid.png", 
    type="cairo",
    units="mm", 
    width=400, 
    height=170, 
    pointsize=10, 
    res=600)
#########
# STANDARDIZE AXES
par(mfrow=c(2,1), mar=c(0.6, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll.BT), 
        c(SB109$biasAll.BT), c(SB110$biasAll.BT), c(SB111$biasAll.BT), c(SB112$biasAll.BT), c(SB113$biasAll.BT), c(SB114$biasAll.BT), c(SB115$biasAll.BT),
        c(SB102$biasAll.BT), c(SB103$biasAll.BT), c(SB104$biasAll.BT), c(SB105$biasAll.BT), c(SB106$biasAll.BT), c(SB107$biasAll.BT), c(SB108$biasAll.BT),
        
        c(SB116$biasAll.BT), c(SB117$biasAll.BT), c(SB118$biasAll.BT), c(SB119$biasAll.BT), c(SB120$biasAll.BT), c(SB121$biasAll.BT),  
        c(SB122$biasAll.BT), c(SB123$biasAll.BT), c(SB124$biasAll.BT), c(SB125$biasAll.BT), c(SB126$biasAll.BT), c(SB127$biasAll.BT),
        c(SB128$biasAll.BT), c(SB129$biasAll.BT), c(SB130$biasAll.BT), c(SB131$biasAll.BT), c(SB132$biasAll.BT), c(SB133$biasAll.BT), 
        
        c(SB137$biasAll.BT), c(SB138$biasAll.BT), c(SB139$biasAll.BT),
        c(SB134$biasAll.BT), c(SB135$biasAll.BT), c(SB136$biasAll.BT), 
        
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=FALSE, cex=0.7 , ylim=c(-1.25,2),panel.first=grid(NA,NULL, lty=1) )
mtext('Back-transformed Annual Bias',side=2,line=1,cex=1.2)
abline(h=0, col='black')
mtext('Sandbar shark: Complete Data', line=-1.4, cex=1.2)

par( mar=c(1.3, 2.1, 0.1, 0.3))
vioplot(SB101$RMSE.BT, 
        SB109$RMSE.BT, SB110$RMSE.BT, SB111$RMSE.BT, SB112$RMSE.BT, SB113$RMSE.BT, SB114$RMSE.BT, SB115$RMSE.BT,
        SB102$RMSE.BT, SB103$RMSE.BT, SB104$RMSE.BT, SB105$RMSE.BT, SB106$RMSE.BT, SB107$RMSE.BT, SB108$RMSE.BT,
        SB116$RMSE.BT, SB117$RMSE.BT, SB118$RMSE.BT, SB119$RMSE.BT, SB120$RMSE.BT, SB121$RMSE.BT,
        SB122$RMSE.BT, SB123$RMSE.BT, SB124$RMSE.BT, SB125$RMSE.BT, SB126$RMSE.BT, SB127$RMSE.BT,
        SB128$RMSE.BT, SB129$RMSE.BT, SB130$RMSE.BT, SB131$RMSE.BT, SB132$RMSE.BT, SB133$RMSE.BT,
        
        SB137$RMSE.BT, SB138$RMSE.BT, SB139$RMSE.BT,
        SB134$RMSE.BT, SB135$RMSE.BT, SB136$RMSE.BT,
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        cex=0.7, ylim=c(0,0.5),panel.first=grid(NA,NULL, lty=1),
        names=c('SB101','SB102','SB103','SB104','SB105','SB106','SB107','SB108','SB109','SB110','SB111','SB112','SB113','SB114','SB115',
                'SB116','SB117','SB118','SB119','SB120','SB121','SB122','SB123','SB124','SB125','SB126','SB127','SB128','SB129',
                'SB130','SB131','SB132','SB133','SB134','SB135','SB136','SB137','SB138','SB139') )
mtext('Back-transformed RMSE',side=2,line=1,cex=1.2)
# abline(h=mean(SB101$RMSE.BT), col='white')

##############

dev.off()




png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Backtransformed_MissingData_FULL_grid.png", 
    type="cairo",
    units="mm", 
    width=400, 
    height=170, 
    pointsize=10, 
    res=600)
###########
# STANDARDIZE AXES
par(mfrow=c(2,1), mar=c(0.6, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB1$biasAll.BT), 
        c(SB9$biasAll.BT), c(SB10$biasAll.BT), c(SB11$biasAll.BT), c(SB12$biasAll.BT), c(SB13$biasAll.BT), c(SB14$biasAll.BT), c(SB15$biasAll.BT),
        c(SB2$biasAll.BT), c(SB3$biasAll.BT), c(SB4$biasAll.BT), c(SB5$biasAll.BT), c(SB6$biasAll.BT), c(SB7$biasAll.BT), c(SB8$biasAll.BT),
        c(SB16$biasAll.BT), c(SB17$biasAll.BT), c(SB18$biasAll.BT), c(SB19$biasAll.BT), c(SB20$biasAll.BT), c(SB21$biasAll.BT), 
        c(SB22$biasAll.BT), c(SB23$biasAll.BT), c(SB24$biasAll.BT), c(SB25$biasAll.BT), c(SB26$biasAll.BT), c(SB27$biasAll.BT),
        c(SB28$biasAll.BT), c(SB29$biasAll.BT), c(SB30$biasAll.BT), c(SB31$biasAll.BT), c(SB32$biasAll.BT), c(SB33$biasAll.BT), 
        c(SB37$biasAll.BT), c(SB38$biasAll.BT), c(SB39$biasAll.BT),
        c(SB34$biasAll.BT), c(SB35$biasAll.BT), c(SB36$biasAll.BT),
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=FALSE,  cex=0.7, ylim=c(-4,3.5), panel.first=grid(NA,NULL, lty=1))
mtext('Back-transformed Annual Bias',side=2,line=1,cex=1.2)
abline(h=0, col='black')
mtext('Sandbar shark: Missing Data', line=-1.4, cex=1.2)

par( mar=c(1.3, 2.1, 0.1, 0.3))
vioplot(SB1$RMSE.BT, 
        SB9$RMSE.BT,  SB10$RMSE.BT, SB11$RMSE.BT, SB12$RMSE.BT, SB13$RMSE.BT, SB14$RMSE.BT, SB15$RMSE.BT,
        SB2$RMSE.BT, SB3$RMSE.BT, SB4$RMSE.BT, SB5$RMSE.BT, SB6$RMSE.BT, SB7$RMSE.BT, SB8$RMSE.BT,
        SB16$RMSE.BT, SB17$RMSE.BT, SB18$RMSE.BT, SB19$RMSE.BT, SB20$RMSE.BT, SB21$RMSE.BT,
        SB22$RMSE.BT, SB23$RMSE.BT, SB24$RMSE.BT, SB25$RMSE.BT, SB26$RMSE.BT, SB27$RMSE.BT,
        SB28$RMSE.BT, SB29$RMSE.BT, SB30$RMSE.BT, SB31$RMSE.BT, SB32$RMSE.BT, SB33$RMSE.BT,
        SB37$RMSE.BT, SB38$RMSE.BT, SB39$RMSE.BT,
        SB34$RMSE.BT, SB35$RMSE.BT, SB36$RMSE.BT, 
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        ylim=c(0,2),cex=0.7, panel.first=grid(NA,NULL, lty=1), 
        names=c('SB1','SB2','SB3','SB4','SB5','SB6','SB7','SB8','SB9','SB10','SB11','SB12','SB13','SB14','SB15',
                'SB16','SB17','SB18','SB19','SB20','SB21','SB22','SB23','SB24','SB25','SB26','SB27','SB28','SB29',
                'SB30','SB31','SB32','SB33','SB34','SB35','SB36','SB37','SB38','SB39') )
mtext('Back-transformed RMSE',side=2,line=1,cex=1.2)
# abline(h=mean(SB101$RMSE), col='white')

##############

dev.off()




png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\FullVSMissingData_LgFont_2.png", 
    type="cairo",
    units="mm", 
    width=170, 
    height=150, 
    pointsize=8, 
    res=600)
#####
# STANDARDIZE AXES
par(mfrow=c(2,1), mar=c(0.6, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0), cex=1.0)
vioplot(c(SB101$biasAll), 
        c(SB1$biasAll), 
        
        c(SB109$biasAll, SB110$biasAll, SB111$biasAll, SB112$biasAll, SB113$biasAll, SB114$biasAll, SB115$biasAll, 
          SB102$biasAll, SB103$biasAll, SB104$biasAll, SB105$biasAll, SB106$biasAll, SB107$biasAll, SB108$biasAll),
        c(SB9$biasAll,  SB10$biasAll, SB11$biasAll, SB12$biasAll, SB13$biasAll, SB14$biasAll, SB15$biasAll, 
          SB2$biasAll, SB3$biasAll, SB4$biasAll, SB5$biasAll, SB6$biasAll, SB7$biasAll, SB8$biasAll),
        
        c(SB116$biasAll, SB117$biasAll, SB118$biasAll, SB119$biasAll, SB120$biasAll, SB121$biasAll, 
          SB122$biasAll, SB123$biasAll, SB124$biasAll, SB125$biasAll, SB126$biasAll, SB127$biasAll),
        c(SB16$biasAll, SB17$biasAll, SB18$biasAll, SB19$biasAll, SB20$biasAll, SB21$biasAll, 
          SB22$biasAll, SB23$biasAll, SB24$biasAll, SB25$biasAll, SB26$biasAll, SB27$biasAll),
        
        c(SB128$biasAll, SB129$biasAll, SB130$biasAll, SB131$biasAll, SB132$biasAll, SB133$biasAll, 
          SB134$biasAll, SB135$biasAll, SB136$biasAll, SB137$biasAll, SB138$biasAll, SB139$biasAll),
        c(SB28$biasAll, SB28$biasAll, SB30$biasAll, SB131$biasAll, SB32$biasAll, SB33$biasAll, 
          SB34$biasAll, SB35$biasAll, SB36$biasAll, SB37$biasAll, SB38$biasAll, SB39$biasAll),
        col=c('firebrick','firebrick1', 'deepskyblue4','deepskyblue','olivedrab','olivedrab1','darkorchid4','darkorchid'),
        names=FALSE) #, ylim=c(-4,4))
mtext('Annual Bias',side=2,line=1.1,cex=1.1)
title('Complete vs. Missing Data', line=-1, font.main=1, cex=0.9)
abline(h=0, col='black')
legend('bottomleft',c('Complete data','Missing data'),fill=c('darkgrey','lightgrey'), bty='n', cex=1.0, border=FALSE)


par(mar=c(1.6, 2.1, 0.1, 0.3))
vioplot(SB101$RMSE, 
        SB1$RMSE,
        
        c(SB109$RMSE, SB110$RMSE, SB111$RMSE, SB112$RMSE, SB113$RMSE, SB114$RMSE, SB115$RMSE, 
          SB102$RMSE, SB103$RMSE, SB104$RMSE, SB105$RMSE, SB106$RMSE, SB107$RMSE, SB108$RMSE),
        c(SB9$RMSE, SB10$RMSE, SB11$RMSE, SB12$RMSE, SB13$RMSE, SB14$RMSE, SB15$RMSE, 
          SB2$RMSE, SB3$RMSE, SB4$RMSE, SB5$RMSE, SB6$RMSE, SB7$RMSE, SB8$RMSE), 
        
        c(SB116$RMSE, SB117$RMSE, SB118$RMSE, SB119$RMSE, SB120$RMSE, SB121$RMSE, 
          SB122$RMSE, SB123$RMSE, SB124$RMSE, SB125$RMSE, SB126$RMSE, SB127$RMSE),
        c(SB16$RMSE, SB17$RMSE, SB18$RMSE, SB19$RMSE, SB20$RMSE, SB21$RMSE, 
          SB22$RMSE, SB23$RMSE, SB24$RMSE, SB25$RMSE, SB26$RMSE, SB27$RMSE),
        
        c(SB128$RMSE, SB129$RMSE, SB130$RMSE, SB131$RMSE, SB132$RMSE, SB133$RMSE, 
          SB134$RMSE, SB135$RMSE, SB136$RMSE, SB137$RMSE, SB138$RMSE, SB139$RMSE),
        c(SB28$RMSE, SB29$RMSE, SB30$RMSE, SB31$RMSE, SB32$RMSE, SB33$RMSE, 
          SB34$RMSE, SB35$RMSE, SB36$RMSE, SB37$RMSE, SB38$RMSE, SB39$RMSE),
        
        col=c('firebrick','firebrick1', 'deepskyblue4','deepskyblue','olivedrab','olivedrab1','darkorchid4','darkorchid'),
        # names=c('              No change','','              Change 1 q','','              Change 2 q','','               Change 3 q',''), 
        names=c('         No change','','      Change 1 q','','       Change 2 q','','     Change 3 q',''), 
        ylim=c(0,2))
mtext('RMSE',side=2,line=1.1,cex=1.1)
# legend('bottomright',c('Complete data','Missing data'),fill=c('darkgrey','lightgrey'), bty='n', cex=0.8)
# legend('bottomright',c('Complete data','Missing data'),fill=c('darkgrey','lightgrey'), bty='n', cex=0.8)
#####
dev.off()




png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\MS1\\Backtransformed+FullVSMissingData_LgFont_2.png", 
    type="cairo",
    units="mm", 
    width=170, 
    height=150, 
    pointsize=8, 
    res=600)
#####
# STANDARDIZE AXES
par(mfrow=c(2,1), mar=c(0.6, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0), cex=1.0)
vioplot(c(SB101$biasAll.BT), 
        c(SB1$biasAll.BT), 
        
        c(SB109$biasAll.BT, SB110$biasAll.BT, SB111$biasAll.BT, SB112$biasAll.BT, SB113$biasAll.BT, SB114$biasAll.BT, SB115$biasAll.BT, 
          SB102$biasAll.BT, SB103$biasAll.BT, SB104$biasAll.BT, SB105$biasAll.BT, SB106$biasAll.BT, SB107$biasAll.BT, SB108$biasAll.BT),
        c(SB9$biasAll.BT,  SB10$biasAll.BT, SB11$biasAll.BT, SB12$biasAll.BT, SB13$biasAll.BT, SB14$biasAll.BT, SB15$biasAll.BT, 
          SB2$biasAll.BT, SB3$biasAll.BT, SB4$biasAll.BT, SB5$biasAll.BT, SB6$biasAll.BT, SB7$biasAll.BT, SB8$biasAll.BT),
        
        c(SB116$biasAll.BT, SB117$biasAll.BT, SB118$biasAll.BT, SB119$biasAll.BT, SB120$biasAll.BT, SB121$biasAll.BT, 
          SB122$biasAll.BT, SB123$biasAll.BT, SB124$biasAll.BT, SB125$biasAll.BT, SB126$biasAll.BT, SB127$biasAll.BT),
        c(SB16$biasAll.BT, SB17$biasAll.BT, SB18$biasAll.BT, SB19$biasAll.BT, SB20$biasAll.BT, SB21$biasAll.BT, 
          SB22$biasAll.BT, SB23$biasAll.BT, SB24$biasAll.BT, SB25$biasAll.BT, SB26$biasAll.BT, SB27$biasAll.BT),
        
        c(SB128$biasAll.BT, SB129$biasAll.BT, SB130$biasAll.BT, SB131$biasAll.BT, SB132$biasAll.BT, SB133$biasAll.BT, 
          SB134$biasAll.BT, SB135$biasAll.BT, SB136$biasAll.BT, SB137$biasAll.BT, SB138$biasAll.BT, SB139$biasAll.BT),
        c(SB28$biasAll.BT, SB28$biasAll.BT, SB30$biasAll.BT, SB131$biasAll.BT, SB32$biasAll.BT, SB33$biasAll.BT, 
          SB34$biasAll.BT, SB35$biasAll.BT, SB36$biasAll.BT, SB37$biasAll.BT, SB38$biasAll.BT, SB39$biasAll.BT),
        col=c('firebrick','firebrick1', 'deepskyblue4','deepskyblue','olivedrab','olivedrab1','darkorchid4','darkorchid'),
        names=FALSE) #, ylim=c(-4,4))
mtext('Back-transformed Annual Bias',side=2,line=1.1,cex=1.1)
title('Complete vs. Missing Data', line=-1, font.main=1, cex=0.9)
abline(h=0, col='black')
legend('bottomleft',c('Complete data','Missing data'),fill=c('darkgrey','lightgrey'), bty='n', cex=1.0, border=FALSE)


par(mar=c(1.6, 2.1, 0.1, 0.3))
vioplot(SB101$RMSE.BT, 
        SB1$RMSE.BT,
        
        c(SB109$RMSE.BT, SB110$RMSE.BT, SB111$RMSE.BT, SB112$RMSE.BT, SB113$RMSE.BT, SB114$RMSE.BT, SB115$RMSE.BT, 
          SB102$RMSE.BT, SB103$RMSE.BT, SB104$RMSE.BT, SB105$RMSE.BT, SB106$RMSE.BT, SB107$RMSE.BT, SB108$RMSE.BT),
        c(SB9$RMSE.BT, SB10$RMSE.BT, SB11$RMSE.BT, SB12$RMSE.BT, SB13$RMSE.BT, SB14$RMSE.BT, SB15$RMSE.BT, 
          SB2$RMSE.BT, SB3$RMSE.BT, SB4$RMSE.BT, SB5$RMSE.BT, SB6$RMSE.BT, SB7$RMSE.BT, SB8$RMSE.BT), 
        
        c(SB116$RMSE.BT, SB117$RMSE.BT, SB118$RMSE.BT, SB119$RMSE.BT, SB120$RMSE.BT, SB121$RMSE.BT, 
          SB122$RMSE.BT, SB123$RMSE.BT, SB124$RMSE.BT, SB125$RMSE.BT, SB126$RMSE.BT, SB127$RMSE.BT),
        c(SB16$RMSE.BT, SB17$RMSE.BT, SB18$RMSE.BT, SB19$RMSE.BT, SB20$RMSE.BT, SB21$RMSE.BT, 
          SB22$RMSE.BT, SB23$RMSE.BT, SB24$RMSE.BT, SB25$RMSE.BT, SB26$RMSE.BT, SB27$RMSE.BT),
        
        c(SB128$RMSE.BT, SB129$RMSE.BT, SB130$RMSE.BT, SB131$RMSE.BT, SB132$RMSE.BT, SB133$RMSE.BT, 
          SB134$RMSE.BT, SB135$RMSE.BT, SB136$RMSE.BT, SB137$RMSE.BT, SB138$RMSE.BT, SB139$RMSE.BT),
        c(SB28$RMSE.BT, SB29$RMSE.BT, SB30$RMSE.BT, SB31$RMSE.BT, SB32$RMSE.BT, SB33$RMSE.BT, 
          SB34$RMSE.BT, SB35$RMSE.BT, SB36$RMSE.BT, SB37$RMSE.BT, SB38$RMSE.BT, SB39$RMSE.BT),
        
        col=c('firebrick','firebrick1', 'deepskyblue4','deepskyblue','olivedrab','olivedrab1','darkorchid4','darkorchid'),
        # names=c('              No change','','              Change 1 q','','              Change 2 q','','               Change 3 q',''), 
        names=c('         No change','','      Change 1 q','','       Change 2 q','','     Change 3 q',''), 
        ylim=c(0,2))
mtext('Back-transformed RMSE',side=2,line=1.1,cex=1.1)
# legend('bottomright',c('Complete data','Missing data'),fill=c('darkgrey','lightgrey'), bty='n', cex=0.8)
# legend('bottomright',c('Complete data','Missing data'),fill=c('darkgrey','lightgrey'), bty='n', cex=0.8)
#####
dev.off()





















png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Backtransformed_FullVSMissingData_3_3.png", 
    type="cairo",
    units="mm", 
    width=250, 
    height=250, 
    pointsize=16, 
    res=600)
#####
# STANDARDIZE AXES
par(mfrow=c(3,1), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll.BT), 
        c(SB1$biasAll.BT), 
        c(SB109$biasAll.BT, SB110$biasAll.BT, SB111$biasAll.BT, SB112$biasAll.BT, SB113$biasAll.BT, SB114$biasAll.BT, SB115$biasAll.BT),
        c(SB9$biasAll.BT,  SB10$biasAll.BT, SB11$biasAll.BT, SB12$biasAll.BT, SB13$biasAll.BT, SB14$biasAll.BT, SB15$biasAll.BT),
        c(SB102$biasAll.BT, SB103$biasAll.BT, SB104$biasAll.BT, SB105$biasAll.BT, SB106$biasAll.BT, SB107$biasAll.BT, SB108$biasAll.BT),
        c(SB2$biasAll.BT, SB3$biasAll.BT, SB4$biasAll.BT, SB5$biasAll.BT, SB6$biasAll.BT, SB7$biasAll.BT, SB8$biasAll.BT),
        c(SB116$biasAll.BT, SB117$biasAll.BT, SB118$biasAll.BT),
        c(SB16$biasAll.BT, SB17$biasAll.BT, SB18$biasAll.BT),
        c(SB119$biasAll.BT, SB120$biasAll.BT, SB121$biasAll.BT), 
        c(SB19$biasAll.BT, SB20$biasAll.BT, SB21$biasAll.BT), 
        c(SB122$biasAll.BT, SB123$biasAll.BT, SB124$biasAll.BT), 
        c(SB22$biasAll.BT, SB23$biasAll.BT, SB24$biasAll.BT), 
        c(SB125$biasAll.BT, SB126$biasAll.BT, SB127$biasAll.BT),
        c(SB25$biasAll.BT, SB26$biasAll.BT, SB27$biasAll.BT),
        
        c(SB128$biasAll.BT, SB129$biasAll.BT, SB130$biasAll.BT),
        c(SB28$biasAll.BT, SB28$biasAll.BT, SB30$biasAll.BT),
        c(SB131$biasAll.BT, SB132$biasAll.BT, SB133$biasAll.BT), 
        c(SB131$biasAll.BT, SB32$biasAll.BT, SB33$biasAll.BT), 
        c(SB134$biasAll.BT, SB135$biasAll.BT, SB136$biasAll.BT), 
        c(SB34$biasAll.BT, SB35$biasAll.BT, SB36$biasAll.BT), 
        c(SB137$biasAll.BT, SB138$biasAll.BT, SB139$biasAll.BT),
        c(SB37$biasAll.BT, SB38$biasAll.BT, SB39$biasAll.BT),
        col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
              'darkslategray4','darkslategray1','mediumorchid1','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
              'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple'),
        # col=c('lightgrey','darkgrey','olivedrab1','olivedrab3','skyblue1','deepskyblue','darkseagreen1','darkseagreen4',
        #       'darkslategray1','darkslategray4','orchid','mediumorchid1','darkorchid1','darkorchid4','springgreen','springgreen4',
        #       'steelblue1','steelblue4','mediumpurple1','mediumpurple4','purple','purple4'),
        names=FALSE) #, ylim=c(-4,3))
mtext('Back-transformed Annual Bias',side=2,line=1.1,cex=0.8)
title('Complete vs. Missing Data', line=-0.9)
abline(h=0, col='black')


vioplot(SB101$RMSE.BT, 
        SB1$RMSE.BT,
        c(SB109$RMSE.BT, SB110$RMSE.BT, SB111$RMSE.BT, SB112$RMSE.BT, SB113$RMSE.BT, SB114$RMSE.BT, SB115$RMSE.BT),
        c(SB9$RMSE.BT, SB10$RMSE.BT, SB11$RMSE.BT, SB12$RMSE.BT, SB13$RMSE.BT, SB14$RMSE.BT, SB15$RMSE.BT),
        c(SB102$RMSE.BT, SB103$RMSE.BT, SB104$RMSE.BT, SB105$RMSE.BT, SB106$RMSE.BT, SB107$RMSE.BT, SB108$RMSE.BT),
        c(SB2$RMSE.BT, SB3$RMSE.BT, SB4$RMSE.BT, SB5$RMSE.BT, SB6$RMSE.BT, SB7$RMSE.BT, SB8$RMSE.BT), 
        
        c(SB116$RMSE.BT, SB117$RMSE.BT, SB118$RMSE.BT),
        c(SB16$RMSE.BT, SB17$RMSE.BT, SB18$RMSE.BT),
        c(SB119$RMSE.BT, SB120$RMSE.BT, SB121$RMSE.BT), 
        c(SB19$RMSE.BT, SB20$RMSE.BT, SB21$RMSE.BT),
        c(SB122$RMSE.BT, SB123$RMSE.BT, SB124$RMSE.BT), 
        c(SB22$RMSE.BT, SB23$RMSE.BT, SB24$RMSE.BT),
        c(SB125$RMSE.BT, SB126$RMSE.BT, SB127$RMSE.BT),
        c(SB25$RMSE.BT, SB26$RMSE.BT, SB27$RMSE.BT),
        
        c(SB128$RMSE.BT, SB129$RMSE.BT, SB130$RMSE.BT),
        c(SB28$RMSE.BT, SB29$RMSE.BT, SB30$RMSE.BT),
        c(SB131$RMSE.BT, SB132$RMSE.BT, SB133$RMSE.BT), 
        c(SB31$RMSE.BT, SB32$RMSE.BT, SB33$RMSE.BT),
        c(SB134$RMSE.BT, SB135$RMSE.BT, SB136$RMSE.BT), 
        c(SB34$RMSE.BT, SB35$RMSE.BT, SB36$RMSE.BT),
        c(SB137$RMSE.BT, SB138$RMSE.BT, SB139$RMSE.BT),
        c(SB37$RMSE.BT, SB38$RMSE.BT, SB39$RMSE.BT),
        col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
              'darkslategray4','darkslategray1','mediumorchid1','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
              'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple'),
        # col=c('lightgrey','darkgrey','olivedrab1','olivedrab3','skyblue1','deepskyblue','darkseagreen1','darkseagreen4',
        #       'darkslategray1','darkslategray4','orchid','mediumorchid1','darkorchid1','darkorchid4','springgreen','springgreen4',
        #       'steelblue1','steelblue4','mediumpurple1','mediumpurple4','purple','purple4'),
        names=FALSE, ylim=c(0,2))
mtext('Back-transformed RMSE',side=2,line=1.1,cex=0.8)


vioplot(apply(SB101$FitRatios, 1, mean),
        apply(SB1$FitRatios, 1, mean),
        c(apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), apply(SB112$FitRatios, 1, mean), 
          apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean)), 
        c(apply(SB9$FitRatios, 1, mean), apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), apply(SB12$FitRatios, 1, mean), 
          apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean)), 
        c(apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), apply(SB105$FitRatios, 1, mean), 
          apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean)), 
        c(apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), apply(SB5$FitRatios, 1, mean), 
          apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean)), 
        
        c(apply(SB116$FitRatios, 1, mean), apply(SB117$FitRatios, 1, mean), apply(SB118$FitRatios, 1, mean)), 
        c(apply(SB16$FitRatios, 1, mean), apply(SB17$FitRatios, 1, mean), apply(SB18$FitRatios, 1, mean)), 
        c(apply(SB119$FitRatios, 1, mean), apply(SB120$FitRatios, 1, mean), apply(SB121$FitRatios, 1, mean)), 
        c(apply(SB19$FitRatios, 1, mean), apply(SB20$FitRatios, 1, mean), apply(SB21$FitRatios, 1, mean)), 
        c(apply(SB122$FitRatios, 1, mean), apply(SB123$FitRatios, 1, mean), apply(SB124$FitRatios, 1, mean)), 
        c(apply(SB22$FitRatios, 1, mean), apply(SB23$FitRatios, 1, mean), apply(SB24$FitRatios, 1, mean)), 
        c(apply(SB125$FitRatios, 1, mean), apply(SB126$FitRatios, 1, mean), apply(SB127$FitRatios, 1, mean)), 
        c(apply(SB25$FitRatios, 1, mean), apply(SB26$FitRatios, 1, mean), apply(SB27$FitRatios, 1, mean)), 
        
        c(apply(SB128$FitRatios, 1, mean), apply(SB129$FitRatios, 1, mean), apply(SB130$FitRatios, 1, mean)), 
        c(apply(SB28$FitRatios, 1, mean), apply(SB29$FitRatios, 1, mean), apply(SB30$FitRatios, 1, mean)), 
        c(apply(SB131$FitRatios, 1, mean), apply(SB132$FitRatios, 1, mean), apply(SB133$FitRatios, 1, mean)), 
        c(apply(SB31$FitRatios, 1, mean), apply(SB32$FitRatios, 1, mean), apply(SB33$FitRatios, 1, mean)), 
        c(apply(SB134$FitRatios, 1, mean), apply(SB135$FitRatios, 1, mean), apply(SB136$FitRatios, 1, mean)), 
        c(apply(SB34$FitRatios, 1, mean), apply(SB35$FitRatios, 1, mean), apply(SB36$FitRatios, 1, mean)), 
        c(apply(SB137$FitRatios, 1, mean), apply(SB138$FitRatios, 1, mean), apply(SB139$FitRatios, 1, mean)), 
        c(apply(SB37$FitRatios, 1, mean), apply(SB38$FitRatios, 1, mean), apply(SB39$FitRatios, 1, mean)), 
        
        col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
              'darkslategray4','darkslategray1','mediumorchid1','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
              'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple'),
        # col=c('lightgrey','darkgrey','olivedrab1','olivedrab3','skyblue1','deepskyblue','darkseagreen1','darkseagreen4',
        #       'darkslategray1','darkslategray4','orchid','mediumorchid1','darkorchid1','darkorchid4','springgreen','springgreen4',
        #       'steelblue1','steelblue4','mediumpurple1','mediumpurple4','purple','purple4'),
        names=c(
          'No
          change','','I1','','D1','','I2','','D2','','I-D','','D-I','','I3','','D3','','I-D-I','','D-I-D',''), 
        ylim=c(0,0.8))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
abline(h=0.6, col='grey', lty=2)
#####
dev.off()




################################


png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\CompleteData_ChangeQ3.png", 
    type="cairo",
    units="mm", 
    width=600, 
    height=250, 
    pointsize=16, 
    res=600)
#####
# STANDARDIZE AXES
par(mfrow=c(3,1), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll), 
        c(SB109$biasAll), c(SB110$biasAll), c(SB111$biasAll), c(SB112$biasAll), c(SB113$biasAll), c(SB114$biasAll), c(SB115$biasAll),
        c(SB102$biasAll), c(SB103$biasAll), c(SB104$biasAll), c(SB105$biasAll), c(SB106$biasAll), c(SB107$biasAll), c(SB108$biasAll),
        
        c(SB116$biasAll), c(SB117$biasAll), c(SB118$biasAll), c(SB119$biasAll), c(SB120$biasAll), c(SB121$biasAll),  
        c(SB122$biasAll), c(SB123$biasAll), c(SB124$biasAll), c(SB125$biasAll), c(SB126$biasAll), c(SB127$biasAll),
        c(SB128$biasAll), c(SB129$biasAll), c(SB130$biasAll), c(SB131$biasAll), c(SB132$biasAll), c(SB133$biasAll), 
        c(SB134$biasAll), c(SB135$biasAll), c(SB136$biasAll), c(SB137$biasAll), c(SB138$biasAll), c(SB139$biasAll),
        
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=FALSE, ylim=c(-1,1.5) )
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
abline(h=0, col='black')
title('Complete Data', line=-0.9)

vioplot(SB101$RMSE, 
        SB109$RMSE, SB110$RMSE, SB111$RMSE, SB112$RMSE, SB113$RMSE, SB114$RMSE, SB115$RMSE,
        SB102$RMSE, SB103$RMSE, SB104$RMSE, SB105$RMSE, SB106$RMSE, SB107$RMSE, SB108$RMSE,
        SB116$RMSE, SB117$RMSE, SB118$RMSE, SB119$RMSE, SB120$RMSE, SB121$RMSE,
        SB122$RMSE, SB123$RMSE, SB124$RMSE, SB125$RMSE, SB126$RMSE, SB127$RMSE,
        SB128$RMSE, SB129$RMSE, SB130$RMSE, SB131$RMSE, SB132$RMSE, SB133$RMSE,
        SB134$RMSE, SB135$RMSE, SB136$RMSE, SB137$RMSE, SB138$RMSE, SB139$RMSE,
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=FALSE, ylim=c(0,0.5))
mtext('RMSE',side=2,line=1.1,cex=0.8)
# abline(h=mean(SB101$RMSE), col='white')


vioplot(apply(SB101$FitRatios, 1, mean), 
        apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), 
        apply(SB112$FitRatios, 1, mean), apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean),
        apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), 
        apply(SB105$FitRatios, 1, mean), apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean), 
        
        apply(SB116$FitRatios, 1, mean), apply(SB117$FitRatios, 1, mean), apply(SB118$FitRatios, 1, mean), 
        apply(SB119$FitRatios, 1, mean), apply(SB120$FitRatios, 1, mean), apply(SB121$FitRatios, 1, mean), 
        
        apply(SB122$FitRatios, 1, mean), apply(SB123$FitRatios, 1, mean), apply(SB124$FitRatios, 1, mean), 
        apply(SB125$FitRatios, 1, mean), apply(SB126$FitRatios, 1, mean), apply(SB127$FitRatios, 1, mean),  
        
        apply(SB128$FitRatios, 1, mean), apply(SB129$FitRatios, 1, mean), apply(SB130$FitRatios, 1, mean), 
        apply(SB131$FitRatios, 1, mean), apply(SB132$FitRatios, 1, mean), apply(SB133$FitRatios, 1, mean), 
        
        apply(SB134$FitRatios, 1, mean), apply(SB135$FitRatios, 1, mean), apply(SB136$FitRatios, 1, mean), 
        apply(SB137$FitRatios, 1, mean), apply(SB138$FitRatios, 1, mean), apply(SB139$FitRatios, 1, mean), 
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=c('None','I-S1','I-S2','I-S3','I-S4','I-S5','I-S6','I-S7','D-S1','D-S2','D-S3','D-S4','D-S5','D-S6','D-S7',
                'I2','I2','I2','D2','D2','D2','I1-D1','I1-D1','I1-D1','D1-I1','D1-I1','D1-I1', 
                'I3','I3','I3','D3','D3','D3','I-D-I','I-D-I','I-D-I','D-I-D','D-I-D','D-I-D'), ylim=c(0,0.6) )
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
# abline(h=mean(apply(SB101$FitRatios, 1, mean)), col='white')
#####

dev.off()



boxplot(SB101$RMSE, 
        SB109$RMSE, SB110$RMSE, SB111$RMSE, SB112$RMSE, SB113$RMSE, SB114$RMSE, SB115$RMSE,
        SB102$RMSE, SB103$RMSE, SB104$RMSE, SB105$RMSE, SB106$RMSE, SB107$RMSE, SB108$RMSE,
        SB116$RMSE, SB117$RMSE, SB118$RMSE, SB119$RMSE, SB120$RMSE, SB121$RMSE,
        SB122$RMSE, SB123$RMSE, SB124$RMSE, SB125$RMSE, SB126$RMSE, SB127$RMSE,
        SB128$RMSE, SB129$RMSE, SB130$RMSE, SB131$RMSE, SB132$RMSE, SB133$RMSE,
        SB134$RMSE, SB135$RMSE, SB136$RMSE, SB137$RMSE, SB138$RMSE, SB139$RMSE,
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=FALSE, ylim=c(0,0.60))

boxplot(apply(SB101$FitRatios, 1, mean), 
        apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), 
        apply(SB112$FitRatios, 1, mean), apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean),
        apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), 
        apply(SB105$FitRatios, 1, mean), apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean), 
        
        apply(SB116$FitRatios, 1, mean), apply(SB117$FitRatios, 1, mean), apply(SB118$FitRatios, 1, mean), 
        apply(SB119$FitRatios, 1, mean), apply(SB120$FitRatios, 1, mean), apply(SB121$FitRatios, 1, mean), 
        
        apply(SB122$FitRatios, 1, mean), apply(SB123$FitRatios, 1, mean), apply(SB124$FitRatios, 1, mean), 
        apply(SB125$FitRatios, 1, mean), apply(SB126$FitRatios, 1, mean), apply(SB127$FitRatios, 1, mean),  
        
        apply(SB128$FitRatios, 1, mean), apply(SB129$FitRatios, 1, mean), apply(SB130$FitRatios, 1, mean), 
        apply(SB131$FitRatios, 1, mean), apply(SB132$FitRatios, 1, mean), apply(SB133$FitRatios, 1, mean), 
        
        apply(SB134$FitRatios, 1, mean), apply(SB135$FitRatios, 1, mean), apply(SB136$FitRatios, 1, mean), 
        apply(SB137$FitRatios, 1, mean), apply(SB138$FitRatios, 1, mean), apply(SB139$FitRatios, 1, mean), 
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=c('None','I-S1','I-S2','I-S3','I-S4','I-S5','I-S6','I-S7','D-S1','D-S2','D-S3','D-S4','D-S5','D-S6','D-S7',
                'I2','I2','I2','D2','D2','D2','I1-D1','I1-D1','I1-D1','D1-I1','D1-I1','D1-I1', 
                'I3','I3','I3','D3','D3','D3','I-D-I','I-D-I','I-D-I','D-I-D','D-I-D','D-I-D'), ylim=c(0,0.7))








png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\MissingData_ChangeQ3.png", 
    type="cairo",
    units="mm", 
    width=600, 
    height=250, 
    pointsize=16, 
    res=600)
#####
# STANDARDIZE AXES
par(mfrow=c(3,1), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB1$biasAll), 
        c(SB9$biasAll), c(SB10$biasAll), c(SB11$biasAll), c(SB12$biasAll), c(SB13$biasAll), c(SB14$biasAll), c(SB15$biasAll),
        c(SB2$biasAll), c(SB3$biasAll), c(SB4$biasAll), c(SB5$biasAll), c(SB6$biasAll), c(SB7$biasAll), c(SB8$biasAll),
        c(SB16$biasAll), c(SB17$biasAll), c(SB18$biasAll), c(SB19$biasAll), c(SB20$biasAll), c(SB21$biasAll), 
        c(SB22$biasAll), c(SB23$biasAll), c(SB24$biasAll), c(SB25$biasAll), c(SB26$biasAll), c(SB27$biasAll),
        c(SB28$biasAll), c(SB29$biasAll), c(SB30$biasAll), c(SB31$biasAll), c(SB32$biasAll), c(SB33$biasAll), 
        c(SB34$biasAll), c(SB35$biasAll), c(SB36$biasAll), c(SB37$biasAll), c(SB38$biasAll), c(SB39$biasAll),
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=FALSE) #, ylim=c(-3.5,4) )
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
abline(h=0, col='black')
title('Missing Data', line=-0.9)

vioplot(SB1$RMSE, 
        SB9$RMSE,  SB10$RMSE, SB11$RMSE, SB12$RMSE, SB13$RMSE, SB14$RMSE, SB15$RMSE,
        SB2$RMSE, SB3$RMSE, SB4$RMSE, SB5$RMSE, SB6$RMSE, SB7$RMSE, SB8$RMSE,
        SB16$RMSE, SB17$RMSE, SB18$RMSE, SB19$RMSE, SB20$RMSE, SB21$RMSE,
        SB22$RMSE, SB23$RMSE, SB24$RMSE, SB25$RMSE, SB26$RMSE, SB27$RMSE,
        SB28$RMSE, SB29$RMSE, SB30$RMSE, SB31$RMSE, SB32$RMSE, SB33$RMSE,
        SB34$RMSE, SB35$RMSE, SB36$RMSE, SB37$RMSE, SB38$RMSE, SB39$RMSE,
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=FALSE, ylim=c(0,2))
mtext('RMSE',side=2,line=1.1,cex=0.8)
# abline(h=mean(SB101$RMSE), col='white')


vioplot(apply(SB1$FitRatios, 1, mean), 
        apply(SB9$FitRatios, 1, mean),  apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), 
        apply(SB12$FitRatios, 1, mean), apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean),
        apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), 
        apply(SB5$FitRatios, 1, mean), apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean), 
        
        apply(SB16$FitRatios, 1, mean), apply(SB17$FitRatios, 1, mean), apply(SB18$FitRatios, 1, mean), 
        apply(SB19$FitRatios, 1, mean), apply(SB20$FitRatios, 1, mean), apply(SB21$FitRatios, 1, mean), 
        apply(SB22$FitRatios, 1, mean), apply(SB23$FitRatios, 1, mean), apply(SB24$FitRatios, 1, mean), 
        apply(SB25$FitRatios, 1, mean), apply(SB26$FitRatios, 1, mean), apply(SB27$FitRatios, 1, mean), 
        
        apply(SB28$FitRatios, 1, mean), apply(SB29$FitRatios, 1, mean), apply(SB30$FitRatios, 1, mean), 
        apply(SB31$FitRatios, 1, mean), apply(SB32$FitRatios, 1, mean), apply(SB33$FitRatios, 1, mean), 
        apply(SB34$FitRatios, 1, mean), apply(SB35$FitRatios, 1, mean), apply(SB36$FitRatios, 1, mean), 
        apply(SB37$FitRatios, 1, mean), apply(SB38$FitRatios, 1, mean), apply(SB39$FitRatios, 1, mean), 
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=c('None','I-S1','I-S2','I-S3','I-S4','I-S5','I-S6','I-S7','D-S1','D-S2','D-S3','D-S4','D-S5','D-S6','D-S7',
                'I2','I2','I2','D2','D2','D2','I1-D1','I1-D1','I1-D1','D1-I1','D1-I1','D1-I1', 
                'I3','I3','I3','D3','D3','D3','I-D-I','I-D-I','I-D-I','D-I-D','D-I-D','D-I-D'), ylim=c(0,1) )
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
# abline(h=mean(apply(SB1$FitRatios, 1, mean)), col='white')
#####

dev.off()







png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\FullVSMissingData_3_2.png", 
    type="cairo",
    units="mm", 
    width=250, 
    height=250, 
    pointsize=16, 
    res=600)
#####
# STANDARDIZE AXES
par(mfrow=c(3,2), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll), 
        c(SB109$biasAll, SB110$biasAll, SB111$biasAll, SB112$biasAll, SB113$biasAll, SB114$biasAll, SB115$biasAll),
        c(SB102$biasAll, SB103$biasAll, SB104$biasAll, SB105$biasAll, SB106$biasAll, SB107$biasAll, SB108$biasAll),
        c(SB116$biasAll, SB117$biasAll, SB118$biasAll), c(SB119$biasAll, SB120$biasAll, SB121$biasAll), 
        c(SB122$biasAll, SB123$biasAll, SB124$biasAll), c(SB125$biasAll, SB126$biasAll, SB127$biasAll),
        c(SB128$biasAll, SB129$biasAll, SB130$biasAll), c(SB131$biasAll, SB132$biasAll, SB133$biasAll), 
        c(SB134$biasAll, SB135$biasAll, SB136$biasAll), c(SB137$biasAll, SB138$biasAll, SB139$biasAll),
        col=c('grey', 'darkolivegreen1', 'deepskyblue', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3',
              'springgreen3', 'steelblue3','mediumpurple3','purple'),
        names=FALSE, ylim=c(-4,3.5))
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
abline(h=0, col='black')
title('Complete Data', line=-0.9)
vioplot(c(SB1$biasAll), 
        c(SB9$biasAll,  SB10$biasAll, SB11$biasAll, SB12$biasAll, SB13$biasAll, SB14$biasAll, SB15$biasAll),
        c(SB2$biasAll, SB3$biasAll, SB4$biasAll, SB5$biasAll, SB6$biasAll, SB7$biasAll, SB8$biasAll),
        c(SB16$biasAll, SB17$biasAll, SB18$biasAll), c(SB19$biasAll, SB20$biasAll, SB21$biasAll),
        c(SB22$biasAll, SB23$biasAll, SB24$biasAll), c(SB25$biasAll, SB26$biasAll, SB27$biasAll),
        c(SB28$biasAll, SB29$biasAll, SB30$biasAll), c(SB31$biasAll, SB32$biasAll, SB33$biasAll),
        c(SB34$biasAll, SB35$biasAll, SB36$biasAll), c(SB37$biasAll, SB38$biasAll, SB39$biasAll),
        col=c('grey', 'darkolivegreen1', 'deepskyblue', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3',
              'springgreen3', 'steelblue3','mediumpurple3','purple'),
        names=FALSE, ylim=c(-4,3.5))
abline(h=0, col='black')
title('Missing Data', line=-0.9)

vioplot(SB101$RMSE, 
        c(SB109$RMSE, SB110$RMSE, SB111$RMSE, SB112$RMSE, SB113$RMSE, SB114$RMSE, SB115$RMSE),
        c(SB102$RMSE, SB103$RMSE, SB104$RMSE, SB105$RMSE, SB106$RMSE, SB107$RMSE, SB108$RMSE),
        c(SB116$RMSE, SB117$RMSE, SB118$RMSE), c(SB119$RMSE, SB120$RMSE, SB121$RMSE), 
        c(SB122$RMSE, SB123$RMSE, SB124$RMSE), c(SB125$RMSE, SB126$RMSE, SB127$RMSE),
        c(SB128$RMSE, SB129$RMSE, SB130$RMSE), c(SB131$RMSE, SB132$RMSE, SB133$RMSE), 
        c(SB134$RMSE, SB135$RMSE, SB136$RMSE), c(SB137$RMSE, SB138$RMSE, SB139$RMSE),
        col=c('grey', 'darkolivegreen1', 'deepskyblue', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3',
              'springgreen3', 'steelblue3','mediumpurple3','purple'),
        names=FALSE, ylim=c(0,2))
mtext('RMSE',side=2,line=1.1,cex=0.8)
vioplot(SB1$RMSE, 
        c(SB9$RMSE, SB10$RMSE, SB11$RMSE, SB12$RMSE, SB13$RMSE, SB14$RMSE, SB15$RMSE),
        c(SB2$RMSE, SB3$RMSE, SB4$RMSE, SB5$RMSE, SB6$RMSE, SB7$RMSE, SB8$RMSE),
        c(SB16$RMSE, SB17$RMSE, SB18$RMSE), c(SB19$RMSE, SB20$RMSE, SB21$RMSE), 
        c(SB22$RMSE, SB23$RMSE, SB24$RMSE), c(SB25$RMSE, SB26$RMSE, SB27$RMSE),
        c(SB28$RMSE, SB29$RMSE, SB30$RMSE), c(SB31$RMSE, SB32$RMSE, SB33$RMSE), 
        c(SB34$RMSE, SB35$RMSE, SB36$RMSE), c(SB37$RMSE, SB38$RMSE, SB39$RMSE),
        col=c('grey', 'darkolivegreen1', 'deepskyblue', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3',
              'springgreen3', 'steelblue3','mediumpurple3','purple'),
        names=FALSE, ylim=c(0,2))


vioplot(apply(SB101$FitRatios, 1, mean),
        c(apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), apply(SB112$FitRatios, 1, mean), 
          apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean)), 
        c(apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), apply(SB105$FitRatios, 1, mean), 
          apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean)), 
        
        c(apply(SB116$FitRatios, 1, mean), apply(SB117$FitRatios, 1, mean), apply(SB118$FitRatios, 1, mean)), 
        c(apply(SB119$FitRatios, 1, mean), apply(SB120$FitRatios, 1, mean), apply(SB121$FitRatios, 1, mean)), 
        
        c(apply(SB122$FitRatios, 1, mean), apply(SB123$FitRatios, 1, mean), apply(SB124$FitRatios, 1, mean)), 
        c(apply(SB125$FitRatios, 1, mean), apply(SB126$FitRatios, 1, mean), apply(SB127$FitRatios, 1, mean)), 
        
        c(apply(SB128$FitRatios, 1, mean), apply(SB129$FitRatios, 1, mean), apply(SB130$FitRatios, 1, mean)), 
        c(apply(SB131$FitRatios, 1, mean), apply(SB132$FitRatios, 1, mean), apply(SB133$FitRatios, 1, mean)), 
        
        c(apply(SB134$FitRatios, 1, mean), apply(SB135$FitRatios, 1, mean), apply(SB136$FitRatios, 1, mean)), 
        c(apply(SB137$FitRatios, 1, mean), apply(SB138$FitRatios, 1, mean), apply(SB139$FitRatios, 1, mean)),
        col=c('grey', 'darkolivegreen1', 'deepskyblue', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3',
              'springgreen3', 'steelblue3','mediumpurple3','purple'),
        names=c(
'No
change','I1','D1','I2','D2','I-D','D-I','I3','D3','I-D-I','D-I-D'), ylim=c(0,1))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)

vioplot(apply(SB1$FitRatios, 1, mean),
        c(apply(SB9$FitRatios, 1, mean), apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), apply(SB12$FitRatios, 1, mean), 
          apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean)), 
        c(apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), apply(SB5$FitRatios, 1, mean), 
          apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean)), 
        
        c(apply(SB16$FitRatios, 1, mean), apply(SB17$FitRatios, 1, mean), apply(SB18$FitRatios, 1, mean)), 
        c(apply(SB19$FitRatios, 1, mean), apply(SB20$FitRatios, 1, mean), apply(SB21$FitRatios, 1, mean)), 
        
        c(apply(SB22$FitRatios, 1, mean), apply(SB23$FitRatios, 1, mean), apply(SB24$FitRatios, 1, mean)), 
        c(apply(SB25$FitRatios, 1, mean), apply(SB26$FitRatios, 1, mean), apply(SB27$FitRatios, 1, mean)), 
        
        c(apply(SB28$FitRatios, 1, mean), apply(SB29$FitRatios, 1, mean), apply(SB30$FitRatios, 1, mean)), 
        c(apply(SB31$FitRatios, 1, mean), apply(SB32$FitRatios, 1, mean), apply(SB33$FitRatios, 1, mean)), 
        
        c(apply(SB34$FitRatios, 1, mean), apply(SB35$FitRatios, 1, mean), apply(SB36$FitRatios, 1, mean)), 
        c(apply(SB37$FitRatios, 1, mean), apply(SB38$FitRatios, 1, mean), apply(SB39$FitRatios, 1, mean)),
        col=c('grey', 'darkolivegreen1', 'deepskyblue', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3',
              'springgreen3', 'steelblue3','mediumpurple3','purple'),
        names=c(
'No
change','I1','D1','I2','D2','I-D','D-I','I3','D3','I-D-I','D-I-D'), ylim=c(0,1))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
#####

dev.off()







png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\FullVSMissingData_3_3.png", 
    type="cairo",
    units="mm", 
    width=250, 
    height=250, 
    pointsize=16, 
    res=600)
#####
# STANDARDIZE AXES
par(mfrow=c(3,1), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll), 
        c(SB1$biasAll), 
        c(SB109$biasAll, SB110$biasAll, SB111$biasAll, SB112$biasAll, SB113$biasAll, SB114$biasAll, SB115$biasAll),
        c(SB9$biasAll,  SB10$biasAll, SB11$biasAll, SB12$biasAll, SB13$biasAll, SB14$biasAll, SB15$biasAll),
        c(SB102$biasAll, SB103$biasAll, SB104$biasAll, SB105$biasAll, SB106$biasAll, SB107$biasAll, SB108$biasAll),
        c(SB2$biasAll, SB3$biasAll, SB4$biasAll, SB5$biasAll, SB6$biasAll, SB7$biasAll, SB8$biasAll),
        c(SB116$biasAll, SB117$biasAll, SB118$biasAll),
        c(SB16$biasAll, SB17$biasAll, SB18$biasAll),
        c(SB119$biasAll, SB120$biasAll, SB121$biasAll), 
        c(SB19$biasAll, SB20$biasAll, SB21$biasAll), 
        c(SB122$biasAll, SB123$biasAll, SB124$biasAll), 
        c(SB22$biasAll, SB23$biasAll, SB24$biasAll), 
        c(SB125$biasAll, SB126$biasAll, SB127$biasAll),
        c(SB25$biasAll, SB26$biasAll, SB27$biasAll),
        
        c(SB128$biasAll, SB129$biasAll, SB130$biasAll),
        c(SB28$biasAll, SB28$biasAll, SB30$biasAll),
        c(SB131$biasAll, SB132$biasAll, SB133$biasAll), 
        c(SB131$biasAll, SB32$biasAll, SB33$biasAll), 
        c(SB134$biasAll, SB135$biasAll, SB136$biasAll), 
        c(SB34$biasAll, SB35$biasAll, SB36$biasAll), 
        c(SB137$biasAll, SB138$biasAll, SB139$biasAll),
        c(SB37$biasAll, SB38$biasAll, SB39$biasAll),
        col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
              'darkslategray4','darkslategray1','mediumorchid1','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
              'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple'),
        # col=c('lightgrey','darkgrey','olivedrab1','olivedrab3','skyblue1','deepskyblue','darkseagreen1','darkseagreen4',
        #       'darkslategray1','darkslategray4','orchid','mediumorchid1','darkorchid1','darkorchid4','springgreen','springgreen4',
        #       'steelblue1','steelblue4','mediumpurple1','mediumpurple4','purple','purple4'),
        names=FALSE, ylim=c(-4,5))
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
title('Complete vs. Missing Data', line=-0.9)
abline(h=0, col='black')


vioplot(SB101$RMSE, 
        SB1$RMSE,
        c(SB109$RMSE, SB110$RMSE, SB111$RMSE, SB112$RMSE, SB113$RMSE, SB114$RMSE, SB115$RMSE),
        c(SB9$RMSE, SB10$RMSE, SB11$RMSE, SB12$RMSE, SB13$RMSE, SB14$RMSE, SB15$RMSE),
        c(SB102$RMSE, SB103$RMSE, SB104$RMSE, SB105$RMSE, SB106$RMSE, SB107$RMSE, SB108$RMSE),
        c(SB2$RMSE, SB3$RMSE, SB4$RMSE, SB5$RMSE, SB6$RMSE, SB7$RMSE, SB8$RMSE), 
        
        c(SB116$RMSE, SB117$RMSE, SB118$RMSE),
        c(SB16$RMSE, SB17$RMSE, SB18$RMSE),
        c(SB119$RMSE, SB120$RMSE, SB121$RMSE), 
        c(SB19$RMSE, SB20$RMSE, SB21$RMSE),
        c(SB122$RMSE, SB123$RMSE, SB124$RMSE), 
        c(SB22$RMSE, SB23$RMSE, SB24$RMSE),
        c(SB125$RMSE, SB126$RMSE, SB127$RMSE),
        c(SB25$RMSE, SB26$RMSE, SB27$RMSE),
        
        c(SB128$RMSE, SB129$RMSE, SB130$RMSE),
        c(SB28$RMSE, SB29$RMSE, SB30$RMSE),
        c(SB131$RMSE, SB132$RMSE, SB133$RMSE), 
        c(SB31$RMSE, SB32$RMSE, SB33$RMSE),
        c(SB134$RMSE, SB135$RMSE, SB136$RMSE), 
        c(SB34$RMSE, SB35$RMSE, SB36$RMSE),
        c(SB137$RMSE, SB138$RMSE, SB139$RMSE),
        c(SB37$RMSE, SB38$RMSE, SB39$RMSE),
        col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
              'darkslategray4','darkslategray1','mediumorchid1','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
              'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple'),
        # col=c('lightgrey','darkgrey','olivedrab1','olivedrab3','skyblue1','deepskyblue','darkseagreen1','darkseagreen4',
        #       'darkslategray1','darkslategray4','orchid','mediumorchid1','darkorchid1','darkorchid4','springgreen','springgreen4',
        #       'steelblue1','steelblue4','mediumpurple1','mediumpurple4','purple','purple4'),
        names=FALSE, ylim=c(0,2))
mtext('RMSE',side=2,line=1.1,cex=0.8)


vioplot(apply(SB101$FitRatios, 1, mean),
        apply(SB1$FitRatios, 1, mean),
        c(apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), apply(SB112$FitRatios, 1, mean), 
          apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean)), 
        c(apply(SB9$FitRatios, 1, mean), apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), apply(SB12$FitRatios, 1, mean), 
          apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean)), 
        c(apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), apply(SB105$FitRatios, 1, mean), 
          apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean)), 
        c(apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), apply(SB5$FitRatios, 1, mean), 
          apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean)), 
        
        c(apply(SB116$FitRatios, 1, mean), apply(SB117$FitRatios, 1, mean), apply(SB118$FitRatios, 1, mean)), 
        c(apply(SB16$FitRatios, 1, mean), apply(SB17$FitRatios, 1, mean), apply(SB18$FitRatios, 1, mean)), 
        c(apply(SB119$FitRatios, 1, mean), apply(SB120$FitRatios, 1, mean), apply(SB121$FitRatios, 1, mean)), 
        c(apply(SB19$FitRatios, 1, mean), apply(SB20$FitRatios, 1, mean), apply(SB21$FitRatios, 1, mean)), 
        c(apply(SB122$FitRatios, 1, mean), apply(SB123$FitRatios, 1, mean), apply(SB124$FitRatios, 1, mean)), 
        c(apply(SB22$FitRatios, 1, mean), apply(SB23$FitRatios, 1, mean), apply(SB24$FitRatios, 1, mean)), 
        c(apply(SB125$FitRatios, 1, mean), apply(SB126$FitRatios, 1, mean), apply(SB127$FitRatios, 1, mean)), 
        c(apply(SB25$FitRatios, 1, mean), apply(SB26$FitRatios, 1, mean), apply(SB27$FitRatios, 1, mean)), 
        
        c(apply(SB128$FitRatios, 1, mean), apply(SB129$FitRatios, 1, mean), apply(SB130$FitRatios, 1, mean)), 
        c(apply(SB28$FitRatios, 1, mean), apply(SB29$FitRatios, 1, mean), apply(SB30$FitRatios, 1, mean)), 
        c(apply(SB131$FitRatios, 1, mean), apply(SB132$FitRatios, 1, mean), apply(SB133$FitRatios, 1, mean)), 
        c(apply(SB31$FitRatios, 1, mean), apply(SB32$FitRatios, 1, mean), apply(SB33$FitRatios, 1, mean)), 
        c(apply(SB134$FitRatios, 1, mean), apply(SB135$FitRatios, 1, mean), apply(SB136$FitRatios, 1, mean)), 
        c(apply(SB34$FitRatios, 1, mean), apply(SB35$FitRatios, 1, mean), apply(SB36$FitRatios, 1, mean)), 
        c(apply(SB137$FitRatios, 1, mean), apply(SB138$FitRatios, 1, mean), apply(SB139$FitRatios, 1, mean)), 
        c(apply(SB37$FitRatios, 1, mean), apply(SB38$FitRatios, 1, mean), apply(SB39$FitRatios, 1, mean)), 
        
        col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
              'darkslategray4','darkslategray1','mediumorchid1','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
              'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple'),
        # col=c('lightgrey','darkgrey','olivedrab1','olivedrab3','skyblue1','deepskyblue','darkseagreen1','darkseagreen4',
        #       'darkslategray1','darkslategray4','orchid','mediumorchid1','darkorchid1','darkorchid4','springgreen','springgreen4',
        #       'steelblue1','steelblue4','mediumpurple1','mediumpurple4','purple','purple4'),
        names=c(
'No
change','','I1','','D1','','I2','','D2','','I-D','','D-I','','I3','','D3','','I-D-I','','D-I-D',''), ylim=c(0,1))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
abline(h=0.6, col='grey', lty=2)
#####
dev.off()










# Plot backtransformed ---------------------------------------------------

png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Backtransformed_CompleteData_ChangeQ3.png", 
    type="cairo",
    units="mm", 
    width=600, 
    height=250, 
    pointsize=16, 
    res=600)
#####
# STANDARDIZE AXES
par(mfrow=c(3,1), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll.BT), 
        c(SB109$biasAll.BT), c(SB110$biasAll.BT), c(SB111$biasAll.BT), c(SB112$biasAll.BT), c(SB113$biasAll.BT), c(SB114$biasAll.BT), c(SB115$biasAll.BT),
        c(SB102$biasAll.BT), c(SB103$biasAll.BT), c(SB104$biasAll.BT), c(SB105$biasAll.BT), c(SB106$biasAll.BT), c(SB107$biasAll.BT), c(SB108$biasAll.BT),
        
        c(SB116$biasAll.BT), c(SB117$biasAll.BT), c(SB118$biasAll.BT), c(SB119$biasAll.BT), c(SB120$biasAll.BT), c(SB121$biasAll.BT),  
        c(SB122$biasAll.BT), c(SB123$biasAll.BT), c(SB124$biasAll.BT), c(SB125$biasAll.BT), c(SB126$biasAll.BT), c(SB127$biasAll.BT),
        c(SB128$biasAll.BT), c(SB129$biasAll.BT), c(SB130$biasAll.BT), c(SB131$biasAll.BT), c(SB132$biasAll.BT), c(SB133$biasAll.BT), 
        c(SB134$biasAll.BT), c(SB135$biasAll.BT), c(SB136$biasAll.BT), c(SB137$biasAll.BT), c(SB138$biasAll.BT), c(SB139$biasAll.BT),
        
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=FALSE) #, ylim=c(-1,1) )
mtext('Back-transformed Annual Bias',side=2,line=1.1,cex=0.8)
abline(h=0, col='black')
title('Complete Data', line=-0.9)

vioplot(SB101$RMSE.BT, 
        SB109$RMSE.BT, SB110$RMSE.BT, SB111$RMSE.BT, SB112$RMSE.BT, SB113$RMSE.BT, SB114$RMSE.BT, SB115$RMSE.BT,
        SB102$RMSE.BT, SB103$RMSE.BT, SB104$RMSE.BT, SB105$RMSE.BT, SB106$RMSE.BT, SB107$RMSE.BT, SB108$RMSE.BT,
        SB116$RMSE.BT, SB117$RMSE.BT, SB118$RMSE.BT, SB119$RMSE.BT, SB120$RMSE.BT, SB121$RMSE.BT,
        SB122$RMSE.BT, SB123$RMSE.BT, SB124$RMSE.BT, SB125$RMSE.BT, SB126$RMSE.BT, SB127$RMSE.BT,
        SB128$RMSE.BT, SB129$RMSE.BT, SB130$RMSE.BT, SB131$RMSE.BT, SB132$RMSE.BT, SB133$RMSE.BT,
        SB134$RMSE.BT, SB135$RMSE.BT, SB136$RMSE.BT, SB137$RMSE.BT, SB138$RMSE.BT, SB139$RMSE.BT,
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=FALSE, ylim=c(0,2))
mtext('Back-transformed RMSE',side=2,line=1.1,cex=0.8)
# abline(h=mean(SB101$RMSE.BT), col='white')


vioplot(apply(SB101$FitRatios, 1, mean), 
        apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), 
        apply(SB112$FitRatios, 1, mean), apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean),
        apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), 
        apply(SB105$FitRatios, 1, mean), apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean), 
        
        apply(SB116$FitRatios, 1, mean), apply(SB117$FitRatios, 1, mean), apply(SB118$FitRatios, 1, mean), 
        apply(SB119$FitRatios, 1, mean), apply(SB120$FitRatios, 1, mean), apply(SB121$FitRatios, 1, mean), 
        
        apply(SB122$FitRatios, 1, mean), apply(SB123$FitRatios, 1, mean), apply(SB124$FitRatios, 1, mean), 
        apply(SB125$FitRatios, 1, mean), apply(SB126$FitRatios, 1, mean), apply(SB127$FitRatios, 1, mean),  
        
        apply(SB128$FitRatios, 1, mean), apply(SB129$FitRatios, 1, mean), apply(SB130$FitRatios, 1, mean), 
        apply(SB131$FitRatios, 1, mean), apply(SB132$FitRatios, 1, mean), apply(SB133$FitRatios, 1, mean), 
        
        apply(SB134$FitRatios, 1, mean), apply(SB135$FitRatios, 1, mean), apply(SB136$FitRatios, 1, mean), 
        apply(SB137$FitRatios, 1, mean), apply(SB138$FitRatios, 1, mean), apply(SB139$FitRatios, 1, mean), 
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=c('None','I-S1','I-S2','I-S3','I-S4','I-S5','I-S6','I-S7','D-S1','D-S2','D-S3','D-S4','D-S5','D-S6','D-S7',
                'I2','I2','I2','D2','D2','D2','I1-D1','I1-D1','I1-D1','D1-I1','D1-I1','D1-I1', 
                'I3','I3','I3','D3','D3','D3','I-D-I','I-D-I','I-D-I','D-I-D','D-I-D','D-I-D'), ylim=c(0,1) )
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
# abline(h=mean(apply(SB101$FitRatios, 1, mean)), col='white')
########

dev.off()



boxplot(SB101$RMSE.BT, 
        SB109$RMSE.BT, SB110$RMSE.BT, SB111$RMSE.BT, SB112$RMSE.BT, SB113$RMSE.BT, SB114$RMSE.BT, SB115$RMSE.BT,
        SB102$RMSE.BT, SB103$RMSE.BT, SB104$RMSE.BT, SB105$RMSE.BT, SB106$RMSE.BT, SB107$RMSE.BT, SB108$RMSE.BT,
        SB116$RMSE.BT, SB117$RMSE.BT, SB118$RMSE.BT, SB119$RMSE.BT, SB120$RMSE.BT, SB121$RMSE.BT,
        SB122$RMSE.BT, SB123$RMSE.BT, SB124$RMSE.BT, SB125$RMSE.BT, SB126$RMSE.BT, SB127$RMSE.BT,
        SB128$RMSE.BT, SB129$RMSE.BT, SB130$RMSE.BT, SB131$RMSE.BT, SB132$RMSE.BT, SB133$RMSE.BT,
        SB134$RMSE.BT, SB135$RMSE.BT, SB136$RMSE.BT, SB137$RMSE.BT, SB138$RMSE.BT, SB139$RMSE.BT,
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=FALSE, ylim=c(0,0.60))

boxplot(apply(SB101$FitRatios, 1, mean), 
        apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), 
        apply(SB112$FitRatios, 1, mean), apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean),
        apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), 
        apply(SB105$FitRatios, 1, mean), apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean), 
        
        apply(SB116$FitRatios, 1, mean), apply(SB117$FitRatios, 1, mean), apply(SB118$FitRatios, 1, mean), 
        apply(SB119$FitRatios, 1, mean), apply(SB120$FitRatios, 1, mean), apply(SB121$FitRatios, 1, mean), 
        
        apply(SB122$FitRatios, 1, mean), apply(SB123$FitRatios, 1, mean), apply(SB124$FitRatios, 1, mean), 
        apply(SB125$FitRatios, 1, mean), apply(SB126$FitRatios, 1, mean), apply(SB127$FitRatios, 1, mean),  
        
        apply(SB128$FitRatios, 1, mean), apply(SB129$FitRatios, 1, mean), apply(SB130$FitRatios, 1, mean), 
        apply(SB131$FitRatios, 1, mean), apply(SB132$FitRatios, 1, mean), apply(SB133$FitRatios, 1, mean), 
        
        apply(SB134$FitRatios, 1, mean), apply(SB135$FitRatios, 1, mean), apply(SB136$FitRatios, 1, mean), 
        apply(SB137$FitRatios, 1, mean), apply(SB138$FitRatios, 1, mean), apply(SB139$FitRatios, 1, mean), 
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=c('None','I-S1','I-S2','I-S3','I-S4','I-S5','I-S6','I-S7','D-S1','D-S2','D-S3','D-S4','D-S5','D-S6','D-S7',
                'I2','I2','I2','D2','D2','D2','I1-D1','I1-D1','I1-D1','D1-I1','D1-I1','D1-I1', 
                'I3','I3','I3','D3','D3','D3','I-D-I','I-D-I','I-D-I','D-I-D','D-I-D','D-I-D'), ylim=c(0,0.7))








png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Backtransformed_MissingData_ChangeQ3.png", 
    type="cairo",
    units="mm", 
    width=600, 
    height=250, 
    pointsize=16, 
    res=600)
##########
# STANDARDIZE AXES
par(mfrow=c(3,1), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB1$biasAll.BT), 
        c(SB9$biasAll.BT), c(SB10$biasAll.BT), c(SB11$biasAll.BT), c(SB12$biasAll.BT), c(SB13$biasAll.BT), c(SB14$biasAll.BT), c(SB15$biasAll.BT),
        c(SB2$biasAll.BT), c(SB3$biasAll.BT), c(SB4$biasAll.BT), c(SB5$biasAll.BT), c(SB6$biasAll.BT), c(SB7$biasAll.BT), c(SB8$biasAll.BT),
        c(SB16$biasAll.BT), c(SB17$biasAll.BT), c(SB18$biasAll.BT), c(SB19$biasAll.BT), c(SB20$biasAll.BT), c(SB21$biasAll.BT), 
        c(SB22$biasAll.BT), c(SB23$biasAll.BT), c(SB24$biasAll.BT), c(SB25$biasAll.BT), c(SB26$biasAll.BT), c(SB27$biasAll.BT),
        c(SB28$biasAll.BT), c(SB29$biasAll.BT), c(SB30$biasAll.BT), c(SB31$biasAll.BT), c(SB32$biasAll.BT), c(SB33$biasAll.BT), 
        c(SB34$biasAll.BT), c(SB35$biasAll.BT), c(SB36$biasAll.BT), c(SB37$biasAll.BT), c(SB38$biasAll.BT), c(SB39$biasAll.BT),
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=FALSE)#, ylim=c(-3.5,3.25) )
mtext('Back-transformed Annual Bias',side=2,line=1.1,cex=0.8)
abline(h=0, col='black')
title('Missing Data', line=-0.9)

vioplot(SB1$RMSE.BT, 
        SB9$RMSE.BT,  SB10$RMSE.BT, SB11$RMSE.BT, SB12$RMSE.BT, SB13$RMSE.BT, SB14$RMSE.BT, SB15$RMSE.BT,
        SB2$RMSE.BT, SB3$RMSE.BT, SB4$RMSE.BT, SB5$RMSE.BT, SB6$RMSE.BT, SB7$RMSE.BT, SB8$RMSE.BT,
        SB16$RMSE.BT, SB17$RMSE.BT, SB18$RMSE.BT, SB19$RMSE.BT, SB20$RMSE.BT, SB21$RMSE.BT,
        SB22$RMSE.BT, SB23$RMSE.BT, SB24$RMSE.BT, SB25$RMSE.BT, SB26$RMSE.BT, SB27$RMSE.BT,
        SB28$RMSE.BT, SB29$RMSE.BT, SB30$RMSE.BT, SB31$RMSE.BT, SB32$RMSE.BT, SB33$RMSE.BT,
        SB34$RMSE.BT, SB35$RMSE.BT, SB36$RMSE.BT, SB37$RMSE.BT, SB38$RMSE.BT, SB39$RMSE.BT,
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=FALSE, ylim=c(0,2))
mtext('Back-transformed RMSE',side=2,line=1.1,cex=0.8)
# abline(h=mean(SB101$RMSE), col='white')


vioplot(apply(SB1$FitRatios, 1, mean), 
        apply(SB9$FitRatios, 1, mean),  apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), 
        apply(SB12$FitRatios, 1, mean), apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean),
        apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), 
        apply(SB5$FitRatios, 1, mean), apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean), 
        
        apply(SB16$FitRatios, 1, mean), apply(SB17$FitRatios, 1, mean), apply(SB18$FitRatios, 1, mean), 
        apply(SB19$FitRatios, 1, mean), apply(SB20$FitRatios, 1, mean), apply(SB21$FitRatios, 1, mean), 
        apply(SB22$FitRatios, 1, mean), apply(SB23$FitRatios, 1, mean), apply(SB24$FitRatios, 1, mean), 
        apply(SB25$FitRatios, 1, mean), apply(SB26$FitRatios, 1, mean), apply(SB27$FitRatios, 1, mean), 
        
        apply(SB28$FitRatios, 1, mean), apply(SB29$FitRatios, 1, mean), apply(SB30$FitRatios, 1, mean), 
        apply(SB31$FitRatios, 1, mean), apply(SB32$FitRatios, 1, mean), apply(SB33$FitRatios, 1, mean), 
        apply(SB34$FitRatios, 1, mean), apply(SB35$FitRatios, 1, mean), apply(SB36$FitRatios, 1, mean), 
        apply(SB37$FitRatios, 1, mean), apply(SB38$FitRatios, 1, mean), apply(SB39$FitRatios, 1, mean), 
        col=c('grey',
              'olivedrab1','darkolivegreen1','darkolivegreen2','darkolivegreen3','olivedrab3','olivedrab','darkolivegreen',
              'skyblue1','skyblue2','skyblue','deepskyblue','deepskyblue2','deepskyblue3','deepskyblue4',
              'darkseagreen1','darkseagreen2','darkseagreen3','darkslategray1','darkslategray3','darkslategray4',
              'orchid1','mediumorchid1','mediumorchid3','darkorchid1','darkorchid3','darkorchid4',
              'springgreen','springgreen3','springgreen4', 'steelblue1','steelblue3','steelblue4',
              'mediumpurple1','mediumpurple3','mediumpurple4','purple','purple3','purple4'),
        names=c('None','I-S1','I-S2','I-S3','I-S4','I-S5','I-S6','I-S7','D-S1','D-S2','D-S3','D-S4','D-S5','D-S6','D-S7',
                'I2','I2','I2','D2','D2','D2','I1-D1','I1-D1','I1-D1','D1-I1','D1-I1','D1-I1', 
                'I3','I3','I3','D3','D3','D3','I-D-I','I-D-I','I-D-I','D-I-D','D-I-D','D-I-D'), ylim=c(0,1) )
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
# abline(h=mean(apply(SB1$FitRatios, 1, mean)), col='white')
#############

dev.off()







png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Backtransformed_FullVSMissingData_3_2.png", 
    type="cairo",
    units="mm", 
    width=250, 
    height=250, 
    pointsize=16, 
    res=600)
#############
# STANDARDIZE AXES
par(mfrow=c(3,2), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll.BT), 
        c(SB109$biasAll.BT, SB110$biasAll.BT, SB111$biasAll.BT, SB112$biasAll.BT, SB113$biasAll.BT, SB114$biasAll.BT, SB115$biasAll.BT),
        c(SB102$biasAll.BT, SB103$biasAll.BT, SB104$biasAll.BT, SB105$biasAll.BT, SB106$biasAll.BT, SB107$biasAll.BT, SB108$biasAll.BT),
        c(SB116$biasAll.BT, SB117$biasAll.BT, SB118$biasAll.BT), c(SB119$biasAll.BT, SB120$biasAll.BT, SB121$biasAll.BT), 
        c(SB122$biasAll.BT, SB123$biasAll.BT, SB124$biasAll.BT), c(SB125$biasAll.BT, SB126$biasAll.BT, SB127$biasAll.BT),
        c(SB128$biasAll.BT, SB129$biasAll.BT, SB130$biasAll.BT), c(SB131$biasAll.BT, SB132$biasAll.BT, SB133$biasAll.BT), 
        c(SB134$biasAll.BT, SB135$biasAll.BT, SB136$biasAll.BT), c(SB137$biasAll.BT, SB138$biasAll.BT, SB139$biasAll.BT),
        col=c('grey', 'darkolivegreen1', 'deepskyblue', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3',
              'springgreen3', 'steelblue3','mediumpurple3','purple'),
        names=FALSE, ylim=c(-4,3.5))
mtext('Back-transformed Annual Bias',side=2,line=1.1,cex=0.8)
abline(h=0, col='black')
title('Complete Data', line=-0.9)
vioplot(c(SB1$biasAll.BT), 
        c(SB9$biasAll.BT,  SB10$biasAll.BT, SB11$biasAll.BT, SB12$biasAll.BT, SB13$biasAll.BT, SB14$biasAll.BT, SB15$biasAll.BT),
        c(SB2$biasAll.BT, SB3$biasAll.BT, SB4$biasAll.BT, SB5$biasAll.BT, SB6$biasAll.BT, SB7$biasAll.BT, SB8$biasAll.BT),
        c(SB16$biasAll.BT, SB17$biasAll.BT, SB18$biasAll.BT), c(SB19$biasAll.BT, SB20$biasAll.BT, SB21$biasAll.BT),
        c(SB22$biasAll.BT, SB23$biasAll.BT, SB24$biasAll.BT), c(SB25$biasAll.BT, SB26$biasAll.BT, SB27$biasAll.BT),
        c(SB28$biasAll.BT, SB29$biasAll.BT, SB30$biasAll.BT), c(SB31$biasAll.BT, SB32$biasAll.BT, SB33$biasAll.BT),
        c(SB34$biasAll.BT, SB35$biasAll.BT, SB36$biasAll.BT), c(SB37$biasAll.BT, SB38$biasAll.BT, SB39$biasAll.BT),
        col=c('grey', 'darkolivegreen1', 'deepskyblue', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3',
              'springgreen3', 'steelblue3','mediumpurple3','purple'),
        names=FALSE, ylim=c(-4,3.5))
abline(h=0, col='black')
title('Missing Data', line=-0.9)

vioplot(SB101$RMSE.BT, 
        c(SB109$RMSE.BT, SB110$RMSE.BT, SB111$RMSE.BT, SB112$RMSE.BT, SB113$RMSE.BT, SB114$RMSE.BT, SB115$RMSE.BT),
        c(SB102$RMSE.BT, SB103$RMSE.BT, SB104$RMSE.BT, SB105$RMSE.BT, SB106$RMSE.BT, SB107$RMSE.BT, SB108$RMSE.BT),
        c(SB116$RMSE.BT, SB117$RMSE.BT, SB118$RMSE.BT), c(SB119$RMSE.BT, SB120$RMSE.BT, SB121$RMSE.BT), 
        c(SB122$RMSE.BT, SB123$RMSE.BT, SB124$RMSE.BT), c(SB125$RMSE.BT, SB126$RMSE.BT, SB127$RMSE.BT),
        c(SB128$RMSE.BT, SB129$RMSE.BT, SB130$RMSE.BT), c(SB131$RMSE.BT, SB132$RMSE.BT, SB133$RMSE.BT), 
        c(SB134$RMSE.BT, SB135$RMSE.BT, SB136$RMSE.BT), c(SB137$RMSE.BT, SB138$RMSE.BT, SB139$RMSE.BT),
        col=c('grey', 'darkolivegreen1', 'deepskyblue', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3',
              'springgreen3', 'steelblue3','mediumpurple3','purple'),
        names=FALSE, ylim=c(0,2))
mtext('Back-transformed RMSE',side=2,line=1.1,cex=0.8)
vioplot(SB1$RMSE.BT, 
        c(SB9$RMSE.BT, SB10$RMSE.BT, SB11$RMSE.BT, SB12$RMSE.BT, SB13$RMSE.BT, SB14$RMSE.BT, SB15$RMSE.BT),
        c(SB2$RMSE.BT, SB3$RMSE.BT, SB4$RMSE.BT, SB5$RMSE.BT, SB6$RMSE.BT, SB7$RMSE.BT, SB8$RMSE.BT),
        c(SB16$RMSE.BT, SB17$RMSE.BT, SB18$RMSE.BT), c(SB19$RMSE.BT, SB20$RMSE.BT, SB21$RMSE.BT), 
        c(SB22$RMSE.BT, SB23$RMSE.BT, SB24$RMSE.BT), c(SB25$RMSE.BT, SB26$RMSE.BT, SB27$RMSE.BT),
        c(SB28$RMSE.BT, SB29$RMSE.BT, SB30$RMSE.BT), c(SB31$RMSE.BT, SB32$RMSE.BT, SB33$RMSE.BT), 
        c(SB34$RMSE.BT, SB35$RMSE.BT, SB36$RMSE.BT), c(SB37$RMSE.BT, SB38$RMSE.BT, SB39$RMSE.BT),
        col=c('grey', 'darkolivegreen1', 'deepskyblue', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3',
              'springgreen3', 'steelblue3','mediumpurple3','purple'),
        names=FALSE, ylim=c(0,2))


vioplot(apply(SB101$FitRatios, 1, mean),
        c(apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), apply(SB112$FitRatios, 1, mean), 
          apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean)), 
        c(apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), apply(SB105$FitRatios, 1, mean), 
          apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean)), 
        
        c(apply(SB116$FitRatios, 1, mean), apply(SB117$FitRatios, 1, mean), apply(SB118$FitRatios, 1, mean)), 
        c(apply(SB119$FitRatios, 1, mean), apply(SB120$FitRatios, 1, mean), apply(SB121$FitRatios, 1, mean)), 
        
        c(apply(SB122$FitRatios, 1, mean), apply(SB123$FitRatios, 1, mean), apply(SB124$FitRatios, 1, mean)), 
        c(apply(SB125$FitRatios, 1, mean), apply(SB126$FitRatios, 1, mean), apply(SB127$FitRatios, 1, mean)), 
        
        c(apply(SB128$FitRatios, 1, mean), apply(SB129$FitRatios, 1, mean), apply(SB130$FitRatios, 1, mean)), 
        c(apply(SB131$FitRatios, 1, mean), apply(SB132$FitRatios, 1, mean), apply(SB133$FitRatios, 1, mean)), 
        
        c(apply(SB134$FitRatios, 1, mean), apply(SB135$FitRatios, 1, mean), apply(SB136$FitRatios, 1, mean)), 
        c(apply(SB137$FitRatios, 1, mean), apply(SB138$FitRatios, 1, mean), apply(SB139$FitRatios, 1, mean)),
        col=c('grey', 'darkolivegreen1', 'deepskyblue', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3',
              'springgreen3', 'steelblue3','mediumpurple3','purple'),
        names=c(
          'No
          change','I1','D1','I2','D2','I-D','D-I','I3','D3','I-D-I','D-I-D'), ylim=c(0,1))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)

vioplot(apply(SB1$FitRatios, 1, mean),
        c(apply(SB9$FitRatios, 1, mean), apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), apply(SB12$FitRatios, 1, mean), 
          apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean)), 
        c(apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), apply(SB5$FitRatios, 1, mean), 
          apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean)), 
        
        c(apply(SB16$FitRatios, 1, mean), apply(SB17$FitRatios, 1, mean), apply(SB18$FitRatios, 1, mean)), 
        c(apply(SB19$FitRatios, 1, mean), apply(SB20$FitRatios, 1, mean), apply(SB21$FitRatios, 1, mean)), 
        
        c(apply(SB22$FitRatios, 1, mean), apply(SB23$FitRatios, 1, mean), apply(SB24$FitRatios, 1, mean)), 
        c(apply(SB25$FitRatios, 1, mean), apply(SB26$FitRatios, 1, mean), apply(SB27$FitRatios, 1, mean)), 
        
        c(apply(SB28$FitRatios, 1, mean), apply(SB29$FitRatios, 1, mean), apply(SB30$FitRatios, 1, mean)), 
        c(apply(SB31$FitRatios, 1, mean), apply(SB32$FitRatios, 1, mean), apply(SB33$FitRatios, 1, mean)), 
        
        c(apply(SB34$FitRatios, 1, mean), apply(SB35$FitRatios, 1, mean), apply(SB36$FitRatios, 1, mean)), 
        c(apply(SB37$FitRatios, 1, mean), apply(SB38$FitRatios, 1, mean), apply(SB39$FitRatios, 1, mean)),
        col=c('grey', 'darkolivegreen1', 'deepskyblue', 'darkslategray3','darkseagreen2','mediumorchid1','darkorchid3',
              'springgreen3', 'steelblue3','mediumpurple3','purple'),
        names=c(
          'No
          change','I1','D1','I2','D2','I-D','D-I','I3','D3','I-D-I','D-I-D'), ylim=c(0,1))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
##########

dev.off()







png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Backtransformed_FullVSMissingData_3_3.png", 
    type="cairo",
    units="mm", 
    width=250, 
    height=250, 
    pointsize=16, 
    res=600)
#############
# STANDARDIZE AXES
par(mfrow=c(3,1), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
vioplot(c(SB101$biasAll.BT), 
        c(SB1$biasAll.BT), 
        c(SB109$biasAll.BT, SB110$biasAll.BT, SB111$biasAll.BT, SB112$biasAll.BT, SB113$biasAll.BT, SB114$biasAll.BT, SB115$biasAll.BT),
        c(SB9$biasAll.BT,  SB10$biasAll.BT, SB11$biasAll.BT, SB12$biasAll.BT, SB13$biasAll.BT, SB14$biasAll.BT, SB15$biasAll.BT),
        c(SB102$biasAll.BT, SB103$biasAll.BT, SB104$biasAll.BT, SB105$biasAll.BT, SB106$biasAll.BT, SB107$biasAll.BT, SB108$biasAll.BT),
        c(SB2$biasAll.BT, SB3$biasAll.BT, SB4$biasAll.BT, SB5$biasAll.BT, SB6$biasAll.BT, SB7$biasAll.BT, SB8$biasAll.BT),
        c(SB116$biasAll.BT, SB117$biasAll.BT, SB118$biasAll.BT),
        c(SB16$biasAll.BT, SB17$biasAll.BT, SB18$biasAll.BT),
        c(SB119$biasAll.BT, SB120$biasAll.BT, SB121$biasAll.BT), 
        c(SB19$biasAll.BT, SB20$biasAll.BT, SB21$biasAll.BT), 
        c(SB122$biasAll.BT, SB123$biasAll.BT, SB124$biasAll.BT), 
        c(SB22$biasAll.BT, SB23$biasAll.BT, SB24$biasAll.BT), 
        c(SB125$biasAll.BT, SB126$biasAll.BT, SB127$biasAll.BT),
        c(SB25$biasAll.BT, SB26$biasAll.BT, SB27$biasAll.BT),
        
        c(SB128$biasAll.BT, SB129$biasAll.BT, SB130$biasAll.BT),
        c(SB28$biasAll.BT, SB28$biasAll.BT, SB30$biasAll.BT),
        c(SB131$biasAll.BT, SB132$biasAll.BT, SB133$biasAll.BT), 
        c(SB131$biasAll.BT, SB32$biasAll.BT, SB33$biasAll.BT), 
        c(SB134$biasAll.BT, SB135$biasAll.BT, SB136$biasAll.BT), 
        c(SB34$biasAll.BT, SB35$biasAll.BT, SB36$biasAll.BT), 
        c(SB137$biasAll.BT, SB138$biasAll.BT, SB139$biasAll.BT),
        c(SB37$biasAll.BT, SB38$biasAll.BT, SB39$biasAll.BT),
        col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
              'darkslategray4','darkslategray1','mediumorchid1','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
              'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple'),
        # col=c('lightgrey','darkgrey','olivedrab1','olivedrab3','skyblue1','deepskyblue','darkseagreen1','darkseagreen4',
        #       'darkslategray1','darkslategray4','orchid','mediumorchid1','darkorchid1','darkorchid4','springgreen','springgreen4',
        #       'steelblue1','steelblue4','mediumpurple1','mediumpurple4','purple','purple4'),
        names=FALSE, ylim=c(-4,3.5))
mtext('Back-transformed Annual Bias',side=2,line=1.1,cex=0.8)
title('Complete vs. Missing Data', line=-0.9)
abline(h=0, col='black')


vioplot(SB101$RMSE.BT, 
        SB1$RMSE.BT,
        c(SB109$RMSE.BT, SB110$RMSE.BT, SB111$RMSE.BT, SB112$RMSE.BT, SB113$RMSE.BT, SB114$RMSE.BT, SB115$RMSE.BT),
        c(SB9$RMSE.BT, SB10$RMSE.BT, SB11$RMSE.BT, SB12$RMSE.BT, SB13$RMSE.BT, SB14$RMSE.BT, SB15$RMSE.BT),
        c(SB102$RMSE.BT, SB103$RMSE.BT, SB104$RMSE.BT, SB105$RMSE.BT, SB106$RMSE.BT, SB107$RMSE.BT, SB108$RMSE.BT),
        c(SB2$RMSE.BT, SB3$RMSE.BT, SB4$RMSE.BT, SB5$RMSE.BT, SB6$RMSE.BT, SB7$RMSE.BT, SB8$RMSE.BT), 
        
        c(SB116$RMSE.BT, SB117$RMSE.BT, SB118$RMSE.BT),
        c(SB16$RMSE.BT, SB17$RMSE.BT, SB18$RMSE.BT),
        c(SB119$RMSE.BT, SB120$RMSE.BT, SB121$RMSE.BT), 
        c(SB19$RMSE.BT, SB20$RMSE.BT, SB21$RMSE.BT),
        c(SB122$RMSE.BT, SB123$RMSE.BT, SB124$RMSE.BT), 
        c(SB22$RMSE.BT, SB23$RMSE.BT, SB24$RMSE.BT),
        c(SB125$RMSE.BT, SB126$RMSE.BT, SB127$RMSE.BT),
        c(SB25$RMSE.BT, SB26$RMSE.BT, SB27$RMSE.BT),
        
        c(SB128$RMSE.BT, SB129$RMSE.BT, SB130$RMSE.BT),
        c(SB28$RMSE.BT, SB29$RMSE.BT, SB30$RMSE.BT),
        c(SB131$RMSE.BT, SB132$RMSE.BT, SB133$RMSE.BT), 
        c(SB31$RMSE.BT, SB32$RMSE.BT, SB33$RMSE.BT),
        c(SB134$RMSE.BT, SB135$RMSE.BT, SB136$RMSE.BT), 
        c(SB34$RMSE.BT, SB35$RMSE.BT, SB36$RMSE.BT),
        c(SB137$RMSE.BT, SB138$RMSE.BT, SB139$RMSE.BT),
        c(SB37$RMSE.BT, SB38$RMSE.BT, SB39$RMSE.BT),
        col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
              'darkslategray4','darkslategray1','mediumorchid1','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
              'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple'),
        # col=c('lightgrey','darkgrey','olivedrab1','olivedrab3','skyblue1','deepskyblue','darkseagreen1','darkseagreen4',
        #       'darkslategray1','darkslategray4','orchid','mediumorchid1','darkorchid1','darkorchid4','springgreen','springgreen4',
        #       'steelblue1','steelblue4','mediumpurple1','mediumpurple4','purple','purple4'),
        names=FALSE, ylim=c(0,2))
mtext('Back-transformed RMSE',side=2,line=1.1,cex=0.8)


vioplot(apply(SB101$FitRatios, 1, mean),
        apply(SB1$FitRatios, 1, mean),
        c(apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), apply(SB112$FitRatios, 1, mean), 
          apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean)), 
        c(apply(SB9$FitRatios, 1, mean), apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), apply(SB12$FitRatios, 1, mean), 
          apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean)), 
        c(apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), apply(SB105$FitRatios, 1, mean), 
          apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean)), 
        c(apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), apply(SB5$FitRatios, 1, mean), 
          apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean)), 
        
        c(apply(SB116$FitRatios, 1, mean), apply(SB117$FitRatios, 1, mean), apply(SB118$FitRatios, 1, mean)), 
        c(apply(SB16$FitRatios, 1, mean), apply(SB17$FitRatios, 1, mean), apply(SB18$FitRatios, 1, mean)), 
        c(apply(SB119$FitRatios, 1, mean), apply(SB120$FitRatios, 1, mean), apply(SB121$FitRatios, 1, mean)), 
        c(apply(SB19$FitRatios, 1, mean), apply(SB20$FitRatios, 1, mean), apply(SB21$FitRatios, 1, mean)), 
        c(apply(SB122$FitRatios, 1, mean), apply(SB123$FitRatios, 1, mean), apply(SB124$FitRatios, 1, mean)), 
        c(apply(SB22$FitRatios, 1, mean), apply(SB23$FitRatios, 1, mean), apply(SB24$FitRatios, 1, mean)), 
        c(apply(SB125$FitRatios, 1, mean), apply(SB126$FitRatios, 1, mean), apply(SB127$FitRatios, 1, mean)), 
        c(apply(SB25$FitRatios, 1, mean), apply(SB26$FitRatios, 1, mean), apply(SB27$FitRatios, 1, mean)), 
        
        c(apply(SB128$FitRatios, 1, mean), apply(SB129$FitRatios, 1, mean), apply(SB130$FitRatios, 1, mean)), 
        c(apply(SB28$FitRatios, 1, mean), apply(SB29$FitRatios, 1, mean), apply(SB30$FitRatios, 1, mean)), 
        c(apply(SB131$FitRatios, 1, mean), apply(SB132$FitRatios, 1, mean), apply(SB133$FitRatios, 1, mean)), 
        c(apply(SB31$FitRatios, 1, mean), apply(SB32$FitRatios, 1, mean), apply(SB33$FitRatios, 1, mean)), 
        c(apply(SB134$FitRatios, 1, mean), apply(SB135$FitRatios, 1, mean), apply(SB136$FitRatios, 1, mean)), 
        c(apply(SB34$FitRatios, 1, mean), apply(SB35$FitRatios, 1, mean), apply(SB36$FitRatios, 1, mean)), 
        c(apply(SB137$FitRatios, 1, mean), apply(SB138$FitRatios, 1, mean), apply(SB139$FitRatios, 1, mean)), 
        c(apply(SB37$FitRatios, 1, mean), apply(SB38$FitRatios, 1, mean), apply(SB39$FitRatios, 1, mean)), 
        
        col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
              'darkslategray4','darkslategray1','mediumorchid1','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
              'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple'),
        # col=c('lightgrey','darkgrey','olivedrab1','olivedrab3','skyblue1','deepskyblue','darkseagreen1','darkseagreen4',
        #       'darkslategray1','darkslategray4','orchid','mediumorchid1','darkorchid1','darkorchid4','springgreen','springgreen4',
        #       'steelblue1','steelblue4','mediumpurple1','mediumpurple4','purple','purple4'),
        names=c(
          'No
          change','','I1','','D1','','I2','','D2','','I-D','','D-I','','I3','','D3','','I-D-I','','D-I-D',''), ylim=c(0,1))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
# abline(h=0.6)
#############
dev.off()







png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Backtransformed_FullVSMissingData_3_1_LgFont.png", 
    type="cairo",
    units="mm", 
    width=250, 
    height=250, 
    pointsize=16, 
    res=600)
##########
# STANDARDIZE AXES
par(mfrow=c(3,1), mar=c(0.6, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0), cex=1.0)
vioplot(c(SB101$biasAll.BT), 
        c(SB1$biasAll.BT), 
        
        c(SB109$biasAll.BT, SB110$biasAll.BT, SB111$biasAll.BT, SB112$biasAll.BT, SB113$biasAll.BT, SB114$biasAll.BT, SB115$biasAll.BT, 
          SB102$biasAll.BT, SB103$biasAll.BT, SB104$biasAll.BT, SB105$biasAll.BT, SB106$biasAll.BT, SB107$biasAll.BT, SB108$biasAll.BT),
        c(SB9$biasAll.BT,  SB10$biasAll.BT, SB11$biasAll.BT, SB12$biasAll.BT, SB13$biasAll.BT, SB14$biasAll.BT, SB15$biasAll.BT, 
          SB2$biasAll.BT, SB3$biasAll.BT, SB4$biasAll.BT, SB5$biasAll.BT, SB6$biasAll.BT, SB7$biasAll.BT, SB8$biasAll.BT),
        
        c(SB116$biasAll.BT, SB117$biasAll.BT, SB118$biasAll.BT, SB119$biasAll.BT, SB120$biasAll.BT, SB121$biasAll.BT, 
          SB122$biasAll.BT, SB123$biasAll.BT, SB124$biasAll.BT, SB125$biasAll.BT, SB126$biasAll.BT, SB127$biasAll.BT),
        c(SB16$biasAll.BT, SB17$biasAll.BT, SB18$biasAll.BT, SB19$biasAll.BT, SB20$biasAll.BT, SB21$biasAll.BT, 
          SB22$biasAll.BT, SB23$biasAll.BT, SB24$biasAll.BT, SB25$biasAll.BT, SB26$biasAll.BT, SB27$biasAll.BT),
        
        c(SB128$biasAll.BT, SB129$biasAll.BT, SB130$biasAll.BT, SB131$biasAll.BT, SB132$biasAll.BT, SB133$biasAll.BT, 
          SB134$biasAll.BT, SB135$biasAll.BT, SB136$biasAll.BT, SB137$biasAll.BT, SB138$biasAll.BT, SB139$biasAll.BT),
        c(SB28$biasAll.BT, SB28$biasAll.BT, SB30$biasAll.BT, SB131$biasAll.BT, SB32$biasAll.BT, SB33$biasAll.BT, 
          SB34$biasAll.BT, SB35$biasAll.BT, SB36$biasAll.BT, SB37$biasAll.BT, SB38$biasAll.BT, SB39$biasAll.BT),
        col=c('firebrick','firebrick1', 'deepskyblue4','deepskyblue','olivedrab','olivedrab1','darkorchid4','darkorchid'),
        names=FALSE, ylim=c(-4,3.5))
mtext('Back-transformed Annual Bias',side=2,line=1.1,cex=1.0)
title('Complete vs. Missing Data', line=-0.9, font.main=1, cex=0.9)
abline(h=0, col='black')
legend('bottomleft',c('Complete data','Missing data'),fill=c('darkgrey','lightgrey'), bty='n', cex=0.8, border=FALSE)


vioplot(SB101$RMSE.BT, 
        SB1$RMSE.BT,
        
        c(SB109$RMSE.BT, SB110$RMSE.BT, SB111$RMSE.BT, SB112$RMSE.BT, SB113$RMSE.BT, SB114$RMSE.BT, SB115$RMSE.BT, 
          SB102$RMSE.BT, SB103$RMSE.BT, SB104$RMSE.BT, SB105$RMSE.BT, SB106$RMSE.BT, SB107$RMSE.BT, SB108$RMSE.BT),
        c(SB9$RMSE.BT, SB10$RMSE.BT, SB11$RMSE.BT, SB12$RMSE.BT, SB13$RMSE.BT, SB14$RMSE.BT, SB15$RMSE.BT, 
          SB2$RMSE.BT, SB3$RMSE.BT, SB4$RMSE.BT, SB5$RMSE.BT, SB6$RMSE.BT, SB7$RMSE.BT, SB8$RMSE.BT), 
        
        c(SB116$RMSE.BT, SB117$RMSE.BT, SB118$RMSE.BT, SB119$RMSE.BT, SB120$RMSE.BT, SB121$RMSE.BT, 
          SB122$RMSE.BT, SB123$RMSE.BT, SB124$RMSE.BT, SB125$RMSE.BT, SB126$RMSE.BT, SB127$RMSE.BT),
        c(SB16$RMSE.BT, SB17$RMSE.BT, SB18$RMSE.BT, SB19$RMSE.BT, SB20$RMSE.BT, SB21$RMSE.BT, 
          SB22$RMSE.BT, SB23$RMSE.BT, SB24$RMSE.BT, SB25$RMSE.BT, SB26$RMSE.BT, SB27$RMSE.BT),
        
        c(SB128$RMSE.BT, SB129$RMSE.BT, SB130$RMSE.BT, SB131$RMSE.BT, SB132$RMSE.BT, SB133$RMSE.BT, 
          SB134$RMSE.BT, SB135$RMSE.BT, SB136$RMSE.BT, SB137$RMSE.BT, SB138$RMSE.BT, SB139$RMSE.BT),
        c(SB28$RMSE.BT, SB29$RMSE.BT, SB30$RMSE.BT, SB31$RMSE.BT, SB32$RMSE.BT, SB33$RMSE.BT, 
          SB34$RMSE.BT, SB35$RMSE.BT, SB36$RMSE.BT, SB37$RMSE.BT, SB38$RMSE.BT, SB39$RMSE.BT),
        col=c('firebrick','firebrick1', 'deepskyblue4','deepskyblue','olivedrab','olivedrab1','darkorchid4','darkorchid'),
        names=FALSE, ylim=c(0,2))
# mtext('RMSE',side=2,line=1.1,cex=0.8)
mtext('Back-transformed RMSE',side=2,line=1.1,cex=1.0)

par(mar=c(1.6, 2.1, 0.1, 0.3))
vioplot(apply(SB101$FitRatios, 1, mean),
        apply(SB1$FitRatios, 1, mean),
        
        c(apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), apply(SB112$FitRatios, 1, mean), 
          apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean), 
          apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), apply(SB105$FitRatios, 1, mean), 
          apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean)), 
        c(apply(SB9$FitRatios, 1, mean), apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), apply(SB12$FitRatios, 1, mean), 
          apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean), 
          apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), apply(SB5$FitRatios, 1, mean), 
          apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean)), 
        
        c(apply(SB116$FitRatios, 1, mean), apply(SB117$FitRatios, 1, mean), apply(SB118$FitRatios, 1, mean), 
          apply(SB119$FitRatios, 1, mean), apply(SB120$FitRatios, 1, mean), apply(SB121$FitRatios, 1, mean), 
          apply(SB122$FitRatios, 1, mean), apply(SB123$FitRatios, 1, mean), apply(SB124$FitRatios, 1, mean), 
          apply(SB125$FitRatios, 1, mean), apply(SB126$FitRatios, 1, mean), apply(SB127$FitRatios, 1, mean)), 
        c(apply(SB16$FitRatios, 1, mean), apply(SB17$FitRatios, 1, mean), apply(SB18$FitRatios, 1, mean), 
          apply(SB19$FitRatios, 1, mean), apply(SB20$FitRatios, 1, mean), apply(SB21$FitRatios, 1, mean), 
          apply(SB22$FitRatios, 1, mean), apply(SB23$FitRatios, 1, mean), apply(SB24$FitRatios, 1, mean), 
          apply(SB25$FitRatios, 1, mean), apply(SB26$FitRatios, 1, mean), apply(SB27$FitRatios, 1, mean)), 
        
        c(apply(SB128$FitRatios, 1, mean), apply(SB129$FitRatios, 1, mean), apply(SB130$FitRatios, 1, mean), 
          apply(SB131$FitRatios, 1, mean), apply(SB132$FitRatios, 1, mean), apply(SB133$FitRatios, 1, mean), 
          apply(SB134$FitRatios, 1, mean), apply(SB135$FitRatios, 1, mean), apply(SB136$FitRatios, 1, mean), 
          apply(SB137$FitRatios, 1, mean), apply(SB138$FitRatios, 1, mean), apply(SB139$FitRatios, 1, mean)), 
        c(apply(SB28$FitRatios, 1, mean), apply(SB29$FitRatios, 1, mean), apply(SB30$FitRatios, 1, mean), 
          apply(SB31$FitRatios, 1, mean), apply(SB32$FitRatios, 1, mean), apply(SB33$FitRatios, 1, mean), 
          apply(SB34$FitRatios, 1, mean), apply(SB35$FitRatios, 1, mean), apply(SB36$FitRatios, 1, mean), 
          apply(SB37$FitRatios, 1, mean), apply(SB38$FitRatios, 1, mean), apply(SB39$FitRatios, 1, mean)), 
        
        col=c('firebrick','firebrick1', 'deepskyblue4','deepskyblue','olivedrab','olivedrab1','darkorchid4','darkorchid'),
        # names=c('              No change','','              Change 1 q','','              Change 2 q','','               Change 3 q',''), 
        names=c('         No change','','      Change 1 q','','       Change 2 q','','     Change 3 q',''), 
        ylim=c(0,1))
# mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)
mtext('Mean Fit Ratio',side=2,line=1.1,cex=1.0)
# legend('bottomright',c('Complete data','Missing data'),fill=c('darkgrey','lightgrey'), bty='n', cex=0.8)
# legend('bottomright',c('Complete data','Missing data'),fill=c('darkgrey','lightgrey'), bty='n', cex=0.8)
#############

dev.off()


# STOP --------------------------------------------------------------------------------------------------

#### TRY WITH BOXPLOTS ####

# col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
#       'darkslategray4','darkslategray1','mediumorchid1','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
#       'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple')


par(mfrow=c(3,1), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
boxplot(c(SB101$biasAll.BT), 
        c(SB1$biasAll.BT), 
        c(SB109$biasAll.BT, SB110$biasAll.BT, SB111$biasAll.BT, SB112$biasAll.BT, SB113$biasAll.BT, SB114$biasAll.BT, SB115$biasAll.BT),
        c(SB9$biasAll.BT,  SB10$biasAll.BT, SB11$biasAll.BT, SB12$biasAll.BT, SB13$biasAll.BT, SB14$biasAll.BT, SB15$biasAll.BT),
        c(SB102$biasAll.BT, SB103$biasAll.BT, SB104$biasAll.BT, SB105$biasAll.BT, SB106$biasAll.BT, SB107$biasAll.BT, SB108$biasAll.BT),
        c(SB2$biasAll.BT, SB3$biasAll.BT, SB4$biasAll.BT, SB5$biasAll.BT, SB6$biasAll.BT, SB7$biasAll.BT, SB8$biasAll.BT),
        c(SB116$biasAll.BT, SB117$biasAll.BT, SB118$biasAll.BT),
        c(SB16$biasAll.BT, SB17$biasAll.BT, SB18$biasAll.BT),
        c(SB119$biasAll.BT, SB120$biasAll.BT, SB121$biasAll.BT), 
        c(SB19$biasAll.BT, SB20$biasAll.BT, SB21$biasAll.BT), 
        c(SB122$biasAll.BT, SB123$biasAll.BT, SB124$biasAll.BT), 
        c(SB22$biasAll.BT, SB23$biasAll.BT, SB24$biasAll.BT), 
        c(SB125$biasAll.BT, SB126$biasAll.BT, SB127$biasAll.BT),
        c(SB25$biasAll.BT, SB26$biasAll.BT, SB27$biasAll.BT),
        
        c(SB128$biasAll.BT, SB129$biasAll.BT, SB130$biasAll.BT),
        c(SB28$biasAll.BT, SB28$biasAll.BT, SB30$biasAll.BT),
        c(SB131$biasAll.BT, SB132$biasAll.BT, SB133$biasAll.BT), 
        c(SB131$biasAll.BT, SB32$biasAll.BT, SB33$biasAll.BT), 
        c(SB134$biasAll.BT, SB135$biasAll.BT, SB136$biasAll.BT), 
        c(SB34$biasAll.BT, SB35$biasAll.BT, SB36$biasAll.BT), 
        c(SB137$biasAll.BT, SB138$biasAll.BT, SB139$biasAll.BT),
        c(SB37$biasAll.BT, SB38$biasAll.BT, SB39$biasAll.BT),
        col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
              'darkslategray4','darkslategray1','orchid4','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
              'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple'),
        # col=c('lightgrey','darkgrey','olivedrab1','olivedrab3','skyblue1','deepskyblue','darkseagreen1','darkseagreen4',
        #       'darkslategray1','darkslategray4','orchid','orchid4','darkorchid1','darkorchid4','springgreen','springgreen4',
        #       'steelblue1','steelblue4','mediumpurple1','mediumpurple4','purple','purple4'),
        names=FALSE, ylim=c(-5,3))
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
title('Complete vs. Missing Data', line=-0.9)
abline(h=0, col='black')


boxplot(SB101$RMSE, 
        SB1$RMSE,
        c(SB109$RMSE, SB110$RMSE, SB111$RMSE, SB112$RMSE, SB113$RMSE, SB114$RMSE, SB115$RMSE),
        c(SB9$RMSE, SB10$RMSE, SB11$RMSE, SB12$RMSE, SB13$RMSE, SB14$RMSE, SB15$RMSE),
        c(SB102$RMSE, SB103$RMSE, SB104$RMSE, SB105$RMSE, SB106$RMSE, SB107$RMSE, SB108$RMSE),
        c(SB2$RMSE, SB3$RMSE, SB4$RMSE, SB5$RMSE, SB6$RMSE, SB7$RMSE, SB8$RMSE), 
        
        c(SB116$RMSE, SB117$RMSE, SB118$RMSE),
        c(SB16$RMSE, SB17$RMSE, SB18$RMSE),
        c(SB119$RMSE, SB120$RMSE, SB121$RMSE), 
        c(SB19$RMSE, SB20$RMSE, SB21$RMSE),
        c(SB122$RMSE, SB123$RMSE, SB124$RMSE), 
        c(SB22$RMSE, SB23$RMSE, SB24$RMSE),
        c(SB125$RMSE, SB126$RMSE, SB127$RMSE),
        c(SB25$RMSE, SB26$RMSE, SB27$RMSE),
        
        c(SB128$RMSE, SB129$RMSE, SB130$RMSE),
        c(SB28$RMSE, SB29$RMSE, SB30$RMSE),
        c(SB131$RMSE, SB132$RMSE, SB133$RMSE), 
        c(SB31$RMSE, SB32$RMSE, SB33$RMSE),
        c(SB134$RMSE, SB135$RMSE, SB136$RMSE), 
        c(SB34$RMSE, SB35$RMSE, SB36$RMSE),
        c(SB137$RMSE, SB138$RMSE, SB139$RMSE),
        c(SB37$RMSE, SB38$RMSE, SB39$RMSE),
        col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
              'darkslategray4','darkslategray1','orchid4','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
              'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple'),
        # col=c('lightgrey','darkgrey','olivedrab1','olivedrab3','skyblue1','deepskyblue','darkseagreen1','darkseagreen4',
        #       'darkslategray1','darkslategray4','orchid','orchid4','darkorchid1','darkorchid4','springgreen','springgreen4',
        #       'steelblue1','steelblue4','mediumpurple1','mediumpurple4','purple','purple4'),
        names=FALSE, ylim=c(0,2))
mtext('RMSE',side=2,line=1.1,cex=0.8)


boxplot(apply(SB101$FitRatios, 1, mean),
        apply(SB1$FitRatios, 1, mean),
        c(apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), apply(SB112$FitRatios, 1, mean), 
          apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean)), 
        c(apply(SB9$FitRatios, 1, mean), apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), apply(SB12$FitRatios, 1, mean), 
          apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean)), 
        c(apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), apply(SB105$FitRatios, 1, mean), 
          apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean)), 
        c(apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), apply(SB5$FitRatios, 1, mean), 
          apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean)), 
        
        c(apply(SB116$FitRatios, 1, mean), apply(SB117$FitRatios, 1, mean), apply(SB118$FitRatios, 1, mean)), 
        c(apply(SB16$FitRatios, 1, mean), apply(SB17$FitRatios, 1, mean), apply(SB18$FitRatios, 1, mean)), 
        c(apply(SB119$FitRatios, 1, mean), apply(SB120$FitRatios, 1, mean), apply(SB121$FitRatios, 1, mean)), 
        c(apply(SB19$FitRatios, 1, mean), apply(SB20$FitRatios, 1, mean), apply(SB21$FitRatios, 1, mean)), 
        c(apply(SB122$FitRatios, 1, mean), apply(SB123$FitRatios, 1, mean), apply(SB124$FitRatios, 1, mean)), 
        c(apply(SB22$FitRatios, 1, mean), apply(SB23$FitRatios, 1, mean), apply(SB24$FitRatios, 1, mean)), 
        c(apply(SB125$FitRatios, 1, mean), apply(SB126$FitRatios, 1, mean), apply(SB127$FitRatios, 1, mean)), 
        c(apply(SB25$FitRatios, 1, mean), apply(SB26$FitRatios, 1, mean), apply(SB27$FitRatios, 1, mean)), 
        
        c(apply(SB128$FitRatios, 1, mean), apply(SB129$FitRatios, 1, mean), apply(SB130$FitRatios, 1, mean)), 
        c(apply(SB28$FitRatios, 1, mean), apply(SB29$FitRatios, 1, mean), apply(SB30$FitRatios, 1, mean)), 
        c(apply(SB131$FitRatios, 1, mean), apply(SB132$FitRatios, 1, mean), apply(SB133$FitRatios, 1, mean)), 
        c(apply(SB31$FitRatios, 1, mean), apply(SB32$FitRatios, 1, mean), apply(SB33$FitRatios, 1, mean)), 
        c(apply(SB134$FitRatios, 1, mean), apply(SB135$FitRatios, 1, mean), apply(SB136$FitRatios, 1, mean)), 
        c(apply(SB34$FitRatios, 1, mean), apply(SB35$FitRatios, 1, mean), apply(SB36$FitRatios, 1, mean)), 
        c(apply(SB137$FitRatios, 1, mean), apply(SB138$FitRatios, 1, mean), apply(SB139$FitRatios, 1, mean)), 
        c(apply(SB37$FitRatios, 1, mean), apply(SB38$FitRatios, 1, mean), apply(SB39$FitRatios, 1, mean)), 
        
        col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
              'darkslategray4','darkslategray1','orchid4','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
              'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple'),
        # col=c('lightgrey','darkgrey','olivedrab1','olivedrab3','skyblue1','deepskyblue','darkseagreen1','darkseagreen4',
        #       'darkslategray1','darkslategray4','orchid','orchid4','darkorchid1','darkorchid4','springgreen','springgreen4',
        #       'steelblue1','steelblue4','mediumpurple1','mediumpurple4','purple','purple4'),
        names=c(
          'No
change','','I1','','D1','','I2','','D2','','I-D','','D-I','','I3','','D3','','I-D-I','','D-I-D',''), ylim=c(0,3))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)

par(mfrow=c(3,1), mar=c(1.3, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.8, 0.2, 0))
boxplot(c(SB101$biasAll.BT), 
        c(SB1$biasAll.BT), 
        c(SB109$biasAll.BT, SB110$biasAll.BT, SB111$biasAll.BT, SB112$biasAll.BT, SB113$biasAll.BT, SB114$biasAll.BT, SB115$biasAll.BT),
        c(SB9$biasAll.BT,  SB10$biasAll.BT, SB11$biasAll.BT, SB12$biasAll.BT, SB13$biasAll.BT, SB14$biasAll.BT, SB15$biasAll.BT),
        c(SB102$biasAll.BT, SB103$biasAll.BT, SB104$biasAll.BT, SB105$biasAll.BT, SB106$biasAll.BT, SB107$biasAll.BT, SB108$biasAll.BT),
        c(SB2$biasAll.BT, SB3$biasAll.BT, SB4$biasAll.BT, SB5$biasAll.BT, SB6$biasAll.BT, SB7$biasAll.BT, SB8$biasAll.BT),
        c(SB116$biasAll.BT, SB117$biasAll.BT, SB118$biasAll.BT),
        c(SB16$biasAll.BT, SB17$biasAll.BT, SB18$biasAll.BT),
        c(SB119$biasAll.BT, SB120$biasAll.BT, SB121$biasAll.BT), 
        c(SB19$biasAll.BT, SB20$biasAll.BT, SB21$biasAll.BT), 
        c(SB122$biasAll.BT, SB123$biasAll.BT, SB124$biasAll.BT), 
        c(SB22$biasAll.BT, SB23$biasAll.BT, SB24$biasAll.BT), 
        c(SB125$biasAll.BT, SB126$biasAll.BT, SB127$biasAll.BT),
        c(SB25$biasAll.BT, SB26$biasAll.BT, SB27$biasAll.BT),
        
        c(SB128$biasAll.BT, SB129$biasAll.BT, SB130$biasAll.BT),
        c(SB28$biasAll.BT, SB28$biasAll.BT, SB30$biasAll.BT),
        c(SB131$biasAll.BT, SB132$biasAll.BT, SB133$biasAll.BT), 
        c(SB131$biasAll.BT, SB32$biasAll.BT, SB33$biasAll.BT), 
        c(SB134$biasAll.BT, SB135$biasAll.BT, SB136$biasAll.BT), 
        c(SB34$biasAll.BT, SB35$biasAll.BT, SB36$biasAll.BT), 
        c(SB137$biasAll.BT, SB138$biasAll.BT, SB139$biasAll.BT),
        c(SB37$biasAll.BT, SB38$biasAll.BT, SB39$biasAll.BT),
        col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
              'darkslategray4','darkslategray1','mediumorchid1','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
              'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple'),
        # col=c('lightgrey','darkgrey','olivedrab1','olivedrab3','skyblue1','deepskyblue','darkseagreen1','darkseagreen4',
        #       'darkslategray1','darkslategray4','orchid','mediumorchid1','darkorchid1','darkorchid4','springgreen','springgreen4',
        #       'steelblue1','steelblue4','mediumpurple1','mediumpurple4','purple','purple4'),
        names=FALSE, ylim=c(-5,3))
mtext('Annual Bias',side=2,line=1.1,cex=0.8)
title('Complete vs. Missing Data', line=-0.9)
abline(h=0, col='black')


boxplot(SB101$RMSE, 
        SB1$RMSE,
        c(SB109$RMSE, SB110$RMSE, SB111$RMSE, SB112$RMSE, SB113$RMSE, SB114$RMSE, SB115$RMSE),
        c(SB9$RMSE, SB10$RMSE, SB11$RMSE, SB12$RMSE, SB13$RMSE, SB14$RMSE, SB15$RMSE),
        c(SB102$RMSE, SB103$RMSE, SB104$RMSE, SB105$RMSE, SB106$RMSE, SB107$RMSE, SB108$RMSE),
        c(SB2$RMSE, SB3$RMSE, SB4$RMSE, SB5$RMSE, SB6$RMSE, SB7$RMSE, SB8$RMSE), 
        
        c(SB116$RMSE, SB117$RMSE, SB118$RMSE),
        c(SB16$RMSE, SB17$RMSE, SB18$RMSE),
        c(SB119$RMSE, SB120$RMSE, SB121$RMSE), 
        c(SB19$RMSE, SB20$RMSE, SB21$RMSE),
        c(SB122$RMSE, SB123$RMSE, SB124$RMSE), 
        c(SB22$RMSE, SB23$RMSE, SB24$RMSE),
        c(SB125$RMSE, SB126$RMSE, SB127$RMSE),
        c(SB25$RMSE, SB26$RMSE, SB27$RMSE),
        
        c(SB128$RMSE, SB129$RMSE, SB130$RMSE),
        c(SB28$RMSE, SB29$RMSE, SB30$RMSE),
        c(SB131$RMSE, SB132$RMSE, SB133$RMSE), 
        c(SB31$RMSE, SB32$RMSE, SB33$RMSE),
        c(SB134$RMSE, SB135$RMSE, SB136$RMSE), 
        c(SB34$RMSE, SB35$RMSE, SB36$RMSE),
        c(SB137$RMSE, SB138$RMSE, SB139$RMSE),
        c(SB37$RMSE, SB38$RMSE, SB39$RMSE),
        col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
              'darkslategray4','darkslategray1','mediumorchid1','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
              'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple'),
        # col=c('lightgrey','darkgrey','olivedrab1','olivedrab3','skyblue1','deepskyblue','darkseagreen1','darkseagreen4',
        #       'darkslategray1','darkslategray4','orchid','mediumorchid1','darkorchid1','darkorchid4','springgreen','springgreen4',
        #       'steelblue1','steelblue4','mediumpurple1','mediumpurple4','purple','purple4'),
        names=FALSE, ylim=c(0,2))
mtext('RMSE',side=2,line=1.1,cex=0.8)


boxplot(apply(SB101$FitRatios, 1, mean),
        apply(SB1$FitRatios, 1, mean),
        c(apply(SB109$FitRatios, 1, mean), apply(SB110$FitRatios, 1, mean), apply(SB111$FitRatios, 1, mean), apply(SB112$FitRatios, 1, mean), 
          apply(SB113$FitRatios, 1, mean), apply(SB114$FitRatios, 1, mean), apply(SB115$FitRatios, 1, mean)), 
        c(apply(SB9$FitRatios, 1, mean), apply(SB10$FitRatios, 1, mean), apply(SB11$FitRatios, 1, mean), apply(SB12$FitRatios, 1, mean), 
          apply(SB13$FitRatios, 1, mean), apply(SB14$FitRatios, 1, mean), apply(SB15$FitRatios, 1, mean)), 
        c(apply(SB102$FitRatios, 1, mean), apply(SB103$FitRatios, 1, mean), apply(SB104$FitRatios, 1, mean), apply(SB105$FitRatios, 1, mean), 
          apply(SB106$FitRatios, 1, mean), apply(SB107$FitRatios, 1, mean), apply(SB108$FitRatios, 1, mean)), 
        c(apply(SB2$FitRatios, 1, mean), apply(SB3$FitRatios, 1, mean), apply(SB4$FitRatios, 1, mean), apply(SB5$FitRatios, 1, mean), 
          apply(SB6$FitRatios, 1, mean), apply(SB7$FitRatios, 1, mean), apply(SB8$FitRatios, 1, mean)), 
        
        c(apply(SB116$FitRatios, 1, mean), apply(SB117$FitRatios, 1, mean), apply(SB118$FitRatios, 1, mean)), 
        c(apply(SB16$FitRatios, 1, mean), apply(SB17$FitRatios, 1, mean), apply(SB18$FitRatios, 1, mean)), 
        c(apply(SB119$FitRatios, 1, mean), apply(SB120$FitRatios, 1, mean), apply(SB121$FitRatios, 1, mean)), 
        c(apply(SB19$FitRatios, 1, mean), apply(SB20$FitRatios, 1, mean), apply(SB21$FitRatios, 1, mean)), 
        c(apply(SB122$FitRatios, 1, mean), apply(SB123$FitRatios, 1, mean), apply(SB124$FitRatios, 1, mean)), 
        c(apply(SB22$FitRatios, 1, mean), apply(SB23$FitRatios, 1, mean), apply(SB24$FitRatios, 1, mean)), 
        c(apply(SB125$FitRatios, 1, mean), apply(SB126$FitRatios, 1, mean), apply(SB127$FitRatios, 1, mean)), 
        c(apply(SB25$FitRatios, 1, mean), apply(SB26$FitRatios, 1, mean), apply(SB27$FitRatios, 1, mean)), 
        
        c(apply(SB128$FitRatios, 1, mean), apply(SB129$FitRatios, 1, mean), apply(SB130$FitRatios, 1, mean)), 
        c(apply(SB28$FitRatios, 1, mean), apply(SB29$FitRatios, 1, mean), apply(SB30$FitRatios, 1, mean)), 
        c(apply(SB131$FitRatios, 1, mean), apply(SB132$FitRatios, 1, mean), apply(SB133$FitRatios, 1, mean)), 
        c(apply(SB31$FitRatios, 1, mean), apply(SB32$FitRatios, 1, mean), apply(SB33$FitRatios, 1, mean)), 
        c(apply(SB134$FitRatios, 1, mean), apply(SB135$FitRatios, 1, mean), apply(SB136$FitRatios, 1, mean)), 
        c(apply(SB34$FitRatios, 1, mean), apply(SB35$FitRatios, 1, mean), apply(SB36$FitRatios, 1, mean)), 
        c(apply(SB137$FitRatios, 1, mean), apply(SB138$FitRatios, 1, mean), apply(SB139$FitRatios, 1, mean)), 
        c(apply(SB37$FitRatios, 1, mean), apply(SB38$FitRatios, 1, mean), apply(SB39$FitRatios, 1, mean)), 
        
        col=c('darkgrey','lightgrey','olivedrab3','olivedrab1','deepskyblue','skyblue1','darkseagreen4', 'darkseagreen1',
              'darkslategray4','darkslategray1','mediumorchid1','orchid','darkorchid4','darkorchid1','springgreen4','springgreen',
              'steelblue4','steelblue1','mediumpurple4','mediumpurple1','purple4','purple'),
        # col=c('lightgrey','darkgrey','olivedrab1','olivedrab3','skyblue1','deepskyblue','darkseagreen1','darkseagreen4',
        #       'darkslategray1','darkslategray4','orchid','mediumorchid1','darkorchid1','darkorchid4','springgreen','springgreen4',
        #       'steelblue1','steelblue4','mediumpurple1','mediumpurple4','purple','purple4'),
        names=c(
          'No
          change','','I1','','D1','','I2','','D2','','I-D','','D-I','','I3','','D3','','I-D-I','','D-I-D',''), ylim=c(0,2.5))
mtext('Mean Fit Ratio',side=2,line=1.1,cex=0.8)






############### RMSE vs. Fit Ratio #############
col=c('firebrick','firebrick1', 'deepskyblue4','deepskyblue','olivedrab','olivedrab1','darkorchid4','darkorchid')
load( file="D:\\vspace1\\DFA_Simulation/Atl_SN/NEW_ANALYSES_2019/BigSNSimulationResults.RData")



png(filename="D:\\vspace1\\DFA_Simulation\\SB\\Plots\\FRvsRMSE_BothSpecies.png", 
    type="cairo",
    units="mm", 
    width=250, 
    height=200, 
    pointsize=16, 
    res=600)
# STANDARDIZE AXES
par(mfrow=c(1,1), mar=c(2.1, 2.1, 0.1, 0.3),tcl = -0.1, mgp = c(0.9, 0.2, 0))
plot(SB101$RMSE, apply(SB101$FitRatios,1,mean),  ylab="Fit Ratios", xlab="RMSE", xlim=c(0,2), ylim=c(0,1), col='firebrick', pch=1)
points(SB1$RMSE, apply(SB1$FitRatios,1,mean), pch=0, col='firebrick1')

for(i in c(1:6)){
  points(BigSNSimulationResults[[i]]$RMSE, apply( BigSNSimulationResults[[i]]$FitRatios,1,mean), col='orange', pch=2)
  j=i+24
  points(BigSNSimulationResults[[j]]$RMSE, apply( BigSNSimulationResults[[j]]$FitRatios,1,mean), col='orange', pch=6)
}
for(i in c(49:51)){
  points(BigSNSimulationResults[[i]]$RMSE, apply( BigSNSimulationResults[[i]]$FitRatios,1,mean), col='orange', pch=2)
  j=i+12
  points(BigSNSimulationResults[[j]]$RMSE, apply( BigSNSimulationResults[[j]]$FitRatios,1,mean), col='orange', pch=6)
}

for(i in c(7:12)){
  points(BigSNSimulationResults[[i]]$RMSE, apply( BigSNSimulationResults[[i]]$FitRatios,1,mean), col='blue', pch=2)
  j=i+24
  points(BigSNSimulationResults[[j]]$RMSE, apply( BigSNSimulationResults[[j]]$FitRatios,1,mean), col='blue', pch=6)
}
for(i in c(52:54)){
  points(BigSNSimulationResults[[i]]$RMSE, apply( BigSNSimulationResults[[i]]$FitRatios,1,mean), col='blue', pch=2)
  j=i+12
  points(BigSNSimulationResults[[j]]$RMSE, apply( BigSNSimulationResults[[j]]$FitRatios,1,mean), col='blue', pch=6)
}

for(i in c(13:18)){
  points(BigSNSimulationResults[[i]]$RMSE, apply( BigSNSimulationResults[[i]]$FitRatios,1,mean), col='forestgreen', pch=2)
  j=i+24
  points(BigSNSimulationResults[[j]]$RMSE, apply( BigSNSimulationResults[[j]]$FitRatios,1,mean), col='forestgreen', pch=6)
}
for(i in c(55:57)){
  points(BigSNSimulationResults[[i]]$RMSE, apply( BigSNSimulationResults[[i]]$FitRatios,1,mean), col='forestgreen', pch=2)
  j=i+12
  points(BigSNSimulationResults[[j]]$RMSE, apply( BigSNSimulationResults[[j]]$FitRatios,1,mean), col='forestgreen', pch=6)
}


for(i in c(19:22)){
  points(BigSNSimulationResults[[i]]$RMSE, apply( BigSNSimulationResults[[i]]$FitRatios,1,mean), col='deeppink', pch=2)
  j=i+24
  points(BigSNSimulationResults[[j]]$RMSE, apply( BigSNSimulationResults[[j]]$FitRatios,1,mean), col='deeppink', pch=6)
}
for(i in c(58:60)){
  points(BigSNSimulationResults[[i]]$RMSE, apply( BigSNSimulationResults[[i]]$FitRatios,1,mean), col='deeppink', pch=2)
  j=i+12
  points(BigSNSimulationResults[[j]]$RMSE, apply( BigSNSimulationResults[[j]]$FitRatios,1,mean), col='deeppink', pch=6)
}


for(i in 2:15){
  j=i+100
  points(get(paste('SB',j,sep=""))$RMSE, apply( get(paste('SB',j,sep=''))$FitRatios ,1,mean), col="deepskyblue4", pch=1)
  points(get(paste('SB',i,sep=""))$RMSE, apply( get(paste('SB',i,sep=''))$FitRatios ,1,mean), col="deepskyblue", pch=0)
}
for(i in 16:27){
  j=i+100
  points(get(paste('SB',j,sep=""))$RMSE, apply( get(paste('SB',j,sep=''))$FitRatios ,1,mean), col="olivedrab", pch=1)
  points(get(paste('SB',i,sep=""))$RMSE, apply( get(paste('SB',i,sep=''))$FitRatios ,1,mean), col="olivedrab1", pch=0)
}
for(i in 28:39){
  j=i+100
  points(get(paste('SB',j,sep=""))$RMSE, apply( get(paste('SB',j,sep=''))$FitRatios ,1,mean), col="darkorchid4", pch=1)
  points(get(paste('SB',i,sep=""))$RMSE, apply( get(paste('SB',i,sep=''))$FitRatios ,1,mean), col="darkorchid", pch=0)
}
points(SB101$RMSE, apply(SB101$FitRatios,1,mean), pch=0, col='firebrick')
points(SB1$RMSE, apply(SB1$FitRatios,1,mean), pch=0, col='firebrick1')

dev.off()


abline(v=0.5)
abline(h=0.5)
abline(h=0.6, col='darkgrey')
abline(h=0.55, col='darkgrey')




FRs = cbind("Fit Ratio" = apply(SB1$FitRatios,1,mean), "RMSE" = SB1$RMSE)
FRs = rbind(FRs,  cbind(apply(SB101$FitRatios,1,mean),SB101$RMSE))
for(i in 2:39){
  j=i+100
  FRs = rbind(FRs, cbind(apply(get(paste("SB",i,sep=""))$FitRatios,1,mean), get(paste("SB",i,sep=""))$RMSE) )
  FRs = rbind(FRs, cbind(apply(get(paste("SB",j,sep=""))$FitRatios,1,mean), get(paste("SB",j,sep=""))$RMSE) )
}
for(i in 1:72){
  FRs = rbind(FRs, cbind(apply(BigSNSimulationResults[[i]]$FitRatios,1,mean), BigSNSimulationResults[[1]]$RMSE) )
}

write.csv(FRs, "D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Compare_MissData\\FRs2.csv")
write.csv(FRs, "D:\\vspace1\\DFA_Simulation\\SB\\Plots\\Compare_MissData\\FRs.csv")
summary(FRs)

summary(apply(SB29$FitRatios,1,mean))

######################## NOTES #################################
names(SB1)







####################### PLOT ########################################
names(Trial24)
par(mfrow=c(4,1), mar=c(0.6, 1.1, 0.6, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
vioplot(ATrial1$RMSE, ATrial2$RMSE, ATrial3$RMSE, ATrial4$RMSE, ATrial5$RMSE, ATrial6$RMSE, ATrial7$RMSE, ATrial8$RMSE, 
        ATrial9$RMSE, ATrial10$RMSE, ATrial11$RMSE, ATrial12$RMSE, ATrial13$RMSE, ATrial14$RMSE, ATrial15$RMSE, ATrial16$RMSE, 
        ATrial17$RMSE, ATrial18$RMSE, ATrial19$RMSE, ATrial20$RMSE, ATrial21$RMSE, ATrial22$RMSE, ATrial23$RMSE, ATrial24$RMSE,
        ylim=c(0,2), col='deepskyblue')
abline(h=1)
vioplot(BTrial1$RMSE, BTrial2$RMSE, BTrial3$RMSE, BTrial4$RMSE, BTrial5$RMSE, BTrial6$RMSE, BTrial7$RMSE, BTrial8$RMSE, 
        BTrial9$RMSE, BTrial10$RMSE, BTrial11$RMSE, BTrial12$RMSE, BTrial13$RMSE, BTrial14$RMSE, BTrial15$RMSE, BTrial16$RMSE, 
        BTrial17$RMSE, BTrial18$RMSE, BTrial19$RMSE, BTrial20$RMSE, BTrial21$RMSE, BTrial22$RMSE, BTrial23$RMSE, BTrial24$RMSE,
        ylim=c(0,2), col='darkolivegreen1')
abline(h=1)
vioplot(ATrial1$RMSE, ATrial2$RMSE, ATrial3$RMSE, CTrial4$RMSE, CTrial5$RMSE, CTrial6$RMSE, ATrial7$RMSE, ATrial8$RMSE, 
        ATrial9$RMSE, CTrial10$RMSE, CTrial11$RMSE, CTrial12$RMSE, ATrial13$RMSE, ATrial14$RMSE, ATrial15$RMSE, CTrial16$RMSE, 
        CTrial17$RMSE, CTrial18$RMSE, ATrial19$RMSE, ATrial20$RMSE, ATrial21$RMSE, CTrial22$RMSE, CTrial23$RMSE, CTrial24$RMSE,
        ylim=c(0,2), col=c('deepskyblue4', rep('deepskyblue', 3), rep('deepskyblue4',3), rep('deepskyblue', 3), rep('deepskyblue4',3), 
                           rep('deepskyblue', 3), rep('deepskyblue4',3), rep('deepskyblue', 3), rep('deepskyblue4',3)))
abline(h=1)
vioplot(BTrial1$RMSE, BTrial2$RMSE, BTrial3$RMSE, DTrial4$RMSE, DTrial5$RMSE, DTrial6$RMSE, BTrial7$RMSE, BTrial8$RMSE, 
        BTrial9$RMSE, DTrial10$RMSE, DTrial11$RMSE, DTrial12$RMSE, BTrial13$RMSE, BTrial14$RMSE, BTrial15$RMSE, DTrial16$RMSE, 
        DTrial17$RMSE, DTrial18$RMSE, BTrial19$RMSE, BTrial20$RMSE, BTrial21$RMSE, DTrial22$RMSE, DTrial23$RMSE, DTrial24$RMSE,
        ylim=c(0,2), col=c('darkolivegreen', rep('darkolivegreen1', 3), rep('darkolivegreen',3), rep('darkolivegreen1', 3), rep('darkolivegreen',3), 
                           rep('darkolivegreen1', 3), rep('darkolivegreen',3), rep('darkolivegreen1', 3), rep('darkolivegreen',3)))
abline(h=1)

### hack vioplot to plot in multiple colors ###
library(sm)




par(mfrow=c(4,1), mar=c(0.6, 1.1, 1.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))
vioplot(ATrial1$RMSE, ATrial2$RMSE, ATrial3$RMSE, ATrial4$RMSE, ATrial5$RMSE, ATrial6$RMSE, 
        CTrial4$RMSE, CTrial5$RMSE, CTrial6$RMSE, 
        BTrial1$RMSE, BTrial2$RMSE, BTrial3$RMSE, BTrial4$RMSE, BTrial5$RMSE, BTrial6$RMSE,
        DTrial4$RMSE, DTrial5$RMSE, DTrial6$RMSE,
        ylim=c(0,2), col=c(rep('deepskyblue',4),rep('lightskyblue',3),rep('cyan',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3),rep('green',3)),
        names=FALSE)
title('RMSE: Constant F', line=0.4)
vioplot(ATrial7$RMSE, ATrial8$RMSE, ATrial9$RMSE, ATrial10$RMSE, ATrial11$RMSE, ATrial12$RMSE, 
        CTrial10$RMSE, CTrial11$RMSE, CTrial12$RMSE,  
        BTrial7$RMSE, BTrial8$RMSE, BTrial9$RMSE, BTrial10$RMSE, BTrial11$RMSE, BTrial12$RMSE, 
        DTrial10$RMSE, DTrial11$RMSE, DTrial12$RMSE,
        ylim=c(0,2), col=c(rep('deepskyblue',4),rep('lightskyblue',3),rep('cyan',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3),rep('green',3)),
        names=FALSE)
title('RMSE: Decreasing F', line=0.4)

vioplot(ATrial13$RMSE, ATrial14$RMSE, ATrial15$RMSE, ATrial16$RMSE, ATrial17$RMSE, ATrial18$RMSE, 
        CTrial16$RMSE, CTrial17$RMSE, CTrial18$RMSE, 
        BTrial13$RMSE, BTrial14$RMSE, BTrial15$RMSE, BTrial16$RMSE, BTrial17$RMSE, BTrial18$RMSE, 
        DTrial16$RMSE, DTrial17$RMSE, DTrial18$RMSE,
        ylim=c(0,2), col=c(rep('deepskyblue',4),rep('lightskyblue',3),rep('cyan',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3),rep('green',3)),
        names=FALSE)
title('RMSE: Increasing F', line=0.4)
vioplot(ATrial19$RMSE, ATrial20$RMSE, ATrial21$RMSE, ATrial22$RMSE, ATrial23$RMSE, ATrial24$RMSE, 
        CTrial22$RMSE, CTrial23$RMSE, CTrial24$RMSE,  
        BTrial19$RMSE, BTrial20$RMSE, BTrial21$RMSE, BTrial22$RMSE, BTrial23$RMSE, BTrial24$RMSE, 
        DTrial22$RMSE, DTrial23$RMSE, DTrial24$RMSE,
        ylim=c(0,2), col=c(rep('deepskyblue',4),rep('lightskyblue',3),rep('cyan',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3),rep('green',3)),
        names=FALSE)
title('RMSE: Decrease-Increasing F', line=0.4)
legend('top', c('Constant q, 3 surveys', 'Knife-edge change in q, 3 surveys', 'Gradual change in q, 3 surveys',
                'Constant q, 4 surveys', 'Knife-edge change in q, 4 surveys', 'Gradual change in q, 4 surveys'), 
       fill=c("deepskyblue",'lightskyblue','cyan','darkolivegreen3','darkolivegreen1','green'),ncol=3)


######################## Fit Ratios ##############
ATrial1$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial1\\FitRatios.csv')
ATrial2$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial2\\FitRatios.csv')
ATrial3$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial3\\FitRatios.csv')
ATrial4$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial4\\FitRatios.csv')
ATrial5$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial5\\FitRatios.csv')
ATrial6$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial6\\FitRatios.csv')
ATrial7$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial7\\FitRatios.csv')
ATrial8$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial8\\FitRatios.csv')
ATrial9$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial9\\FitRatios.csv')
ATrial10$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial10\\FitRatios.csv')
ATrial11$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial11\\FitRatios.csv')
ATrial12$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial12\\FitRatios.csv')
ATrial13$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial13\\FitRatios.csv')
ATrial14$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial14\\FitRatios.csv')
ATrial15$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial15\\FitRatios.csv')
ATrial16$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial16\\FitRatios.csv')
ATrial17$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial17\\FitRatios.csv')
ATrial18$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial18\\FitRatios.csv')
ATrial19$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial19\\FitRatios.csv')
ATrial20$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial20\\FitRatios.csv')
ATrial21$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial21\\FitRatios.csv')
ATrial22$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial22\\FitRatios.csv')
ATrial23$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial23\\FitRatios.csv')
ATrial24$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Base\\Trial24\\FitRatios.csv')

BTrial1$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial1\\FitRatios.csv')
BTrial2$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial2\\FitRatios.csv')
BTrial3$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial3\\FitRatios.csv')
BTrial4$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial4\\FitRatios.csv')
BTrial5$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial5\\FitRatios.csv')
BTrial6$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial6\\FitRatios.csv')
BTrial7$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial7\\FitRatios.csv')
BTrial8$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial8\\FitRatios.csv')
BTrial9$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial9\\FitRatios.csv')
BTrial10$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial10\\FitRatios.csv')
BTrial11$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial11\\FitRatios.csv')
BTrial12$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial12\\FitRatios.csv')
BTrial13$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial13\\FitRatios.csv')
BTrial14$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial14\\FitRatios.csv')
BTrial15$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial15\\FitRatios.csv')
BTrial16$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial16\\FitRatios.csv')
BTrial17$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial17\\FitRatios.csv')
BTrial18$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial18\\FitRatios.csv')
BTrial19$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial19\\FitRatios.csv')
BTrial20$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial20\\FitRatios.csv')
BTrial21$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial21\\FitRatios.csv')
BTrial22$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial22\\FitRatios.csv')
BTrial23$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial23\\FitRatios.csv')
BTrial24$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\FourIndices\\Trial24\\FitRatios.csv')


CTrial4$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Gradual\\Trial4\\FitRatios.csv')
CTrial5$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Gradual\\Trial5\\FitRatios.csv')
CTrial6$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Gradual\\Trial6\\FitRatios.csv')
CTrial10$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Gradual\\Trial10\\FitRatios.csv')
CTrial11$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Gradual\\Trial11\\FitRatios.csv')
CTrial12$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Gradual\\Trial12\\FitRatios.csv')
CTrial16$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Gradual\\Trial16\\FitRatios.csv')
CTrial17$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Gradual\\Trial17\\FitRatios.csv')
CTrial18$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Gradual\\Trial18\\FitRatios.csv')
CTrial22$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Gradual\\Trial22\\FitRatios.csv')
CTrial23$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Gradual\\Trial23\\FitRatios.csv')
CTrial24$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Gradual\\Trial24\\FitRatios.csv')

DTrial4$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndicesGradual\\Trial4\\FitRatios.csv')
DTrial5$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndicesGradual\\Trial5\\FitRatios.csv')
DTrial6$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndicesGradual\\Trial6\\FitRatios.csv')
DTrial10$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndicesGradual\\Trial10\\FitRatios.csv')
DTrial11$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndicesGradual\\Trial11\\FitRatios.csv')
DTrial12$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndicesGradual\\Trial12\\FitRatios.csv')
DTrial16$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndicesGradual\\Trial16\\FitRatios.csv')
DTrial17$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndicesGradual\\Trial17\\FitRatios.csv')
DTrial18$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndicesGradual\\Trial18\\FitRatios.csv')
DTrial22$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndicesGradual\\Trial22\\FitRatios.csv')
DTrial23$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndicesGradual\\Trial23\\FitRatios.csv')
DTrial24$FitRatios = read.csv('N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndicesGradual\\Trial24\\FitRatios.csv')



ATrial1$AbsBias = apply(abs(ATrial1$biasAll),1,sum)
ATrial2$AbsBias = apply(abs(ATrial2$biasAll),1,sum)
ATrial3$AbsBias = apply(abs(ATrial3$biasAll),1,sum)
ATrial4$AbsBias = apply(abs(ATrial4$biasAll),1,sum)
ATrial5$AbsBias = apply(abs(ATrial5$biasAll),1,sum)
ATrial6$AbsBias = apply(abs(ATrial6$biasAll),1,sum)
ATrial7$AbsBias = apply(abs(ATrial7$biasAll),1,sum)
ATrial8$AbsBias = apply(abs(ATrial8$biasAll),1,sum)
ATrial9$AbsBias = apply(abs(ATrial9$biasAll),1,sum)
ATrial10$AbsBias = apply(abs(ATrial10$biasAll),1,sum)
ATrial11$AbsBias = apply(abs(ATrial11$biasAll),1,sum)
ATrial12$AbsBias = apply(abs(ATrial12$biasAll),1,sum)
ATrial13$AbsBias = apply(abs(ATrial13$biasAll),1,sum)
ATrial14$AbsBias = apply(abs(ATrial14$biasAll),1,sum)
ATrial15$AbsBias = apply(abs(ATrial15$biasAll),1,sum)
ATrial16$AbsBias = apply(abs(ATrial16$biasAll),1,sum)
ATrial17$AbsBias = apply(abs(ATrial17$biasAll),1,sum)
ATrial18$AbsBias = apply(abs(ATrial18$biasAll),1,sum)
ATrial19$AbsBias = apply(abs(ATrial19$biasAll),1,sum)
ATrial20$AbsBias = apply(abs(ATrial20$biasAll),1,sum)
ATrial21$AbsBias = apply(abs(ATrial21$biasAll),1,sum)
ATrial22$AbsBias = apply(abs(ATrial22$biasAll),1,sum)
ATrial23$AbsBias = apply(abs(ATrial23$biasAll),1,sum)
ATrial24$AbsBias = apply(abs(ATrial24$biasAll),1,sum)


BTrial1$AbsBias = apply(abs(BTrial1$biasAll),1,sum)
BTrial2$AbsBias = apply(abs(BTrial2$biasAll),1,sum)
BTrial3$AbsBias = apply(abs(BTrial3$biasAll),1,sum)
BTrial4$AbsBias = apply(abs(BTrial4$biasAll),1,sum)
BTrial5$AbsBias = apply(abs(BTrial5$biasAll),1,sum)
BTrial6$AbsBias = apply(abs(BTrial6$biasAll),1,sum)
BTrial7$AbsBias = apply(abs(BTrial7$biasAll),1,sum)
BTrial8$AbsBias = apply(abs(BTrial8$biasAll),1,sum)
BTrial9$AbsBias = apply(abs(BTrial9$biasAll),1,sum)
BTrial10$AbsBias = apply(abs(BTrial10$biasAll),1,sum)
BTrial11$AbsBias = apply(abs(BTrial11$biasAll),1,sum)
BTrial12$AbsBias = apply(abs(BTrial12$biasAll),1,sum)
BTrial13$AbsBias = apply(abs(BTrial13$biasAll),1,sum)
BTrial14$AbsBias = apply(abs(BTrial14$biasAll),1,sum)
BTrial15$AbsBias = apply(abs(BTrial15$biasAll),1,sum)
BTrial16$AbsBias = apply(abs(BTrial16$biasAll),1,sum)
BTrial17$AbsBias = apply(abs(BTrial17$biasAll),1,sum)
BTrial18$AbsBias = apply(abs(BTrial18$biasAll),1,sum)
BTrial19$AbsBias = apply(abs(BTrial19$biasAll),1,sum)
BTrial20$AbsBias = apply(abs(BTrial20$biasAll),1,sum)
BTrial21$AbsBias = apply(abs(BTrial21$biasAll),1,sum)
BTrial22$AbsBias = apply(abs(BTrial22$biasAll),1,sum)
BTrial23$AbsBias = apply(abs(BTrial23$biasAll),1,sum)
BTrial24$AbsBias = apply(abs(BTrial24$biasAll),1,sum)

# CTrial1$AbsBias = apply(abs(CTrial1$biasAll),1,sum)
# CTrial2$AbsBias = apply(abs(CTrial2$biasAll),1,sum)
# CTrial3$AbsBias = apply(abs(CTrial3$biasAll),1,sum)
CTrial4$AbsBias = apply(abs(CTrial4$biasAll),1,sum)
CTrial5$AbsBias = apply(abs(CTrial5$biasAll),1,sum)
CTrial6$AbsBias = apply(abs(CTrial6$biasAll),1,sum)
# CTrial7$AbsBias = apply(abs(CTrial7$biasAll),1,sum)
# CTrial8$AbsBias = apply(abs(CTrial8$biasAll),1,sum)
# CTrial9$AbsBias = apply(abs(CTrial9$biasAll),1,sum)
CTrial10$AbsBias = apply(abs(CTrial10$biasAll),1,sum)
CTrial11$AbsBias = apply(abs(CTrial11$biasAll),1,sum)
CTrial12$AbsBias = apply(abs(CTrial12$biasAll),1,sum)
# CTrial13$AbsBias = apply(abs(CTrial13$biasAll),1,sum)
# CTrial14$AbsBias = apply(abs(CTrial14$biasAll),1,sum)
# CTrial15$AbsBias = apply(abs(CTrial15$biasAll),1,sum)
CTrial16$AbsBias = apply(abs(CTrial16$biasAll),1,sum)
CTrial17$AbsBias = apply(abs(CTrial17$biasAll),1,sum)
CTrial18$AbsBias = apply(abs(CTrial18$biasAll),1,sum)
# CTrial19$AbsBias = apply(abs(CTrial19$biasAll),1,sum)
# CTrial20$AbsBias = apply(abs(CTrial20$biasAll),1,sum)
# CTrial21$AbsBias = apply(abs(CTrial21$biasAll),1,sum)
CTrial22$AbsBias = apply(abs(CTrial22$biasAll),1,sum)
CTrial23$AbsBias = apply(abs(CTrial23$biasAll),1,sum)
CTrial24$AbsBias = apply(abs(CTrial24$biasAll),1,sum)



# DTrial1$AbsBias = apply(abs(DTrial1$biasAll),1,sum)
# DTrial2$AbsBias = apply(abs(DTrial2$biasAll),1,sum)
# DTrial3$AbsBias = apply(abs(DTrial3$biasAll),1,sum)
DTrial4$AbsBias = apply(abs(DTrial4$biasAll),1,sum)
DTrial5$AbsBias = apply(abs(DTrial5$biasAll),1,sum)
DTrial6$AbsBias = apply(abs(DTrial6$biasAll),1,sum)
# DTrial7$AbsBias = apply(abs(DTrial7$biasAll),1,sum)
# DTrial8$AbsBias = apply(abs(DTrial8$biasAll),1,sum)
# DTrial9$AbsBias = apply(abs(DTrial9$biasAll),1,sum)
DTrial10$AbsBias = apply(abs(DTrial10$biasAll),1,sum)
DTrial11$AbsBias = apply(abs(DTrial11$biasAll),1,sum)
DTrial12$AbsBias = apply(abs(DTrial12$biasAll),1,sum)
# DTrial13$AbsBias = apply(abs(DTrial13$biasAll),1,sum)
# DTrial14$AbsBias = apply(abs(DTrial14$biasAll),1,sum)
# DTrial15$AbsBias = apply(abs(DTrial15$biasAll),1,sum)
DTrial16$AbsBias = apply(abs(DTrial16$biasAll),1,sum)
DTrial17$AbsBias = apply(abs(DTrial17$biasAll),1,sum)
DTrial18$AbsBias = apply(abs(DTrial18$biasAll),1,sum)
# DTrial19$AbsBias = apply(abs(DTrial19$biasAll),1,sum)
# DTrial20$AbsBias = apply(abs(DTrial20$biasAll),1,sum)
# DTrial21$AbsBias = apply(abs(DTrial21$biasAll),1,sum)
DTrial22$AbsBias = apply(abs(DTrial22$biasAll),1,sum)
DTrial23$AbsBias = apply(abs(DTrial23$biasAll),1,sum)
DTrial24$AbsBias = apply(abs(DTrial24$biasAll),1,sum)

par(mfrow=c(4,1), mar=c(0.6, 1.1, 1.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))

vioplot(ATrial1$AbsBias, ATrial2$AbsBias, ATrial3$AbsBias, ATrial4$AbsBias, ATrial5$AbsBias, ATrial6$AbsBias, 
        CTrial4$AbsBias, CTrial5$AbsBias, CTrial6$AbsBias, 
        BTrial1$AbsBias, BTrial2$AbsBias, BTrial3$AbsBias, BTrial4$AbsBias, BTrial5$AbsBias, BTrial6$AbsBias,
        DTrial4$AbsBias, DTrial5$AbsBias, DTrial6$AbsBias,
        ylim=c(0,50),
        col=c(rep('deepskyblue',4),rep('lightskyblue',3),rep('cyan',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3),rep('green',3)),
        names=FALSE)
title('No Trend', line=-1)
vioplot(ATrial7$AbsBias, ATrial8$AbsBias, ATrial9$AbsBias, ATrial10$AbsBias, ATrial11$AbsBias, ATrial12$AbsBias, 
        CTrial10$AbsBias, CTrial11$AbsBias, CTrial12$AbsBias,  
        BTrial7$AbsBias, BTrial8$AbsBias, BTrial9$AbsBias, BTrial10$AbsBias, BTrial11$AbsBias, BTrial12$AbsBias, 
        DTrial10$AbsBias, DTrial11$AbsBias, DTrial12$AbsBias,
        ylim=c(0,10),
        col=c(rep('deepskyblue',4),rep('lightskyblue',3),rep('cyan',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3),rep('green',3)),
        names=FALSE)
title('Decreasing Trend', line=-1)

vioplot(ATrial13$AbsBias, ATrial14$AbsBias, ATrial15$AbsBias, ATrial16$AbsBias, ATrial17$AbsBias, ATrial18$AbsBias, 
        CTrial16$AbsBias, CTrial17$AbsBias, CTrial18$AbsBias, 
        BTrial13$AbsBias, BTrial14$AbsBias, BTrial15$AbsBias, BTrial16$AbsBias, BTrial17$AbsBias, BTrial18$AbsBias, 
        DTrial16$AbsBias, DTrial17$AbsBias, DTrial18$AbsBias,
        ylim=c(0,10), 
        col=c(rep('deepskyblue',4),rep('lightskyblue',3),rep('cyan',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3),rep('green',3)),
        names=FALSE)
title('Increasing Trend', line=-1)

vioplot(ATrial19$AbsBias, ATrial20$AbsBias, ATrial21$AbsBias, ATrial22$AbsBias, ATrial23$AbsBias, ATrial24$AbsBias, 
        CTrial22$AbsBias, CTrial23$AbsBias, CTrial24$AbsBias,  
        BTrial19$AbsBias, BTrial20$AbsBias, BTrial21$AbsBias, BTrial22$AbsBias, BTrial23$AbsBias, BTrial24$AbsBias, 
        DTrial22$AbsBias, DTrial23$AbsBias, DTrial24$AbsBias,
        ylim=c(0,40),
        col=c(rep('deepskyblue',4),rep('lightskyblue',3),rep('cyan',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3),rep('green',3)),
        names=FALSE)
title('Decrease - Increasing Trend', line=-1)







legend('top', c('Constant q, 3 surveys', 'Knife-edge change in q, 3 surveys', 'Gradual change in q, 3 surveys',
                'Constant q, 4 surveys', 'Knife-edge change in q, 4 surveys', 'Gradual change in q, 4 surveys'), 
       fill=c("deepskyblue",'lightskyblue','cyan','darkolivegreen3','darkolivegreen1','green'),ncol=3)








tiff(filename="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Plots\\AllBiasRMSE_FitRatios.tiff", 
     type="cairo",
     units="mm", 
     width=500, 
     height=167, 
     pointsize=16, 
     res=600)
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
# title('Annual Bias', line=-1.4)

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
        ylim=c(-5,5), col=c(rep('deepskyblue',2),'lightskyblue','darkolivegreen3','darkolivegreen1'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
title('Annual Bias', line=-1.4)
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
        ylim=c(0,2), col=c(rep('deepskyblue',2),'lightskyblue','darkolivegreen3','darkolivegreen1'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
title('RMSE', line=-1.4)



vioplot(c(apply(ATrial1$FitRatios[,2:4], 1, FUN=mean), apply(ATrial2$FitRatios[,2:4], 1, FUN=mean), 
          apply(ATrial3$FitRatios[,2:4], 1, FUN=mean), apply(ATrial4$FitRatios[,2:4], 1, FUN=mean),
          apply(ATrial5$FitRatios[,2:4], 1, FUN=mean), apply(ATrial6$FitRatios[,2:4], 1, FUN=mean), 
          apply(CTrial4$FitRatios[,2:4], 1, FUN=mean), apply(CTrial5$FitRatios[,2:4], 1, FUN=mean),
          apply(CTrial6$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial1$FitRatios[,2:4], 1, FUN=mean), apply(BTrial2$FitRatios[,2:4], 1, FUN=mean), 
          apply(BTrial3$FitRatios[,2:4], 1, FUN=mean), apply(BTrial4$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial5$FitRatios[,2:4], 1, FUN=mean), apply(BTrial6$FitRatios[,2:4], 1, FUN=mean), 
          apply(DTrial4$FitRatios[,2:4], 1, FUN=mean), apply(CTrial5$FitRatios[,2:4], 1, FUN=mean),
          apply(DTrial6$FitRatios[,2:4], 1, FUN=mean)),
        
        c(apply(ATrial7$FitRatios[,2:4], 1, FUN=mean), apply(ATrial8$FitRatios[,2:4], 1, FUN=mean), 
          apply(ATrial9$FitRatios[,2:4], 1, FUN=mean), apply(ATrial10$FitRatios[,2:4], 1, FUN=mean),
          apply(ATrial11$FitRatios[,2:4], 1, FUN=mean), apply(ATrial12$FitRatios[,2:4], 1, FUN=mean), 
          apply(CTrial10$FitRatios[,2:4], 1, FUN=mean), apply(CTrial11$FitRatios[,2:4], 1, FUN=mean),
          apply(CTrial12$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial7$FitRatios[,2:4], 1, FUN=mean), apply(BTrial8$FitRatios[,2:4], 1, FUN=mean), 
          apply(BTrial9$FitRatios[,2:4], 1, FUN=mean), apply(BTrial10$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial11$FitRatios[,2:4], 1, FUN=mean), apply(BTrial12$FitRatios[,2:4], 1, FUN=mean), 
          apply(DTrial10$FitRatios[,2:4], 1, FUN=mean), apply(CTrial11$FitRatios[,2:4], 1, FUN=mean),
          apply(DTrial12$FitRatios[,2:4], 1, FUN=mean)),
        
        c(apply(ATrial13$FitRatios[,2:4], 1, FUN=mean), apply(ATrial14$FitRatios[,2:4], 1, FUN=mean), 
          apply(ATrial15$FitRatios[,2:4], 1, FUN=mean), apply(ATrial16$FitRatios[,2:4], 1, FUN=mean),
          apply(ATrial17$FitRatios[,2:4], 1, FUN=mean), apply(ATrial18$FitRatios[,2:4], 1, FUN=mean), 
          apply(CTrial16$FitRatios[,2:4], 1, FUN=mean), apply(CTrial17$FitRatios[,2:4], 1, FUN=mean),
          apply(CTrial18$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial13$FitRatios[,2:4], 1, FUN=mean), apply(BTrial14$FitRatios[,2:4], 1, FUN=mean), 
          apply(BTrial15$FitRatios[,2:4], 1, FUN=mean), apply(BTrial16$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial17$FitRatios[,2:4], 1, FUN=mean), apply(BTrial18$FitRatios[,2:4], 1, FUN=mean), 
          apply(DTrial16$FitRatios[,2:4], 1, FUN=mean), apply(CTrial17$FitRatios[,2:4], 1, FUN=mean),
          apply(DTrial18$FitRatios[,2:4], 1, FUN=mean)),
        
        c(apply(ATrial19$FitRatios[,2:4], 1, FUN=mean), apply(ATrial20$FitRatios[,2:4], 1, FUN=mean), 
          apply(ATrial21$FitRatios[,2:4], 1, FUN=mean), apply(ATrial22$FitRatios[,2:4], 1, FUN=mean),
          apply(ATrial23$FitRatios[,2:4], 1, FUN=mean), apply(ATrial24$FitRatios[,2:4], 1, FUN=mean), 
          apply(CTrial22$FitRatios[,2:4], 1, FUN=mean), apply(CTrial23$FitRatios[,2:4], 1, FUN=mean),
          apply(CTrial24$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial19$FitRatios[,2:4], 1, FUN=mean), apply(BTrial20$FitRatios[,2:4], 1, FUN=mean), 
          apply(BTrial21$FitRatios[,2:4], 1, FUN=mean), apply(BTrial22$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial23$FitRatios[,2:4], 1, FUN=mean), apply(BTrial24$FitRatios[,2:4], 1, FUN=mean), 
          apply(DTrial22$FitRatios[,2:4], 1, FUN=mean), apply(CTrial23$FitRatios[,2:4], 1, FUN=mean),
          apply(DTrial24$FitRatios[,2:4], 1, FUN=mean)),
        
        ylim=c(0,1), col=c(rep('deepskyblue',2),'lightskyblue','darkolivegreen3','darkolivegreen1'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
title('Fit Ratio', line=-1.4)

dev.off()

par(mfrow=c(1,1))
vioplot(apply(ATrial1$FitRatios[,2:4], 1, FUN=mean), ylim=c(0,1), col='deepskyblue', names=c('Trial Base 1'))



tiff(filename="N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Plots\\FitRatios.tiff", 
     type="cairo",
     units="mm", 
     width=250, 
     height=250, 
     pointsize=16, 
     res=600)
par(mfrow=c(1,1), mar=c(1.6, 1.6, 1.1, 0.1),tcl = -0.1, mgp = c(2, 0.6, 0))


vioplot(c(apply(ATrial1$FitRatios[,2:4], 1, FUN=mean), apply(ATrial2$FitRatios[,2:4], 1, FUN=mean), 
          apply(ATrial3$FitRatios[,2:4], 1, FUN=mean), apply(ATrial4$FitRatios[,2:4], 1, FUN=mean),
          apply(ATrial5$FitRatios[,2:4], 1, FUN=mean), apply(ATrial6$FitRatios[,2:4], 1, FUN=mean), 
          apply(CTrial4$FitRatios[,2:4], 1, FUN=mean), apply(CTrial5$FitRatios[,2:4], 1, FUN=mean),
          apply(CTrial6$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial1$FitRatios[,2:4], 1, FUN=mean), apply(BTrial2$FitRatios[,2:4], 1, FUN=mean), 
          apply(BTrial3$FitRatios[,2:4], 1, FUN=mean), apply(BTrial4$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial5$FitRatios[,2:4], 1, FUN=mean), apply(BTrial6$FitRatios[,2:4], 1, FUN=mean), 
          apply(DTrial4$FitRatios[,2:4], 1, FUN=mean), apply(CTrial5$FitRatios[,2:4], 1, FUN=mean),
          apply(DTrial6$FitRatios[,2:4], 1, FUN=mean)),
        
        c(apply(ATrial7$FitRatios[,2:4], 1, FUN=mean), apply(ATrial8$FitRatios[,2:4], 1, FUN=mean), 
          apply(ATrial9$FitRatios[,2:4], 1, FUN=mean), apply(ATrial10$FitRatios[,2:4], 1, FUN=mean),
          apply(ATrial11$FitRatios[,2:4], 1, FUN=mean), apply(ATrial12$FitRatios[,2:4], 1, FUN=mean), 
          apply(CTrial10$FitRatios[,2:4], 1, FUN=mean), apply(CTrial11$FitRatios[,2:4], 1, FUN=mean),
          apply(CTrial12$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial7$FitRatios[,2:4], 1, FUN=mean), apply(BTrial8$FitRatios[,2:4], 1, FUN=mean), 
          apply(BTrial9$FitRatios[,2:4], 1, FUN=mean), apply(BTrial10$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial11$FitRatios[,2:4], 1, FUN=mean), apply(BTrial12$FitRatios[,2:4], 1, FUN=mean), 
          apply(DTrial10$FitRatios[,2:4], 1, FUN=mean), apply(CTrial11$FitRatios[,2:4], 1, FUN=mean),
          apply(DTrial12$FitRatios[,2:4], 1, FUN=mean)),
        
        c(apply(ATrial13$FitRatios[,2:4], 1, FUN=mean), apply(ATrial14$FitRatios[,2:4], 1, FUN=mean), 
          apply(ATrial15$FitRatios[,2:4], 1, FUN=mean), apply(ATrial16$FitRatios[,2:4], 1, FUN=mean),
          apply(ATrial17$FitRatios[,2:4], 1, FUN=mean), apply(ATrial18$FitRatios[,2:4], 1, FUN=mean), 
          apply(CTrial16$FitRatios[,2:4], 1, FUN=mean), apply(CTrial17$FitRatios[,2:4], 1, FUN=mean),
          apply(CTrial18$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial13$FitRatios[,2:4], 1, FUN=mean), apply(BTrial14$FitRatios[,2:4], 1, FUN=mean), 
          apply(BTrial15$FitRatios[,2:4], 1, FUN=mean), apply(BTrial16$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial17$FitRatios[,2:4], 1, FUN=mean), apply(BTrial18$FitRatios[,2:4], 1, FUN=mean), 
          apply(DTrial16$FitRatios[,2:4], 1, FUN=mean), apply(CTrial17$FitRatios[,2:4], 1, FUN=mean),
          apply(DTrial18$FitRatios[,2:4], 1, FUN=mean)),
        
        c(apply(ATrial19$FitRatios[,2:4], 1, FUN=mean), apply(ATrial20$FitRatios[,2:4], 1, FUN=mean), 
          apply(ATrial21$FitRatios[,2:4], 1, FUN=mean), apply(ATrial22$FitRatios[,2:4], 1, FUN=mean),
          apply(ATrial23$FitRatios[,2:4], 1, FUN=mean), apply(ATrial24$FitRatios[,2:4], 1, FUN=mean), 
          apply(CTrial22$FitRatios[,2:4], 1, FUN=mean), apply(CTrial23$FitRatios[,2:4], 1, FUN=mean),
          apply(CTrial24$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial19$FitRatios[,2:4], 1, FUN=mean), apply(BTrial20$FitRatios[,2:4], 1, FUN=mean), 
          apply(BTrial21$FitRatios[,2:4], 1, FUN=mean), apply(BTrial22$FitRatios[,2:4], 1, FUN=mean),
          apply(BTrial23$FitRatios[,2:4], 1, FUN=mean), apply(BTrial24$FitRatios[,2:4], 1, FUN=mean), 
          apply(DTrial22$FitRatios[,2:4], 1, FUN=mean), apply(CTrial23$FitRatios[,2:4], 1, FUN=mean),
          apply(DTrial24$FitRatios[,2:4], 1, FUN=mean)),
        
        ylim=c(0,1), col=c(rep('deepskyblue',2),'lightskyblue','darkolivegreen3','darkolivegreen1'),
        names=c('constant F','increasing F','decreasing F', 'incr-decrease F'))
title('Fit Ratio', line=-1.4)
dev.off()


##### SAME SCALE ############
setwd("N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Plots")

png(filename="SN_BiasRMSE.png", 
    type="cairo",
    units="mm", 
    width=300, 
    height=375, 
    pointsize=16, 
    res=600)

par(mfrow=c(4,2), mar=c(0.6,1.6,1.1,0), tcl=-0.25, mgp=c(2, 0.4, 0))
vioplot(c(ATrial1$biasAll), c(ATrial2$biasAll), c( ATrial3$biasAll), c( ATrial4$biasAll), c( ATrial5$biasAll), c( ATrial6$biasAll), 
        c(CTrial4$biasAll), c( CTrial5$biasAll), c( CTrial6$biasAll), 
        c(BTrial1$biasAll), c( BTrial2$biasAll), c( BTrial3$biasAll), c( BTrial4$biasAll), c( BTrial5$biasAll), c( BTrial6$biasAll),
        c(DTrial4$biasAll), c( DTrial5$biasAll), c( DTrial6$biasAll),
        ylim=c(-5,5),
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('No Trend', line=-1)
title('Annual Bias', line=0.4, cex=1.5)


vioplot(ATrial1$RMSE, ATrial2$RMSE, ATrial3$RMSE, ATrial4$RMSE, ATrial5$RMSE, ATrial6$RMSE, 
        CTrial4$RMSE, CTrial5$RMSE, CTrial6$RMSE, 
        BTrial1$RMSE, BTrial2$RMSE, BTrial3$RMSE, BTrial4$RMSE, BTrial5$RMSE, BTrial6$RMSE,
        DTrial4$RMSE, DTrial5$RMSE, DTrial6$RMSE,
        ylim=c(0,2), 
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('No Trend', line=-1)
title('RMSE', line=0.4, cex=1.5)



par(mar=c(0.6,1.6,0,0), tcl=-0.25, mgp=c(2, 0.4, 0))
vioplot(c(ATrial7$biasAll), c( ATrial8$biasAll), c( ATrial9$biasAll), c( ATrial10$biasAll), c( ATrial11$biasAll), c( ATrial12$biasAll), 
        c(CTrial10$biasAll), c( CTrial11$biasAll), c( CTrial12$biasAll),  
        c(BTrial7$biasAll), c( BTrial8$biasAll), c( BTrial9$biasAll), c( BTrial10$biasAll), c( BTrial11$biasAll), c( BTrial12$biasAll), 
        c(DTrial10$biasAll), c( DTrial11$biasAll), c( DTrial12$biasAll),
        ylim=c(-5,5),
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Decreasing Trend', line=-1)

vioplot(ATrial7$RMSE, ATrial8$RMSE, ATrial9$RMSE, ATrial10$RMSE, ATrial11$RMSE, ATrial12$RMSE, 
        CTrial10$RMSE, CTrial11$RMSE, CTrial12$RMSE,  
        BTrial7$RMSE, BTrial8$RMSE, BTrial9$RMSE, BTrial10$RMSE, BTrial11$RMSE, BTrial12$RMSE, 
        DTrial10$RMSE, DTrial11$RMSE, DTrial12$RMSE,
        ylim=c(0,2), 
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Decreasing Trend', line=-1)




vioplot(c(ATrial13$biasAll),c(ATrial14$biasAll),c( ATrial15$biasAll),c( ATrial16$biasAll),c( ATrial17$biasAll),c( ATrial18$biasAll), 
        c(CTrial16$biasAll),c( CTrial17$biasAll),c( CTrial18$biasAll), 
        c(BTrial13$biasAll),c( BTrial14$biasAll),c( BTrial15$biasAll),c( BTrial16$biasAll),c( BTrial17$biasAll),c( BTrial18$biasAll), 
        c(DTrial16$biasAll),c( DTrial17$biasAll),c( DTrial18$biasAll),
        ylim=c(-5,5), 
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Increasing Trend', line=-1)


vioplot(ATrial13$RMSE, ATrial14$RMSE, ATrial15$RMSE, ATrial16$RMSE, ATrial17$RMSE, ATrial18$RMSE, 
        CTrial16$RMSE, CTrial17$RMSE, CTrial18$RMSE, 
        BTrial13$RMSE, BTrial14$RMSE, BTrial15$RMSE, BTrial16$RMSE, BTrial17$RMSE, BTrial18$RMSE, 
        DTrial16$RMSE, DTrial17$RMSE, DTrial18$RMSE,
        ylim=c(0,2), 
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Increasing Trend', line=-1)

vioplot(c(ATrial19$biasAll),c( ATrial20$biasAll),c( ATrial21$biasAll),c( ATrial22$biasAll),c( ATrial23$biasAll),c( ATrial24$biasAll), 
        c(CTrial22$biasAll),c( CTrial23$biasAll),c( CTrial24$biasAll),  
        c(BTrial19$biasAll),c( BTrial20$biasAll),c( BTrial21$biasAll),c( BTrial22$biasAll),c( BTrial23$biasAll),c(BTrial24$biasAll), 
        c(DTrial22$biasAll),c( DTrial23$biasAll),c( DTrial24$biasAll),
        ylim=c(-5,5),
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Decrease - Increasing Trend', line=-1)


vioplot(ATrial19$RMSE, ATrial20$RMSE, ATrial21$RMSE, ATrial22$RMSE, ATrial23$RMSE, ATrial24$RMSE, 
        CTrial22$RMSE, CTrial23$RMSE, CTrial24$RMSE,  
        BTrial19$RMSE, BTrial20$RMSE, BTrial21$RMSE, BTrial22$RMSE, BTrial23$RMSE, BTrial24$RMSE, 
        DTrial22$RMSE, DTrial23$RMSE, DTrial24$RMSE,
        ylim=c(0,2), 
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Decrease - Increasing Trend', line=-1)

dev.off()









##### CHANGE SCALE ############
setwd("N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Plots")

png(filename="SN_BiasRMSE_Scale.png", 
    type="cairo",
    units="mm", 
    width=300, 
    height=375, 
    pointsize=16, 
    res=600)

par(mfrow=c(4,2), mar=c(0.6,1.6,1.1,0), tcl=-0.25, mgp=c(2, 0.4, 0))
vioplot(c(ATrial1$biasAll), c(ATrial2$biasAll), c( ATrial3$biasAll), c( ATrial4$biasAll), c( ATrial5$biasAll), c( ATrial6$biasAll), 
        c(CTrial4$biasAll), c( CTrial5$biasAll), c( CTrial6$biasAll), 
        c(BTrial1$biasAll), c( BTrial2$biasAll), c( BTrial3$biasAll), c( BTrial4$biasAll), c( BTrial5$biasAll), c( BTrial6$biasAll),
        c(DTrial4$biasAll), c( DTrial5$biasAll), c( DTrial6$biasAll),
        ylim=c(-5,5),
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('No Trend', line=-1)
title('Annual Bias', line=0.4, cex=1.5)


vioplot(ATrial1$RMSE, ATrial2$RMSE, ATrial3$RMSE, ATrial4$RMSE, ATrial5$RMSE, ATrial6$RMSE, 
        CTrial4$RMSE, CTrial5$RMSE, CTrial6$RMSE, 
        BTrial1$RMSE, BTrial2$RMSE, BTrial3$RMSE, BTrial4$RMSE, BTrial5$RMSE, BTrial6$RMSE,
        DTrial4$RMSE, DTrial5$RMSE, DTrial6$RMSE,
        ylim=c(0,2), 
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('No Trend', line=-1)
title('RMSE', line=0.4, cex=1.5)



par(mar=c(0.6,1.6,0,0), tcl=-0.25, mgp=c(2, 0.4, 0))
vioplot(c(ATrial7$biasAll), c( ATrial8$biasAll), c( ATrial9$biasAll), c( ATrial10$biasAll), c( ATrial11$biasAll), c( ATrial12$biasAll), 
        c(CTrial10$biasAll), c( CTrial11$biasAll), c( CTrial12$biasAll),  
        c(BTrial7$biasAll), c( BTrial8$biasAll), c( BTrial9$biasAll), c( BTrial10$biasAll), c( BTrial11$biasAll), c( BTrial12$biasAll), 
        c(DTrial10$biasAll), c( DTrial11$biasAll), c( DTrial12$biasAll),
        ylim=c(-1,1),
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Decreasing Trend', line=-1)

vioplot(ATrial7$RMSE, ATrial8$RMSE, ATrial9$RMSE, ATrial10$RMSE, ATrial11$RMSE, ATrial12$RMSE, 
        CTrial10$RMSE, CTrial11$RMSE, CTrial12$RMSE,  
        BTrial7$RMSE, BTrial8$RMSE, BTrial9$RMSE, BTrial10$RMSE, BTrial11$RMSE, BTrial12$RMSE, 
        DTrial10$RMSE, DTrial11$RMSE, DTrial12$RMSE,
        ylim=c(0,0.3), 
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Decreasing Trend', line=-1)




vioplot(c(ATrial13$biasAll),c(ATrial14$biasAll),c( ATrial15$biasAll),c( ATrial16$biasAll),c( ATrial17$biasAll),c( ATrial18$biasAll), 
        c(CTrial16$biasAll),c( CTrial17$biasAll),c( CTrial18$biasAll), 
        c(BTrial13$biasAll),c( BTrial14$biasAll),c( BTrial15$biasAll),c( BTrial16$biasAll),c( BTrial17$biasAll),c( BTrial18$biasAll), 
        c(DTrial16$biasAll),c( DTrial17$biasAll),c( DTrial18$biasAll),
        ylim=c(-1,1), 
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Increasing Trend', line=-1)


vioplot(ATrial13$RMSE, ATrial14$RMSE, ATrial15$RMSE, ATrial16$RMSE, ATrial17$RMSE, ATrial18$RMSE, 
        CTrial16$RMSE, CTrial17$RMSE, CTrial18$RMSE, 
        BTrial13$RMSE, BTrial14$RMSE, BTrial15$RMSE, BTrial16$RMSE, BTrial17$RMSE, BTrial18$RMSE, 
        DTrial16$RMSE, DTrial17$RMSE, DTrial18$RMSE,
        ylim=c(0,0.3), 
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Increasing Trend', line=-1)

vioplot(c(ATrial19$biasAll),c( ATrial20$biasAll),c( ATrial21$biasAll),c( ATrial22$biasAll),c( ATrial23$biasAll),c( ATrial24$biasAll), 
        c(CTrial22$biasAll),c( CTrial23$biasAll),c( CTrial24$biasAll),  
        c(BTrial19$biasAll),c( BTrial20$biasAll),c( BTrial21$biasAll),c( BTrial22$biasAll),c( BTrial23$biasAll),c(BTrial24$biasAll), 
        c(DTrial22$biasAll),c( DTrial23$biasAll),c( DTrial24$biasAll),
        ylim=c(-5,5),
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Decrease - Increasing Trend', line=-1)


vioplot(ATrial19$RMSE, ATrial20$RMSE, ATrial21$RMSE, ATrial22$RMSE, ATrial23$RMSE, ATrial24$RMSE, 
        CTrial22$RMSE, CTrial23$RMSE, CTrial24$RMSE,  
        BTrial19$RMSE, BTrial20$RMSE, BTrial21$RMSE, BTrial22$RMSE, BTrial23$RMSE, BTrial24$RMSE, 
        DTrial22$RMSE, DTrial23$RMSE, DTrial24$RMSE,
        ylim=c(0,1.5), 
        col=c(rep('deepskyblue3',4),rep('deepskyblue',3),rep('lightskyblue1',3),rep('darkolivegreen4',3),rep('darkolivegreen3',3),rep('darkolivegreen1',3)),
        names=FALSE)
title('Decrease - Increasing Trend', line=-1)

dev.off()



##### Legend ############
setwd("N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\Plots")

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

setwd("N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\FourIndicesGradual\\Trial22")

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




setwd("N:\\Documents\\DFA_Simulation\\DFA_Simulation\\Atl_SN\\\\Plots")

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
legend('topright', c("DFA predicted trend +/- 90% CIs",'Rescaled simulated abundance'), lwd=2, col=c('black','darkorchid'))

dev.off()