## Script to run the TB model as compiled C code - adapted to do first pass of COR modelling

## Calibrate model to: TB incidence (by HIV)
##                     TB prevalence
##                     TB mortality (by HIV)
## by varying reactivation rate, proportion primary disease, beta, HIV RRs 
## (should give us range of disease due to transmission which will be important in determining indirect impact of strategy)
## use an importance resampling approach 

## then run the model assuming 15% COR positivity which accounts of 70% of disease and treatment efficacy of 80% (put uncertainty on these?)
## prevent 56% of HIV- adult TB cases in each year (reduce incidence flows accordingly) and but these individuals back into the latent box
## repeat for n years
## also count numbers screened and treated?
## repeat for QFT+?

## Set working directory
setwd("C:/Users/TOM SUMMER/Documents/COR_modelling/COR")

## Load packages
library(deSolve)
library(reshape2)
library(ggplot2)

## Compile and load the C code
system("R CMD SHLIB TB_model_v4.c") # Compile
dyn.load("TB_model_v4.dll") # Load
dyn.unload("TB_model_v4.dll") # Unload - need to do this before recompiling

##############################################################################################################################

## Load UN population data
UN_pop_age <- as.data.frame(read.table("SA_pop_age.txt",header=TRUE)) # Load UN Population data
UN_pop_age_low <- as.data.frame(read.table("SA_pop_age_low.txt",header=TRUE)) # Load UN Population data
UN_pop_age_high <- as.data.frame(read.table("SA_pop_age_high.txt",header=TRUE)) # Load UN Population data
# add total to data
UN_pop_age_t <- cbind(UN_pop_age,rowSums(UN_pop_age[,2:18]))
colnames(UN_pop_age_t) <- c(colnames(UN_pop_age),"Total")
UN_pop_age_low_t <- cbind(UN_pop_age_low,rowSums(UN_pop_age_low[,2:18]))
colnames(UN_pop_age_low_t) <- c(colnames(UN_pop_age_low),"Total")
UN_pop_age_high_t <- cbind(UN_pop_age_high,rowSums(UN_pop_age_high[,2:18]))
colnames(UN_pop_age_high_t) <- c(colnames(UN_pop_age_high),"Total")

# Load TB burden data and calculate variances etc for use in likelihood
TB_data <- as.data.frame(read.table("TB_burden.txt",header=TRUE))
# Variances in data approximated from WHO ranges assuming they represent 95% range of normal distribution
v_prev<-(((TB_data[,"prev_hi"]-TB_data[,"prev_lo"])/2)/1.96)^2
v_inc<-(((TB_data[,"inc_hi"]-TB_data[,"inc_lo"])/2)/1.96)^2
v_mort<-(((TB_data[,"mort_hi"]-TB_data[,"mort_lo"])/2)/1.96)^2

# Set up age structure
ages <- c(4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,79,100) # upper end of age classes
num_ages <- length(ages) # calculates the number of age classes

# Set up the forcing functions for birth and death - all from 1970 onwards

## Fertility
# Uses crude birth rate (per 1000 population) from UN pop for South Africa which has values for 5 year periods
# values post 2010 based on medium fertiliy
birth_rate <- cbind(seq(1972.5,2047.5,5),
                    c(38,36,34,31,27,25,24,22,21,20,18,17,17,16,15,14))

## Survival
Survive_age <- as.data.frame(read.table("SA_survival_age.txt",header=TRUE)) # Load survival proportions calculated from life tables
# Proportion surviving from age 0 to 1 - used to determine entry to first age group
s_birth <- cbind(Survive_age$Year,Survive_age$X1)
# and to other ages
s5 <- cbind(Survive_age$Year,Survive_age$X5)
s10 <- cbind(Survive_age$Year,Survive_age$X10)
s15 <- cbind(Survive_age$Year,Survive_age$X15)
s20 <- cbind(Survive_age$Year,Survive_age$X20)
s25 <- cbind(Survive_age$Year,Survive_age$X25)
s30 <- cbind(Survive_age$Year,Survive_age$X30)
s35 <- cbind(Survive_age$Year,Survive_age$X35)
s40 <- cbind(Survive_age$Year,Survive_age$X40)
s45 <- cbind(Survive_age$Year,Survive_age$X45)
s50 <- cbind(Survive_age$Year,Survive_age$X50)
s55 <- cbind(Survive_age$Year,Survive_age$X55)
s60 <- cbind(Survive_age$Year,Survive_age$X60)
s65 <- cbind(Survive_age$Year,Survive_age$X65)
s70 <- cbind(Survive_age$Year,Survive_age$X70)
s75 <- cbind(Survive_age$Year,Survive_age$X75)
s80 <- cbind(Survive_age$Year,Survive_age$X80)

# HIV Incidence by age and year - based on AIM output, but ignoring childhood infections (this is what Carel does in TIME)
HIV_Inc_age <- as.data.frame(read.table("HIV_Inc_age.txt",header=TRUE)) # Load HIV incidence data taken from AIM                                       # Data from AIM is rate per 1000 
#HIV_Inc_age[,2:18]=HIV_Inc_age[,2:18]*0

h0 <- 0*cbind(HIV_Inc_age$Year,HIV_Inc_age$X0/1000)
h5 <- 0*cbind(HIV_Inc_age$Year,HIV_Inc_age$X5/1000)
h10 <- 0*cbind(HIV_Inc_age$Year,HIV_Inc_age$X10/1000)
h15 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X15/1000)
h20 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X20/1000)
h25 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X25/1000)
h30 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X30/1000)
h35 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X35/1000)
h40 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X40/1000)
h45 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X45/1000)
h50 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X50/1000)
h55 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X55/1000)
h60 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X60/1000)
h65 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X65/1000)
h70 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X70/1000)
h75 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X75/1000)
h80 <- cbind(HIV_Inc_age$Year,HIV_Inc_age$X80/1000)

# ART coverage - based on AIM, use CD4 eligibility threshold and % of those in need on ART
ART_data <- as.data.frame(read.table("ART_data.txt",header=TRUE)) # Load data
# Create forcing function of threshold category
Athresh <- cbind(ART_data[,"Year"],ART_data[,"CD4_cat"])
# Create forcing functions which account for threshold and coverage
A50 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>50,1,0)/100)
A99 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>99,1,0)/100)
A199 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>199,1,0)/100)
A249 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>249,1,0)/100)
A349 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>349,1,0)/100)
A500 <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>499,1,0)/100)
Ahigh <- cbind(ART_data[,"Year"],ART_data[,"Percent"]*ifelse(ART_data[,"CD4_t"]>500,1,0)/100)

# Pop adjust - to turn off population adjust for TB/HIV deaths from 2015 onwards
pop_ad <- cbind(c(2014,2015,2016),c(1,0,0))

# BCG coverage - currently assume 90% at all times
BCG_cov <- cbind(c(1972,1973,2050),c(0.9,0.9,0.9))

# Case detection rate - generalised logistic function, don't think this quite matches TIME 
k <- cbind(seq(1970,2050),(30 + ((120-30)/((1+exp(-0.5*(seq(1970,2050)-2004)))^(1/2))))/100)

# DST coverage among new and previously treated cases
dst_n <- cbind(seq(1970,2050),(0 + ((95-0)/((1+exp(-1*(seq(1970,2050)-1993)))^(1/2))))/100)
dst_p <- cbind(seq(1970,2050),(0 + ((95-0)/((1+exp(-1*(seq(1970,2050)-1993)))^(1/2))))/100)

##### COR screen and treat intervention - currently turned off
# Screen 20% of HIV- adult population per year
# Avert 0.7*treatment_eff of the disease in that 20% - move these people to latent. zero for children
COR_prevent <- cbind(c(0,0,0),c(0,0,0))
# Also diagnose 20% of prevalent cases through excluding active disease (symptom screen with assumed sensitivity of 70% - HIV- population, WHO review)
# http://www.who.int/tb/Review2Accuracyofscreeningtests.pdf
COR_diagnose <- cbind(c(0,0,0),c(0,0,0))
 
# Combine forcing functions into a list
force <- list(birth_rate,s_birth,s5,s10,s15,s20,s25,s30,s35,s40,s45,s50,s55,s60,s65,s70,s75,s80,
              h0,h5,h10,h15,h20,h25,h30,h35,h40,h45,h50,h55,h60,h65,h70,h75,h80,
              Ahigh,A500,A349,A249,A199,A99,A50,Athresh,
              BCG_cov,pop_ad,k,dst_n,dst_p,
              COR_prevent,COR_diagnose)

# Fit model using importance resampling approach

# Number of runs to do
n_run = 100
# Years to run model for - always start at 1970
times_run <- seq(1970,2035 , by=1)

# Create arrays for storing outputs and likelihood 
I = mat.or.vec(length(times_run),n_run)
M = I
P = I    
LTBI = I
per_HIV = I

# Calculate likelihood at each time point (1990-2013) for each output
L_inc = mat.or.vec(n_run,24)
L_mort = L_inc
L_prev = L_inc
L = rep(0,n_run)

# Sampled parameters
beta <- rep(0,n_run)
a_a <- beta
v <- beta

ptm <- proc.time()

for (runs in 1:n_run){
  
  # sample parameters - just beta for now
  
  beta[runs] = runif(1,17,20)
  a_a[runs] = runif(1,0.08,0.15)
  v[runs] = runif(1,0.0001,0.0025)
  print(c(runs,beta[runs],a_a[runs],v[runs]))
  
  # Fitness of MDR, used to calculate parameter for superinfections
  # Both of these are passed into "parms" together with e, the MDR acquisition rate (we set this to zero to exclude MDR in equilibirum phase)

  fit_cost=0.7
  g = fit_cost/(1+fit_cost) # superinfections 
  e = 0.01
  
  # proportion primary (a), proportion smear pos (sig) and mortality rates (muN and mu_I) take different values for 
  # adults (>15) (_a), 0-4 (_0), 5-9 (_5) and 10-14 (_10)
  
  # create parameter vector to pass to model
  parms <- c(age1 = 1/5, age2 = 1/21, beta = beta[runs], 
            a_a = a_a[runs], a0 = 0.26432, a5 = 0.14056, a10 = 0.056,  
            p = 0.65, v = v[runs], 
            sig_a = 0.45, sig0 = 0.0684, sig5 = 0.0414, sig10 = 0.0846, rel_inf = 0.25, theta = 0.02, r = 0.25, 
            mu_N = 0.25, mu_N0 = 0.426, mu_I = 0.35, mu_I0 = 0.59, fit_cost = fit_cost, e = e, g=g, l_s = 0.83, l_m = 0.7, d = 0.8, tau_s = 0.76, tau_m = 0.5,
            eff_n = 0.61, eff_p = 0.45, 
            muN_H = 0.45, muI_H = 0.6, RR1a = 2, RR2a = 1.288, RR1v = 3, RR2v = 3, RR1p = 0.5, RR2p = 1.1,
            ART_TB1 = 0.7, ART_TB2 = 0.5, ART_TB3 = 0.35, ART_mort1 = 0.5, ART_mort2 = 0.4, ART_mort3 = 0.3,
            BCG_eff = 0.39,
            sig_H = 0.35,r_H=0.15)

  # Run the model for eq.

  # First run the model from 1970 pop with 1970 birth/death rates and care and control parameters with 100 TB cases (no MDR or HIV) 
  # for 100 years to get stable age structure and disease state

  # Times to run model for
  times <- seq(0,100, by=1)

  # Initial conditions - all susceptible
  temp <- c()
  for (i in 1:num_ages){temp[i]<-UN_pop_age[21,i+1]}
  xstart <- c(S=c(temp),
              Lsn=rep(0,num_ages),Lsp=rep(0,num_ages),Lmn=rep(0,num_ages),Lmp=rep(0,num_ages),
              Nsn=rep(0,num_ages),Nsp=rep(0,num_ages),Nmn=rep(0,num_ages),Nmp=rep(0,num_ages),
              Isn=c(rep(0,5),100,rep(0,11)),Isp=rep(0,num_ages),Imn=rep(0,num_ages),Imp=rep(0,num_ages),
              S_H=rep(0,num_ages*7),
              Lsn_H=rep(0,num_ages*7),Lsp_H=rep(0,num_ages*7),Lmn_H=rep(0,num_ages*7),Lmp_H=rep(0,num_ages*7),
              Nsn_H=rep(0,num_ages*7),Nsp_H=rep(0,num_ages*7),Nmn_H=rep(0,num_ages*7),Nmp_H=rep(0,num_ages*7),
              Isn_H=rep(0,num_ages*7),Isp_H=rep(0,num_ages*7),Imn_H=rep(0,num_ages*7),Imp_H=rep(0,num_ages*7),
              S_A=rep(0,num_ages*7*3),
              Lsn_A=rep(0,num_ages*7*3),Lsp_A=rep(0,num_ages*7*3),Lmn_A=rep(0,num_ages*7*3),Lmp_A=rep(0,num_ages*7*3),
              Nsn_A=rep(0,num_ages*7*3),Nsp_A=rep(0,num_ages*7*3),Nmn_A=rep(0,num_ages*7*3),Nmp_A=rep(0,num_ages*7*3),
              Isn_A=rep(0,num_ages*7*3),Isp_A=rep(0,num_ages*7*3),Imn_A=rep(0,num_ages*7*3),Imp_A=rep(0,num_ages*7*3))

  # For initialisation run turn off MDR by setting e = 0
  parms["e"]=0

  # Run the model
  time_eq <- system.time(out_eq <- ode(y=xstart, times, func = "derivsc",
                         parms = parms, dllname = "TB_model_v4",initforc = "forcc",
                         forcings=force, initfunc = "parmsc", nout = 44,
                         outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                         "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                         "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                         "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50",
                         "v1","v2","v3","v4","v5","v6","v7","ART_tot","ART_need","ART_new","ART_on","TB_deaths",
                         "Cases_neg","Cases_pos","Cases_ART"), method = rkMethod("rk34f")))

  # Adjust pop down to 1970 values and reassign initial conditions - model can now be run from 1970 with TB and HIV
  temp <- out_eq[dim(out_eq)[1],2:6410]
  temp <- temp/(sum(temp)/22502) # 22502 is total pop from UN estimates in 1970)
  xstart <- temp

  # Reset e to allow MDR
  parms["e"]=e

  # Run the model
  time_run <-system.time(out <- ode(y=xstart, times_run, func = "derivsc",
                         parms = parms, dllname = "TB_model_v4",initforc = "forcc",
                         forcings=force, initfunc = "parmsc", nout = 44,
                         outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                        "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                        "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                        "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50",
                        "v1","v2","v3","v4","v5","v6","v7","ART_tot","ART_need","ART_new","ART_on","TB_deaths",
                        "Cases_neg","Cases_pos","Cases_ART"), method = rkMethod("rk34f")))
 
  # Store things we need for likelihood
  
  I[,runs] = 100000*(out[,"Cases_neg"]+out[,"Cases_pos"]+out[,"Cases_ART"])/out[,"Total"]  # TB incidence
  M[,runs] = 100000*out[,"TB_deaths"]/out[,"Total"]                                        # TB mortality
  P[,runs] = 100000*(out[,"Total_DS"]+out[,"Total_MDR"])/out[,"Total"]                     # TB prevalence  
  LTBI[,runs] = 100*(out[,"Total_L"]+out[,"Total_DS"]+out[,"Total_MDR"])/out[,"Total"]     # LTBI
  per_HIV[,runs] = 100*(out[,"Cases_pos"]+out[,"Cases_ART"])/(out[,"Cases_neg"]+out[,"Cases_pos"]+out[,"Cases_ART"]) # % incident TB in HIV+
    
  # Calculate likelihood at each time point for each output
  L_inc[runs,]<-((2*pi*v_inc)^(-1/2))*exp(-(1/(2*v_inc))*((I[,runs][21:44]-TB_data[,"inc_mid"])^2))
  L_mort[runs,]<-((2*pi*v_mort)^(-1/2))*exp(-(1/(2*v_mort))*((M[,runs][21:44]-TB_data[,"mort_mid"])^2))
  L_prev[runs,]<-((2*pi*v_prev)^(-1/2))*exp(-(1/(2*v_prev))*((P[,runs][21:44]-TB_data[,"prev_mid"])^2))
  
  # Total likelihood is based on 2004 and 2013 data points (last 10 years) (product of likelihood for each data point)
  L[runs]<-prod(L_inc[runs,15:24])*prod(L_mort[runs,15:24])*prod(L_prev[runs,15:24])
  
}

proc.time() - ptm
  
# Resample parameter sets based on likelihood to get a set of good parameters
N_resamp<-10000
t<-sample(seq(1:n_run),N_resamp,replace=TRUE,prob=L/sum(L))
unique_t<-unique(t)

## Now rerun for each unique accepted sample with COR screening turned on

##### COR screen and treat intervention - turn it on
# Screen 20% of HIV- adult population per year
# Avert 0.7*treatment_eff of the disease in that 20% - move these people to latent. zero for children

Ic = mat.or.vec(length(times_run),length(unique_t))
Mc = Ic
Pc = Ic   
LTBIc = Ic

Ic20 = Ic
Mc20 = Ic
Pc20 = Ic
LTBIc20 = Ic

Ic40 = Ic
Mc40 = Ic
Pc40 = Ic
LTBIc40 = Ic

Is <- Ic
Ms <- Ic
Ps <- Ic
LTBIs <- Ic

# sort the parameter sets to run into ascending order
t_to_run <- sort(unique_t)
for (runs in 1:length(unique_t)){
  
  # % screen = 30%
  p_screen <- 0.3
  # sensitivity of COR = 70%
  sens_COR <- 0.7
  # treatment efficacy = 70% (60-80)
  v_eff<-(((0.8-0.6)/2)/1.96)
  treat_eff <- rnorm(1,0.7,v_eff)
  # sensitivity of COR for active disease = 85% (80-90)
  v_active<-(((0.9-0.8)/2)/1.96)
  sens_active <- rnorm(1,0.85,v_eff)
  # If just ACF would use symptom screen - 60% sensitivity
  # http://www.who.int/tb/Review2Accuracyofscreeningtests.pdf
  sens_sym <- 0.6
  
  COR_prevent <- cbind(c(1970,2020,2021),c(0,0,p_screen*sens_COR*treat_eff))
  COR_diagnose <- cbind(c(1970,2020,2021),c(0,0,p_screen*sens_active))

  # Combine forcing functions into a list
  force <- list(birth_rate,s_birth,s5,s10,s15,s20,s25,s30,s35,s40,s45,s50,s55,s60,s65,s70,s75,s80,
                h0,h5,h10,h15,h20,h25,h30,h35,h40,h45,h50,h55,h60,h65,h70,h75,h80,
                Ahigh,A500,A349,A249,A199,A99,A50,Athresh,
                BCG_cov,pop_ad,k,dst_n,dst_p,
                COR_prevent,COR_diagnose)
  
  j <- t_to_run[runs]
  
  # create parameter vector to pass to model
  parms <- c(age1 = 1/5, age2 = 1/21, beta = beta[j], 
             a_a = a_a[j], a0 = 0.26432, a5 = 0.14056, a10 = 0.056,  
             p = 0.65, v = v[j], 
             sig_a = 0.45, sig0 = 0.0684, sig5 = 0.0414, sig10 = 0.0846, rel_inf = 0.25, theta = 0.02, r = 0.25, 
             mu_N = 0.25, mu_N0 = 0.426, mu_I = 0.35, mu_I0 = 0.59, fit_cost = fit_cost, e = e, g=g, l_s = 0.83, l_m = 0.7, d = 0.8, tau_s = 0.76, tau_m = 0.5,
             eff_n = 0.61, eff_p = 0.45, 
             muN_H = 0.45, muI_H = 0.6, RR1a = 2, RR2a = 1.288, RR1v = 3, RR2v = 3, RR1p = 0.5, RR2p = 1.1,
             ART_TB1 = 0.7, ART_TB2 = 0.5, ART_TB3 = 0.35, ART_mort1 = 0.5, ART_mort2 = 0.4, ART_mort3 = 0.3,
             BCG_eff = 0.39,
             sig_H = 0.35,r_H=0.15)
  
  # Run the model for eq.
  
  # First run the model from 1970 pop with 1970 birth/death rates and care and control parameters with 100 TB cases (no MDR or HIV) 
  # for 100 years to get stable age structure and disease state
  
  # Times to run model for
  times <- seq(0,100, by=1)
  
  # Initial conditions - all susceptible
  temp <- c()
  for (i in 1:num_ages){temp[i]<-UN_pop_age[21,i+1]}
  xstart <- c(S=c(temp),
              Lsn=rep(0,num_ages),Lsp=rep(0,num_ages),Lmn=rep(0,num_ages),Lmp=rep(0,num_ages),
              Nsn=rep(0,num_ages),Nsp=rep(0,num_ages),Nmn=rep(0,num_ages),Nmp=rep(0,num_ages),
              Isn=c(rep(0,5),100,rep(0,11)),Isp=rep(0,num_ages),Imn=rep(0,num_ages),Imp=rep(0,num_ages),
              S_H=rep(0,num_ages*7),
              Lsn_H=rep(0,num_ages*7),Lsp_H=rep(0,num_ages*7),Lmn_H=rep(0,num_ages*7),Lmp_H=rep(0,num_ages*7),
              Nsn_H=rep(0,num_ages*7),Nsp_H=rep(0,num_ages*7),Nmn_H=rep(0,num_ages*7),Nmp_H=rep(0,num_ages*7),
              Isn_H=rep(0,num_ages*7),Isp_H=rep(0,num_ages*7),Imn_H=rep(0,num_ages*7),Imp_H=rep(0,num_ages*7),
              S_A=rep(0,num_ages*7*3),
              Lsn_A=rep(0,num_ages*7*3),Lsp_A=rep(0,num_ages*7*3),Lmn_A=rep(0,num_ages*7*3),Lmp_A=rep(0,num_ages*7*3),
              Nsn_A=rep(0,num_ages*7*3),Nsp_A=rep(0,num_ages*7*3),Nmn_A=rep(0,num_ages*7*3),Nmp_A=rep(0,num_ages*7*3),
              Isn_A=rep(0,num_ages*7*3),Isp_A=rep(0,num_ages*7*3),Imn_A=rep(0,num_ages*7*3),Imp_A=rep(0,num_ages*7*3))
  
  # For initialisation run turn off MDR by setting e = 0
  parms["e"]=0
  
  # Run the model
  time_eq <- system.time(out_eq <- ode(y=xstart, times, func = "derivsc",
                                       parms = parms, dllname = "TB_model_v4",initforc = "forcc",
                                       forcings=force, initfunc = "parmsc", nout = 44,
                                       outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                                                    "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                                                    "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                                                    "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50",
                                                    "v1","v2","v3","v4","v5","v6","v7","ART_tot","ART_need","ART_new","ART_on","TB_deaths",
                                                    "Cases_neg","Cases_pos","Cases_ART"), method = rkMethod("rk34f")))
  
  # Adjust pop down to 1970 values and reassign initial conditions - model can now be run from 1970 with TB and HIV
  temp <- out_eq[dim(out_eq)[1],2:6410]
  temp <- temp/(sum(temp)/22502) # 22502 is total pop from UN estimates in 1970)
  xstart <- temp
  
  # Reset e to allow MDR
  parms["e"]=e
  
  # Run the model
  time_run <-system.time(out <- ode(y=xstart, times_run, func = "derivsc",
                                    parms = parms, dllname = "TB_model_v4",initforc = "forcc",
                                    forcings=force, initfunc = "parmsc", nout = 44,
                                    outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                                                 "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                                                 "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                                                 "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50",
                                                 "v1","v2","v3","v4","v5","v6","v7","ART_tot","ART_need","ART_new","ART_on","TB_deaths",
                                                 "Cases_neg","Cases_pos","Cases_ART"), method = rkMethod("rk34f")))
  
  Ic[,runs] = 100000*(out[,"Cases_neg"]+out[,"Cases_pos"]+out[,"Cases_ART"])/out[,"Total"]  # TB incidence
  Mc[,runs] = 100000*out[,"TB_deaths"]/out[,"Total"]                                        # TB mortality
  Pc[,runs] = 100000*(out[,"Total_DS"]+out[,"Total_MDR"])/out[,"Total"]                     # TB prevalence  
  LTBIc[,runs] = 100*(out[,"Total_L"]+out[,"Total_DS"]+out[,"Total_MDR"])/out[,"Total"]     # LTBI
  
  # Repeat with just symptom screen i.e. set COR prevent = 0
  
  COR_prevent <- cbind(c(1970,2020,2021),c(0,0,0))
  # Also diagnose 20% of prevalent cases through excluding active disease (symptom screen with assumed sensitivity of 70% - HIV- population, WHO review)
  # http://www.who.int/tb/Review2Accuracyofscreeningtests.pdf
  COR_diagnose <- cbind(c(1970,2020,2021),c(0,0,p_screen*sens_sym))
  
  # Combine forcing functions into a list
  force <- list(birth_rate,s_birth,s5,s10,s15,s20,s25,s30,s35,s40,s45,s50,s55,s60,s65,s70,s75,s80,
                h0,h5,h10,h15,h20,h25,h30,h35,h40,h45,h50,h55,h60,h65,h70,h75,h80,
                Ahigh,A500,A349,A249,A199,A99,A50,Athresh,
                BCG_cov,pop_ad,k,dst_n,dst_p,
                COR_prevent,COR_diagnose)
  
  # Run the model
  time_run <-system.time(out <- ode(y=xstart, times_run, func = "derivsc",
                                    parms = parms, dllname = "TB_model_v4",initforc = "forcc",
                                    forcings=force, initfunc = "parmsc", nout = 44,
                                    outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                                                 "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                                                 "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                                                 "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50",
                                                 "v1","v2","v3","v4","v5","v6","v7","ART_tot","ART_need","ART_new","ART_on","TB_deaths",
                                                 "Cases_neg","Cases_pos","Cases_ART"), method = rkMethod("rk34f")))
  
  Is[,runs] = 100000*(out[,"Cases_neg"]+out[,"Cases_pos"]+out[,"Cases_ART"])/out[,"Total"]  # TB incidence
  Ms[,runs] = 100000*out[,"TB_deaths"]/out[,"Total"]                                        # TB mortality
  Ps[,runs] = 100000*(out[,"Total_DS"]+out[,"Total_MDR"])/out[,"Total"]                     # TB prevalence  
  LTBIs[,runs] = 100*(out[,"Total_L"]+out[,"Total_DS"]+out[,"Total_MDR"])/out[,"Total"]     # LTBI
  
  # And repeat for 20% coverage
  p_screen <- 0.2
  
  COR_prevent <- cbind(c(1970,2020,2021),c(0,0,p_screen*sens_COR*treat_eff))
  # Also diagnose 20% of prevalent cases through excluding active disease (symptom screen with assumed sensitivity of 70% - HIV- population, WHO review)
  # http://www.who.int/tb/Review2Accuracyofscreeningtests.pdf
  COR_diagnose <- cbind(c(1970,2020,2021),c(0,0,p_screen*sens_active))
  
  # Combine forcing functions into a list
  force <- list(birth_rate,s_birth,s5,s10,s15,s20,s25,s30,s35,s40,s45,s50,s55,s60,s65,s70,s75,s80,
                h0,h5,h10,h15,h20,h25,h30,h35,h40,h45,h50,h55,h60,h65,h70,h75,h80,
                Ahigh,A500,A349,A249,A199,A99,A50,Athresh,
                BCG_cov,pop_ad,k,dst_n,dst_p,
                COR_prevent,COR_diagnose)
  
  # Run the model
  time_run <-system.time(out <- ode(y=xstart, times_run, func = "derivsc",
                                    parms = parms, dllname = "TB_model_v4",initforc = "forcc",
                                    forcings=force, initfunc = "parmsc", nout = 44,
                                    outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                                                 "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                                                 "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                                                 "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50",
                                                 "v1","v2","v3","v4","v5","v6","v7","ART_tot","ART_need","ART_new","ART_on","TB_deaths",
                                                 "Cases_neg","Cases_pos","Cases_ART"), method = rkMethod("rk34f")))
  
  Ic20[,runs] = 100000*(out[,"Cases_neg"]+out[,"Cases_pos"]+out[,"Cases_ART"])/out[,"Total"]  # TB incidence
  Mc20[,runs] = 100000*out[,"TB_deaths"]/out[,"Total"]                                        # TB mortality
  Pc20[,runs] = 100000*(out[,"Total_DS"]+out[,"Total_MDR"])/out[,"Total"]                     # TB prevalence  
  LTBIc20[,runs] = 100*(out[,"Total_L"]+out[,"Total_DS"]+out[,"Total_MDR"])/out[,"Total"]     # LTBI
  
  # And repeat for 40% coverage
  p_screen <- 0.4
  
  COR_prevent <- cbind(c(1970,2020,2021),c(0,0,p_screen*sens_COR*treat_eff))
  # Also diagnose 20% of prevalent cases through excluding active disease (symptom screen with assumed sensitivity of 70% - HIV- population, WHO review)
  # http://www.who.int/tb/Review2Accuracyofscreeningtests.pdf
  COR_diagnose <- cbind(c(1970,2020,2021),c(0,0,p_screen*sens_active))
  
  # Combine forcing functions into a list
  force <- list(birth_rate,s_birth,s5,s10,s15,s20,s25,s30,s35,s40,s45,s50,s55,s60,s65,s70,s75,s80,
                h0,h5,h10,h15,h20,h25,h30,h35,h40,h45,h50,h55,h60,h65,h70,h75,h80,
                Ahigh,A500,A349,A249,A199,A99,A50,Athresh,
                BCG_cov,pop_ad,k,dst_n,dst_p,
                COR_prevent,COR_diagnose)
  
  # Run the model
  time_run <-system.time(out <- ode(y=xstart, times_run, func = "derivsc",
                                    parms = parms, dllname = "TB_model_v4",initforc = "forcc",
                                    forcings=force, initfunc = "parmsc", nout = 44,
                                    outnames = c("Total","Total_S","Total_Ls","Total_Lm","Total_L","Total_Ns","Total_Nm",
                                                 "Total_N","Total_Is","Total_Im","Total_I","Total_DS","Total_MDR","FS","FM",
                                                 "CD4500","CD4350_500","CD4250_349","CD4200_249","CD4100_199","CD450_99","CD450",
                                                 "ART500","ART350_500","ART250_349","ART200_249","ART100_199","ART50_99","ART50",
                                                 "v1","v2","v3","v4","v5","v6","v7","ART_tot","ART_need","ART_new","ART_on","TB_deaths",
                                                 "Cases_neg","Cases_pos","Cases_ART"), method = rkMethod("rk34f")))
  
  Ic40[,runs] = 100000*(out[,"Cases_neg"]+out[,"Cases_pos"]+out[,"Cases_ART"])/out[,"Total"]  # TB incidence
  Mc40[,runs] = 100000*out[,"TB_deaths"]/out[,"Total"]                                        # TB mortality
  Pc40[,runs] = 100000*(out[,"Total_DS"]+out[,"Total_MDR"])/out[,"Total"]                     # TB prevalence  
  LTBIc40[,runs] = 100*(out[,"Total_L"]+out[,"Total_DS"]+out[,"Total_MDR"])/out[,"Total"]     # LTBI
  
}

# Then duplicate each unique run the correct number of times  

I_c <- Ic[,rep(1:dim(Ic)[2],unname(table(t)))]
M_c <- Mc[,rep(1:dim(Mc)[2],unname(table(t)))]
P_c <- Pc[,rep(1:dim(Pc)[2],unname(table(t)))]
LTBI_c <- LTBIc[,rep(1:dim(LTBIc)[2],unname(table(t)))]

I_c20 <- Ic20[,rep(1:dim(Ic20)[2],unname(table(t)))]
M_c20 <- Mc20[,rep(1:dim(Mc20)[2],unname(table(t)))]
P_c20 <- Pc20[,rep(1:dim(Pc20)[2],unname(table(t)))]
LTBI_c20 <- LTBIc20[,rep(1:dim(LTBIc20)[2],unname(table(t)))]

I_c40 <- Ic40[,rep(1:dim(Ic40)[2],unname(table(t)))]
M_c40 <- Mc40[,rep(1:dim(Mc40)[2],unname(table(t)))]
P_c40 <- Pc40[,rep(1:dim(Pc40)[2],unname(table(t)))]
LTBI_c40 <- LTBIc40[,rep(1:dim(LTBIc40)[2],unname(table(t)))]

I_s <- Is[,rep(1:dim(Is)[2],unname(table(t)))]
M_s <- Ms[,rep(1:dim(Ms)[2],unname(table(t)))]
P_s <- Ps[,rep(1:dim(Ps)[2],unname(table(t)))]
LTBI_s <- LTBIs[,rep(1:dim(LTBIs)[2],unname(table(t)))]

## Calculate % reduction in TB incidence and mortality from 2020 level by year and generate boxplot

## By coverage
I_red_30 <- as.data.frame(cbind(100*(I[51,sort(t)] - I_c[56,])/I[51,sort(t)],
               100*(I[51,sort(t)] - I_c[61,])/I[51,sort(t)],
               100*(I[51,sort(t)] - I_c[66,])/I[51,sort(t)]))
colnames(I_red_30) <- c("2025","2030","2035")
I_red_30 <- melt(I_red_30)
I_red_30 <- cbind(I_red_30,"Incidence","30%")
colnames(I_red_30) <- c("Year","Reduction","Var","Coverage")

M_red_30 <- as.data.frame(cbind(100*(M[51,sort(t)] - M_c[56,])/M[51,sort(t)],
               100*(M[51,sort(t)] - M_c[61,])/M[51,sort(t)],
               100*(M[51,sort(t)] - M_c[66,])/M[51,sort(t)]))
colnames(M_red_30) <- c("2025","2030","2035")
M_red_30 <- melt(M_red_30)
M_red_30 <- cbind(M_red_30,"Mortality","30%")
colnames(M_red_30) <- c("Year","Reduction","Var","Coverage")

I_red_20 <- as.data.frame(cbind(100*(I[51,sort(t)] - I_c20[56,])/I[51,sort(t)],
               100*(I[51,sort(t)] - I_c20[61,])/I[51,sort(t)],
               100*(I[51,sort(t)] - I_c20[66,])/I[51,sort(t)]))
colnames(I_red_20) <- c("2025","2030","2035")
I_red_20 <- melt(I_red_20)
I_red_20 <- cbind(I_red_20,"Incidence","20%")
colnames(I_red_20) <- c("Year","Reduction","Var","Coverage")

M_red_20 <- as.data.frame(cbind(100*(M[51,sort(t)] - M_c20[56,])/M[51,sort(t)],
               100*(M[51,sort(t)] - M_c20[61,])/M[51,sort(t)],
               100*(M[51,sort(t)] - M_c20[66,])/M[51,sort(t)]))
colnames(M_red_20) <- c("2025","2030","2035")
M_red_20 <- melt(M_red_20)
M_red_20 <- cbind(M_red_20,"Mortality","20%")
colnames(M_red_20) <- c("Year","Reduction","Var","Coverage")

I_red_40 <- as.data.frame(cbind(100*(I[51,sort(t)] - I_c40[56,])/I[51,sort(t)],
               100*(I[51,sort(t)] - I_c40[61,])/I[51,sort(t)],
               100*(I[51,sort(t)] - I_c40[66,])/I[51,sort(t)]))
colnames(I_red_40) <- c("2025","2030","2035")
I_red_40 <- melt(I_red_40)
I_red_40 <- cbind(I_red_40,"Incidence","40%")
colnames(I_red_40) <- c("Year","Reduction","Var","Coverage")

M_red_40 <- as.data.frame(cbind(100*(M[51,sort(t)] - M_c40[56,])/M[51,sort(t)],
               100*(M[51,sort(t)] - M_c40[61,])/M[51,sort(t)],
               100*(M[51,sort(t)] - M_c40[66,])/M[51,sort(t)]))
colnames(M_red_40) <- c("2025","2030","2035")
M_red_40 <- melt(M_red_40)
M_red_40 <- cbind(M_red_40,"Mortality","40%")
colnames(M_red_40) <- c("Year","Reduction","Var","Coverage")

temp <- rbind(I_red_20,M_red_20,I_red_30,M_red_30,I_red_40,M_red_40)
cov_boxplot <- ggplot(aes(y = Reduction, x = Coverage, fill = Year), data = temp) + geom_boxplot()+
facet_wrap(~Var)+
scale_y_continuous(limits=c(0,25),expand = c(0, 0),"% Reduction from 2020") +
theme_bw()

### By component

I_red_all <- as.data.frame(cbind(100*(I[51,sort(t)] - I_c[56,])/I[51,sort(t)],
                                100*(I[51,sort(t)] - I_c[61,])/I[51,sort(t)],
                                100*(I[51,sort(t)] - I_c[66,])/I[51,sort(t)]))
colnames(I_red_all) <- c("2025","2030","2035")
I_red_all <- melt(I_red_all)
I_red_all <- cbind(I_red_all,"Incidence","COR+ACF")
colnames(I_red_all) <- c("Year","Reduction","Var","Strategy")

M_red_all <- as.data.frame(cbind(100*(M[51,sort(t)] - M_c[56,])/M[51,sort(t)],
                                100*(M[51,sort(t)] - M_c[61,])/M[51,sort(t)],
                                100*(M[51,sort(t)] - M_c[66,])/M[51,sort(t)]))
colnames(M_red_all) <- c("2025","2030","2035")
M_red_all <- melt(M_red_all)
M_red_all <- cbind(M_red_all,"Mortality","COR+ACF")
colnames(M_red_all) <- c("Year","Reduction","Var","Strategy")

I_red_s <- as.data.frame(cbind(100*(I[51,sort(t)] - I_s[56,])/I[51,sort(t)],
                  100*(I[51,sort(t)] - I_s[61,])/I[51,sort(t)],
                  100*(I[51,sort(t)] - I_s[66,])/I[51,sort(t)]))
colnames(I_red_s) <- c("2025","2030","2035")
I_red_s <- melt(I_red_s)
I_red_s <- cbind(I_red_s,"Incidence","ACF")
colnames(I_red_s) <- c("Year","Reduction","Var","Strategy")

M_red_s <- as.data.frame(cbind(100*(M[51,sort(t)] - M_s[56,])/M[51,sort(t)],
                  100*(M[51,sort(t)] - M_s[61,])/M[51,sort(t)],
                  100*(M[51,sort(t)] - M_s[66,])/M[51,sort(t)]))
colnames(M_red_s) <- c("2025","2030","2035")
M_red_s <- melt(M_red_s)
M_red_s <- cbind(M_red_s,"Mortality","ACF")
colnames(M_red_s) <- c("Year","Reduction","Var","Strategy")

temp <- rbind(I_red_all,M_red_all,I_red_s,M_red_s)
comp_boxplot <- ggplot(aes(y = Reduction, x = Strategy, fill = Year), data = temp) + geom_boxplot()+
  facet_wrap(~Var)+
  scale_y_continuous(limits=c(0,20),expand = c(0, 0),"% Reduction from 2020") +
  theme_bw()

# Get quantiles of incidence, prevalence and mortality

# Base case
I_m=apply(I[,t],1,function(x) quantile(x,probs=c(0.5)))
I_l9=apply(I[,t],1,function(x) quantile(x,probs=c(0.025)))
I_u9=apply(I[,t],1,function(x) quantile(x,probs=c(0.975)))

M_m=apply(M[,t],1,function(x) quantile(x,probs=c(0.5)))
M_l9=apply(M[,t],1,function(x) quantile(x,probs=c(0.025)))
M_u9=apply(M[,t],1,function(x) quantile(x,probs=c(0.975)))

P_m=apply(P[,t],1,function(x) quantile(x,probs=c(0.5)))
P_l9=apply(P[,t],1,function(x) quantile(x,probs=c(0.025)))
P_u9=apply(P[,t],1,function(x) quantile(x,probs=c(0.975)))

L_m=apply(LTBI[,t],1,function(x) quantile(x,probs=c(0.5)))
L_l9=apply(LTBI[,t],1,function(x) quantile(x,probs=c(0.025)))
L_u9=apply(LTBI[,t],1,function(x) quantile(x,probs=c(0.975)))

H_m=apply(per_HIV[,t],1,function(x) quantile(x,probs=c(0.5)))
H_l9=apply(per_HIV[,t],1,function(x) quantile(x,probs=c(0.025)))
H_u9=apply(per_HIV[,t],1,function(x) quantile(x,probs=c(0.975)))

# Intervention
I_c_m=apply(I_c,1,function(x) quantile(x,probs=c(0.5)))
I_c_l9=apply(I_c,1,function(x) quantile(x,probs=c(0.025)))
I_c_u9=apply(I_c,1,function(x) quantile(x,probs=c(0.975)))

M_c_m=apply(M_c,1,function(x) quantile(x,probs=c(0.5)))
M_c_l9=apply(M_c,1,function(x) quantile(x,probs=c(0.025)))
M_c_u9=apply(M_c,1,function(x) quantile(x,probs=c(0.975)))

P_c_m=apply(P_c,1,function(x) quantile(x,probs=c(0.5)))
P_c_l9=apply(P_c,1,function(x) quantile(x,probs=c(0.025)))
P_c_u9=apply(P_c,1,function(x) quantile(x,probs=c(0.975)))

L_c_m=apply(LTBI_c,1,function(x) quantile(x,probs=c(0.5)))
L_c_l9=apply(LTBI_c,1,function(x) quantile(x,probs=c(0.025)))
L_c_u9=apply(LTBI_c,1,function(x) quantile(x,probs=c(0.975)))

###################################################
# Plot model outputs vs data

# Convert data to long format
TB_data_plot <- melt(TB_data,id.vars="year")
TB_data_plot <- as.data.frame(cbind(TB_data_plot$year,(do.call(rbind, strsplit(as.character(TB_data_plot$variable), "_"))),TB_data_plot$value))
colnames(TB_data_plot) <- c("Year","var","level","value")
TB_data_plot <- transform(TB_data_plot, var = c("Incidence (/100k)","Mortality (/100k)","Prevalence (/100k)")[as.factor(var)])
TB_data_plot <- droplevels(TB_data_plot)

# Convert model outputs to long format - basecase
TB_model_plot <- as.data.frame(cbind(out[,"time"],I_m,M_m,P_m,L_m,H_m))
colnames(TB_model_plot) <- c("Year","Incidence (/100k)","Mortality (/100k)","Prevalence (/100k)","LTBI (%)","HIV+ TB (%)")
TB_model_plot1 <- cbind(melt(TB_model_plot,id.vars="Year"),"mid")
colnames(TB_model_plot1) <- c("Year","var","value","level")

TB_model_plot <- as.data.frame(cbind(out[,"time"],I_l9,M_l9,P_l9,L_l9,H_l9))
colnames(TB_model_plot) <- c("Year","Incidence (/100k)","Mortality (/100k)","Prevalence (/100k)","LTBI (%)","HIV+ TB (%)")
TB_model_plot2 <- cbind(melt(TB_model_plot,id.vars="Year"),"lo")
colnames(TB_model_plot2) <- c("Year","var","value","level")

TB_model_plot <- as.data.frame(cbind(out[,"time"],I_u9,M_u9,P_u9,L_u9,H_u9))
colnames(TB_model_plot) <- c("Year","Incidence (/100k)","Mortality (/100k)","Prevalence (/100k)","LTBI (%)","HIV+ TB (%)")
TB_model_plot3 <- cbind(melt(TB_model_plot,id.vars="Year"),"hi")
colnames(TB_model_plot3) <- c("Year","var","value","level")

TB_model_plot_base <- rbind(TB_model_plot1,TB_model_plot2,TB_model_plot3)

# Add dummy values to specify scales
TB_model_plot_base <- rbind(TB_model_plot_base,c(0,"Incidence (/100k)",1100,"mid"))
TB_model_plot_base <- rbind(TB_model_plot_base,c(0,"Incidence (/100k)",0,"mid"))
TB_model_plot_base <- rbind(TB_model_plot_base,c(0,"Mortality (/100k)",400,"mid"))
TB_model_plot_base <- rbind(TB_model_plot_base,c(0,"Mortality (/100k)",0,"mid"))
TB_model_plot_base <- rbind(TB_model_plot_base,c(0,"Prevalence (/100k)",1300,"mid"))
TB_model_plot_base <- rbind(TB_model_plot_base,c(0,"Prevalence (/100k)",0,"mid"))
TB_model_plot_base <- rbind(TB_model_plot_base,c(0,"LTBI (%)",100,"mid"))
TB_model_plot_base <- rbind(TB_model_plot_base,c(0,"LTBI (%)",0,"mid"))
TB_model_plot_base <- rbind(TB_model_plot_base,c(0,"HIV+ TB (%)",100,"mid"))
TB_model_plot_base <- rbind(TB_model_plot_base,c(0,"HIV+ TB (%)",0,"mid"))

# Convert model outputs to long format - COR
TB_model_plot <- as.data.frame(cbind(out[,"time"],I_c_m,M_c_m,P_c_m,L_c_m))
colnames(TB_model_plot) <- c("Year","Incidence (/100k)","Mortality (/100k)","Prevalence (/100k)","LTBI (%)")
TB_model_plot1 <- cbind(melt(TB_model_plot,id.vars="Year"),"mid")
colnames(TB_model_plot1) <- c("Year","var","value","level")

TB_model_plot <- as.data.frame(cbind(out[,"time"],I_c_l9,M_c_l9,P_c_l9,L_c_l9))
colnames(TB_model_plot) <- c("Year","Incidence (/100k)","Mortality (/100k)","Prevalence (/100k)","LTBI (%)")
TB_model_plot2 <- cbind(melt(TB_model_plot,id.vars="Year"),"lo")
colnames(TB_model_plot2) <- c("Year","var","value","level")

TB_model_plot <- as.data.frame(cbind(out[,"time"],I_c_u9,M_c_u9,P_c_u9,L_c_u9))
colnames(TB_model_plot) <- c("Year","Incidence (/100k)","Mortality (/100k)","Prevalence (/100k)","LTBI (%)")
TB_model_plot3 <- cbind(melt(TB_model_plot,id.vars="Year"),"hi")
colnames(TB_model_plot3) <- c("Year","var","value","level")

TB_model_plot_COR <- rbind(TB_model_plot1,TB_model_plot2,TB_model_plot3)

# Create plot
fit_plot <- ggplot(TB_data_plot, aes(x=as.numeric(as.character(Year)),y=as.numeric(as.character(value)),linetype=level))+ 
  facet_wrap(~var,scales="free_y")+
  geom_line()+
  geom_line(data=TB_model_plot_base,aes(x=as.numeric(as.character(Year)),y=as.numeric(as.character(value)),linetype=level),colour="red")+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(limits=c(2004,2015),expand = c(0, 0),"Year")+
  scale_linetype_manual(values = c("solid","solid","dashed"))+
  theme_bw() + theme(legend.position = "none")

## Create plot of incidence and mortality trends showing COR trajectories too 
TB_d_plot <- TB_data_plot[TB_data_plot$var %in% c("Incidence (/100k)","Mortality (/100k)"),]
TB_m_plot_base <- TB_model_plot_base[TB_model_plot_base$var %in% c("Incidence (/100k)","Mortality (/100k)"),]
TB_m_plot_COR <- TB_model_plot_COR[TB_model_plot_COR$var %in% c("Incidence (/100k)","Mortality (/100k)"),]

int_plot <- ggplot(TB_d_plot, aes(x=as.numeric(as.character(Year)),y=as.numeric(as.character(value)),linetype=level))+ 
  facet_wrap(~var,scales="free_y")+
  geom_line()+
  geom_line(data=TB_m_plot_base,aes(x=as.numeric(as.character(Year)),y=as.numeric(as.character(value)),linetype=level),colour="red")+
  geom_line(data=TB_m_plot_COR,aes(x=as.numeric(as.character(Year)),y=as.numeric(as.character(value)),linetype=level),colour="blue")+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(limits=c(2004,2035),expand = c(0, 0),"Year")+
  scale_linetype_manual(values = c("solid","solid","dashed"))+
  theme_bw() + theme(legend.position = "none")


