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

##### COR screen and treat intervention
# Screen 20% of HIV- adult population per year
# Avert 0.7*treatment_eff of the disease in that 20% - move these people to latent. zero for children
COR_prevent <- cbind(c(1970,2015,2016),c(0,0,0.5*0.7*0.8))
# Also diagnose 20% of prevalent cases through excluding active disease (symptom screen with assumed sensitivity of 70% - HIV- population, WHO review)
# http://www.who.int/tb/Review2Accuracyofscreeningtests.pdf
COR_diagnose <- cbind(c(1970,2015,2016),c(0,0,0.2*0.7))
 
# Combine forcing functions into a list
force <- list(birth_rate,s_birth,s5,s10,s15,s20,s25,s30,s35,s40,s45,s50,s55,s60,s65,s70,s75,s80,
              h0,h5,h10,h15,h20,h25,h30,h35,h40,h45,h50,h55,h60,h65,h70,h75,h80,
              Ahigh,A500,A349,A249,A199,A99,A50,Athresh,
              BCG_cov,pop_ad,k,dst_n,dst_p,
              COR_prevent,COR_diagnose)

# Fit model using importance resampling approach

# Number of runs to do
n_run = 10
# Years to run model for - always start at 1970
times_run <- seq(1970,2050 , by=1)

# Create arrays for storing outputs and likelihood 
I = mat.or.vec(length(times_run),n_run)
Ih = I
M = I
P = I                                                    

# Calculate likelihood at each time point (1990-2013) for each output
L_inc = mat.or.vec(n_run,24)
L_HIV_inc = L_inc
L_mort = L_inc
L_prev = L_inc
L = rep(0,n_run)

# Sampled parameters
beta <- rep(0,n_run)

for (runs in 1:n_run){
  
  # sample parameters - just beta for now
  
  beta[runs] = runif(1,19,20)
  print(c(runs,beta[runs]))

  # Fitness of MDR, used to calculate parameter for superinfections
  # Both of these are passed into "parms" together with e, the MDR acquisition rate (we set this to zero to exclude MDR in equilibirum phase)

  fit_cost=0.7
  g = fit_cost/(1+fit_cost) # superinfections 
  e = 0.01
  
  # proportion primary (a), proportion smear pos (sig) and mortality rates (muN and mu_I) take different values for 
  # adults (>15) (_a), 0-4 (_0), 5-9 (_5) and 10-14 (_10)
  
  # create parameter vector to pass to model
  parms <- c(age1 = 1/5, age2 = 1/21, beta = 19, 
            a_a = 0.14, a0 = 0.26432, a5 = 0.14056, a10 = 0.056,  
            p = 0.65, v = 0.001, 
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
    
  # Calculate likelihood at each time point for each output
  L_inc[runs,]<-((2*pi*v_inc)^(-1/2))*exp(-(1/(2*v_inc))*((I[,runs][21:44]-TB_data[,"inc_mid"])^2))
  L_mort[runs,]<-((2*pi*v_mort)^(-1/2))*exp(-(1/(2*v_mort))*((M[,runs][21:44]-TB_data[,"mort_mid"])^2))
  L_prev[runs,]<-((2*pi*v_prev)^(-1/2))*exp(-(1/(2*v_prev))*((P[,runs][21:44]-TB_data[,"prev_mid"])^2))
  
  # Total likelihood is based on 2004 and 2013 data points (last 10 years) (product of likelihood for each data point)
  L[runs]<-prod(L_inc[runs,15:24])*prod(L_mort[runs,15:24])*prod(L_prev[runs,15:24])
  
}
  
###################################################
# Plot model outputs vs data

# Convert data to long format
TB_data_plot <- melt(TB_data,id.vars="year")
TB_data_plot <- as.data.frame(cbind(TB_data_plot$year,(do.call(rbind, strsplit(as.character(TB_data_plot$variable), "_"))),TB_data_plot$value))
colnames(TB_data_plot) <- c("Year","var","level","value")
TB_data_plot <- transform(TB_data_plot, var = c("Incidence","Mortality","Prevalence")[as.factor(var)])
TB_data_plot <- droplevels(TB_data_plot)

# Convert best fitting model output to long format
best_fit <- which.max(L) 
TB_model_plot <- as.data.frame(cbind(out[,"time"],I[,best_fit],M[,best_fit],P[,best_fit]))
colnames(TB_model_plot) <- c("Year","Incidence","Mortality","Prevalence")
TB_model_plot <- cbind(melt(TB_model_plot,id.vars="Year"),"mid")
colnames(TB_model_plot) <- c("Year","var","value","level")

# Convert best fitting model output to long format
best_fit <- which.max(L) 
TB_model_plotn <- as.data.frame(cbind(out[,"time"],In[,best_fit],Mn[,best_fit],Pn[,best_fit]))
colnames(TB_model_plotn) <- c("Year","Incidence","Mortality","Prevalence")
TB_model_plotn <- cbind(melt(TB_model_plotn,id.vars="Year"),"mid")
colnames(TB_model_plotn) <- c("Year","var","value","level")

# Convert best fitting model output to long format
best_fit <- which.max(L) 
TB_model_plotc <- as.data.frame(cbind(out[,"time"],Ic[,best_fit],Mc[,best_fit],Pc[,best_fit]))
colnames(TB_model_plotc) <- c("Year","Incidence","Mortality","Prevalence")
TB_model_plotc <- cbind(melt(TB_model_plotc,id.vars="Year"),"mid")
colnames(TB_model_plotc) <- c("Year","var","value","level")

fit_plot <- ggplot(TB_data_plot, aes(x=as.numeric(as.character(Year)),y=as.numeric(as.character(value)),linetype=level))+ 
  facet_wrap(~var,scales="free_y")+
  geom_line()+
  geom_line(data=TB_model_plot,aes(x=as.numeric(as.character(Year)),y=as.numeric(as.character(value))),linetype="dashed",colour="red")+
  geom_line(data=TB_model_plotn,aes(x=as.numeric(as.character(Year)),y=as.numeric(as.character(value))),linetype="solid",colour="red")+
  geom_line(data=TB_model_plotc,aes(x=as.numeric(as.character(Year)),y=as.numeric(as.character(value))),linetype="dotted",colour="red")+
  scale_y_continuous(expand = c(0, 0),"Rate/100,000") +
  scale_x_continuous(limits=c(2000,2050),expand = c(0, 0),"Year")+
  scale_linetype_manual(values = c("solid","solid","dashed"))+
  theme_bw() + theme(legend.position = "none")


# Plot pop against UN data ###################

# convert UN data to long format
temp_data <- melt(UN_pop_age_t,id="Year")
temp_data_l <- melt(UN_pop_age_low_t,id="Year")
temp_data_h <- melt(UN_pop_age_high_t,id="Year")

# sum up model outputs over age groups and turn into long format
tot<-mat.or.vec(81,17)
for(i in 1:17){
  tot[,i] <- apply(out,1,function(x) sum(x[seq(i+1,6410,17)]))
}

temp_model <- as.data.frame(cbind(seq(1970,2050),tot,out[,"Total"]))
colnames(temp_model) <- colnames(UN_pop_age_t)
temp_model <- melt(temp_model,id="Year")

# and plot
plot_pop <- ggplot(temp_model,aes(x=Year,y=value))+
  geom_line(colour="red")+
  geom_line(data=temp_data,aes(x=Year,y=value),colour="black")+
  geom_line(data=temp_data_l,aes(x=Year,y=value),colour="black",linetype="dashed")+
  geom_line(data=temp_data_h,aes(x=Year,y=value),colour="black",linetype="dashed")+
  facet_wrap(~variable,scales="free")+
  xlim(c(1970,2100))

# Resample parameter sets based on likelihood to get a set of good parameters
N_resamp<-200000
t<-sample(seq(1:N_runs),N_resamp,replace=TRUE,prob=L/sum(L))
unique_t<-unique(t)

# Calculate prevalence of LTBI
LTBI <- 100*Infected[,t]/Total[,t]

# Get quantiles of incidence, mort and LTBI
I_m=apply(I[,t],1,function(x) quantile(x,probs=c(0.5)))
I_l9=apply(I[,t],1,function(x) quantile(x,probs=c(0.025)))
I_u9=apply(I[,t],1,function(x) quantile(x,probs=c(0.975)))

Ih_m=apply(Ih[,t],1,function(x) quantile(x,probs=c(0.5)))
Ih_l9=apply(Ih[,t],1,function(x) quantile(x,probs=c(0.025)))
Ih_u9=apply(Ih[,t],1,function(x) quantile(x,probs=c(0.975)))

M_m=apply(M[,t],1,function(x) quantile(x,probs=c(0.5)))
M_l9=apply(M[,t],1,function(x) quantile(x,probs=c(0.025)))
M_u9=apply(M[,t],1,function(x) quantile(x,probs=c(0.975)))

Mh_m=apply(Mh[,t],1,function(x) quantile(x,probs=c(0.5)))
Mh_l9=apply(Mh[,t],1,function(x) quantile(x,probs=c(0.025)))
Mh_u9=apply(Mh[,t],1,function(x) quantile(x,probs=c(0.975)))

LTBI_m=apply(LTBI[,t],1,function(x) quantile(x,probs=c(0.5)))
LTBI_l9=apply(LTBI[,t],1,function(x) quantile(x,probs=c(0.025)))
LTBI_u9=apply(LTBI[,t],1,function(x) quantile(x,probs=c(0.975)))

# Plots

model_data <- as.data.frame(rbind(cbind(seq(1990, 2025, by = 1),I_m,I_l9,I_u9,"Inc","Model"),
                                  cbind(seq(1990, 2025, by = 1),Ih_m,Ih_l9,Ih_u9,"Inc_HIV","Model"),
                                  cbind(seq(1990, 2025, by = 1),M_m,M_l9,M_u9,"Mort","Model"),
                                  cbind(seq(1990, 2025, by = 1),Mh_m,Mh_l9,Mh_u9,"Mort_HIV","Model"),
                                  cbind(TB_incid,"Inc","Data"),
                                  cbind(TB_HIV_incid,"Inc_HIV","Data"),
                                  cbind(TB_mort,"Mort","Data"),
                                  cbind(TB_HIV_mort,"Mort_HIV","Data")))
colnames(model_data) <- c("Time","Rate","Lower","Upper","Var","Type")

fit_plot <- ggplot(model_data, aes(x=as.numeric(as.character(Time)),y=as.numeric(as.character(Rate)),color=Type))+
  # set axis scales and titles
  scale_y_continuous(expand = c(0, 0),"Rate/100,000") +
  scale_x_continuous(limits=c(1990,2015),expand = c(0, 0),"Year") +  
  facet_wrap(~Var,scales="free_y")+
  # line showing median
  geom_line()+
  # ribbon showing 95% CI
  geom_ribbon(aes(ymin=as.numeric(as.character(Lower)), ymax=as.numeric(as.character(Upper)),fill=Type),alpha=0.3,colour=NA)+
  # set theme (to remove grey background, position legend within plot and remove legend title)
  theme_bw() + theme(legend.position = c(0.15,0.9)) + theme(legend.title=element_blank()) + 
  # this line removes the lines from the legend elements created by the geom_line, geom_point and geom_errobar so that only ribbon shows
  # also reverses order so that Placebo comes first (defualt is alphabetical)
  guides(colour=FALSE)+guides(fill=FALSE)+guides(linetype = guide_legend(reverse=TRUE))+
  # this removes the grid line
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

LTBI_data <- as.data.frame(cbind(seq(1990,2025,by=1),LTBI_m,LTBI_l9,LTBI_u9))
colnames(LTBI_data) <- c("Time","Rate","Lower","Upper")

LTBI_plot <- ggplot(LTBI_data, aes(x=as.numeric(as.character(Time)),y=as.numeric(as.character(Rate))))+
  # set axis scales and titles
  scale_y_continuous(limits=c(0,100),expand = c(0, 0),"LTBI prevalence (%)") +
  scale_x_continuous(limits=c(1990,2015),expand = c(0, 0),"Year") +  
  # line showing median
  geom_line()+
  # ribbon showing 95% CI
  geom_ribbon(aes(ymin=as.numeric(as.character(Lower)), ymax=as.numeric(as.character(Upper))),alpha=0.3,colour=NA)+
  # set theme (to remove grey background, position legend within plot and remove legend title)
  theme_bw() + theme(legend.position = c(0.15,0.9)) + theme(legend.title=element_blank()) + 
  # this line removes the lines from the legend elements created by the geom_line, geom_point and geom_errobar so that only ribbon shows
  # also reverses order so that Placebo comes first (defualt is alphabetical)
  guides(colour=FALSE)+guides(fill=FALSE)+guides(linetype = guide_legend(reverse=TRUE))+
  # this removes the grid line
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


 

# Plot TB prevalence ################################
plot(out[,"time"],100*out[,"Total_DS"]/out[,"Total"],ylim=c(0,1))


# Arrange some outputs to take out to excel
cbind(out[,"time"],1000*out[,"Total_DS"],1000*out[,"Total_MDR"]) # Prev in numbers
cbind(out[,"time"],100*out[,"Total_DS"]/out[,"Total"],100*out[,"Total_MDR"]/out[,"Total"],
                   100*(out[,"Total_DS"]+out[,"Total_MDR"])/out[,"Total"]) # Prev in %

cbind(out[,"time"],1000*out[,"Total"]) # Total population

cbind(out[,"time"],100*out[,"TB_deaths"]/out[,"Total"]) # Mortality in %

cbind(out[,"time"],100*out[,"Total_L"]/out[,"Total"]) # LTBI in %

cbind(out[,"time"],1000*out[,"Cases_neg"],1000*out[,"Cases_pos"],1000*out[,"Cases_ART"]) # TB cases


# distribution CD4 no ART
temp <- rbind(out[,"time"],1000*out[,"CD4500"],1000*out[,"CD4350_500"],1000*out[,"CD4250_349"],1000*out[,"CD4200_249"],
                   1000*out[,"CD4100_199"],1000*out[,"CD450_99"],1000*out[,"CD450"])
write.table(temp,file="CD4_no_ART.txt",sep=" ")

# distribution CD4 with ART
temp <- rbind(out[,"time"],1000*out[,"ART500"],1000*out[,"ART350_500"],1000*out[,"ART250_349"],1000*out[,"ART200_249"],
              1000*out[,"ART100_199"],1000*out[,"ART50_99"],1000*out[,"ART50"])
write.table(temp,file="CD4_ART.txt",sep=" ")

# HIV prevalence 15+
cbind(out[,"time"],100*(1000*out[,"CD4500"]+1000*out[,"CD4350_500"]+1000*out[,"CD4250_349"]+1000*out[,"CD4200_249"]+
      1000*out[,"CD4100_199"]+1000*out[,"CD450_99"]+1000*out[,"CD450"]+1000*out[,"ART500"]+1000*out[,"ART350_500"]+
      1000*out[,"ART250_349"]+1000*out[,"ART200_249"]+1000*out[,"ART100_199"]+1000*out[,"ART50_99"]+1000*out[,"ART50"])/(1000*rowSums(tot[,4:17])))

# Number of HIV positives - we currently ignore childhood HIV
cbind(out[,"time"],1000*out[,"CD4500"]+1000*out[,"CD4350_500"]+1000*out[,"CD4250_349"]+1000*out[,"CD4200_249"]+
                          1000*out[,"CD4100_199"]+1000*out[,"CD450_99"]+1000*out[,"CD450"]+1000*out[,"ART500"]+1000*out[,"ART350_500"]+
                          1000*out[,"ART250_349"]+1000*out[,"ART200_249"]+1000*out[,"ART100_199"]+1000*out[,"ART50_99"]+1000*out[,"ART50"])

