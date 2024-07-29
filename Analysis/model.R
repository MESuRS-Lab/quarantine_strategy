#calling packages
library(truncnorm)
library(pse)
library(reshape)
library(gridExtra)
library(grid)
library(abind)
library(matrixStats)
library(dplyr)
library(forcats)
library(metR)
library(scales)
library(ggplot2)
library(odin)
library(dde)

#I MODELS

##I.1 model of the reference scenario
model_generator_reference <- odin::odin({
  #parameters definition
  ##time
  dt <- 1 #day
  update(time) <- (step + 1)*dt
  
  ##fixed parameters
  N_tot <- user(50000) #population of the city considered
  prop_hcw <- user(1.5/100) #health care workers proportion in the population
  delta <- user(1/5.4) #usual hospitals discharge rate
  
  ##biological parameters
  T_E <- user(4.6) #days - average duration of the incubation phase
  T_I <- user(9.8) #days - average duration of the infection (contagious) phase
  phi <- user(0.2) #asymptomatic cases proportion; 1 - phi = severe or mild cases proportion
  psi <- user(0.07) #severe proportion; 1 - psi = asymptomatic or mild cases proportion
  kappa <- user(0.32) #infectiousness reduction factor for asymptomatic subjects
  
  ##sanitary context-dependent parameters
  R_C <- user(1.1) #basic reproduction number in the community
  R_1_ref <- user(2) #basic reproduction number in the usual hospital
  
  ##organizational parameters
  patients_perHCW1 <- user(3) #patients to health care workers ratio in hospitals
  T_shift1 <- user(1) #day - Working shift duration in the usual hospital
  eps <- user(0.2) #factor representing how much  symptomatic HCWs respect self quarantine (0 = they respect it perfectly)
  
  
  #calculation of intermediate parameters
  kappa_mild <- kappa*(phi/(1-psi)) + (1-phi - psi)/(1-psi)  #infectiousness reduction factor for asymptomatic or mild subjects
  
  ##basic reproduction number (p = patient, hcw = health care workers -> phcw = from patients to health care workers)
  ###usual hospital
  R_pp_1 <- R_1_ref 
  R_phcw_1 <- R_1_ref
  R_hcwhcw_1 <- R_1_ref/3 #??because HCWs work 8h a days namely 1/3 of the day
  R_hcwp_1 <- 2*R_1_ref/3 #??
  
  ###community
  R_pp_C <- R_C
  R_phcw_C <- R_C
  R_hcwhcw_C <- R_C
  R_hcwp_C <- R_C
  
  ##transmission rates
  ###affecting HCWs
  beta_hcwhcw_1 <- R_hcwhcw_1/T_I
  beta_phcw_1 <- R_phcw_1/T_I
  beta_hcwhcw_C <- R_hcwhcw_C/T_I
  beta_phcw_C <- R_phcw_C/T_I
  
  ###acffecting patients
  beta_pp_1 <- R_pp_1/T_I
  beta_hcwp_1 <- R_hcwp_1/T_I
  beta_hcwp_C <- R_hcwp_C/T_I
  beta_pp_C <- R_pp_C/T_I
  
  ##infection forces
  ###affecting HCWs
  FI_E_hcw_1 <- beta_hcwhcw_1*kappa*A_hcw_1/N_hcw_1 + beta_phcw_1*(kappa_mild*AM_p_1 + I_p_1)/N_p_1
  FI_E_hcw_C <- beta_hcwhcw_C*((1 - eps)*kappa*A_hcw_C/N_hcw_C +
                                 eps*(kappa*A_hcw_C + IM_hcw_L)/(N_hcw_C + IM_hcw_L)) +
    beta_phcw_C*kappa_mild*AM_p_C/N_p_C
  
  ###acffecting patients
  FI_E_p_1 <- beta_pp_1*(kappa_mild*AM_p_1 + I_p_1)/N_p_1 + beta_hcwp_1*kappa*A_hcw_1/N_hcw_1
  FI_E_p_C <- beta_hcwp_C*((1 - eps)*kappa*A_hcw_C/N_hcw_C + eps*(kappa*A_hcw_C + IM_hcw_L )/(N_hcw_C + IM_hcw_L)) + beta_pp_C*kappa_mild*AM_p_C/N_p_C
  
  ##probability calculation
  
  ###becoming exposed
  ####for HCWs
  p_E_hcw_1 <- 1 - exp(-FI_E_hcw_1*dt)
  p_E_hcw_C <- 1 - exp(-FI_E_hcw_C*dt)
  ####for patients
  p_E_p_1 <- 1 - exp(-FI_E_p_1*dt)
  p_E_p_C <- 1 - exp(-FI_E_p_C*dt)
  ###becoming infectious
  p_IMA <- 1 - exp(-dt/T_E)
  ###recovery
  p_R <- 1 - exp(-dt/T_I)
  
  #compartments upfate
  ##HCWs
  ###HCWs in hospitals
  ####total
  N_hcw_1 <- S_hcw_1 + E_hcw_1 + A_hcw_1 + R_hcw_1
  ####contaminations
  new_E_hcw_1 <- rbinom(S_hcw_1, p_E_hcw_1)### to define
  new_IMA_hcw_1 <- rbinom(E_hcw_1, p_IMA)
  ####symptomatic or asymptomatic
  new_IM_new_A_hcw_1[] <- rmultinom(new_IMA_hcw_1, p_hcw_1)### to define
  p_hcw_1[1] <- 1 - phi
  p_hcw_1[2] <- phi
  dim(p_hcw_1) <- 2
  dim(new_IM_new_A_hcw_1) <- 2
  
  new_IM_hcw_1 <- new_IM_new_A_hcw_1[1]
  new_A_hcw_1 <- new_IM_new_A_hcw_1[2]
  ####recovery
  new_R_hcw_1 <- rbinom(A_hcw_1, p_R)
  ####provisional compartments update
  S_hcw_1_next <- max(S_hcw_1 - new_E_hcw_1,0)
  E_hcw_1_next <- max(E_hcw_1 + new_E_hcw_1 - new_IMA_hcw_1,0)
  A_hcw_1_next <- max( A_hcw_1 + new_A_hcw_1 - new_R_hcw_1,0)
  R_hcw_1_next <- R_hcw_1 + new_R_hcw_1
  
  ###HCWs in community
  ####total
  N_hcw_C <- S_hcw_C + E_hcw_C + A_hcw_C + R_hcw_C
  ####contaminations
  new_E_hcw_C <- rbinom(S_hcw_C, p_E_hcw_C)
  new_IMA_hcw_C <- rbinom(E_hcw_C, p_IMA)
  ####symptomatic or asymptomatic
  new_IM_new_A_hcw_C[] <- rmultinom(new_IMA_hcw_C, p_hcw_C)
  p_hcw_C[1] <- 1-phi
  p_hcw_C[2] <- phi
  dim(p_hcw_C) <- 2
  dim(new_IM_new_A_hcw_C) <- 2
  
  new_IM_hcw_C <- new_IM_new_A_hcw_C[1]
  new_A_hcw_C <- new_IM_new_A_hcw_C[2]
  ####recovery
  new_R_hcw_C <- rbinom(A_hcw_C, p_R)
  ####provisional compartments update
  S_hcw_C_next <- max(S_hcw_C - new_E_hcw_C,0) ## - rebalancing1to2
  E_hcw_C_next <- max(E_hcw_C + new_E_hcw_C - new_IMA_hcw_C,0)
  A_hcw_C_next <- max(A_hcw_C + new_A_hcw_C - new_R_hcw_C,0)
  R_hcw_C_next <- R_hcw_C + new_R_hcw_C +  new_IM_to_R_hcw_L
  
  
  ###symptomatic HCWs confined in the community (L stands for lockdown)
  ####contaminations
  new_IM_to_R_hcw_L <- rbinom(IM_hcw_L,p_R)
  IM_hcw_L_next <- max(IM_hcw_L +  new_IM_hcw_1 +  new_IM_hcw_C - new_IM_to_R_hcw_L,0)
  
  ###staff rotation
  #initial  Cgoinin1_fact <- (N_hcw_working - IM_hcw_L_next)/(N_hcw - N_hcw_working) #normalisation du flux entrant dans H1
  Cgoinin1_fact <- min(1,(N_hcw_working)/(S_hcw_C_next+E_hcw_C_next+A_hcw_C_next+R_hcw_C_next)) #adjust staff to maintain initial hospital capacity
  
  #Cgoinin1_fact <- min(1,(N_hcw_working)/(S_hcw_C_next+E_hcw_C_next+A_hcw_C_next+R_hcw_C_next)) #on veut un nb fixe de patients dans H1
  S_hcw_C_goingin1 <- round((S_hcw_C_next*Cgoinin1_fact - S_hcw_1_next)*dt/T_shift1)
  E_hcw_C_goingin1 <- round((E_hcw_C_next*Cgoinin1_fact - E_hcw_1_next)*dt/T_shift1)
  A_hcw_C_goingin1 <- round((A_hcw_C_next*Cgoinin1_fact - A_hcw_1_next)*dt/T_shift1)
  R_hcw_C_goingin1 <- round((R_hcw_C_next*Cgoinin1_fact - R_hcw_1_next)*dt/T_shift1)
  
  ###final compartments update
  ####symptomatic HCWs confined in the community
  update(IM_hcw_L) <- IM_hcw_L_next
  ####HCWs in the community
  update(S_hcw_C) <- max(S_hcw_C_next - S_hcw_C_goingin1,0)
  update(E_hcw_C) <- max(E_hcw_C_next - E_hcw_C_goingin1,0)
  update(A_hcw_C) <- max(A_hcw_C_next - A_hcw_C_goingin1,0)
  update(R_hcw_C) <- max(R_hcw_C_next - R_hcw_C_goingin1,0)
  ####HCWs in the usual hospital
  update(S_hcw_1) <- max(S_hcw_1_next + S_hcw_C_goingin1,0)
  update(E_hcw_1) <- max(E_hcw_1_next + E_hcw_C_goingin1,0)
  update(A_hcw_1) <- max(A_hcw_1_next + A_hcw_C_goingin1,0)
  update(R_hcw_1) <- max(R_hcw_1_next + R_hcw_C_goingin1,0)
  
  
  ##Patients
  ###patients in the usual hospital
  ####total
  N_p_1 <- S_p_1 + E_p_1 + I_p_1 + AM_p_1 + R_p_1
  ####contaminations
  new_E_p_1 <- rbinom(S_p_1, p_E_p_1)### to define
  new_IMA_p_1 <- rbinom(E_p_1, p_IMA)
  ####severe or not severe (= asymptomatic or mild)
  new_I_new_AM_p_1[] <- rmultinom(new_IMA_p_1, p_p_1)
  p_p_1[1] <- psi
  p_p_1[2] <- 1 - psi
  dim(p_p_1) <- 2
  dim(new_I_new_AM_p_1) <- 2
  
  new_I_p_1 <- new_I_new_AM_p_1[1]
  new_AM_p_1 <- new_I_new_AM_p_1[2]
  ####recovery
  new_AM_to_R_p_1 <- rbinom(AM_p_1, p_R)
  new_I_to_R_p_1 <- rbinom(I_p_1, p_R)
  
  ####provisional compartments update
  S_p_1_next <- max(S_p_1 - new_E_p_1,0)
  E_p_1_next <- max(E_p_1 + new_E_p_1 - new_IMA_p_1,0)
  I_p_1_next <- max(I_p_1 + new_I_p_C + new_I_p_1 - new_I_to_R_p_1,0)
  AM_p_1_next <- max(AM_p_1 + new_AM_p_1  - new_AM_to_R_p_1,0)
  R_p_1_next <- max(R_p_1 + new_AM_to_R_p_1,0)
  
  ###patients in the community
  ####total
  N_p_C <- S_p_C + E_p_C + AM_p_C + R_p_C
  ####contaminations
  new_E_p_C <- rbinom(S_p_C, p_E_p_C)### to define
  new_IMA_p_C <- rbinom(E_p_C, p_IMA)
  ####severe or not severe (= asymptomatic or mild)
  new_I_new_AM_p_C[] <- rmultinom(new_IMA_p_C, p_p_C)### to define
  p_p_C[1] <- psi
  p_p_C[2] <- 1 - psi
  dim(p_p_C) <- 2
  dim(new_I_new_AM_p_C) <- 2
  
  new_I_p_C <- new_I_new_AM_p_C[1]
  new_AM_p_C <- new_I_new_AM_p_C[2]
  ####recovery
  new_AM_to_R_p_C <- rbinom(AM_p_C, p_R)
  ####provisional compartments update
  S_p_C_next <- max(S_p_C - new_E_p_C,0)
  E_p_C_next <- max(E_p_C + new_E_p_C - new_IMA_p_C,0)
  AM_p_C_next <- max(AM_p_C + new_AM_p_C - new_AM_to_R_p_C,0)
  R_p_C_next <- max(R_p_C + new_I_to_R_p_1 +  new_AM_to_R_p_C,0)
  
  ###hospital admission
  #gamma <- delta*(N_p_1_ini - I_p_1_next)/(S_p_C_next+E_p_C_next+AM_p_C_next+R_p_C_next)#mise ? jour taux d'admission pour conserver effectif de l'h?pital
  gamma <- delta*(patients_perHCW1*(S_hcw_1_next + E_hcw_1_next + A_hcw_1_next + R_hcw_1_next)-I_p_1_next)/(S_p_C_next+E_p_C_next+AM_p_C_next+R_p_C_next) #patients are admitted provided there is enough staff 
  #gamma <- delta*N_p_1_ini/(S_p_C_next+E_p_C_next+AM_p_C_next+R_p_C_next)#mise ? jour taux d'admission sans refus de soin
  
  S_p_admitted_in_1 <- round(gamma*S_p_C_next - delta*S_p_1_next)
  E_p_admitted_in_1 <- round(gamma*E_p_C_next - delta*E_p_1_next)
  AM_p_admitted_in_1 <- round(gamma*AM_p_C_next - delta*AM_p_1_next)
  R_p_admitted_in_1 <- round(gamma*R_p_C_next - delta*R_p_1_next)
  
  ###final compartments update
  ####usual hospital
  update(S_p_1) <- max(S_p_1_next + S_p_admitted_in_1,0)
  update(E_p_1) <- max(E_p_1_next + E_p_admitted_in_1,0)
  update(I_p_1) <- I_p_1_next
  update(AM_p_1) <- max(AM_p_1_next + AM_p_admitted_in_1,0)
  update(R_p_1) <- max(R_p_1_next + R_p_admitted_in_1,0)
  ####community
  update(S_p_C) <- max(S_p_C_next - S_p_admitted_in_1,0)
  update(E_p_C) <- max(E_p_C_next - E_p_admitted_in_1,0)
  update(AM_p_C) <- max(AM_p_C_next - AM_p_admitted_in_1,0)
  update(R_p_C) <- max(R_p_C_next - R_p_admitted_in_1,0)
  
  
  #initial condition
  ##time
  initial(time) <- 0
  
  ##initial totals
  N_p <- round((1-prop_hcw)*N_tot)
  
  N_hcw <- N_tot - N_p
  N_hcw_working <- round(N_hcw/3)
  
  N_p_1_ini <- round(N_hcw_working*patients_perHCW1) ##capacity of the usual hospital
  
  ##start of the epidemic
  E_p_ini <- user(100)
  E_p_C_ini <- min(ceiling(((N_p - N_p_1_ini)/N_p)*E_p_ini),E_p_ini-2)
  E_p_1_ini <- E_p_ini - E_p_C_ini
  E_hcw_ini <- max(ceiling(prop_hcw*E_p_ini),1)
  E_hcw_1_ini <- round(E_hcw_ini/3)
  E_hcw_C_ini <- E_hcw_ini - E_hcw_1_ini
  
  ##compartments initialization
  ###HCWs
  ####community
  initial(S_hcw_C) <- N_hcw - N_hcw_working - E_hcw_C_ini
  initial(E_hcw_C) <- E_hcw_C_ini
  initial(A_hcw_C) <- 0
  initial(R_hcw_C) <- 0
  ####usual hospital
  initial(S_hcw_1) <- N_hcw_working - E_hcw_1_ini
  initial(E_hcw_1) <- E_hcw_1_ini
  initial(A_hcw_1) <- 0
  initial(R_hcw_1) <- 0
  ####symptmatic HCWs isolated in the community
  initial(IM_hcw_L) <- 0
  
  ###Patients
  ####community
  initial(S_p_C) <- N_p - N_p_1_ini - E_p_C_ini
  initial(E_p_C) <- E_p_C_ini
  initial(AM_p_C) <- 0
  initial(R_p_C) <- 0
  ####usual hospital
  initial(S_p_1) <- N_p_1_ini - E_p_1_ini
  initial(E_p_1) <- E_p_1_ini
  initial(AM_p_1) <- 0
  initial(I_p_1) <- 0
  initial(R_p_1) <- 0
  
  
  #epidemic indicators
  ##burden and infection peak
  initial(Infections_hcw) <- 0
  initial(max_infections_hcw) <- 0
  initial(Infections_tot) <- 0
  initial(max_infections_tot) <- 0
  initial(cumul_prop_inf_hcw) <- 0
  initial(cumul_prop_inf_p) <- 0
  
  update(Infections_hcw) <- Infections_hcw + new_IMA_hcw_1 + new_IMA_hcw_C
  update(max_infections_hcw) <- max(max_infections_hcw, IM_hcw_L +  A_hcw_C + A_hcw_1)
  update(Infections_tot) <- Infections_tot + new_IMA_p_1 + new_IMA_p_C + new_IMA_hcw_1 + new_IMA_hcw_C
  update(max_infections_tot) <- max(max_infections_tot, AM_p_1 + I_p_1 + AM_p_C + IM_hcw_L +  A_hcw_C + A_hcw_1)
  update(cumul_prop_inf_hcw) <- 100*Infections_hcw/N_hcw
  update(cumul_prop_inf_p) <- 100*(Infections_tot -Infections_hcw )/N_p
  
  ##nosocomial burden
  initial(new_cases_among_p_1) <- E_p_1_ini
  initial(total_hospitalizations) <- N_p_1_ini
  
  update(new_cases_among_p_1) <- new_cases_among_p_1 + new_E_p_1
  update(total_hospitalizations) <- total_hospitalizations + gamma*S_p_C_next
  
  
  
  ##risk for health care workers in the usual hospital
  ###average risk on the whole period
  initial(total_cases_among_hcw_1) <- E_hcw_1_ini
  initial(total_man_day_hcw_1) <- N_hcw_working
  
  update(total_cases_among_hcw_1) <- total_cases_among_hcw_1 + new_E_hcw_1
  update(total_man_day_hcw_1) <- total_man_day_hcw_1 + S_hcw_1
  
  initial(risk_hcw_1) <- 100*E_hcw_1_ini/N_hcw_working
  update(risk_hcw_1) <- 100*total_cases_among_hcw_1/total_man_day_hcw_1
  ### rolling risk on five week-long period
  
  initial(rolling_total_cases_among_hcw_1) <- E_hcw_1_ini
  initial(rolling_total_man_day_hcw_1) <- N_hcw_working
  
  is_accounted_day_hcw_1 <- 1 - ceiling((max(time-1,8)%%35)/35)
  
  update(rolling_total_cases_among_hcw_1) <- (1-is_accounted_day_hcw_1)*(new_E_hcw_1+rolling_total_cases_among_hcw_1) +is_accounted_day_hcw_1*new_E_hcw_1
  ##might be N instead of S
  update(rolling_total_man_day_hcw_1) <- (1-is_accounted_day_hcw_1)*(S_hcw_1+rolling_total_man_day_hcw_1) +is_accounted_day_hcw_1*S_hcw_1
  
  initial(rolling_risk_hcw_1) <-0
  update(rolling_risk_hcw_1) <- 100*rolling_total_cases_among_hcw_1/rolling_total_man_day_hcw_1
  
  
  ##totals to make sure compartments evolution is coherent (mes stands for measured)
  initial(mes_N_hcw_1) <- N_hcw_working
  initial(mes_N_hcw_C) <- N_hcw - N_hcw_working
  initial(mes_N_hcw_L) <- 0
  initial(mes_N_hcw) <- N_hcw
  
  initial(mes_N_p_1) <- N_p_1_ini
  initial(mes_N_p_C) <- N_p - N_p_1_ini
  initial(mes_N_p) <- N_p
  
  update(mes_N_hcw_1) <- N_hcw_1
  update(mes_N_hcw_C) <- N_hcw_C
  update(mes_N_hcw_L) <- IM_hcw_L
  update(mes_N_hcw) <- N_hcw_1 + N_hcw_C + IM_hcw_L
  
  update(mes_N_p_1) <- N_p_1
  update(mes_N_p_C) <- N_p_C
  update(mes_N_p) <- N_p_1 + N_p_C
  
  
  ##compartment propotion
  ###Total
  initial(Prop_inf_hcw) <- 0
  initial(Prop_inf_p) <- 0
  
  ###HCWs
  ####community
  initial(S_hcw_C_prop) <- 100*(N_hcw - N_hcw_working - E_hcw_C_ini)/(N_hcw - N_hcw_working)
  initial(E_hcw_C_prop) <- 100*E_hcw_C_ini/(N_hcw - N_hcw_working)
  initial(A_hcw_C_prop) <- 0
  initial(R_hcw_C_prop) <- 0
  ####usual hospitals
  initial(S_hcw_1_prop) <- 100*(N_hcw_working - E_hcw_1_ini)/N_hcw_working
  initial(E_hcw_1_prop) <- 100*E_hcw_1_ini/N_hcw_working
  initial(A_hcw_1_prop) <- 0
  initial(R_hcw_1_prop) <- 0
  ####isolated symptomatic HCWs
  
  ###Patients
  ####community
  initial(S_p_C_prop) <- 100*(N_p - N_p_1_ini - E_p_C_ini) / (N_p - N_p_1_ini)
  initial(E_p_C_prop) <- 100*E_p_C_ini/ (N_p - N_p_1_ini)
  initial(AM_p_C_prop) <- 0
  initial(R_p_C_prop) <- 0
  ####usual hospitals
  initial(S_p_1_prop) <- 100*(N_p_1_ini - E_p_1_ini)/N_p_1_ini
  initial(E_p_1_prop) <- 100*E_p_1_ini/N_p_1_ini
  initial(AM_p_1_prop) <- 0
  initial(R_p_1_prop) <- 0
  
  ###Infectious people
  update(Prop_inf_hcw) <- 100*(A_hcw_C+A_hcw_1+IM_hcw_L)/N_hcw
  update(Prop_inf_p) <- 100*(AM_p_C+AM_p_1+I_p_1)/N_p
  
  ###HCWs
  ####community
  update(S_hcw_C_prop) <- 100*S_hcw_C/N_hcw_C
  update(E_hcw_C_prop) <- 100*E_hcw_C/N_hcw_C
  update(A_hcw_C_prop) <- 100*A_hcw_C/N_hcw_C
  update(R_hcw_C_prop) <- 100*R_hcw_C/N_hcw_C
  ####usual hospital
  update(S_hcw_1_prop) <- 100*S_hcw_1/N_hcw_1
  update(E_hcw_1_prop) <- 100*E_hcw_1/N_hcw_1
  update(A_hcw_1_prop) <- 100*A_hcw_1/N_hcw_1
  update(R_hcw_1_prop) <- 100*R_hcw_1/N_hcw_1
  ####isolated symptomatic HCWs
  
  ###Patients
  ####community
  update(S_p_C_prop) <- 100*S_p_C/N_p_C
  update(E_p_C_prop) <- 100*E_p_C/N_p_C
  update(AM_p_C_prop) <- 100*AM_p_C/N_p_C
  update(R_p_C_prop) <- 100*R_p_C/N_p_C
  ####usual hospital
  update(S_p_1_prop) <- 100*S_p_1/N_p_1
  update(E_p_1_prop) <- 100*E_p_1/N_p_1
  update(AM_p_1_prop) <- 100*AM_p_1/N_p_1
  update(R_p_1_prop) <- 100*R_p_1/N_p_1
  
  ##acquisition routes
  initial(FI_on_p_1_prop_from_hcw) <- 0
  initial(FI_on_p_1_prop_from_p) <- 0
  initial(FI_on_p_C_prop_from_hcw_C) <- 0
  initial(FI_on_p_C_prop_from_hcw_L) <- 0
  initial(FI_on_p_C_prop_from_p) <- 0
  
  update(FI_on_p_1_prop_from_hcw) <- (beta_hcwp_1*kappa*A_hcw_1/N_hcw_1) /max(.00000001,FI_E_p_1)
  update(FI_on_p_1_prop_from_p) <- (beta_pp_1*(kappa_mild*AM_p_1 + I_p_1)/N_p_1) /max(.00000001,FI_E_p_1  )
  update(FI_on_p_C_prop_from_hcw_C) <- (beta_hcwp_C*((1 - eps)*kappa*A_hcw_C/N_hcw_C + eps*(kappa*A_hcw_C)/(N_hcw_C + IM_hcw_L))) /max(.00000001,FI_E_p_C)
  update(FI_on_p_C_prop_from_hcw_L) <- beta_hcwhcw_C*(eps*(IM_hcw_L)/(N_hcw_C + IM_hcw_L)) /max(.00000001,FI_E_p_C  )
  update(FI_on_p_C_prop_from_p) <-(beta_pp_C*kappa_mild*AM_p_C/N_p_C)/max(.00000001,FI_E_p_C  )
  
}, verbose = FALSE)

##I.2 Model of the quarantine scenario

model_generator_quarantine <- odin::odin({
  #parameters definition
  ##time
  dt <- 1 #day
  update(time) <- (step + 1)*dt
  
  ##fixed parameters
  N_tot <- user(50000) #population of the city considered
  prop_hcw <- user(1.5/100) #health care workers proportion in the population
  delta <- user(1/5.4) #usual hospitals discharge rate
  
  ##biological parameters
  T_E <- user(4.6) #days - average duration of the incubation phase
  T_I <- user(9.8) #days - average duration of the infection (contagious) phase
  phi <- user(0.2) #asymptomatic cases proportion; 1 - phi = severe or mild cases proportion
  psi <- user(0.07) #severe proportion; 1 - psi = asymptomatic or mild cases proportion
  kappa <- user(0.32) #infectiousness reduction factor for asymptomatic subjects
  
  ##sanitary context-dependent
  R_C <- user(1.1) #basic reproduction number in the community
  R_1_quarantine <- user(1.5) #basic reproduction number in the usual hospital
  R_2 <- user(2.5) #basic reproduction number in the quarantine hospital
  
  ##organizational parameters
  patients_perHCW1 <- user(3) #patients to health care workers ratio in hospitals
  T_shift1 <- user(1) #day - Working shift duration in the usual hospital
  patients_perHCW <- patients_perHCW1 #this line exist in case you xant to distinguish patients to hcw ratio between the two hospitals
  T_shift2 <- user(7) #days - Working shift duration in the quarantine hospital
  eps <- user(0.2) #factor representing how much symptomatic HCWs respect self quarantine (0 = they respect it perfectly)
  eps_leavingH2 <- user(0.5) #factor representing how much asymptomatic HCWs leaving the quarantine hospitals respect self quarantine (0 = they respect it perfectly)
  
  #calculation of intermediate parameters
  kappa_mild <- kappa*(phi/(1-psi)) + (1-phi - psi)/(1-psi) #infectiousness reduction factor for asymptomatic or mild subjects
  
  ##basic reproduction number (p = patient, hcw = health care workers -> phcw = from patients to health care workers)
  ###usual hospital
  R_pp_1 <- R_1_quarantine
  R_phcw_1 <- R_1_quarantine
  R_hcwhcw_1 <- R_1_quarantine/3
  R_hcwp_1 <- 2*R_1_quarantine/3
  
  ###community
  R_pp_C <- R_C
  R_phcw_C <- R_C
  R_hcwhcw_C <- R_C
  R_hcwp_C <- R_C
  
  ###quarantine hospitals
  R_phcw_2 <- R_2
  R_hcwhcw_2 <- R_2
  
  ##transmission rates
  ###affecting HCWs
  beta_hcwhcw_1 <- R_hcwhcw_1/T_I
  beta_phcw_1 <- R_phcw_1/T_I
  beta_hcwhcw_2 <- R_hcwhcw_2/T_I
  beta_phcw_2 <- R_phcw_2/T_I
  beta_hcwhcw_C <- R_hcwhcw_C/T_I
  beta_phcw_C <- R_phcw_C/T_I
  
  ###acffecting patients
  beta_pp_1 <- R_pp_1/T_I
  beta_hcwp_1 <- R_hcwp_1/T_I
  beta_hcwp_C <- R_hcwp_C/T_I
  beta_pp_C <- R_pp_C/T_I
  
  ##infection forces
  ###infection sources inside the community
  contacts <- kappa*A_hcw_C/N_hcw_C*(1 - eps)*(1 - eps_leavingH2) +
    (kappa*A_hcw_C + kappa*A_hcw_L)/(N_hcw_C + N_hcw_L - IM_hcw_L)*(1 - eps)*eps_leavingH2 +
    (kappa*A_hcw_C + IM_hcw_L)/(N_hcw_C + IM_hcw_L)*eps*(1 - eps_leavingH2) +
    (kappa*A_hcw_C + kappa*A_hcw_L + IM_hcw_L)/(N_hcw_C + N_hcw_L)*eps*eps_leavingH2
  
  ###affecting HCWs
  FI_E_hcw_1 <- beta_hcwhcw_1*kappa*A_hcw_1/N_hcw_1 + beta_phcw_1*kappa_mild*AM_p_1/N_p_1
  FI_E_hcw_2 <- beta_hcwhcw_2*kappa*A_hcw_2/max(1,N_hcw_2) + beta_phcw_2*I_p_2/(I_p_2+max(1,N_hcw_2))
  FI_E_hcw_C <- beta_hcwhcw_C*contacts + beta_phcw_C*kappa_mild*AM_p_C/N_p_C
  FI_E_hcw_L <- eps*FI_E_hcw_C
  
  ###affectant patients
  FI_E_p_1 <- beta_pp_1*kappa_mild*AM_p_1/N_p_1 + beta_hcwp_1*kappa*A_hcw_1/N_hcw_1
  FI_E_p_C <- beta_hcwp_C*contacts + beta_pp_C*kappa_mild*AM_p_C/N_p_C
  
  ##probability calculation
  
  ###becoming exposed
  ####for HCWs
  p_E_hcw_1 <- 1 - exp(-FI_E_hcw_1*dt)
  p_E_hcw_2 <- 1 - exp(-FI_E_hcw_2*dt)
  p_E_hcw_C <- 1 - exp(-FI_E_hcw_C*dt)
  p_E_hcw_L <- 1 - exp(-FI_E_hcw_L*dt)
  ####for patients
  p_E_p_1 <-1 - exp(-FI_E_p_1*dt)
  p_E_p_C <- 1 - exp(-FI_E_p_C*dt)
  
  ###becoming infectious
  p_IMA <- 1 - exp(-dt/T_E)
  
  ###recovery
  p_R <- 1 - exp(-dt/T_I)
  
  
  #update des compartiments
  
  #compartments upfate
  ##HCWs
  ###HCWs in the usual hospital
  ####total
  N_hcw_1 <- S_hcw_1 + E_hcw_1 + A_hcw_1 + R_hcw_1
  ####contaminations
  new_E_hcw_1 <- rbinom(S_hcw_1, p_E_hcw_1)### to define
  new_IMA_hcw_1 <- rbinom(E_hcw_1, p_IMA)
  ####symptomatic or asymptomatic
  new_IM_new_A_hcw_1[] <- rmultinom(new_IMA_hcw_1, p_hcw_1)### to define
  p_hcw_1[1] <- 1-phi
  p_hcw_1[2] <- phi
  dim(p_hcw_1) <- 2
  dim(new_IM_new_A_hcw_1) <- 2
  
  new_IM_hcw_1 <- new_IM_new_A_hcw_1[1]
  new_A_hcw_1 <- new_IM_new_A_hcw_1[2]
  ####recovery
  new_R_hcw_1 <- rbinom(A_hcw_1, p_R)
  
  ####provisional compartments update
  S_hcw_1_next <- max(S_hcw_1 - new_E_hcw_1,0) ## - rebalancing1to2
  E_hcw_1_next <- max(E_hcw_1 + new_E_hcw_1 - new_IMA_hcw_1,0)
  A_hcw_1_next <- max( A_hcw_1 + new_A_hcw_1 - new_R_hcw_1,0)
  R_hcw_1_next <- R_hcw_1 + new_R_hcw_1
  
  
  ###HCWs in the quarantine hospital
  ####staff
  N_hcw_2 <- S_hcw_2 + E_hcw_2 + A_hcw_2 + R_hcw_2
  ####contaminations
  new_E_hcw_2 <- rbinom(S_hcw_2, p_E_hcw_2)### to define
  new_IMA_hcw_2 <- rbinom(E_hcw_2, p_IMA)
  ####symptomatic or asymptomatic
  new_IM_new_A_hcw_2[] <- rmultinom(new_IMA_hcw_2, p_hcw_2)### to define
  p_hcw_2[1] <- 1-phi
  p_hcw_2[2] <- phi
  dim(p_hcw_2) <- 2
  dim(new_IM_new_A_hcw_2) <- 2
  
  new_IM_hcw_2 <- new_IM_new_A_hcw_2[1]
  new_A_hcw_2 <- new_IM_new_A_hcw_2[2]
  ####recovery
  new_R_hcw_2 <- rbinom(A_hcw_2, p_R)
  ####provisional compartments update
  S_hcw_2_next <- max(S_hcw_2 - new_E_hcw_2,0)
  E_hcw_2_next <- max(E_hcw_2 + new_E_hcw_2 - new_IMA_hcw_2,0)
  A_hcw_2_next <- max(A_hcw_2 + new_A_hcw_2 - new_R_hcw_2,0)
  R_hcw_2_next <- R_hcw_2 + new_R_hcw_2
  N_hcw_2_current <- S_hcw_2_next + E_hcw_2_next + A_hcw_2_next + R_hcw_2_next
  
  ###HCWs in community
  ####total
  N_hcw_C <- S_hcw_C + E_hcw_C + A_hcw_C + R_hcw_C
  ####contaminations
  new_E_hcw_C <- rbinom(S_hcw_C, p_E_hcw_C)### to define
  new_IMA_hcw_C <- rbinom(E_hcw_C, p_IMA)
  ####symptomatic or asymptomatic
  new_IM_new_A_hcw_C[] <- rmultinom(new_IMA_hcw_C, p_hcw_C)### to define
  p_hcw_C[1] <- 1-phi
  p_hcw_C[2] <- phi
  dim(p_hcw_C) <- 2
  dim(new_IM_new_A_hcw_C) <- 2
  
  new_IM_hcw_C <- new_IM_new_A_hcw_C[1]
  new_A_hcw_C <- new_IM_new_A_hcw_C[2]
  ####recovery
  new_R_hcw_C <- rbinom(A_hcw_C, p_R)
  ####provisional compartments update
  S_hcw_C_next <- max(S_hcw_C - new_E_hcw_C,0)
  E_hcw_C_next <- max(E_hcw_C + new_E_hcw_C - new_IMA_hcw_C,0)
  A_hcw_C_next <- max(A_hcw_C + new_A_hcw_C - new_R_hcw_C,0)
  R_hcw_C_next <- R_hcw_C + new_R_hcw_C +  new_IM_to_R_hcw_L
  N_hcw_C_current <- S_hcw_C_next + E_hcw_C_next + A_hcw_C_next + R_hcw_C_next
  
  
  ###symptomatic HCWs confined in the community (L stands for lockdown)
  ####contaminations
  N_hcw_L <- S_hcw_L + E_hcw_L + IM_hcw_L + A_hcw_L + R_hcw_L
  ####contaminations
  new_E_hcw_L <- rbinom(S_hcw_L, p_E_hcw_L)### to define
  new_IMA_hcw_L <- rbinom(E_hcw_L, p_IMA)
  ####symptomatic or asymptomatic
  new_IM_new_A_hcw_L[] <- rmultinom(new_IMA_hcw_L, p_hcw_L)### to define
  p_hcw_L[1] <- 1-phi
  p_hcw_L[2] <- phi
  dim(p_hcw_L) <- 2
  dim(new_IM_new_A_hcw_L) <- 2
  
  new_IM_hcw_L <- new_IM_new_A_hcw_L[1]
  new_A_hcw_L <- new_IM_new_A_hcw_L[2]
  ####recovery
  new_IM_to_R_hcw_L <- rbinom(IM_hcw_L,p_R)
  new_A_to_R_hcw_L <- rbinom(A_hcw_L, p_R)
  ####provisional compartments update
  S_hcw_L_next <- max(S_hcw_L - new_E_hcw_L,0)
  E_hcw_L_next <- max(E_hcw_L + new_E_hcw_L - new_IMA_hcw_L,0)
  IM_hcw_L_next <- max(IM_hcw_L + new_IM_hcw_L + new_IM_hcw_1 + new_IM_hcw_2 + new_IM_hcw_C - new_IM_to_R_hcw_L,0)
  A_hcw_L_next <- max(A_hcw_L + new_A_hcw_L - new_A_to_R_hcw_L,0)
  R_hcw_L_next <- R_hcw_L + new_A_to_R_hcw_L
  
  N_hcw_L_current <- S_hcw_L_next + E_hcw_L_next + IM_hcw_L_next +A_hcw_L_next +R_hcw_L_next
  
  ###staff rotation
  is_shift_day <- 1 - ceiling(max(time,8)%%T_shift2/T_shift2) #quarantine hospital rotation when it is shift day
  
  #Cgoinin1_fact <- (N_hcw_working - N_hcw_L_current - N_hcw_2_current)/(N_hcw - N_hcw_working - N_hcw_2_current)
  Cgoinin1_fact <- min(1,(N_hcw_working-N_hcw_2_current)/(S_hcw_C_next + E_hcw_C_next + R_hcw_C_next)) #attempt to adjust staff to maintain initial hospital capacity
  
  S_hcw_C_goingin1 <- round((S_hcw_C_next*Cgoinin1_fact - S_hcw_1_next)*dt/T_shift1)
  E_hcw_C_goingin1 <- round((E_hcw_C_next*Cgoinin1_fact - E_hcw_1_next)*dt/T_shift1)
  A_hcw_C_goingin1 <- round((A_hcw_C_next*Cgoinin1_fact - A_hcw_1_next)*dt/T_shift1)
  R_hcw_C_goingin1 <- round((R_hcw_C_next*Cgoinin1_fact - R_hcw_1_next)*dt/T_shift1)
  
  #Need_hcw_2 <- max(min(round(ceiling(I_p_2/reorg_step)*reorg_step/patients_perHCW),N_hcw_max),N_hcw_2_min) #besoin d'effectif dans H2
  Need_hcw_2 <- max(min(round(I_p_2/patients_perHCW),N_hcw_max),N_hcw_2_min) #staff requirements in the quarantine hospital
  S_hcw_C_goingin2 <- is_shift_day*round((Need_hcw_2/N_hcw_C_current)*S_hcw_C_next)
  E_hcw_C_goingin2 <- is_shift_day*round((Need_hcw_2/N_hcw_C_current)*E_hcw_C_next)
  A_hcw_C_goingin2 <- is_shift_day*round((Need_hcw_2/N_hcw_C_current)*A_hcw_C_next)
  R_hcw_C_goingin2 <- is_shift_day*round((Need_hcw_2/N_hcw_C_current)*R_hcw_C_next)
  
  
  ###final compartments update
  ####symptomatic HCWs confined in the community
  update(S_hcw_L) <- (1-is_shift_day)*S_hcw_L_next + is_shift_day*S_hcw_2_next
  update(E_hcw_L) <- (1-is_shift_day)*E_hcw_L_next + is_shift_day*E_hcw_2_next
  update(IM_hcw_L) <- IM_hcw_L_next
  update(A_hcw_L) <- (1-is_shift_day)*A_hcw_L_next + is_shift_day*A_hcw_2_next
  update(R_hcw_L) <- (1-is_shift_day)*R_hcw_L_next + is_shift_day*R_hcw_2_next
  ####HCWs in the community
  update(S_hcw_C) <- max(S_hcw_C_next + is_shift_day*S_hcw_L_next - S_hcw_C_goingin2 - S_hcw_C_goingin1,0)
  update(E_hcw_C) <- max(E_hcw_C_next + is_shift_day*E_hcw_L_next - E_hcw_C_goingin2 - E_hcw_C_goingin1,0)
  update(A_hcw_C) <- max(A_hcw_C_next + is_shift_day*A_hcw_L_next - A_hcw_C_goingin2 - A_hcw_C_goingin1 ,0)
  update(R_hcw_C) <- max(R_hcw_C_next + is_shift_day*R_hcw_L_next - R_hcw_C_goingin2 - R_hcw_C_goingin1,0)
  ####HCWs in the usual hospital
  update(S_hcw_1) <- max(S_hcw_1_next  + S_hcw_C_goingin1,0)
  update(E_hcw_1) <- max(E_hcw_1_next  + E_hcw_C_goingin1,0)
  update(A_hcw_1) <- max(A_hcw_1_next  + A_hcw_C_goingin1,0)
  update(R_hcw_1) <- max(R_hcw_1_next  + R_hcw_C_goingin1,0)
  ####HCWs in the quarantine hospital
  update(S_hcw_2) <- (1-is_shift_day)*S_hcw_2_next + S_hcw_C_goingin2
  update(E_hcw_2) <- (1-is_shift_day)*E_hcw_2_next + E_hcw_C_goingin2
  update(A_hcw_2) <- (1-is_shift_day)*A_hcw_2_next + A_hcw_C_goingin2
  update(R_hcw_2) <- (1-is_shift_day)*R_hcw_2_next + R_hcw_C_goingin2
  
  
  
  ##Patients
  ###patients in the usual hospital
  ####total
  N_p_1 <- S_p_1 + E_p_1 + AM_p_1 + R_p_1
  ####contaminations
  new_E_p_1 <- rbinom(S_p_1, p_E_p_1)### to define
  new_IMA_p_1 <- rbinom(E_p_1, p_IMA)
  ####severe or not severe (= asymptomatic or mild)
  new_I_new_AM_p_1[] <- rmultinom(new_IMA_p_1, p_p_1)### to define
  p_p_1[1] <- psi
  p_p_1[2] <- 1 - psi
  dim(p_p_1) <- 2
  dim(new_I_new_AM_p_1) <- 2
  
  new_I_p_1 <- new_I_new_AM_p_1[1]
  new_AM_p_1 <- new_I_new_AM_p_1[2]
  ####recovery
  new_R_p_1 <- rbinom(AM_p_1, p_R)
  ####provisional compartments update 
  S_p_1_next <- max(S_p_1 - new_E_p_1,0)
  E_p_1_next <- max(E_p_1 + new_E_p_1 - new_IMA_p_1,0)
  AM_p_1_next <- max(AM_p_1 + new_AM_p_1  - new_R_p_1,0)
  R_p_1_next <- max(R_p_1 + new_R_p_1,0)
  
  
  ###patients in the quarantine hospitals
  ####recovery
  new_R_p_2 <- rbinom(I_p_2, p_R)
  ####compartment update
  I_p_2_next <- I_p_2 + new_I_p_1 + new_I_p_C - new_R_p_2
  update(I_p_2) <- I_p_2_next
  
  ###patients in the community
  ####total
  N_p_C <- S_p_C + E_p_C + AM_p_C + R_p_C
  ####contaminations
  new_E_p_C <- rbinom(S_p_C, p_E_p_C)
  new_IMA_p_C <- rbinom(E_p_C, p_IMA)
  ####severe or not severe (= asymptomatic or mild)
  new_I_new_AM_p_C[] <- rmultinom(new_IMA_p_C, p_p_C)
  p_p_C[1] <- psi
  p_p_C[2] <- 1 - psi
  dim(p_p_C) <- 2
  dim(new_I_new_AM_p_C) <- 2
  
  new_I_p_C <- new_I_new_AM_p_C[1]
  new_AM_p_C <- new_I_new_AM_p_C[2]
  ####recovery
  new_R_p_C <- rbinom(AM_p_C, p_R)
  ####provisional compartments update
  S_p_C_next <- max(S_p_C - new_E_p_C,0)
  E_p_C_next <- max(E_p_C + new_E_p_C - new_IMA_p_C,0)
  AM_p_C_next <- max(AM_p_C + new_AM_p_C - new_R_p_C,0)
  R_p_C_next <- max(R_p_C + new_R_p_C + new_R_p_2,0)
  
  ###hospital admission
  #gamma <- delta*(N_p_1_ini)/(N_p - N_p_1_ini - I_p_2_next)# initial ( garder fixe le nb de patients)
  gamma <- delta*(patients_perHCW1*(S_hcw_1_next + E_hcw_1_next + A_hcw_1_next + R_hcw_1_next))/(S_p_C_next +E_p_C_next +AM_p_C_next +R_p_C_next)  #patients are admitted according to staff
  #gamma <- delta*(N_p_1_ini)/(S_p_C_next +E_p_C_next +AM_p_C_next +R_p_C_next)# pas de refus de soin
  S_p_admitted_in_1 <- round(gamma*S_p_C_next - delta*S_p_1_next)
  E_p_admitted_in_1 <- round(gamma*E_p_C_next - delta*E_p_1_next)
  AM_p_admitted_in_1 <- round(gamma*AM_p_C_next - delta*AM_p_1_next)
  R_p_admitted_in_1 <- round(gamma*R_p_C_next - delta*R_p_1_next)
  
  ###final compartments update
  ####usual hospital
  update(S_p_1) <- max(S_p_1_next + S_p_admitted_in_1,0)
  update(E_p_1) <- max(E_p_1_next + E_p_admitted_in_1,0)
  update(AM_p_1) <- max(AM_p_1_next + AM_p_admitted_in_1,0)
  update(R_p_1) <- max(R_p_1_next + R_p_admitted_in_1,0)
  ####community
  update(S_p_C) <- max(S_p_C_next - S_p_admitted_in_1,0)
  update(E_p_C) <- max(E_p_C_next - E_p_admitted_in_1,0)
  update(AM_p_C) <- max(AM_p_C_next - AM_p_admitted_in_1,0)
  update(R_p_C) <- max(R_p_C_next - R_p_admitted_in_1,0)
  
  #initial condition
  ##time
  initial(time) <- 0
  
  
  ##initial totals
  N_p <- round((1-prop_hcw)*N_tot)
  
  N_hcw <- N_tot - N_p
  N_hcw_working <- round(N_hcw/3)
  N_hcw_max <- round(N_hcw_working/2)
  N_hcw_2_min <-3# max(ceiling(N_hcw_working*5/100),10) #we expect a minimum number of HCWs in the quarantine regarless of the number of patients
  
  N_p_1_ini <- round(N_hcw_working*patients_perHCW1) #capacity of the usual hospital
  
  ##start of the epidemic
  E_p_ini <- user(100)
  E_p_C_ini <- min(ceiling(((N_p - N_p_1_ini)/N_p)*E_p_ini),E_p_ini-2)
  E_p_1_ini <- E_p_ini - E_p_C_ini
  E_hcw_ini <- max(ceiling(prop_hcw*E_p_ini),1)
  E_hcw_1_ini <- round(E_hcw_ini/3)
  E_hcw_C_ini <- E_hcw_ini - E_hcw_1_ini
  
  ##compartments initialization
  ###HCWs
  ####community
  initial(S_hcw_C) <- N_hcw - N_hcw_working - E_hcw_C_ini
  initial(E_hcw_C) <- E_hcw_C_ini
  initial(A_hcw_C) <- 0
  initial(R_hcw_C) <- 0
  ####usual hospital
  initial(S_hcw_1) <- N_hcw_working - N_hcw_2_min - E_hcw_1_ini
  initial(E_hcw_1) <- E_hcw_1_ini
  initial(A_hcw_1) <- 0
  initial(R_hcw_1) <- 0
  ####quarantine hospital
  initial(S_hcw_2) <- N_hcw_2_min
  initial(E_hcw_2) <- 0
  initial(A_hcw_2) <- 0
  initial(R_hcw_2) <- 0
  ####symptmatic HCWs isolated in the community
  initial(S_hcw_L) <- 0
  initial(E_hcw_L) <- 0
  initial(IM_hcw_L) <- 0
  initial(A_hcw_L) <- 0
  initial(R_hcw_L) <- 0
  
  ###Patients
  ####community
  initial(S_p_C) <- N_p - N_p_1_ini - E_p_C_ini
  initial(E_p_C) <- E_p_C_ini
  initial(AM_p_C) <- 0
  initial(R_p_C) <- 0
  ####usual hospital
  initial(S_p_1) <- N_p_1_ini - E_p_1_ini
  initial(E_p_1) <- E_p_1_ini
  initial(AM_p_1) <- 0
  initial(R_p_1) <- 0
  ####quarantine hospital
  initial(I_p_2) <- 0
  
  
  #epidemic indicators
  ##burden and infection peak
  initial(Infections_hcw) <- 0
  initial(cumul_prop_inf_hcw)<- 0
  initial(cumul_prop_inf_p)<- 0
  initial(max_infections_hcw) <- 0
  initial(Infections_tot) <- 0
  initial(max_infections_tot) <- 0
  
  update(Infections_hcw) <- Infections_hcw + new_IMA_hcw_1 + new_IMA_hcw_2 + new_IMA_hcw_C + new_IMA_hcw_L
  update(cumul_prop_inf_hcw) <- 100*Infections_hcw/N_hcw
  update(max_infections_hcw) <- max(max_infections_hcw, IM_hcw_L + A_hcw_L + A_hcw_C + A_hcw_1 + A_hcw_2)
  update(Infections_tot) <- Infections_tot + new_IMA_p_1 + new_IMA_p_C + new_IMA_hcw_1 + new_IMA_hcw_2 + new_IMA_hcw_C + new_IMA_hcw_L
  update(cumul_prop_inf_p) <- 100*(Infections_tot -Infections_hcw )/N_p
  update(max_infections_tot) <- max(max_infections_tot, AM_p_1 + I_p_2 + AM_p_C + IM_hcw_L + A_hcw_L + A_hcw_C + A_hcw_1 + A_hcw_2)
  
  
  ##totals to make sure compartments evolution is coherent (mes stands for measured)
  update(mes_N_hcw_1) <- N_hcw_1
  update(mes_N_hcw_2) <- N_hcw_2
  update(mes_N_hcw_C) <- N_hcw_C
  update(mes_N_hcw_L) <- N_hcw_L
  update(mes_N_hcw) <- N_hcw_1 + N_hcw_2 + N_hcw_C + N_hcw_L
  
  update(mes_N_p_1) <- N_p_1
  update(mes_N_p_2) <- I_p_2
  update(mes_N_p_C) <- N_p_C
  update(mes_N_p) <- N_p_1 + I_p_2 + N_p_C
  
  initial(mes_N_hcw_1) <- N_hcw_working - N_hcw_2_min
  initial(mes_N_hcw_2) <- N_hcw_2_min
  initial(mes_N_hcw_C) <- N_hcw - N_hcw_working
  initial(mes_N_hcw_L) <- 0
  initial(mes_N_hcw) <- N_hcw
  
  initial(mes_N_p_1) <- N_p_1_ini
  initial(mes_N_p_2) <- 0
  initial(mes_N_p_C) <- N_p - N_p_1_ini
  initial(mes_N_p) <- N_p
  
  ##nosocomial burden
  initial(new_cases_among_p_1) <- E_p_1_ini
  initial(total_hospitalizations) <- N_p_1_ini
  
  update(new_cases_among_p_1) <- new_cases_among_p_1 + new_E_p_1
  update(total_hospitalizations) <- total_hospitalizations + gamma*S_p_C
  
  
  ##risk for health care workers in the usual hospital
  ###average risk over the whole epidemic
  initial(total_cases_among_hcw_1) <- E_hcw_1_ini
  initial(total_man_day_hcw_1) <- N_hcw_working - N_hcw_2_min
  
  update(total_cases_among_hcw_1) <- total_cases_among_hcw_1 + new_E_hcw_1
  update(total_man_day_hcw_1) <- total_man_day_hcw_1 + S_hcw_1
  
  initial(risk_hcw_1) <- 100*E_hcw_1_ini/(N_hcw_working - N_hcw_2_min)
  update(risk_hcw_1) <- 100*total_cases_among_hcw_1/total_man_day_hcw_1
  
  ###rolling risk over five week-long period
  initial(rolling_total_cases_among_hcw_1) <- E_hcw_1_ini
  initial(rolling_total_man_day_hcw_1) <- N_hcw_working - N_hcw_2_min
  
  is_accounted_day_hcw_1 <- 1 - ceiling((max(time-1,8)%%35)/35)
  
  update(rolling_total_cases_among_hcw_1) <- (1-is_accounted_day_hcw_1)*(new_E_hcw_1+rolling_total_cases_among_hcw_1) +is_accounted_day_hcw_1*new_E_hcw_1
  update(rolling_total_man_day_hcw_1) <- (1-is_accounted_day_hcw_1)*(S_hcw_1+rolling_total_man_day_hcw_1) +is_accounted_day_hcw_1*S_hcw_1
  
  initial(rolling_risk_hcw_1) <-0
  update(rolling_risk_hcw_1) <- 100*rolling_total_cases_among_hcw_1/rolling_total_man_day_hcw_1
  
  ##risk for health care workers in the quarantine hospital
  ###average risk over the whole epidemic
  initial(total_cases_among_hcw_2) <- 0
  initial(total_man_day_hcw_2) <- N_hcw_2_min
  
  update(total_cases_among_hcw_2) <- total_cases_among_hcw_2 + new_E_hcw_2
  update(total_man_day_hcw_2) <- total_man_day_hcw_2 + S_hcw_2
  
  initial(risk_hcw_2) <- 0
  update(risk_hcw_2) <- 100*total_cases_among_hcw_2/total_man_day_hcw_2
  
  initial(total_cases_among_hcw_2_shift) <- 0
  initial(total_man_day_hcw_2_shift) <- N_hcw_2_min
  
  ###risk for each shift
  ####for a given shift, risk is updated daily so that the risk calculated
  ####the last day of the shift is equal to the risk over the whole shift
  
  is_accounted_day <- 1 - ceiling((max(time-1,8)%%T_shift2)/T_shift2)
  
  update(total_cases_among_hcw_2_shift) <- (1-is_accounted_day)*(new_E_hcw_2+total_cases_among_hcw_2_shift) +is_accounted_day*new_E_hcw_2
  update(total_man_day_hcw_2_shift) <- (1-is_accounted_day)*total_man_day_hcw_2_shift +is_accounted_day*S_hcw_2
  
  
  initial(risk_hcw_2_shift) <- 0
  update(risk_hcw_2_shift) <- 100*total_cases_among_hcw_2_shift/total_man_day_hcw_2_shift
  
  
  
  ##compartment proportion
  ###Total infections
  initial(Prop_inf_hcw) <- 0
  initial(Prop_inf_p) <- 0
  ###HCWs
  ####community
  initial(S_hcw_C_prop) <- 100*(N_hcw - N_hcw_working - E_hcw_C_ini)/(N_hcw - N_hcw_working)
  initial(E_hcw_C_prop) <- 100*E_hcw_C_ini/(N_hcw - N_hcw_working)
  initial(A_hcw_C_prop) <- 0
  initial(R_hcw_C_prop) <- 0
  ####usual hospital
  initial(S_hcw_1_prop) <- 100*(N_hcw_working - N_hcw_2_min - E_hcw_1_ini)/(N_hcw_working - N_hcw_2_min)
  initial(E_hcw_1_prop) <- 100*E_hcw_1_ini/N_hcw_working
  initial(A_hcw_1_prop) <- 0
  initial(R_hcw_1_prop) <- 0
  ####quarantine hospital
  initial(S_hcw_2_prop) <- 100
  initial(E_hcw_2_prop) <- 0
  initial(A_hcw_2_prop) <- 0
  initial(R_hcw_2_prop) <- 0
  ###Patients
  ####community
  initial(S_p_C_prop) <- 100*(N_p - N_p_1_ini - E_p_C_ini) / (N_p - N_p_1_ini)
  initial(E_p_C_prop) <- 100*E_p_C_ini/ (N_p - N_p_1_ini)
  initial(AM_p_C_prop) <- 0
  initial(R_p_C_prop) <- 0
  ####usual hospital
  initial(S_p_1_prop) <- 100*(N_p_1_ini - E_p_1_ini)/N_p_1_ini
  initial(E_p_1_prop) <- 100*E_p_1_ini/N_p_1_ini
  initial(AM_p_1_prop) <- 0
  initial(R_p_1_prop) <- 0
  
  ###Total
  update(Prop_inf_hcw) <- 100*(A_hcw_C+A_hcw_1+A_hcw_2+A_hcw_L+IM_hcw_L)/N_hcw
  update(Prop_inf_p) <- 100*(AM_p_C+AM_p_1+I_p_2)/N_p 
  ###HCWs
  ####community
  update(S_hcw_C_prop) <- 100*S_hcw_C/N_hcw_C
  update(E_hcw_C_prop) <- 100*E_hcw_C/N_hcw_C
  update(A_hcw_C_prop) <- 100*A_hcw_C/N_hcw_C
  update(R_hcw_C_prop) <- 100*R_hcw_C/N_hcw_C
  ####usual hospital
  update(S_hcw_1_prop) <- 100*S_hcw_1/N_hcw_1
  update(E_hcw_1_prop) <- 100*E_hcw_1/N_hcw_1
  update(A_hcw_1_prop) <- 100*A_hcw_1/N_hcw_1
  update(R_hcw_1_prop) <- 100*R_hcw_1/N_hcw_1
  ####quarantine hospital
  update(S_hcw_2_prop) <- 100*S_hcw_2/max(1,N_hcw_2)
  update(E_hcw_2_prop) <- 100*E_hcw_2/max(1,N_hcw_2)
  update(A_hcw_2_prop) <- 100*A_hcw_2/max(1,N_hcw_2)
  update(R_hcw_2_prop) <- 100*R_hcw_2/max(1,N_hcw_2)
  ###Patients
  ####community
  update(S_p_C_prop) <- 100*S_p_C/N_p_C
  update(E_p_C_prop) <- 100*E_p_C/N_p_C
  update(AM_p_C_prop) <- 100*AM_p_C/N_p_C
  update(R_p_C_prop) <- 100*R_p_C/N_p_C
  ####usual hospital
  update(S_p_1_prop) <- 100*S_p_1/N_p_1
  update(E_p_1_prop) <- 100*E_p_1/N_p_1
  update(AM_p_1_prop) <- 100*AM_p_1/N_p_1
  update(R_p_1_prop) <- 100*R_p_1/N_p_1
  
  ##acquisition routes
  initial(FI_on_p_1_prop_from_hcw) <- 0
  initial(FI_on_p_1_prop_from_p) <- 0
  initial(FI_on_hcw_2_prop_from_hcw) <- 0
  initial(FI_on_hcw_2_prop_from_p) <- 0
  initial(FI_on_p_C_prop_from_hcw_C) <- 0
  initial(FI_on_p_C_prop_from_hcw_L) <- 0
  initial(FI_on_p_C_prop_from_hcw_AL) <- 0
  initial(FI_on_p_C_prop_from_p) <- 0
  
  update(FI_on_p_1_prop_from_hcw) <- (beta_hcwp_1*kappa*A_hcw_1/N_hcw_1) /max(.00000001,FI_E_p_1)
  update(FI_on_p_1_prop_from_p) <- (beta_pp_1*kappa_mild*AM_p_1/N_p_1) /max(.00000001,FI_E_p_1)
  update(FI_on_hcw_2_prop_from_hcw) <- beta_hcwhcw_2*kappa*A_hcw_2/N_hcw_2 /max(.00000001,FI_E_hcw_2)
  update(FI_on_hcw_2_prop_from_p) <- (beta_phcw_2*I_p_2/(I_p_2+N_hcw_2)) /max(.00000001,FI_E_hcw_2)
  update(FI_on_p_C_prop_from_hcw_C) <- beta_hcwp_C*(kappa*A_hcw_C/N_hcw_C*(1 - eps)*(1 - eps_leavingH2) + kappa*A_hcw_C/(N_hcw_C + N_hcw_L - IM_hcw_L)*(1 - eps)*eps_leavingH2 +
                                                      kappa*A_hcw_C/(N_hcw_C + IM_hcw_L)*eps*(1 - eps_leavingH2) + kappa*A_hcw_C/(N_hcw_C + N_hcw_L)*eps*eps_leavingH2
  ) / max(.00000001,FI_E_p_C)
  update(FI_on_p_C_prop_from_hcw_L) <- beta_hcwp_C* (IM_hcw_L/(N_hcw_C + IM_hcw_L)*eps*(1 - eps_leavingH2) +
                                                       IM_hcw_L/(N_hcw_C + N_hcw_L)*eps*eps_leavingH2) / max(.00000001,FI_E_p_C)
  update(FI_on_p_C_prop_from_hcw_AL) <- beta_hcwp_C*( kappa*A_hcw_L/(N_hcw_C + N_hcw_L - IM_hcw_L)*(1 - eps)*eps_leavingH2 +
                                                        kappa*A_hcw_L/(N_hcw_C + N_hcw_L)*eps*eps_leavingH2 ) / max(.00000001,FI_E_p_C)
  update(FI_on_p_C_prop_from_p) <- (beta_pp_C*kappa_mild*AM_p_C/N_p_C)/max(.00000001,FI_E_p_C  )
}, verbose = FALSE)

#II GENERATE SIMULATIONS FOR EPIDEMIC CURVES

##baseline values
paper_pars <- list(T_E = 4.6, T_I = 9.8, phi = 0.2, psi = 0.07,
                   kappa = 0.32 ,R_C = 1.05, R_1_ref = 1.15,
                   R_1_quarantine = 1.1 , R_2 = 1.2, 
                   patients_perHCW = 3,
                   T_shift2 = 7 , eps = 0.15,eps_leavingH2 = 0.4)



date_lim <- 365 #days - simulations ending date
nb_sim <-  1000 #number of simulations

##function generating a given number of simulations of a given scenario up to a given ending date for a given set of parameters
model <- function(type,parameters,date_lim = 300,nb_sim = 1000){
  if (type=="ref"){
    model <- model_generator_reference$new(user=parameters)
  }else{model <- model_generator_quarantine$new(user=parameters)}
  set.seed(1)
  
  res <- model$run(0:date_lim, replicate = nb_sim)
  res_200 <- model$transform_variables(res)[-1]
  
  curve <- lapply(res_200,rowQuantiles, probs = c(0.025,0.5,0.975))
  names(curve) <- paste(names(curve),type, sep = "_")
  low <- lapply(curve, FUN = `[`,TRUE, 1)
  mean <- lapply(curve, FUN = `[`,TRUE, 2)
  up <- lapply(curve, FUN = `[`,TRUE, 3)
  
  low <- do.call(cbind, low)
  mean <- do.call(cbind, mean)
  up <- do.call(cbind, up)
  
  df_low <- as.data.frame(low)
  df_mean <- as.data.frame(mean)
  df_up <-as.data.frame(up)
  return(list(df_low,df_mean,df_up))
}






