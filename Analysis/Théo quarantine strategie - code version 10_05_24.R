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



###build the dataset containing simulations according to the baseline values, chosen ending date and number of simulations
df_all <- model("ref",paper_pars, date_lim,nb_sim)
df_quarantine <- model("quarantine",paper_pars, date_lim,nb_sim)
for (i in 1:3){
  df_all[[i]] <- cbind(df_all[[i]],df_quarantine[[i]])
}



###function for plotting epidemic curves
plot_curves <- function(df,columns,curve_colors,title, labels, ylab = "Proportion (%)"){
  n <- length(df[[1]]$time_ref)
  df_for_plot <- NULL
  for (col in columns){
    df_col <- data.frame(date = df[[1]]$time_ref, nb_5p = df[[1]][[col]], nb = df[[2]][[col]], nb_95p = df[[3]][[col]], comp = rep(col, n))
    df_for_plot <- rbind(df_for_plot,df_col)
  }
  df_for_plot <- df_for_plot %>% arrange(factor(comp, levels = columns))
  
  ggplot(data=df_for_plot, aes(x=date, y=nb, ymin=nb_5p, ymax=nb_95p, fill=comp, linetype = comp)) +
    geom_line() +
    geom_ribbon(alpha=0.5) +
    xlab("Time (days)") +  guides(linetype=FALSE) +
    theme(legend.position="top") + #c(0.8, 0.5)) +
    scale_fill_manual(values =curve_colors, name = "", limits = columns,labels=labels)+
    ylab(ylab)+ theme(text = element_text(size = 20)) + theme(plot.title=element_text(hjust=0.5))
  #+ggtitle(title)
}
####useful color palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

####usual epidemic curves
burden_reduction_prop_hcw <- c(paste(c("cumul_prop_inf_hcw"),"ref",sep="_"),paste("cumul_prop_inf_hcw","quarantine",sep="_"))
plot_curves(df_all, burden_reduction_prop_hcw, c("blue","red"), "Cumulative proportion of infections among HCWs",
            c("Reference strategy","Quarantine strategy"),"Cumulative proportion of infections among HCWs (%)")
burden_reduction_prop_non_hcw <- c(paste(c("cumul_prop_inf_p"),"ref",sep="_"),paste("cumul_prop_inf_p","quarantine",sep="_"))
plot_curves(df_all, burden_reduction_prop_non_hcw, c("blue","red"), "Cumulative proportion of infections among non HCWs",
            c("Reference strategy","Quarantine strategy"),"Cumulative proportion of infections among non HCWs (%)")

prop_infection_hcw <- c(paste(c("Prop_inf_hcw"),"ref",sep="_"),paste("Prop_inf_hcw","quarantine",sep="_"))
plot_curves(df_all, prop_infection_hcw, c("blue","red"), "Proportion of infectious individuals",
            c("Reference strategy","Quarantine strategy"),"Proportion of infectious HCWs (%)")
prop_infection_non_hcw <- c(paste(c("Prop_inf_p"),"ref",sep="_"),paste("Prop_inf_p","quarantine",sep="_"))
plot_curves(df_all, prop_infection_non_hcw, c("blue","red"), "Proportion of infectious individuals",
            c("Reference strategy","Quarantine strategy"),"Proportion of infectious non HCWs (%)")

virus_exposition <- c(paste(c("S_p_C_prop"),"ref",sep="_"),paste("S_p_C_prop","quarantine",sep="_"))
plot_curves(df_all, virus_exposition, c("blue","red"),  "Susceptible patients among community",c("Reference strategy","Quarantine strategy"),"Proportion of susceptible non HCWs among community (%)")

FI_on_p_C_evo_ref <- paste(c("FI_on_p_C_prop_from_p","FI_on_p_C_prop_from_hcw_C","FI_on_p_C_prop_from_hcw_L"),"ref",sep="_") #,"FI_on_p_C_prop_from_hcw_AL"
plot_curves(df_all, FI_on_p_C_evo_ref, c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
            "Virus acquisition routes among community without quarantine strategy",
            c( "From patients","From resting HCWs", "From strictly\nisolated HCWs"),
            "Acquisition route frequency")# with reference strategy")

FI_on_p_C_evo_quar <- paste(c("FI_on_p_C_prop_from_p","FI_on_p_C_prop_from_hcw_C","FI_on_p_C_prop_from_hcw_L","FI_on_p_C_prop_from_hcw_AL"),"quarantine",sep="_") #,"FI_on_p_C_prop_from_hcw_AL"
plot_curves(df_all, FI_on_p_C_evo_quar, c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7"),
            "Virus acquisition routes among community with quarantine strategy",
            c("From patients","From resting HCWs",  "From strictly\nisolated HCWs",  "From HCWs isolated after working\nin the quarantine hospital"),
            "Acquisition route frequency")# with quarantine strategy")

Nb_of_patients_ref <- paste(c( "mes_N_p_1", "I_p_1"),"ref",sep="_") 
plot_curves(df_all, Nb_of_patients_ref, c( "black","red"),
            "",
            c( "In usual hospital","Of which symptomatic patients"),
            "Number of patients")

Nb_of_patients_quarantine <- paste(c("mes_N_p_1","mes_N_p_2"),"quarantine",sep="_") 
plot_curves(df_all, Nb_of_patients_quarantine, c( "black", "red"),
            "Number of hospitalized patients with quarantine strategy",
            c( "In usual hospital", "In quarantine hospital"),
            "Number of patients")

Nb_of_hcw_hospitals <- c("mes_N_hcw_1_ref","mes_N_hcw_1_quarantine","mes_N_hcw_2_quarantine")
plot_curves(df_all, Nb_of_hcw_hospitals, c( "black", "grey","red"),
            "",
            c( "In usual hospital\nno quarantine","In usual hospital\nwith quarantine" ,"In quarantine hospital"),
            "Number of hcws")
Nb_of_hcw_C <- c("mes_N_hcw_L_ref","mes_N_hcw_L_quarantine")
plot_curves(df_all, Nb_of_hcw_C, c( "black", "blue"),
            "Number of  hcw with quarantine strategy",
            c( "Reference strategy","Quarantine strategy"),
            "Number of hcw in the community" )

#Risk_for_hcw <- c("risk_hcw_1_ref","risk_hcw_1_quarantine","risk_hcw_2_shift_quarantine") 
#plot_curves(df_all, Risk_for_hcw, c( "black", "grey","red"),
#            "Risk for hcw",
#            c( "In usual hospital\nwith quarantine", "In usual hospital\nwithout quarantine", "In quarantine hospital (shift)"),
#            "Incidence (%)")
#Risk_for_hcw_mean <- c("risk_hcw_1_ref","risk_hcw_1_quarantine","risk_hcw_2_quarantine") 
#plot_curves(df_all, Risk_for_hcw_mean, c( "black", "grey","red"),
#            "Risk for hcw",
#            c( "In usual hospital\nwith quarantine", "In usual hospital\nwithout quarantine", "In quarantine hospital"),
#            "Incidence (%)")

df_risk <- data.frame(date = 1:42, nb_5p = df_all[[1]]$risk_hcw_2_shift_quarantine[7*(1:42)+3]
                      , nb = df_all[[2]]$risk_hcw_2_shift_quarantine[7*(1:42)+3]
                      , nb_95p = df_all[[3]]$risk_hcw_2_shift_quarantine[7*(1:42)+3])


ggplot(data=df_risk, aes(x=date, y=nb, ymin=nb_5p, ymax=nb_95p)) +
  geom_point() +
  geom_ribbon(alpha=0.5) +
  xlab("Number of weeks") +  guides(linetype=FALSE) +
  theme(legend.position="top") + #c(0.8, 0.5)) +
 
  ylab("7-days shift risk (%)")+ theme(text = element_text(size = 20)) + theme(plot.title=element_text(hjust=0.5))



#III ANALYSIS

###print the table with risk for HCWs values

####function calcultating the epidemic peak week index
maximazing_week <- function(vector){
  n_week <- length(vector)/7
  weekly_mean <- mean(vector[1:7])
  for (i in 1:n_week-1){
    weekly_mean <- c(weekly_mean, mean(vector[(7*i+1):(7*(i+1))]))
  }
  return(median(which(weekly_mean == max(weekly_mean))))
}

#### calculating risk with prediction band
peak_weak_quarantine_strategy1 <- maximazing_week(df_all[[1]]$Prop_inf_p_quarantine)
peak_weak_quarantine_strategy2 <- maximazing_week(df_all[[2]]$Prop_inf_p_quarantine)
peak_weak_quarantine_strategy3 <- maximazing_week(df_all[[3]]$Prop_inf_p_quarantine)


med_risk_shift_2_5 <- mean(df_all[[1]]$risk_hcw_2_shift_quarantine[(7*((peak_weak_quarantine_strategy1-3):(peak_weak_quarantine_strategy1+1)))+3])
med_risk_shift_2_50 <- mean(df_all[[2]]$risk_hcw_2_shift_quarantine[(7*((peak_weak_quarantine_strategy2-3):(peak_weak_quarantine_strategy2+1)))+3])
med_risk_shift_2_95 <- mean(df_all[[3]]$risk_hcw_2_shift_quarantine[(7*((peak_weak_quarantine_strategy3-3):(peak_weak_quarantine_strategy3+1)))+3])

med_risk_hcw_1_quar_5 <- mean(df_all[[1]]$rolling_risk_hcw_1_quarantine[((7*(peak_weak_quarantine_strategy1-3)):(7*(peak_weak_quarantine_strategy1+2)-1))])
med_risk_hcw_1_quar_50 <- mean(df_all[[2]]$rolling_risk_hcw_1_quarantine[((7*(peak_weak_quarantine_strategy2-3)):(7*(peak_weak_quarantine_strategy2+2)-1))])
med_risk_hcw_1_quar_95 <- mean(df_all[[3]]$rolling_risk_hcw_1_quarantine[((7*(peak_weak_quarantine_strategy3-3)):(7*(peak_weak_quarantine_strategy3+2)-1))])

peak_weak_ref_strategy1 <- maximazing_week(df_all[[1]]$Prop_inf_p_ref)
peak_weak_ref_strategy2 <- maximazing_week(df_all[[2]]$Prop_inf_p_ref)
peak_weak_ref_strategy3 <- maximazing_week(df_all[[3]]$Prop_inf_p_ref)

med_risk_hcw_1_ref_5 <- mean(df_all[[1]]$rolling_risk_hcw_1_ref[((7*(peak_weak_ref_strategy1-3)):(7*(peak_weak_ref_strategy1+2)-1))])
med_risk_hcw_1_ref_50 <- mean(df_all[[2]]$rolling_risk_hcw_1_ref[((7*(peak_weak_ref_strategy2-3)):(7*(peak_weak_ref_strategy2+2)-1))])
med_risk_hcw_1_ref_95 <- mean(df_all[[3]]$rolling_risk_hcw_1_ref[((7*(peak_weak_ref_strategy3-3)):(7*(peak_weak_ref_strategy3+2)-1))])


#### gathering results and printing them
risk_hcw_1_ref <- c(med_risk_hcw_1_ref_5 ,med_risk_hcw_1_ref_50 ,med_risk_hcw_1_ref_95 )
risk_hcw_1_quarantine <- c(med_risk_hcw_1_quar_5 ,med_risk_hcw_1_quar_50 ,med_risk_hcw_1_quar_95)
weekly_risk_hcw_2_quarantine <- c(med_risk_shift_2_5,med_risk_shift_2_50,med_risk_shift_2_95)
risk <- t(matrix(c(risk_hcw_1_ref,risk_hcw_1_quarantine,weekly_risk_hcw_2_quarantine), ncol = 3))
rownames(risk) <- c("risk_hcw_1_ref","risk_hcw_1_quarantine","weekly_risk_hcw_2_quarantine")
colnames(risk) <- c( "2.5%", "50%","97.5%")
print(risk)

##indicators for the heatmap at the end of the paper

###function calculating epidemic indicators on one simulation given one set of parameters
oneRun <- function(pars){
  model <- model_generator_quarantine$new(user = pars)
  
  res <- model$run(0:date_lim)
  res_matrix <- model$transform_variables(res)[-1]
  infections_tot <- res_matrix$Infections_tot[date_lim]
  max_infections_tot <- res_matrix$max_infections_tot[date_lim]
  infections_hcw <- res_matrix$Infections_hcw[date_lim]
  max_infections_hcw <- res_matrix$max_infections_hcw[date_lim]
  
  nosocomial_burden <- res_matrix$new_cases_among_p_1[date_lim]
  total_hospitalizations <- res_matrix$total_hospitalizations[date_lim]
  
  
  model_reference <- model_generator_reference$new(user = pars)
  
  res_reference <- model_reference$run(0:date_lim)
  res_matrix_reference <- model_reference$transform_variables(res_reference)[-1]
  infections_tot_reference <- max(res_matrix_reference$Infections_tot[date_lim],1)
  max_infections_tot_reference<- max(res_matrix_reference$max_infections_tot[date_lim],1)
  infections_hcw_reference <- max(res_matrix_reference$Infections_hcw[date_lim],1)
  max_infections_hcw_reference<- max(res_matrix_reference$max_infections_hcw[date_lim],1)
  
  nosocomial_burden_reference <- max(res_matrix_reference$new_cases_among_p_1[date_lim],1)
  total_hospitalizations_reference <- max(res_matrix_reference$total_hospitalizations[date_lim],1)
  
  
  infections_tot_reduction <- -round(100*(infections_tot - infections_tot_reference) / infections_tot_reference, digits = 3)
  max_infections_tot_reduction <- -round(100*(max_infections_tot - max_infections_tot_reference) / max_infections_tot_reference, digits = 3)
  infections_hcw_reduction <- -round(100*(infections_hcw - infections_hcw_reference) / infections_hcw_reference, digits = 3)
  max_infections_hcw_reduction <- -round(100*(max_infections_hcw - max_infections_hcw_reference) / max_infections_hcw_reference, digits = 3)
  nosocomial_burden_reduction <- -round(100*(nosocomial_burden/total_hospitalizations - nosocomial_burden_reference/total_hospitalizations_reference) / (nosocomial_burden_reference/total_hospitalizations_reference), digits = 3)
  
  indicators <- c(max(-50,as.numeric(infections_tot_reduction)), max(-50,as.numeric(max_infections_tot_reduction)),
                  max(-50,as.numeric(infections_hcw_reduction)),max(-50,as.numeric(max_infections_hcw_reduction)),
                  max(-50,as.numeric(nosocomial_burden_reduction)))
  return(indicators)
}


###function doing the same but adding absolute indicators and paving the way for prediction bands far all indicators
oneRun_withPB <- function(pars){
  model <- model_generator_quarantine$new(user = pars)
  
  res <- model$run(0:date_lim)
  res_matrix <- model$transform_variables(res)[-1]
  
  model_reference <- model_generator_reference$new(user = pars)
  res_reference <- model_reference$run(0:date_lim)
  res_matrix_reference <- model_reference$transform_variables(res_reference)[-1]
   
 
  infections_tot_reduction <- -round(100*(res_matrix$Infections_tot[date_lim] - max(res_matrix_reference$Infections_tot[date_lim],1)) / max(res_matrix_reference$Infections_tot[date_lim],1), digits = 3)
  max_infections_tot_reduction <- -round(100*(res_matrix$max_infections_tot[date_lim] - max(res_matrix_reference$max_infections_tot[date_lim],1)) / max(res_matrix_reference$max_infections_tot[date_lim],1), digits = 3)
  infections_hcw_reduction <- -round(100*(res_matrix$Infections_hcw[date_lim] - max(res_matrix_reference$Infections_hcw[date_lim],1)) / max(res_matrix_reference$Infections_hcw[date_lim],1), digits = 3)
  max_infections_hcw_reduction <- -round(100*(res_matrix$max_infections_hcw[date_lim] - max(res_matrix_reference$max_infections_hcw[date_lim],1)) / max(res_matrix_reference$max_infections_hcw[date_lim],1), digits = 3)
  nosocomial_burden_reduction <- -round(100*(max(res_matrix$new_cases_among_p_1[date_lim],1)/max(res_matrix$total_hospitalizations[date_lim],1) - max(res_matrix_reference$new_cases_among_p_1[date_lim],1)/max(res_matrix_reference$total_hospitalizations[date_lim],1)) / (max(res_matrix_reference$new_cases_among_p_1[date_lim],1)/max(res_matrix_reference$total_hospitalizations[date_lim],1)), digits = 3)
 
  return(list(infections_tot_ref= max(res_matrix_reference$Infections_tot[date_lim],1),
              
                max_infections_tot_ref= max(res_matrix_reference$max_infections_tot[date_lim],1),
              
                infections_hcw_ref= max(res_matrix_reference$Infections_hcw[date_lim],1),
              
                max_infections_hcw_ref= max(res_matrix_reference$max_infections_hcw[date_lim],1),
              
                nosocomial_burden_ref= max(res_matrix_reference$new_cases_among_p_1[date_lim],1),
              
                total_hospitalizations_ref= max(res_matrix_reference$total_hospitalizations[date_lim],1),
              
                infections_tot_quar= res_matrix$Infections_tot[date_lim],
              
                max_infections_tot_quar= res_matrix$max_infections_tot[date_lim],
              
                infections_hcw_quar= res_matrix$Infections_hcw[date_lim],
              
                max_infections_hcw_quar= res_matrix$max_infections_hcw[date_lim],
              
                nosocomial_burden_quar= res_matrix$new_cases_among_p_1[date_lim],
              
                total_hospitalizations_quar= res_matrix$total_hospitalizations[date_lim],
              
                risk_p_1_ref= 100*max(res_matrix_reference$new_cases_among_p_1[date_lim],1)/max(res_matrix_reference$total_hospitalizations[date_lim],1),
              
                risk_p_1_quar= 100*max(res_matrix$new_cases_among_p_1[date_lim],1)/max(res_matrix$total_hospitalizations[date_lim],1),
              
                risk_hcw_1_ref= res_matrix_reference$risk_hcw_1[date_lim],
              
                risk_hcw_1_quar=  res_matrix$risk_hcw_1[date_lim],
              
                risk_hcw_2_quar= res_matrix$risk_hcw_2[date_lim],
              
                infections_reduc= max(-50,as.numeric(infections_tot_reduction)),
              
                infection_peak_reduc= max(-50,as.numeric(max_infections_tot_reduction)),
              
                HCW_infections_reduc= max(-50,as.numeric(infections_hcw_reduction)),
              
                HCW_infection_peak_reduc= max(-50,as.numeric(max_infections_hcw_reduction)),
              
                nosocomial_burden_reduc= max(-50,as.numeric(nosocomial_burden_reduction))))
  
  #return(c(values, indicators))
}


###function calculating median epidemic indicators of a given number of simulations given one set of parameters
severalRun <- function(n_sim,pars){
  set.seed(1)
  mat <- matrix(oneRun(pars), ncol = 5, dimnames = list(NULL,indicators_name))
  for (l in 1:n_sim-1){
    mat <- rbind(mat,oneRun(pars))
  }
  return(colQuantiles(mat, probs = 0.5))
  #mean <- lapply(mat,rowQuantiles, probs = 0.5)
  #return(mean)
}

###function doing the same but adding absolute indicators and prediction bands far all indicators
severalRun_withPB<- function(n_sim,pars){
  set.seed(1)
  row_names <- names(oneRun_withPB(pars))
  mat <- matrix(as.numeric(format(unlist(oneRun_withPB(pars)), digits = 2)), ncol = 22, dimnames = list(NULL,row_names))
  for (l in 1:n_sim-1){
    mat <- rbind(mat,as.numeric(format(unlist(oneRun_withPB(pars)), digits = 2, scientific = FALSE)))
  }
  return(colQuantiles(mat, probs = c(0.025,0.5,0.975)))
  
}

###call the latter function to display indicators for the paper
severalRun_withPB(100, paper_pars)

##########
##########


##Sensitivity analysis
###parameter space explored
#distrib <- list(T_E = c("T_E","qtruncnorm"),
                #T_E_arg = list(a=0, b=Inf, mean=5, sd=2),
                #T_I = c("T_I","qtruncnorm"),
                #T_I_arg = list(a=0, b=Inf, mean=9.5, sd=2),
                #phi = c("phi","qtruncnorm"),
                #phi_arg = list(a=0, b=1, mean=0.2, sd=0.02),
                #psi = c("psi","qtruncnorm"),
                #psi_arg = list(a=0, b=1, mean=0.07, sd=0.02),
                #kappa = c("kappa","qtruncnorm"),
                #kappa_arg = list(a=0, b=1, mean=0.32, sd=0.07),
                #R_C = c("R_C","qtruncnorm"),
                #R_C_arg = list(a=0, b=Inf, mean=1.1, sd=0.4),
                #R_1_ref = c("R_1_ref","qtruncnorm"),
                #R_1_ref_arg = list( a=0, b=Inf, mean=2, sd=0.3),
                #R_1_quarantine = c("R_1_quarantine","qtruncnorm"),
                #R_1_quarantine = list(a=0, b=Inf, mean=1.5, sd=0.3),
                #R_2 = c("R_2","qtruncnorm"),
                #R_2_arg = list(a=0, b=Inf, mean=2.5, sd=1),
                #T_shift2 = c("T_shift2","qunif"),
                #T_shift2_arg = list(min=7, max=14),
                #eps = c("eps","qunif"),
                #eps_arg = list(min=0, max=0.4),
                #eps_leavingH2 = c("eps_leavingH2","qunif"),
                #eps_leavingH2_arg = list(min=0.4, max=0.7))
#n_pars = length(distrib)/2

distrib <- list(T_E = c("T_E","qtruncnorm"),
                T_E_arg = list(a=0, b=Inf, mean=5, sd=2),
                T_I = c("T_I","qtruncnorm"),
                T_I_arg = list(a=0, b=Inf, mean=9.5, sd=2),
                phi = c("phi","qtruncnorm"),
                phi_arg = list(a=0, b=1, mean=0.2, sd=0.02),
                psi = c("psi","qtruncnorm"),
                psi_arg = list(a=0, b=1, mean=0.07, sd=0.02),
                kappa = c("kappa","qtruncnorm"),
                kappa_arg = list(a=0, b=1, mean=0.32, sd=0.07),
                R_C = c("R_C","qtruncnorm"),
                R_C_arg = list(a=0, b=Inf, mean=1.05, sd=0.15),
                R_1_ref = c("R_1_ref","qtruncnorm"),
                R_1_ref_arg = list( a=0, b=Inf, mean=1.15, sd=0.1),
                R_1_quarantine = c("R_1_quarantine","qtruncnorm"),
                R_1_quarantine_arg = list(a=0, b=Inf, mean=1.1, sd=0.1),
                R_2 = c("R_2","qtruncnorm"),
                R_2_arg = list(a=0, b=Inf, mean=1.2, sd=0.15),
                T_shift2 = c("T_shift2","qunif"),
                T_shift2_arg = list(min=7, max=14),
                eps = c("eps","qunif"),
                eps_arg = list(min=0, max=0.3),
                eps_leavingH2 = c("eps_leavingH2","qunif"),
                eps_leavingH2_arg = list(min=0.2, max=0.6))
n_pars = length(distrib)/2
###list allowing for maths symbols in the PRCC graph legend
name_pars_maths <- list(T_E = "T[E]",
                T_I = "T[I]",
                phi = "phi",
                psi = "psi","",
                kappa = "kappa",
                R_C = "R[C]",
                R_1_ref = "R[1]^ref",
                R_1_quarantine = "R[1]^quarantine",
                R_2 = "R[2]",
                patients_perHCW = "patients[perHCW]",
                T_shift2 = "T[2]^shift",
                eps = "epsilon",
                eps_leavingH2 = "epsilon[leavingH2]")

###function that might be useful
mean_pars <- function(distrib){
  mean_pars <- NULL
  for (i in 1:n_pars){
    arg <- distrib[[2*i]]
    if (is.null(arg$mean)){
      mean_par <- .5*( arg$min + arg$max)
    }else {mean_par <- arg$mean}
    name_par <- distrib[[2*i-1]][1]
    mean_pars[[name_par]] <- mean_par
  }
  return(mean_pars)
}

###creating a set of argument desciing the parameter space explored for the sensitivity analysis
factors <- NULL
q <- NULL
q.arg <- NULL
for (i in 1:n_pars){
  factors <- c(factors, distrib[[2*i-1]][1])
  q <- c(q, distrib[[2*i-1]][2])
  q.arg <- c(q.arg,list(distrib[[2*i]]) )
}

###function calculating epidemic indicators value one simulation for each set of parameters stored in a given parameters dataframe
modelRun <- function (my.data){
  N <- nrow(my.data)
  n_indics <- length(oneRun(my.data[1,]))
  res <- NULL
  for (i in 1:N){
    res <- c(res,oneRun(my.data[i,]))
  }
  return(array(res,dim=c(n_indics,N)))
}

###exploring parameter space
res.names <- "burden reduction"  #c("pop infections reduction","pop infection peak reduction","HCW infections reduction","HCW infection peak reduction")
myLHS <- LHS(modelRun, N=150, factors, q, q.arg, res.names, nboot=100, repetitions = 200)
#myLHS <- LHS(modelRun, N=20, factors, q, q.arg, res.names, nboot=10, repetitions = 30)
###plotting PRCC
plot_beautiful_prcc <- function(myLHS, index.res = 1){
  prcc <- myLHS$prcc[[index.res]]$PRCC
  prcc_data <- data.frame(PRCC = row.names(prcc), group = c( rep('Biological', 5), rep('Sanitary context', 4), rep('Organizational', 3)),
                          value =prcc[,1] ,  abs_value = abs(prcc[,1]), sign =sign(prcc[,1]), min_c_i = prcc[,4] , max_c_i = prcc[,5] )
  prcc_data = prcc_data %>% arrange(group, abs_value)
  data = prcc_data %>% arrange(group, abs_value)
  
  # Set a number of 'empty bar' to add at the end of each group
  empty_bar <- 1
  to_add <- data.frame( matrix(NA, empty_bar*nlevels(as.factor(data$group)), ncol(data)) )
  colnames(to_add) <- colnames(data)
  to_add$group <- rep(levels(as.factor(data$group)), each=empty_bar)
  data <- rbind(data, to_add)
  data <- data %>% arrange(group)
  target <- c("Biological", "Sanitary context", "Organizational")
  data <- data %>% arrange(factor(group, levels = target))
  data$id <- seq(1, nrow(data))
  ## prepare a data frame for base lines
  #base_data <- data %>%
  #  group_by(group) %>%
  #  summarize(start=min(id), end=max(id) - empty_bar) %>%
  #  rowwise() %>%
  #  mutate(title=mean(c(start, end)))
  
  # Make the plot
  ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    
    geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity") +scale_fill_manual("Group", values=c("Biological"="#A8D08D","Sanitary context"="#ACB9CA","Organizational"  ="#F7CAAC"))+
    
    geom_errorbar( aes(x=as.factor(id), ymin=min_c_i, ymax=max_c_i), stat= "identity" ,width=0.4, colour="#A6A6A6", alpha=0.9, size=1.3) +
    ylim(-1.1,.7)+
    theme_minimal() + theme(panel.grid.major.x = element_blank(), plot.title =element_text(hjust = 0.5),axis.text.x = element_blank(), axis.title.x = element_blank())+ylab(label="PRCC")+
    geom_segment(x = 0, y = 0.5, xend = 17.5, yend = 0.5,
                 colour = "red", lty = "dotted", alpha=1, size=1 , inherit.aes = FALSE ) +
    geom_segment(x = 0, y = -0.5, xend = 17.5, yend = -0.5,
                 colour = "red", alpha=1, lty = "dotted", size=1 , inherit.aes = FALSE ) +
    
    geom_text(x = 1,
          y=(data$max_c_i[1]+0.12)*(data$value[1]>0)+(data$min_c_i[1]-0.12)*(data$value[1]<=0),
          hjust = 0.5,#abs(max_c_i)<=abs(min_c_i),
          label=parse_format()(name_pars_maths[[data$PRCC[1]]]), vjust = 0.5,
      color="black",
      fontface="bold",
      alpha=0.6,
      size=9,
      angle= 45,
      na.rm = TRUE,
      inherit.aes = FALSE ) +
    geom_text(x = 2,
              y=(data$max_c_i[2]+0.12)*(data$value[2]>0)+(data$min_c_i[2]-0.12)*(data$value[2]<=0),
              hjust = 0.5,#abs(max_c_i)<=abs(min_c_i),
              label=parse_format()(name_pars_maths[[data$PRCC[2]]]), vjust = 0.5,
              color="black",
              fontface="bold",
              alpha=0.6,
              size=9,
              angle= 45,
              na.rm = TRUE,
              inherit.aes = FALSE ) +
    geom_text(x = 3,
              y=(data$max_c_i[3]+0.12)*(data$value[3]>0)+(data$min_c_i[3]-0.12)*(data$value[3]<=0),
              hjust = 0.5,#abs(max_c_i)<=abs(min_c_i),
              label=parse_format()(name_pars_maths[[data$PRCC[3]]]), vjust = 0.5,
              color="black",
              fontface="bold",
              alpha=0.6,
              size=9,
              angle= 45,
              na.rm = TRUE,
              inherit.aes = FALSE ) +
    geom_text(x = 4,
              y=(data$max_c_i[4]+0.12)*(data$value[4]>0)+(data$min_c_i[4]-0.12)*(data$value[4]<=0),
              hjust = 0.5,#abs(max_c_i)<=abs(min_c_i),
              label=parse_format()(name_pars_maths[[data$PRCC[4]]]), vjust = 0.5,
              color="black",
              fontface="bold",
              alpha=0.6,
              size=9,
              angle= 45,
              na.rm = TRUE,
              inherit.aes = FALSE ) +
    geom_text(x = 5,
              y=(data$max_c_i[5]+0.12)*(data$value[5]>0)+(data$min_c_i[5]-0.12)*(data$value[5]<=0),
              hjust = 0.5,#abs(max_c_i)<=abs(min_c_i),
              label=parse_format()(name_pars_maths[[data$PRCC[5]]]), vjust = 0.5,
              color="black",
              fontface="bold",
              alpha=0.6,
              size=9,
              angle= 45,
              na.rm = TRUE,
              inherit.aes = FALSE ) +
    geom_text(x = 7,
              y=(data$max_c_i[7]+0.12)*(data$value[7]>0)+(data$min_c_i[7]-0.12)*(data$value[7]<=0),
              hjust = 0.5,#abs(max_c_i)<=abs(min_c_i),
              label=parse_format() (name_pars_maths[[data$PRCC[7]]]), vjust = 0.5,
              color="black",
              fontface="bold",
              alpha=0.6,
              size=9,
              angle= 45,
              na.rm = TRUE,
              inherit.aes = FALSE ) +
    geom_text(x = 8,
              y=(data$max_c_i[8]+0.12)*(data$value[8]>0)+(data$min_c_i[8]-0.12)*(data$value[8]<=0),
              hjust = 0.5,#abs(max_c_i)<=abs(min_c_i),
              label=parse_format() (name_pars_maths[[data$PRCC[8]]]), vjust = 0.5,
              color="black",
              fontface="bold",
              alpha=0.6,
              size=9,
              angle= 45,
              na.rm = TRUE,
              inherit.aes = FALSE ) +
    geom_text(x = 9,
              y=(data$max_c_i[9]+0.12)*(data$value[9]>0)+(data$min_c_i[9]-0.12)*(data$value[9]<=0),
              hjust = 0.5,#abs(max_c_i)<=abs(min_c_i),
              label=parse_format() (name_pars_maths[[data$PRCC[9]]]), vjust = 0.5,
              color="black",
              fontface="bold",
              alpha=0.6,
              size=9,
              angle= 45,
              na.rm = TRUE,
              inherit.aes = FALSE ) +
    geom_text(x = 10,
              y=(data$max_c_i[10]+0.12)*(data$value[10]>0)+(data$min_c_i[10]-0.12)*(data$value[10]<=0),
              hjust = 0.5,#abs(max_c_i)<=abs(min_c_i),
              label=parse_format() (name_pars_maths[[data$PRCC[10]]]), vjust = 0.5,
              color="black",
              fontface="bold",
              alpha=0.6,
              size=9,
              angle= 45,
              na.rm = TRUE,
              inherit.aes = FALSE ) +

    geom_text(x = 12,
              y=(data$max_c_i[12]+0.12)*(data$value[12]>0)+(data$min_c_i[12]-0.12)*(data$value[12]<=0),
              hjust = 0.5,#abs(max_c_i)<=abs(min_c_i),
              label=parse_format() (name_pars_maths[[data$PRCC[12]]]), vjust = 0.5,
              color="black",
              fontface="bold",
              alpha=0.6,
              size=9,
              angle= 45,
              na.rm = TRUE,
              inherit.aes = FALSE ) +
    geom_text(x = 13,
              y=(data$max_c_i[13]+0.12)*(data$value[13]>0)+(data$min_c_i[13]-0.12)*(data$value[13]<=0),
              hjust = 0.5,#abs(max_c_i)<=abs(min_c_i),
              label=parse_format() (name_pars_maths[[data$PRCC[13]]]), vjust = 0.5,
              color="black",
              fontface="bold",
              alpha=0.6,
              size=9,
              angle= 45,
              na.rm = TRUE,
              inherit.aes = FALSE ) +
    geom_text(x = 14,
              y=(data$max_c_i[14]+0.12)*(data$value[14]>0)+(data$min_c_i[14]-0.12)*(data$value[14]<=0),
              hjust = 0.5,#abs(max_c_i)<=abs(min_c_i),
              label=parse_format() (name_pars_maths[[data$PRCC[14]]]), vjust = 0.5,
              color="black",
              fontface="bold",
              alpha=0.6,
              size=9,
              angle= 45,
              na.rm = TRUE,
              inherit.aes = FALSE ) + theme(text = element_text(size = 20))
} # index.res is the indicator index of which the user choses to display the PRCC graph
plot_beautiful_prcc(myLHS) #remplace plotprcc(myLHS,index.res = 1)

###graph that allows to evaluate the convergence of the calculated PRCCs
plotcv(myLHS)

###additional lines to evaluate the convergence
#####newLHS <- LHS(modelRun, N=200, factors, q, q.arg, res.names, nboot=100, repetitions = 200)
#####mySbma <- sbma(myLHS,newLHS)



##HEATMAPS
###indicators chosen for the built heatmap
indicators_name <- c("infections_reduc","infection_peak_reduc","HCW_infections_reduc","HCW_infection_peak_reduc","nosocomial_burden_reduc")

default_pars <- paper_pars
build_heatmap_dataset <- function(pars,par_x, x_min, x_max, ncol, par_y, y_min, y_max, nrow, n_sim){
  x_range <- seq(x_min, x_max, length.out = ncol)
  y_range <- seq(y_min, y_max, length.out = nrow)
  reference_value <- NULL
  pars_built <- pars
  dataset <- NULL
  for (ind in indicators_name) {dataset[[ind]] <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c(par_x, par_y,ind))}
  for (l in 1:nrow){
    for (k in 1:ncol){
      x_par_value <- x_range[k]
      y_par_value <- y_range[l]
      pars_built[[par_x]] <- x_par_value
      pars_built[[par_y]] <- y_par_value
      data_value <- severalRun(n_sim,pars_built)
      
      for (ind in indicators_name){
        dataset[[ind]][nrow(dataset[[ind]]) + 1,] <- c(x_par_value,y_par_value,data_value[[ind]])
      }
    }
  }
  dataset[["min"]] <- NULL
  dataset[["max"]] <- NULL
  for (ind in indicators_name) {
    dataset[["min"]][[ind]] <- min(dataset[[ind]][[3]],na.rm = TRUE)
    dataset[["max"]][[ind]] <- max(dataset[[ind]][[3]],na.rm = TRUE)
  }
  dataset[["def"]] <- list( par_x = par_x, par_y = par_y)
  return(dataset)
}

###heatmpa resolution = nb_dots*nb_dots
nb_dots = 15
###number of independent realisations for each point
n_sim=100

###data creation !!! very long compilation
R_C_x_eps <- build_heatmap_dataset(paper_pars,"R_C",0.85,1.3,nb_dots,"eps",0.05,0.3,nb_dots,n_sim)


###indicateurs we chose to plot
####indicators_name_plot <- c("HCW_infection_peak_reduc")

####generate associated title names
indicators_title <- function(indicators_name_plot){
  indicators_title <- NULL
  for (ind in indicators_name_plot){
    if (ind=="infections_reduc"){
      indicators_title[[ind]] <- "Total infections\nreduction (%)"
    }else if (ind=="infection_peak_reduc"){
      indicators_title[[ind]] <- "Infection peak\nreduction (%)"
    }else if (ind=="HCW_infections_reduc"){
      indicators_title[[ind]] <-"HCWs infection\nreduction (%)"
    }else if (ind=="HCW_infection_peak_reduc"){
      indicators_title[[ind]] <-"HCWs infection\npeak reduction (%)"
    }else {
      indicators_title[[ind]] <-"Nosocomial burden\nreduction (%)"
    }
  }
  return(indicators_title)
}

###adjusting dimensions when arranging several graphs
opt_dims <- function(N){
  ncol <- floor(sqrt(N))
  nrow <- ncol
  if (N>ncol*nrow){ncol <- ncol + 1}
  if (N>ncol*nrow){nrow <- nrow + 1}
  return(c(ncol,nrow))
}

indicators_name_plot <- c("infections_reduc","infection_peak_reduc",
                          "HCW_infections_reduc","HCW_infection_peak_reduc")#,"nosocomial_burden_reduc")


plot_heatmap<- function(ind){#dataset){
  dataset <- R_C_x_eps
  min <- 100
  max <- 0
  
  if (dataset$min[[ind]]<min){min <- dataset$min[[ind]]}
  if (dataset$max[[ind]]>max){max <- dataset$max[[ind]]}
  
  lims = c(0,80)#min,max)
  par_x <- dataset$def$par_x
  par_y <- dataset$def$par_y
  
  ggplot(dataset[[ind]]) +
    geom_tile(aes_string(x = par_x, y = par_y, fill = ind)) +
    scale_fill_gradient2(low = "white",  mid = "darkolivegreen1",  high = "green", midpoint = (min+max)/2, lim = lims) +
    #scale_fill_gradient2(low = "white",  mid = "#56B4E9",  high = "#0072B2", midpoint = (min+max)/2, lim = lims) +
    
    labs( fill = indicators_title(indicators_name_plot)[[ind]] ) +
    geom_contour(aes_string(x = par_x, y = par_y, z = ind), colour = "orange",
                 breaks = c(-20,-10,0,10,15,20,30,40,50),show.legend = FALSE) +
    geom_text_contour(aes_string(x = par_x, y = par_y, z = ind),
                      breaks = c(-20,-10,0,10,15,20,30,40,50) ,skip = 0, size = 6) +
    
    xlab(parse_format()("R[C]")) + ylab(parse_format()("epsilon"))+
    
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title = element_text(size = 20), legend.title  = element_text(size = 15),
          axis.text = element_text(size = 15), legend.text  = element_text(size = 12))
  #arrange_pars <- NULL
  
  
}


###display heatmap
plot_heatmap("infections_reduc")#R_C_x_eps)
plot_heatmap("infection_peak_reduc")
plot_heatmap("HCW_infections_reduc")
plot_heatmap("HCW_infection_peak_reduc")

