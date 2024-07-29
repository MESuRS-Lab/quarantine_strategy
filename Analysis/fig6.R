
library(ggplot2)
library(dplyr)

source(here::here("Analysis", "model.R"))

##Sensitivity analysis
###parameter space explored
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
                R_C_arg = list(a=0, b=Inf, mean=1.05, sd=0.105),
                R_1_ref = c("R_1_ref","qtruncnorm"),
                R_1_ref_arg = list( a=0, b=Inf, mean=1.15, sd=0.115),
                R_1_quarantine = c("R_1_quarantine","qtruncnorm"),
                R_1_quarantine_arg = list(a=0, b=Inf, mean=1.1, sd=0.11),
                R_2 = c("R_2","qtruncnorm"),
                R_2_arg = list(a=0, b=Inf, mean=1.2, sd=0.12),
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
                        psi = "psi",
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

suppressWarnings({myLHS <- LHS(modelRun, N=300, factors, q, q.arg, res.names,
                               nboot=100, repetitions = 500)})

###plotting PRCC

myLHS
index.res = 1
# index.res is the indicator index of which the user chooses to display the PRCC graph

prcc <- myLHS$prcc[[index.res]]$PRCC
prcc_data <- data.frame(PRCC = row.names(prcc), group = c( rep('Biological', 5), rep('Sanitary context', 4), rep('Organizational', 3)),
                        value =prcc[,1] ,  abs_value = abs(prcc[,1]), sign =sign(prcc[,1]), min_c_i = prcc[,4] , max_c_i = prcc[,5] )
prcc_data = prcc_data %>% arrange(group, abs_value)
data = prcc_data %>% arrange(group, abs_value)

data <- data %>%
  mutate(group = factor(group, levels = c("Biological", "Sanitary context", "Organizational"))) %>%
  arrange(group) %>%
  mutate(PRCC = factor(PRCC, levels = c("T_E", "T_I", "phi", "kappa", "psi",
                                        "R_2", "R_1_ref", "R_1_quarantine", "R_C",
                                        "T_shift2", "eps_leavingH2", "eps")))
data$id <- seq(1, nrow(data))

# Make the plot

ggplot(data) +
  geom_bar(aes(PRCC, value, fill = group), stat= "identity") +
  geom_errorbar(aes(x=PRCC, ymin=min_c_i, ymax=max_c_i),
                stat= "identity" ,width=0.4, colour="#A6A6A6", alpha=0.9, linewidth=1.3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = -0.5, linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual("Group", values=c("Biological"="#A8D08D","Sanitary context"="#ACB9CA","Organizational"  ="#F7CAAC"))+
  scale_y_continuous(breaks = seq(-1,1,0.25)) +
  scale_x_discrete(labels = c(bquote(T[E]),
                              bquote(T[I]),
                              bquote(phi),
                              bquote(kappa),
                              bquote(psi),
                              bquote(beta[2]),
                              bquote(beta[1]^ref),
                              bquote(beta[1]^q),
                              bquote(beta[C]),
                              bquote(T[2]^shift),
                              bquote(epsilon[H2]),
                              bquote(epsilon))) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  labs(x = "\n", y = "Correlation coefficient") +
  guides(fill = "none") +
  coord_cartesian(clip = "off", ylim = c(-1,1)) +
  annotate("text", x = 3, y = -1.44, label = "Biological parameters") +
  annotate("text", x = 7.5, y = -1.44, label = "Sanitary context\nparameters") +
  annotate("text", x = 11, y = -1.44, label = "Organizational\nparameters") +
  annotate("segment", x = 1, xend = 5, y = -1.33, yend = -1.33) +
  annotate("segment", x = 6, xend = 9, y = -1.33, yend = -1.33) +
  annotate("segment", x = 10, xend = 12, y = -1.33, yend = -1.33)

ggsave(here::here("Figures", "fig6.png"))

###additional lines to evaluate the convergence
###graph that allows to evaluate the convergence of the calculated PRCCs
plotcv(myLHS)

#####newLHS <- LHS(modelRun, N=200, factors, q, q.arg, res.names, nboot=100, repetitions = 200)
#####mySbma <- sbma(myLHS,newLHS)

