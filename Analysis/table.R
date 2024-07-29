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
