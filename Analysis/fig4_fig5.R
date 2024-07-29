

library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(cowplot)

df_all = readRDS(here::here("Results", "df_all.rds"))

source(here::here("Analysis", "model.R"))

##HEATMAPS
###indicators chosen for the built heatmap
indicators_name <- c("infections_reduc","infection_peak_reduc","HCW_infections_reduc","HCW_infection_peak_reduc","nosocomial_burden_reduc")

paper_pars <- list(T_E = 4.6, T_I = 9.8, phi = 0.2, psi = 0.07,
                   kappa = 0.32 ,R_C = 1.05, R_1_ref = 1.15,
                   R_1_quarantine = 1.1 , R_2 = 1.2, 
                   patients_perHCW = 3,
                   T_shift2 = 7 , eps = 0.15,eps_leavingH2 = 0.4)

default_pars <- paper_pars

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

###call the latter function to display indicators for the paper
severalRun_withPB(100, paper_pars)

###heatmap resolution = nb_dots*nb_dots
nb_dots = 15
###number of independent realisations for each point
n_sim=100

###data creation !!! very long compilation
suppressWarnings({
  R_C_x_eps <- build_heatmap_dataset(paper_pars,"R_C",0.85,1.3,nb_dots,"eps",0.05,0.3,nb_dots,n_sim)
})


col_low = "white"
col_mid = "slateblue1"
col_high = "slateblue4"
par_x <- R_C_x_eps$def$par_x
par_y <- R_C_x_eps$def$par_y

pa_data = left_join(R_C_x_eps[["infections_reduc"]],
                    R_C_x_eps[["HCW_infections_reduc"]]) %>%
  dplyr::rename("All infections" = "infections_reduc") %>%
  dplyr::rename("HCWs infections" = "HCW_infections_reduc") %>%
  melt(id = c("R_C", "eps"))

lims = c(0,80)

pa = ggplot(pa_data) +
  geom_tile(aes(x = R_C, y = eps, fill = value)) +
  scale_fill_gradient2(low = col_low,  mid = col_mid,  high = col_high, midpoint = (lims[1]+lims[2])/2, lim = lims) +
  facet_grid(.~variable) +
  geom_contour(aes(x = R_C, y = eps, z = value), colour = "grey40",
               breaks = c(-20,-10,0,10,20,30,40,50), show.legend = FALSE, linewidth = 1) +
  geom_text_contour(aes(x = R_C, y = eps, z = value),
                    breaks = c(-20,-10,0,10,20,30,40,50) ,skip = 0, size = 6, fontface = "bold") +
  xlab(parse_format()("R[C]")) +
  ylab(parse_format()("epsilon")) +
  labs(fill = "Total infections\nreduction (%)") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        legend.title  = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text  = element_text(size = 12),
        strip.text = element_text(size = 11))

col_low = "white"
col_mid = "tan1"
col_high = "tan4"
par_x <- R_C_x_eps$def$par_x
par_y <- R_C_x_eps$def$par_y

pb_data = left_join(R_C_x_eps[["infection_peak_reduc"]],
                    R_C_x_eps[["HCW_infection_peak_reduc"]]) %>%
  dplyr::rename("All infections" = "infection_peak_reduc") %>%
  dplyr::rename("HCWs infections" = "HCW_infection_peak_reduc") %>%
  melt(id = c("R_C", "eps"))

lims = c(0,80)

pb = ggplot(pb_data) +
  geom_tile(aes(x = R_C, y = eps, fill = value)) +
  scale_fill_gradient2(low = col_low,  mid = col_mid,  high = col_high, midpoint = (lims[1]+lims[2])/2, lim = lims) +
  facet_grid(.~variable) +
  geom_contour(aes(x = R_C, y = eps, z = value), colour = "grey40",
               breaks = c(-20,-10,0,10,20,30,40,50), show.legend = FALSE, linewidth = 1) +
  geom_text_contour(aes(x = R_C, y = eps, z = value),
                    breaks = c(-20,-10,0,10,20,30,40,50) ,skip = 0, size = 6, fontface = "bold") +
  xlab(parse_format()("R[C]")) +
  ylab(parse_format()("epsilon")) +
  labs(fill = "Infection peak\nreduction (%)") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        legend.title  = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text  = element_text(size = 12),
        strip.text = element_text(size = 11))

plot_grid(pa, pb, ncol = 1, align = "v", labels = c("a)", "b)"), hjust = 0)

ggsave(here::here("Figures", "fig4.png"))


col_low = "white"
col_mid = "violetred2"
col_high = "violetred4"
par_x <- R_C_x_eps$def$par_x
par_y <- R_C_x_eps$def$par_y

plot_data = R_C_x_eps[["nosocomial_burden_reduc"]]

min <- R_C_x_eps$min[["nosocomial_burden_reduc"]]
max <- R_C_x_eps$max[["nosocomial_burden_reduc"]]
lims = c(min,max)

ggplot(plot_data) +
  geom_tile(aes(x = R_C, y = eps, fill = nosocomial_burden_reduc)) +
  scale_fill_gradient2(low = col_low,  mid = col_mid,  high = col_high, midpoint = (lims[1]+lims[2])/2, lim = lims) +
  geom_contour(aes(x = R_C, y = eps, z = nosocomial_burden_reduc), colour = "grey40",
               breaks = c(86, 90, 94), show.legend = FALSE, linewidth = 1) +
  geom_text_contour(aes(x = R_C, y = eps, z = nosocomial_burden_reduc),
                    breaks = c(86, 90, 94) ,skip = 0, size = 6, fontface = "bold") +
  xlab(parse_format()("R[C]")) +
  ylab(parse_format()("epsilon")) +
  labs(fill = "Nosocomial burden\nreduction (%)") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        legend.title  = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text  = element_text(size = 12),
        strip.text = element_text(size = 11))

ggsave(here::here("Figures", "fig5.png"))
