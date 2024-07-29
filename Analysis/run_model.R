
source(here::here("Analysis", "model.R"))

###build the dataset containing simulations according to the baseline values, chosen ending date and number of simulations
df_all <- model("ref",paper_pars, date_lim,nb_sim)
df_quarantine <- model("quarantine",paper_pars, date_lim,nb_sim)
for (i in 1:3){
  df_all[[i]] <- cbind(df_all[[i]],df_quarantine[[i]])
}


saveRDS(df_all, here::here("Results", "df_all.rds"))
