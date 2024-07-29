
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)

df_all = readRDS(here::here("Results", "df_all.rds"))


###function for plotting epidemic curves
plot_curves <- function(df, columns, curve_colors, labels, ylab = "Proportion (%)",
                        compare_peak = F, compare_attack = F){
  
  n <- length(df[[1]]$time_ref)
  df_for_plot <- NULL
  for (col in columns){
    df_col <- data.frame(date = df[[1]]$time_ref, nb_5p = df[[1]][[col]], nb = df[[2]][[col]], nb_95p = df[[3]][[col]], comp = rep(col, n))
    df_for_plot <- rbind(df_for_plot,df_col)
  }
  df_for_plot <- df_for_plot %>% arrange(factor(comp, levels = columns))
  
  if(compare_peak){
    df_for_plot %>%
      group_by(comp) %>%
      filter(nb==max(nb)) %>%
      summarise(date = median(date), nb_5p = median(nb_5p), nb = median(nb), nb_95p = median(nb_95p)) %>%
      print
  }
  
  if(compare_attack){
    df_for_plot %>%
      group_by(comp) %>%
      filter(date==max(date)) %>%
      group_by(date) %>%
      summarise(nb=max(nb)-min(nb), nb_5p=max(nb_5p)-min(nb_5p), nb_95p=max(nb_95p)-min(nb_95p)) %>%
      print
  }
  
  ggplot(data=df_for_plot, aes(x=date, y=nb, ymin=nb_5p, ymax=nb_95p, fill=comp)) +
    geom_line(aes(colour=comp), linewidth=1) +
    geom_ribbon(alpha=0.5) +
    labs(x = "Time (days)", y = ylab) + 
    guides(linetype="none") +
    theme_bw() +
    scale_fill_manual(values = curve_colors, name = "", limits = columns, labels = labels)+
    scale_colour_manual(values = curve_colors, name = "", limits = columns, labels = labels)+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  
}

col_pal = c("grey50", "royalblue3")

####usual epidemic curves
prop_infection_non_hcw <- c(paste(c("Prop_inf_p"),"ref",sep="_"),paste("Prop_inf_p","quarantine",sep="_"))
pa = plot_curves(df_all, prop_infection_non_hcw, col_pal,
                 c("Reference strategy","Quarantine strategy"), "Infectious non-HCWs (%)", compare_peak = T)

prop_infection_hcw <- c(paste(c("Prop_inf_hcw"),"ref",sep="_"),paste("Prop_inf_hcw","quarantine",sep="_"))
pb = plot_curves(df_all, prop_infection_hcw, col_pal,
                 c("Reference strategy","Quarantine strategy"), "Infectious HCWs (%)", compare_peak = T)

burden_reduction_prop_non_hcw <- c(paste(c("cumul_prop_inf_p"),"ref",sep="_"),paste("cumul_prop_inf_p","quarantine",sep="_"))
pc = plot_curves(df_all, burden_reduction_prop_non_hcw, col_pal,
                 c("Reference strategy","Quarantine strategy"), "Attack rate for non-HCWs (%)", compare_attack = T)

burden_reduction_prop_hcw <- c(paste(c("cumul_prop_inf_hcw"),"ref",sep="_"),paste("cumul_prop_inf_hcw","quarantine",sep="_"))
pd = plot_curves(df_all, burden_reduction_prop_hcw, col_pal,
                 c("Reference strategy","Quarantine strategy"), "Attack rate for HCWs (%)", compare_attack = T)


# if legend_position = "bottom", legend is 3rd component of guide_box
legend_for_plot = get_plot_component(pa + theme(legend.position = "bottom"), "guide-box", return_all = TRUE)[[3]]

plot_grid(plot_grid(pa+theme(legend.position = "none"), pb+theme(legend.position = "none"),
                    pc+theme(legend.position = "none"), pd+theme(legend.position = "none"),
                    labels = c("a)", "b)", "c)", "d)"), hjust = 0),
          legend_for_plot,
          rel_heights = c(1,0.05), ncol = 1)

ggsave(here::here("Figures", "fig2.png"), height = 7)



# virus_exposition <- c(paste(c("S_p_C_prop"),"ref",sep="_"),paste("S_p_C_prop","quarantine",sep="_"))
# plot_curves(df_all, virus_exposition, c("blue","red"),
#             c("Reference strategy","Quarantine strategy"),"Proportion of susceptible non-HCWs among community (%)")


col_pal2 <- brewer.pal(4, "Dark2")

FI_on_p_C_evo_ref <- paste(c("FI_on_p_C_prop_from_p","FI_on_p_C_prop_from_hcw_C","FI_on_p_C_prop_from_hcw_L"),"ref",sep="_") #,"FI_on_p_C_prop_from_hcw_AL"

n <- length(df_all[[1]]$time_ref)
df_for_plot <- NULL
for (col in FI_on_p_C_evo_ref){
  df_col <- data.frame(date = df_all[[1]]$time_ref, nb_5p = df_all[[1]][[col]], nb = df_all[[2]][[col]], nb_95p = df_all[[3]][[col]],
                       comp = rep(col, n), strategy = "Reference strategy")
  df_for_plot <- rbind(df_for_plot,df_col)
}
FI_on_p_C_evo_quar <- paste(c("FI_on_p_C_prop_from_p","FI_on_p_C_prop_from_hcw_C","FI_on_p_C_prop_from_hcw_L","FI_on_p_C_prop_from_hcw_AL"),
                            "quarantine",sep="_") #,"FI_on_p_C_prop_from_hcw_AL"
for (col in FI_on_p_C_evo_quar){
  df_col <- data.frame(date = df_all[[1]]$time_ref, nb_5p = df_all[[1]][[col]], nb = df_all[[2]][[col]], nb_95p = df_all[[3]][[col]],
                       comp = rep(col, n), strategy = "Quarantine hospital strategy")
  df_for_plot <- rbind(df_for_plot,df_col)
}

df_for_plot = df_for_plot %>%
  mutate(comp = replace(comp, grepl("FI_on_p_C_prop_from_p", comp), "From non-HCWs")) %>%
  mutate(comp = replace(comp, grepl("FI_on_p_C_prop_from_hcw_C", comp), "From resting HCWs")) %>%
  mutate(comp = replace(comp, grepl("FI_on_p_C_prop_from_hcw_L", comp), "From symptomatic\nisolated HCWs")) %>%
  mutate(comp = replace(comp, grepl("FI_on_p_C_prop_from_hcw_AL", comp), "From HCWs isolated after working\nin the quarantine hospital"))

df_for_plot <- df_for_plot %>%
  mutate(comp = factor(comp,
                        levels = c("From non-HCWs", "From resting HCWs",  "From symptomatic\nisolated HCWs", 
                                   "From HCWs isolated after working\nin the quarantine hospital")))
  

ggplot(data=df_for_plot, aes(x=date, y=nb, ymin=nb_5p, ymax=nb_95p, fill=comp)) +
  geom_line(aes(colour=comp), linewidth=1) +
  geom_ribbon(alpha=0.5) +
  
  geom_hline(aes(yintercept = df_for_plot %>% filter(strategy=="Reference strategy" & comp == "From non-HCWs") %>%
                   group_by(comp) %>% summarise(nb = median(nb)) %>% select(nb) %>% pull, colour=comp),
             data = data.frame(strategy = "Reference strategy", comp = "From non-HCWs"), linetype = "dashed", linewidth = 1) +
  
  geom_hline(aes(yintercept = df_for_plot %>% filter(strategy=="Reference strategy" & comp == "From resting HCWs") %>%
                   group_by(comp) %>% summarise(nb = median(nb)) %>% select(nb) %>% pull, colour=comp),
             data = data.frame(strategy = "Reference strategy", comp = "From resting HCWs"), linetype = "dashed", linewidth = 1) +
  
  geom_hline(aes(yintercept = df_for_plot %>% filter(strategy=="Reference strategy" & comp == "From symptomatic\nisolated HCWs") %>%
                   group_by(comp) %>% summarise(nb = median(nb)) %>% select(nb) %>% pull, colour=comp),
             data = data.frame(strategy = "Reference strategy", comp = "From symptomatic\nisolated HCWs"), linetype = "dashed", linewidth = 1) +

  
  geom_hline(aes(yintercept = df_for_plot %>% filter(strategy=="Quarantine hospital strategy" & comp == "From non-HCWs") %>%
                   group_by(comp) %>% summarise(nb = median(nb)) %>% select(nb) %>% pull, colour=comp),
             data = data.frame(strategy = "Quarantine hospital strategy", comp = "From non-HCWs"), linetype = "dashed", linewidth = 1) +
  
  geom_hline(aes(yintercept = df_for_plot %>% filter(strategy=="Quarantine hospital strategy" & comp == "From resting HCWs") %>%
                   group_by(comp) %>% summarise(nb = median(nb)) %>% select(nb) %>% pull, colour=comp),
             data = data.frame(strategy = "Quarantine hospital strategy", comp = "From resting HCWs"), linetype = "dashed", linewidth = 1) +
  
  geom_hline(aes(yintercept = df_for_plot %>% filter(strategy=="Quarantine hospital strategy" & comp == "From symptomatic\nisolated HCWs") %>%
                   group_by(comp) %>% summarise(nb = median(nb)) %>% select(nb) %>% pull, colour=comp),
             data = data.frame(strategy = "Quarantine hospital strategy", comp = "From symptomatic\nisolated HCWs"), linetype = "dashed", linewidth = 1) +

  geom_hline(aes(yintercept = df_for_plot %>% filter(strategy=="Quarantine hospital strategy" & comp == "From HCWs isolated after working\nin the quarantine hospital") %>%
                   group_by(comp) %>% summarise(nb = median(nb)) %>% select(nb) %>% pull, colour=comp),
             data = data.frame(strategy = "Quarantine hospital strategy", comp = "From HCWs isolated after working\nin the quarantine hospital"), linetype = "dashed", linewidth = 1) +
  
    
  facet_grid(.~factor(strategy, levels = c("Reference strategy", "Quarantine hospital strategy"))) +
  labs(x = "Time (days)", y = "Acquisition route frequency", colour = "", fill = "") + 
  guides(linetype="none") +
  theme_bw() +
  scale_fill_discrete(type = col_pal2) +
  scale_colour_discrete(type = col_pal2) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11),
        legend.position = "bottom")

ggsave(here::here("Figures", "fig3.png"))


# supplementary stuff

Nb_of_patients_ref <- paste(c("mes_N_p_1", "I_p_1"),"ref",sep="_") 
plot_curves(df_all, Nb_of_patients_ref, c("black","red"),
            c( "In usual hospital","Of which symptomatic patients"),
            "Number of patients")

ggsave(here::here("Figures", "suppfigx1.png"))

Nb_of_patients_quarantine <- paste(c("mes_N_p_1","mes_N_p_2"),"quarantine",sep="_") 
plot_curves(df_all, Nb_of_patients_quarantine, c("black", "red"),
            c("In usual hospital", "In quarantine hospital"),
            "Number of patients")

ggsave(here::here("Figures", "suppfigx2.png"))

Nb_of_hcw_hospitals <- c("mes_N_hcw_1_ref","mes_N_hcw_1_quarantine","mes_N_hcw_2_quarantine")
plot_curves(df_all, Nb_of_hcw_hospitals, c("black", "grey","red"),
            c( "In usual hospital\nno quarantine", "In usual hospital\nwith quarantine", "In quarantine hospital"),
            "Number of HCWs")

ggsave(here::here("Figures", "suppfigx3.png"))

Nb_of_hcw_C <- c("mes_N_hcw_L_ref","mes_N_hcw_L_quarantine")
plot_curves(df_all, Nb_of_hcw_C, c( "black", "blue"),
            c( "Reference strategy","Quarantine strategy"),
            "Number of HCWs in the community" )

ggsave(here::here("Figures", "suppfigx4.png"))

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

df_risk <- data.frame(date = 1:51, nb_5p = df_all[[1]]$risk_hcw_2_shift_quarantine[7*(1:51)+3]
                      , nb = df_all[[2]]$risk_hcw_2_shift_quarantine[7*(1:51)+3]
                      , nb_95p = df_all[[3]]$risk_hcw_2_shift_quarantine[7*(1:51)+3])

# add number of contaminations to show that large bands at the end are because of small number left susceptible

## CHECK: was plotted up to 42 weeks, increase to 52 for total period coverage ????

ggplot(data=df_risk, aes(x=date, y=nb, ymin=nb_5p, ymax=nb_95p)) +
  geom_point() +
  geom_ribbon(alpha=0.5) +
  labs(x = "Time (weeks)", y = "7-days shift risk (%)") +
  guides(linetype="none") +
  theme_bw() +
  theme(legend.position="top",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

ggsave(here::here("Figures", "suppfigx5.png"))

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

