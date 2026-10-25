# R script for STB Project ####
# CO2 emission ####

require(tidybayes)
require(brms)
require(rstan)
require(bayesplot)

library(vegan)
library(tidyverse)
library(ggpubr)

library(ggplot2)
library(patchwork)
library(cowplot)




#*********************####
## Load data ####
data_co2 <- read_csv("respiration/data_respiration.csv") |>
  mutate(time_cal = measurement_point * 0.5,
         date = as.factor(date))

metadata_soil <- read_csv("metadata_soil.csv")


## Bayesian model ####

model_respiration <- bf(CO2 ~ C * (1 - moisture) * (1 - exp(-K * (1 - moisture) * (time_cal + D))),
                        C ~ 0 + exp_id:date, 
                        K ~ 0 + exp_id:date, 
                        D ~ 1 + date + (1 | exp_id),
                        nl = TRUE, family = gaussian()
)

p <- c(
  set_prior("normal(8000,4000)", class = "b", nlpar = "C", lb = 0),
  set_prior("normal(0.075,0.015)", class = "b", nlpar = "K", lb = 0, ub = 0.2),
  set_prior("normal(0,1)", class = "b", nlpar = "D", lb = 0, ub = 6),
  set_prior("normal(0,500)", class = "sigma", lb = 0, ub = 600)
)

fit_respiration <- brm(model_respiration, data_co2,
                       prior = p,
                       cores = 10, chains = 10,
                       warmup = 2000, iter = 3000,
                       control = list(adapt_delta = 0.99, max_treedepth = 12),
                       drop_unused_levels = TRUE,
                       init = 0,
                       seed = 123,
                       file = "respiration/models_posteriors/fit_respiration_final"
)



### graphical posterior predictive checks ####
#print(fit_respiration, probs = c(0.025, 0.5, 0.975))

#get_variables(fit_respiration)

# Supplementary

p1 <- pp_check(fit_respiration, ndraws = 100) +
  labs(title = "Posterior Predictive Checks: Respiration Model",
       x = "Value",
       y = "Density")+
  theme_bw(base_size = 15)

p1

p2 <- pp_check(fit_respiration, type = "scatter_avg") +
  labs(title = "Observed vs Predicted Average",
       x = "Average Predicted Value",
       y = "Observed Value") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  theme_bw(base_size = 15)

p2

p3 <- pp_check(fit_respiration, type = "stat", stat = "mean") +
  labs(title = "Posterior Predictive Check: Mean") +
  theme_bw(base_size = 15)

p3

p4 <- pp_check(fit_respiration, type = "stat", stat = "sd") +
  labs(title = "Posterior Predictive Check: SD") +
  theme_bw(base_size = 15)

p4

p <- plot_grid(
  p1,p2,p3,p4,
  ncol=2,
  align="hv",
  axis="tblr",
  rel_widths=c(1, 1, 1, 1)
)

p

ggsave("Figure_SI_posterior_predictive_checks_fit_respiration.jpeg", dpi=600, width = 34, height = 25, units = "cm")




### posterior draws ####

fit_respiration <- read_rds("respiration/models_posteriors/fit_respiration_final.rds")

## posterior draws D
fix_D <- as_draws_df(fit_respiration) %>%
  select(.draw, starts_with("b_D"))

ran_D <-  fit_respiration %>%
  spread_draws(r_exp_id__D[exp_id, term]) %>%
  filter(term == "Intercept") %>%
  select(.draw, exp_id, r_D = r_exp_id__D)

draws_D <- expand_grid(
  .draw = unique(fix_D$.draw),
  exp_id = unique(ran_D$exp_id),
  date = c(0,3,6,9,15,21,27)) %>%
  left_join(fix_D, by=".draw") %>%
  left_join(ran_D, by=c(".draw","exp_id")) %>%
  mutate(D = b_D_Intercept + r_D +
           case_when(
             date==0  ~ 0,
             date==3  ~ b_D_date3,
             date==6  ~ b_D_date6,
             date==9  ~ b_D_date9,
             date==15 ~ b_D_date15,
             date==21 ~ b_D_date21,
             date==27 ~ b_D_date27)) %>%
  dplyr::select(.draw,exp_id,date,D)

# posterior draws respiration
draws_respiration <- fit_respiration %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  ungroup() %>%
  mutate(exp_id = str_extract(.variable, "(?<=exp_id)[^:]+")) %>%
  mutate(date = str_extract(.variable, "(?<=date)\\w+")) %>%
  mutate(variable = str_extract(.variable, "(?<=b_)[^_]+")) %>%
  select(-.variable) %>%
  pivot_wider(names_from = "variable", values_from = .value) %>%
  mutate(date = as.numeric(date)) %>%
  select(-D) %>%
  drop_na(exp_id)




## Emission rate ####
df_respiration <- draws_respiration %>%
  left_join(metadata_soil %>%
              select(exp_id, group, season, moisture) %>% distinct(), by = c("exp_id"))%>%
  left_join(draws_D, by=c(".draw","exp_id","date")) %>%
  mutate(initial_rate = C*K*(1-moisture)^2) %>%
  rename(Group = group, Season = season, Day = date, Exp_id = exp_id) %>%
  filter(!(Group %in% c("SW2","SW3") & Day %in% c("9","27"))) %>%
  mutate(Soil_id = gsub("(S)[CW]", "\\1", Group))%>%
  mutate(Site_type = case_when(Soil_id == "S1" ~ "Agriculture",
                               Soil_id == "S2" ~ "Grassland",
                               Soil_id == "S3" ~ "Forest"))%>%
  mutate(Soil = case_when(
    Group == "SC1" ~ "Cold - Agriculture",
    Group == "SC2" ~ "Cold - Grassland",
    Group == "SC3" ~ "Cold - Forest",
    Group == "SW1" ~ "Warm - Agriculture",
    Group == "SW2" ~ "Warm - Grassland",
    Group == "SW3" ~ "Warm - Forest")) 

saveRDS(df_respiration, "respiration/df_respiration.rds")




## Respiration summary ####
df_respiration <- readRDS("respiration/df_respiration.rds")

summary_respiration <- df_respiration %>%
  group_by(Exp_id, Day, Group) %>%
  summarise(
    ci_low_C = quantile(C, 0.025),
    ci_high_C = quantile(C, 0.975),
    C = median(C),
    ci_low_K = quantile(K, 0.025),
    ci_high_K = quantile(K, 0.975),
    K = median(K),
    ci_low_D = quantile(D, 0.025),
    ci_high_D = quantile(D, 0.975),
    D = median(D),
    ci_low_R0 = quantile(initial_rate, 0.025),
    ci_high_R0 = quantile(initial_rate, 0.975),
    initial_rate = median(initial_rate),
    .groups = "drop"
  ) %>%
  filter(!(Group %in% c("SW2","SW3") & Day %in% c(9,27)))

saveRDS(summary_respiration, "respiration/summary_respiration.rds")


#*********************####
# Correlation with redox activity ####

## Load data ####

summary_respiration <- readRDS("respiration/summary_respiration.rds")

df_redox <- readRDS("logos/df_redox.rds")

summary_redox <- df_redox %>%
  group_by(Exp_id, Day, Group) %>%
  summarise(
    median_fluorescence_se = sd(median_fluorescence)/sqrt(n()),
    median_fluorescence = mean(median_fluorescence),
    rsg_density_se = sd(density)/sqrt(n()),
    rsg_density = mean(density),
    syb_density_se = sd(syb_density)/sqrt(n()),
    syb_density = mean(syb_density),
    cumulative_fluo_se = sd(cumulative_fluo)/sqrt(n()),
    cumulative_fluo = mean(cumulative_fluo),
    mean_cumulative_fluo_se = sd(mean_cumulative_fluo)/sqrt(n()),
    mean_cumulative_fluo = mean(mean_cumulative_fluo),
    .groups="drop"
  )

## combine data
summary_respiration_redox <- summary_respiration %>%
  left_join(summary_redox,
            by = c("Exp_id","Day","Group")) %>%
  mutate(Soil_id = gsub("(S)[CW]", "\\1", Group)) %>%
  mutate(Site_type = case_when(Soil_id == "S1" ~ "Agriculture",
                               Soil_id == "S2" ~ "Grassland",
                               Soil_id == "S3" ~ "Forest")) %>%
  mutate(Season = case_when(
    Group == "SC1" ~ "Cold", Group == "SC2" ~ "Cold", Group == "SC3" ~ "Cold",
    Group == "SW1" ~ "Warm", Group == "SW2" ~ "Warm", Group == "SW3" ~ "Warm"))  %>%
  mutate(Soil = case_when(
    Group == "SC1" ~ "Cold - Agriculture", Group == "SC2" ~ "Cold - Grassland",
    Group == "SC3" ~ "Cold - Forest", Group == "SW1" ~ "Warm - Agriculture",
    Group == "SW2" ~ "Warm - Grassland", Group == "SW3" ~ "Warm - Forest"))





## plot correlation ####

season_shape <- c(21,24)
site_color <- c("Agriculture"="#b35806","Grassland"="#762a83","Forest"="#01665e")
day_size <- c("0"=1.8,"3"=2.2,"6"=2.6,"9"=3.0,
              "15"=3.4,"21"=3.8,"27"=5)
group_order <- c("Agriculture","Grassland","Forest")

### Community median redox activity ~ CO2 emission ####
p4.1 <- ggplot(summary_respiration_redox, aes(x = median_fluorescence, y = initial_rate/15)) +
  geom_point(aes(color = Site_type, shape = Season, size = factor(Day)), alpha = 0.7) +
  geom_errorbar(aes(ymin = ci_low_R0/15, ymax = ci_high_R0/15, color = Site_type), width = 0) +
  geom_errorbarh(aes(xmin = median_fluorescence - median_fluorescence_se, 
                     xmax = median_fluorescence + median_fluorescence_se, color = Site_type), height = 0) +
  geom_smooth(method = "lm", se = FALSE, color = "grey10", alpha = 0.2, linetype = 2) +
  geom_ribbon(stat = "smooth", method = "lm", alpha=0.05, color="grey30", linetype = 2)+
  #stat_regline_equation(label.y = 87, size = 5)+
  stat_cor(method = "spearman", size = 5, label.y = 85, p.accuracy = 0.01,
           aes(label = paste(gsub("italic(R)", "italic(rho)", after_stat(r.label), fixed = TRUE), 
                             gsub("italic(p)", "italic(P)", after_stat(p.label), fixed = TRUE), 
                             sep = "~`,`~")))+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))+
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)))+
  #facet_wrap(~Season)+
  labs(y = expression("CO"[2]*" emission ("~ppm~h^{-1}*" g"^{-1}*")"),
       x = "Fluorescence intensity (a.u.)",
       title = "Median redox activity")+
  scale_fill_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
                    values = c("#b35806","#762a83","#01665e"))+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  scale_size_manual(name="Day",values=day_size)+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black'),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size=13, color='black'), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 15, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)),
         size = guide_legend(title = "Day"))

p4.1

### Community active-cell abundance ~ CO2 emission ####
p4.2 <- ggplot(summary_respiration_redox, aes(x = rsg_density, y = initial_rate/15)) +
  geom_point(aes(color = Site_type, shape = Season, size = factor(Day)), alpha = 0.7) +
  geom_errorbar(aes(ymin = ci_low_R0/15, ymax = ci_high_R0/15, color = Site_type), width = 0) +
  geom_errorbarh(aes(xmin = rsg_density - rsg_density_se, 
                     xmax = rsg_density + rsg_density_se, color = Site_type), height = 0) +
  geom_smooth(method = "lm", se = FALSE, color = "grey10", alpha = 0.2, linetype = 1) +
  geom_ribbon(stat = "smooth", method = "lm", alpha=0.05, color="grey30", linetype = 2)+
  #stat_regline_equation(label.y = 87, size = 5)+
  stat_cor(method = "spearman", size = 5, label.y = 85, p.accuracy = 0.001,
           aes(label = paste(gsub("italic(R)", "italic(rho)", after_stat(r.label), fixed = TRUE), 
                             gsub("italic(p)", "italic(P)", after_stat(p.label), fixed = TRUE), 
                             sep = "~`,`~")))+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))+
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)))+
  #facet_wrap(~Site_type)+
  labs(y = expression("CO"[2]*" emission ("~ppm~h^{-1}*" g"^{-1}*")"),
       x = expression("Cell abundance (g"^{-1}*" soil)"),
       title = "Active cell abundance")+
  scale_fill_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
                    values = c("#b35806","#762a83","#01665e"))+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  scale_size_manual(name="Day",values=day_size)+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black'),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size=13, color='black'), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 15, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)),
         size = guide_legend(title = "Day"))

p4.2

### Community cumulative activity ~ CO2 emission ####
p4.3 <- ggplot(summary_respiration_redox, aes(x = cumulative_fluo, y = initial_rate/15)) +
  geom_point(aes(color = Site_type, shape = Season, size = factor(Day)), alpha = 0.7) +
  geom_errorbar(aes(ymin = ci_low_R0/15, ymax = ci_high_R0/15, color = Site_type), width = 0) +
  geom_errorbarh(aes(xmin = cumulative_fluo - cumulative_fluo_se, 
                     xmax = cumulative_fluo + cumulative_fluo_se, color = Site_type), height = 0) +
  geom_smooth(method = "lm", se = FALSE, color = "grey10", alpha = 0.2, linetype = 1) +
  geom_ribbon(stat = "smooth", method = "lm", alpha=0.05, color="grey30", linetype = 2)+
  #stat_regline_equation(label.y = 87, size = 5)+
  stat_cor(method = "spearman", size = 5, label.y = 85, p.accuracy = 0.001,
           aes(label = paste(gsub("italic(R)", "italic(rho)", after_stat(r.label), fixed = TRUE), 
                             gsub("italic(p)", "italic(P)", after_stat(p.label), fixed = TRUE), 
                             sep = "~`,`~")))+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))+
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)))+
  #facet_wrap(~Site_type)+
  labs(y = expression("CO"[2]*" emission ("~ppm~h^{-1}*" g"^{-1}*")"),
       x = expression("Fluorescence intensity (g"^{-1}*" soil)"),
       title = "Cumulative activity")+
  scale_fill_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
                    values = c("#b35806","#762a83","#01665e"))+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  scale_size_manual(name="Day",values=day_size)+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black'),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size=13, color='black'), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 15, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)),
         size = guide_legend(title = "Day"))

p4.3

### cowplot ####

f1.d <- plot_grid(
  p4.1, p4.2, p4.3,
  ncol=3, nrow = 1,
  align="hv",
  axis="tblr",
  rel_widths=c(1,1,1),
  rel_heights=c(1,1,1)
)

f1.d




#*********************####
# Supplementary ####

### plot respiration emission ####
p4.6 <- ggplot(df_respiration, aes(x = Day, y = initial_rate/15, colour = Site_type, fill = Site_type, group = Group)) +
  stat_pointinterval(.width = c(0.66, 0.95), point_interval = "median_hdci",
                     point_size = 3.5 , point_alpha = 0.6,
                     interval_alpha = 0.6,
                     position = position_dodge(width = 1)) +
  geom_line(aes(linetype = Season), stat = "summary", fun = median, linewidth = 1.5, alpha = 0.6,
            position = position_dodge(width = 1)) +
  #facet_wrap(~group, scales = "free_y")+
  scale_x_continuous(name = "Day after disturbance", breaks = seq(0, 30, 3), expand = c(0, 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.3)))+
  labs(x = "Day after disturbance",
       y = expression("CO"[2]*" emission ("~ppm~h^{-1}*" g"^{-1}*")")) +
  scale_color_manual(name = c("Site"),
                     limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  scale_fill_manual(name = c("Site"),
                    limits = c("Agriculture", "Grassland", "Forest"),
                    values = c("#b35806","#762a83","#01665e"))+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black'),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size=13, color='black'), 
        axis.text.y=element_text(size=11, color='black')) +
  theme(strip.text = element_text(size = 15, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "right")+
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)),
         linetype = guide_legend(title = "Season"))

p4.6

ggsave("Figure_SI_co2_emission.jpeg", dpi=600, width = 18, height = 10, units = "cm")


## VPA ####
vpa <- varpart(
  summary_respiration_redox$initial_rate,
  summary_respiration_redox %>%
    select(median_fluorescence),
  summary_respiration_redox %>%
    select(rsg_density)
)

plot(vpa)

summary(vpa)


vpa <- varpart(
  summary_respiration_redox$initial_rate,
  summary_respiration_redox %>%
    select(mean_cumulative_fluo),
  summary_respiration_redox %>%
    select(syb_density)
)

plot(vpa)

summary(vpa)



## plot correlation ####
p4.4 <- ggplot(summary_respiration_redox, aes(x = mean_cumulative_fluo, y = initial_rate/15)) +
  geom_point(aes(color = Site_type, shape = Season, size = as.factor(Day)),alpha = 0.7) +
  geom_errorbar(aes(ymin = ci_low_R0/15, ymax = ci_high_R0/15, color = Site_type), width = 0) +
  geom_errorbarh(aes(xmin = mean_cumulative_fluo - mean_cumulative_fluo_se, 
                     xmax = mean_cumulative_fluo + mean_cumulative_fluo_se, color = Site_type), height = 0) +
  geom_smooth(method = "lm", se = FALSE, color = "grey10", alpha = 0.2, linetype = 2) +
  geom_ribbon(stat = "smooth", method = "lm", alpha=0.05, color="grey30", linetype = 2)+
  #stat_regline_equation(label.y = 87, size = 5)+
  stat_cor(method = "spearman", size = 5, label.y = 85, p.accuracy = 0.01,
           aes(label = paste(gsub("italic(R)", "italic(rho)", after_stat(r.label), fixed = TRUE), 
                             gsub("italic(p)", "italic(P)", after_stat(p.label), fixed = TRUE), 
                             sep = "~`,`~")))+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))+
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)))+
  labs(y = expression("CO"[2]*" emission ("~ppm~h^{-1}*" g"^{-1}*")"),
       x = "Fluorescence intensity (a.u.)",
       title = "Mean redox activity")+
  scale_fill_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
                    values = c("#b35806","#762a83","#01665e"))+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  scale_size_manual(name="Day",values=day_size)+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black'),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size=13, color='black'), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 15, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)),
         size = guide_legend(title = "Day"))

p4.4

p4.5 <- ggplot(summary_respiration_redox, aes(x = syb_density, y = initial_rate/15)) +
  geom_point(aes(color = Site_type, shape = Season, size = as.factor(Day)),alpha = 0.7) +
  geom_errorbar(aes(ymin = ci_low_R0/15, ymax = ci_high_R0/15, color = Site_type), width = 0) +
  geom_errorbarh(aes(xmin = syb_density - syb_density_se, 
                     xmax = syb_density + syb_density_se, color = Site_type), height = 0) +
  geom_smooth(method = "lm", se = FALSE, color = "grey10", alpha = 0.2, linetype = 2) +
  geom_ribbon(stat = "smooth", method = "lm", alpha=0.05, color="grey30", linetype = 2)+
  #stat_regline_equation(label.y = 87, size = 5)+
  stat_cor(method = "spearman", size = 5, label.y = 85, p.accuracy = 0.001,
           aes(label = paste(gsub("italic(R)", "italic(rho)", after_stat(r.label), fixed = TRUE), 
                             gsub("italic(p)", "italic(P)", after_stat(p.label), fixed = TRUE), 
                             sep = "~`,`~")))+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))+
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)))+
  #facet_wrap(~Site_type)+
  labs(y = expression("CO"[2]*" emission ("~ppm~h^{-1}*" g"^{-1}*")"),
       x = expression("Cell abundance (g"^{-1}*" soil)"),
       title = "Total cell abundance")+
  scale_fill_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
                    values = c("#b35806","#762a83","#01665e"))+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  scale_size_manual(name="Day",values=day_size)+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black'),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size=13, color='black'), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 15, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)),
         size = guide_legend(title = "Day"))

p4.5


p <- plot_grid(
  p4.4, p4.5,
  ncol=2, nrow = 1,
  align="hv",
  axis="tblr",
  rel_widths=c(1,1),
  rel_heights=c(1,1)
)

p

ggsave("Figure_SI_co2_vs_total_abundance.jpeg", dpi=600, width = 21, height = 10, units = "cm")

