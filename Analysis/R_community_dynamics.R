# R script for STB Project ####
# Multi-dimensional and multi-level responses ####

#********************************************####

library(readr)
library(tidyverse)

library(vegan)
library(phyloseq)

library(lmerTest)
library(performance)
library(emmeans)
library(ggpubr)

library(patchwork)
library(ggplot2)
library(cowplot)
library(ggalluvial)



# Load data ####

df_cell_count <- readRDS("logos/df_cell_count.rds")

df_redox <- readRDS("logos/df_redox.rds")

metadata_dna = readRDS(file = "phyloseq/metadata_dna.rds")




#*********************####
## Legend ####

season_shape <- c("Cold"=21, "Warm"=24)
site_color <- c("Agriculture"="#b35806","Grassland"="#762a83","Forest"="#01665e")
day_size <- c("0"=1.8,"3"=2.2,"6"=2.6,"9"=3.0,"15"=3.4,"21"=3.8,"27"=5)
season_linetype <- c("Cold"="solid","Warm"="dashed")

legend_df <- expand.grid(
  Site_type=names(site_color),Season=names(season_shape),Day=names(day_size)
) %>% mutate(x=1,xend=2,y=1,yend=1)

p_legend_v <- ggplot(legend_df,aes(x=x,y=y,xend=xend,yend=yend))+
  geom_segment(aes(linetype=Season),colour="black",linewidth=1)+
  geom_point(aes(x=1.5,colour=Site_type,fill=Site_type,size=Day,shape=Season),stroke=1)+
  scale_colour_manual(name="Site",values=site_color)+
  scale_fill_manual(name="Site",values=site_color)+
  scale_size_manual(name="Day",values=day_size)+
  scale_shape_manual(name="Season",values=season_shape)+
  scale_linetype_manual(name="Season",values=season_linetype)+
  guides(
    colour=guide_legend(order=1,override.aes=list(shape=21,size=5,fill=site_color,linetype=0)),
    fill="none",
    shape=guide_legend(order=3,override.aes=list(size=5,fill="black",colour="black")),
    linetype=guide_legend(order=3,override.aes=list(shape=c(21,24), size=5, fill=NA, colour="black")),
    size=guide_legend(order=4,override.aes=list(shape=21,colour="black",fill="black",linetype=0))
  )+
  theme_void()+
  theme(
    legend.position="right",
    legend.box="vertical",
    legend.text=element_text(size=14),
    legend.title=element_text(size=15),
    legend.key=element_blank(),
    legend.key.width=grid::unit(1.5,"cm")
  )

p_legend_v

legend_day <- get_legend(p_legend_v + guides(size=guide_legend(order=1),
                                          color="none",shape="none",linetype="none"))


# horizontal
p_legend_h <- ggplot(legend_df,aes(x=x,y=y,xend=xend,yend=yend))+
  geom_segment(aes(colour=Site_type, linetype=Season),linewidth=1)+
  geom_point(aes(x=1.5,colour=Site_type,fill=Site_type,size=Day),stroke=1)+
  scale_colour_manual(name="Site",values=site_color)+
  scale_fill_manual(name="Site",values=site_color)+
  scale_size_manual(name="Day",values=day_size)+
  scale_linetype_manual(name="Season",values=season_linetype)+
  guides(
    colour=guide_legend(order=1,override.aes=list(shape=21,size=5,fill=site_color,linetype=0)),
    fill="none",
    linetype=guide_legend(order=2,override.aes=list(colour="black",linewidth=1)),
    size="none"
  )+
  theme_void()+
  theme(
    legend.position="bottom",
    legend.box="horizontal",
    legend.text=element_text(size=15),
    legend.title=element_text(size=15),
    legend.key=element_blank(),
    legend.key.width=grid::unit(1,"cm")
  )

p_legend_h

legend_h <- cowplot::get_legend(p_legend_h)





#********************************************####

# calculate short-term (Day 3) and long-term (Day 27) resistance
calc_resistance <- function(df, value_col, group_cols = c("Group","Type"),
                            baseline = 0, early = 3, late = 27) {
  value_col <- rlang::ensym(value_col)
  day0_ref <- df %>%
    filter(Day == baseline) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(Day0 = mean(!!value_col, na.rm = TRUE), .groups = "drop")
  df %>%
    filter(Day %in% c(early, late)) %>%
    select(all_of(c(group_cols, "Replicate")), Day, value = !!value_col) %>%
    pivot_wider(names_from = Day, values_from = value, names_prefix = "Day") %>%
    left_join(day0_ref, by = group_cols) %>%
    mutate(
      resistance_Day3 = ifelse(Day0 > 0, 1 / (1 + abs(Day3 - Day0) / Day0), NA_real_),
      resistance_Day27 = ifelse(Day0 > 0, 1 / (1 + abs(Day27 - Day0) / Day0), NA_real_)
    ) %>%
    rename(value_Day0 = Day0, value_Day3 = Day3, value_Day27 = Day27)
}



# Cell abundance ####

## Dynamics of active cells ####
df_rsg <- df_cell_count %>%
  filter(Type == "Active") %>%
  mutate(density_rsg = density) %>%
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

### plot dynamics ####
p1.1 <- ggplot(df_rsg) +
  aes(x = Day, y = density, color = Site_type, fill = Site_type, group = interaction(Site_type, Season)) +
  stat_summary(aes(color = Site_type, linetype = Season), fun = mean,
               geom = "line", alpha = 0.7, linewidth = 1.5)+
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.15, color=NA) +
  scale_color_manual(name = c("Site"),
                     limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  scale_fill_manual(name = c("Site"),
                    limits = c("Agriculture", "Grassland", "Forest"),
                    values = c("#b35806","#762a83","#01665e"))+
  scale_x_continuous(name = "Day after disturbance", breaks = seq(0, 30, 3), expand = c(0, 0.5)) +
  scale_y_log10(name = expression("Cell abundance (g "^{-1}*" dry soil)"),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
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
        legend.position = "none")

p1.1



## Resistance ####

resist_rsg_count <- calc_resistance(
  df = df_rsg,
  value_col = density_rsg,
  group_cols = c("Group", "Type", "Replicate"))%>%
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


### plot ####
p1.2 <- ggplot(resist_rsg_count, aes(x = Season, y = resistance_Day3))+
  geom_boxplot(aes(group = Season), width = 0.4, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(shape = Season, colour = Site_type), width = 0.15, alpha = 0.7, size = 3) +
  scale_y_continuous(limits = c(0.3, 1.01))+
  scale_color_manual(name = c("Site"),
                     limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  scale_fill_manual(name = c("Site"),
                    limits = c("Agriculture", "Grassland", "Forest"),
                    values = c("#b35806","#762a83","#01665e"))+
  labs(y = "Resistance (Day 3)")+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_blank(),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size = 15, color = "black"), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 15, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 5)))

p1.2


p1.3 <- ggplot(resist_rsg_count, aes(x = Season, y = resistance_Day27))+
  geom_boxplot(aes(group = Season), width = 0.4, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(shape = Season, colour = Site_type), width = 0.15, alpha = 0.7, size = 3) +
  scale_y_continuous(limits = c(0.3, 1))+
  scale_color_manual(name = c("Site"),
                     limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  scale_fill_manual(name = c("Site"),
                    limits = c("Agriculture", "Grassland", "Forest"),
                    values = c("#b35806","#762a83","#01665e"))+
  labs(y = "Resistance (Day 27)")+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_blank(),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size = 15, color = "black"), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 15, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 5)))

p1.3


### linear model ####

# Day 3
model <- lm(
  resistance_Day3 ~ Season * Site_type,
  data = resist_rsg_count
)

anova(model)

r2(model)

pairs(emmeans(model, ~ Site_type | Season), adjust = "tukey")

# Day 27
model <- lm(
  resistance_Day27 ~ Season * Site_type,
  data = resist_rsg_count
)

anova(model)

r2(model)

pairs(emmeans(model, ~ Site_type | Season), adjust = "tukey")


## Temporal stability ####

# normalized temporal cv
stability_rsg_count <- df_cell_count %>%
  filter(Type %in% c("Active")) %>%
  group_by(Group, Replicate, Type) %>%
  summarise(
    mean = mean(density, na.rm = TRUE),
    sd = sd(density, na.rm = TRUE),
    CV = sd / mean,
    stability = ifelse(is.finite(CV), 1 / (1 + CV), NA_real_),
    mean_3_27 = mean(density[Day >= 3 & Day <= 27], na.rm = TRUE),
    sd_3_27 = sd(density[Day >= 3 & Day <= 27], na.rm = TRUE),
    CV_3_27 = sd_3_27 / mean_3_27,
    stability_3_27 = ifelse(is.finite(CV_3_27), 1 / (1 + CV_3_27), NA_real_),
    .groups = "drop"
  )

stability_rsg_count <- stability_rsg_count %>%
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

### plot ####
p1.4 <- ggplot(stability_rsg_count %>% filter(Type %in% c("Active")), 
               aes(x = Season, y = stability))+
  geom_boxplot(aes(group = Season), width = 0.4, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(shape = Season, colour = Site_type), width = 0.15, alpha = 0.7, size = 3) +
  scale_y_continuous(limits = c(0.5, 1))+
  scale_color_manual(name = c("Site"),
                     limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  scale_fill_manual(name = c("Site"),
                    limits = c("Agriculture", "Grassland", "Forest"),
                    values = c("#b35806","#762a83","#01665e"))+
  labs(y = "Temporal stability")+
  #stat_compare_means(method = "wilcox.test", comparisons = list(c("Cold", "Warm")),
  #                   label = "p.signif", label.x.npc = 0.4, size = 6, label.y.npc = 0.999) +
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_blank(),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size = 15, color = "black"), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 15, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 5)))

p1.4

### linear model ####

model <- lm(
  stability ~ Season * Site_type,
  data = stability_rsg_count
)

anova(model)

r2(model)

pairs(emmeans(model, ~ Site_type | Season), adjust = "tukey")







#*********************####
# Redox activity ####

## Dynamic ####

### plot ####
p2.1 <- ggplot(df_redox) +
  aes(x = Day, y = median_fluorescence, color = Site_type, group = interaction(Site_type,Season)) +
  stat_summary(aes(color = Site_type, linetype = Season, group = interaction(Site_type,Season)), fun = mean,
               geom = "line", alpha = 0.8, linewidth = 1.5)+
  stat_summary(aes(fill = Site_type, group = interaction(Site_type,Season)), fun.data = mean_se,
               geom = "ribbon", alpha = 0.15, color=NA) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))+
  scale_x_continuous(name = "Day after disturbance", breaks = seq(0, 30, 3), expand = c(0, 0.5)) +
  labs(y = "Median redox activity (a.u.)", x = "Day after disturbance") +
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
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 15, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")

p2.1 


## Resistance ####
resist_median_redox <- calc_resistance(
  df = df_redox,
  value_col = median_fluorescence,
  group_cols = c("Group", "Type", "Replicate"))%>%
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
  


### plot ####
p2.2 <- ggplot(resist_median_redox, aes(x = Season, y = resistance_Day3))+
  geom_boxplot(aes(group = Season), width = 0.4, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(shape = Season, colour = Site_type), width = 0.15, alpha = 0.7, size = 3) +
  scale_y_continuous(limits = c(0.3, 1))+
  scale_color_manual(name = c("Site"),
                     limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  scale_fill_manual(name = c("Site"),
                    limits = c("Agriculture", "Grassland", "Forest"),
                    values = c("#b35806","#762a83","#01665e"))+
  labs(y = "Resistance (Day 3)")+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_blank(),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size = 15, color = "black"), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 15, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 5)))

p2.2


p2.3 <- ggplot(resist_median_redox, aes(x = Season, y = resistance_Day27))+
  geom_boxplot(aes(group = Season), width = 0.4, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(shape = Season, colour = Site_type), width = 0.15, alpha = 0.7, size = 3) +
  scale_y_continuous(limits = c(0.3, 1))+
  scale_color_manual(name = c("Site"),
                     limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  scale_fill_manual(name = c("Site"),
                    limits = c("Agriculture", "Grassland", "Forest"),
                    values = c("#b35806","#762a83","#01665e"))+
  labs(y = "Resistance (Day 27)")+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_blank(),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size = 15, color = "black"), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 15, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 5)))

p2.3


### linear model ####

# Day3
model <- lm(
  resistance_Day3 ~ Season * Site_type,
  data = resist_median_redox
)

anova(model)

r2(model)

pairs(emmeans(model, ~ Site_type | Season), adjust = "tukey")

# Day27
model <- lm(
  resistance_Day27 ~ Season * Site_type,
  data = resist_median_redox
)

anova(model)

r2(model)

pairs(emmeans(model, ~ Site_type | Season), adjust = "tukey")


## Temporal stability ####

# normalized temporal cv
stability_median_redox <- df_redox %>%
  group_by(Group, Replicate, Type) %>%
  summarise(
    mean = mean(median_fluorescence, na.rm = TRUE),
    sd = sd(median_fluorescence, na.rm = TRUE),
    CV = sd / mean,
    stability = ifelse(is.finite(CV), 1 / (1 + CV), NA_real_),
    mean_3_27 = mean(median_fluorescence[Day >= 3 & Day <= 27], na.rm = TRUE),
    sd_3_27 = sd(median_fluorescence[Day >= 3 & Day <= 27], na.rm = TRUE),
    CV_3_27 = sd_3_27 / mean_3_27,
    stability_3_27 = ifelse(is.finite(CV_3_27), 1 / (1 + CV_3_27), NA_real_),
    .groups = "drop"
  )

stability_median_redox <- stability_median_redox %>%
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


p2.4 <- ggplot(stability_median_redox, aes(x = Season, y = stability))+
  geom_boxplot(aes(group = Season), width = 0.4, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(shape = Season, colour = Site_type), width = 0.15, alpha = 0.7, size = 3) +
  scale_y_continuous(limits = c(0.5, 1))+
  scale_color_manual(name = c("Site"),
                     limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  scale_fill_manual(name = c("Site"),
                    limits = c("Agriculture", "Grassland", "Forest"),
                    values = c("#b35806","#762a83","#01665e"))+
  labs(y = "Temporal stability")+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_blank(),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size = 15, color = "black"), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 15, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 5)))

p2.4

### linear model ####

model <- lm(
  stability ~ Season * Site_type,
  data = stability_median_redox
)

anova(model)

r2(model)

pairs(emmeans(model, ~ Site_type | Season), adjust = "tukey")






#*********************####
# Composition trajectory ####

ps <- readRDS("phyloseq/ps_otu_dna_GCN_corrected_absolute.rds")

## Dynamic ####
## PCoA
season_shape <- c(21,24)
site_color <- c("Agriculture"="#b35806","Grassland"="#762a83","Forest"="#01665e")
day_size <- c("0"=1.8,"3"=2.2,"6"=2.6,"9"=3.0,
              "15"=3.4,"21"=3.8,"27"=5)
group_order <- c("Agriculture","Grassland","Forest")

plot_list <- list(); pcoa_list <- list(); distance_list <- list()

for(g in group_order){
  ps_sub <- subset_samples(ps,Site_type==g) %>% prune_taxa(taxa_sums(.)>0,.)
  if(!ape::is.binary(phy_tree(ps_sub))) phy_tree(ps_sub) <- ape::multi2di(phy_tree(ps_sub))
  if(!is.null(phy_tree(ps_sub)$root.edge)) phy_tree(ps_sub)$root.edge <- NULL
  
  dist_sub <- phyloseq::distance(ps_sub,method="wunifrac",normalized=TRUE,fast=TRUE)
  
  meta_sub <- data.frame(sample_data(ps_sub)) %>%
    mutate(Sample_id=as.character(Sample_id),Day=factor(Day,levels=c("0","3","6","9","15","21","27")),
           Season=factor(Season),Trajectory=interaction(Group,Replicate,sep="_"))
  
  distance_list[[g]] <- as.matrix(dist_sub)
  ord_sub <- ape::pcoa(as.matrix(dist_sub))
  
  ord_df <- data.frame(ord_sub$vectors) %>%
    tibble::rownames_to_column("Sample_id") %>%
    select(Sample_id,Axis.1,Axis.2) %>%
    left_join(meta_sub,by="Sample_id") %>%
    arrange(Trajectory,Day)
  
  pcoa_list[[g]] <- ord_df
  
  Xaxis <- paste0("PCoA1 (",round(ord_sub$values$Relative_eig[1]*100,2),"%)")
  Yaxis <- paste0("PCoA2 (",round(ord_sub$values$Relative_eig[2]*100,2),"%)")
  
  plot_list[[g]] <- ggplot(ord_df,aes(Axis.1,Axis.2))+
    geom_point(aes(fill=Site_type,colour=Site_type,shape=Season,size=Day),
               stroke=1.1,alpha=.7)+
    scale_fill_manual(values=site_color,limits=names(site_color))+
    scale_color_manual(values=site_color,limits=names(site_color))+
    scale_size_manual(name="Day",values=day_size)+
    scale_shape_manual(name="Season",values=season_shape)+
    scale_linetype_manual(name="Season",values=c("Cold"="solid","Warm"="dashed"))+
    geom_hline(yintercept=0,linetype=2,colour="grey75")+
    geom_vline(xintercept=0,linetype=2,colour="grey75")+
    #coord_fixed(ratio=1)+
    labs(title=g,x=Xaxis,y=Yaxis,colour="Site type",fill="Site type")+
    theme_bw()+
    theme(panel.background=element_rect(fill="white",colour="black"),
          plot.title=element_text(size=17,hjust=.5),
          axis.title=element_text(size=17,colour="black"),
          axis.text.x=element_text(size=13,colour="black"),
          axis.text.y=element_text(size=13,colour="black"),
          legend.position = "none")
}


p3.1 <- cowplot::plot_grid(
  plot_list[[1]],plot_list[[2]],plot_list[[3]],legend_day,
  ncol=4, align="hv", axis="tblr", rel_widths=c(1,1,1,0.5)
)

p3.1





## Succession ####

meta <- metadata_dna %>%
  mutate(
    Sample_id=as.character(Sample_id),
    Day=as.numeric(as.character(Day)),
    Replicate=factor(Replicate),
    Trajectory=interaction(Group,Replicate,sep="_")
  )

wunifrac_dist <- readRDS(file = "phyloseq/wunifrac_dist_dna_gcn_absolute.rds")

wu_mat <- as.matrix(wunifrac_dist)

df_wu <- meta %>%
  arrange(Trajectory, Day) %>%
  group_by(Trajectory) %>%
  mutate(
    previous_sample = lag(Sample_id),
    delta_day = Day - lag(Day),
    step_distance = mapply(
      \(x, y) ifelse(is.na(y), NA,
                     wu_mat[x, y]),
      Sample_id, previous_sample),
    succession_rate = step_distance / delta_day,
    interval = paste0(lag(Day), "-", Day),
    Replicate = as.numeric(Replicate)
  ) %>%
  filter(!is.na(succession_rate)) %>%
  ungroup()

# stat
adonis_res <- vegan::adonis2(
  wunifrac_dist ~ Site_type * Season * factor(Day),
  data = meta,
  permutations = 999,
  strata = meta$Trajectory,
  by = "terms"
)

adonis_res




## Resistance and Temporal stability ####

resist_community <- meta %>%
  select(Sample_id, Trajectory, Day, Replicate, Group) %>%
  group_by(Trajectory) %>%
  summarise(
    Group = first(Group),
    Replicate = first(Replicate),
    Sample_Day0 = Sample_id[Day == 0][1],
    Sample_Day3 = Sample_id[Day == 3][1],
    Sample_Day27 = Sample_id[Day == 27][1],
    .groups = "drop"
  ) %>%
  mutate(
    dist_Day3 = mapply(\(x,y) wu_mat[x,y], Sample_Day0, Sample_Day3),
    dist_Day27 = mapply(\(x,y) wu_mat[x,y], Sample_Day0, Sample_Day27),
    resistance_Day3 = 1 / (1 + dist_Day3),
    resistance_Day27 = 1 / (1 + dist_Day27)
  )



# Trajectory temporal variability 
stability_traj <- df_wu %>%
  group_by(Group, Replicate, Site_type, Season, Trajectory) %>%
  summarise(
    mean = mean(step_distance, na.rm = TRUE),
    sd = sd(step_distance, na.rm = TRUE),
    CV = sd / mean,
    stability = ifelse(is.finite(CV), 1 / (1 + CV), NA_real_),
    mean_3_27 = mean(step_distance[Day >= 6 & Day <= 27], na.rm = TRUE),
    sd_3_27 = sd(step_distance[Day >= 6 & Day <= 27], na.rm = TRUE),
    CV_3_27 = sd_3_27 / mean_3_27,
    stability_3_27 = ifelse(is.finite(CV_3_27), 1 / (1 + CV_3_27), NA_real_),
    .groups = "drop"
  )





#*********************####

## Fig.1 cowplot ####

f1.a <- plot_grid(
  p1.1,p1.2,p1.3,p1.4,
  ncol=4,
  align="hv",
  axis="lb",
  rel_widths=c(0.6, 0.2, 0.2, 0.2)
)

f1.a


f1.b <- plot_grid(
  p2.1,p2.2,p2.3,p2.4,
  ncol=4,
  align="hv",
  axis="lb",
  rel_widths=c(0.6, 0.2, 0.2, 0.2)
)

f1.b



p3.1 <- cowplot::plot_grid(
  plot_list[[1]],plot_list[[2]],plot_list[[3]],legend_day,
  ncol=4, align="hv", axis="tblr",rel_widths=c(1,1,1,0.5)
)

p3.1










#********************************************####

# Change by interval (community) ####

## Community level
calc_change <- function(df, value_col, time_col = Day, group_cols = c("Group", "Type", "Replicate"),
                             log_transform = FALSE) {
  value_col <- rlang::ensym(value_col); time_col <- rlang::ensym(time_col)
  df %>% arrange(!!!rlang::syms(group_cols), !!time_col) %>%
    group_by(!!!rlang::syms(group_cols)) %>%
    mutate(
      time_start = lag(!!time_col), time_end = !!time_col,
      interval = ifelse(is.na(time_start), NA, paste0(time_start, "-", time_end)),
      value_start = lag(!!value_col), value_end = !!value_col,
      dt = time_end - time_start,
      change = value_end - value_start,
      #abs_change = abs(change),
      relative_change = change / value_start,
      relative_change_rate = relative_change / dt,
      relative_change_abs_rate = abs(relative_change_rate),
      log2FC = log2(value_end / value_start),
      rate = ((if (log_transform) log10(value_end) else value_end) -
                (if (log_transform) log10(value_start) else value_start)) / dt,
      abs_rate = abs(rate)
    ) %>%
    ungroup()
}

change_syb_count <- calc_change(
  df = df_syb,
  value_col = density_syb,
  time_col = Day,
  group_cols = c("Group", "Type", "Replicate")
)

change_rsg_count <- calc_change(
  df = df_rsg,
  value_col = density_rsg,
  time_col = Day,
  group_cols = c("Group", "Type", "Replicate")
)

change_redox_median <- calc_change(
  df = df_redox,
  value_col = median_fluorescence,
  time_col = Day,
  group_cols = c("Group", "Type", "Replicate")
)

### combine ####
change_combined <- change_syb_count %>%
  dplyr::select(Group, Season, Replicate, interval, 
                relative_change_rate_syb = relative_change_rate, 
                relative_change_syb = relative_change, log2FC_syb = log2FC) %>%
  inner_join(
    change_redox_median %>%
      dplyr::select(Group, Season, Replicate, interval, 
                    relative_change_rate_redox = relative_change_rate, 
                    relative_change_redox = relative_change, log2FC_redox = log2FC),
    by = c("Group", "Season", "Replicate", "interval")
  )%>%
  inner_join(
    change_rsg_count %>%
      dplyr::select(Group, Season, Replicate, interval, 
                    relative_change_rate_rsg = relative_change_rate,
                    relative_change_rsg = relative_change, log2FC_rsg = log2FC),
    by = c("Group", "Season", "Replicate", "interval")
  )%>%
  inner_join(
    df_wu %>%
      dplyr::select(Group, Season, Replicate, interval, relative_change_rate_com = succession_rate,
                    relative_change_com = step_distance),
    by = c("Group", "Season", "Replicate", "interval")
  )%>%
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
    Group == "SW3" ~ "Warm - Forest"))%>%
  filter(!is.na(relative_change_rate_syb), 
         !is.na(relative_change_rate_redox), 
         !is.na(relative_change_rate_com)) %>%
  mutate(Trajectory=interaction(Group,Replicate,sep="_"))


saveRDS(change_combined, "phyloseq/community_change_combined.rds")


### corr activity ~ community size ####
change_combined$interval <- factor(change_combined$interval, levels = c("0-3", "3-6", "6-9", "9-15", "15-21", "21-27"))

p5.1 <- ggplot(change_combined, aes(x = log2FC_rsg, y = log2FC_redox)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_point(color = "grey30", size = 3, alpha = 0.3) +
  geom_smooth(aes(color = Site_type, linetype = Season), 
              method = "lm", se = FALSE, alpha = 0.9, linewidth = 1.5) +
  stat_regline_equation(label.y.npc = 1, size = 5)+
  stat_cor(method = "pearson", size = 5, label.y.npc = 0.9, p.accuracy = 0.001,
           aes(label = paste(after_stat(r.label), 
                             gsub("italic(p)", "italic(P)", after_stat(p.label), fixed = TRUE), 
                             sep = "~`,`~")))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  #facet_wrap(~interval)+
  labs(x = "Change in\nactive cell abundance (by interval)",
       y = "Change in\nredox activity (by interval)",
       title = "Community level")+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black', hjust = 0.5),
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
         size = guide_legend(title = "Time interval (Day)"))

p5.1


model <- lmer(log2FC_redox ~ log2FC_rsg * Season * Site_type,
            data = change_combined)

anova(model)

r2(model)

summary(model)

emtrends(model, ~ Season * Site_type, var = "log2FC_rsg")





# Stability summary ####
# community stability index

p5.2 <- ggplot(filter(stability_com),
               aes(x = redox_stability, y = rsg_stability)) +
  geom_point(color = "grey30", size = 3, alpha = 0.3) +
  geom_smooth(color = "grey5", method = "lm", se = FALSE, alpha = 0.9, linewidth = 1.5) +
  stat_regline_equation(label.y.npc = 1, size = 5)+
  stat_cor(method = "pearson", size = 5, label.y.npc = 0.9, p.accuracy = 0.001,
           aes(label = paste(after_stat(r.label), 
                             gsub("italic(p)", "italic(P)", after_stat(p.label), fixed = TRUE), 
                             sep = "~`,`~")))+
  #facet_wrap(~Day) +
  theme_bw()+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  labs(x = "Temporal stability\n(active cell abundance)",
       y = "Temporal stability\n(redox activity)",
       title = "Community level")+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black', hjust = 0.5),
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
         size = guide_legend(title = "Time interval (Day)"))

p5.2



# @ R_population_dynamics.R for population changes











