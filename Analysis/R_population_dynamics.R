# R script for STB Project ####
# Multi-dimensional population responses ####

library(readr)
library(tidyverse)
library(dplyr)
library(tidyr)

library(vegan)
library(phyloseq)
library(lmerTest)
library(performance)
library(emmeans)
library(ggpubr)
library(dtw)

library(patchwork)
library(ggplot2)
library(cowplot)
library(ggalluvial)


#********************************************####

# Load data ####

df_cell_count <- readRDS("logos/df_cell_count.rds")

df_redox <- readRDS("logos/df_redox.rds")

dna_gcn_rel <- readRDS(file = "phyloseq/otu_dna_GCN_corrected_relative.rds")

metadata_dna = readRDS(file = "phyloseq/metadata_dna.rds")

metadata_dna <- metadata_dna %>%
  mutate(Microcosm = as.factor(paste(Group, Replicate, sep = "_")))


pred.GCN <- readRDS("phyloseq/DNA_OTU_GCN_prediction.rds")

pred.GCN <- pred.GCN$tab %>%
  dplyr::transmute(
    OTU=label,
    GCN=x,
    prob=probs
  )


#********************************************####

# Change by interval (population) ####
# rRNA vs rDNA relative ####

## Population level


## Load data
dna_gcn_rel <- readRDS(
  "phyloseq/otu_dna_GCN_corrected_relative.rds"
)
cdna_gcn_rel <- readRDS(
  "phyloseq/otu_cdna_GCN_corrected_relative.rds"
)
metadata_dna <- readRDS(
  "phyloseq/metadata_dna.rds"
)

otu_dna  <- as.matrix(dna_gcn_rel)
otu_cdna <- as.matrix(cdna_gcn_rel)

## Check and align taxa and samples
stopifnot(
  !is.null(rownames(otu_dna)),
  !is.null(colnames(otu_dna)),
  !is.null(rownames(otu_cdna)),
  !is.null(colnames(otu_cdna))
)

stopifnot(setequal(rownames(otu_cdna), rownames(otu_dna)))
stopifnot(setequal(colnames(otu_cdna), colnames(otu_dna)))

#otu_cdna <- otu_cdna[rownames(otu_dna), colnames(otu_dna), drop = FALSE]


## Align metadata
metadata_dna <- as.data.frame(metadata_dna)

metadata_dna <- metadata_dna %>%
  mutate(Microcosm = as.factor(paste(Group, Replicate, sep = "_")))


stopifnot(all(colnames(otu_dna) %in% rownames(metadata_dna)))

metadata_dna <- metadata_dna[colnames(otu_dna), , drop = FALSE]

stopifnot(all(rownames(metadata_dna) == colnames(otu_dna)))

## Verify one sample per microcosm per day
dup_samples <- metadata_dna %>%
  count(Microcosm, Day) %>%
  filter(n > 1)

if (nrow(dup_samples) > 0) {
  stop("Some Microcosm x Day combinations contain multiple samples.
       Check whether Replicate identifies the true microcosm.")
}

## Prevalence filter: retain OTUs detected in >= 2 samples
prev_cdna <- rowSums(otu_cdna > 0, na.rm = TRUE)
prev_dna  <- rowSums(otu_dna > 0, na.rm = TRUE)

keep_otu <- prev_cdna >= 2 & prev_dna >= 2

otu_cdna <- otu_cdna[keep_otu, , drop = FALSE]
otu_dna  <- otu_dna[keep_otu, , drop = FALSE]

stopifnot(nrow(otu_dna) > 1)

# to long format
otu_to_long <- function(mat, value_name) {
  as.data.frame(t(mat), check.names = FALSE) %>%
    tibble::rownames_to_column("Sample_id") %>%
    pivot_longer(-Sample_id, names_to = "OTU", values_to = value_name, values_drop_na = FALSE)
}

otu_dna_long <- otu_to_long(otu_dna, "rDNA") %>%
  left_join(metadata_dna,
            by = "Sample_id")

otu_cdna_long <- otu_to_long(otu_cdna, "rRNA") %>%
  select(Sample_id, OTU, rRNA = rRNA)

otu_long <- inner_join(otu_dna_long, otu_cdna_long,
                       by = c("Sample_id", "OTU")) %>%
  mutate(Trajectory = as.character(Microcosm)) %>%
  arrange(OTU, Trajectory, Day)%>%
  mutate(Day = as.numeric(as.character(Day)))




## Change by interval ####
## calc_otu_interval
calc_otu_interval <- function(df) {
  pseudocount = 0
  df %>%
    group_by(OTU, Trajectory) %>%
    arrange(Day, .by_group = TRUE) %>%
    mutate(
      Day_start = lag(Day), Day_end = Day,
      interval = if_else(is.na(Day_start), NA_character_,
                         paste0(Day_start, "-", Day_end)),
      dt = Day_end - Day_start,
      rDNA_start = lag(rDNA),
      rDNA_end = rDNA,
      rRNA_start = lag(rRNA),
      rRNA_end = rRNA,
      log2FC_rDNA = log2((rDNA + pseudocount) /
                           (lag(rDNA) + pseudocount)),
      log2FC_rRNA = log2((rRNA + pseudocount) /
                           (lag(rRNA) + pseudocount))
    ) %>%
    ungroup() %>%
    filter(!is.na(interval), dt > 0)
}

# calculation
otu_change_interval <- calc_otu_interval(otu_long)

# 
otu_change_interval_res <- otu_change_interval %>% 
  filter(!is.na(log2FC_rDNA), !is.na(log2FC_rRNA),
         rDNA_start > 0, rDNA_end > 0, rRNA_start > 0, rRNA_end > 0)%>%
  mutate(log_rDNA_start = log2(rDNA_start),
         log_rRNA_start = log2(rRNA_start),)




### linear mixed-effects model ####
model <- lmer(log2FC_rRNA ~ log2FC_rDNA * Season * Site_type +
                (1 | Trajectory) + (1 | interval) +
                (0 + log2FC_rDNA | Trajectory:interval),
              data = otu_change_interval_res, REML = TRUE)

VarCorr(model)

anova(model)

r2(model)

summary(model)

emtrends(model, ~ Season * Site_type, var = "Stability_all_rDNA")




p6.1 <- ggplot(otu_change_interval_res,
               aes(log2FC_rDNA, log2FC_rRNA)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hex(bins = 40) +
  scale_fill_gradient(trans = "log10", low = "grey90", high = "grey30", name = "OTU density")+
  geom_smooth(aes(color = Site_type, linetype = Season), 
              method = "lm", se = FALSE, alpha = 0.9, linewidth = 1.5) +
  stat_regline_equation(label.y.npc = 1, size = 5)+
  stat_cor(method = "pearson", size = 5, label.y.npc = 0.9, p.accuracy = 0.001,
           aes(label = paste(after_stat(r.label), 
                             gsub("italic(p)", "italic(P)", after_stat(p.label), fixed = TRUE), 
                             sep = "~`,`~")))+
  labs(x = "Change in\nrDNA abundance (by interval)",
       y = "Change in\nrRNA abundance (by interval)",
       title = "Population level")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  #facet_wrap(~interval)+
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

p6.1






## Stability ####
otu_cdna <- otu_cdna[keep_otu, , drop = FALSE]
otu_dna  <- otu_dna[keep_otu, , drop = FALSE]

table_DNA  <- as.matrix(dna_gcn_rel)
table_RNA <- as.matrix(cdna_gcn_rel)

## Prevalence filter: retain OTUs detected in >= 2 samples
prev_cdna <- rowSums(table_RNA > 0, na.rm = TRUE)
prev_dna  <- rowSums(table_DNA > 0, na.rm = TRUE)

keep_otu <- prev_cdna >= 2 & prev_dna >= 2

table_RNA <- table_RNA[keep_otu, , drop = FALSE]
table_DNA  <- table_DNA[keep_otu, , drop = FALSE]


# Check that all metadata samples exist in both tables
metadata_dna <- as.data.frame(metadata_dna)

metadata_dna <- metadata_dna %>%
  mutate(Microcosm = as.factor(paste(Group, Replicate, sep = "_")))

stopifnot(all(metadata_dna$Sample_id %in% colnames(table_DNA)))
stopifnot(all(metadata_dna$Sample_id %in% colnames(table_RNA)))


# rDNA table
df_DNA <- as.data.frame(otu_dna) %>%
  rownames_to_column("OTU") %>%
  pivot_longer(
    cols = -OTU,
    names_to = "Sample_id",
    values_to = "Abundance"
  ) %>%
  left_join(metadata_dna, by = "Sample_id") %>%
  mutate(Data_type = "rDNA")

df_RNA <- as.data.frame(otu_cdna) %>%
  rownames_to_column("OTU") %>%
  pivot_longer(
    cols = -OTU,
    names_to = "Sample_id",
    values_to = "Abundance"
  ) %>%
  left_join(metadata_dna, by = "Sample_id") %>%
  mutate(Data_type = "rRNA")

# Combine rDNA and rRNA
df_otu <- bind_rows(df_DNA, df_RNA) %>%
  filter(!is.na(Microcosm), !is.na(Day)) %>%
  mutate(Abundance = ifelse(is.nan(Abundance), NA_real_, Abundance))

# calculate stability
stability_otu <- df_otu %>%
  filter(Day %in% c(0, 3, 6, 9, 15, 21, 27)) %>%
  group_by(Data_type, Group, Microcosm, OTU) %>%
  summarise(
    n_nonNA = sum(!is.na(Abundance)),
    n_detected = sum(Abundance > 0, na.rm = TRUE),
    mean = mean(Abundance, na.rm = TRUE),
    sd = sd(Abundance, na.rm = TRUE),
    CV = sd / mean,
    Stability = ifelse(is.finite(CV), 1 / (1 + CV), NA_real_),
    .groups = "drop"
  )

stability_otu_res <- stability_otu %>%
  pivot_wider(
    id_cols = c(Group, Microcosm, OTU),
    names_from = Data_type,
    values_from = c(
      n_nonNA, n_detected,
      mean, sd, CV, Stability
    ),
    names_glue = "{.value}_{Data_type}"
  ) %>%
  left_join(
    metadata_dna %>%
      select(Group, Season, Site_type) %>%
      distinct(),
    by = "Group"
  )


# filter, only OTUs that have 7 time points are used
stability_otu_res <- stability_otu_res %>%
  filter(n_nonNA_rRNA == n_detected_rRNA)%>%
  filter(n_nonNA_rDNA == n_detected_rDNA)


## plot and filter
p6.2 <- ggplot(stability_otu_res,
               aes(Stability_rDNA, Stability_rRNA)) +
  geom_hex(bins = 40) +
  scale_fill_gradient(trans = "log10", low = "grey90", high = "grey30", name = "OTU density")+
  geom_smooth(aes(color = Site_type, linetype = Season), 
              method = "lm", se = FALSE, alpha = 0.9, linewidth = 1.5) +
  stat_regline_equation(label.y.npc = 0.99, size = 5)+
  stat_cor(method = "pearson", size = 5, label.y.npc = 0.88, p.accuracy = 0.001,
           aes(label = paste(after_stat(r.label), 
                             gsub("italic(p)", "italic(P)", after_stat(p.label), fixed = TRUE), 
                             sep = "~`,`~"))) +
  labs(x = "Temporal stability (rDNA)",
       y = "Temporal stability (rRNA)",
       title = "Population level")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
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

p6.2

### linear mixed-effects model ####
options(contrasts = c("contr.sum", "contr.poly"))

model <- lmer(Stability_rRNA ~ Stability_rDNA * Season * Site_type
              + (1 | OTU) + (1 | Microcosm) + (0 + Stability_rDNA | Microcosm),
            data = stability_otu_res, REML = TRUE)

anova(model)

performance::r2(model)

summary(model)

emtrends(model, ~ Season * Site_type, var = "Stability_rDNA")







## Trajectory asynchrony ####


target_days <- c(0, 3, 6, 9, 15, 21, 27)

otu_summary <- otu_long %>%
  filter(Day %in% target_days) %>%
  group_by(Group, Trajectory, OTU) %>%
  summarise(
    n_nonNA_rRNA = sum(is.finite(rRNA)),
    n_nonNA_rDNA = sum(is.finite(rDNA)),
    n_detected_rRNA = sum(is.finite(rRNA) & rRNA > 0),
    n_detected_rDNA = sum(is.finite(rDNA) & rDNA > 0),
    n_days = n_distinct(Day),
    .groups = "drop"
  )

otu_filter <- otu_summary %>%
  filter(
    n_days == length(target_days),
    n_nonNA_rRNA == n_detected_rRNA,
    n_nonNA_rDNA == n_detected_rDNA
  ) %>%
  select(Group, Trajectory, OTU)

otu_complete <- otu_long %>%
  semi_join(otu_filter, by = c("Group", "Trajectory", "OTU")) %>%
  filter(Day %in% target_days) %>%
  arrange(Group, Trajectory, OTU, Day)


# Optional: report how many OTUs remain
otu_filter_summary <- otu_complete %>%
  distinct(Group, Trajectory, OTU) %>%
  count(Group, Trajectory, name = "n_OTUs")

print(otu_filter_summary)


# Log-transform
otu_complete <- otu_complete %>%
  mutate(
    log_rDNA = log10(rDNA),
    log_rRNA = log10(rRNA)
  )





### calc_asynchrony ####
calc_asynchrony <- function(dat, window_size = 1) {
  
  dat <- dat %>% arrange(Day)
  
  # Initial abundance and activity
  initial_rDNA <- dat$rDNA[dat$Day == 0][1]
  initial_rRNA <- dat$rRNA[dat$Day == 0][1]
  initial_16S_ratio <- initial_rRNA / initial_rDNA
  
  # Remove paired missing observations
  dat <- dat %>%
    filter(
      is.finite(rDNA), is.finite(rRNA),
      rDNA > 0, rRNA > 0
    )
  
  day <- dat$Day
  x <- log10(dat$rDNA)
  y <- log10(dat$rRNA)
  
  n_days <- length(day)
  
  empty_result <- tibble(
    n_days = n_days,
    n_intervals = NA_integer_,
    initial_rDNA = initial_rDNA,
    initial_rRNA = initial_rRNA,
    initial_16S_ratio = initial_16S_ratio,
    directional_asynchrony = NA_real_,
    magnitude_asynchrony = NA_real_,
    temporal_asynchrony_days = NA_real_,
    mean_signed_offset_days = NA_real_,
    DTW_distance = NA_real_,
    DTW_path_length = NA_integer_
  )
  
  # Validate trajectory
  valid <- (
    n_days >= 3 &&
      all(day %in% target_days) &&
      !anyDuplicated(day) &&
      sd(x) > 0 &&
      sd(y) > 0
  )
  
  if (!valid) return(empty_result)
  
  dt <- diff(day)
  
  # Standardize trajectories separately
  x_z <- as.numeric(scale(x))
  y_z <- as.numeric(scale(y))
  
  # Changes per day
  vx <- diff(x_z) / dt
  vy <- diff(y_z) / dt
  
  # A. Directional asynchrony
  valid_direction <- (vx != 0 & vy != 0)
  
  directional_asynchrony <- if (any(valid_direction)) {
    mean(vx[valid_direction] * vy[valid_direction] < 0)
  } else {
    NA_real_
  }
  
  n_intervals <- sum(valid_direction)
  
  # B. Magnitude asynchrony
  magnitude_asynchrony <- mean(abs(abs(vx) - abs(vy)))
  
  # C. Temporal asynchrony and DTW distance
  fit <- tryCatch(
    dtw(
      x_z, y_z,
      window.type = "sakoechiba",
      window.size = window_size,
      distance.only = FALSE
    ),
    error = function(e) NULL
  )
  
  temporal_asynchrony_days <- NA_real_
  mean_signed_offset_days <- NA_real_
  dtw_distance <- NA_real_
  dtw_path_length <- NA_integer_
  
  if (
    !is.null(fit) &&
    is.finite(fit$distance) &&
    length(fit$index1) > 0 &&
    length(fit$index1) == length(fit$index2)
  ) {
    
    dtw_path_length <- length(fit$index1)
    dtw_distance <- fit$distance / dtw_path_length
    
    day_dna <- day[fit$index1]
    day_rna <- day[fit$index2]
    
    time_offset <- day_rna - day_dna
    
    temporal_asynchrony_days <- mean(abs(time_offset))
    mean_signed_offset_days <- mean(time_offset)
  }
  
  tibble(
    n_days = n_days,
    n_intervals = n_intervals,
    initial_rDNA = initial_rDNA,
    initial_rRNA = initial_rRNA,
    initial_16S_ratio = initial_16S_ratio,
    directional_asynchrony = directional_asynchrony,
    magnitude_asynchrony = magnitude_asynchrony,
    temporal_asynchrony_days = temporal_asynchrony_days,
    mean_signed_offset_days = mean_signed_offset_days,
    DTW_distance = dtw_distance,
    DTW_path_length = dtw_path_length
  )
}

# calculation
otu_asynchrony <- otu_complete  %>%
  group_by(Group, Trajectory, OTU) %>%
  group_modify(
    ~ calc_asynchrony(.x, window_size = 0)
  ) %>%
  ungroup() %>%
  arrange(Group, OTU, Trajectory)

otu_asynchrony <- otu_asynchrony %>%
  left_join(pred.GCN %>%
              dplyr::select(OTU, GCN),
            by = "OTU") %>%
  mutate(DTW_scaled = (DTW_distance - min(DTW_distance, na.rm = TRUE)) /
           (max(DTW_distance, na.rm = TRUE) - min(DTW_distance, na.rm = TRUE))) %>%
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

saveRDS(otu_asynchrony, "phyloseq/otu_interval_7_asynchrony.rds")


ggplot(otu_asynchrony, aes(x = log(initial_rDNA), y = DTW_scaled)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  geom_hex(bins = 40) +
  geom_smooth(aes(color=Site_type, linetype = Season), method = "lm", se = FALSE, 
              alpha = 0.9, linewidth = 1.5) +
  #geom_ribbon(stat = "smooth", method = "lm", alpha=0.2, color="grey30", linetype = 2)+
  stat_regline_equation(label.y.npc = 0.15, size = 5)+
  stat_cor(method = "pearson", size = 5, label.y.npc = 0.05, p.accuracy = 0.001,
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")))


ggplot(otu_asynchrony, aes(x = log(initial_16S_ratio), y = DTW_scaled)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  geom_hex(bins = 40) +
  geom_smooth(aes(color=Site_type, linetype = Season), method = "lm", se = FALSE, 
              alpha = 0.9, linewidth = 1.5) +
  #geom_ribbon(stat = "smooth", method = "lm", alpha=0.2, color="grey30", linetype = 2)+
  stat_regline_equation(label.y.npc = 0.15, size = 5)+
  stat_cor(method = "pearson", size = 5, label.y.npc = 0.05, p.accuracy = 0.001,
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")))



### plot Temporal trajectory asynchrony (DTW distance)  ####

library(scales)
library(ggridges)

otu_asynchrony$Soil <- factor(
  otu_asynchrony$Soil, 
  levels = c("Cold - Agriculture","Cold - Grassland", "Cold - Forest",
             "Warm - Agriculture", "Warm - Grassland", "Warm - Forest"))

p7.2 <- ggplot(otu_asynchrony,
  aes(x = DTW_scaled, y = Soil, fill = Soil)) +
  geom_vline(
    #xintercept = mean(otu_asynchrony$DTW_scaled, na.rm = TRUE),
    xintercept = 0.5,
    linetype = "dashed",colour = "black", linewidth = 0.6) +
  geom_density_ridges(
    scale = 1.5,
    rel_min_height = 0.02,
    alpha = 0.5,
    color = "white",
    linewidth = 0.3,
    bandwidth = 0.06) +
  labs(x = "Trajectory asynchrony", y = NULL)  +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=15),
        axis.title.x =element_text(size=13, color='black'),
        axis.title.y=element_text(size=13, color='black'),
        axis.text.x =element_text(size = 13, color = "black"), 
        axis.text.y=element_text(size=11, color='black')) +
  theme(strip.text = element_text(size = 15, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")

p7.2



#### statistic ####
library(sandwich)
library(lmtest)

# linear model
model <- lm(
  DTW_scaled ~ Season * Site_type,
  data = otu_asynchrony
)

# Model results
summary(model)

car::Anova(model, type = 2, white.adjust = "hc3")


# Pairwise comparisons
emmeans(model, pairwise ~ Season * Site_type, adjust = "tukey")




### visualization ####

selected_otus <- c("043d0dd93e8348f867854edfb1f54b43")


plot_data <- otu_complete %>%
  filter(
    OTU %in% selected_otus,
    Day %in% target_days
  ) %>%
  pivot_longer(
    cols = c(rDNA, rRNA),
    names_to = "Metric",
    values_to = "Abundance"
  ) %>%
  mutate(log_abundance = log10(Abundance)) %>%
  group_by(Group, Soil, Trajectory, OTU, Metric) %>%
  mutate(
    z_abundance = if (
      all(is.finite(log_abundance)) &&
      sd(log_abundance) > 0
    ) {
      as.numeric(scale(log_abundance))
    } else {
      NA_real_
    }
  ) %>%
  ungroup()

#Group means and SE across Microcosms
summary_data <- plot_data %>%
  group_by(Group, Soil, OTU, Day, Metric) %>%
  summarise(
    n = sum(is.finite(z_abundance)),
    mean_z = if (n > 0) mean(z_abundance, na.rm = TRUE) else NA_real_,
    se_z = if (n > 1) sd(z_abundance, na.rm = TRUE) / sqrt(n) else NA_real_,
    .groups = "drop"
  )


# plot
summary_data$Soil <- factor(
  summary_data$Soil, 
  levels = c("Cold - Agriculture","Cold - Grassland", "Cold - Forest",
             "Warm - Agriculture", "Warm - Grassland", "Warm - Forest"))

p7.1 <- ggplot(summary_data,
               aes(x = Day, y = mean_z, color = Metric, fill = Metric)) +
  geom_line(alpha = 0.55, linewidth = 0.6) +
  geom_ribbon(
    aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
    alpha = 0.05, color = NA, na.rm = TRUE) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  facet_wrap(~ Soil, ncol = 3) +
  scale_x_continuous(
    breaks = target_days,
    expand = expansion(mult = 0.05)) +
  scale_color_manual(values = c(rDNA = "#0072B2", rRNA = "#D55E00")) +
  scale_fill_manual(values = c(rDNA = "#0072B2", rRNA = "#D55E00")) +
  labs(x = "Day after disturbance", y = "Standardized relative abundance",
    color = NULL, fill = NULL, title = "Selected OTU (Desertimonas)") +
  theme_bw(base_size = 15) +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=15),
        axis.title.x =element_text(size=15, color='black', hjust = 0.5),
        axis.title.y=element_text(size=15, color='black'),
        axis.text.x =element_text(size=13, color='black'), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 13, color = "black"))+
  theme(legend.text=element_text(size=11, color='black'),
        legend.title=element_text(size=11, color='black'),
        legend.position = c(0.3, 0.75),
        legend.justification = c("right", "top"),
        legend.spacing.y = unit(0.01, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.background = element_rect(
          fill = NA, color = NA),
        legend.key = element_rect(fill = NA, color = NA))

p7.1


## Fig.3 cowplot ####

f3 <- plot_grid(
  p5.2,p5.1,p6.2,p6.1,
  ncol=2,
  align="hv",
  axis="lb",
  rel_widths=c(1, 1, 1, 1),
  labels = "auto",
  label_size = 22
)

f3


f3 <- plot_grid(
  f3, legend_h,
  ncol=1,
  align="hv",
  axis="lb",
  rel_widths=c(1, 0.9),
  rel_heights=c(1, 0.08)
)

f3

ggsave("Figure_3_01.jpeg", dpi=600, width = 26, height = 20, units = "cm")


f3.b <- plot_grid(
  p7.1, p7.2,
  ncol=2,
  align="hv",
  axis="b",
  rel_widths=c(0.75, 0.4),
  labels = c("e", "f"),
  label_size = 22
)

f3.b


ggsave("Figure_3_02.jpeg", dpi=600, width = 26, height = 10, units = "cm")










#********************************************####
# Population transitions ####

dna_gcn_abs <- readRDS("phyloseq/otu_dna_GCN_corrected_absolute.rds")

ratio_mat <- readRDS(file = "phyloseq/cdna_dna_ratio_mat.rds")
ratio_mat_clean <- ratio_mat
ratio_mat_clean[is.nan(ratio_mat_clean)] <- NA
ratio_mat_clean[is.infinite(ratio_mat_clean)] <- NA


table_a <- dna_gcn_abs

table_b <- ratio_mat_clean

keep_otu <- rowSums(table_a > 0, na.rm = TRUE) >= 2

table_a <- table_a[keep_otu, , drop = FALSE]

table_b <- table_b[keep_otu, , drop = FALSE]


## Calculation ####
### all dates
## check sample order always !!!
data.frame(
  table_a  = colnames(table_a),
  table_b = colnames(table_b)
) |> head(20)

all(rownames(table_a) == rownames(table_b))
all(colnames(table_a) == colnames(table_b))


groups <- unique(metadata_dna$Group)

transition_list <- list()

## Calculate DNA abundance and ratio dynamics
for(g in groups){
  
  ## Metadata and samples for one group
  meta_g <- metadata_dna[metadata_dna$Group == g, , drop = FALSE]
  meta_g$Day <- as.numeric(as.character(meta_g$Day))
  day_groups <- split(rownames(meta_g), meta_g$Day)
  
  ## Mean DNA abundance and ratio across biological replicates
  dna_day_mean <- sapply(
    day_groups, \(x) rowMeans(table_a[, x, drop = FALSE], na.rm = TRUE))
  
  ratio_day_mean <- sapply(
    day_groups, \(x) rowMeans(table_b[, x, drop = FALSE], na.rm = TRUE))
  
  ## Convert NaN generated by all-NA rows to NA
  dna_day_mean[is.nan(dna_day_mean)] <- NA
  ratio_day_mean[is.nan(ratio_day_mean)] <- NA
  
  ## Check whether Day 0 and Day 3 are available
  if(!all(c("0", "3") %in% colnames(dna_day_mean))){
    stop(paste("Day 0 or Day 3 is missing in Group:", g))
  }
  
  day0_samples <- rownames(meta_g)[meta_g$Day == 0]
  day3_samples <- rownames(meta_g)[meta_g$Day == 3]
  
  ## Day 0 and Day 3 values
  DNA_Day0 <- dna_day_mean[, "0"]
  DNA_Day3 <- dna_day_mean[, "3"]
  
  Ratio_Day0 <- ratio_day_mean[, "0"]
  Ratio_Day3 <- ratio_day_mean[, "3"]
  
  ## Absolute changes
  DNA_Day3_change <- ifelse(is.na(DNA_Day3), 0, DNA_Day3)
  Ratio_Day3_change <- ifelse(is.na(Ratio_Day3), 0, Ratio_Day3)
  
  DNA_change_D3 <- DNA_Day3_change - DNA_Day0
  Ratio_change_D3 <- Ratio_Day3_change - Ratio_Day0
  
  ## Calculate log2 fold changes
  
  pseudo_DNA <- min(dna_day_mean[dna_day_mean > 0],na.rm = TRUE) / 4
  pseudo_Ratio <- min(ratio_day_mean[ratio_day_mean > 0],na.rm = TRUE) / 4
  
  DNA_logFC_D3 <- log2((DNA_Day3 + pseudo_DNA) / (DNA_Day0 + pseudo_DNA))
  Ratio_logFC_D3 <- log2((Ratio_Day3 + pseudo_Ratio) / (Ratio_Day0 + pseudo_Ratio))
  
  DNA_resistance_Day3 = ifelse(DNA_Day0 > 0, 1 / (1 + abs(DNA_Day3 - DNA_Day0) / DNA_Day0), NA_real_)
  Ratio_resistance_Day3 = ifelse(Ratio_Day0 > 0, 1 / (1 + abs(Ratio_Day3 - Ratio_Day0) / Ratio_Day0), NA_real_)
  
  ## Statistical comparison between Day 0 and Day 3
  # Welch t-test
  OTU <- rownames(dna_day_mean)
  
  DNA_p <- sapply(seq_len(nrow(dna_day_mean)), function(i){
    x <- table_a[OTU[i], day0_samples]
    y <- table_a[OTU[i], day3_samples]
    
    x0 <- ifelse(is.finite(x), x, 0); y0 <- ifelse(is.finite(y), y, 0)
    
    if(all(x0 == 0) && sum(y0 > 0) >= 2) return(0.0001)## give a pvalue when Day0 is always zero, indicate activated
    if(all(y0 == 0) && sum(x0 > 0) >= 2) return(0.0001)## give a pvalue when Day3 is always zero, indicate loss
    
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    if(length(x) < 2 | length(y) < 2) return(NA_real_)
    x <- log2(x + pseudo_DNA); y <- log2(y + pseudo_DNA)
    
    if(var(c(x, y)) < .Machine$double.eps) return(NA_real_)
    tryCatch(t.test(x, y, var.equal = FALSE)$p.value, error = function(e) NA_real_)
  })
  
  Ratio_p <- sapply(seq_len(nrow(ratio_day_mean)), function(i){
    x <- table_b[OTU[i], day0_samples]
    y <- table_b[OTU[i], day3_samples]
    
    x0 <- ifelse(is.finite(x), x, 0); y0 <- ifelse(is.finite(y), y, 0)
    if(all(x0 == 0) && sum(y0 > 0) >= 2) return(0.0001)## give a pvalue when Day0 is always zero, indicate activated
    if(all(y0 == 0) && sum(x0 > 0) >= 2) return(0.0001)## give a pvalue when Day3 is always zero, indicate loss
    
    #x[!is.finite(x)] <- 0; y[!is.finite(y)] <- 0
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    if(length(x) < 2 | length(y) < 2) return(NA_real_)
    x <- log2(x + pseudo_Ratio); y <- log2(y + pseudo_Ratio)
    
    if(var(c(x, y)) < .Machine$double.eps) return(NA_real_)
    tryCatch(t.test(x, y, var.equal = FALSE)$p.value, error = function(e) NA_real_)
  })
  
  DNA_significance <- case_when(
    is.na(DNA_p) ~ "NA",
    DNA_p < 0.001 ~ "***",
    DNA_p < 0.01 ~ "**",
    DNA_p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  
  Ratio_significance <- case_when(
    is.na(Ratio_p) ~ "NA",
    Ratio_p < 0.001 ~ "***",
    Ratio_p < 0.01 ~ "**",
    Ratio_p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  
  ## Calculate temporal mean and stability
  ## Replace NA with zero only for stability calculations
  DNA_stability_data <- dna_day_mean
  Ratio_stability_data <- ratio_day_mean
  
  DNA_n_timepoints <- rowSums(DNA_stability_data > 0, na.rm = TRUE)
  Ratio_n_timepoints <- rowSums(Ratio_stability_data > 0, na.rm = TRUE)
  
  ### All timepoints
  # DNA NA = truly absent → 0; Ratio NA = undefined
  DNA_stability_data[is.na(DNA_stability_data)] <- 0
  
  ## Mean values and variation across all timepoints
  DNA_mean <- rowMeans(DNA_stability_data)
  DNA_sd <- apply(DNA_stability_data, 1, sd)
  
  Ratio_mean <- rowMeans(Ratio_stability_data, na.rm = TRUE)
  Ratio_sd <- apply(Ratio_stability_data, 1, sd, na.rm = TRUE)
  
  ## Normalize each ASV by its temporal mean
  DNA_CV <- DNA_sd / DNA_mean
  Ratio_CV <- Ratio_sd / Ratio_mean
  
  DNA_stability <- 1 / (1 + DNA_CV)
  Ratio_stability <- 1 / (1 + Ratio_CV)
  
  DNA_stability[!is.finite(DNA_stability)] <- NA
  Ratio_stability[!is.finite(Ratio_stability)] <- NA
  
  # minimum observation criterion
  Ratio_stability[Ratio_n_timepoints < 6] <- NA
  DNA_stability[DNA_n_timepoints < 6] <- NA
  
  
  ### Stability during Day 3–27
  day_3_27 <- colnames(dna_day_mean)[as.numeric(colnames(dna_day_mean)) >= 3]
  
  DNA_stability_3_27_data <- dna_day_mean[, day_3_27, drop = FALSE]
  Ratio_stability_3_27_data <- ratio_day_mean[, day_3_27, drop = FALSE]
  
  # DNA: NA = absent → 0
  DNA_stability_3_27_data[is.na(DNA_stability_3_27_data)] <- 0
  DNA_mean_3_27 <- rowMeans(DNA_stability_3_27_data)
  DNA_sd_3_27 <- apply(DNA_stability_3_27_data, 1, sd)
  
  # Ratio: NA = undefined → ignored
  Ratio_mean_3_27 <- rowMeans(Ratio_stability_3_27_data, na.rm = TRUE)
  Ratio_sd_3_27 <- apply(Ratio_stability_3_27_data, 1, sd, na.rm = TRUE)
  
  DNA_CV_3_27 <- DNA_sd_3_27 / DNA_mean_3_27
  Ratio_CV_3_27 <- Ratio_sd_3_27 / Ratio_mean_3_27
  
  DNA_stability_3_27 <- 1 / (1 + DNA_CV_3_27)
  Ratio_stability_3_27 <- 1 / (1 + Ratio_CV_3_27)
  
  DNA_stability_3_27[!is.finite(DNA_stability_3_27)] <- NA
  Ratio_stability_3_27[!is.finite(Ratio_stability_3_27)] <- NA
  
  DNA_n_timepoints_3_27 <- rowSums(DNA_stability_3_27_data > 0, na.rm = TRUE)
  Ratio_n_timepoints_3_27 <- rowSums(Ratio_stability_3_27_data > 0, na.rm = TRUE)
  
  # minimum observation criterion
  DNA_stability_3_27[DNA_n_timepoints_3_27 < 6] <- NA
  Ratio_stability_3_27[Ratio_n_timepoints_3_27 < 6] <- NA
  
  ## 7. Save results
  transition_list[[g]] <- data.frame(
    OTU = rownames(dna_day_mean),
    Group = g,
    
    DNA_Day0 = DNA_Day0,
    DNA_Day3 = DNA_Day3,
    DNA_mean = DNA_mean,
    DNA_change_D3 = DNA_change_D3,
    DNA_logFC_D3 = DNA_logFC_D3,
    DNA_resistance_Day3 = DNA_resistance_Day3,
    DNA_p = DNA_p,
    DNA_significance = DNA_significance,
    DNA_sd = DNA_sd,
    DNA_stability = DNA_stability,
    DNA_n_timepoints = DNA_n_timepoints,
    DNA_stability_3_27 = DNA_stability_3_27,
    DNA_mean_3_27 = DNA_mean_3_27,
    DNA_sd_3_27 = DNA_sd_3_27,
    DNA_n_timepoints_3_27 = DNA_n_timepoints_3_27,
    
    Ratio_Day0 = Ratio_Day0,
    Ratio_Day3 = Ratio_Day3,
    Ratio_mean = Ratio_mean,
    Ratio_change_D3 = Ratio_change_D3,
    Ratio_logFC_D3 = Ratio_logFC_D3,
    Ratio_resistance_Day3 = Ratio_resistance_Day3,
    Ratio_p = Ratio_p,
    Ratio_significance = Ratio_significance,
    Ratio_sd = Ratio_sd,
    Ratio_stability = Ratio_stability,
    Ratio_n_timepoints = Ratio_n_timepoints,
    Ratio_stability_3_27 = Ratio_stability_3_27,
    Ratio_mean_3_27 = Ratio_mean_3_27,
    Ratio_sd_3_27 = Ratio_sd_3_27,
    Ratio_n_timepoints_3_27 = Ratio_n_timepoints_3_27,
    
    stringsAsFactors = FALSE
  )
}

## Merge all groups
transition_results <- bind_rows(transition_list)

saveRDS(transition_results, "phyloseq/res_stability_16Sratio_abundance_abs.rds")


# data filter
transition_results <- readRDS("phyloseq/res_stability_16Sratio_abundance_abs.rds")

## remove zero ratios
transition_logFC <- transition_results %>%
  filter(Ratio_Day0 != 0 & Ratio_Day3 != 0)%>%
  filter(DNA_Day0 != 0 & DNA_Day3 != 0)

transition_logFC <- transition_logFC  %>%
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

transition_logFC$Soil <- factor(
  transition_logFC$Soil, 
  levels = c("Cold - Agriculture","Cold - Grassland", "Cold - Forest",
             "Warm - Agriculture", "Warm - Grassland", "Warm - Forest"))




#*********************####
# Changes in activity and abundance at Day 3 ####

## Ratio-state transitions day 0-3 ####

# Day 0 ratio classes and Day 3 ratio-state transitions

ratio_levels <- c("NA", "Zero", "Low", "Moderate", "High")

transition_results <- transition_results %>%
  mutate(
    Ratio_class_Day0 = case_when(
      is.na(Ratio_Day0) ~ "NA",
      Ratio_Day0 == 0 ~ "Zero",
      Ratio_Day0 < 0.25 ~ "Low",
      Ratio_Day0 < 1 ~ "Moderate",
      Ratio_Day0 >= 1 ~ "High"
    ),
    Ratio_class_Day3 = case_when(
      is.na(Ratio_Day3) ~ "NA",
      Ratio_Day3 == 0 ~ "Zero",
      Ratio_Day3 < 0.25 ~ "Low",
      Ratio_Day3 < 1 ~ "Moderate",
      Ratio_Day3 >= 1 ~ "High"
    ),
    Ratio_class_Day0 = factor(Ratio_class_Day0, levels = ratio_levels),
    Ratio_class_Day3 = factor(Ratio_class_Day3, levels = ratio_levels)
  )


### 16S-ratio transited taxa ####

# summary
significant_ratio_OTU <- transition_results %>%
  filter(!is.na(Ratio_p), Ratio_p < 0.05) %>%
  mutate(
    Ratio_change = case_when(
      Ratio_change_D3 > 0 ~ "Activated",
      Ratio_change_D3 < 0 ~ "Deactivated",
      TRUE ~ NA_character_),
    Ratio_change = factor(
      Ratio_change,
      levels = c("Deactivated", "Activated")),
    Ratio_change_class = paste(
      Ratio_class_Day0, "→", Ratio_class_Day3)
  ) %>%
  filter(!is.na(Ratio_change)) %>%
  left_join(pred.GCN %>%
              dplyr::select(OTU, GCN),
            by = "OTU")


#### count and taxonomy ####

tax_split <- readRDS("phyloseq/taxonomy_split_all_OTUs.rds") %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% tibble::rownames_to_column("OTU")

significant_ratio_OTU <- significant_ratio_OTU %>%
  left_join(tax_split, by = "OTU")

tax_level <- "Class"

tax_summary_sig_ratio <- significant_ratio_OTU %>%
  filter(!is.na(Ratio_change)) %>%
  count(Group, Ratio_change, Taxon = .data[[tax_level]], name = "n") %>%
  mutate(
    Taxon = case_when(
      is.na(Taxon) ~ "Others",
      grepl("Incertae[_ ]Sedis", Taxon, ignore.case = TRUE) ~ "Others",
      TRUE ~ Taxon
    )
  ) %>%
  group_by(Group, Ratio_change, Taxon) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Group, Ratio_change) %>%
  mutate(
    Total_n = sum(n),
    prop = n / Total_n,
    percent = prop * 100
  ) %>%
  ungroup() %>%
  arrange(Group, Ratio_change, desc(prop))

top_taxa <- tax_summary_sig_ratio %>%
  group_by(Taxon) %>%
  summarise(total_n = sum(n), .groups = "drop") %>%
  slice_max(total_n, n = 12) %>%
  pull(Taxon)

tax_summary_plot <- tax_summary_sig_ratio %>%
  mutate(Taxon = if_else(Taxon %in% top_taxa,Taxon,"Others")) %>%
  group_by(Group, Ratio_change, Taxon) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Group, Ratio_change) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()%>%
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

tax_summary_plot$Soil <- factor(
  tax_summary_plot$Soil, 
  levels = c("Cold - Agriculture","Cold - Grassland", "Cold - Forest",
             "Warm - Agriculture", "Warm - Grassland", "Warm - Forest"))

p1.1 <- ggplot(tax_summary_plot, aes(Soil, n, fill = Taxon)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.15) +
  facet_wrap(~Ratio_change, nrow = 1) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(),
    expand = expansion(mult = c(0, 0.05))) +
  scale_fill_brewer(palette = "Paired", name = tax_level) +
  labs(x = NULL, y = "Number of significant OTUs", fill = tax_level) +
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black'),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size=12, color='black', angle = 50, hjust = 1), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 17, color = "black"))+
  theme(legend.text=element_text(size=11, color='black'),
        legend.title=element_text(size=13, color='black'),
        legend.key = element_blank(),
        legend.position = "right")+
  guides(colour = guide_legend(nrow = 3, override.aes = list(size = 4, alpha = 1)))

p1.1



#### summary abundance in total count ####
ratio_sig_abundance <- significant_ratio_OTU %>%
  filter(!is.na(Ratio_change)) %>%
  group_by(Group, Ratio_change) %>%
  summarise(Abundance = sum(DNA_Day0, na.rm = TRUE), .groups = "drop") %>%
  ungroup()%>%
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


p1.2 <- ggplot(ratio_sig_abundance, aes(x = Season, y = Abundance, group = Site_type)) +
  geom_line(aes(colour = Site_type), linewidth = 0.8, alpha = 0.7) +
  geom_point(aes(colour = Site_type), size = 4) +
  facet_wrap(~Ratio_change, nrow = 1) +
  #scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::breaks_width(0.03),
  #                   expand = expansion(mult = c(0.05, 0.50))) +
  scale_color_manual(limits = c("Agriculture", "Grassland", "Forest"),
                     values = c("#b35806","#762a83","#01665e"))+
  labs(x = NULL, y = "Absolute abundance at Day 0", colour = "Site") +
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black'),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size=13, color='black'), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 17, color = "black"))+
  theme(legend.text=element_text(size=13, color='black'),
        legend.title=element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom")+
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4, alpha = 1)))

p1.2


#### summary transition class ####

ratio_levels <- c("NA", "Zero", "Low", "Moderate", "High")

ratio_sig_class <- significant_ratio_OTU %>%
  count(Group, Ratio_change, Ratio_class_Day0, Ratio_class_Day3, name = "n") %>%
  filter(!is.na(Ratio_class_Day0), !is.na(Ratio_class_Day3))%>%
  group_by(Group) %>%
  group_by(Group, Ratio_change) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()%>%
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


# plot by proportion
ratio_sig_class$Soil <- factor(
  ratio_sig_class$Soil, 
  levels = c("Cold - Agriculture","Cold - Grassland", "Cold - Forest",
             "Warm - Agriculture", "Warm - Grassland", "Warm - Forest"))

p1.3 <- ggplot(ratio_sig_class, aes(axis1 = Ratio_class_Day0, axis2 = Ratio_class_Day3, y = prop)) +
  geom_alluvium(aes(fill = Ratio_class_Day0), width = 0.14, alpha = 0.7) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.14, colour = "black") +
  #geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5, angle = 90) +
  scale_x_discrete(limits = c("Day 0", "Day 3"), expand = c(0.12, 0.12)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
  facet_wrap(~Soil+Ratio_change, nrow = 2, scales = "fixed") +
  scale_fill_manual(values = c("NA" = "#E0E0E0", Zero = "#8daad3", Low = "#4575b4", Moderate = "#fdae61", High = "#d73027"), drop = FALSE) +
  theme_bw(base_size = 14) +
  labs(x = NULL, y = "Proportion", fill = "Class of 16S ratio") +
  theme(panel.grid = element_blank(), 
        axis.title.x =element_text(size=17, color='black'),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size=15, color='black'), 
        axis.text.y=element_text(size=13, color='black'))+
  theme(strip.text = element_text(size = 15, color = "black"))+
  theme(legend.text=element_text(size=13, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "bottom")

p1.3






#*********************####
## Fig.4 cowplot ####
p.a <- plot_grid(
  p1.1,p1.2,
  ncol=2, nrow = 1,
  #align="hv",
  rel_widths=c(1,0.55),
  rel_heights=c(1,1)
)

p.a <- plot_grid(
  p.a,p1.3,
  ncol=1, nrow = 2,
  align="hv",
  #axis="tblr",
  rel_widths=c(1,1),
  rel_heights=c(1,1)
)

p.a

ggsave("Figure_4_ratio_transitions.jpeg", dpi=600, width = 33, height = 26, units = "cm")






