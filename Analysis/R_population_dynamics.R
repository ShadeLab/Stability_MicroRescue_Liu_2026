# R script for STB Project ####
# Multi-dimensional population responses ####

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

#********************************************####

# Change by interval (population) ####

## Population level

# filter extreme rare OTUs
metadata_dna = readRDS(file = "phyloseq/metadata_dna.rds")

dna_gcn_abs <- readRDS("phyloseq/otu_dna_GCN_corrected_absolute.rds")
my_dna_table <- dna_gcn_abs

ratio_mat <- readRDS(file = "phyloseq/cdna_dna_ratio_mat.rds")
ratio_mat_clean <- ratio_mat
ratio_mat_clean[is.nan(ratio_mat_clean)] <- NA
ratio_mat_clean[is.infinite(ratio_mat_clean)] <- NA
my_ratio_table <- ratio_mat_clean

keep_otu <- rowSums(my_dna_table > 0, na.rm = TRUE) >= 2

my_dna_table <- my_dna_table[keep_otu, , drop = FALSE]

my_ratio_table <- my_ratio_table[keep_otu, , drop = FALSE]


## Calculation ####
### all dates
## check sample order always !!!
data.frame(
  DNA  = colnames(my_dna_table),
  my_ratio_table = colnames(my_ratio_table)
) |> head(20)

all(rownames(my_dna_table) == rownames(my_ratio_table))
all(colnames(my_dna_table) == colnames(my_ratio_table))

groups <- unique(metadata_dna$Group)

otu_interval_list <- list()

for(g in groups){
  meta_g <- metadata_dna[metadata_dna$Group == g, , drop = FALSE]
  meta_g$Day <- as.numeric(as.character(meta_g$Day))
  day_groups <- split(rownames(meta_g), meta_g$Day)
  days <- sort(as.numeric(names(day_groups)))
  
  dna_mean <- sapply(day_groups, \(x) rowMeans(my_dna_table[, x, drop = FALSE], na.rm = TRUE))
  ratio_mean <- sapply(day_groups, \(x) rowMeans(my_ratio_table[, x, drop = FALSE], na.rm = TRUE))
  dna_mean[is.nan(dna_mean)] <- NA
  ratio_mean[is.nan(ratio_mean)] <- NA
  colnames(dna_mean) <- colnames(ratio_mean) <- as.character(days)
  
  ## Pseudocounts calculated within each Group
  #pseudo_DNA <- min(dna_mean[dna_mean > 0], na.rm = TRUE) / 4
  #pseudo_Ratio <- min(ratio_mean[ratio_mean > 0], na.rm = TRUE) / 4
  
  pseudo_DNA <- 0
  pseudo_Ratio <- 0
  
  ## Calculate all adjacent intervals
  for(i in seq_len(length(days) - 1)){
    d1 <- as.character(days[i]); d2 <- as.character(days[i + 1])
    DNA_start <- dna_mean[, d1]; 
    DNA_end <- dna_mean[, d2]
    Ratio_start <- ratio_mean[, d1]; 
    Ratio_end <- ratio_mean[, d2]
    
    otu_interval_list[[paste(g, d1, d2, sep = "_")]] <- data.frame(
      Group = g, OTU = rownames(dna_mean), Day_start = days[i], Day_end = days[i + 1],
      interval = paste0(days[i], "-", days[i + 1]), Duration = days[i + 1] - days[i],
      DNA_start = DNA_start, DNA_end = DNA_end,
      DNA_change = DNA_end - DNA_start,
      DNA_relative_change = (DNA_end - DNA_start) / (DNA_start + pseudo_DNA),
      DNA_log2FC = log2((DNA_end + pseudo_DNA) / (DNA_start + pseudo_DNA)),
      Ratio_start = Ratio_start, Ratio_end = Ratio_end,
      Ratio_change = Ratio_end - Ratio_start,
      Ratio_relative_change = (Ratio_end - Ratio_start) / (Ratio_start + pseudo_Ratio),
      Ratio_log2FC = log2((Ratio_end + pseudo_Ratio) / (Ratio_start + pseudo_Ratio)),
      DNA_detect_n_start = rowSums(my_dna_table[rownames(dna_mean), day_groups[[d1]], drop = FALSE] > 0, na.rm = TRUE),
      DNA_detect_n_end = rowSums(my_dna_table[rownames(dna_mean), day_groups[[d2]], drop = FALSE] > 0, na.rm = TRUE),
      Ratio_n_start = rowSums(!is.na(my_ratio_table[rownames(dna_mean), day_groups[[d1]], drop = FALSE])),
      Ratio_n_end = rowSums(!is.na(my_ratio_table[rownames(dna_mean), day_groups[[d2]], drop = FALSE])),
      stringsAsFactors = FALSE
    )
  }
}

otu_interval_res <- dplyr::bind_rows(otu_interval_list)

saveRDS(otu_interval_res, file = "phyloseq/OTU_DNA_Ratio_interval_summary_absolute.rds")




## Corr abundance ~ activity ####

otu_interval_res <- readRDS("phyloseq/OTU_DNA_Ratio_interval_summary_absolute.rds")

# remove values of zero
df_otu_interval <- otu_interval_res %>%
  dplyr::filter(
    is.finite(DNA_start), is.finite(DNA_end),
    is.finite(Ratio_start), is.finite(Ratio_end),
    DNA_start > 0, DNA_end > 0,
    Ratio_start > 0, Ratio_end > 0
  )%>%
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


df_otu_interval$interval <- factor(df_otu_interval$interval, levels = c("0-3", "3-6", "6-9", "9-15", "15-21", "21-27"))

df_otu_interval$Soil <- factor(
  df_otu_interval$Soil, 
  levels = c("Cold - Agriculture","Cold - Grassland", "Cold - Forest",
             "Warm - Agriculture", "Warm - Grassland", "Warm - Forest"))

p6.1 <- ggplot(df_otu_interval, aes(x = DNA_log2FC, y = Ratio_log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hex(bins = 40) +
  geom_smooth(aes(color=Site_type, linetype = Season), method = "lm", se = FALSE, 
              alpha = 0.9, linewidth = 1.5) +
  #geom_ribbon(stat = "smooth", method = "lm", alpha=0.2, color="grey30", linetype = 2)+
  stat_regline_equation(label.y.npc = 0.99, size = 5)+
  stat_cor(method = "pearson", size = 5, label.y.npc = 0.88,
           aes(label = paste(after_stat(r.label), 
                             gsub("italic(p)", "italic(P)", after_stat(p.label), fixed = TRUE), 
                             sep = "~`,`~")))+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))+
  labs(x = "Change in\nOTU absolute abundance",
       y = "Change in\nOTU 16S ratio",
       title = "Population level")+
  #facet_wrap(~interval, nrow = 2) +
  scale_fill_gradient(trans = "log10", low = "grey90", high = "grey30", name = "OTU density")+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
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
        legend.position = "right")+
  guides(color = guide_legend(override.aes = list(size = 5)))

p6.1

### linear model ####

model <- lmer(Ratio_log2FC ~ DNA_log2FC * Season * Site_type + (1 | OTU),
            data = df_otu_interval)

anova(model)

emtrends(model, ~ Season * Site_type, var = "DNA_log2FC")






## Lag effect by interval ####

# interval order
interval_info <- df_otu_interval %>%
  distinct(interval, Day_start, Day_end) %>%
  arrange(Day_start)

interval_levels <- interval_info$interval

# all lagged interval pairs
lagged_pairs <- df_otu_interval %>%
  select(Group, OTU, Day_start, Day_end, interval, Ratio_log2FC, DNA_log2FC) %>%
  rename(interval_a = interval, Day_start_a = Day_start, Day_end_a = Day_end,
         Ratio_a = Ratio_log2FC, DNA_a = DNA_log2FC) %>%
  inner_join(
    df_otu_interval %>%
      select(Group, OTU, Day_start, Day_end, interval, Ratio_log2FC, DNA_log2FC) %>%
      rename(interval_b = interval, Day_start_b = Day_start, Day_end_b = Day_end,
             Ratio_b = Ratio_log2FC, DNA_b = DNA_log2FC),
    by = c("Group", "OTU")
  ) %>%
  filter(Day_start_b >= Day_end_a) %>%
  mutate(
    interval_a = factor(interval_a, levels = interval_levels),
    interval_b = factor(interval_b, levels = interval_levels)
  )

# Four lagged relationships
lagged_long <- bind_rows(
  lagged_pairs %>% transmute(Group, OTU, interval_a, interval_b,
                             Relationship = "Activity → Abundance", X = Ratio_a, Y = DNA_b),
  lagged_pairs %>% transmute(Group, OTU, interval_a, interval_b,
                             Relationship = "Abundance → Activity", X = DNA_a, Y = Ratio_b),
  lagged_pairs %>% transmute(Group, OTU, interval_a, interval_b,
                             Relationship = "Activity → Activity", X = Ratio_a, Y = Ratio_b),
  lagged_pairs %>% transmute(Group, OTU, interval_a, interval_b,
                             Relationship = "Abundance → Abundance", X = DNA_a, Y = DNA_b)
)


# Correlation and P value within each Group
cor_by_group <- lagged_long %>%
  group_by(Relationship, interval_a, interval_b, Group) %>%
  summarise(
    n = sum(is.finite(X) & is.finite(Y)),
    cor = if (n >= 10 &&
              sd(X, na.rm = TRUE) > 0 &&
              sd(Y, na.rm = TRUE) > 0)
      cor(X, Y, method = "pearson", use = "complete.obs") else NA_real_,
    p = if (n >= 10 &&
            sd(X, na.rm = TRUE) > 0 &&
            sd(Y, na.rm = TRUE) > 0)
      cor.test(X, Y, method = "pearson")$p.value else NA_real_,
    .groups = "drop"
  )%>%
  mutate(
    p_adj = p.adjust(p, method = "BH")
  )

# General trend across soil groups
cor_heatmap <- cor_by_group %>%
  group_by(Relationship, interval_a, interval_b) %>%
  summarise(
    cor = median(cor, na.rm = TRUE),
    n_group = sum(!is.na(cor)),
    n_sig = sum(p_adj < 0.05, na.rm = TRUE),
    consistent = n_sig == 6,
    .groups = "drop"
  )

cor_heatmap$Relationship <- factor(
  cor_heatmap$Relationship,
  levels = c(
    "Activity → Abundance",
    "Abundance → Activity",
    "Activity → Activity",
    "Abundance → Abundance"
  )
)

p6.2 <- ggplot(cor_heatmap, aes(interval_b, interval_a, fill = cor)) +
  geom_tile(color = "grey85", linewidth = 0.6) +
  geom_text(
    data = filter(cor_heatmap, consistent),
    aes(label = sprintf("%.2f", cor)),
    fontface = "bold",
    size = 4
  ) +
  facet_wrap(~Relationship, ncol = 2) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1),
    name = "Pearson's R"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(
    x = "Subsequent interval",
    y = "Preceding interval"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, size = 13, color='black'),
    axis.text.y = element_text(size = 13, color='black'),
    axis.title = element_text(size = 17, color='black'),
    strip.text = element_text(size = 15, color='black'),
    legend.title = element_text(size = 15, color='black'),
    legend.text = element_text(size = 12, color='black'),
    panel.spacing = unit(0.8, "lines")
  )

p6.2




#*********************####
## Fig.2 cowplot ####

fig.2 <- plot_grid(
  p5.1, p6.1,
  ncol=2, nrow = 1,
  align="h",
  rel_widths=c(1,1.3),
  rel_heights=c(1,1)
)

fig.2

fig.2 <- plot_grid(
  fig.2, p6.2,
  ncol=1, nrow = 2,
  align="hv",
  rel_widths=c(1,1),
  rel_heights=c(1,1)
)

fig.2

ggsave("Figure_2.jpeg", dpi=600, width = 28, height = 24, units = "cm")






#*********************####
# Supplementary ####

## SI: absolute abundance across periods ####
p6.3 <- ggplot(df_otu_interval, aes(x = DNA_log2FC, y = Ratio_log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hex(bins = 40) +
  geom_smooth(aes(color=Site_type, linetype = Season), method = "lm", se = FALSE, 
              alpha = 0.9, linewidth = 1.5) +
  #geom_ribbon(stat = "smooth", method = "lm", alpha=0.2, color="grey30", linetype = 2)+
  stat_regline_equation(label.y.npc = 0.15, size = 5)+
  stat_cor(method = "pearson", size = 5, label.y.npc = 0.05, 
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  labs(x = "Change in\nOTU absolute abundance (cell count)",
       y = "Change in\nOTU 16S ratio (rRNA:rDNA)",
       title = "OTU absolute abundance vs 16S ratio")+
  facet_wrap(~interval, nrow = 2) +
  scale_fill_gradient(trans = "log10", low = "grey90", high = "grey30", name = "OTU density")+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
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
        legend.position = "right")+
  guides(color = guide_legend(override.aes = list(size = 5)))

p6.3



## SI: rRNA vs rDNA ####

dna_gcn_rel = readRDS(file="phyloseq/otu_dna_GCN_corrected_relative.rds")
cdna_gcn_rel = readRDS(file="phyloseq/otu_cdna_GCN_corrected_relative.rds")

metadata_dna = readRDS(file="phyloseq/metadata_dna.rds")

otu_dna  <- as.matrix(dna_gcn_rel)
otu_cdna <- as.matrix(cdna_gcn_rel)


## filter prevalence
prev_cdna <- rowSums(otu_cdna > 0, na.rm = TRUE)
prev_dna  <- rowSums(otu_dna > 0, na.rm = TRUE)

keep_otu <- (prev_cdna >= 2) & (prev_dna >= 2)

otu_cdna <- otu_cdna[keep_otu, , drop = FALSE]
otu_dna  <- otu_dna[keep_otu, , drop = FALSE]

# align samples and taxa
stopifnot(all(rownames(otu_cdna) == rownames(otu_dna)))
stopifnot(all(colnames(otu_cdna) %in% rownames(metadata_dna)))

metadata_dna <- metadata_dna[colnames(otu_dna), , drop = FALSE]
stopifnot(all(rownames(metadata_dna) == colnames(otu_cdna)))

metadata_dna$Day <- as.numeric(as.character(metadata_dna$Day))
groups <- unique(metadata_dna$Group)
otu_interval_list <- list()

for (g in groups) {
  meta_g <- metadata_dna[metadata_dna$Group == g, , drop = FALSE]
  day_groups <- split(rownames(meta_g), meta_g$Day)
  days <- sort(as.numeric(names(day_groups)))
  
  cdna_mean <- sapply(day_groups, \(s) rowMeans(otu_cdna[, s, drop = FALSE], na.rm = TRUE))
  dna_mean  <- sapply(day_groups, \(s) rowMeans(otu_dna[, s, drop = FALSE], na.rm = TRUE))
  
  if (is.null(dim(cdna_mean))) cdna_mean <- matrix(cdna_mean, ncol = 1)
  if (is.null(dim(dna_mean))) dna_mean <- matrix(dna_mean, ncol = 1)
  colnames(cdna_mean) <- colnames(dna_mean) <- as.character(days)
  cdna_mean[is.nan(cdna_mean)] <- NA
  dna_mean[is.nan(dna_mean)] <- NA
  
  # Same pseudocount for cDNA and DNA
  all_positive <- c(cdna_mean[cdna_mean > 0],
                    dna_mean[dna_mean > 0])
  
  #pseudocount <- min(all_positive, na.rm = TRUE) / 4
  pseudocount <- 0
  
  for (i in seq_len(length(days) - 1)) {
    d1 <- as.character(days[i])
    d2 <- as.character(days[i + 1])
    C1 <- cdna_mean[, d1]; C2 <- cdna_mean[, d2]
    D1 <- dna_mean[, d1]; D2 <- dna_mean[, d2]
    
    otu_interval_list[[paste(g, d1, d2, sep = "_")]] <- data.frame(
      Group = g, OTU = rownames(dna_mean),
      Day_start = days[i], Day_end = days[i + 1],
      interval = paste0(days[i], "-", days[i + 1]),
      Duration = days[i + 1] - days[i],
      cDNA_start = C1, cDNA_end = C2,
      DNA_start = D1, DNA_end = D2,
      delta_rRNA = log2(C2 + pseudocount) - log2(C1 + pseudocount),
      delta_rDNA = log2(D2 + pseudocount) - log2(D1 + pseudocount),
      cDNA_log2FC = log2((C2 + pseudocount) / (C1 + pseudocount)),
      DNA_log2FC = log2((D2 + pseudocount) / (D1 + pseudocount)),
      # ratio: no pseudocount
      Ratio_start = C1 / D1,
      Ratio_end   = C2 / D2,
      log_ratio_start = log2(C1 / D1),
      log_ratio_end   = log2(C2 / D2),
      Ratio_log2FC = log2(C2 / D2) - log2(C1 / D1),
      stringsAsFactors = FALSE
    )
  }
}

otu_interval_rel <- bind_rows(otu_interval_list)

otu_interval_rel <- otu_interval_rel %>%
  dplyr::filter(
    is.finite(DNA_start), is.finite(DNA_end),
    is.finite(cDNA_start), is.finite(cDNA_end),
    DNA_start > 0, DNA_end > 0,
    cDNA_start > 0, cDNA_end > 0
  )%>%
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

otu_interval_rel$interval <- factor(otu_interval_rel$interval, levels = c("0-3", "3-6", "6-9", "9-15", "15-21", "21-27"))


## change in rDNA and rRNA
ggplot(otu_interval_rel, aes(x = DNA_log2FC, y = cDNA_log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hex(bins = 40) +
  geom_smooth(aes(color=Site_type, linetype = Season), method = "lm", se = FALSE, 
              alpha = 0.9, linewidth = 1.5) +
  #geom_ribbon(stat = "smooth", method = "lm", alpha=0.2, color="grey30", linetype = 2)+
  stat_regline_equation(label.y.npc = 0.15, size = 5)+
  stat_cor(method = "pearson", size = 5, label.y.npc = 0.05, 
           aes(label = paste(after_stat(r.label), sep = "~`,`~")))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  labs(x = "Change in\nrelative abundance (rDNA)",
       y = "Change in\nrelative abundance (rRNA)",
       title = "rDNA vs rRNA")+
  facet_wrap(~interval, nrow = 2) +
  scale_fill_gradient(trans = "log10", low = "grey90", high = "grey30", name = "OTU density")+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
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
        legend.position = "right")+
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave("Figure_SI_corr_changes_rRNArDNA_across_periods.jpeg", dpi=600, width = 27, height = 16, units = "cm")



## change in rDNA and ratio
p6.4 <- ggplot(otu_interval_rel, aes(x = DNA_log2FC, y = Ratio_log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hex(bins = 40) +
  geom_smooth(aes(color=Site_type, linetype = Season), method = "lm", se = FALSE, 
              alpha = 0.9, linewidth = 1.5) +
  #geom_ribbon(stat = "smooth", method = "lm", alpha=0.2, color="grey30", linetype = 2)+
  stat_regline_equation(label.y.npc = 0.15, size = 5)+
  stat_cor(method = "pearson", size = 5, label.y.npc = 0.05, 
           aes(label = paste(after_stat(r.label), after_stat(p.label),  sep = "~`,`~")))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  labs(x = "Change in\nOTU relative abundance (rDNA)",
       y = "Change in\nOTU 16S ratio (rRNA:rDNA)",
       title = "OTU relative abundance vs 16S ratio")+
  facet_wrap(~interval, nrow = 2) +
  scale_fill_gradient(trans = "log10", low = "grey90", high = "grey30", name = "OTU density")+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
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
        legend.position = "right")+
  guides(color = guide_legend(override.aes = list(size = 5)))

p6.4

fig.s2 <- plot_grid(
  p6.3,p6.4, 
  ncol=1, nrow = 2,
  align="hv",
  rel_widths=c(1,1),
  rel_heights=c(1,1),
  labels = "auto",
  label_size = 18
)

fig.s2

ggsave("Figure_SI_corr_changes_rDNA16Sratio_across_periods.jpeg", dpi=600, width = 26, height = 30, units = "cm")






#********************************************####
# Stability of populations ####

## Calculation ####
groups <- unique(metadata_dna$Group)

stability_list <- list()

## Calculate DNA abundance and ratio dynamics
for(g in groups){
  
  ## Metadata and samples for one group
  meta_g <- metadata_dna[metadata_dna$Group == g, , drop = FALSE]
  meta_g$Day <- as.numeric(as.character(meta_g$Day))
  day_groups <- split(rownames(meta_g), meta_g$Day)
  
  ## Mean DNA abundance and ratio across biological replicates
  dna_day_mean <- sapply(
    day_groups, \(x) rowMeans(my_dna_table[, x, drop = FALSE], na.rm = TRUE))
  
  ratio_day_mean <- sapply(
    day_groups, \(x) rowMeans(my_ratio_table[, x, drop = FALSE], na.rm = TRUE))
  
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
    x <- my_dna_table[OTU[i], day0_samples]
    y <- my_dna_table[OTU[i], day3_samples]
    
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
    x <- my_ratio_table[OTU[i], day0_samples]
    y <- my_ratio_table[OTU[i], day3_samples]
    
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
  stability_list[[g]] <- data.frame(
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
stability_results <- bind_rows(stability_list)

saveRDS(stability_results, "phyloseq/res_stability_16Sratio_abundance_abs.rds")


# data filter
stability_results <- readRDS("phyloseq/res_stability_16Sratio_abundance_abs.rds")

## remove zero ratios
stability_logFC <- stability_results %>%
  filter(Ratio_Day0 != 0 & Ratio_Day3 != 0)%>%
  filter(DNA_Day0 != 0 & DNA_Day3 != 0)

stability_logFC <- stability_logFC  %>%
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

stability_logFC$Soil <- factor(
  stability_logFC$Soil, 
  levels = c("Cold - Agriculture","Cold - Grassland", "Cold - Forest",
             "Warm - Agriculture", "Warm - Grassland", "Warm - Forest"))





#********************************************####
# Changes in activity and abundance at Day 3 ####

## Ratio-state transitions day 0-3 ####

# Day 0 ratio classes and Day 3 ratio-state transitions

ratio_levels <- c("NA", "Zero", "Low", "Moderate", "High")

stability_results <- stability_results %>%
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
significant_ratio_OTU <- stability_results %>%
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
## Fig.3 cowplot ####
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

ggsave("Figure_3_ratio_transitions.jpeg", dpi=600, width = 33, height = 26, units = "cm")




#*********************####
# Supplementary ####

## SI: Abundance transited taxa ####

# summary
significant_abun_OTU <- stability_results %>%
  filter(!is.na(DNA_p), DNA_p < 0.05) %>%
  mutate(
    Abun_change = case_when(
      DNA_change_D3 > 0 ~ "Increased",
      DNA_change_D3 < 0 ~ "Decreased",
      TRUE ~ NA_character_),
    Abun_change = factor(
      Abun_change,
      levels = c("Decreased", "Increased"))) %>%
  filter(!is.na(Abun_change)) %>%
  left_join(pred.GCN %>%
              dplyr::select(OTU, GCN),
            by = "OTU")


#### count and taxonomy ####
tax_split <- readRDS("phyloseq/taxonomy_split_all_OTUs.rds") %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% tibble::rownames_to_column("OTU")

significant_abun_OTU <- significant_abun_OTU %>%
  left_join(tax_split, by = "OTU")

tax_level <- "Class"

tax_summary_sig_abun <- significant_abun_OTU %>%
  filter(!is.na(Abun_change)) %>%
  count(Group, Abun_change, Taxon = .data[[tax_level]], name = "n") %>%
  mutate(
    Taxon = case_when(
      is.na(Taxon) ~ "Others",
      grepl("Incertae[_ ]Sedis", Taxon, ignore.case = TRUE) ~ "Others",
      TRUE ~ Taxon
    )
  ) %>%
  group_by(Group, Abun_change, Taxon) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Group, Abun_change) %>%
  mutate(
    Total_n = sum(n),
    prop = n / Total_n,
    percent = prop * 100
  ) %>%
  ungroup() %>%
  arrange(Group, Abun_change, desc(prop))

top_taxa <- tax_summary_sig_abun %>%
  group_by(Taxon) %>%
  summarise(total_n = sum(n), .groups = "drop") %>%
  slice_max(total_n, n = 12) %>%
  pull(Taxon)

tax_summary_plot <- tax_summary_sig_abun %>%
  mutate(Taxon = if_else(Taxon %in% top_taxa,Taxon,"Others")) %>%
  group_by(Group, Abun_change, Taxon) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Group, Abun_change) %>%
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

p2.1 <- ggplot(tax_summary_plot, aes(Soil, n, fill = Taxon)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.15) +
  facet_wrap(~Abun_change, nrow = 1) +
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

p2.1

sig_abun_proportion <- tax_summary_sig_abun %>%
  group_by(Group, Abun_change) %>%
  summarise(n_OTU = sum(n), .groups = "drop") %>%
  mutate(Soil = recode(Group, SC1="Cold - Agriculture", SC2="Cold - Grassland", SC3="Cold - Forest", 
                       SW1="Warm - Agriculture", SW2="Warm - Grassland", SW3="Warm - Forest"))

#### summary abundance in total count ####
abun_sig_abundance <- significant_abun_OTU %>%
  filter(!is.na(Abun_change)) %>%
  group_by(Group, Abun_change) %>%
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

p2.2 <- ggplot(abun_sig_abundance, aes(x = Season, y = Abundance, group = Site_type)) +
  geom_line(aes(colour = Site_type), linewidth = 0.8, alpha = 0.7) +
  geom_point(aes(colour = Site_type), size = 4) +
  facet_wrap(~Abun_change, nrow = 1) +
  #scale_y_continuous(labels = scales::percent_format(accuracy = 1),
  #                   breaks = scales::breaks_width(0.03),
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

p2.2


#### cowplot ####
# Abundance transited taxa
p.b <- plot_grid(
  p2.1,p2.2,
  ncol=2, nrow = 1,
  #align="hv",
  rel_widths=c(1,0.6),
  rel_heights=c(1,1),
  labels = "auto",
  label_size = 18
)

p.b

ggsave("Figure_SI_abundance_transited_taxa.jpeg", dpi=600, width = 31, height = 12, units = "cm")



## SI: Corr initial states to short-term changes ####

#### 16S ratio
p1.1 <- ggplot(stability_logFC, aes(x = log2(Ratio_Day0), y = Ratio_logFC_D3)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  geom_hex(bins = 40) +
  geom_smooth(aes(color=Site_type, linetype = Season), method = "lm", se = FALSE, 
              alpha = 0.9, linewidth = 1.5) +
  #geom_ribbon(stat = "smooth", method = "lm", alpha=0.2, color="grey30", linetype = 2)+
  stat_regline_equation(label.y.npc = 0.15, size = 5)+
  stat_cor(method = "pearson", size = 5, label.y.npc = 0.05, p.accuracy = 0.001,
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  labs(y = "Sensitivity of 16S ratio",
       x = "16S ratio at Day 0")+
  #facet_wrap(~Soil, nrow = 2) +
  scale_fill_gradient(trans = "log10", low = "grey90", high = "grey30", name = "OTU density")+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
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
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 5)))

p1.1

p1.2 <- ggplot(stability_logFC, aes(x = log2(DNA_Day0), y = Ratio_logFC_D3)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  #geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hex(bins = 40) +
  geom_smooth(aes(color=Site_type, linetype = Season), method = "lm", se = FALSE, 
              alpha = 0.9, linewidth = 1.5) +
  #geom_ribbon(stat = "smooth", method = "lm", alpha=0.2, color="grey30", linetype = 2)+
  stat_regline_equation(label.y.npc = 0.15, size = 5)+
  stat_cor(method = "pearson", size = 5, label.y.npc = 0.05, p.accuracy = 0.001,
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  labs(y = "Sensitivity of 16S ratio",
       x = "Absolute abundance at Day 0")+
  #facet_wrap(~Soil, nrow = 2) +
  scale_fill_gradient(trans = "log10", low = "grey90", high = "grey30", name = "OTU density")+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
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
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 5)))

p1.2

p1.3 <- ggplot(stability_logFC, aes(x = log2(Ratio_Day0), y = DNA_logFC_D3)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  geom_hex(bins = 40) +
  geom_smooth(aes(color=Site_type, linetype = Season), method = "lm", se = FALSE, 
              alpha = 0.9, linewidth = 1.5) +
  #geom_ribbon(stat = "smooth", method = "lm", alpha=0.2, color="grey30", linetype = 2)+
  stat_regline_equation(label.y.npc = 0.15, size = 5)+
  stat_cor(method = "pearson", size = 5, label.y.npc = 0.05, p.accuracy = 0.001,
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  labs(y = "Sensitivity of absolute abundance",
       x = "16S ratio at Day 0")+
  #facet_wrap(~Soil, nrow = 2) +
  scale_fill_gradient(trans = "log10", low = "grey90", high = "grey30", name = "OTU density")+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
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
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 5)))

p1.3

p1.4 <- ggplot(stability_logFC, aes(x = log2(DNA_Day0), y = DNA_logFC_D3)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  #geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hex(bins = 40) +
  geom_smooth(aes(color=Site_type, linetype = Season), method = "lm", se = FALSE, 
              alpha = 0.9, linewidth = 1.5) +
  #geom_ribbon(stat = "smooth", method = "lm", alpha=0.2, color="grey30", linetype = 2)+
  stat_regline_equation(label.y.npc = 0.15, size = 5)+
  stat_cor(method = "pearson", size = 5, label.y.npc = 0.05, p.accuracy = 0.001,
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  labs(y = "Sensitivity of absolute abundance",
       x = "Absolute abundance at Day 0")+
  #facet_wrap(~Soil, nrow = 2) +
  scale_fill_gradient(trans = "log10", low = "grey90", high = "grey30", name = "OTU density")+
  scale_color_manual(name = "Site", limits = c("Agriculture", "Grassland", "Forest"),
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
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 5)))

p1.4
### cowplot ####

legend_3 <- get_legend(p1.1 + theme(legend.position = "bottom"))

p <- plot_grid(
  p1.1,p1.3,
  p1.2,p1.4,
  ncol=2, nrow = 2,
  align="hv",
  axis="tblr",
  rel_widths=c(1,1,1,1),
  rel_heights=c(1,1,1,1)
)

p

p <- plot_grid(
  legend_3, p,
  ncol=1, nrow = 2,
  align="hv",
  axis="tblr",
  rel_widths=c(1,1),
  rel_heights=c(0.1,1)
)

p

ggsave("Figure_SI_corr_initial_states_responses.jpeg", dpi=600, width = 26, height = 24, units = "cm")





## SI: compare class - ratio ####

significant_ratio_OTU <- significant_ratio_OTU %>%
  mutate(Soil = case_when(
    Group == "SC1" ~ "Cold - Agriculture", Group == "SC2" ~ "Cold - Grassland",
    Group == "SC3" ~ "Cold - Forest", Group == "SW1" ~ "Warm - Agriculture",
    Group == "SW2" ~ "Warm - Grassland", Group == "SW3" ~ "Warm - Forest"))

significant_ratio_OTU$Soil <- factor(
  significant_ratio_OTU$Soil, 
  levels = c("Cold - Agriculture","Cold - Grassland", "Cold - Forest",
             "Warm - Agriculture", "Warm - Grassland", "Warm - Forest"))

p1.1 <- ggplot(significant_ratio_OTU, aes(x = Ratio_change, y = Ratio_Day0)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.15) +
  geom_jitter(aes(color = Ratio_change), width = 0.15, alpha = 0.3) +
  facet_wrap(~Soil, nrow = 1) +
  labs(x = NULL, y = "16S ratio at Day 0", colour = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.4)))+
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Deactivated", "Activated")),
                     method.args = list(paired = FALSE, exact = FALSE),
                     label = "p.signif", label.x.npc = 0.4, size = 6, label.y.npc = 0.6) +
  scale_colour_brewer(palette = "Set1")+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black'),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size=13, color='black'), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 13, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")
p1.1

p1.2 <- ggplot(significant_ratio_OTU, aes(x = Ratio_change, y = DNA_Day0)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.15) +
  geom_jitter(aes(color = Ratio_change), width = 0.15, alpha = 0.3) +
  facet_wrap(~Soil, nrow = 1) +
  labs(x = NULL, y = "Absolute abundance at Day 0", colour = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.4)))+
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Deactivated", "Activated")),
                     method.args = list(paired = FALSE, exact = FALSE),
                     label = "p.signif", label.x.npc = 0.4, size = 6, label.y.npc = 0.6) +
  scale_colour_brewer(palette = "Set1")+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black'),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size=13, color='black'), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 13, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")
p1.2

p1.3 <- ggplot(significant_ratio_OTU, aes(x = Ratio_change, y = GCN)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.15) +
  geom_jitter(aes(color = Ratio_change), width = 0.15, alpha = 0.3) +
  facet_wrap(~Soil, nrow = 1) +
  labs(x = NULL, y = "16S rrn copy number", colour = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.4)))+
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Deactivated", "Activated")),
                     method.args = list(paired = FALSE, exact = FALSE),
                     label = "p.signif", label.x.npc = 0.4, size = 6, label.y.npc = 0.6) +
  scale_colour_brewer(palette = "Set1")+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black'),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size=13, color='black'), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 13, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")
p1.3

#### cowplot ####
p <- plot_grid(
  p1.1,p1.2,p1.3,
  ncol=1, nrow = 3,
  align="hv",
  axis="tblr",
  rel_widths=c(1,1,1),
  rel_heights=c(1,1,1)
)

p

ggsave("Figure_SI_compare_ratio_categories.jpeg", dpi=600, width = 34, height = 24, units = "cm")



## SI: compare class - abundance ####

significant_abun_OTU <- significant_abun_OTU %>%
  mutate(Soil = case_when(
    Group == "SC1" ~ "Cold - Agriculture", Group == "SC2" ~ "Cold - Grassland",
    Group == "SC3" ~ "Cold - Forest", Group == "SW1" ~ "Warm - Agriculture",
    Group == "SW2" ~ "Warm - Grassland", Group == "SW3" ~ "Warm - Forest"))

significant_abun_OTU$Soil <- factor(
  significant_abun_OTU$Soil, 
  levels = c("Cold - Agriculture","Cold - Grassland", "Cold - Forest",
             "Warm - Agriculture", "Warm - Grassland", "Warm - Forest"))

p1.1 <- ggplot(significant_abun_OTU, aes(x = Abun_change, y = Ratio_Day0)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.15) +
  geom_jitter(aes(color = Abun_change), width = 0.15, alpha = 0.3) +
  facet_wrap(~Soil, nrow = 1) +
  labs(x = NULL, y = "16S ratio at Day 0", colour = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.4)))+
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Decreased", "Increased")),
                     method.args = list(paired = FALSE, exact = FALSE),
                     label = "p.signif", label.x.npc = 0.4, size = 6, label.y.npc = 0.6) +
  scale_colour_brewer(palette = "Set1")+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black'),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size=13, color='black'), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 13, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")
p1.1

p1.2 <- ggplot(significant_abun_OTU, aes(x = Abun_change, y = DNA_Day0)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.15) +
  geom_jitter(aes(color = Abun_change), width = 0.15, alpha = 0.3) +
  facet_wrap(~Soil, nrow = 1) +
  labs(x = NULL, y = "Absolute abundance at Day 0", colour = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.4)))+
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Decreased", "Increased")),
                     method.args = list(paired = FALSE, exact = FALSE),
                     label = "p.signif", label.x.npc = 0.4, size = 6, label.y.npc = 0.6) +
  scale_colour_brewer(palette = "Set1")+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black'),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size=13, color='black'), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 13, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")
p1.2

p1.3 <- ggplot(significant_abun_OTU, aes(x = Abun_change, y = GCN)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.15) +
  geom_jitter(aes(color = Abun_change), width = 0.15, alpha = 0.3) +
  facet_wrap(~Soil, nrow = 1) +
  labs(x = NULL, y = "16S rrn copy number", colour = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.4)))+
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Decreased", "Increased")),
                     method.args = list(paired = FALSE, exact = FALSE),
                     label = "p.signif", label.x.npc = 0.4, size = 6, label.y.npc = 0.6) +
  scale_colour_brewer(palette = "Set1")+
  theme_bw() +
  theme(panel.background=element_rect(fill='white', color='black'),
        plot.title = element_text(size=17),
        axis.title.x =element_text(size=17, color='black'),
        axis.title.y=element_text(size=17, color='black'),
        axis.text.x =element_text(size=13, color='black'), 
        axis.text.y=element_text(size=13, color='black')) +
  theme(strip.text = element_text(size = 13, color = "black"))+
  theme(legend.text=element_text(size=15, color='black'),
        legend.title=element_text(size=15, color='black'),
        legend.key = element_blank(),
        legend.position = "none")
p1.3

#### cowplot ####
p <- plot_grid(
  p1.1,p1.2,p1.3,
  ncol=1, nrow = 3,
  align="hv",
  axis="tblr",
  rel_widths=c(1,1,1),
  rel_heights=c(1,1,1)
)

p

ggsave("Figure_SI_compare_abundance_categories.jpeg", dpi=600, width = 34, height = 24, units = "cm")























