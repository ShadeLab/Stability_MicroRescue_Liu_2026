# R script for STB Project ####
# Indicative role of activity transitions ####

library(tidyverse)
library(seminr)


#********************************************####
# PLS-PM ####
# Partial Least Squares Path Modeling

## Data collection ####

all_metrics <- readRDS(all_metrics, "Community/Analysis/phyloseq/all_metrics_stability.rds")


## Multicollinearity analysis ####

inital_vars <- c("density_syb", "density_rsg", "median_fluorescence", "intensity_per_area", "mean_cumulative_fluo", "mean_fluorescence")

transition_act_vars <- c("Activated_count", "Activated_abun"
                         #"Deactivated_count", "Deactivated_abun"
)

transition_abun_vars <- c("Increased_count", "Increased_abun"
                          #"Decreased_count", "Decreased_abun"
)

stab_all_vars <- c("syb_stability", "rsg_stability", "redox_stability", "com_stability",
                   "syb_stability_3_27", "rsg_stability_3_27", "redox_stability_3_27", "com_stability_3_27")

resis_later_vars <- c("syb_resistance_Day27", "rsg_resistance_Day27", 
                      "redox_resistance_Day27", "com_resistance_Day27")


cor_matrix <- cor(all_metrics[, inital_vars], use = "complete.obs", method = "spearman")
cor_matrix







## Models ####


### model 1 ####

measurement_model <- constructs(
  composite("Initial_abundance", 
            c("density_rsg"),
            weights = mode_A),
  composite("Initial_activity",
            c("median_fluorescence"),
            weights = mode_A),
  composite("Transition_activity", 
            c("Activated_count", "Activated_count_pct", "Activated_abun"),
            weights = mode_A),
  composite("Transition_abundance", 
            c("Increased_count", "Increased_count_pct", "Increased_abun"),
            weights = mode_A),
  composite("Overall_community_stability", 
            c("syb_stability", "rsg_stability", "redox_stability", "com_stability"),
            weights = mode_A)
)

structural_model <- relationships(
  paths(from = "Initial_activity", to = "Transition_activity"),
  paths(from = "Initial_activity", to = "Transition_abundance"),
  paths(from = "Initial_activity", to = "Overall_community_stability"),
  paths(from = "Initial_abundance", to = "Transition_activity"),
  paths(from = "Initial_abundance", to = "Transition_abundance"),
  paths(from = "Initial_abundance", to = "Overall_community_stability"),
  paths(from = "Transition_abundance", to = "Overall_community_stability"),
  paths(from = "Transition_activity", to = "Overall_community_stability")
)

pls_model <- estimate_pls(
  data = all_metrics,
  measurement_model = measurement_model,
  structural_model = structural_model
)

summary(pls_model)

#plot_scores(pls_model)
#summary(pls_model)$validity

# Bootstrap
boot_model <- bootstrap_model(
  seminr_model = pls_model,
  nboot = 10000,
  cores = 6,
  seed = 123
)

boot_summary <- summary(boot_model, alpha = 0.05)

boot_summary$bootstrapped_paths


### model 2 ####
# Later community stability (composition)

measurement_model <- constructs(
  composite("Initial_abundance", 
            c("density_rsg"),
            weights = mode_A),
  composite("Initial_activity",
            c("median_fluorescence"),
            weights = mode_A),
  composite("Resistance_activity", 
            c("redox_resistance_Day3"),
            weights = mode_A),
  composite("Resistance_abundance", 
            c("rsg_resistance_Day3"),
            weights = mode_A),
  composite("Stability_later",
            c("com_stability_3_27"),
            weights = mode_A)
)

structural_model <- relationships(
  paths(from = "Initial_activity", to = "Resistance_activity"),
  paths(from = "Initial_activity", to = "Resistance_abundance"),
  paths(from = "Initial_activity", to = "Stability_later"),
  paths(from = "Initial_abundance", to = "Resistance_activity"),
  paths(from = "Initial_abundance", to = "Resistance_abundance"),
  paths(from = "Initial_abundance", to = "Stability_later"),
  paths(from = "Resistance_abundance", to = "Stability_later"),
  paths(from = "Resistance_activity", to = "Stability_later")
)

pls_model <- estimate_pls(
  data = all_metrics,
  measurement_model = measurement_model,
  structural_model = structural_model
)


summary(pls_model)


#plot_scores(pls_model)
#summary(pls_model)$validity

# Bootstrap
boot_model <- bootstrap_model(
  seminr_model = pls_model,
  nboot = 10000,          
  cores = 6,
  seed = 123
)

boot_summary <- summary(boot_model, alpha = 0.05)

boot_summary$bootstrapped_paths


### model 3 ####
# Later community stability (active-cell abundance)

measurement_model <- constructs(
  composite("Initial_abundance", 
            c("density_rsg"),
            weights = mode_A),
  composite("Initial_activity",
            c("median_fluorescence"),
            weights = mode_A),
  composite("Resistance_activity", 
            c("redox_resistance_Day3"),
            weights = mode_A),
  composite("Resistance_abundance", 
            c("rsg_resistance_Day3"),
            weights = mode_A),
  composite("Stability_later",
            c("rsg_stability_3_27"),
            weights = mode_A)
)

structural_model <- relationships(
  paths(from = "Initial_activity", to = "Resistance_activity"),
  paths(from = "Initial_activity", to = "Resistance_abundance"),
  paths(from = "Initial_activity", to = "Stability_later"),
  paths(from = "Initial_abundance", to = "Resistance_activity"),
  paths(from = "Initial_abundance", to = "Resistance_abundance"),
  paths(from = "Initial_abundance", to = "Stability_later"),
  paths(from = "Resistance_abundance", to = "Stability_later"),
  paths(from = "Resistance_activity", to = "Stability_later")
)

pls_model <- estimate_pls(
  data = all_metrics,
  measurement_model = measurement_model,
  structural_model = structural_model
)

summary(pls_model)
#plot_scores(pls_model)
#summary(pls_model)$validity

# Bootstrap
boot_model <- bootstrap_model(
  seminr_model = pls_model,
  nboot = 10000,          
  cores = 6,
  seed = 123
)

boot_summary <- summary(boot_model, alpha = 0.05)

boot_summary$bootstrapped_paths



### model 4 ####
# Later stability of community activity

measurement_model <- constructs(
  composite("Initial_abundance", 
            c("density_rsg"),
            weights = mode_A),
  composite("Initial_activity",
            c("median_fluorescence"),
            weights = mode_A),
  composite("Resistance_activity", 
            c("redox_resistance_Day3"),
            weights = mode_A),
  composite("Resistance_abundance", 
            c("rsg_resistance_Day3"),
            weights = mode_A),
  composite("Stability_later",
            c("redox_stability_3_27"),
            weights = mode_A)
)

structural_model <- relationships(
  paths(from = "Initial_activity", to = "Resistance_activity"),
  paths(from = "Initial_activity", to = "Resistance_abundance"),
  paths(from = "Initial_activity", to = "Stability_later"),
  paths(from = "Initial_abundance", to = "Resistance_activity"),
  paths(from = "Initial_abundance", to = "Resistance_abundance"),
  paths(from = "Initial_abundance", to = "Stability_later"),
  paths(from = "Resistance_abundance", to = "Stability_later"),
  paths(from = "Resistance_activity", to = "Stability_later")
)

pls_model <- estimate_pls(
  data = all_metrics,
  measurement_model = measurement_model,
  structural_model = structural_model
)

summary(pls_model)
#plot_scores(pls_model)
#summary(pls_model)$validity

# Bootstrap
boot_model <- bootstrap_model(
  seminr_model = pls_model,
  nboot = 10000,          
  cores = 6,
  seed = 123
)

boot_summary <- summary(boot_model, alpha = 0.05)

boot_summary$bootstrapped_paths





### model 5 ####
# Population level - abundance

cols <- c("DNA_stability_3_27", "Ratio_stability_3_27")

all_metrics_population <- stability_logFC[
  rowSums(is.finite(as.matrix(stability_logFC[cols]))) == length(cols), 
]

measurement_model <- constructs(
  composite("Initial_abundance", 
            c("DNA_Day0"),
            weights = mode_A),
  composite("Initial_activity",
            c("Ratio_Day0"),
            weights = mode_A),
  composite("Resistance_activity", 
            c("Ratio_resistance_Day3"),
            weights = mode_A),
  composite("Resistance_abundance", 
            c("DNA_resistance_Day3"),
            weights = mode_A),
  composite("Stability_later",
            c("DNA_stability_3_27"),
            weights = mode_A)
)

structural_model <- relationships(
  paths(from = "Initial_activity", to = "Resistance_activity"),
  paths(from = "Initial_activity", to = "Resistance_abundance"),
  paths(from = "Initial_activity", to = "Stability_later"),
  paths(from = "Initial_abundance", to = "Resistance_activity"),
  paths(from = "Initial_abundance", to = "Resistance_abundance"),
  paths(from = "Initial_abundance", to = "Stability_later"),
  paths(from = "Resistance_abundance", to = "Stability_later"),
  paths(from = "Resistance_activity", to = "Stability_later")
)


pls_model <- estimate_pls(
  data = all_metrics_population,
  measurement_model = measurement_model,
  structural_model = structural_model
)

summary(pls_model)
#plot_scores(pls_model)

# Bootstrap
boot_model <- bootstrap_model(
  seminr_model = pls_model,
  nboot = 10000,          
  cores = 6,
  seed = 123
)

boot_summary <- summary(boot_model, alpha = 0.05)

boot_summary$bootstrapped_paths


### model 6 ####
# Population level - activity

cols <- c("DNA_stability_3_27", "Ratio_stability_3_27")

all_metrics_population <- stability_logFC[
  rowSums(is.finite(as.matrix(stability_logFC[cols]))) == length(cols), 
]

measurement_model <- constructs(
  composite("Initial_abundance", 
            c("DNA_Day0"),
            weights = mode_A),
  composite("Initial_activity",
            c("Ratio_Day0"),
            weights = mode_A),
  composite("Resistance_activity", 
            c("Ratio_resistance_Day3"),
            weights = mode_A),
  composite("Resistance_abundance", 
            c("DNA_resistance_Day3"),
            weights = mode_A),
  composite("Stability_later",
            c("Ratio_stability_3_27"),
            weights = mode_A)
)

structural_model <- relationships(
  paths(from = "Initial_activity", to = "Resistance_activity"),
  paths(from = "Initial_activity", to = "Resistance_abundance"),
  paths(from = "Initial_activity", to = "Stability_later"),
  paths(from = "Initial_abundance", to = "Resistance_activity"),
  paths(from = "Initial_abundance", to = "Resistance_abundance"),
  paths(from = "Initial_abundance", to = "Stability_later"),
  paths(from = "Resistance_abundance", to = "Stability_later"),
  paths(from = "Resistance_activity", to = "Stability_later")
)


pls_model <- estimate_pls(
  data = all_metrics_population,
  measurement_model = measurement_model,
  structural_model = structural_model
)

summary(pls_model)
#plot_scores(pls_model)

# Bootstrap
boot_model <- bootstrap_model(
  seminr_model = pls_model,
  nboot = 10000,          
  cores = 6,
  seed = 123
)

boot_summary <- summary(boot_model, alpha = 0.05)

boot_summary$bootstrapped_paths

