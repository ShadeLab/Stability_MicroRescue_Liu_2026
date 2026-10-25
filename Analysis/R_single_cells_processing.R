# R script for STB Project ####
# Single-cell data processing ####

library(readr)
library(tidyverse)
library(ggpubr)

#********************************************####
# Image data ####

## Load data ####

# Load csv
if (!file.exists("logos/logos_syb_rsg_raw.rds")) {
  data_logos_raw <- list.files(
    "logos/data/",
    pattern = "_data\\.csv$", full.names = TRUE, recursive = TRUE
  ) |>
    map_df(~ read_csv(.))
  write_rds(data_logos_raw, "logos/logos_syb_rsg_raw.rds")
} else {
  data_logos_raw <- read_rds("logos/logos_syb_rsg_raw.rds")
}


# Join metadata
if (!file.exists("logos/logos_syb_rsg.rds")) {
  logos_syb_rsg <- read_rds("logos/logos_syb_rsg_raw.rds") |>
    inner_join(read_csv("logos/STB_metadata_logos.csv"),
               by = "experiment" )
  write_rds(logos_syb_rsg, "logos/logos_syb_rsg.rds")
} else {
  logos_syb_rsg <- read_rds("logos/logos_syb_rsg.rds")
}


## Filtering ####

logos_filtered <- logos_syb_rsg |>
  filter(
    flag_ring == FALSE,
    euler_number == 1,
    mean_distance_5_neighbors > 14,
    axis_major_length > 1,
    axis_major_length < 10,
    axis_minor_length > 0.5,
    axis_minor_length < 5,
    equivalent_diameter_area > 0.5,
    equivalent_diameter_area < 10,
    area < 15,
    intensity_cv > 0.1,
    intensity_cv < 0.3,
    solidity > 0.75,
    intensity_total > 0
  )

saveRDS(logos_filtered, "logos/logos_syb_rsg_filtered.rds")



#*********************####
# Community cell abundance ####
# Calculate cell abundance of communities.

logos_filtered <- readRDS("logos/logos_syb_rsg_filtered.rds")

concentration_factor <- 14342

df_cell_count <- logos_filtered %>%
  filter(Type %in% c("Total", "Active")) %>%
  group_by(Sample_id, Exp_id, Day, Replicate, Microcosm, Dye, Site_id, Season, 
           Moisture, Dilution, Original_temp, Group, Type) %>%
  summarise(counts = n()) %>%
  mutate(density = counts * concentration_factor * Dilution / (1 - Moisture)) %>%
  ungroup()%>%
  mutate(Sample = as.factor(paste(Exp_id, Day, Replicate, sep = "_")))

saveRDS(df_cell_count, "logos/df_cell_count.rds")


#*********************####
# Community redox activity ####

# calculate median_fluorescence, cumulative and mean total_fluorescence
df_redox <- logos_filtered %>%
  filter(Type %in% c("Active")) %>%
  group_by(Exp_id, Day, Replicate, Microcosm, Dye, Site_id, Season, 
           Moisture, Dilution, Original_temp, Group, Type) %>%
  summarise(counts = n(),
            intensity_sum = sum(intensity_total),
            area_sum =  sum(area),
            mean_fluorescence = intensity_sum/counts,
            median_fluorescence = median(intensity_total, na.rm = TRUE),
            intensity_per_area = intensity_sum / area_sum,
            .groups = "drop") %>%
  mutate(cumulative_fluo = mean_fluorescence * counts * concentration_factor * Dilution / (1 - Moisture)) %>%
  mutate(density = counts * concentration_factor * Dilution / (1 - Moisture))%>%
  mutate(Soil_id = gsub("(S)[CW]", "\\1", Group))%>%
  mutate(Site_type = case_when(Soil_id == "S1" ~ "Agriculture",Soil_id == "S2" ~ "Grassland", Soil_id == "S3" ~ "Forest")) %>%
  mutate(Sample = as.factor(paste(Exp_id, Day, Replicate, sep = "_")))

syb_counts <- df_cell_count %>%
  filter(Type == "Total") %>%
  dplyr::select(Sample, syb_counts = counts, syb_density = density)

df_redox <- df_redox %>%
  filter(Type %in% c("Active"))%>%
  left_join(syb_counts, by = c("Sample")) %>%
  mutate(mean_cumulative_fluo = cumulative_fluo / syb_density)

saveRDS(df_redox, "logos/df_redox.rds")

