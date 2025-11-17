#!/usr/bin/env Rscript
# Generate all main manuscript figures for ca_ref profile type

library(tidyverse)
library(cowplot)
library(here)

theme_set(theme_minimal())
theme_set(theme_cowplot(font_size = 10))

# Parameters
profile_type <- "ca_ref"
mcsa_id_case <- 2

cat("=== Generating figures for profile type:", profile_type, "===\n\n")

# Load data
cat("Loading data...\n")
dataset <- read_csv(here("data", paste0("final_dataset_", profile_type, ".csv")),
                    show_col_types = FALSE) %>%
  mutate(mcsa_id = as.character(mcsa_id))

profiles <- read_csv(here("data", paste0("final_dataset_profiles_", profile_type, ".csv")),
                     show_col_types = FALSE) %>%
  mutate(mcsa_id = as.character(mcsa_id))

data_gof <- read_csv(here("data", paste0("final_dataset_gof_", profile_type, ".csv")),
                     show_col_types = FALSE) %>%
  mutate(mcsa_id = as.character(mcsa_id))

dataset_and_gof <- dataset %>%
  inner_join(data_gof, by = "mcsa_id")

cat("Data loaded successfully.\n\n")

# Generate figures
cat("Generating figures...\n")

# Figure 1: Observed profiles and patterns
cat("  - fig_obs_", profile_type, ".pdf\n", sep = "")
source(here("scripts", "plot_obs.R"))
plot <- plot_obs(profiles, data_gof, mcsa_id_case)
ggsave(here("figures", "fig_obs.pdf"),
       plot, width = 6.5, height = 4.1)

# Figure 2: Models vs. observations
cat("  - fig_fit_vs_obs_", profile_type, ".pdf\n", sep = "")
source(here("scripts", "plot_fits_vs_obs.R"))
plot <- plot_fits_vs_obs(profiles, data_gof, mcsa_id_case)
ggsave(here("figures", "fig_fit_vs_obs.pdf"),
       plot, width = 6.5, height = 6.8)


# Figure 3: Split s1 and s2 contributions
cat("  - fig_split_", profile_type, ".pdf\n", sep = "")
source(here("scripts", "plot_split.R"))
plot <- plot_split(profiles, data_gof, mcsa_id_case)
ggsave(here("figures", "fig_split.pdf"),
       plot, width = 6.5, height = 4.1)



# Figure 4: Shapley map and examples
cat("  - fig_shapley_map_", profile_type, ".pdf\n", sep = "")
source(here("scripts", "plot_shapley_map.R"))
plot <- plot_shapley_map(
  data_gof = data_gof,
  dataset = dataset,
  mcsa_id_examples = c(858, 15, 2, 252, 908),
  mcsa_id_outliers = c(749, 258)
)
ggsave(here("figures", "fig_shapley_map.pdf"),
       plot, width = 5.5, height = 4.5)


# Figure 5: Active site divergence
cat("  - fig_asite_", profile_type, ".pdf\n", sep = "")
source(here("scripts", "plot_asite.R"))
plot <- plot_asite(profiles, data_gof)
ggsave(here("figures", "fig_asite.pdf"),
       plot, width = 3.5, height = 7)

cat("\n=== All figures generated successfully ===\n")
cat("Output location: figures/fig_*_", profile_type, ".pdf\n", sep = "")
