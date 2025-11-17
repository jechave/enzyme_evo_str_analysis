#!/usr/bin/env Rscript
# Generate all case-by-case figures

library(tidyverse)
library(here)
library(cowplot)

cat("=== Generating All Case-by-Case Figures ===\n\n")

# Load data
cat("Loading data...\n")
dataset <- read_csv(here("data", "final_dataset_ca_ref.csv"),
                    show_col_types = FALSE) %>%
  mutate(mcsa_id = as.character(mcsa_id))

profiles <- read_csv(here("data", "final_dataset_profiles_ca_ref.csv"),
                     show_col_types = FALSE) %>%
  mutate(mcsa_id = as.character(mcsa_id))

data_gof <- read_csv(here("data", "final_dataset_gof_ca_ref.csv"),
                     show_col_types = FALSE) %>%
  mutate(mcsa_id = as.character(mcsa_id))

# Source plotting functions
cat("Loading plotting functions...\n")
source(here("supplementary/scripts", "plot_case.R"))

# Get list of all mcsa_ids
d <- data_gof %>%
  arrange(as.numeric(mcsa_id))
mcsa_ids <- d$mcsa_id

cat(sprintf("Generating %d case figures...\n", length(mcsa_ids)))

# Generate each figure
for (id in mcsa_ids) {
  cat(sprintf("  Processing MCSA ID %s...\n", id))

  # Generate the plot
  plot <- plot_case(profiles, data_gof, id)

  # Save with zero-padded ID for proper sorting
  output_file <- here("supplementary/figures",
                      sprintf("case_%s.pdf", id))
  ggsave(output_file, plot, width = 6.5, height = 9)
}

cat("\n=== Done ===\n")
cat(sprintf("Generated %d case figures in figures/\n", length(mcsa_ids)))
