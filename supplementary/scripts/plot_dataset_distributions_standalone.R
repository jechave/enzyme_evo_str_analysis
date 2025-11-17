#!/usr/bin/env Rscript
# Dataset distributions figure

library(tidyverse)
library(here)
library(cowplot)

cat("=== Generating Dataset Distributions Figure ===\n\n")

# Load dataset
cat("Loading dataset...\n")
dataset <- read_csv(here("data", "final_dataset_ca_ref.csv"),
                    show_col_types = FALSE) %>%
  mutate(mcsa_id = as.character(mcsa_id))

# Source the plotting function
source(here("supplementary/scripts", "plot_dataset_distributions.R"))

cat("Creating plot...\n")
plot <- plot_dataset_distributions(dataset)

# Save
output_file <- here("supplementary/figures",
                    "dataset_distributions.pdf")
ggsave(output_file, plot, width = 6.5, height = 9)

cat("\n=== Done ===\n")
cat("Saved:", output_file, "\n")
