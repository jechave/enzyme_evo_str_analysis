#!/usr/bin/env Rscript
# Correlation between flexibility and distance figure

library(tidyverse)
library(here)
library(cowplot)
library(ggpubr)

cat("=== Generating Flex-Dist Correlation Figure ===\n\n")

# Load profiles data
cat("Loading profiles data...\n")
profiles <- read_csv(here("data", "final_dataset_profiles_ca_ref.csv"),
                     show_col_types = FALSE) %>%
  mutate(mcsa_id = as.character(mcsa_id))

cat("Calculating statistics...\n")

# First calculate the statistics
stats_summary <- profiles %>%
  group_by(mcsa_id) %>%
  summarise(rho_rmsf_dactive = cor(lrmsf, dactive, method = "spearman")) %>%
  summarise(
    mean_rho = mean(rho_rmsf_dactive, na.rm = TRUE),
    sd_rho = sd(rho_rmsf_dactive, na.rm = TRUE),
    min_rho = min(rho_rmsf_dactive, na.rm = TRUE),
    max_rho = max(rho_rmsf_dactive, na.rm = TRUE)
  )

cat("Creating plots...\n")

# Create the plot with annotations
p1 <- profiles %>%
  group_by(mcsa_id) %>%
  summarise(rho_rmsf_dactive = cor(lrmsf, dactive, method = "spearman")) %>%
  ggplot(aes(rho_rmsf_dactive)) +
  geom_histogram(fill = "grey", color = "black", bins = 10) +
  xlab(bquote(rho(nlRMSF, d))) +
  ylab("Number of cases") +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    hjust = -0.1,
    vjust = 1.5,
    label = sprintf(
      "Range: [%.2f, %.2f]\nMean: %.2f\nSD: %.2f",
      stats_summary$min_rho,
      stats_summary$max_rho,
      stats_summary$mean_rho,
      stats_summary$sd_rho
    )
  )  +
  theme_cowplot()

p2 <- profiles %>%
  group_by(mcsa_id) %>%
  summarise(rho_rmsf_dactive = cor(lrmsf, dactive, method = "spearman"),
            nlrmsf_asite = mean(nlrmsf[near(dactive, 0)])) %>%
  ggplot(aes(nlrmsf_asite, rho_rmsf_dactive)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor() +
  ylim(-0.1, 1) +
  xlab("Active site flexibility (nlRMSF)") +
  ylab(bquote(rho(nlRMSF, d)))  +
  theme_cowplot()

final_plot <- plot_grid(p1, p2, ncol = 1,
                        labels = c("(a)", "(b)"),
                        hjust = 0,
                        vjust = 0) +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"))

# Save
output_file <- here("supplementary/figures",
                    "flex_dist_correlation.pdf")
ggsave(output_file, final_plot, width = 6.5, height = 5.5)

cat("\n=== Done ===\n")
cat("Saved:", output_file, "\n")
