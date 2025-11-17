#!/usr/bin/env Rscript
# Active site constraint correlations figure

library(tidyverse)
library(here)
library(cowplot)
library(ggpubr)

cat("=== Generating Active Site Correlations Figure ===\n\n")

# Load data
cat("Loading data...\n")
profiles <- read_csv(here("data", "final_dataset_profiles_ca_ref.csv"),
                     show_col_types = FALSE) %>%
  mutate(mcsa_id = as.character(mcsa_id))

data_gof <- read_csv(here("data", "final_dataset_gof_ca_ref.csv"),
                     show_col_types = FALSE) %>%
  mutate(mcsa_id = as.character(mcsa_id))

cat("Calculating family-averaged active site metrics...\n")

# Calculate family-averaged active site metrics (same as in supplement tables)
family_values <- profiles %>%
  filter(mcsa_id != "749") %>%
  group_by(mcsa_id) %>%
  mutate(
    s1 = lrmsd.fit12.c1 - mean(lrmsd.fit12.c1, na.rm = TRUE),
    s2 = lrmsd.fit12.c2 - mean(lrmsd.fit12.c2, na.rm = TRUE),
    nlrmsf = lrmsf - mean(lrmsf, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(near(dactive, 0)) %>%
  mutate(
    abs_sum = abs(s1) + abs(s2),
    prop_s1 = abs(s1) / abs_sum,
    prop_s2 = abs(s2) / abs_sum
  ) %>%
  group_by(mcsa_id) %>%
  summarize(
    s1 = mean(s1, na.rm = TRUE),
    s2 = mean(s2, na.rm = TRUE),
    rs1 = mean(prop_s1, na.rm = TRUE),
    rs2 = mean(prop_s2, na.rm = TRUE),
    nlrmsf = mean(nlrmsf, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(data_gof %>% select(mcsa_id, shapley_lrmsf, shapley_dactive), by = "mcsa_id")

cat("Creating plots...\n")

# Helper function to convert p-value to significance stars
p_to_stars <- function(p) {
  if (p < 0.0001) return("****")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

# Calculate correlations and p-values
cor_test_1 <- cor.test(family_values$nlrmsf, family_values$s1, method = "pearson")
cor_test_2 <- cor.test(family_values$shapley_dactive, family_values$s2, method = "pearson")
cor_test_3 <- cor.test(family_values$nlrmsf, family_values$rs1, method = "pearson")
cor_test_4 <- cor.test(family_values$shapley_dactive, family_values$rs1, method = "pearson")

# Create labels
label_1 <- sprintf("R = %.2f, %s", cor_test_1$estimate, p_to_stars(cor_test_1$p.value))
label_2 <- sprintf("R = %.2f, %s", cor_test_2$estimate, p_to_stars(cor_test_2$p.value))
label_3 <- sprintf("R = %.2f, %s", cor_test_3$estimate, p_to_stars(cor_test_3$p.value))
label_4 <- sprintf("R = %.2f, %s", cor_test_4$estimate, p_to_stars(cor_test_4$p.value))

# Calculate shared y-axis limits for panels a and b
y_min_ab <- min(c(family_values$s1, family_values$s2))
y_max_ab <- max(c(family_values$s1, family_values$s2))
y_range_ab <- y_max_ab - y_min_ab
y_limits_ab <- c(y_min_ab - 0.05 * y_range_ab, y_max_ab + 0.15 * y_range_ab)

# Panel (a): s1 vs nlRMSF
p1 <- ggplot(family_values, aes(nlrmsf, s1)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = -Inf, y = Inf, label = label_1,
           hjust = -0.1, vjust = 1.5, size = 4) +
  scale_y_continuous(limits = y_limits_ab) +
  xlab("Active site nlRMSF") +
  ylab(expression(paste("Active site ", s[1]))) +
  theme_cowplot() +
  theme(plot.margin = unit(c(0.5, 0.25, 0.25, 0.25), "cm"),
        panel.grid.major = element_line(color = "gray90"))

# Panel (b): s2 vs shapley_dactive (SC(s2))
p2 <- ggplot(family_values, aes(shapley_dactive, s2)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = -Inf, y = Inf, label = label_2,
           hjust = -0.1, vjust = 1.5, size = 4) +
  scale_y_continuous(limits = y_limits_ab) +
  xlab(expression(paste("Whole-protein SC(", s[2], ")"))) +
  ylab(expression(paste("Active site ", s[2]))) +
  theme_cowplot() +
  theme(plot.margin = unit(c(0.5, 0.25, 0.25, 0.25), "cm"),
        panel.grid.major = element_line(color = "gray90"))

# Panel (c): rs1 vs nlRMSF
p3 <- ggplot(family_values, aes(nlrmsf, rs1)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = -Inf, y = Inf, label = label_3,
           hjust = -0.1, vjust = 1.5, size = 4) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  xlab("Active site nlRMSF") +
  ylab(expression(paste("Active site ", rs[1]))) +
  theme_cowplot() +
  theme(plot.margin = unit(c(0.5, 0.25, 0.25, 0.25), "cm"),
        panel.grid.major = element_line(color = "gray90"))

# Panel (d): rs1 vs shapley_dactive
p4 <- ggplot(family_values, aes(shapley_dactive, rs1)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = -Inf, y = Inf, label = label_4,
           hjust = -0.1, vjust = 1.5, size = 4) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  xlab(expression(paste("Whole-protein SC(", s[2], ")"))) +
  ylab(expression(paste("Active site ", rs[1]))) +
  theme_cowplot() +
  theme(plot.margin = unit(c(0.5, 0.25, 0.25, 0.25), "cm"),
        panel.grid.major = element_line(color = "gray90"))

# Combine into 2x2 grid
final_plot <- plot_grid(p1, p2, p3, p4,
                        ncol = 2,
                        labels = c("(a)", "(b)", "(c)", "(d)"),
                        hjust = 0,
                        vjust = 1) +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"))

# Save
output_file <- here("supplementary/figures", "asite_correlations.pdf")
ggsave(output_file, final_plot, width = 8, height = 7)

cat("\n=== Done ===\n")
cat("Saved:", output_file, "\n")
