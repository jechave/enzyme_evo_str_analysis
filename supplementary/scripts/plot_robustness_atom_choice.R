#!/usr/bin/env Rscript
# Robustness figure: CA vs. CB atom choice
# 4 rows × 3 columns layout (12 subpanels total)

library(tidyverse)
library(here)
library(ggrepel)
library(cowplot)

# Color palette (consistent with supplementary material)
COLOR_PALETTE <- c(fit1 = "#35B779FF", fit2 = "#FDE725FF", fit12 = "#31688EFF")

cat("=== Generating Robustness Figure: Atom Choice ===\n\n")

# Load GOF data
cat("Loading GOF data...\n")
gof_files <- list.files(here("data"),
                        pattern = "^final_dataset_gof_.*\\.csv$",
                        full.names = TRUE)

gof_combined <- gof_files %>%
  set_names(basename) %>%
  map_dfr(~read_csv(.x, show_col_types = FALSE), .id = "file") %>%
  mutate(
    profile_type = str_extract(file, "(?<=final_dataset_gof_)[^.]+"),
    mcsa_id = as.character(mcsa_id)
  ) %>%
  filter(profile_type %in% c("ca_ref", "cb_ref"))

# Add RSC calculation
gof_combined <- gof_combined %>%
  mutate(rsc_dist = shapley_dactive / (shapley_lrmsf + shapley_dactive))

# Pivot to comparison format
comparison_data <- gof_combined %>%
  select(mcsa_id, profile_type,
         dev.expl.fit1, dev.expl.fit2, dev.expl.fit12,
         rmse_lrmsf.fit1, rmse_lrmsf.fit2, rmse_lrmsf.fit12,
         rmse_dactive.fit1, rmse_dactive.fit2, rmse_dactive.fit12,
         shapley_lrmsf, shapley_dactive, rsc_dist) %>%
  pivot_wider(
    names_from = profile_type,
    values_from = c(dev.expl.fit1, dev.expl.fit2, dev.expl.fit12,
                    rmse_lrmsf.fit1, rmse_lrmsf.fit2, rmse_lrmsf.fit12,
                    rmse_dactive.fit1, rmse_dactive.fit2, rmse_dactive.fit12,
                    shapley_lrmsf, shapley_dactive, rsc_dist)
  )

cat("Data prepared.\n\n")

# Helper function to create individual subpanel
make_subpanel <- function(data, x_col, y_col, title, lims = NULL, color = "black") {
  # Extract x and y values
  x_vals <- data[[x_col]]
  y_vals <- data[[y_col]]

  # Calculate statistics
  cor_val <- cor(x_vals, y_vals, use = "complete.obs")
  mean_x <- mean(x_vals, na.rm = TRUE)
  sd_x <- sd(x_vals, na.rm = TRUE)
  mean_y <- mean(y_vals, na.rm = TRUE)
  sd_y <- sd(y_vals, na.rm = TRUE)

  # Get limits (use provided if available, otherwise calculate)
  if (is.null(lims)) {
    all_vals <- c(x_vals, y_vals)
    lims <- c(min(all_vals, na.rm = TRUE), max(all_vals, na.rm = TRUE))
  }

  # Prepare plot data
  plot_data <- data %>%
    select(mcsa_id, x = all_of(x_col), y = all_of(y_col))

  # Create annotation text lines separately for better control
  ann_line1 <- sprintf("r = %.2f", cor_val)
  ann_line2 <- sprintf("CA = %.2f ± %.2f", mean_x, sd_x)
  ann_line3 <- sprintf("CB = %.2f ± %.2f", mean_y, sd_y)

  # Position for annotation (with some padding from edges)
  x_pos <- lims[1] + 0.02 * (lims[2] - lims[1])
  y_pos <- lims[2] - 0.02 * (lims[2] - lims[1])

  ggplot(plot_data, aes(x = x, y = y)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey50", linewidth = 0.3) +
    geom_point(size = 1.2, alpha = 0.75, color = color) +
    coord_fixed(xlim = lims, ylim = lims) +
    labs(x = "CA", y = "CB", title = title) +
    annotate("text", x = x_pos, y = y_pos, label = ann_line1,
             hjust = 0, vjust = 1, size = 2.8) +
    annotate("text", x = x_pos, y = y_pos - 0.06 * (lims[2] - lims[1]), label = ann_line2,
             hjust = 0, vjust = 1, size = 2.8) +
    annotate("text", x = x_pos, y = y_pos - 0.12 * (lims[2] - lims[1]), label = ann_line3,
             hjust = 0, vjust = 1, size = 2.8) +
    theme_cowplot(font_size = 10) +
    theme(
      plot.margin = unit(c(4, 2, 4, 2), "pt"),
      plot.title = element_text(face = "plain", size = 9)
    )
}

cat("Generating subpanels...\n")

# Panel (a): Deviance Explained - 3 subpanels horizontally with shared limits
lims_a <- range(c(comparison_data$dev.expl.fit1_ca_ref, comparison_data$dev.expl.fit1_cb_ref,
                  comparison_data$dev.expl.fit2_ca_ref, comparison_data$dev.expl.fit2_cb_ref,
                  comparison_data$dev.expl.fit12_ca_ref, comparison_data$dev.expl.fit12_cb_ref), na.rm = TRUE)
p_a1 <- make_subpanel(comparison_data, "dev.expl.fit1_ca_ref", "dev.expl.fit1_cb_ref", expression("Dev.Expl. " * M[1]), lims_a, COLOR_PALETTE["fit1"])
p_a2 <- make_subpanel(comparison_data, "dev.expl.fit2_ca_ref", "dev.expl.fit2_cb_ref", expression("Dev.Expl. " * M[2]), lims_a, COLOR_PALETTE["fit2"])
p_a3 <- make_subpanel(comparison_data, "dev.expl.fit12_ca_ref", "dev.expl.fit12_cb_ref", expression("Dev.Expl. " * M[12]), lims_a, COLOR_PALETTE["fit12"])
panel_a <- plot_grid(p_a1, p_a2, p_a3, ncol = 3)

# Panel (b): RMSE(flex) - 3 subpanels horizontally with shared limits (colored by model)
lims_b <- range(c(comparison_data$rmse_lrmsf.fit1_ca_ref, comparison_data$rmse_lrmsf.fit1_cb_ref,
                  comparison_data$rmse_lrmsf.fit2_ca_ref, comparison_data$rmse_lrmsf.fit2_cb_ref,
                  comparison_data$rmse_lrmsf.fit12_ca_ref, comparison_data$rmse_lrmsf.fit12_cb_ref), na.rm = TRUE)
p_b1 <- make_subpanel(comparison_data, "rmse_lrmsf.fit1_ca_ref", "rmse_lrmsf.fit1_cb_ref", expression("RMSE(flex) " * M[1]), lims_b, COLOR_PALETTE["fit1"])
p_b2 <- make_subpanel(comparison_data, "rmse_lrmsf.fit2_ca_ref", "rmse_lrmsf.fit2_cb_ref", expression("RMSE(flex) " * M[2]), lims_b, COLOR_PALETTE["fit2"])
p_b3 <- make_subpanel(comparison_data, "rmse_lrmsf.fit12_ca_ref", "rmse_lrmsf.fit12_cb_ref", expression("RMSE(flex) " * M[12]), lims_b, COLOR_PALETTE["fit12"])
panel_b <- plot_grid(p_b1, p_b2, p_b3, ncol = 3)

# Panel (c): RMSE(dist) - 3 subpanels horizontally with shared limits (colored by model)
lims_c <- range(c(comparison_data$rmse_dactive.fit1_ca_ref, comparison_data$rmse_dactive.fit1_cb_ref,
                  comparison_data$rmse_dactive.fit2_ca_ref, comparison_data$rmse_dactive.fit2_cb_ref,
                  comparison_data$rmse_dactive.fit12_ca_ref, comparison_data$rmse_dactive.fit12_cb_ref), na.rm = TRUE)
p_c1 <- make_subpanel(comparison_data, "rmse_dactive.fit1_ca_ref", "rmse_dactive.fit1_cb_ref", expression("RMSE(dist) " * M[1]), lims_c, COLOR_PALETTE["fit1"])
p_c2 <- make_subpanel(comparison_data, "rmse_dactive.fit2_ca_ref", "rmse_dactive.fit2_cb_ref", expression("RMSE(dist) " * M[2]), lims_c, COLOR_PALETTE["fit2"])
p_c3 <- make_subpanel(comparison_data, "rmse_dactive.fit12_ca_ref", "rmse_dactive.fit12_cb_ref", expression("RMSE(dist) " * M[12]), lims_c, COLOR_PALETTE["fit12"])
panel_c <- plot_grid(p_c1, p_c2, p_c3, ncol = 3)

# Panel (d): Shapley values - first 2 subpanels share limits, RSC has its own
# SC(flex)=green, SC(dist) and RSC(dist)=yellow
lims_d_shapley <- range(c(comparison_data$shapley_lrmsf_ca_ref, comparison_data$shapley_lrmsf_cb_ref,
                          comparison_data$shapley_dactive_ca_ref, comparison_data$shapley_dactive_cb_ref), na.rm = TRUE)
lims_d_rsc <- range(c(comparison_data$rsc_dist_ca_ref, comparison_data$rsc_dist_cb_ref), na.rm = TRUE)
p_d1 <- make_subpanel(comparison_data, "shapley_lrmsf_ca_ref", "shapley_lrmsf_cb_ref", "SC(flex)", lims_d_shapley, COLOR_PALETTE["fit1"])
p_d2 <- make_subpanel(comparison_data, "shapley_dactive_ca_ref", "shapley_dactive_cb_ref", "SC(dist)", lims_d_shapley, COLOR_PALETTE["fit2"])
p_d3 <- make_subpanel(comparison_data, "rsc_dist_ca_ref", "rsc_dist_cb_ref", "RSC(dist)", lims_d_rsc, COLOR_PALETTE["fit2"])
panel_d <- plot_grid(p_d1, p_d2, p_d3, ncol = 3)

cat("Combining panels...\n")

# Combine all 4 panels vertically with spacing between rows
final_plot <- plot_grid(panel_a, panel_b, panel_c, panel_d,
                        ncol = 1,
                        labels = c("(a)", "(b)", "(c)", "(d)"),
                        label_size = 10,
                        rel_heights = c(1, 1, 1, 1),
                        align = "v")

# Save
output_file <- here("supplementary/figures",
                    "robustness_atom_choice.pdf")
ggsave(output_file, final_plot, width = 7.5, height = 10)

cat("\n=== Done ===\n")
cat("Saved:", output_file, "\n")
cat("4 rows × 3 columns (12 subpanels total)\n")
