#!/usr/bin/env Rscript
# Generate robustness table comparing GOF metrics across methodological choices
# Output: LaTeX table in manuscript/tables/table_robustness.tex

library(tidyverse)
library(here)

cat("=== Generating Robustness Table ===\n\n")

# Load GOF data for all three methods
cat("Loading GOF data...\n")
gof_ca_ref <- read_csv(here("data", "final_dataset_gof_ca_ref.csv"),
                       show_col_types = FALSE) %>%
  mutate(mcsa_id = as.character(mcsa_id),
         method = "ca_ref")

gof_ca_hom <- read_csv(here("data", "final_dataset_gof_ca_homolog.csv"),
                       show_col_types = FALSE) %>%
  mutate(mcsa_id = as.character(mcsa_id),
         method = "ca_hom")

gof_cb_ref <- read_csv(here("data", "final_dataset_gof_cb_ref.csv"),
                       show_col_types = FALSE) %>%
  mutate(mcsa_id = as.character(mcsa_id),
         method = "CB")

# Combine all data
gof_all <- bind_rows(gof_ca_ref, gof_ca_hom, gof_cb_ref)

# Calculate RSC
gof_all <- gof_all %>%
  mutate(rsc_dist = shapley_dactive / (shapley_lrmsf + shapley_dactive))

cat("Calculating summary statistics...\n")

# Define metrics and their labels
metrics <- tribble(
  ~metric_name, ~metric_label, ~metric_type,
  "dev.expl.fit1", "M\\textsubscript{1}", "Dev.Expl.",
  "dev.expl.fit2", "M\\textsubscript{2}", "Dev.Expl.",
  "dev.expl.fit12", "M\\textsubscript{12}", "Dev.Expl.",
  "rmse_lrmsf.fit1", "M\\textsubscript{1}", "RMSE(flex)",
  "rmse_lrmsf.fit2", "M\\textsubscript{2}", "RMSE(flex)",
  "rmse_lrmsf.fit12", "M\\textsubscript{12}", "RMSE(flex)",
  "rmse_dactive.fit1", "M\\textsubscript{1}", "RMSE(dist)",
  "rmse_dactive.fit2", "M\\textsubscript{2}", "RMSE(dist)",
  "rmse_dactive.fit12", "M\\textsubscript{12}", "RMSE(dist)",
  "shapley_lrmsf", "SC(flex)", "Shapley",
  "shapley_dactive", "SC(dist)", "Shapley",
  "rsc_dist", "RSC(dist)", "Shapley"
)

# Calculate mean ± SD for each metric and method
summary_list <- list()
for (i in 1:nrow(metrics)) {
  metric_col <- metrics$metric_name[i]

  stats <- gof_all %>%
    group_by(method) %>%
    summarise(
      mean_val = mean(.data[[metric_col]], na.rm = TRUE),
      sd_val = sd(.data[[metric_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(formatted = sprintf("%.2f$\\pm$%.2f", mean_val, sd_val)) %>%
    select(method, formatted) %>%
    pivot_wider(names_from = method, values_from = formatted)

  summary_list[[i]] <- bind_cols(
    metrics[i, c("metric_type", "metric_label")],
    stats
  )
}

summary_stats <- bind_rows(summary_list)

cat("Building LaTeX table...\n")

# Start building LaTeX table
latex_lines <- c(
  "\\begin{table}[htbp]",
  "\\small",
  "\\centering",
  "\\setlength{\\tabcolsep}{6pt}",
  "\\renewcommand{\\arraystretch}{1.2}",
  "\\begin{tabular}{llrrr}",
  "\\toprule",
  "Metric & Model & ca\\_ref & ca\\_hom & CB \\\\",
  "\\midrule"
)

# Add rows grouped by metric type
current_type <- ""
for (i in 1:nrow(summary_stats)) {
  row <- summary_stats[i, ]

  # Add line spacing between metric types
  if (row$metric_type != current_type && current_type != "") {
    latex_lines <- c(latex_lines, "\\addlinespace")
  }
  current_type <- row$metric_type

  # Build table row
  metric_col <- if (i == 1 || summary_stats$metric_type[i] != summary_stats$metric_type[i-1]) {
    paste0("\\textbf{", row$metric_type, "}")
  } else {
    ""
  }

  table_row <- sprintf("%s & %s & %s & %s & %s \\\\",
                      metric_col,
                      row$metric_label,
                      row$ca_ref,
                      row$ca_hom,
                      row$CB)
  latex_lines <- c(latex_lines, table_row)
}

# Complete the table
latex_lines <- c(latex_lines,
  "\\bottomrule",
  "\\end{tabular}",
  paste0("\\caption{Robustness of goodness-of-fit metrics to methodological choices. ",
         "Values show mean $\\pm$ standard deviation across 53 enzyme families. ",
         "ca\\_ref: reference protein with CA atoms; ",
         "ca\\_hom: closest homolog with CA atoms; ",
         "CB: reference protein with CB atoms. ",
         "High consistency across methods indicates robustness. ",
         "See Supplementary Figures for detailed scatter plots and correlations.}"),
  "\\label{tab:robustness}",
  "\\end{table}"
)

# Write to file
output_file <- here("tables", "table_robustness.tex")
writeLines(latex_lines, output_file)

cat("\n=== Done ===\n")
cat("Table saved to:", output_file, "\n")
cat(sprintf("Generated %d rows\n", nrow(summary_stats)))
