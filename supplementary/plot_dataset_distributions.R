plot_dataset_distributions <- function(dataset) {
  # Check and load required libraries
  check_and_load_libraries()

  # Define variables to plot and their labels
  numeric_vars <- c("nprot", "nresidues", "seqid_mean", "rmsd_mean", "nresidues_asite", "rmsd_mean_0")
  numeric_labels <- c("Number of proteins", "Reference protein size", "Id%", "RMSD", "Active site size", "Active site RMSD")

  categorical_vars <- c("cath_class", "cath_architecture", "ec_class")
  categorical_labels <- c("CATH class", "CATH architecture", "EC class")

  # Create plots with adjusted font sizes and labels
  numeric_plots <- map2(numeric_vars, numeric_labels,
                        ~create_histogram(dataset, .x, .y, font_size = 8, annotation_size = 2.5))
  categorical_plots <- map2(categorical_vars, categorical_labels,
                            ~create_barplot(dataset, .x, .y, font_size = 8))

  # Combine bar plots into a single column with increased spacing
  combined_barplot <- plot_grid(plotlist = categorical_plots, ncol = 1, align = 'v', axis = 'lr', rel_heights = c(1, 1.2, 1))

  # Combine histograms into a grid with appropriate number of rows and columns
  n_cols <- ceiling(length(numeric_plots) / 2)
  combined_histogram <- plot_grid(plotlist = numeric_plots, ncol = n_cols, rel_widths = rep(1, n_cols))

  # Combine all plots into a single figure with more space between barplots and histograms
  final_plot <- plot_grid(
    combined_barplot,
    combined_histogram,
    ncol = 1,
    rel_heights = c(1, 1.5)  # Adjusted the ratio for visual balance
  )

  # Add margins to the final combined plot
  final_plot_with_margins <- final_plot +
    theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"))

  return(final_plot_with_margins)
}

# Helper function to create a histogram for numeric variables with annotations
create_histogram <- function(data, variable, label, font_size = 8, annotation_size = 2) {
  stats <- data %>%
    summarise(
      mean = mean(.data[[variable]], na.rm = TRUE),
      sd = sd(.data[[variable]], na.rm = TRUE),
      min = min(.data[[variable]], na.rm = TRUE),
      max = max(.data[[variable]], na.rm = TRUE)
    )
  annotation <- sprintf(
    "Range: %.2f - %.2f\nMean: %.2f\nSD: %.2f",
    stats$min, stats$max, stats$mean, stats$sd
  )
  p <- ggplot(data, aes(x = .data[[variable]])) +
    geom_histogram(fill = "grey", color = "black") +
    labs(x = label, y = "Count") +
    theme_cowplot(font_size = font_size) +
    theme(plot.margin = margin(5, 5, 5, 5))
  y_range <- ggplot_build(p)$layout$panel_params[[1]]$y.range
  y_pos <- y_range[2] * 0.9
  p + annotate("text", x = Inf, y = y_pos, label = annotation,
               hjust = 1, vjust = 1, size = annotation_size)
}
# Helper function to create a bar plot for categorical variables
create_barplot <- function(data, variable, label, font_size = 8) {
  data %>%
    count(.data[[variable]]) %>%
    ggplot(aes(x = reorder(.data[[variable]], n), y = n)) +
    geom_bar(stat = "identity", fill = "grey") +
    coord_flip() +
    labs(x = label, y = "Count") +
    theme_cowplot(font_size = font_size) +
    theme(plot.margin = margin(5, 5, 5, 5))
}

# Function to check and load required libraries
check_and_load_libraries <- function() {
  required_packages <- c("tidyverse", "cowplot", "ggplot2")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop(paste("The following required packages are not installed:",
               paste(missing_packages, collapse = ", "),
               "\nPlease install these packages before running this function."))
  }

  invisible(sapply(required_packages, library, character.only = TRUE))
}
