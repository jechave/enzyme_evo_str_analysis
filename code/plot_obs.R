library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# Constants for styling
ANNOTATION_SIZE <- 2.5

# Main function
plot_obs <- function(profiles, data_gof, mcsa_id) {
  all_stats <- calculate_all_stats(profiles)
  plot_data <- prepare_plot_data(profiles, data_gof, mcsa_id)
  active_site_data <- filter(plot_data, dactive == 0)
  correlations <- calculate_correlations(plot_data)
  rmsd_stats <- calculate_rmsd_stats(plot_data)
  plots <- list(
    create_nlrmsd_plot(plot_data, active_site_data, rmsd_stats, title = "Divergence profile"),
    create_nlrmsd_vs_nlrmsf_plot(plot_data, active_site_data, correlations$nlrmsf, title = "Flexibility trend"),
    create_nlrmsd_vs_dactive_plot(plot_data, active_site_data, correlations$dactive, title = "Distance trend"),
    create_histogram(all_stats, "rmsd_ratio", "RMSD ratio (Q90/Q10)", title = "Divergence variation"),
    create_histogram(all_stats, "cor_lrmsd_lrmsf", expression(rho(nlRMSD, nlRMSF)), title = "Flexibility correlation"),
    create_histogram(all_stats, "cor_lrmsd_dactive", expression(rho(nlRMSD, d)), title = "Distance correlation")
  )
  create_final_plot(plots, mcsa_id)
}

# Helper functions
calculate_all_stats <- function(profiles) {
  profiles %>%
    group_by(mcsa_id) %>%
    summarise(
      rmsd_q90 = quantile(exp(lrmsd), 0.9),
      rmsd_q10 = quantile(exp(lrmsd), 0.1),
      rmsd_ratio = rmsd_q90 / rmsd_q10,
      cor_lrmsd_lrmsf = cor(lrmsd, lrmsf, method = "spearman"),
      cor_lrmsd_dactive = cor(lrmsd, dactive, method = "spearman")
    )
}

prepare_plot_data <- function(profiles, data_gof, mcsa_id) {
  profiles %>%
    filter(mcsa_id == !!mcsa_id) %>%
    mutate(
      nlrmsd = lrmsd - mean(lrmsd),
      nlrmsf = lrmsf - mean(lrmsf)
    ) %>%
    left_join(filter(data_gof, mcsa_id == !!mcsa_id))
}

calculate_correlations <- function(plot_data) {
  list(
    nlrmsf = cor.test(plot_data$nlrmsd, plot_data$nlrmsf, method = "spearman"),
    dactive = cor.test(plot_data$nlrmsd, plot_data$dactive, method = "spearman")
  )
}

calculate_rmsd_stats <- function(plot_data) {
  rmsd_q90 <- quantile(exp(plot_data$lrmsd), 0.9)
  rmsd_q10 <- quantile(exp(plot_data$lrmsd), 0.1)
  list(
    q90 = rmsd_q90,
    q10 = rmsd_q10,
    ratio = rmsd_q90 / rmsd_q10
  )
}

format_pvalue <- function(p) {
  if (p < 0.001) return("p < 0.001")
  sprintf("p = %.3f", p)
}

create_nlrmsd_plot <- function(plot_data, active_site_data, rmsd_stats, title = NULL) {
  y_limits <- range(plot_data$nlrmsd)
  y_limits[2] <- y_limits[2] * 1.2

  label_text <- sprintf("RMSD ratio (Q90/Q10) = %.1f / %.1f = %.1f",
                        rmsd_stats$q90, rmsd_stats$q10, rmsd_stats$ratio)

  p <- ggplot(plot_data) +
    geom_line(aes(x = ref_site, y = nlrmsd), color = "black", linewidth = 0.6) +
    geom_vline(data = active_site_data, aes(xintercept = ref_site),
               color = "red", linetype = "dashed", alpha = 0.3) +
    annotate("label",
             x = -Inf, y = Inf,
             label = label_text,
             hjust = -0.1, vjust = 1.1, size = ANNOTATION_SIZE,
             label.padding = unit(0.15, "lines"),
             label.r = unit(0, "lines"),
             fill = "white",
             color = "black",
             label.size = 0.2) +
    labs(x = "Residue", y = "nlRMSD") +
    coord_cartesian(ylim = y_limits) +
    theme_cowplot(font_size = 8)

  if (!is.null(title)) {
    p <- p + ggtitle(title) +
      theme( plot.title = element_text(face = "plain", color = "#4a4a4a"))
  }
  p
}
create_nlrmsd_vs_nlrmsf_plot <- function(plot_data, active_site_data, correlation, title = NULL) {
  rho_value <- sprintf("%.2f", correlation$estimate)
  p_value <- format_pvalue(correlation$p.value)

  rho_label <- bquote(rho == .(rho_value) * ", " * .(p_value))
  p <- ggplot(plot_data) +
    geom_point(aes(x = nlrmsf, y = nlrmsd), color = "black", alpha = 0.1, size = 0.5) +
    geom_point(data = active_site_data, aes(x = nlrmsf, y = nlrmsd),
               color = "red", size = 2, shape = 21, fill = "white") +
    geom_smooth(aes(x = nlrmsf, y = nlrmsd), method = "loess", se = TRUE,
                color = "black", linetype = "dotted", linewidth = 0.8) +
    geom_label(
      aes(x = -Inf, y = Inf),
      label = as.expression(rho_label),
      hjust = -0.1, vjust = 1.1, size = ANNOTATION_SIZE,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0, "lines"),
      fill = "white",
      color = "black",
      label.size = 0.2,
      inherit.aes = FALSE
    ) +
    labs(x = "nlRMSF", y = "nlRMSD") +
    theme_cowplot(font_size = 8)

  if (!is.null(title)) {
    p <- p + ggtitle(title) +
      theme( plot.title = element_text(face = "plain", color = "#4a4a4a"))
  }
  p
}

create_nlrmsd_vs_dactive_plot <- function(plot_data, active_site_data, correlation, title = NULL) {
  rho_value <- sprintf("%.2f", correlation$estimate)
  p_value <- format_pvalue(correlation$p.value)

  rho_label <- bquote(rho == .(rho_value) * ", " * .(p_value))
  p <- ggplot(plot_data) +
    geom_point(aes(x = dactive, y = nlrmsd), color = "black", alpha = 0.1, size = 0.5) +
    geom_point(data = active_site_data, aes(x = dactive, y = nlrmsd),
               color = "red", size = 2, shape = 21, fill = "white") +
    geom_smooth(aes(x = dactive, y = nlrmsd), method = "loess", se = TRUE,
                color = "black", linetype = "dotted", linewidth = 0.8) +
    geom_label(
      aes(x = -Inf, y = Inf),
      label = as.expression(rho_label),
      hjust = -0.1, vjust = 1.1, size = ANNOTATION_SIZE,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0, "lines"),
      fill = "white",
      color = "black",
      label.size = 0.2,
      inherit.aes = FALSE
    ) +
    labs(x = "d", y = "nlRMSD") +
    theme_cowplot(font_size = 8)

  if (!is.null(title)) {
    p <- p + ggtitle(title) +
      theme( plot.title = element_text(face = "plain", color = "#4a4a4a"))
  }
  p
}

create_histogram <- function(data, x_var, x_lab, title = NULL) {
  mean_value <- mean(data[[x_var]], na.rm = TRUE)
  sd_value <- sd(data[[x_var]], na.rm = TRUE)

  stat_label <- bquote(mu == .(sprintf("%.2f", mean_value)) %+-% .(sprintf("%.2f", sd_value)))

  p <- ggplot(data, aes(x = !!sym(x_var))) +
    geom_histogram(bins = 10, fill = "grey", color = "black") +
    labs(x = x_lab, y = "Count") +
    geom_label(
      aes(x = -Inf, y = Inf),
      label = as.expression(stat_label),
      hjust = -0.1, vjust = 1.1, size = ANNOTATION_SIZE,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0, "lines"),
      fill = "white",
      color = "black",
      label.size = 0.2,
      inherit.aes = FALSE
    ) +
    theme_cowplot(font_size = 8)

  if (!is.null(title)) {
    p <- p + ggtitle(title)  +
      theme( plot.title = element_text(face = "plain", size = 9, color = "#4a4a4a"))
  }

  # hist_data <- ggplot_build(p)$data[[1]]
  # y_max <- max(hist_data$count) * 1.1
  # p + coord_cartesian(ylim = c(0, y_max))

  p +  scale_y_continuous(expand = expansion(mult = c(0, 0.35)))
}

combine_plots <- function(plots, row_index) {
  plot_grid(
    plots[[1 + (row_index - 1) * 3]],
    plots[[2 + (row_index - 1) * 3]],
    plots[[3 + (row_index - 1) * 3]],
    ncol = 3,
    rel_widths = c(1.5, 1, 1),
    labels = if(row_index == 1) c("(a)", "(b)", "(c)") else c("(d)", "(e)", "(f)"),
    label_size = 9
  )
}

create_final_plot <- function(plots, mcsa_id) {

  # Create main title
  main_title <- ggdraw() +
    draw_label(
      paste0("Observed structural divergence patterns"),
      fontface = "plain",
      size = 12,
      color = "#4a4a4a"
    )

  # Create the first row title
  first_row_title <- ggdraw() +
    draw_label(
      paste0("Example family: M-CSA ID = ", mcsa_id),
      fontface = "plain",
      size = 10,
      color = "#4a4a4a"
    )

  # Create the second row title
  second_row_title <- ggdraw() +
    draw_label(
      "Aggregate results across all families",
      fontface = "plain",
      size = 10,
      color = "#4a4a4a"
    )

  # Combine all elements using a single plot_grid
  final_plot <- plot_grid(
    main_title,
    first_row_title,
    combine_plots(plots, 1),
    ggdraw(),  # empty space for separation
    second_row_title,
    combine_plots(plots, 2),
    ncol = 1,
    rel_heights = c(0.3, 0.3, 1.5, 0.2, 0.3, 1.5),
    align = 'v'
  )

  final_plot <- final_plot +
    theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"))

  return(final_plot)
}
