library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# Constants for styling
ANNOTATION_SIZE <- 2.5

# Main function
plot_case_obs <- function(profiles, data_gof, mcsa_id) {
  plot_data <- prepare_plot_data(profiles, data_gof, mcsa_id)
  active_site_data <- filter(plot_data, dactive == 0)
  correlations <- calculate_correlations(plot_data)
  rmsd_stats <- calculate_rmsd_stats(plot_data)

  plots <- list(
    create_nlrmsd_plot(plot_data, active_site_data, rmsd_stats, title = "Divergence profile"),
    create_nlrmsd_vs_nlrmsf_plot(plot_data, active_site_data, correlations$nlrmsf, title = "Flexibility trend"),
    create_nlrmsd_vs_dactive_plot(plot_data, active_site_data, correlations$dactive, title = "Distance trend")
  )

  create_final_plot(plots, mcsa_id)
}

# Helper functions
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
      theme(plot.title = element_text(face = "plain", color = "#4a4a4a"))
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
      theme(plot.title = element_text(face = "plain", color = "#4a4a4a"))
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
      theme(plot.title = element_text(face = "plain", color = "#4a4a4a"))
  }
  p
}

create_final_plot <- function(plots, mcsa_id) {
  # Combine plots
  final_plot <- plot_grid(
    plots[[1]], plots[[2]], plots[[3]],
    ncol = 3,
    rel_widths = c(1.5, 1, 1),
    labels = c("(a)", "(b)", "(c)"),
    label_size = 9,
    hjust = 0
  )

  return(final_plot)
}
