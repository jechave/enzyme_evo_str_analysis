# Constants
ANNOTATION_SIZE <- 2.5
CUSTOM_PALETTE <- c(Observed = "#000000", M12 = "#31688EFF", s1 = "#35B779FF", s2 = "#FDE725FF")

# Main function to create the plot
plot_case_split <- function(profiles, data_gof, mcsa_id) {
  prepared_data <- prepare_data(profiles, data_gof, mcsa_id)
  plots <- create_plots(prepared_data)
  create_layout(plots, mcsa_id)
}

# Helper function to prepare data
prepare_data <- function(profiles, data_gof, mcsa_id) {
  gof_data <- data_gof %>% filter(mcsa_id == !!mcsa_id)
  plot_data <- profiles %>%
    filter(mcsa_id == !!mcsa_id) %>%
    mutate(
      Observed = lrmsd - lrmsd.fit12.c0,
      M12 = lrmsd.fit12 - lrmsd.fit12.c0,
      s1 = lrmsd.fit12.c1,
      s2 = lrmsd.fit12.c2,
      nlrmsf = lrmsf - mean(lrmsf),
      nlrmsd = lrmsd - mean(lrmsd),
      nlrmsd.fit12 = lrmsd.fit12 - mean(lrmsd.fit12),
      nlrmsd_s1 = lrmsd.fit12.c1 - mean(lrmsd.fit12.c1),
      nlrmsd_s2 = lrmsd.fit12.c2 - mean(lrmsd.fit12.c2)
    ) %>%
    left_join(gof_data, by = "mcsa_id")
  list(plot_data = plot_data, gof_data = gof_data)
}

# Helper functions for creating individual plots
create_plots <- function(prepared_data) {
  plot_data <- prepared_data$plot_data
  gof_data <- prepared_data$gof_data
  p1 <- create_panel1(plot_data, gof_data)
  p2 <- create_panel2(plot_data)
  p3 <- create_panel3(plot_data)
  list(p1 = p1, p2 = p2, p3 = p3)
}

create_panel1 <- function(plot_data, gof_data) {
  active_site_data <- plot_data %>% filter(dactive == 0)
  y_limits <- range(plot_data$Observed, plot_data$M12, plot_data$s1, plot_data$s2, na.rm = TRUE)
  y_limits[2] <- y_limits[2] * 1.2
  ggplot(plot_data) +
    geom_vline(data = active_site_data, aes(xintercept = ref_site),
               color = "red", linetype = "dashed", alpha = 0.3) +
    geom_line(aes(x = ref_site, y = M12, color = "M12"), linewidth = 0.6) +
    geom_line(aes(x = ref_site, y = s1, color = "s1"), linewidth = 0.6) +
    geom_line(aes(x = ref_site, y = s2, color = "s2"), linewidth = 0.6) +
    geom_label(
      aes(x = -Inf, y = Inf),
      label = sprintf("Expl. Dev: %.2f, s1: %.2f, s2: %.2f", gof_data$dev.expl.fit12, gof_data$shapley_lrmsf, gof_data$shapley_dactive),
      hjust = -0.1, vjust = 1.1, size = ANNOTATION_SIZE,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0, "lines"),
      fill = "white",
      color = "black",
      label.size = 0.2,
      inherit.aes = FALSE
    ) +
    labs(x = "Residue", y = "shifted lRMSD", title = "Profile") +
    coord_cartesian(ylim = y_limits) +
    theme_cowplot(font_size = 8) +
    scale_color_manual(values = CUSTOM_PALETTE) +
    theme(legend.position = "none", plot.title = element_text(face = "plain", color = "#4a4a4a"))
}

create_panel2 <- function(plot_data) {
  ggplot(plot_data) +
    geom_smooth(aes(x = nlrmsf, y = nlrmsd, color = "Observed"), method = "loess", se = FALSE, linetype = "dotted", linewidth = 0.8) +
    geom_smooth(aes(x = nlrmsf, y = nlrmsd.fit12, color = "M12"), method = "loess", se = FALSE, linewidth = 0.6) +
    geom_smooth(aes(x = nlrmsf, y = nlrmsd_s1, color = "s1"), method = "loess", se = FALSE, linewidth = 0.6) +
    geom_smooth(aes(x = nlrmsf, y = nlrmsd_s2, color = "s2"), method = "loess", se = FALSE, linewidth = 0.6) +
    labs(x = "nlRMSF", y = "nlRMSD", title = "Flexibility trend") +
    theme_cowplot(font_size = 8) +
    scale_color_manual(values = CUSTOM_PALETTE) +
    theme(legend.position = "none", plot.title = element_text(face = "plain", color = "#4a4a4a"))
}

create_panel3 <- function(plot_data) {
  ggplot(plot_data) +
    geom_smooth(aes(x = dactive, y = nlrmsd, color = "Observed"), method = "loess", se = FALSE, linetype = "dotted", linewidth = 0.8) +
    geom_smooth(aes(x = dactive, y = nlrmsd.fit12, color = "M12"), method = "loess", se = FALSE, linewidth = 0.6) +
    geom_smooth(aes(x = dactive, y = nlrmsd_s1, color = "s1"), method = "loess", se = FALSE, linewidth = 0.6) +
    geom_smooth(aes(x = dactive, y = nlrmsd_s2, color = "s2"), method = "loess", se = FALSE, linewidth = 0.6) +
    labs(x = "d", y = "nlRMSD", title = "Distance trend") +
    theme_cowplot(font_size = 8) +
    scale_color_manual(values = CUSTOM_PALETTE) +
    theme(legend.position = "none", plot.title = element_text(face = "plain", color = "#4a4a4a"))
}

# Helper function to create the final layout
create_layout <- function(plots, mcsa_id) {

  first_row_plots <- plot_grid(
    plots$p1, plots$p2, plots$p3,
    ncol = 3,
    rel_widths = c(1.5, 1, 1),
    labels = c("(g)", "(h)", "(i)"),
    label_size = 9,
    hjust = 0
  )

  legend_data <- data.frame(
    x = 1:3,
    y = rep(1, 3),
    label = factor(c("M12", "s1", "s2"), levels = c("M12", "s1", "s2"))
  )

  legend <- ggplot(legend_data) +
    geom_line(aes(x = x, y = y, color = label), size = 1) +
    scale_color_manual(
      values = CUSTOM_PALETTE[c("M12", "s1", "s2")],
      labels = c(expression(M[12]), expression(s[1]), expression(s[2]))
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8)
    ) +
    guides(color = guide_legend(override.aes = list(size = 1)))

  final_plot <- plot_grid(
    first_row_plots,
    legend,
    ncol = 1,
    rel_heights = c( 1.5, 0.2),
    align = 'v'
  )



  return(final_plot)
}
