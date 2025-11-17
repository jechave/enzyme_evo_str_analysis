# Constants
ANNOTATION_SIZE <- 2.5
CUSTOM_PALETTE <- c(Observed = "#000000", M12 = "#31688EFF", s1 = "#35B779FF", s2 = "#FDE725FF")

# Main function to create the plot
plot_split <- function(profiles, data_gof, mcsa_id) {
  prepared_data <- prepare_data(profiles, data_gof, mcsa_id)
  plots <- create_plots(prepared_data)
  create_layout(plots, mcsa_id)
}

# Helper function to prepare data
prepare_data <- function(profiles, data_gof, mcsa_id) {
  data_gof <- data_gof %>%
    left_join(profiles %>% group_by(mcsa_id) %>% summarize(
      lrmsf_sd = sd(lrmsf),
      dactive_mean = mean(dactive)
    ), by = "mcsa_id")

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

  list(plot_data = plot_data, gof_data = gof_data, data_gof = data_gof)
}

# Helper functions for creating individual plots
create_plots <- function(prepared_data) {
  plot_data <- prepared_data$plot_data
  gof_data <- prepared_data$gof_data
  data_gof <- prepared_data$data_gof

  p1 <- create_panel1(plot_data, gof_data)
  p2 <- create_panel2(plot_data)
  p3 <- create_panel3(plot_data)
  p4 <- create_panel4(data_gof)
  p5 <- create_panel5(data_gof)
  p6 <- create_panel6(data_gof)

  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6)
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
      #label = sprintf("Expl. Dev: %.2f, SC(s1): %.2f, SC(s2): %.2f", gof_data$dev.expl.fit12, gof_data$shapley_lrmsf, gof_data$shapley_dactive),
      label = as.expression(bquote(paste("Expl. Dev: ",
                                         .(sprintf("%.2f", gof_data$dev.expl.fit12)),
                                         ", SC(s"[1], "): ",
                                         .(sprintf("%.2f", gof_data$shapley_lrmsf)),
                                         ", SC(s"[2], "): ",
                                         .(sprintf("%.2f", gof_data$shapley_dactive))
      ))),
      hjust = -0.05, vjust = 1.1, size = 2.5,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0, "lines"),
      fill = "white",
      color = "black",
      label.size = 0.2,
      inherit.aes = FALSE
    ) +
    labs(x = "Residue", y = "lRMSD", title = "Profile") +
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

create_panel4 <- function(data_gof) {
  p4_data <- data_gof %>%
    arrange(prop_dactive) %>%
    mutate(mcsa_id = factor(mcsa_id, levels = mcsa_id)) %>%
    pivot_longer(cols = c(prop_lrmsf, prop_dactive), names_to = "component", values_to = "value") %>%
    mutate(component = ifelse(component == "prop_lrmsf", "s1", "s2"))
  ggplot(p4_data, aes(x = mcsa_id, y = value, fill = component)) +
    geom_col(position = "stack", color = CUSTOM_PALETTE["M12"]) +
    scale_fill_manual(values = CUSTOM_PALETTE) +
    labs(x = "M-CSA ID", y = "RSC",
         title = expression("Relative contrib. of s"[1]*" and s"[2])) +
           theme_cowplot(font_size = 8) +  # Changed to 7 as requested
           theme(axis.text.x = element_text(angle = 90, hjust = 0.0, vjust = 1.0, size = 5),
                 legend.position = "none", plot.title = element_text(face = "plain", color = "#4a4a4a")) +
           scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
           scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                              labels = scales::number_format(accuracy = 0.1))  # One decimal place
}

create_panel5 <- function(data_gof) {
  p <- ggplot(data_gof, aes(x = lrmsf_sd, y = shapley_lrmsf)) +
    geom_point(color = CUSTOM_PALETTE["s1"]) +
    geom_smooth(method = "lm", se = FALSE, color = CUSTOM_PALETTE["s1"]) +
    #labs(x = "sd(lRMSF)", y = "SC(s1)",
    labs(x = "sd(lRMSF)", y = expression("SC(s"[1]*")"),
         #title = "s1 contrib.") +
         title = expression("s"[1]*" contrib.")) +
    theme_cowplot(font_size = 8) +
    scale_y_continuous(expand = expansion(mult = c(0.0, 0.2), add = c(0, 0)))

  add_stats_annotation(p, data_gof$lrmsf_sd, data_gof$shapley_lrmsf)
}

create_panel6 <- function(data_gof) {
  p <- ggplot(data_gof, aes(x = dactive_mean, y = shapley_dactive)) +
    geom_point(color = CUSTOM_PALETTE["s2"]) +
    geom_smooth(method = "lm", se = FALSE, color = CUSTOM_PALETTE["s2"]) +
    #labs(x = expression("mean(d)"), y = expression("SC(s2)"),
    labs(x = "sd(lRMSF)", y = expression("SC(s"[2]*")"),
         #title = "s2 contrib.") +
         title = expression("s"[2]*" contrib.")) +
    theme_cowplot(font_size = 8) +
    scale_y_continuous(expand = expansion(mult = c(0.0, 0.2), add = c(0, 0)))

  add_stats_annotation(p, data_gof$dactive_mean, data_gof$shapley_dactive)
}

add_stats_annotation <- function(plot, x, y) {
  lm_fit <- lm(y ~ x)
  r_squared <- summary(lm_fit)$r.squared
  p_value <- summary(lm_fit)$coefficients[2, 4]

  p_value_text <- if(p_value < 0.001) "p < 0.001" else sprintf("p = %.3f", p_value)

  plot +
    geom_label(
      aes(x = -Inf, y = Inf),
      label = sprintf("R² = %.2f, %s", r_squared, p_value_text),
      hjust = -0.1, vjust = 1.1, size = ANNOTATION_SIZE,
      label.padding = unit(0.15, "lines"),
      label.r = unit(0, "lines"),
      fill = "white",
      color = "black",
      label.size = 0.2,
      inherit.aes = FALSE
    ) +
    theme(plot.title = element_text(face = "plain", color = "#4a4a4a"))
}

# Helper function to create the final layout
create_layout <- function(plots, mcsa_id) {
  #main_title <- ggdraw() + draw_label("Decomposition into s1 and s2 contributions",
  main_title <- ggdraw() + draw_label(expression(paste("Decomposition into s"[1], " and s"[2], " contributions")),
                                      fontface = "plain", size = 12, color = "#4a4a4a")

  first_row_title <- ggdraw() + draw_label(paste0("Example family: M-CSA ID = ", mcsa_id),
                                           fontface = "plain", size = 10, color = "#4a4a4a")

  second_row_title <- ggdraw() + draw_label("Aggregate results across families",
                                            fontface = "plain", size = 10, color = "#4a4a4a")

  first_row_plots <- plot_grid(plots$p1, plots$p2, plots$p3, ncol = 3,
                               rel_widths = c(1.5, 1, 1), labels = c("(a)", "(b)", "(c)"),
                               label_size = 9)

  second_row_plots <- plot_grid(plots$p4, plots$p5, plots$p6, ncol = 3,
                                rel_widths = c(1.5, 1, 1), labels = c("(d)", "(e)", "(f)"),
                                label_size = 9)



  ## legend for first row of plots

  legend_data <- data.frame(x = 1:3, y = rep(1, 3),
                            label = factor(c("M12", "s1", "s2"),
                                           levels = c("M12", "s1", "s2")))

  legend_first_row <- ggplot(legend_data) +
    geom_line(aes(x = x, y = y, color = label), size = 1) +
    scale_color_manual(values = CUSTOM_PALETTE[c("M12", "s1", "s2")],
                       labels = c(expression(M[12]), expression(s[1]), expression(s[2]))) +
    theme_void() +
    theme(legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 8)) +
    guides(color = guide_legend(override.aes = list(size = 1)))


  # legend for second row of plots
  legend_data <- data.frame(x = 1:2, y = rep(1, 2),
                            label = factor(c("s1", "s2"), levels = c("s1", "s2")))

  legend_second_row <- ggplot(legend_data) +
    geom_line(aes(x = x, y = y, color = label), size = 1) +
    scale_color_manual(values = CUSTOM_PALETTE[c("s1", "s2")],
                       labels = c("s1 contrib.", "s2 contrib.")) +
    theme_void() +
    theme(legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 8)) +
    guides(color = guide_legend(override.aes = list(size = 1)))


  ## Only one legend

  legend <- legend_first_row

  plot_grid(
    main_title,
    first_row_title,
    first_row_plots,
    ggdraw(),
    second_row_title,
    second_row_plots,
    legend,
    ncol = 1, rel_heights = c(0.3, 0.3, 1.5, 0.2, 0.3, 1.5, 0.2), align = 'v'
  ) +
    theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"))


  # plot_grid(
  #   main_title,
  #   first_row_title,
  #   first_row_plots,
  #   legend_first_row,
  #   ggdraw(),
  #   second_row_title,
  #   second_row_plots,
  #   legend_second_row,
  #   ncol = 1, rel_heights = c(0.3, 0.3, 1.5, 0.2, 0.2, 0.3, 1.5, 0.2), align = 'v'
  # ) +
  #   theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"))
}
