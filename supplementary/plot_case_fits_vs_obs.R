library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)

# Constants
MODEL_NAMES <- c("fit1", "fit2", "fit12")
MODEL_LABELS <- setNames(c(expression(M[1]), expression(M[2]), expression(M[12])), MODEL_NAMES)
COLOR_PALETTE <- c(Observed = "#000000", fit1 = "#35B779FF", fit2 = "#FDE725FF", fit12 = "#31688EFF")
ANNOTATION_SIZE <- 2.5

plot_case_fit_vs_obs <- function(profiles, data_gof, mcsa_id) {
  plot_data <- prepare_plot_data(profiles, data_gof, mcsa_id)
  gof_metrics <- data.frame(
    model = MODEL_NAMES,
    dev_expl = c(plot_data$dev.expl.fit1[1], plot_data$dev.expl.fit2[1], plot_data$dev.expl.fit12[1]),
    aic = c(plot_data$aic.fit1[1], plot_data$aic.fit2[1], plot_data$aic.fit12[1]),
    rmse_lrmsf = c(plot_data$rmse_lrmsf.fit1[1], plot_data$rmse_lrmsf.fit2[1], plot_data$rmse_lrmsf.fit12[1]),
    rmse_dactive = c(plot_data$rmse_dactive.fit1[1], plot_data$rmse_dactive.fit2[1], plot_data$rmse_dactive.fit12[1])
  )

  # Get y-axis limits
  y_limits <- get_y_limits(plot_data)

  # Create faceted plots with standardized y-axis
  p1 <- create_faceted_p1(plot_data, gof_metrics, y_limits)
  p2 <- create_faceted_p2(plot_data, gof_metrics, y_limits)
  p3 <- create_faceted_p3(plot_data, gof_metrics, y_limits)

  # Combine main plots
  main_plots <- plot_grid(p1, p2, p3, ncol = 3,
                          labels = c("(d)", "(e)", "(f)"), label_size = 9,
                          rel_widths = c(1.5, 1, 1),
                          hjust = 0)

  # Create legend
  legend <- create_legend()



  # Combine all elements using a single plot_grid
  final_plot <- plot_grid(
    main_plots,
    legend,
    ncol = 1,
    rel_heights = c(3.0, 0.2),
    align = 'v'
  )



  return(final_plot)
}



# Function to prepare data for plotting
prepare_plot_data <- function(profiles, data_gof, mcsa_id) {
  gof_data <- data_gof %>%
    filter(mcsa_id == !!mcsa_id)
  profiles %>%
    filter(mcsa_id == !!mcsa_id) %>%
    mutate(
      nlrmsd = lrmsd - mean(lrmsd),
      nlrmsf = lrmsf - mean(lrmsf),
      nlrmsd.fit1 = lrmsd.fit1 - mean(lrmsd.fit1),
      nlrmsd.fit2 = lrmsd.fit2 - mean(lrmsd.fit2),
      nlrmsd.fit12 = lrmsd.fit12 - mean(lrmsd.fit12)
    ) %>%
    left_join(gof_data)
}

# Helper function to get y-axis limits
get_y_limits <- function(plot_data) {
  y_min <- min(c(plot_data$nlrmsd, plot_data$nlrmsd.fit1, plot_data$nlrmsd.fit2, plot_data$nlrmsd.fit12), na.rm = TRUE)
  y_max <- max(c(plot_data$nlrmsd, plot_data$nlrmsd.fit1, plot_data$nlrmsd.fit2, plot_data$nlrmsd.fit12), na.rm = TRUE)
  y_range <- y_max - y_min
  return(c(y_min, y_max + 0.1 * y_range))  # Add 10% to the top
}

# Updated create_faceted_p1 function
create_faceted_p1 <- function(plot_data, gof_metrics, y_limits) {
  long_data <- plot_data %>%
    pivot_longer(cols = c(nlrmsd.fit1, nlrmsd.fit2, nlrmsd.fit12),
                 names_to = "model",
                 values_to = "nlrmsd_fit",
                 names_prefix = "nlrmsd.") %>%
    mutate(model = factor(model, levels = MODEL_NAMES))
  active_sites <- plot_data %>%
    filter(dactive == 0) %>%
    select(ref_site) %>%
    distinct()
  gof_metrics <- gof_metrics %>%
    mutate(model = factor(model, levels = MODEL_NAMES)) %>%
    arrange(model)
  ggplot(long_data, aes(x = ref_site)) +
    geom_vline(data = active_sites,
               aes(xintercept = ref_site),
               color = "red", alpha = 0.3, linetype = "dashed") +
    geom_line(aes(y = nlrmsd), color = COLOR_PALETTE["Observed"], linewidth = 0.6) +
    geom_line(aes(y = nlrmsd_fit, color = model), linewidth = 1) +
    facet_wrap(~model, ncol = 1, labeller = labeller(model = MODEL_LABELS)) +
    geom_label(data = gof_metrics,
               aes(x = -Inf, y = Inf,
                   label = sprintf("AIC: %.2f, Expl. Dev: %.2f", aic, dev_expl)),
               hjust = -0.1, vjust = 1.1, size = ANNOTATION_SIZE,
               label.padding = unit(0.15, "lines"),
               label.r = unit(0, "lines"),
               fill = "white",
               color = "black",
               label.size = 0.2) +  # Finer frame
    labs(x = "Residue", y = "nlRMSD", title = "Profiles") +
    theme_cowplot(font_size = 8) +
    scale_color_manual(values = COLOR_PALETTE[-1]) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing = unit(0.2, "lines"),
          plot.title = element_text(face = "plain", color = "#4a4a4a")) +
    coord_cartesian(ylim = y_limits) +
    scale_x_continuous(
      limits = c(0, max(long_data$ref_site)),
      expand = expansion(mult = c(0.05, 0), add = c(0, 0))
    )
}

# Updated create_faceted_p2 function
create_faceted_p2 <- function(plot_data, gof_metrics, y_limits) {
  long_data <- plot_data %>%
    pivot_longer(cols = c(nlrmsd.fit1, nlrmsd.fit2, nlrmsd.fit12),
                 names_to = "model",
                 values_to = "nlrmsd_fit",
                 names_prefix = "nlrmsd.") %>%
    mutate(model = factor(model, levels = MODEL_NAMES))
  gof_metrics <- gof_metrics %>%
    mutate(model = factor(model, levels = MODEL_NAMES)) %>%
    arrange(model)
  ggplot(long_data, aes(x = nlrmsf)) +
    geom_smooth(aes(y = nlrmsd), method = "loess", se = FALSE,
                color = COLOR_PALETTE["Observed"], linetype = "dotted", linewidth = 0.6) +
    geom_smooth(aes(y = nlrmsd_fit, color = model), method = "loess", se = FALSE, linewidth = 1) +
    facet_wrap(~model, ncol = 1, labeller = labeller(model = MODEL_LABELS)) +
    geom_label(data = gof_metrics,
               aes(x = -Inf, y = Inf,
                   label = sprintf("RMSE: %.2f", rmse_lrmsf)),
               hjust = -0.1, vjust = 1.1, size = ANNOTATION_SIZE,
               label.padding = unit(0.15, "lines"),
               label.r = unit(0, "lines"),
               fill = "white",
               color = "black",
               label.size = 0.2) +  # Finer frame
    labs(x = "nlRMSF", y = "nlRMSD", title = "Flexibility trends") +
    theme_cowplot(font_size = 8) +
    scale_color_manual(values = COLOR_PALETTE[-1]) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing = unit(0.2, "lines"),
          plot.title = element_text(face = "plain", color = "#4a4a4a")) +
    coord_cartesian(ylim = y_limits) +
    scale_x_continuous(
      limits = c(min(long_data$nlrmsf), max(long_data$nlrmsf)),
      expand = expansion(mult = c(0.05, 0), add = c(0, 0))
    )
}

# Updated create_faceted_p3 function
create_faceted_p3 <- function(plot_data, gof_metrics, y_limits) {
  long_data <- plot_data %>%
    pivot_longer(cols = c(nlrmsd.fit1, nlrmsd.fit2, nlrmsd.fit12),
                 names_to = "model",
                 values_to = "nlrmsd_fit",
                 names_prefix = "nlrmsd.") %>%
    mutate(model = factor(model, levels = MODEL_NAMES))
  gof_metrics <- gof_metrics %>%
    mutate(model = factor(model, levels = MODEL_NAMES)) %>%
    arrange(model)
  ggplot(long_data, aes(x = dactive)) +
    geom_smooth(aes(y = nlrmsd), method = "loess", se = FALSE,
                color = COLOR_PALETTE["Observed"], linetype = "dotted", linewidth = 0.6) +
    geom_smooth(aes(y = nlrmsd_fit, color = model), method = "loess", se = FALSE, linewidth = 1) +
    facet_wrap(~model, ncol = 1, labeller = labeller(model = MODEL_LABELS)) +
    geom_label(data = gof_metrics,
               aes(x = -Inf, y = Inf,
                   label = sprintf("RMSE: %.2f", rmse_dactive)),
               hjust = -0.1, vjust = 1.1, size = ANNOTATION_SIZE,
               label.padding = unit(0.15, "lines"),
               label.r = unit(0, "lines"),
               fill = "white",
               color = "black",
               label.size = 0.2) +  # Finer frame
    labs(x = "d", y = "nlRMSD", title = "Distance trends") +
    theme_cowplot(font_size = 8) +
    scale_color_manual(values = COLOR_PALETTE[-1]) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing = unit(0.2, "lines"),
          plot.title = element_text(face = "plain", color = "#4a4a4a")) +
    coord_cartesian(ylim = y_limits) +
    scale_x_continuous(
      limits = c(min(long_data$dactive), max(long_data$dactive)),
      expand = expansion(mult = c(0.05, 0), add = c(0, 0))
    )
}

create_legend <- function() {
  legend_data <- data.frame(
    x = 1:5,
    y = rep(1, 5),
    label = factor(c("Observed", MODEL_NAMES, "Active site"),
                   levels = c("Observed", MODEL_NAMES, "Active site"))
  )
  # Create a named vector for labels, using expressions for subscripts
  label_mapping <- c(
    "Observed" = "Obs.",
    setNames(MODEL_LABELS, MODEL_NAMES),
    "Active site" = "Active site"
  )
  ggplot(legend_data) +
    geom_line(aes(x = x, y = y, color = label, linetype = label), size = 0.5) +
    scale_color_manual(
      values = c(COLOR_PALETTE, "Active site" = "red"),
      labels = label_mapping,
      breaks = c("Observed", MODEL_NAMES, "Active site")
    ) +
    scale_linetype_manual(
      values = c(Observed = "solid", fit1 = "solid", fit2 = "solid", fit12 = "solid", "Active site" = "dashed"),
      labels = label_mapping,
      breaks = c("Observed", MODEL_NAMES, "Active site")
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.3, "cm"),
      legend.key.width = unit(0.6, "cm"),
      legend.spacing.x = unit(0.1, "cm")
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 0.5)),
      linetype = guide_legend(override.aes = list(size = 0.5))
    )
}

