library(tidyverse)
library(ggpubr)
library(ggforce)
library(cowplot)
library(ggrepel)

plot_asite <- function(profiles, data_gof) {
  # Filter out outlier case 749 only
  profiles <- profiles %>% filter(mcsa_id != "749")
  data_gof <- data_gof %>% filter(mcsa_id != "749")

  custom_palette <- c(fit0 = "#000000",
                      fit12 = "#31688EFF",
                      fit1 = "#35B779FF",
                      fit2 = "#FDE725FF")

  # Contribution colors from plot_split.R
  contrib_palette <- c(s1 = "#35B779FF", s2 = "#FDE725FF", s1_plus_s2 = "#31688EFF")

  # Colors for null, observed, and predicted
  panel_a_palette <- c(null = "#969696", observed = "#000000", predicted = "#31688EFF")

  nunique <- function(x) length(unique(x))

  # Panel (a): Scatter plot of observed vs M12-predicted active site divergence
  # FAMILY-AVERAGED VERSION: Average across active site residues within each family
  panel_a_data <- profiles %>%
    left_join(data_gof) %>%
    group_by(mcsa_id) %>%
    mutate(
      nlrmsd = lrmsd - mean(lrmsd, na.rm = TRUE),
      nlrmsd.fit12 = lrmsd.fit12 - mean(lrmsd.fit12, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    filter(near(dactive, 0)) %>%
    # FAMILY AVERAGING: Average divergence across active site residues for each family
    group_by(mcsa_id) %>%
    summarize(
      observed = mean(nlrmsd, na.rm = TRUE),
      predicted = mean(nlrmsd.fit12, na.rm = TRUE),
      .groups = "drop"
    )

  # Calculate statistics for annotation
  cor_result <- cor.test(panel_a_data$observed, panel_a_data$predicted, method = "pearson")
  obs_mean <- mean(panel_a_data$observed)
  obs_sd <- sd(panel_a_data$observed)
  pred_mean <- mean(panel_a_data$predicted)
  pred_sd <- sd(panel_a_data$predicted)

  # Create annotation with proper formatting (separate lines to avoid parsing issues)
  annotation_text <- sprintf(
    "R = %.2f\nObs: %.2f ± %.2f\nM12: %.2f ± %.2f",
    cor_result$estimate, obs_mean, obs_sd, pred_mean, pred_sd
  )

  # Identify outlier points (observed > -0.1, i.e., close to 0 or positive)
  outlier_points <- panel_a_data %>% filter(observed > -0.1)
  outlier_ids <- outlier_points$mcsa_id

  panel_a <- ggplot(panel_a_data, aes(x = predicted, y = observed)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 1.5, alpha = 0.6, color = "#31688EFF") +
    geom_text_repel(data = outlier_points, aes(label = mcsa_id),
                    size = 2.5, color = "black",
                    box.padding = 0.2, point.padding = 0.2,
                    min.segment.length = 0, segment.size = 0.3) +
    annotate("label", x = -1.2, y = 0.3, label = annotation_text,
             hjust = 0, vjust = 1, size = 2.75, lineheight = 1,
             fill = "white", label.size = 0.2) +
    scale_x_continuous(limits = c(-1.2, 0.3), breaks = c(-1, -0.5, 0)) +
    scale_y_continuous(limits = c(-1.2, 0.3), breaks = c(-1, -0.5, 0)) +
    theme_cowplot(font_size = 9) +
    labs(x = expression(M[12] ~ nlRMSD),
         y = "Observed nlRMSD",
         title = "Active site divergence") +
    theme(plot.title = element_text(face = "plain", color = "#4a4a4a"),
          panel.grid.major = element_line(color = "gray90"))

  # Panel (b): Constraint contributions for active site residues
  # FAMILY-AVERAGED VERSION: Average across active site residues within each family
  panel_b_data <- profiles %>%
    group_by(mcsa_id) %>%
    mutate(
      s1 = lrmsd.fit12.c1 - mean(lrmsd.fit12.c1, na.rm = TRUE),
      s2 = lrmsd.fit12.c2 - mean(lrmsd.fit12.c2, na.rm = TRUE),
      s1_plus_s2 = s1 + s2
    ) %>%
    ungroup() %>%
    filter(near(dactive, 0)) %>%
    select(mcsa_id, ref_pdb_chain_id, ref_pdb_site, s1, s2, s1_plus_s2) %>%
    pivot_longer(
      cols = c(s1, s2, s1_plus_s2),
      names_to = "component",
      values_to = "value"
    ) %>%
    mutate(component = factor(component, levels = c("s1_plus_s2", "s1", "s2"))) %>%
    # FAMILY AVERAGING: Average component values across active site residues for each family
    group_by(mcsa_id, component) %>%
    summarize(value = mean(value, na.rm = TRUE), .groups = "drop")

  panel_b <- ggplot(panel_b_data, aes(component, value, color = component)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_sina(size = 0.7, alpha = 0.7) +
    stat_summary(fun.data = function(x) {
      n <- length(x)
      sd <- sd(x)
      mean <- mean(x)
      return(data.frame(y = mean, ymin = mean-sd, ymax = mean+sd))
    },
    geom = "pointrange", size = 0.3, color = "blue") +
    theme_cowplot(font_size = 9) +
    stat_compare_means(method = "t.test",
                       paired = TRUE,
                       comparisons = list(c("s1", "s1_plus_s2"),
                                          c("s2", "s1_plus_s2"),
                                          c("s1", "s2")),
                       aes(label = "p.signif"),
                       geom = "label",
                       vjust = 0,
                       size = 2.3,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, Inf),
                                          symbols = c("****", "***", "**", "ns"))) +
    scale_color_manual(values = contrib_palette) +
    scale_x_discrete(labels = c(s1 = expression(s[1]),
                                s2 = expression(s[2]),
                                s1_plus_s2 = expression(s[1] + s[2])),
                     expand = expansion(add = c(0.6, 0.9))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(y = "Component value", x = "Component",
         title = "Constraint contributions") +
    theme(legend.position = "none",
          plot.title = element_text(face = "plain", color = "#4a4a4a"))

  # Panel (c): Relative contributions of s1 and s2 to active site conservation
  # Calculate proportion of MAGNITUDE for each residue, then average per family
  # Using absolute values to show which constraint dominates in magnitude
  # Sort by active-site RSC(s2)
  panel_c_data <- profiles %>%
    group_by(mcsa_id) %>%
    mutate(
      s1 = lrmsd.fit12.c1 - mean(lrmsd.fit12.c1, na.rm = TRUE),
      s2 = lrmsd.fit12.c2 - mean(lrmsd.fit12.c2, na.rm = TRUE)
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
      prop_s1_avg = mean(prop_s1, na.rm = TRUE),
      prop_s2_avg = mean(prop_s2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(prop_s2_avg) %>%
    mutate(mcsa_id = factor(mcsa_id, levels = mcsa_id)) %>%
    pivot_longer(
      cols = c(prop_s1_avg, prop_s2_avg),
      names_to = "component",
      values_to = "value"
    ) %>%
    mutate(component = ifelse(component == "prop_s1_avg", "s1", "s2"))

  panel_c <- ggplot(panel_c_data, aes(x = mcsa_id, y = value, fill = component)) +
    geom_col(position = "stack", color = contrib_palette["s1_plus_s2"], linewidth = 0.3) +
    scale_fill_manual(values = contrib_palette) +
    labs(x = "M-CSA ID", y = "Relative contribution",
         title = "Relative contributions") +
    theme_cowplot(font_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          legend.position = "none",
          plot.title = element_text(face = "plain", color = "#4a4a4a")) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                       labels = scales::number_format(accuracy = 0.1))

  # Combine panels vertically
  combined_plot <- plot_grid(panel_a, panel_b, panel_c,
                             ncol = 1,
                             labels = c("(a)", "(b)", "(c)"),
                             label_size = 10)

  return(combined_plot)
}
