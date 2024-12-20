library(tidyverse)
library(ggpubr)
library(ggforce)
library(cowplot)

plot_asite <- function(profiles, data_gof, data_asite_residues) {
  custom_palette <- c(fit0 = "#000000",
                      fit12 = "#31688EFF",
                      fit1 = "#35B779FF",
                      fit2 = "#FDE725FF")

  nunique <- function(x) length(unique(x))

  p1 <- profiles %>%
    left_join(data_gof) %>%
    filter(shell %in% c(0)) %>%
    mutate(nlrmsd.fit0 = 0) %>%
    pivot_longer(
      cols = c(nlrmsd.fit0, nlrmsd.fit1, nlrmsd.fit12),
      names_to = "model",
      values_to = "nlrmsd_model"
    ) %>%
    mutate(model = str_remove(model, "nlrmsd.")) %>%
    mutate(model = factor(model, levels = c("fit0", "fit1", "fit12"))) %>%
    ggplot(aes(model, nlrmsd - nlrmsd_model, color = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_sina(size = 0.7, alpha = 0.2) +
    stat_summary(fun.data = function(x) {
      n <- length(x)
      sd <- sd(x)
      mean <- mean(x)
      return(data.frame(y = mean, ymin = mean-sd, ymax = mean+sd))
    },
    geom = "pointrange", size = 0.5, color = "blue") +
    stat_summary(fun.data = function(x) {
      n <- length(x)
      sd <- sd(x)
      mean <- mean(x)
      return(data.frame(
        y = mean,
        label = sprintf("%.2f ± %.2f", mean, sd)
      ))
    },
    geom = "label",
    color = "black",
    size = 3,
    position = position_nudge(x = 0.3)) +
    theme_cowplot(font_size = 10) +
    stat_compare_means(method = "t.test",
                       paired = TRUE,
                       comparisons = list(c("fit0", "fit1"),
                                          c("fit1", "fit12"),
                                          c("fit0", "fit12")),
                       aes(label = "p.signif"),
                       geom = "label",
                       vjust = 0.5) +
    scale_color_manual(values = custom_palette) +
    scale_x_discrete(labels = c(fit0 = expression(M[0]),
                                fit1 = expression(M[1]),
                                fit12 = expression(M[12]))) +
    labs(y = "Residual (obs. - base.)",
         x = "Baseline") +
    ggtitle("Active site structural divergence") +
    theme(legend.position = "none",
          # axis.title.x = element_blank(),
          plot.title = element_text(face = "plain", color = "#4a4a4a"))


  return(p1)
}
