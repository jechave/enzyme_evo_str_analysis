library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)
library(purrr)

plot_shapley_map <- function(
    data_gof,
    dataset,
    mcsa_id_examples = NULL,
    mcsa_id_outliers = NULL,
    color_palette = c(fit1 = "#35B779FF",
                      fit2 = "#FDE725FF",
                      fit12 = "#31688EFF"),
    xlim = c(-.02, .6),
    ylim = c(-.02, .6),
    font_size = 12
) {

  # Join datasets
  plot_data <- data_gof %>%
    left_join(dataset) %>%
    mutate(best_model = factor(best_model, levels = c("fit1", "fit2", "fit12")))

  # Create base plot
  p <- plot_data %>%
    ggplot(aes(shapley_lrmsf, shapley_dactive)) +
    geom_abline(color = "grey") +
    geom_point(aes(color = best_model), size = 2.5)

  # Add example labels if provided
  if (!is.null(mcsa_id_examples)) {
    p <- p + geom_label_repel(
      data = filter(plot_data, mcsa_id %in% mcsa_id_examples),
      aes(label = mcsa_id)
    )
  }

  # Add outlier labels if provided
  if (!is.null(mcsa_id_outliers)) {
    p <- p + geom_label_repel(
      data = filter(plot_data, mcsa_id %in% mcsa_id_outliers),
      aes(label = mcsa_id),
      color = "red"
    )
  }

  # Add styling
  p <- p +
    scale_color_manual(
      values = color_palette,
      labels = c(expression(M[1]), expression(M[2]), expression(M[12]))
    ) +
    theme_cowplot(font_size = font_size) +
    coord_fixed(ratio = 1) +
    xlim(xlim) +
    ylim(ylim) +
    labs(
      x = expression("Non-functional constr. SC(s"[1]*")"),
      y = expression("Functional constr. SC(s"[2]*")"),
      color = "Best Model"
    )

  return(p)
}
