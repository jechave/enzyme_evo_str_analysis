library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

plot_case <- function(profiles, data_gof, mcsa_id_case) {
  get_pdb_case <- function(profiles, mcsa_id_case) {
    dat <- profiles %>%
      select(mcsa_id, ref_pdb_chain_id)
    dat$ref_pdb_chain_id[[1]]
  }

  pdb_case <- get_pdb_case(profiles, mcsa_id_case)


  source(here("supplementary", "plot_case_obs.R"))
  plot_obs <- plot_case_obs(profiles, data_gof, mcsa_id_case)

  source(here("supplementary", "plot_case_fits_vs_obs.R"))
  plot_fits_vs_obs <- plot_case_fit_vs_obs(profiles, data_gof, mcsa_id_case)

  source(here("supplementary", "plot_case_split.R"))
  plot_split <- plot_case_split(profiles, data_gof, mcsa_id_case)


  title_plot_obs <- ggdraw() +
    draw_label(
      paste("Observed patterns"),
      fontface = "plain",
      size = 10,
      color = "#4a4a4a"
    )

  title_plot_fits_vs_obs <- ggdraw() +
    draw_label(
      paste("Model predictions vs. observations"),
      fontface = "plain",
      size = 10,
      color = "#4a4a4a"
    )

  title_plot_split <- ggdraw() +
    draw_label(
      paste("Decomposition into flexibility and distance contributions"),
      fontface = "plain",
      size = 10,
      color = "#4a4a4a"
    )

  final_plot <- plot_grid(
    title_plot_obs,
    plot_obs,
    ggdraw(),
    title_plot_fits_vs_obs,
    plot_fits_vs_obs,
    ggdraw(),
    title_plot_split,
    plot_split,
    ncol = 1,
    rel_heights = c(0.2, 1.5, 0.2, 0.2, 3.5, 0.2, 0.2, 1.7)
  ) +
    theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"))

  return(final_plot)

}



