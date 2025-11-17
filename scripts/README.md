# Figure and Table Generation Scripts

This directory contains scripts to generate the main manuscript figures and tables from the final processed data.

## Requirements

R packages:
- `tidyverse`
- `here`
- `cowplot`
- `ggpubr`
- `ggforce`
- `ggrepel`

## Usage

### Generate all main figures

```bash
Rscript make_figures.R
```

This generates all main paper figures in `../figures/`:
- `fig_obs.pdf` - Observed divergence and flexibility profiles
- `fig_fit_vs_obs.pdf` - Model predictions vs observations
- `fig_split.pdf` - Split constraint contributions (s1 and s2)
- `fig_shapley_map.pdf` - Shapley value map across families
- `fig_asite.pdf` - Active site analysis

### Generate robustness table

```bash
Rscript make_table_robustness.R
```

This generates `../tables/table_robustness.tex` comparing model performance across three analysis variants (ca_ref, ca_homolog, cb_ref).

## Data Dependencies

All scripts read from `../data/`:
- `final_dataset_ca_ref.csv` - Family-level properties (primary analysis)
- `final_dataset_gof_ca_ref.csv` - Model goodness-of-fit metrics
- `final_dataset_profiles_ca_ref.csv` - Residue-level profiles

The robustness table also uses the corresponding files for ca_homolog and cb_ref variants.

## Plotting Functions

Individual plotting functions are defined in:
- `plot_obs.R`
- `plot_fits_vs_obs.R`
- `plot_split.R`
- `plot_shapley_map.R`
- `plot_asite.R`

These are sourced by `make_figures.R` and should not be run directly.
