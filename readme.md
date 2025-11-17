# Reproducibility Package

This repository contains the data and code to reproduce all figures and tables from our paper.

## Repository Structure

```
.
├── data/           # Final processed datasets
├── scripts/        # Figure and table generation scripts
├── figures/        # Main paper figures (generated)
├── tables/         # Tables (generated)
└── supplementary/  # Supplementary materials
    ├── scripts/    # Supplementary figure scripts
    ├── figures/    # Supplementary figures (generated)
    └── Makefile    # Build system for supplement
```

## Data Files

All data files required to reproduce the figures are in `data/`. Three analysis variants are provided (ca_ref is the primary analysis):

- **ca_ref**: Reference protein structure, C-alpha atoms (primary)
- **ca_homolog**: Closest homolog structure, C-alpha atoms
- **cb_ref**: Reference protein structure, C-beta atoms

### Primary Data Files (ca_ref variant)

- `final_dataset_ca_ref.csv`: Family-level properties and annotations
- `final_dataset_gof_ca_ref.csv`: Goodness-of-fit metrics for models M1, M2, and M12
- `final_dataset_profiles_ca_ref.csv`: Residue-level structural profiles including:
  - Divergence: RMSD, lRMSD, nlRMSD (observed and model-predicted)
  - Flexibility: lRMSF, nlRMSF
  - Distance to active site: dactive
  - Constraint decomposition: s1, s2 (from Shapley analysis)
- `final_dataset_residues.csv`: Active site residue annotations from M-CSA database
- `final_dataset_table.csv`: Summary table of dataset families
- `final_dataset_pdbs.rda`: PDB metadata (R data file)

The same files exist for ca_homolog and cb_ref variants.

## Reproducing Figures and Tables

### Main Figures

Generate all main paper figures:
```bash
Rscript scripts/make_figures.R
```

This creates in `figures/`:
- `fig_obs.pdf` - Observed divergence and flexibility profiles
- `fig_fit_vs_obs.pdf` - Model predictions vs observations
- `fig_split.pdf` - Split constraint contributions (s1 and s2)
- `fig_shapley_map.pdf` - Shapley value map across families
- `fig_asite.pdf` - Active site analysis

### Tables

Generate the robustness table:
```bash
Rscript scripts/make_table_robustness.R
```

This creates `tables/table_robustness.tex` comparing model performance across the three analysis variants.

### Supplementary Materials

The supplement uses a Makefile-based build system:

```bash
cd supplementary
make          # Generate all figures and compile PDF
make figures  # Generate figures only
make clean    # Remove generated files
```

This produces `supplementary/supplement.pdf` with all supplementary figures and tables.

## Requirements

R packages:
- `tidyverse`
- `here`
- `cowplot`
- `ggpubr`
- `ggforce`
- `ggrepel`
- `knitr`
- `rmarkdown`

For the supplement, you also need LaTeX installed.

## Notes

- All scripts use the `here` package for path management, so they work from any working directory within the repository
- The `data/` folder contains final processed datasets ready for analysis
- Data preparation scripts are not included (only final data is provided for reproducibility)
