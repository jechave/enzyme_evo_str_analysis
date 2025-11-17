# Supplementary Material

This directory contains the code and scripts to generate the supplementary material PDF for the enzyme evolution structural divergence paper.

## Quick Start

To build the complete supplementary material:

```bash
make
```

This will:
1. Generate all figures (dataset distributions, robustness analyses, case-by-case profiles)
2. Compile the PDF document from `make_supplement.Rmd`
3. Output: `supplement.pdf`

## Prerequisites

### Required R Packages
- `tidyverse` - data manipulation and visualization
- `here` - path management
- `cowplot` - plot composition
- `patchwork` - combining plots
- `ggpubr` - publication-ready plots
- `rmarkdown` - document compilation

### Data Dependencies
The scripts require processed data files from `manuscript/data/`:
- `final_dataset_ca_ref.csv`
- `final_dataset_profiles_ca_ref.csv`
- `final_dataset_gof_ca_ref.csv`

These must be generated first by running the data processing pipeline in the main project.

## Directory Structure

```
manuscript/supplement/
├── Makefile                    # Build automation
├── README.md                   # This file
├── make_supplement.Rmd         # Main RMarkdown document
├── scripts/                    # R scripts for figure generation
│   ├── plot_dataset_distributions.R
│   ├── plot_flex_dist_correlation.R
│   ├── plot_robustness_*.R
│   ├── plot_all_cases.R
│   └── plot_case*.R
└── figures/                    # Generated figures (created by make)
```

## Make Targets

### `make` or `make all` (default)
Generate all figures and compile the PDF document.

### `make figures`
Generate all figures only, without compiling the document.

### `make clean`
Remove all generated files (figures, PDFs, intermediate LaTeX files).

### `make help`
Display available targets and their descriptions.

## Output Files

After running `make`, the following files are generated:

- **`supplement.pdf`** - Final supplementary material document
- **`supplement.tex`** - Intermediate LaTeX source
- **`supplement.log`** - LaTeX compilation log
- **`figures/*.pdf`** - Individual figure files:
  - `dataset_distributions.pdf` - Dataset summary statistics
  - `flex_dist_correlation.pdf` - Flexibility-distance correlation
  - `robustness_protein_choice.pdf` - Robustness analysis for protein selection
  - `robustness_atom_choice.pdf` - Robustness analysis for atom type selection
  - `case_*.pdf` - 34 case-by-case profile figures

## Troubleshooting

### Missing data files
If you get errors about missing CSV files, ensure you've run the main data processing pipeline first:
```bash
cd ../../
# Run data processing scripts in code/Rmd/
```

### Missing R packages
Install required packages:
```r
install.packages(c("tidyverse", "here", "cowplot", "patchwork", "ggpubr", "rmarkdown"))
```

### Stale figures
If figures don't update after data changes:
```bash
make clean
make
```

### Build fails partway through
The Makefile tracks dependencies. You can rebuild specific figures:
```bash
make figures/dataset_distributions.pdf
```
Or rebuild just the document (if figures are already generated):
```bash
rm supplement.pdf
make supplement.pdf
```

## Notes

- The build process can take several minutes, especially generating all 34 case-by-case figures
- Intermediate files (`.tex`, `.log`) are kept for debugging but can be safely removed
- The Makefile uses dependency tracking - only out-of-date targets are rebuilt
