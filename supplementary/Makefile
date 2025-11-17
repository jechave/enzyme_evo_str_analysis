# Makefile for Supplementary Material
# Manages figure generation and document compilation

# Directories
FIGURES_DIR = figures
SCRIPTS_DIR = scripts
DATA_DIR = ../data

# Data files (dependencies)
DATA_FILES = $(DATA_DIR)/final_dataset_ca_ref.csv \
             $(DATA_DIR)/final_dataset_profiles_ca_ref.csv \
             $(DATA_DIR)/final_dataset_gof_ca_ref.csv

# Figure targets
FIG_DATASET = $(FIGURES_DIR)/dataset_distributions.pdf
FIG_FLEXDIST = $(FIGURES_DIR)/flex_dist_correlation.pdf
FIG_ROB_PROT = $(FIGURES_DIR)/robustness_protein_choice.pdf
FIG_ROB_ATOM = $(FIGURES_DIR)/robustness_atom_choice.pdf
FIG_ASITE = $(FIGURES_DIR)/asite_correlations.pdf

# Case figures - we'll use a marker file since there are 34 of them
CASE_MARKER = $(FIGURES_DIR)/.case_figures_done

# All figures
ALL_FIGURES = $(FIG_DATASET) $(FIG_FLEXDIST) $(FIG_ROB_PROT) $(FIG_ROB_ATOM) $(FIG_ASITE) $(CASE_MARKER)

# Final document
DOCUMENT = supplement.pdf

.PHONY: all figures clean help

# Default target
all: $(DOCUMENT)

# Generate all figures
figures: $(ALL_FIGURES)

# Dataset distributions figure
$(FIG_DATASET): $(SCRIPTS_DIR)/plot_dataset_distributions_standalone.R \
                $(SCRIPTS_DIR)/plot_dataset_distributions.R \
                $(DATA_DIR)/final_dataset_ca_ref.csv
	@echo "Generating dataset distributions figure..."
	Rscript $(SCRIPTS_DIR)/plot_dataset_distributions_standalone.R

# Flex-dist correlation figure
$(FIG_FLEXDIST): $(SCRIPTS_DIR)/plot_flex_dist_correlation.R \
                 $(DATA_DIR)/final_dataset_profiles_ca_ref.csv
	@echo "Generating flex-dist correlation figure..."
	Rscript $(SCRIPTS_DIR)/plot_flex_dist_correlation.R

# Robustness protein choice figure
$(FIG_ROB_PROT): $(SCRIPTS_DIR)/plot_robustness_protein_choice.R \
                 $(DATA_FILES)
	@echo "Generating robustness protein choice figure..."
	Rscript $(SCRIPTS_DIR)/plot_robustness_protein_choice.R

# Robustness atom choice figure
$(FIG_ROB_ATOM): $(SCRIPTS_DIR)/plot_robustness_atom_choice.R \
                 $(DATA_FILES)
	@echo "Generating robustness atom choice figure..."
	Rscript $(SCRIPTS_DIR)/plot_robustness_atom_choice.R

# Active site correlations figure
$(FIG_ASITE): $(SCRIPTS_DIR)/plot_asite_correlations.R \
              $(DATA_FILES)
	@echo "Generating active site correlations figure..."
	Rscript $(SCRIPTS_DIR)/plot_asite_correlations.R

# All case figures (34 files)
$(CASE_MARKER): $(SCRIPTS_DIR)/plot_all_cases.R \
                $(SCRIPTS_DIR)/plot_case.R \
                $(SCRIPTS_DIR)/plot_case_obs.R \
                $(SCRIPTS_DIR)/plot_case_fits_vs_obs.R \
                $(SCRIPTS_DIR)/plot_case_split.R \
                $(DATA_FILES)
	@echo "Generating all case-by-case figures..."
	Rscript $(SCRIPTS_DIR)/plot_all_cases.R
	@touch $(CASE_MARKER)

# Compile the document
$(DOCUMENT): make_supplement.Rmd $(ALL_FIGURES)
	@echo "Compiling supplementary material PDF..."
	Rscript -e "rmarkdown::render('make_supplement.Rmd', output_file = 'supplement.pdf')"

# Clean generated files
clean:
	@echo "Cleaning generated files..."
	rm -f $(FIGURES_DIR)/*.pdf
	rm -f $(CASE_MARKER)
	rm -f supplement.pdf supplement.tex supplement.log make_supplement.tex make_supplement.log
	rm -rf make_supplement_files supplement_files

# Help
help:
	@echo "Supplementary Material Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all      - Generate all figures and compile document (default)"
	@echo "  figures  - Generate all figures only"
	@echo "  clean    - Remove all generated files"
	@echo "  help     - Show this help message"
	@echo ""
	@echo "Individual figures:"
	@echo "  $(FIG_DATASET)"
	@echo "  $(FIG_FLEXDIST)"
	@echo "  $(FIG_ROB_PROT)"
	@echo "  $(FIG_ROB_ATOM)"
	@echo "  case figures (34 files)"
