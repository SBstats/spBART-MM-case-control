# Real Data Analysis

Real-data analysis pipeline for the Probit BART high-dimensional cohort study (multiple myeloma case/control, pooled UCMM + BC cohorts, 5hmC gene-body features). 

---

## Data access (read this first)

**No raw data are included in this repository.** The analysis depends on controlled-access datasets from the UChicago MM (UCMM) Study and the British Columbia MM Case-control Study. 


**All file paths are hard-coded** as `data/<filename>` relative to the script's working directory. If your data live elsewhere, either (a) place the files into a `data/` subfolder alongside the scripts, or (b) search each script for `data/` and edit the paths to match your storage location. There is no environment variable or config file to switch paths centrally.

---

## Pipeline / execution order

Run in numeric order. Each script writes its outputs to a `results*/` (or `dart_results*/`) subdirectory created at the script's working directory.

| Order | Script | Purpose |
|---|---|---|
| 1 | `01_Data_preprocessing_pool_cohort_and_case_control.R` | Pool UCMM + BC cohorts; apply exclusion criteria; build case-control design; write Table 1. Produces `pooled_metadata` and `pooled_genebody_data_filtered_normalized`. |
| 2 | `02_naive_analysis_ElasticNet.R` | Naive elastic-net logistic regression comparator (separate penalties for clinical vs. gene features, stability selection). Retained for reference; not reported in the main text. |
| 3 | `03_real_data_analysis_spBART.R` | **Main analysis.** Semi-parametric probit BART (spBART) with two-stage gene screening, 5-fold CV + Bayesian FDR variable selection, refit on development set, validation on held-out set. Produces the results in Section 3 and Tables 8.2–8.4 of the writeup. |
| 4 | `04_comparative_analysis_ElasticNet_vs_spBART_binaryBMInewCutoff_studyIndicatorInBART_Dec05.R` | Side-by-side comparison of elastic net and spBART gene selection at matched FDR. |
| 5 | `05_gene_set_analysis.Rmd` | Post-hoc functional annotation of the spBART-selected genes. Runs GO (BP/MF/CC), KEGG, MSigDB (C2:CGP, C7), and a curated MM-pathway panel. Outputs are persisted to CSV/RDS; the writeup cites the biological interpretation narratively (Section 3, lines 887–888) but does not reproduce formal enrichment tables. |

---

## Prerequisites

### R version
Tested on R ≥ 4.2 (macOS / Linux). Other versions may work but are untested.

### Modified BART package (required)
The analysis depends on a **modified** BART package (`BART_2.9.9.tar.gz`, shipped in this folder) that exposes stateful Gibbs-sampler functions (`wbart_create`, `wbart_update_response`, `wbart_run_iteration`, `wbart_destroy`). The CRAN BART package will not work.

The installation steps:

1. Remove any pre-existing `BART` package from your library.
2. Install from source: `install.packages("BART_2.9.9.tar.gz", repos = NULL, type = "source")`.
3. Verify the four stateful functions exist in the loaded `BART` namespace.

### CRAN packages (scripts 01–04)
`glmnet`, `caret`, `mclust`, `truncnorm`, `pROC`, `readxl`, `dplyr`, `tidyverse`, `ggplot2`, `gridExtra`, `cowplot`, `viridis`, `MASS`

### Bioconductor packages (script 01 for DESeq2; script 05)
`DESeq2` (script 01, variance-stabilizing transformation), and for script 05: `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `msigdbr`, `AnnotationDbi`, `DOSE`, `pathview`, `GOSemSim`.

Script 05 auto-installs missing Bioconductor packages via `BiocManager::install()`. Scripts 01–04 expect required packages to be pre-installed.

---

## Runtime and compute

Script 03 is the expensive one. On a single modern CPU core:

- Stage-1 screening (univariate probit on ~19,000 genes): minutes.
- Stage-2 GMM + 5-fold CV fits of the spBART Gibbs sampler: few hours.
- Final refit + validation-set evaluation: ~30–60 minutes.



Scripts 01, 02, 04, and 05 complete in minutes to ~1 hour each, depending on hardware and available RAM.

---

## Model specification summary

- Probit link, semi-parametric: `Pr(Y=1 | X, Z) = Φ(f(X) + Z'β)`.
- `X` = screened 5hmC gene signatures + study indicator (into the tree ensemble `f`).
- `Z` = four linear covariates: age, sex, race, BMI (BMI is binary at ≥25 kg/m²; see writeup Table 8.4).
- DART (Dirichlet) sparsity prior on splitting variables with `a=0.5, b=1`, `ρ = D+1`.
- `σ_β² = 10` for the linear coefficients.
- `T = 200` trees, 2,000 burn-in, 5,000 post-burn-in, thinning 5, yielding `Q = 1,000` posterior samples.
- Two-stage gene screening: Stage 1 = univariate probit p < 0.05 adjusting for age + sex; Stage 2 = GMM clustering on `|β̂_j|` (2–4 components, BIC, retain strongest cluster).
- 5-fold stratified CV on a development set of `N_Dev = 500`, held-out validation set of `N_Val = 369`.

---

## Outputs

Each script writes to a `dart_results*/` (or `results*/`) folder created next to the script. Key outputs:

- Script 01: `output/table1_demographics.tex`, `pooled_metadata`, normalized gene matrix.
- Script 03: `probit_gmm_screening_results.RDS`, `cv_gene_selection.RDS`, `cv_performance_metrics.RDS`, `beta_posterior_draws.RDS`, `beta_summary.RDS`, validation-set ROC/calibration plots, posterior PIP distributions, MCMC traceplots.
- Script 04: `comparative_analysis_results.RDS`, gene-overlap figures.
- Script 05: per-enrichment CSV and RDS files (GO, KEGG, MSigDB C2:CGP, MSigDB C7, curated MM pathways). The knitted Rmd is lean by design — the formal tables/figures live in the writeup.

---

## Reproducibility notes

- Stratified train/validation split uses `set.seed(789)` (see `03_real_data_analysis_spBART.R`).
- Elastic-net stability selection seeds are `789 + rep` (rep = 1..100); elastic-net final model seed = 456; spBART CV fold creation seed = 890.
- The exact reproducibility of the manuscript's 8-gene list depends on the random-seed state at the data-partitioning step. 
---

## Known gotchas

- Running any script without first installing `BART_2.9.9.tar.gz` will fail with "could not find function wbart_create" or similar.
- The hard-coded `data/<filename>` pattern means `setwd()` matters. Always run from the `real data analysis/` directory, not from a parent folder.
- Script 03 writes a large number of PDF plots (`posterior_pip_distributions.pdf`, `mcmc_traceplots.pdf`, `roc_curves.pdf`, etc.). These can take several minutes to render after the model fits.

