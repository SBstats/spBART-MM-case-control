# Semi-parametric BART for High-Dimensional Cohort Studies


The repository implements **spBART**, a semi-parametric Bayesian additive regression trees model with a probit link, designed for binary outcomes in high-dimensional cohort studies. It includes simulation studies and real-data application to a pooled multiple-myeloma cohort study.

---

## Repository layout

```
.
├── simulation study/              # Simulation pipeline 
    ├── README.md                  # Detailed usage for the simulation study
    ├── 00–06 *.R / *.sbatch       # Per-method simulation scripts and SLURM jobs
    ├── 05_visualization_*.Rmd     # Compiles per-replicate RDS into manuscript tables
    ├── BART_2.9.9.tar.gz          # Modified BART package source (see below)
    ├── utils/data_generation.R    # True-DGM data generator
    └── results/, logs/            # Output / log directories
└── real data analysis/            # Real-data pipeline
    ├── README.md                  # Detailed usage for the real-data analysis
    ├── 00–05 *.R / *.Rmd          # Preprocessing, validation, spBART fit, comparator, pathway analysis
    └── BART_2.9.9.tar.gz          # Modified BART package source (same tarball as above)
```

Both subdirectories ship with their own `README.md` covering the precise execution order, prerequisites, runtime expectations, and output schemas. **Read those before running any code.**

---

## Modified BART package (required, included as `tar.gz`)

The semi-parametric Gibbs sampler used throughout this repository requires a **modified BART package** that exposes stateful sampler functions. The CRAN BART package will not work.

The modified package source is shipped as `BART_2.9.9.tar.gz` in both the [simulation study/](simulation%20study/) and [real data analysis/](real%20data%20analysis/) folders. **All scripts in this repo are configured to install BART from this local tarball** — there is no remote repository step.

### Installation

Run once, before any analysis script:

```r
install.packages("BART_2.9.9.tar.gz", repos = NULL, type = "source")
```

In each subfolder, `00_setup_BART_package.R` (simulation) and `00_validate_modified_BART_package.Rmd` (real data) handle the install and verify that the four required stateful functions are exported:

- `wbart_create()`
- `wbart_update_response()`
- `wbart_run_iteration()`
- `wbart_destroy()`

If any of these are missing after install, scripts that fit spBART will fail with `could not find function "wbart_create"`.


The modified package is a **fork of the BART R package** by Sparapani et al. (2021). The fork extends the original package by exposing the four stateful sampler functions above so that BART trees can be updated one MCMC iteration at a time inside an outer Gibbs loop — required for the semi-parametric structure where the latent function `f` and the linear coefficients `β` are updated in alternation. No modifications were made to the underlying tree-sampling algorithm; the modifications are restricted to exposing per-iteration entry points to the existing C++ backend.

The DART sparsity prior used throughout (`a = 0.5, b = 1`, `ρ` calibration) follows Linero (2018).

### References

> Sparapani, R., Spanbauer, C., & McCulloch, R. (2021). *Nonparametric Machine Learning and Efficient Computation with Bayesian Additive Regression Trees: The BART R Package.* Journal of Statistical Software, **97**(1), 1–66. <https://doi.org/10.18637/jss.v097.i01>

> Linero, A. R. (2018). *Bayesian Regression Trees for High-Dimensional Prediction and Variable Selection.* Journal of the American Statistical Association, **113**(522), 626–636.

Original BART package: <https://cran.r-project.org/package=BART> (GPL ≥ 2).

---

## Real data — not included

**The real cohort data are not bundled with this repository.** The application uses two controlled-access epigenetic datasets:

- The **UChicago MM (UCMM) Epidemiology Study** (`N₁ = 293` newly diagnosed MM cases)
- The **British Columbia MM Case-control Study** (`N₂ = 576`; 282 cases + 294 controls)

Together these contribute the pooled `N = 869` analytic dataset described in the manuscript.

 The expected file names, locations, and the (hard-coded) data paths used by the scripts are documented in [real data analysis/README.md](real%20data%20analysis/README.md). All paths are of the form `data/<filename>` relative to the script's working directory; you'll need to either place the files into a `data/` subfolder or edit the paths manually.

The simulation study in [simulation study/](simulation%20study/) is fully self-contained — it generates all data internally via `utils/data_generation.R` — and can be run end-to-end without the controlled-access cohort data.

---

## Quick start

1. Install the modified BART package (see above).
2. To reproduce the **simulation tables**, follow [simulation study/README.md](simulation%20study/README.md). On a SLURM cluster, this is roughly:
   ```bash
   cd "simulation study"
   sbatch 01_TrueModel1_regular_BART.sbatch
   sbatch 02_TrueModel1_spBART.sbatch
   sbatch 03_TrueModel2_regular_BART.sbatch
   sbatch 04_TrueModel2_spBART.sbatch
   sbatch 06_TrueModel2_logistic_regression.sbatch
   # After all jobs finish:
   Rscript -e 'rmarkdown::render("05_visualization_comparison.Rmd")'
   ```
3. To reproduce the **real-data analysis** (after obtaining the data), follow [real data analysis/README.md](real%20data%20analysis/README.md):
   ```bash
   cd "real data analysis"
   Rscript 01_Data_preprocessing_pool_cohort_and_case_control.R
   Rscript 03_real_data_analysis_spBART.R
   # Comparator + post-hoc pathway analysis as desired.
   ```

---

