# Simulation Study

Simulation pipeline for the semi-parametric BART model. 

---

## What this folder contains

| File / dir | Purpose |
|---|---|
| `00_setup_BART_package.R` | One-time install + sanity check of the modified BART package |
| `01_TrueModel1_regular_BART.R` | Standard BART (DART) under DGM 1 |
| `02_TrueModel1_spBART.R` | Semi-parametric BART under DGM 1 |
| `03_TrueModel2_regular_BART.R` | Standard BART (DART) under DGM 2 |
| `04_TrueModel2_spBART.R` | Semi-parametric BART under DGM 2 |
| `06_TrueModel2_logistic_regression.R` | Logistic + probit GLM baselines under DGM 2 (Appendix C) |
| `05_visualization_comparison.Rmd` | Compiles per-replicate `.RDS` files into the manuscript tables |
| `*.sbatch` | SLURM batch scripts (one per simulation script) |
| `BART_2.9.9.tar.gz` | Modified BART package source (required) |
| `utils/data_generation.R` | True-DGM data generator used by all simulation scripts |
| `results/` | Output directory for per-replicate `.RDS` files (empty initially) |
| `logs/` | SLURM stdout/stderr logs |

---

## Simulation design 

- **Sample sizes:** `n ∈ {1500, 2500}`
- **Gene-set sizes:** `p ∈ {500, 2000, 3000}`
- **True signal genes:** 10 (always genes 1–10, fixed across replicates)
- **Covariates (linear or split):** age, sex, race, BMI (all standardized / encoded as in Sec 7.1)
- **Replicates:** 500 per (DGM, n, p) cell
- **Two data-generating models:**
  - **DGM 1 — Pure sum of trees.** 20 trees: 10 marginal-gene depth-1, 10 gene-covariate depth-2 interactions. **No linear Z'β.** spBART is misspecified.
  - **DGM 2 — Gene-only trees + linear covariates.** 15 gene-only trees (10 marginal + 5 gene-gene) + linear `Z'β = (0.50, 0.60, -0.40, 0.45)` for (age, sex, race, BMI). spBART is correctly specified.
- **Outcome:** `Y_i ~ Bernoulli(Φ(f(X_i, Z_i)))` (probit link).
- **Methods compared:** standard BART vs spBART (both DGMs); logistic + probit GLM baselines under DGM 2 only (Appendix C).

---

## Model / MCMC settings 

| Parameter | Value | Where |
|---|---|---|
| Number of trees `T` | 200 | scripts 01–04 |
| Burn-in iterations | 2,000 | scripts 01–04 |
| Post-burn-in iterations | 5,000 | scripts 01–04 |
| Thinning | 5 → `Q = 1,000` saved samples | scripts 01–04 |
| DART hyperparameters | `a = 0.5, b = 1` | scripts 01–04 |
| DART concentration `ρ` | `p_genes` (no study indicator in simulation) | scripts 02, 04 |
| Linear-coefficient prior | `β ~ N(0, σ_β² I)` with `σ_β² = 10` | scripts 02, 04 |
| FDR control level | `α = 0.05` (Bayesian FDR on PIPs) | all scripts |
| CV folds | 5 (stratified) | all scripts |

The MCMC schedule above gives total iterations `2,000 + 5,000 = 7,000` per fit.

---

## Pipeline / execution order

Run in numeric order. Each simulation script is independent — they can be parallelized across SLURM jobs.

1. **`00_setup_BART_package.R`** — install the modified BART package and verify the four stateful sampler functions exist. Run once before anything else.
2. **`01–04`** — fit the four BART variants. Each script writes one `.RDS` per replicate to `results/TrueModel<N>/n<n>/p<p>/<method>/rep_<NNN>.RDS`.
3. **`06`** — fit the GLM baselines for DGM 2. Writes to `results/TrueModel2/n<n>/p<p>/{logistic_reg,probit_reg}/rep_<NNN>.RDS`.
4. **`05_visualization_comparison.Rmd`** — knit to produce the manuscript tables. Reads from `results/`.

### Two execution modes

- **SLURM cluster (recommended for full 500-replicate runs):** submit each `.sbatch` as a job array. Each array task is one replicate; each replicate iterates over all (n, p) combinations sequentially. Example:
  ```bash
  sbatch 01_TrueModel1_regular_BART.sbatch   # 500 array tasks
  sbatch 02_TrueModel1_spBART.sbatch
  sbatch 03_TrueModel2_regular_BART.sbatch
  sbatch 04_TrueModel2_spBART.sbatch
  sbatch 06_TrueModel2_logistic_regression.sbatch
  ```
  Adjust the `--account=` and `--qos=` lines at the top of each sbatch file to your cluster's allocation.

- **Local (for debugging or small reps):** invoke the R scripts directly. They read the SLURM array index from the environment variable `SLURM_ARRAY_TASK_ID`; if absent, set it manually:
  ```bash
  SLURM_ARRAY_TASK_ID=1 Rscript 02_TrueModel1_spBART.R
  ```

### Resource expectations (per replicate, all (n, p) cells)

| Script | Wall time | Memory |
|---|---|---|
| 01 / 03 (Regular BART) | ~5 hours | 32 GB |
| 02 / 04 (spBART) | ~6 hours | 32 GB |
| 06 (logistic + probit) | ~30 minutes | 32 GB |



---

## Prerequisites

### R version
Tested on R ≥ 4.2 (Linux HPC and macOS). Other versions may work but are untested.

### Modified BART package (required)
The simulation depends on a **modified** BART package (`BART_2.9.9.tar.gz`, included in this folder) that exposes stateful Gibbs-sampler functions (`wbart_create`, `wbart_update_response`, `wbart_run_iteration`, `wbart_destroy`). The CRAN BART package will **not** work for `02_TrueModel1_spBART.R`, `04_TrueModel2_spBART.R`, or any script that calls `fit_semiparametric_probit_DART()`.

Run `00_setup_BART_package.R` once to handle installation.

### CRAN packages
`BART` (loaded after the modified install above), `glmnet`, `caret`, `mclust`, `truncnorm`, `pROC`, `MASS`, `tidyverse`, `knitr`, `kableExtra`. Used by simulation scripts and `05_visualization_comparison.Rmd`.


---

## Reproducibility notes

- Seeds are deterministic per (DGM, n, p, replicate):
  - DGM 1 scripts use `set.seed(1000 * n + p_genes + rep_id)`.
  - DGM 2 scripts use `set.seed(2000 * n + p_genes + rep_id)`.
  - DGM 1 replicate `k` and DGM 2 replicate `k` therefore run on **different** simulated datasets and are not paired.
- Cutpoints in the true DGMs are fixed population quantiles from `N(0,1)` (see `utils/data_generation.R`), not sample quantiles, so the true regression function is identical across replicates.

---

## Known gotchas

- Running any spBART script (`02`, `04`) without first installing `BART_2.9.9.tar.gz` will fail with `could not find function "wbart_create"` or similar. Always run `00_setup_BART_package.R` first.
- The `.sbatch` files set `--account=` and `--qos=` to empty strings. You **must** fill these in for your cluster, or remove the lines.
- `results/` is created lazily by the simulation scripts; it ships empty.
- The replicate output is one `.RDS` per replicate per (DGM, n, p, method). A full run produces 500 × 6 (n, p) × 5 (methods) ≈ 15,000 small files. Plan disk and inode budgets accordingly.

