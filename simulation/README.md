# Simulation studies

This folder contains the simulation studies for the paper *An interpretable semiparametric Bayesian machine
learning model for high-dimensional circulating epigenetic signatures of multiple myeloma* (Section 5 and
Supplementary Section B). Install the modified BART package first (`Rscript BART_package/00_install_BART.R` from the
repository root; see the main [README](../README.md)).

All studies except the semi-synthetic ones simulate their own data and need no real data (see
[Semi-synthetic scenarios](#semi-synthetic-scenarios) below).

## Studies

The paper numbers the studies 1 to 8. The scripts keep their own numbers, which also name the job files in
`../slurm/` and the results folders.

| Script | Paper | What it studies | Replicates and settings | Results folder |
|---|---|---|---|---|
| `14_study_confounder_adjustment.R` | Study 1 | eight linear covariates, a strong age effect and 30 null genes correlated with age | β_age ∈ {0.5, 1.0, 1.5}, ρ ∈ {0.5, 0.8}, p ∈ {2000, 5000}, n = 1500, plus two scenarios with a U-shaped age effect; 14 scenarios × 100 replicates | `results/14_confounder_adjustment/` |
| `12_study_bmi_grid.R` | Study 2 | detection of a small covariate effect as its size varies | β_BMI ∈ {0, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8}, n = 1500, p ∈ {500, 2000}; 100 replicates | `results/12_bmi_grid/` |
| `09_study_semisynthetic.R` | Study 3 | semi-synthetic data: real 5hmC matrix and covariates, simulated outcome | 8 scenarios × 200 replicates | `results/09_semisynthetic/` |
| `01_TrueModel1_regular_BART.R` | Study 4 | BART (Dirichlet prior) with the cross-validation union procedure under DGM 1, in which the covariates act only through gene–covariate interactions (spBART misspecified) | 500 replicates; n ∈ {1500, 2500}, p ∈ {500, 2000, 3000} | `results/TrueModel1/` |
| `02_TrueModel1_spBART.R` | Study 4 | spBART under DGM 1 | as above | `results/TrueModel1/` |
| `03_TrueModel2_regular_BART.R` | Study 4 | BART under DGM 2, with linear covariate effects (spBART correctly specified) | as above | `results/TrueModel2/` |
| `04_TrueModel2_spBART.R` | Study 4 | spBART under DGM 2 | as above | `results/TrueModel2/` |
| `06_TrueModel2_logistic_regression.R` | Study 4 | unpenalized logistic and probit regressions with all genes and the four covariates, under DGM 2 | as above | `results/TrueModel2/` |
| `15_study_ultrahigh_dimension.R` | Study 5 | very high dimension without screening | n ∈ {1000, 1500}, p ∈ {5000, 10000}; 4 scenarios × 100 replicates | `results/15_ultrahigh_dimension/` |
| `07_study_sigma_beta.R` | Study 6 | sensitivity to the prior variance σ²_β ∈ {0.1, 1, 10, 100, 1000}, and the contributions of the linear term and the trees to prediction | 3 scenarios × 100 replicates | `results/07_sigma_beta/` |
| `11_study_multichain.R` | Study 7 | one long MCMC chain (S0) versus 10 chains with the same total number of iterations (S10eq) and 10, 20 or 50 longer chains (S10, S20, S50) | 6 scenarios × 100 replicates | `results/11_multichain/` |
| `08_study_dgm_extensions.R` | Study 8 | smooth nonlinear f, block-correlated genes, covariates correlated with signal genes, and partial misspecification of the linear term | n = 1500, p ∈ {500, 2000}; 100 replicates | `results/08_dgm_extensions/` |
| `10_study_full_pipeline.R` | – | the complete procedure of the application: screening, 5-fold cross-validation with Bayesian FDR, union of the fold selections and final fit | DGM 1 and DGM 2, n = 1500, p ∈ {2000, 5000}; 4 scenarios × 100 replicates | `results/10_full_pipeline/` |

Depending on the study, spBART is compared with BART (`pbart` with the Dirichlet prior), `stan4bart`, the logistic
elastic net (covariates unpenalized) and model-X knockoffs. The header of each script describes its scenarios,
methods and metrics in detail.

Shared code is in `utils/`: `data_generation.R` (DGMs 1 and 2) and `study_functions.R` (the further data-generating
models, method wrappers, metrics, screening and the cross-validation pipeline). The spBART sampler itself is in
`../R/spBART_sampler.R`.

## Running a study locally

Every script takes one integer argument, the task number. For Scripts 01–06 the task number is the replicate, and
one task runs every (n, p) setting. For the other studies each task is one (scenario, replicate) unit, or a small
group of scenarios for one replicate, as described in the header of the script. Run from the repository root:

```bash
Rscript simulation/14_study_confounder_adjustment.R 1   # task 1 of Study 1
Rscript simulation/04_TrueModel2_spBART.R 1             # replicate 1 of spBART under DGM 2 (Study 4)
```

Every script runs the full settings of the paper. Results are saved as one file per replicate,
`results/<study>/<setting>/rep_XXX.RDS`, each holding the selection, prediction and covariate results of every
method for that replicate. The repository contains no results, so `results/` is created and filled when you run
the studies. A task skips any setting whose result file already exists, so an interrupted run resumes where it
stopped.

## Running every replicate on a SLURM cluster

The full studies need far more computing time than a single machine provides. The job scripts in `../slurm/`
(`sim_01_*.sbatch` to `sim_15_*.sbatch`) run each study as a job array, one array task per task number. The job
files do not name an account, so **you need to provide your own cluster credentials** when you submit. Submit all
jobs, with their dependencies, from the `slurm/` folder:

```bash
cd slurm
SPBART_ACCOUNT=<your_account> SPBART_QOS=<your_qos> bash submit_all.sh --dry-run
SPBART_ACCOUNT=<your_account> SPBART_QOS=<your_qos> bash submit_all.sh
```

`SPBART_QOS` is optional; without it, also set `SPBART_MAX_JOBS=<your queued-job limit>`. See the main
[README](../README.md#running-on-a-slurm-cluster) for submitting single jobs and rerunning missing replicates.

## Semi-synthetic scenarios

Study 3 (`09_study_semisynthetic.R`), and one scenario each in Studies 6 (`07_study_sigma_beta.R`) and 7
(`11_study_multichain.R`), use the real 5hmC matrix and covariates of the application, with a simulated outcome.
They read `data_semisynthetic/semisynthetic_base.RDS`, which is built from the real data by
`../real_data_analysis/00_prepare_analysis_base.R`. **The real data are not provided with this repository**, so
these scenarios can only be run by users who have their own copies of the data (see
[`../real_data_analysis/README.md`](../real_data_analysis/README.md)). The file is derived from patient-level data
and is excluded from version control. The other scenarios of Studies 6 and 7, and all other studies, run without it.

## Settings

Unless a study varies them: T = 200 trees, Dirichlet sparsity prior (a = 0.5, b = 1), β ~ N(0, 10 I) for spBART,
the `pbart` leaf prior and offset for both spBART and BART, 2,000 burn-in iterations followed by 5,000 iterations
thinned by 5 (1,000 saved draws), and Bayesian FDR control at 0.05. Prediction metrics (AUC and Brier score) are
computed on an independent test set of n = 1,000 drawn from the same data-generating model; the semi-synthetic
scenarios use the real validation-set predictors (n = 369) instead.
