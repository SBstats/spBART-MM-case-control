# Real-data application

This folder contains the real-data analysis of the paper *An interpretable semiparametric Bayesian machine learning
model for high-dimensional circulating epigenetic signatures of multiple myeloma* (Section 4): a pooled case–control
analysis of gene-level 5-hydroxymethylcytosine (5hmC) profiles from circulating cell-free DNA. It pools multiple
myeloma cases from the University of Chicago MM Epidemiology (MMEP) study (labelled UCMM in the code) with cases and
controls from the British Columbia case–control study (BC).

## No real data are provided

**This repository does not contain any of the real data used in the paper, nor any results derived from them.** The
data are patient-level and cannot be shared. To run these scripts you must have your own copies of the data and
place them in `data/` (this folder is excluded from version control). Without them, none of the scripts in this
folder will run.

The scripts expect the following files in `data/`:

| File | Contents | Fields used |
|---|---|---|
| `ERRC.ID.xlsx` | map from MMEP patient ID to 5hmC sequencing ID | `ERRC.ID`, `Sequencing.ID` |
| `genebody_377.RDS` | MMEP gene-body 5hmC read counts (genes × sequencing IDs) | gene names as row names |
| `patient_metadata_full_n797.RData` | MMEP patient metadata (data frame `analysis_data`) | `errcid`, `age_diag`, `sex_emr`, `race_composite`, `bmi_dx_emr`, `dx_errc` |
| `Canada_case_control_genebody_count.RDS` | BC gene-body 5hmC read counts (genes × sample IDs) | gene names as row names |
| `MM_Questionnaire_BC_for_BART_projects.RDS` | BC questionnaire data for cases and controls | `StudyID`, `Age`, `Sex`, `EthSelf`, `HtCurrFt`, `HtCurrIn`, `WtCurrLb`, `MM_STATUS` |
| `Canada_case_control_sample_key.csv` | map from BC study ID to 5hmC sample ID | `Study.ID`, `Assigned.ID` |

If your data are organized differently, adapt the data-loading section at the top of each script (Section 2,
"Data pre-processing", of `03_real_data_analysis.R`).

## Analysis steps

All scripts apply the same cohort construction and preprocessing (manuscript Section 2.1):

1. Cohorts: MMEP patients with a multiple myeloma diagnosis (all cases), and BC cases and controls; participants
   with missing covariates are excluded, and race is restricted to European American (White) and African American.
2. Covariates: age, sex, race, and BMI (≥ 25 kg/m² vs < 25; BMI at diagnosis for MMEP, current BMI for BC). A study
   indicator (S = 1 for MMEP, S = 0 for BC) enters the tree ensemble.
3. 5hmC data: genes with fewer than 10 reads in more than 5% of samples are removed, and counts are normalized with
   the DESeq2 variance-stabilizing transformation.
4. Split: a stratified random split (seed 789) into a development set (n = 500) and a held-out validation set
   (n = 369).
5. Screening (development set only): a probit regression of case status on each gene, adjusted for age and sex;
   genes with p < 0.05 are clustered by the size of their coefficient with a Gaussian mixture model (2–4
   components), and the cluster with the largest coefficients is kept.
6. spBART with 5-fold cross-validation (folds: seed 890): in each fold, a fit with age and sex in the linear term and
   the Bayesian FDR rule at 0.05; the union of the genes selected in any fold; and a final fit on the development set
   with the union genes, the study indicator and all four covariates in the linear term. Discrimination (AUC) and the
   Brier score are reported in cross-validation and on the held-out validation set, together with the posterior
   summaries of the covariate coefficients (manuscript Sections 3.3–3.4 and 4). The average marginal effects of the
   covariates are computed by `06_compare_spBART_BART.R`.

spBART is called through `fit_spBART_analysis(analysis = "application")`: no latent offset and a leaf prior s.d. of
0.5·sd(U₀ − Zβ₀)/√T (about 0.5/√T), the settings used for the application in the paper. The simulation studies use
`fit_spBART()` with its defaults, the `pbart` settings (offset Φ⁻¹(Ȳ), leaf prior s.d. 3/(k√T) with k = 2); see
`R/spBART_sampler.R`.

## Scripts

| Script | What it does | Needs | Writes |
|---|---|---|---|
| `03_real_data_analysis.R` | main spBART analysis (steps 1–6; manuscript Section 4) | `data/` | `results_spBART_<date>/`, `spBART_results/` |
| `07_multichain_stability.R` | the pipeline of `03` with four chains per spBART fit (own seed and initial latent variables; 10,000 burn-in iterations, then 10,000 saved every tenth), convergence check (rank-normalized split-R̂ ≤ 1.02 and bulk ESS ≥ 390) before pooling, and the fold selections and union of each chain (Section 4.2 and Supplementary Section C) | `data/` | `results_spBART_multichain_<date>/` |
| `08_BConly_analysis.R` | the pipeline of `03` on the BC case–control participants alone, without the study indicator; every other step is the same (Sections 4.1–4.2 and Supplementary Section D) | `data/` | `results_spBART_BConly_<date>/` |
| `09_supporting_summaries.R` | DeLong and bootstrap confidence intervals for the validation AUC (pooled, within BC, BC-only), probit coefficients of *IL1RAP* and *ST5* adjusted for age and sex (pooled development set and all BC participants), and their levels among BC cases and controls | latest results of `03` and `08` | `summary/` |
| `00_prepare_analysis_base.R` | builds the analysis data once (split, normalized genes, covariates, screened genes) for `05` and for the semi-synthetic simulations | `data/` | `analysis_base/analysis_base.RDS`, `../simulation/data_semisynthetic/semisynthetic_base.RDS` |
| `04_comparative_analysis_ElasticNet_vs_spBART.R` | elastic net with the same screening, compared with spBART's gene selection across FDR levels | `data/` | `results_comparative_analysis_<date>/` |
| `05_BART_application.R` | BART (`pbart` with the Dirichlet prior) run through the same pipeline, with the covariates in the trees instead of a linear term | output of `00` | `results_BART_application_<date>/` |
| `06_compare_spBART_BART.R` | average marginal effects of the covariates over the validation participants (spBART by posterior g-computation; age: +10 years, others: 1 vs 0), as reported in Table 6 of the paper, and the gene overlap and marginal effects of spBART versus BART | latest results of `03` and `05` | `summary/` |

Run order: `03`, `04`, `07` and `08` are independent of the others. Run `00` before `05`, `06` after `03` and `05`,
and `09` after `03` and `08`.

`07` and `08` run `03` with different options (`spBART.n_chains`, `spBART.cohort`), so all three share one
implementation of the pipeline. `03` and `08` also save the partial dependence of the predicted probability on
*IL1RAP* and *ST5* in the final refit (`partial_dependence.RDS`, `partial_dependence.pdf`), when these genes are in
the union.

## Running the scripts

Install the modified BART package first (from the repository root, `Rscript BART_package/00_install_BART.R`; see the
main [README](../README.md)). Each script sets its working directory to this folder, so it can be called from
anywhere:

```bash
Rscript real_data_analysis/03_real_data_analysis.R
Rscript real_data_analysis/07_multichain_stability.R
Rscript real_data_analysis/08_BConly_analysis.R
Rscript real_data_analysis/09_supporting_summaries.R
Rscript real_data_analysis/00_prepare_analysis_base.R
Rscript real_data_analysis/04_comparative_analysis_ElasticNet_vs_spBART.R
Rscript real_data_analysis/05_BART_application.R
Rscript real_data_analysis/06_compare_spBART_BART.R
```

The analyses use long MCMC chains (2,000 burn-in iterations, then 5,000 iterations thinned by 5, for each fit) and
take several hours each; see the job scripts in `../slurm/` for the memory and time we requested.

On a SLURM cluster, the corresponding job scripts are `real_03_spBART_real_data_analysis.sbatch`,
`real_07_multichain_stability.sbatch` (four cores, one chain per core), `real_08_BConly_analysis.sbatch`,
`00b_prepare_analysis_base.sbatch`, `real_04_comparative_ElasticNet_vs_spBART.sbatch` and
`real_05_BART_application.sbatch` in `../slurm/`; `submit_all.sh` submits them with the right dependencies. The job
files do not name an account, so **you need to provide your own cluster credentials**, for example:

```bash
cd slurm
SPBART_ACCOUNT=<your_account> SPBART_QOS=<your_qos> bash submit_all.sh real_03_spBART_real_data_analysis.sbatch
```

See the main [README](../README.md#running-on-a-slurm-cluster) for the full set of options.

## Keep the outputs private

Every output folder of these scripts (`results_*`, `spBART_results/`, `analysis_base/` and `summary/`) contains
quantities derived from the patient-level data, for example the normalized 5hmC matrix and predicted probabilities
for individual participants. The same applies to `../simulation/data_semisynthetic/`. All of these folders are
excluded from version control. Treat them with the same care as the data themselves, and do not share or commit
them.
