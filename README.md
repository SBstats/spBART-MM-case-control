# spBART: semiparametric probit Bayesian additive regression trees

This repository contains the code for the paper

> **An interpretable semiparametric Bayesian machine learning model for high-dimensional circulating epigenetic signatures of multiple myeloma**
> Saurabh Bhandari, Parveen Bhatti, Brian C.-H. Chiu, and Yuan Ji

spBART is a probit model for a binary outcome that combines two parts:

- a **linear term** for a small set of prespecified covariates (age, sex, race and BMI), so that their associations
  are estimated and reported directly, with credible intervals; and
- a **sparse tree ensemble** (BART with the Dirichlet splitting prior of Linero, 2018) for a high-dimensional set of
  features, here gene-level 5-hydroxymethylcytosine (5hmC) signals from circulating cell-free DNA, together with a
  study indicator.

The model is fitted with a Gibbs sampler that cycles through the latent probit variables, the linear coefficients,
the trees and the splitting probabilities. Genes are selected with the Bayesian false discovery rate rule applied to
their posterior inclusion probabilities. The repository contains the real-data application of the paper (a pooled
case–control analysis of cfDNA 5hmC profiles in multiple myeloma) and its simulation studies.

## Data availability

**This repository does not contain any real data.** The application uses patient-level data that cannot be shared,
and none of these data, nor any quantity derived from them, are included here. To run the real-data analysis, you
need your own copies of the data, placed in `real_data_analysis/data/` (see
[`real_data_analysis/README.md`](real_data_analysis/README.md) for the expected files). That folder, and every
folder to which the real-data scripts write their output, is excluded from version control.

The simulation studies generate their own data and run without any real data, with one exception: the semi-synthetic
scenarios (Study 3 of the paper, and one scenario each in Studies 6 and 7) are built from the real 5hmC matrix and
covariates, so they can only be run by users who have the real data.

The repository contains no results. All results are generated when you run the code, and the output folders are
excluded from version control.

## Citation

If you use our code, please cite our paper:

```
Bhandari S, Bhatti P, Chiu BC-H, Ji Y. An interpretable semiparametric Bayesian machine learning model for
high-dimensional circulating epigenetic signatures of multiple myeloma. 2026. Manuscript under review.
```

```bibtex
@unpublished{bhandari2026spbart,
  author = {Bhandari, Saurabh and Bhatti, Parveen and Chiu, Brian C.-H. and Ji, Yuan},
  title  = {An interpretable semiparametric {B}ayesian machine learning model for high-dimensional
            circulating epigenetic signatures of multiple myeloma},
  year   = {2026},
  note   = {Manuscript under review}
}
```

The citation will be updated with the journal reference once the paper is published.

## Repository layout

```
spBART-MM-case-control/
├── BART_package/                   modified BART package used by every analysis
│   ├── BART_2.9.9.1.tar.gz         CRAN BART 2.9.9 plus a stateful sampler for spBART
│   ├── 00_install_BART.R           installs this BART build and every R package used; runs self-tests
│   ├── build_BART_tarball.sh       rebuilds the tarball from CRAN BART 2.9.9 + modified_files/
│   └── modified_files/             the only files that differ from CRAN BART 2.9.9
├── R/spBART_sampler.R              shared spBART sampler, Bayesian FDR gene selection, package check
├── real_data_analysis/             real-data application (Section 4); see real_data_analysis/README.md
├── simulation/                     simulation studies (Section 5 and Supplementary Section B); see simulation/README.md
└── slurm/                          job scripts for a SLURM cluster; logs are written to slurm/logs/
```

## Requirements

- R (version 4.x).
- The modified BART package in `BART_package/` (version 2.9.9.1). The spBART sampler needs its stateful functions;
  the CRAN release of BART does not work.
- CRAN packages: Rcpp, nlme, survival, MASS, truncnorm, caret, pROC, dplyr, readxl, data.table, glmnet, mclust,
  future, future.apply, matrixStats, gtools, ggplot2, stan4bart, corpcor, posterior, BiocManager.
- Bioconductor package: DESeq2 (real-data analysis only).

Install everything, and check the BART build, with:

```bash
Rscript BART_package/00_install_BART.R
```

The installer writes to your personal R library (`R_LIBS_USER` if set, otherwise `~/R/library/<R version>`). It
replaces any BART package already installed there. To keep another BART installation intact, point `R_LIBS_USER` to
a dedicated library first, for example `export R_LIBS_USER=~/R/spBART_library`. Every script that fits BART or
spBART checks at start-up that BART 2.9.9.1 is the version loaded, and stops otherwise.

## Quick start

```bash
# 1. install the modified BART package and all dependencies
Rscript BART_package/00_install_BART.R

# 2. one replicate of a simulation study, with the settings of the paper
Rscript simulation/14_study_confounder_adjustment.R 1

# 3. the real-data analysis (requires your own copies of the data)
Rscript real_data_analysis/03_real_data_analysis.R
```

Folder-specific instructions are in [`simulation/README.md`](simulation/README.md) and
[`real_data_analysis/README.md`](real_data_analysis/README.md).

## Running on a SLURM cluster

The full set of simulations needs many CPU hours, so we ran it on a SLURM cluster. Every job script is in `slurm/`;
submit from that folder. The job scripts are templates: they request time, memory and cores only, and name no
account, partition, user or e-mail address.

**You need to provide your own cluster credentials.** Pass them through environment variables on every submission:

```bash
cd slurm
SPBART_ACCOUNT=<your_account> SPBART_QOS=<your_qos> bash submit_all.sh --dry-run   # show the plan; submit nothing
SPBART_ACCOUNT=<your_account> SPBART_QOS=<your_qos> bash submit_all.sh             # submit every job
```

- `SPBART_ACCOUNT` (required): the SLURM account to charge.
- `SPBART_QOS` (optional): the quality of service to submit under. If your cluster does not use one, leave it out
  and set `SPBART_MAX_JOBS=<n>` to your limit on queued jobs instead.
- `SPBART_MAX_JOBS` (optional with a QOS): overrides the queued-job limit that is otherwise read from the QOS.

`submit_all.sh` submits the package installation first, then the real-data preparation, then all simulation arrays
and real-data jobs; each job starts only after the jobs it depends on finish successfully. The header of
`submit_all.sh` shows the full dependency tree. The script sizes the job arrays to fit your queue limit and creates
the log folders under `slurm/logs/` that the jobs write to. Adapt the `module load R` line in the job files to the R
module available on your cluster.

If compute nodes cannot reach CRAN, run `Rscript ../BART_package/00_install_BART.R` on a login or development node
first, and then submit the remaining jobs by name:

```bash
SPBART_ACCOUNT=<your_account> SPBART_QOS=<your_qos> bash submit_all.sh <job1>.sbatch <job2>.sbatch ...
```

To submit a single job file directly with `sbatch` instead, give your credentials on the command line and create its
log folder first (see the `--output` line of the job file), for example:

```bash
mkdir -p logs/14_confounder_adjustment
sbatch --account=<your_account> --qos=<your_qos> sim_14_study_confounder_adjustment.sbatch
```

To find and rerun missing simulation results (credentials as above):

```bash
cd slurm
Rscript check_missing.R                                  # lists missing results and the job files to rerun
SPBART_ACCOUNT=<your_account> SPBART_QOS=<your_qos> bash rerun_studies.sh --dry-run <job>.sbatch ...   # preview
SPBART_ACCOUNT=<your_account> SPBART_QOS=<your_qos> bash rerun_studies.sh <job>.sbatch ...             # fill gaps only
SPBART_ACCOUNT=<your_account> SPBART_QOS=<your_qos> bash rerun_studies.sh --fresh <job>.sbatch ...     # rerun every replicate
```

`--fresh` moves the old results to `simulation/results_archive/<timestamp>/`; nothing is deleted. Without it,
finished replicates are skipped. `bash status_report.sh` summarizes the state of all submitted jobs.

## Contact

Please open an issue in this repository for questions about the code.
