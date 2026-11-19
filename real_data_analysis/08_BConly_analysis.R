################################################################################
# 08_BConly_analysis.R
#
# Analysis restricted to the BC Case-control study (manuscript Sections 4.1-4.2
# and Appendix D).
#
# Runs the complete pipeline of 03_real_data_analysis.R on the BC Case-control
# participants alone. Every step is the same as in the pooled analysis: gene
# filtering and variance-stabilizing normalization, the stratified split into
# 500 development participants and a held-out validation set, screening
# (probit regression of case status on each gene adjusted for age and sex, then
# the Gaussian mixture on the absolute coefficients), 5-fold cross-validation
# with age and sex in the linear part, the Bayesian FDR rule at 0.05 in each
# fold, the union, and the refit with all four covariates. The only differences
# are the data (no MMEP participants) and the absence of the study indicator,
# which is constant within one study.
#
# Usage (from this folder or from slurm/):
#   Rscript 08_BConly_analysis.R
#
# Writes results_spBART_BConly_<date>/ with the same outputs as 03.
################################################################################

local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("--file=", "", grep("--file=", args, value = TRUE)[1])
  here <- if (is.na(f)) getwd() else dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE)))
  assign("spbart_driver_dir", here, envir = globalenv())
})

options(spBART.cohort = "BC_only")

source(file.path(spbart_driver_dir, "03_real_data_analysis.R"))
