################################################################################
# check_missing.R -- which simulation results are missing, per study and setting.
# Compares simulation/results/ with the replicates every study should
# produce, prints the gaps, and lists the job files to rerun.
#
# Usage (from slurm/):  module load R && Rscript check_missing.R
################################################################################
args_all <- commandArgs(trailingOnly = FALSE)
here <- tryCatch(dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "",
  grep("^--file=", args_all, value = TRUE)[1]), fixed = TRUE))), error = function(e) getwd())
res <- file.path(dirname(here), "simulation", "results")
if (!dir.exists(res)) stop("No results folder at ", res)

have <- function(dir) {
  f <- list.files(dir, pattern = "^rep_[0-9]+\\.RDS$")
  as.integer(sub("rep_([0-9]+)\\.RDS", "\\1", f))
}
show <- function(miss) {
  if (!length(miss)) return("")
  if (length(miss) > 30) paste(paste(head(miss, 30), collapse = " "), "...") else paste(miss, collapse = " ")
}
report <- list()
check <- function(job, dir, reps) {
  miss <- setdiff(seq_len(reps), have(dir))
  report[[length(report) + 1]] <<- data.frame(job = job, setting = sub(paste0(res, "/"), "", dir, fixed = TRUE),
                                              expected = reps, missing = length(miss), ids = show(miss))
}

# Paper simulations (01-04, 06): results/TrueModel{1,2}/n{n}/p{p}/{method}/rep_###.RDS, 500 replicates
main <- list(sim_01_TrueModel1_regular_BART = c("TrueModel1", "regular_BART"),
             sim_02_TrueModel1_spBART       = c("TrueModel1", "spBART"),
             sim_03_TrueModel2_regular_BART = c("TrueModel2", "regular_BART"),
             sim_04_TrueModel2_spBART       = c("TrueModel2", "spBART"),
             sim_06_TrueModel2_logistic_regression = c("TrueModel2", "logistic_reg"),
             sim_06_TrueModel2_logistic_regression = c("TrueModel2", "probit_reg"))
for (j in seq_along(main)) for (n in c(1500, 2500)) for (p in c(500, 2000, 3000))
  check(names(main)[j], file.path(res, main[[j]][1], paste0("n", n), paste0("p", p), main[[j]][2]), 500)

# Studies 07-15: results/<study>/<scenario>/rep_###.RDS
studies <- list(
  sim_07_study_sigma_beta     = list(dir = "07_sigma_beta", reps = 100, n_scen = 3),
  sim_08_study_dgm_extensions = list(dir = "08_dgm_extensions", reps = 100, n_scen = 28),
  sim_09_study_semisynthetic  = list(dir = "09_semisynthetic", reps = 200, n_scen = 8),
  sim_10_study_full_pipeline  = list(dir = "10_full_pipeline", reps = 100, n_scen = 4),
  sim_11_study_multichain     = list(dir = "11_multichain", reps = 100, n_scen = 6),
  sim_12_study_bmi_grid       = list(dir = "12_bmi_grid", reps = 100, n_scen = 14),
  sim_14_study_confounder_adjustment = list(dir = "14_confounder_adjustment", reps = 100, n_scen = 14),
  sim_15_study_ultrahigh_dimension   = list(dir = "15_ultrahigh_dimension", reps = 100, n_scen = 4))
for (j in names(studies)) {
  s <- studies[[j]]; base <- file.path(res, s$dir)
  scen <- if (dir.exists(base)) list.dirs(base, recursive = FALSE) else character(0)
  if (length(scen) < s$n_scen)
    report[[length(report) + 1]] <- data.frame(job = j, setting = paste0(s$dir, " (scenario folders)"), expected = s$n_scen,
                                               missing = s$n_scen - length(scen), ids = "whole scenarios missing")
  for (d in scen) check(j, d, s$reps)
}

r <- do.call(rbind, report)
gaps <- r[r$missing > 0, ]
cat(sprintf("Checked %d settings in %s\n", nrow(r), res))
if (!nrow(gaps)) { cat("Nothing is missing.\n"); quit(status = 0) }
cat(sprintf("\n%d settings have missing results:\n\n", nrow(gaps)))
for (j in unique(gaps$job)) {
  g <- gaps[gaps$job == j, ]
  cat(sprintf("== %s.sbatch: %d missing results in %d settings\n", j, sum(g$missing), nrow(g)))
  for (i in seq_len(nrow(g))) cat(sprintf("   %-48s %3d of %3d missing: %s\n", g$setting[i], g$missing[i], g$expected[i], g$ids[i]))
}
jobs <- paste0(unique(gaps$job), ".sbatch")
cat("\nTo fill only the gaps (finished results are kept and skipped):\n   bash rerun_studies.sh", paste(jobs, collapse = " "), "\n")
cat("To rerun these studies from scratch (old results moved to simulation/results_archive/):\n   bash rerun_studies.sh --fresh", paste(jobs, collapse = " "), "\n")
