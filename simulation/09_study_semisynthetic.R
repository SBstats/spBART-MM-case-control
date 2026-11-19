################################################################################
# 09_study_semisynthetic.R
#
# Semi-synthetic data: the real development-set 5hmC matrix and covariates
# (n = 500) for training and the real validation set (n = 369) for testing;
# only the outcome is simulated, so real gene-gene and gene-covariate
# dependence is kept. 2 x 2 x 2 = 8 scenarios, 200 replicates:
#   gene set : screened (115 genes) | random2000
#   signals  : random (10 genes) | cluster (a gene and its 9 most correlated)
#   f        : trees (DGM 2, empirical-quantile cutpoints) | smooth (DGM 3)
# Methods: spBART, BART, stan4bart, elastic net, knockoffs (second-order).
# Metrics: gene sensitivity / FDR, test AUC / Brier on the real validation set,
# coefficients (spBART, stan4bart), BART's covariate selection (Bayesian FDR on the covariate PIPs),
# and average marginal effects (AMEs) of the 4 covariates on the probability
# scale for spBART, BART (partial dependence) and stan4bart, against their true
# values. Covariate detection (95% interval excludes 0, or BART selection) is the
# main covariate comparison.
# Requires simulation/data_semisynthetic/semisynthetic_base.RDS
# (real_data_analysis/00_prepare_analysis_base.R).
#
# Array: task = (scenario chunk, replicate) = 2 x 200 = 400 tasks.
################################################################################

script_dir <- tryCatch({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE)))
}, error = function(e) getwd())
spbart_root <- dirname(script_dir)
source(file.path(spbart_root, "R", "spBART_sampler.R")); load_spBART_BART()
source(file.path(script_dir, "utils", "data_generation.R"))
source(file.path(script_dir, "utils", "study_functions.R"))

study   <- "09_semisynthetic"
n_reps  <- 200
n_chunk <- 2
mc <- default_mc
methods <- c("spBART", "BART", "stan4bart", "enet", "knockoff")
base <- readRDS(file.path(script_dir, "data_semisynthetic", "semisynthetic_base.RDS"))

grid <- expand.grid(set = c("screened", "random2000"), signals = c("random", "cluster"),
                    f = c("trees", "smooth"), stringsAsFactors = FALSE)
grid$label <- with(grid, paste0("semisyn_", set, "_", signals, "_", f))

id <- task_id()
chunk <- ceiling(id / n_reps); rep <- id - (chunk - 1) * n_reps
todo <- lapply(seq(chunk, nrow(grid), by = n_chunk), function(s) c(s, rep))
out_dir <- study_root(script_dir, study)

for (sr in todo) {
  s <- sr[1]; rep <- sr[2]; g <- grid[s, ]
  f_out <- result_file(out_dir, g$label, rep)
  if (file.exists(f_out)) { cat("Done already:", f_out, "\n"); next }
  seed <- 9e6 + s * 1e4 + rep
  dat <- add_ame_truth(gen_semisyn_data(base, g$set, g$signals, g$f, seed), n_ame = 300)
  res <- run_methods(dat, methods, mc, seed = seed)
  for (m in names(res)) res[[m]]$pip <- NULL
  saveRDS(list(study = study, scenario = as.list(g), rep = rep, seed = seed,
               true_genes = dat$true_genes, ame_truth = dat$ame_truth, seedfix = TRUE, results = res), f_out)
  cat("Saved", f_out, "\n")
}
