################################################################################
# 15_study_ultrahigh_dimension.R  -- very high dimension without screening
#
# DGM 2 (10 signal genes, gene trees + linear effects of age, sex, race, BMI),
# p in {5000, 10000}, n in {1000, 1500}, no screening step, 100 replicates
# (4 scenarios).
# Methods (one full-data fit each; test set n = 1000): spBART, BART, stan4bart,
# elastic net. Metrics: gene sensitivity / FDR, test AUC / Brier, coefficients,
# and average marginal effects of the 4 covariates (BART by partial dependence).
#
# Array: unit = (scenario, replicate) = 4 x 100 = 400 units.
################################################################################

script_dir <- tryCatch({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE)))
}, error = function(e) getwd())
spbart_root <- dirname(script_dir)
source(file.path(spbart_root, "R", "spBART_sampler.R")); load_spBART_BART()
source(file.path(script_dir, "utils", "data_generation.R"))
source(file.path(script_dir, "utils", "study_functions.R"))

study   <- "15_ultrahigh_dimension"
n_reps  <- 100
mc <- default_mc
methods <- c("spBART", "BART", "stan4bart", "enet")
grid <- expand.grid(n = c(1000, 1500), p = c(5000, 10000))
scenarios <- lapply(seq_len(nrow(grid)), function(i)
  scenario(sprintf("DGM2_n%d_p%d", grid$n[i], grid$p[i]), grid$n[i], grid$p[i]))

id <- task_id()
todo <- list(c(ceiling(id / n_reps), id - (ceiling(id / n_reps) - 1) * n_reps))
out_dir <- study_root(script_dir, study)

for (sr in todo) {
  s <- sr[1]; rep <- sr[2]; sc <- scenarios[[s]]
  if (is.null(sc)) { cat("No scenario for unit", id, "\n"); next }
  f_out <- result_file(out_dir, sc$label, rep)
  if (file.exists(f_out)) { cat("Done already:", f_out, "\n"); next }
  seed <- 15e6 + s * 1e4 + rep
  dat <- add_ame_truth(gen_study_data(sc, seed, n_test = 1000), n_ame = 300)
  res <- run_methods(dat, methods, mc, seed = seed)
  for (m in names(res)) res[[m]]$pip <- NULL
  saveRDS(list(study = study, scenario = sc, rep = rep, seed = seed, ame_truth = dat$ame_truth,
               seedfix = TRUE, results = res), f_out)
  cat("Saved", f_out, "\n")
}
