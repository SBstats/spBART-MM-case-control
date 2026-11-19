################################################################################
# 14_study_confounder_adjustment.R  -- a realistic adjustment set and null
# genes that track a strong confounder
#
# 8 linear covariates (age, sex, race, BMI, smoking, education, activity, family
# history) with small effects (0.10-0.40) and a strong age effect beta_age;
# 10 signal genes (DGM 2 trees, independent of the covariates); 30 NULL genes
# correlated with age (rho_proxy) that have no effect on the outcome given age.
#   beta_age  in {0.5, 1.0, 1.5}, rho_proxy in {0.5, 0.8}, p in {2000, 5000},
#   n = 1500, 100 replicates (12 scenarios). rho = 0.8 is far above the
#   gene-age correlations in the real data.
#   Two further scenarios let age act as a U-shape (beta_age = 1.0,
#   rho_proxy = 0.5), which a linear term cannot represent.
# Mechanism: BART builds the age effect from tree splits and can pick up the
# leftover age signal through the proxy genes (false discoveries); spBART's
# linear term absorbs age, so the proxies carry no information given age;
# stan4bart adjusts linearly but has no sparsity prior.
# Methods (one full-data fit each; test set n = 1000): spBART, BART, stan4bart,
# elastic net. Metrics: gene sensitivity / FDR, proxy genes selected, test AUC /
# Brier, coefficients (spBART, stan4bart), and average marginal effects of the
# 8 covariates on the probability scale for every method (BART by partial
# dependence), against their true values.
#
# Array: unit = (scenario, replicate) = 14 x 100 = 1400 units.
################################################################################

script_dir <- tryCatch({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE)))
}, error = function(e) getwd())
spbart_root <- dirname(script_dir)
source(file.path(spbart_root, "R", "spBART_sampler.R")); load_spBART_BART()
source(file.path(script_dir, "utils", "data_generation.R"))
source(file.path(script_dir, "utils", "study_functions.R"))

study   <- "14_confounder_adjustment"
n_reps  <- 100
mc <- default_mc
methods <- c("spBART", "BART", "stan4bart", "enet")
nn <- 1500
pp <- c(2000, 5000)

scenarios <- list()
for (p in pp) {
  for (ba in c(0.5, 1.0, 1.5)) for (rp in c(0.5, 0.8))
    scenarios[[length(scenarios) + 1]] <- conf_scenario(sprintf("conf_age%.1f_rho%.1f_n1500_p%d", ba, rp, p), nn, p, ba, rp)
  scenarios[[length(scenarios) + 1]] <- conf_scenario(sprintf("conf_ushape_age1.0_rho0.5_n1500_p%d", p), nn, p, 1.0, 0.5,
                                                      age_form = "ushape")
}

id <- task_id()
todo <- list(c(ceiling(id / n_reps), id - (ceiling(id / n_reps) - 1) * n_reps))
out_dir <- study_root(script_dir, study)

for (sr in todo) {
  s <- sr[1]; rep <- sr[2]; sc <- scenarios[[s]]
  if (is.null(sc)) { cat("No scenario for unit", id, "\n"); next }
  f_out <- result_file(out_dir, sc$label, rep)
  if (file.exists(f_out)) { cat("Done already:", f_out, "\n"); next }
  seed <- 14e6 + s * 1e4 + rep
  dat <- gen_confounder_data(sc, seed, n_test = 1000, n_ame = 300)
  res <- run_methods(dat, methods, mc, seed = seed)
  for (m in names(res)) res[[m]]$pip <- NULL
  saveRDS(list(study = study, scenario = sc, rep = rep, seed = seed, ame_truth = dat$ame_truth,
               seedfix = TRUE, results = res), f_out)
  cat("Saved", f_out, "\n")
}
