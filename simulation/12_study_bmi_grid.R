################################################################################
# 12_study_bmi_grid.R
#
# DGM 2 with the BMI coefficient varied, everything else fixed:
#   beta_BMI in {0, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8}; n = 1500, p in {500, 2000};
#   100 replicates. beta_BMI = 0 gives each method's false-positive rate.
# BART: covariate selected by Bayesian FDR on the four covariate PIPs.
# spBART: 95% credible interval for beta_BMI excludes 0.
# One full-data fit per method; test set n = 1000.
#
# Array: task = (p, replicate) = 2 x 100 = 200 tasks, each running the 7 values.
################################################################################

script_dir <- tryCatch({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE)))
}, error = function(e) getwd())
spbart_root <- dirname(script_dir)
source(file.path(spbart_root, "R", "spBART_sampler.R")); load_spBART_BART()
source(file.path(script_dir, "utils", "data_generation.R"))
source(file.path(script_dir, "utils", "study_functions.R"))

study   <- "12_bmi_grid"
n_reps  <- 100
bmi_grid <- c(0, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8)
mc <- default_mc
p_values <- c(500, 2000)
nn <- 1500

id <- task_id()
todo <- list(c(ceiling(id / n_reps), id - (ceiling(id / n_reps) - 1) * n_reps))
out_dir <- study_root(script_dir, study)

for (pr in todo) {
  p <- p_values[pr[1]]; rep <- pr[2]
  for (b in bmi_grid) {
    sc <- scenario(sprintf("DGM2_bmi%.2f_n1500_p%d", b, p), nn, p, beta_bmi = b)
    f_out <- result_file(out_dir, sc$label, rep)
    if (file.exists(f_out)) { cat("Done already:", f_out, "\n"); next }
    seed <- 12e6 + (match(b, bmi_grid) + 10 * pr[1]) * 1e4 + rep
    dat <- gen_study_data(sc, seed, n_test = 1000)
    res <- run_methods(dat, c("spBART", "BART"), mc, seed = seed)
    for (m in names(res)) res[[m]]$pip <- NULL
    saveRDS(list(study = study, scenario = sc, beta_bmi = b, rep = rep, seed = seed, results = res), f_out)
    cat("Saved", f_out, "\n")
  }
}
