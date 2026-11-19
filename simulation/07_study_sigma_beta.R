################################################################################
# 07_study_sigma_beta.R
#
# Sensitivity of spBART to the prior variance sigma_beta^2 in {0.1, 1, 10,
# 100, 1000}. The real-data analysis uses sigma_beta^2 = 10.
# Scenarios (100 replicates each):
#   DGM2_n1500_p500, DGM2_n1500_p2000   (simulated, independent genes)
#   semisyn_app_n500_p115               (real development-set X and Z)
# Per sigma_beta^2: beta bias / coverage / width, variance shares of Z'beta and
# f, test AUC/Brier of Z'beta alone, f alone and the full model, selection.
#
# Array: task = (scenario, replicate), 3 x 100 = 300 tasks.
################################################################################

script_dir <- tryCatch({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE)))
}, error = function(e) getwd())
spbart_root <- dirname(script_dir)
source(file.path(spbart_root, "R", "spBART_sampler.R")); load_spBART_BART()
source(file.path(script_dir, "utils", "data_generation.R"))
source(file.path(script_dir, "utils", "study_functions.R"))

study  <- "07_sigma_beta"
grid   <- c(0.1, 1, 10, 100, 1000)
n_reps <- 100
mc     <- default_mc
nn     <- 1500
scenarios <- list(
  list(label = "DGM2_n1500_p500",  sc = scenario("DGM2_n1500_p500",  nn, 500)),
  list(label = "DGM2_n1500_p2000", sc = scenario("DGM2_n1500_p2000", nn, 2000)),
  list(label = "semisyn_app_n500_p115", semisyn = list(set = "screened", signals = "random", f = "trees")))

id <- task_id()
todo <- list(c(ceiling(id / n_reps), id - (ceiling(id / n_reps) - 1) * n_reps))
out_dir <- study_root(script_dir, study)
semisyn_base <- NULL

for (sr in todo) {
  s <- sr[1]; rep <- sr[2]; S <- scenarios[[s]]
  f_out <- result_file(out_dir, S$label, rep)
  if (file.exists(f_out)) { cat("Done already:", f_out, "\n"); next }
  seed <- 7e6 + s * 1e4 + rep
  dat <- if (is.null(S$semisyn)) {
    gen_study_data(S$sc, seed, n_test = 1000)
  } else {
    if (is.null(semisyn_base))
      semisyn_base <- readRDS(file.path(script_dir, "data_semisynthetic", "semisynthetic_base.RDS"))
    gen_semisyn_data(semisyn_base, S$semisyn$set, S$semisyn$signals, S$semisyn$f, seed)
  }
  res <- lapply(grid, function(g) {
    r <- tryCatch(m_spbart(dat, mc, sigma2_beta = g), error = function(e) list(error = conditionMessage(e)))
    r$pip <- NULL
    r
  })
  names(res) <- paste0("sigma2_", grid)
  saveRDS(list(study = study, scenario = S$label, rep = rep, seed = seed, grid = grid, results = res), f_out)
  cat("Saved", f_out, "\n")
}
