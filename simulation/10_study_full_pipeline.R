################################################################################
# 10_study_full_pipeline.R
#
# The complete procedure used in the application: probit + GMM screening
# (adjusted for age and sex), then 5-fold CV with Bayesian FDR (0.05) per fold,
# union of fold selections, final fit on the union, test-set prediction.
# Run for spBART and for BART on the same screened genes.
# Scenarios: DGM 1 and DGM 2 at n = 1500, p in {2000, 5000}; 100 replicates.
# Reports screening recall (true genes surviving screening), then selection
# and prediction of each pipeline.
#
# Array: task = (scenario, replicate) = 4 x 100 = 400 tasks.
################################################################################

script_dir <- tryCatch({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE)))
}, error = function(e) getwd())
spbart_root <- dirname(script_dir)
source(file.path(spbart_root, "R", "spBART_sampler.R")); load_spBART_BART()
source(file.path(script_dir, "utils", "data_generation.R"))
source(file.path(script_dir, "utils", "study_functions.R"))

study  <- "10_full_pipeline"
n_reps <- 100
mc <- default_mc
nn <- 1500
pp <- c(2000, 5000)
scenarios <- list(
  scenario(paste0("DGM2_n1500_p", pp[1]), nn, pp[1]),
  scenario(paste0("DGM2_n1500_p", pp[2]), nn, pp[2]),
  scenario(paste0("DGM1_n1500_p", pp[1]), nn, pp[1], type = "dgm1"),
  scenario(paste0("DGM1_n1500_p", pp[2]), nn, pp[2], type = "dgm1"))

id <- task_id()
todo <- list(c(ceiling(id / n_reps), id - (ceiling(id / n_reps) - 1) * n_reps))
out_dir <- study_root(script_dir, study)

for (sr in todo) {
  s <- sr[1]; rep <- sr[2]; sc <- scenarios[[s]]
  f_out <- result_file(out_dir, sc$label, rep)
  if (file.exists(f_out)) { cat("Done already:", f_out, "\n"); next }
  seed <- 10e6 + s * 1e4 + rep
  dat <- gen_study_data(sc, seed, n_test = 1000)
  t0 <- Sys.time()
  screened <- screen_genes(dat$X, dat$Y, dat$Z[, c("age", "sex")])
  screening <- list(n_screened = length(screened), recall = mean(dat$true_genes %in% screened),
                    minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  res <- list()
  if (length(screened) > 0) {
    for (m in c("spBART", "BART")) {
      res[[m]] <- tryCatch(cv_union_pipeline(m, dat, screened, mc, seed = seed),
                           error = function(e) list(method = m, error = conditionMessage(e)))
    }
  }
  saveRDS(list(study = study, scenario = sc, rep = rep, seed = seed, screening = screening,
               screened = screened, results = res), f_out)
  cat("Saved", f_out, "\n")
}
