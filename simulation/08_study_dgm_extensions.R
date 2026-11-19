################################################################################
# 08_study_dgm_extensions.R
#
# Further data-generating models at n = 1500, p in {500, 2000}, 100 replicates:
#   DGM2                    baseline (independent genes, tree signal)
#   smooth                  DGM 3: smooth nonlinear f
#   corrX_rho{.5,.8}_{sep,same}  block-correlated genes (blocks of 50); signal
#                           genes in separate blocks or 5 + 5 in two blocks
#   corrZ_c{0,.3,.6,.9}_{trees,linear}  covariates correlated with signal genes,
#                           tree or near-linear f (identifiability);
#                           c = 0 with trees is the DGM2 baseline
#   partial                 sex and race linear; age and BMI nonlinear and
#                           interacting with signal genes
# Methods (one full-data fit each, test set n = 1000): spBART, BART,
# stan4bart, elastic net, knockoffs; plus spBART with age and BMI moved to the
# trees in the "partial" scenario.
#
# Array: task = (p, scenario chunk, replicate) = 2 x 4 x 100 = 800 tasks.
################################################################################

script_dir <- tryCatch({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE)))
}, error = function(e) getwd())
spbart_root <- dirname(script_dir)
source(file.path(spbart_root, "R", "spBART_sampler.R")); load_spBART_BART()
source(file.path(script_dir, "utils", "data_generation.R"))
source(file.path(script_dir, "utils", "study_functions.R"))

study   <- "08_dgm_extensions"
n_reps  <- 100
n_chunk <- 4
mc <- default_mc
methods <- c("spBART", "BART", "stan4bart", "enet", "knockoff")

make_scenarios <- function(p, n) {
  L <- function(x) paste0(x, "_n", n, "_p", p)
  c(list(scenario(L("DGM2"), n, p),
         scenario(L("smooth"), n, p, f = "smooth")),
    unlist(lapply(c(0.5, 0.8), function(r) lapply(c("sep", "same"), function(lay)
      scenario(L(paste0("corrX_rho", r, "_", lay)), n, p, rho = r, layout = lay))), recursive = FALSE),
    list(scenario(L("corrZ_c0_linear"), n, p, f = "linear")),
    unlist(lapply(c(0.3, 0.6, 0.9), function(cz) lapply(c("trees", "linear"), function(ff)
      scenario(L(paste0("corrZ_c", cz, "_", ff)), n, p, corrZ = cz, f = ff))), recursive = FALSE),
    list(scenario(L("partial"), n, p, type = "partial")))
}
p_values <- c(500, 2000)       # corrX layouts need p >= 460
n_train  <- 1500

id <- task_id()
per_p <- n_chunk * n_reps
p <- p_values[ceiling(id / per_p)]
rem <- id - (ceiling(id / per_p) - 1) * per_p
chunk <- ceiling(rem / n_reps); rep <- rem - (chunk - 1) * n_reps
scs <- make_scenarios(p, n_train)
todo <- lapply(scs[seq(chunk, length(scs), by = n_chunk)], function(sc) list(sc = sc, rep = rep))
out_dir <- study_root(script_dir, study)

for (t in todo) {
  sc <- t$sc; rep <- t$rep
  f_out <- result_file(out_dir, sc$label, rep)
  if (file.exists(f_out)) { cat("Done already:", f_out, "\n"); next }
  seed <- 8e6 + (match(sc$label, sapply(make_scenarios(sc$p, sc$n), `[[`, "label")) + 100 * (sc$p == 2000)) * 1e4 + rep
  dat <- gen_study_data(sc, seed, n_test = 1000)
  ms <- c(methods, if (sc$type == "partial") "spBART_partialZ")
  res <- run_methods(dat, ms, mc, seed = seed)
  for (m in names(res)) res[[m]]$pip <- NULL
  saveRDS(list(study = study, scenario = sc, rep = rep, seed = seed, results = res), f_out)
  cat("Saved", f_out, "\n")
}
