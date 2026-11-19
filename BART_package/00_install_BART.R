################################################################################
# 00_install_BART.R
#
# Install the modified BART package (version 2.9.9.1) used by ALL spBART
# analyses (simulation and real data), then run three self-tests.
#
# Run ONCE per machine / cluster account, before any analysis:
#   Rscript BART_package/00_install_BART.R
# or on the cluster:
#   cd slurm && sbatch 00_install_BART.sbatch
#
# Never install from inside array jobs: concurrent installs into the same
# library can corrupt it.
################################################################################

script_dir <- tryCatch({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE)))  # Rscript encodes spaces as ~+~
  } else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    getwd()
  }
}, error = function(e) getwd())
spbart_root <- dirname(script_dir)
source(file.path(spbart_root, "R", "spBART_sampler.R"))   # SPBART_BART_VERSION, spbart_user_lib()

tarball <- file.path(script_dir, paste0("BART_", SPBART_BART_VERSION, ".tar.gz"))
if (!file.exists(tarball)) stop("Tarball not found: ", tarball)

lib <- spbart_user_lib()
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(lib, .libPaths())))
cat("Installing", basename(tarball), "into", lib, "\n")

# Every package the simulation and real-data scripts load, installed here ONCE so
# that array tasks never install anything (concurrent installs corrupt libraries)
cran_pkgs <- c("Rcpp", "nlme", "survival",                       # BART
               "MASS", "truncnorm", "caret", "pROC",             # simulations
               "dplyr", "readxl", "data.table", "glmnet", "mclust",
               "future", "future.apply", "matrixStats", "gtools", "ggplot2",  # real data
               "stan4bart", "corpcor", "posterior",             # simulation studies 07-15
               "BiocManager")
missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  cat("Installing CRAN packages:", paste(missing_cran, collapse = ", "), "\n")
  install.packages(missing_cran, lib = lib, repos = "https://cloud.r-project.org")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  cat("Installing Bioconductor package DESeq2\n")
  BiocManager::install("DESeq2", lib = lib, ask = FALSE, update = FALSE)
}
still_missing <- c(cran_pkgs, "DESeq2")[!vapply(c(cran_pkgs, "DESeq2"), requireNamespace,
                                                logical(1), quietly = TRUE)]
if (length(still_missing) > 0) stop("Could not install: ", paste(still_missing, collapse = ", "))
install.packages(tarball, lib = lib, repos = NULL, type = "source",
                 INSTALL_opts = "--no-multiarch")

#-------------------------------------------------------------------------------
# Self-tests (expected results in brackets)
#-------------------------------------------------------------------------------
load_spBART_BART()
cat("Loaded BART", as.character(packageVersion("BART")), "from",
    dirname(find.package("BART")), "\n")

set.seed(1)
n <- 600; p <- 200
X <- matrix(rnorm(n * p), n, p)
y <- as.numeric(2 * (X[, 1] > 0) + 2 * (X[, 2] > 0) - 2 + rnorm(n))
s <- wbart_create(X, y, sparse = TRUE, a = 0.5, b = 1, rho = p, ntree = 50)
wbart_start_dart(s)
for (it in 1:1500) r <- wbart_run_iteration(s)
vc <- replicate(200, wbart_run_iteration(s)$varcount)
null_pip <- mean(rowMeans(vc > 0)[-(1:2)])
wbart_destroy(s)

set.seed(2); s2 <- wbart_create(X, y, ntree = 20)
invisible(wbart_run_iteration(s2)); u_after <- runif(3); wbart_destroy(s2)
set.seed(2); u_noiter <- runif(3)

checks <- c(
  "DART active (mean null-variable PIP < 0.10)"  = null_pip < 0.10,
  "Error s.d. fixed at 1"                        = isTRUE(all.equal(r$sigma, 1)),
  "C++ draws advance R's RNG stream"             = !identical(u_after, u_noiter)
)
print(data.frame(check = names(checks), passed = unname(checks)), row.names = FALSE)
cat(sprintf("Mean null-variable PIP = %.3f\n", null_pip))
if (!all(checks)) stop("Self-test failed: do not run analyses with this build.")
cat("\nModified BART", SPBART_BART_VERSION, "installed and verified.\n")
