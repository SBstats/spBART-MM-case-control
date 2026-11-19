################################################################################
# 06_compare_spBART_BART.R
#
# Compares the real-data fits of 03_real_data_analysis.R (spBART) and
# 05_BART_application.R (BART) from their latest results folders:
#   - overlap of the selected genes (Jaccard index), BART validation AUC
#   - average marginal effects (risk differences) of each covariate over the
#     validation participants, spBART by posterior g-computation
#     (age: +10 years; others: 1 vs 0)
# Writes real_data_analysis/summary/real_data_gene_overlap.csv and
# real_data_analysis/summary/real_data_average_marginal_effects.csv.
#
# Usage (after real_03 and real_05 have run):  Rscript real_data_analysis/06_compare_spBART_BART.R
################################################################################
rd <- tryCatch({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE)))
}, error = function(e) getwd())
out <- file.path(rd, "summary"); dir.create(out, showWarnings = FALSE)
save_csv <- function(df, name) {
  write.csv(df, file.path(out, name), row.names = FALSE)
  cat("wrote", name, "(", nrow(df), "rows )\n")
}
latest <- function(pattern) { d <- sort(list.files(rd, pattern = pattern, full.names = TRUE)); if (length(d)) tail(d, 1) else NA }
sp_dir <- latest("^results_spBART_[0-9]")
bt_dir <- latest("^results_BART_application_[0-9]")
if (is.na(sp_dir) || !file.exists(file.path(sp_dir, "final_model_draws.RDS")))
  stop("No spBART results with final_model_draws.RDS; run 03_real_data_analysis.R first")

# Validation participants: f at their predictors from the refit, and their covariates
fd <- readRDS(file.path(sp_dir, "final_model_draws.RDS"))
if (is.null(fd$f_test_draws) || is.null(fd$Z_test))
  stop("final_model_draws.RDS has no validation-set draws; rerun 03_real_data_analysis.R")
Z <- fd$Z_test
stopifnot(ncol(fd$f_test_draws) == nrow(Z))
ame <- do.call(rbind, lapply(colnames(Z), function(v) {
  Z1 <- Z; Z0 <- Z
  if (v == "age") Z1[, v] <- Z[, v] + 10 else { Z1[, v] <- 1; Z0[, v] <- 0 }
  eta_f <- fd$offset + fd$f_test_draws
  d <- rowMeans(pnorm(eta_f + fd$beta_draws %*% t(Z1)) - pnorm(eta_f + fd$beta_draws %*% t(Z0)))
  data.frame(covariate = v, contrast = if (v == "age") "+10 years" else "1 vs 0", AME = mean(d),
             lower = unname(quantile(d, .025)), upper = unname(quantile(d, .975)), P_positive = mean(d > 0))
}))
ame$method <- "spBART"
if (!is.na(bt_dir)) {
  bt <- readRDS(file.path(bt_dir, "BART_application_results.RDS"))
  ame <- rbind(ame, transform(bt$ame, method = "BART"))
  sp_genes <- setdiff(readRDS(file.path(sp_dir, "selected_gene_names.RDS")), "study_indicator")
  bt_genes <- setdiff(bt$final_selected, "study_indicator")
  save_csv(data.frame(spBART_genes = paste(sp_genes, collapse = ";"), BART_genes = paste(bt_genes, collapse = ";"),
                      n_spBART = length(sp_genes), n_BART = length(bt_genes), n_both = length(intersect(sp_genes, bt_genes)),
                      jaccard = length(intersect(sp_genes, bt_genes)) / max(1, length(union(sp_genes, bt_genes))),
                      BART_validation_AUC = bt$validation["AUC"]),
           "real_data_gene_overlap.csv")
} else cat("No BART application results yet (05_BART_application.R); writing spBART AMEs only\n")
save_csv(ame, "real_data_average_marginal_effects.csv")
