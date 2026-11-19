################################################################################
# 05_BART_application.R
#
# BART with the Dirichlet prior (Linero 2018; pbart, sparse = TRUE) applied to
# the pooled MM data with exactly the spBART pipeline of 03_real_data_analysis.R:
# same split, same 115 screened genes, same folds (seed 890), Bayesian FDR 0.05
# per fold, union of fold selections (study indicator always kept), final fit
# on the development set, validation on the held-out set. The difference is
# that age and sex (CV) and all four covariates (final model) enter the trees
# instead of a linear term.
# Also reports covariate PIPs / FDR selection and average marginal effects
# (risk differences) of each covariate over the validation participants,
# for comparison with spBART.
#
# Requires real_data_analysis/analysis_base/analysis_base.RDS
# (00_prepare_analysis_base.R).
################################################################################

script_path <- tryCatch({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE)))
}, error = function(e) getwd())
setwd(script_path)
spbart_root <- dirname(script_path)
source(file.path(spbart_root, "R", "spBART_sampler.R")); load_spBART_BART()
suppressPackageStartupMessages({ library(caret); library(pROC) })

nskip <- 2000
ndpost <- 1000
keepevery <- 5
out_dir <- file.path(script_path, paste0("results_BART_application_",
                                         format(Sys.Date(), "%Y_%m_%d")))
dir.create(out_dir, showWarnings = FALSE)

b <- readRDS(file.path("analysis_base", "analysis_base.RDS"))
genes <- b$screened_genes
X_cv   <- cbind(b$X_genes_cv[, genes],   study_indicator = b$study_cv)
X_test <- cbind(b$X_genes_test[, genes], study_indicator = b$study_test)
Y_cv <- b$Y_cv; Y_test <- b$Y_test
cat("Development", nrow(X_cv), "x", ncol(X_cv), "; validation", nrow(X_test), "\n")

auc_of <- function(y, p) as.numeric(auc(roc(y, p, quiet = TRUE, levels = c(0, 1), direction = "<")))
fit_pbart <- function(x, y, xt) {
  pbart(x.train = x, y.train = y, x.test = xt, ntree = 200, ndpost = ndpost, nskip = nskip,
        keepevery = keepevery, sparse = TRUE, a = 0.5, b = 1, rho = ncol(x), printevery = 1e6)
}

#-------------------------------------------------------------------------------
# 5-fold CV on the development set (folds identical to 03: seed 890)
#-------------------------------------------------------------------------------
set.seed(890)
folds <- createFolds(Y_cv, k = 5, list = TRUE, returnTrain = FALSE)
fold_sel <- list(); fold_perf <- list(); fold_cov_pip <- list()
for (k in seq_along(folds)) {
  val <- folds[[k]]; tr <- setdiff(seq_along(Y_cv), val)
  x  <- cbind(X_cv[tr, ],  b$Z_cv[tr, ]);  xt <- cbind(X_cv[val, ], b$Z_cv[val, ])
  fit <- fit_pbart(x, Y_cv[tr], xt)
  colnames(fit$varcount) <- colnames(x)
  pip <- colMeans(fit$varcount > 0)
  fold_sel[[k]] <- fdr_gene_selection(pip[colnames(X_cv)], 0.05)$selected_genes
  fold_cov_pip[[k]] <- pip[colnames(b$Z_cv)]
  aucs <- apply(fit$prob.test, 1, function(pr) auc_of(Y_cv[val], pr))
  briers <- apply(fit$prob.test, 1, function(pr) mean((Y_cv[val] - pr)^2))
  fold_perf[[k]] <- c(AUC = mean(aucs), AUC_lo = unname(quantile(aucs, .025)), AUC_hi = unname(quantile(aucs, .975)),
                      Brier = mean(briers), Brier_lo = unname(quantile(briers, .025)), Brier_hi = unname(quantile(briers, .975)))
  cat(sprintf("Fold %d: %d selected, AUC %.3f\n", k, length(fold_sel[[k]]), fold_perf[[k]]["AUC"]))
}
union <- unique(unlist(fold_sel))
if (!"study_indicator" %in% union) union <- c(union, "study_indicator")   # as in 03
cat("Union:", length(union), "variables:", paste(union, collapse = ", "), "\n")

#-------------------------------------------------------------------------------
# Final model on the development set (union + all four covariates), validation
#-------------------------------------------------------------------------------
x  <- cbind(X_cv[, union, drop = FALSE],   b$Z_cv_full)
xt <- cbind(X_test[, union, drop = FALSE], b$Z_test_full)
fit <- fit_pbart(x, Y_cv, xt)
colnames(fit$varcount) <- colnames(x)
pip_full <- colMeans(fit$varcount > 0)
final_sel <- fdr_gene_selection(pip_full[union], 0.05)$selected_genes
cov_pip <- pip_full[colnames(b$Z_cv_full)]
cov_sel <- fdr_gene_selection(cov_pip, 0.05)$selected_genes
aucs <- apply(fit$prob.test, 1, function(pr) auc_of(Y_test, pr))
briers <- apply(fit$prob.test, 1, function(pr) mean((Y_test - pr)^2))
validation <- c(AUC = mean(aucs), AUC_lo = unname(quantile(aucs, .025)), AUC_hi = unname(quantile(aucs, .975)),
                Brier = mean(briers), Brier_lo = unname(quantile(briers, .025)), Brier_hi = unname(quantile(briers, .975)))
cat(sprintf("Validation AUC %.3f (%.3f, %.3f); Brier %.3f\n", validation["AUC"], validation["AUC_lo"],
            validation["AUC_hi"], validation["Brier"]))

#-------------------------------------------------------------------------------
# Average marginal effects of each covariate (risk difference over the
# validation set, as for spBART): binary covariates 1 vs 0, age +10 years
#-------------------------------------------------------------------------------
ame <- lapply(colnames(b$Z_cv_full), function(v) {
  x1 <- xt; x0 <- xt
  if (v == "age") { x1[, v] <- xt[, v] + 10 } else { x1[, v] <- 1; x0[, v] <- 0 }
  p1 <- predict(fit, newdata = x1)$prob.test
  p0 <- predict(fit, newdata = x0)$prob.test
  d <- rowMeans(p1 - p0)
  data.frame(covariate = v, contrast = if (v == "age") "+10 years" else "1 vs 0",
             AME = mean(d), lower = unname(quantile(d, .025)), upper = unname(quantile(d, .975)),
             P_positive = mean(d > 0))
})
ame <- do.call(rbind, ame)
print(ame, row.names = FALSE)

saveRDS(list(fold_selected = fold_sel, fold_performance = do.call(rbind, fold_perf),
             fold_covariate_pip = do.call(rbind, fold_cov_pip), union = union,
             final_pip = pip_full, final_selected = final_sel,
             covariate_pip = cov_pip, covariate_selected = cov_sel,
             validation = validation, ame = ame,
             settings = list(nskip = nskip, ndpost = ndpost, keepevery = keepevery, ntree = 200)),
        file.path(out_dir, "BART_application_results.RDS"))
write.csv(ame, file.path(out_dir, "BART_average_marginal_effects.csv"), row.names = FALSE)
cat("Saved results to", out_dir, "\n")
