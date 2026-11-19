################################################################################
# 07_multichain_stability.R
#
# Stability of the gene selection across Markov chains (manuscript Section 4.2
# and Appendix C).
#
# Runs the complete pipeline of 03_real_data_analysis.R (same data, screening,
# development/validation split, folds, fold-wise Bayesian FDR rule, union and
# refit) with four chains for every spBART fit. Each chain has its own seed and
# its own random initial latent variables, and runs 10,000 burn-in iterations
# followed by 10,000 iterations saved every tenth (1,000 draws per chain). For
# every fit, rank-normalized split-Rhat and bulk effective sample size are
# computed across the chains for the log-likelihood, the total number of tree
# splits, each coefficient and the latent means of three fixed participants
# (criterion: Rhat <= 1.02 and bulk ESS >= 390). The fold selections, the union
# and the refit use the pooled draws of the four chains. Each chain's own fold
# selections and union are also reported.
#
# Usage (from this folder or from slurm/):
#   Rscript 07_multichain_stability.R
# Optional: options(spBART.primary_results_dir = "<results_spBART_<date>>") to
# compare with the union of the primary single-chain run of 03 (Jaccard index).
#
# Writes results_spBART_multichain_<date>/: the outputs of 03, plus
#   multichain_diagnostics.csv     Rhat and bulk ESS for every fit
#   multichain_gene_selection.csv  selection of every gene across folds and chains
#   multichain_chain_unions.csv    union of the fold selections of each chain
#   multichain_summary.RDS         all of the above, with the seeds
################################################################################

local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("--file=", "", grep("--file=", args, value = TRUE)[1])
  here <- if (is.na(f)) getwd() else dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE)))
  assign("spbart_driver_dir", here, envir = globalenv())
})

options(spBART.n_chains = 4L,
        spBART.mc_n_burn = 10000, spBART.mc_n_iter = 1000, spBART.mc_n_thin = 10,
        spBART.mc_seed = 20261020,
        spBART.mc_cores = min(4L, parallel::detectCores()),   # chains in parallel; results do not depend on it
        spBART.pd_genes = character(0))   # partial dependence: primary analysis (03) only

source(file.path(spbart_driver_dir, "03_real_data_analysis.R"))

# ------------------------------------------------------------------------------
# Multi-chain summaries
# ------------------------------------------------------------------------------
alpha_mc <- 0.05
fold_labels <- sprintf("CV fold %d", seq_len(K))

diag_tab <- do.call(rbind, lapply(names(spbart_chain_log), function(lb)
  data.frame(fit = lb, spbart_chain_log[[lb]]$diagnostics)))
diag_tab$all_pass_in_fit <- ave(diag_tab$pass, diag_tab$fit, FUN = all)

# PIPs and selections of each chain in each fold
chain_PIP <- lapply(fold_labels, function(lb) lapply(spbart_chain_log[[lb]]$varcount, function(vc) {
  p <- colMeans(vc > 0); names(p) <- colnames(X_cv); p }))
chain_sel <- lapply(chain_PIP, function(fold) lapply(fold, function(p) fdr_gene_selection(p, alpha_mc)$selected_genes))
chain_unions <- lapply(seq_len(spbart_n_chains), function(ch) unique(unlist(lapply(chain_sel, `[[`, ch))))
pooled_union <- unique(unlist(fold_selected_genes))   # pooled draws (03, before adding the study indicator)

genes_any <- unique(c(pooled_union, unlist(chain_unions)))
all_chain_PIPs <- do.call(rbind, unlist(chain_PIP, recursive = FALSE))
gene_tab <- data.frame(
  gene = genes_any,
  chain_fold_fits = vapply(genes_any, function(g) sum(vapply(unlist(chain_sel, recursive = FALSE), function(s) g %in% s, logical(1))), numeric(1)),
  chain_unions = vapply(genes_any, function(g) sum(vapply(chain_unions, function(u) g %in% u, logical(1))), numeric(1)),
  pooled_folds = vapply(genes_any, function(g) sum(vapply(fold_selected_genes, function(s) g %in% s, logical(1))), numeric(1)),
  in_pooled_union = genes_any %in% pooled_union,
  mean_pooled_fold_PIP = colMeans(cv_PIPs_by_fold)[genes_any],
  chain_PIP_min = apply(all_chain_PIPs[, genes_any, drop = FALSE], 2, min),
  chain_PIP_max = apply(all_chain_PIPs[, genes_any, drop = FALSE], 2, max),
  row.names = NULL)
gene_tab <- gene_tab[order(-gene_tab$chain_fold_fits, -gene_tab$mean_pooled_fold_PIP), ]
gene_tab$chain_fold_fits_of <- K * spbart_n_chains

jaccard <- function(a, b) if (length(union(a, b))) length(intersect(a, b)) / length(union(a, b)) else 1
primary_dir <- getOption("spBART.primary_results_dir", NULL)
primary_union <- NULL
if (!is.null(primary_dir) && file.exists(file.path(primary_dir, "cv_gene_selection.RDS"))) {
  primary_union <- unique(unlist(readRDS(file.path(primary_dir, "cv_gene_selection.RDS"))$fold_selected_genes))
}
union_tab <- data.frame(
  analysis = c("Chains pooled", sprintf("Chain %d alone", seq_len(spbart_n_chains))),
  n_selected = c(length(pooled_union), lengths(chain_unions)),
  jaccard_vs_primary = if (is.null(primary_union)) NA else
    vapply(c(list(pooled_union), chain_unions), jaccard, numeric(1), b = primary_union),
  union = vapply(c(list(pooled_union), chain_unions), paste, character(1), collapse = ", "))

write.csv(diag_tab, result_path("multichain_diagnostics.csv"), row.names = FALSE)
write.csv(gene_tab, result_path("multichain_gene_selection.csv"), row.names = FALSE)
write.csv(union_tab, result_path("multichain_chain_unions.csv"), row.names = FALSE)
saveRDS(list(diagnostics = diag_tab, genes = gene_tab, unions = union_tab,
             chain_selections = chain_sel, chain_PIPs = chain_PIP,
             seeds = lapply(spbart_chain_log, `[[`, "seeds"),
             settings = list(n_chains = spbart_n_chains, n_burn = n_burn, n_iter_per_chain = n_iter_per_chain,
                             n_thin = n_thin, alpha = alpha_mc, rhat_max = 1.02, ess_min = 390)),
        result_path("multichain_summary.RDS"))

message("\n=== Multi-chain convergence (criterion: Rhat <= 1.02, bulk ESS >= 390) ===")
for (lb in unique(diag_tab$fit)) {
  d <- diag_tab[diag_tab$fit == lb, ]
  message(sprintf("  %-25s max Rhat %.3f, min bulk ESS %.0f: %s", lb, max(d$rhat), min(d$ess_bulk),
                  if (all(d$pass)) "pass" else paste("FAIL:", paste(d$quantity[!d$pass], collapse = ", "))))
}
message("\n=== Unions of the fold selections ===")
for (i in seq_len(nrow(union_tab))) message(sprintf("  %-15s %s", union_tab$analysis[i], union_tab$union[i]))
message("\nMulti-chain summaries saved to: ", results_dir)
