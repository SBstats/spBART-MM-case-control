################################################################################
# 09_supporting_summaries.R
#
# Summaries reported alongside the main results, computed from the latest
# results folders of 03_real_data_analysis.R (pooled analysis) and
# 08_BConly_analysis.R (BC Case-control alone). No model is refitted.
#   (1) Validation AUC of the posterior-mean predicted probabilities with
#       DeLong and stratified bootstrap 95% confidence intervals, in the pooled
#       validation set, in its BC Case-control participants and in the
#       validation set of the BC-only analysis
#   (2) For each gene saved by 03 (default IL1RAP and ST5): probit coefficient
#       of MM status on the variance-stabilized level, adjusted for age and
#       sex, in the pooled development set and in all BC Case-control
#       participants (per unit of the variance-stabilized scale, and per
#       standard deviation of the gene within the same participants), with
#       the medians among BC cases and controls and a Wilcoxon test
#   (3) Violin plot of these levels among BC cases and controls
# Writes real_data_analysis/summary/auc_confidence_intervals.csv,
# real_data_analysis/summary/gene_coefficients.csv and
# real_data_analysis/summary/BC_<genes>_violin.png.
#
# Usage (after real_03 and real_08 have run):  Rscript real_data_analysis/09_supporting_summaries.R
################################################################################
suppressPackageStartupMessages({ library(pROC); library(ggplot2) })
rd <- tryCatch({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE)))
}, error = function(e) getwd())
out <- file.path(rd, "summary"); dir.create(out, showWarnings = FALSE)
latest <- function(pattern) { d <- sort(list.files(rd, pattern = pattern, full.names = TRUE)); if (length(d)) tail(d, 1) else NA }
sp_dir <- latest("^results_spBART_[0-9]")
bc_dir <- latest("^results_spBART_BConly_[0-9]")
if (is.na(sp_dir) || !file.exists(file.path(sp_dir, "summary_inputs.RDS")))
  stop("No pooled spBART results with summary_inputs.RDS; run 03_real_data_analysis.R first")
cat("Pooled results:", basename(sp_dir), "\nBC-only results:", if (is.na(bc_dir)) "none" else basename(bc_dir), "\n")

si <- readRDS(file.path(sp_dir, "summary_inputs.RDS"))
meta <- si$metadata
stopifnot(identical(rownames(si$vst), as.character(meta$PooledData_ID)))   # genes and participants aligned
tp <- readRDS(file.path(sp_dir, "test_set_performance.RDS"))
stopifnot(identical(as.numeric(meta$outcome[si$test_indices]), as.numeric(tp$Y_test)))

#-------------------------------------------------------------------------------
# (1) AUC of the posterior-mean predictions, with DeLong and bootstrap intervals
#-------------------------------------------------------------------------------
auc_row <- function(y, p_mean, p_draws, label, n_boot = 2000, seed = 20261031) {
  r <- roc(y, p_mean, quiet = TRUE, levels = c(0, 1), direction = "<")
  dl <- ci.auc(r, method = "delong")
  set.seed(seed)
  bs <- ci.auc(r, method = "bootstrap", boot.n = n_boot, boot.stratified = TRUE)
  post <- apply(p_draws, 1, function(p) as.numeric(auc(roc(y, p, quiet = TRUE, levels = c(0, 1), direction = "<"))))
  data.frame(evaluation_set = label, n = length(y), cases = sum(y == 1),
             AUC_posterior_mean_predictions = as.numeric(auc(r)),
             DeLong_lower = dl[1], DeLong_upper = dl[3], bootstrap_lower = bs[1], bootstrap_upper = bs[3],
             AUC_posterior_mean = mean(post), AUC_CrI_lower = unname(quantile(post, .025)),
             AUC_CrI_upper = unname(quantile(post, .975)), row.names = NULL)
}
y_val <- as.numeric(tp$Y_test)
bc_val <- which(meta$study_indicator[si$test_indices] == 0)          # S = 0: BC Case-control
auc_tab <- rbind(
  auc_row(y_val, tp$predictions_mean, tp$predictions_full, "Pooled validation set"),
  auc_row(y_val[bc_val], tp$predictions_mean[bc_val], tp$predictions_full[, bc_val, drop = FALSE],
          "BC Case-control participants of the pooled validation set"))
if (!is.na(bc_dir)) {
  bo <- readRDS(file.path(bc_dir, "test_set_performance.RDS"))
  auc_tab <- rbind(auc_tab, auc_row(as.numeric(bo$Y_test), bo$predictions_mean, bo$predictions_full,
                                    "Validation set of the BC-only analysis"))
} else cat("No BC-only results yet (08_BConly_analysis.R); skipping its validation set\n")
write.csv(auc_tab, file.path(out, "auc_confidence_intervals.csv"), row.names = FALSE)
print(auc_tab, digits = 3, row.names = FALSE)

#-------------------------------------------------------------------------------
# (2) Probit coefficients of the genes, adjusted for age and sex
#-------------------------------------------------------------------------------
genes <- colnames(si$vst)
if (length(genes) == 0) stop("summary_inputs.RDS contains no genes")
d_all <- data.frame(meta[, c("outcome", "study_indicator", "age", "sex_binary")], si$vst, check.names = FALSE)
sets <- list("Pooled development set" = d_all[si$cv_indices, ],
             "All BC Case-control participants" = d_all[d_all$study_indicator == 0, ])
probit_coef <- function(d, g) {
  m <- glm(reformulate(c(sprintf("`%s`", g), "age", "sex_binary"), "outcome"),
           family = binomial(link = "probit"), data = d)
  cs <- summary(m)$coefficients[2, ]
  data.frame(coef_per_unit = unname(cs[1]), se = unname(cs[2]), p_value = unname(cs[4]),
             gene_sd = sd(d[[g]]), coef_per_sd = unname(cs[1]) * sd(d[[g]]))
}
bc <- sets[["All BC Case-control participants"]]
gene_tab <- do.call(rbind, lapply(genes, function(g) do.call(rbind, lapply(names(sets), function(s) {
  d <- sets[[s]]
  data.frame(gene = g, participants = s, n = nrow(d), cases = sum(d$outcome == 1), probit_coef(d, g),
             BC_median_cases = median(bc[[g]][bc$outcome == 1]), BC_median_controls = median(bc[[g]][bc$outcome == 0]),
             BC_wilcoxon_p = wilcox.test(bc[[g]] ~ bc$outcome)$p.value)
}))))
write.csv(gene_tab, file.path(out, "gene_coefficients.csv"), row.names = FALSE)
print(gene_tab[, c("gene", "participants", "n", "cases", "coef_per_unit", "se", "coef_per_sd")], digits = 3, row.names = FALSE)

#-------------------------------------------------------------------------------
# (3) Levels among BC cases and controls
#-------------------------------------------------------------------------------
long <- do.call(rbind, lapply(genes, function(g) data.frame(gene = g, vst = bc[[g]], status = bc$outcome)))
long$gene <- factor(long$gene, levels = genes)
long$status <- factor(long$status, levels = c(0, 1), labels = c("Controls", "MM cases"))
# points below the 0.1th percentile of a gene are not shown (their number is printed)
lo <- tapply(long$vst, long$gene, quantile, probs = 0.001)
shown <- long$vst >= lo[as.character(long$gene)]
cat("Points below the plotted range:", paste(names(lo), tapply(!shown, long$gene, sum), collapse = ", "), "\n")
pl <- ggplot(long[shown, ], aes(status, vst, fill = status)) +
  geom_violin(trim = TRUE, alpha = 0.6, colour = NA) + geom_boxplot(width = 0.15, outlier.size = 0.4, fill = "white") +
  facet_wrap(~ gene, scales = "free_y") + scale_fill_manual(values = c("#4C72B0", "#DD8452"), guide = "none") +
  labs(x = NULL, y = "5hmC level (variance-stabilized)") + theme_bw(base_size = 11)
png(file.path(out, paste0("BC_", paste(genes, collapse = "_"), "_violin.png")), width = 6, height = 3.4, units = "in", res = 300)
print(pl); invisible(dev.off())
cat("Summaries saved to:", out, "\n")
