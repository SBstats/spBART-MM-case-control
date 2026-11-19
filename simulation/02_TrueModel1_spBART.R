################################################################################
# 02_TrueModel1_spBART.R
#
# PURPOSE:
# Fit Semi-parametric Probit BART (spBART) on data generated from True Model 1.
#
# TRUE MODEL 1: Nonparametric with Gene-Covariate Interactions
# f(X, Z) includes gene-covariate interaction terms (x7*z1, x8*z2, x9*z4),
# making the spBART additive assumption MISSPECIFIED.
# This tests robustness of spBART to model misspecification.
#
# VARIABLE SELECTION PIPELINE:
# 1. No initial screening - all p genes passed directly to BART
# 2. 5-fold CV: fit spBART on all genes, apply FDR control (alpha = 0.05)
# 3. Final selection: union of genes selected across all folds
# 4. Final model: fit spBART on union genes to compute AUC, Brier, RMSE
#
# SIMULATION DESIGN:
# - Sample sizes: n = 1500, 2500
# - Number of genes: p = 500, 2000, 3000
# - True signal genes: 10 (fixed across all p values)
#
# SLURM ARRAY JOB:
# - Each array task (1-500) runs one replicate for ALL n,p combinations
#
# PREREQUISITE: Rscript BART_package/00_install_BART.R
#
# OUTPUTS:
# - results/TrueModel1/n{1500,2500}/p{500,2000,3000}/spBART/rep_{001-500}.RDS
#
################################################################################

# ==============================================================================
# SECTION 1: Get Replicate ID from SLURM Array
# ==============================================================================

rep_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if (is.na(rep_id)) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    rep_id <- as.integer(args[1])
  } else {
    rep_id <- 1
    cat("WARNING: No SLURM_ARRAY_TASK_ID found. Using rep_id = 1 for testing.\n")
  }
}

cat("\n")
cat("================================================================================\n")
cat("  02_TrueModel1_spBART.R: Semi-parametric Probit BART\n")
cat("  True Model 1: Gene-covariate interactions (spBART MISSPECIFIED)\n")
cat("  Variable Selection: 5-Fold CV FDR Union (No Initial Screening)\n")
cat("  Replicate ID:", rep_id, "\n")
cat("================================================================================\n")
cat("Script started at:", format(Sys.time()), "\n\n")

# ==============================================================================
# SECTION 2: Setup
# ==============================================================================

# Determine script directory robustly for both interactive and batch modes
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
spbart_root <- dirname(script_dir)   # spBART/
cat("Script directory:", script_dir, "\n")

source(file.path(script_dir, "utils", "data_generation.R"))

cat("\nLoading required packages...\n")

required_packages <- c("MASS", "truncnorm", "caret", "pROC")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# Shared spBART sampler (fit_spBART, fdr_gene_selection) and the modified BART
# package check (version 2.9.9.1, installed from BART_package/)
source(file.path(spbart_root, "R", "spBART_sampler.R"))
load_spBART_BART()

# ==============================================================================
# SECTION 3: spBART Gibbs Sampler
# ==============================================================================
# fit_spBART() is defined in R/spBART_sampler.R (shared by all analyses).

# ==============================================================================
# SECTION 4: Simulation Settings
# ==============================================================================

sample_sizes <- c(1500, 2500)
gene_counts <- c(500, 2000, 3000)

ntrees <- 200
n_burn <- 2000    # Burn-in iterations
n_iter <- 1000    # Number of posterior samples to save
n_thin <- 5       # Thinning: save every 5th iteration (5000 post-burn / 5 = 1000 saved)
sparse <- TRUE
a <- 0.5
b <- 1
sigma_beta <- sqrt(10)   # beta ~ N(0, 10 I), as in the manuscript
n_test <- 1000           # independent test set for AUC / Brier

alpha_FDR <- 0.05
K_folds <- 5
true_model <- 1

results_root <- "results"

cat("================================================================================\n")
cat("  Simulation Settings\n")
cat("================================================================================\n")
cat("Replicate ID:", rep_id, "\n")
cat("True Model:", true_model, "(spBART MISSPECIFIED)\n")
cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
cat("Gene counts:", paste(gene_counts, collapse = ", "), "\n")
cat("spBART: ntrees =", ntrees, ", n_burn =", n_burn, ", n_iter =", n_iter, ", n_thin =", n_thin, "\n")
cat("FDR threshold:", alpha_FDR, "\n")
cat("No initial screening - all genes passed to BART\n\n")

# ==============================================================================
# SECTION 5: Run Replicate for Each (n, p) Combination
# ==============================================================================

for (n in sample_sizes) {
  for (p_genes in gene_counts) {

    cat("\n")
    cat("################################################################################\n")
    cat("  Running replicate", rep_id, "for n =", n, ", p =", p_genes, "\n")
    cat("################################################################################\n\n")

    results_dir <- file.path(script_dir, results_root, "TrueModel1",
                             paste0("n", n), paste0("p", p_genes), "spBART")
    if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

    output_file <- file.path(results_dir, sprintf("rep_%03d.RDS", rep_id))
    if (file.exists(output_file)) {
      cat("Replicate", rep_id, "already completed for n =", n, ", p =", p_genes, ". Skipping.\n")
      next
    }

    # Seed depends on n, p, and rep_id for reproducibility
    set.seed(1000 * n + p_genes + rep_id)
    rep_start <- Sys.time()

    tryCatch({
      # ========================================================================
      # Step 1: Generate data
      # ========================================================================
      cat("Step 1: Generating simulated data...\n")
      sim_data <- generate_simulation_data(
        n = n, p_genes = p_genes, n_true_genes = 10,
        true_model = true_model, seed = 1000 * n + p_genes + rep_id
      )

      # Independent test set from the same DGM (same fixed trees, fresh draws)
      test_data <- generate_simulation_data(
        n = n_test, p_genes = p_genes, n_true_genes = 10,
        true_model = true_model, seed = 1e8 + 1000 * n + p_genes + rep_id
      )
      X_genes_test <- test_data$X_genes
      Z_cov_test <- test_data$Z_covariates
      Y_test <- test_data$Y

      X_genes <- sim_data$X_genes
      Z_cov <- sim_data$Z_covariates
      Y <- sim_data$Y
      f_true <- sim_data$f_true  # True latent function values

      gene_names <- colnames(X_genes)
      true_genes <- sim_data$true_gene_names
      beta_true <- sim_data$beta_true

      cat("  n =", n, ", p =", ncol(X_genes), ", Y=1:", sum(Y), "\n")
      cat("  True genes:", paste(true_genes, collapse = ", "), "\n")

      # ========================================================================
      # Step 2: 5-Fold CV with FDR control on ALL genes (no screening)
      # ========================================================================
      cat("Step 2: 5-Fold CV with FDR control on", p_genes, "genes (no screening)...\n")

      # Create stratified folds to ensure balanced Y=0/Y=1 in each fold
      fold_assignment <- createFolds(factor(Y), k = K_folds, list = TRUE, returnTrain = FALSE)

      # Verify fold balance
      cat("  Fold balance check:\n")
      for (k in 1:K_folds) {
        fold_k_indices <- fold_assignment[[k]]
        n_cases_k <- sum(Y[fold_k_indices] == 1)
        n_controls_k <- sum(Y[fold_k_indices] == 0)
        cat(sprintf("    Fold %d: %d patients (%d cases, %d controls)\n",
                    k, length(fold_k_indices), n_cases_k, n_controls_k))
      }

      fold_selected_genes <- vector("list", K_folds)
      fold_PIPs <- vector("list", K_folds)

      for (k in 1:K_folds) {
        cat("  Fold", k, "/", K_folds, "...\n")

        val_indices <- fold_assignment[[k]]
        train_indices <- setdiff(1:n, val_indices)

        X_train_k <- X_genes[train_indices, , drop = FALSE]
        Z_train_k <- Z_cov[train_indices, , drop = FALSE]
        Y_train_k <- Y[train_indices]

        fit_k <- fit_spBART(
          X = X_train_k, Z = Z_train_k, Y = Y_train_k,
          n_burn = n_burn, n_iter = n_iter, n_thin = n_thin,
          n_trees = ntrees, sparse = sparse,
          a = a, b = b, rho = p_genes,
          sigma_beta = sigma_beta
        )

        # PIPs
        var_used_k <- fit_k$varcount > 0
        PIP_k <- colMeans(var_used_k)
        names(PIP_k) <- gene_names
        fold_PIPs[[k]] <- PIP_k

        # FDR selection
        fdr_result_k <- fdr_gene_selection(PIP_k, alpha = alpha_FDR)
        fold_selected_genes[[k]] <- fdr_result_k$selected_genes
        cat("    Selected:", length(fdr_result_k$selected_genes), "genes\n")

        # Explicit cleanup to free memory
        rm(fit_k, var_used_k)
        gc()
      }

      # ========================================================================
      # Step 3: Union of selected genes
      # ========================================================================
      cat("Step 3: Computing union of selected genes...\n")
      union_genes <- unique(unlist(fold_selected_genes))
      cat("  Union size:", length(union_genes), "\n")

      gene_fold_frequency <- if (length(union_genes) > 0) {
        sapply(union_genes, function(g) sum(sapply(fold_selected_genes, function(fg) g %in% fg)))
      } else integer(0)

      # ========================================================================
      # Step 4: Compute selection metrics
      # ========================================================================
      cat("Step 4: Computing selection metrics...\n")

      selected_mask <- gene_names %in% union_genes
      true_mask <- gene_names %in% true_genes

      TP <- sum(selected_mask & true_mask)
      FP <- sum(selected_mask & !true_mask)
      TN <- sum(!selected_mask & !true_mask)
      FN <- sum(!selected_mask & true_mask)

      sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
      specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 1)
      FDR_actual <- ifelse((TP + FP) > 0, FP / (TP + FP), 0)
      PPV <- 1 - FDR_actual

      cat("  TP:", TP, "FP:", FP, "FN:", FN, "\n")
      cat("  Sensitivity:", round(sensitivity, 3), "\n")
      cat("  FDR:", round(FDR_actual, 3), "\n")

      # ========================================================================
      # Step 5: Fit final model on union genes and compute AUC, Brier, RMSE
      # ========================================================================
      cat("Step 5: Fitting final model on union genes...\n")

      if (length(union_genes) > 0) {
        # Create design matrix with union genes
        X_union <- X_genes[, union_genes, drop = FALSE]

        # Fit final spBART model on full data
        final_fit <- fit_spBART(
          X = X_union, Z = Z_cov, Y = Y,
          X_test = X_genes_test[, union_genes, drop = FALSE], Z_test = Z_cov_test,
          n_burn = n_burn, n_iter = n_iter, n_thin = n_thin,
          n_trees = ntrees, sparse = sparse,
          a = a, b = b, rho = length(union_genes),
          sigma_beta = sigma_beta
        )

        # Predicted probabilities
        mu_hat <- final_fit$mu_hat
        prob_hat <- final_fit$prob.test.mean   # posterior mean of Phi(.) on the test set

        # AUC
        if (length(unique(Y_test)) == 2) {
          roc_obj <- tryCatch({
            roc(Y_test, prob_hat, quiet = TRUE, levels = c(0, 1), direction = "<")
          }, error = function(e) NULL)
          final_AUC <- if (!is.null(roc_obj)) as.numeric(auc(roc_obj)) else NA
        } else {
          final_AUC <- NA
        }

        # Brier Score
        final_Brier <- mean((Y_test - prob_hat)^2)

        # RMSE of latent function
        final_RMSE <- sqrt(mean((mu_hat - f_true)^2))

        # Beta estimation
        final_beta_draws <- final_fit$beta_draws
        final_beta_hat <- colMeans(final_beta_draws)
        final_beta_sd <- apply(final_beta_draws, 2, sd)
        final_beta_ci_lower <- apply(final_beta_draws, 2, quantile, 0.025)
        final_beta_ci_upper <- apply(final_beta_draws, 2, quantile, 0.975)
        final_beta_coverage <- (beta_true >= final_beta_ci_lower) & (beta_true <= final_beta_ci_upper)

        # Beta estimation metrics (per covariate)
        final_beta_bias <- final_beta_hat - beta_true
        final_beta_mse <- final_beta_bias^2
        final_beta_rmse <- abs(final_beta_bias)
        final_beta_ci_width <- final_beta_ci_upper - final_beta_ci_lower

        # Compute posterior probability of correct sign
        final_beta_sign_prob <- sapply(1:ncol(final_beta_draws), function(j) {
          if (is.na(beta_true[j])) return(NA)  # Model 1 has no true beta
          mean(sign(final_beta_draws[, j]) == sign(beta_true[j]))
        })
        names(final_beta_sign_prob) <- colnames(final_beta_draws)

      } else {
        # No genes selected - use covariates only
        cat("  WARNING: No genes selected. Fitting model with covariates only.\n")

        # For spBART with no genes, just fit a probit model with covariates
        # Use a single-column dummy matrix for X
        # Non-constant dummy column (a constant column would be dropped by BART)
        X_dummy <- matrix(rnorm(n), ncol = 1, dimnames = list(NULL, "dummy"))
        X_dummy_test <- matrix(rnorm(n_test), ncol = 1, dimnames = list(NULL, "dummy"))

        final_fit <- fit_spBART(
          X = X_dummy, Z = Z_cov, Y = Y,
          X_test = X_dummy_test, Z_test = Z_cov_test,
          n_burn = n_burn, n_iter = n_iter, n_thin = n_thin,
          n_trees = ntrees, sparse = sparse,
          a = a, b = b, rho = 1,
          sigma_beta = sigma_beta
        )

        mu_hat <- final_fit$mu_hat
        prob_hat <- final_fit$prob.test.mean   # posterior mean of Phi(.) on the test set

        if (length(unique(Y_test)) == 2) {
          roc_obj <- tryCatch({
            roc(Y_test, prob_hat, quiet = TRUE, levels = c(0, 1), direction = "<")
          }, error = function(e) NULL)
          final_AUC <- if (!is.null(roc_obj)) as.numeric(auc(roc_obj)) else NA
        } else {
          final_AUC <- NA
        }

        final_Brier <- mean((Y_test - prob_hat)^2)
        final_RMSE <- sqrt(mean((mu_hat - f_true)^2))

        final_beta_draws <- final_fit$beta_draws
        final_beta_hat <- colMeans(final_beta_draws)
        final_beta_sd <- apply(final_beta_draws, 2, sd)
        final_beta_ci_lower <- apply(final_beta_draws, 2, quantile, 0.025)
        final_beta_ci_upper <- apply(final_beta_draws, 2, quantile, 0.975)
        final_beta_coverage <- (beta_true >= final_beta_ci_lower) & (beta_true <= final_beta_ci_upper)

        # Beta estimation metrics (per covariate)
        final_beta_bias <- final_beta_hat - beta_true
        final_beta_mse <- final_beta_bias^2
        final_beta_rmse <- abs(final_beta_bias)
        final_beta_ci_width <- final_beta_ci_upper - final_beta_ci_lower

        # Compute posterior probability of correct sign
        final_beta_sign_prob <- sapply(1:ncol(final_beta_draws), function(j) {
          if (is.na(beta_true[j])) return(NA)  # Model 1 has no true beta
          mean(sign(final_beta_draws[, j]) == sign(beta_true[j]))
        })
        names(final_beta_sign_prob) <- colnames(final_beta_draws)
      }

      cat("  Final AUC:", round(final_AUC, 3), "\n")
      cat("  Final Brier:", round(final_Brier, 4), "\n")
      cat("  Final RMSE:", round(final_RMSE, 4), "\n")
      cat("  Beta estimates:", paste(round(final_beta_hat, 3), collapse = ", "), "\n")

      # ========================================================================
      # Step 6: Save results
      # ========================================================================
      rep_end <- Sys.time()
      runtime <- difftime(rep_end, rep_start, units = "mins")

      result <- list(
        # Variable selection
        union_genes = union_genes,
        n_union_genes = length(union_genes),
        gene_fold_frequency = gene_fold_frequency,
        fold_selected_genes = fold_selected_genes,
        fold_PIPs = fold_PIPs,

        # Selection metrics
        TP = TP, FP = FP, TN = TN, FN = FN,
        sensitivity = sensitivity, specificity = specificity,
        FDR_actual = FDR_actual, PPV = PPV,

        # Final model metrics (single values)
        final_AUC = final_AUC,
        final_Brier = final_Brier,
        final_RMSE = final_RMSE,           # latent index vs truth, training data
        n_test = n_test,                   # AUC / Brier are on this test set

        # Beta estimation from final model (all per-covariate vectors)
        beta_hat = final_beta_hat,
        beta_sd = final_beta_sd,
        beta_bias = final_beta_bias,
        beta_mse = final_beta_mse,
        beta_rmse = final_beta_rmse,
        beta_ci_lower = final_beta_ci_lower,
        beta_ci_upper = final_beta_ci_upper,
        beta_ci_width = final_beta_ci_width,
        beta_coverage = final_beta_coverage,
        beta_sign_prob = final_beta_sign_prob,

        # Ground truth
        true_gene_names = true_genes,
        n_true_genes = length(true_genes),
        beta_true = beta_true,

        # Metadata
        sample_size = n,
        p_genes = p_genes,
        replicate_id = rep_id,
        true_model = true_model,
        method = "spBART",
        runtime_mins = as.numeric(runtime),
        seed = 1000 * n + p_genes + rep_id
      )

      saveRDS(result, file = output_file)
      cat("\nCompleted n =", n, ", p =", p_genes, "in", round(as.numeric(runtime), 2), "minutes\n")

    }, error = function(e) {
      cat("ERROR in replicate", rep_id, "for n =", n, ", p =", p_genes, ":", conditionMessage(e), "\n")
    })

    gc()
  }
}

cat("\n")
cat("================================================================================\n")
cat("  Replicate", rep_id, "Complete!\n")
cat("================================================================================\n")
cat("Script completed at:", format(Sys.time()), "\n")
