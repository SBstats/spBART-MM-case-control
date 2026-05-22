################################################################################
# 01_TrueModel1_regular_BART.R
#
# PURPOSE:
# Fit Regular BART with sparsity (DART) on data generated from True Model 1.
#
# TRUE MODEL 1: Nonparametric with Gene-Covariate Interactions
# f(X, Z) includes gene-covariate interaction terms (x7*z1, x8*z2, x9*z4),
# making the semi-parametric assumption of spBART MISSPECIFIED.
# Regular BART should capture these interactions and perform well.
#
# VARIABLE SELECTION PIPELINE:
# 1. No initial screening - all p genes passed directly to BART
# 2. 5-fold CV: fit BART on all genes, apply FDR control (alpha = 0.05)
# 3. Final selection: union of genes selected across all folds
# 4. Final model: fit BART on union genes to compute AUC, Brier, RMSE
#
# SIMULATION DESIGN:
# - Sample sizes: n = 1500, 2500
# - Number of genes: p = 500, 2000, 3000
# - True signal genes: 10 (fixed across all p values)
#
# SLURM ARRAY JOB:
# - Each array task (1-500) runs one replicate for ALL n,p combinations
# - Usage: sbatch 01_TrueModel1_regular_BART.sbatch
#
# OUTPUTS:
# - results/TrueModel1/n{1500,2500}/p{500,2000,3000}/regular_BART/rep_{001-500}.RDS
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
cat("  01_TrueModel1_regular_BART.R: Regular BART with Sparsity (DART)\n")
cat("  True Model 1: Gene-covariate interactions (spBART misspecified)\n")
cat("  Variable Selection: 5-Fold CV FDR Union (No Initial Screening)\n")
cat("  Replicate ID:", rep_id, "\n")
cat("================================================================================\n")
cat("Script started at:", format(Sys.time()), "\n\n")

# ==============================================================================
# SECTION 2: Setup
# ==============================================================================

# Determine script directory robustly for both interactive and batch modes
script_dir <- NULL

# Method 1: SLURM environment variable (best for HPC batch jobs)
if (Sys.getenv("SLURM_SUBMIT_DIR") != "") {
  script_dir <- Sys.getenv("SLURM_SUBMIT_DIR")
}

# Method 2: Works when script is run via source()
if (is.null(script_dir) || script_dir == "" || script_dir == ".") {
  script_dir <- tryCatch({
    dirname(sys.frame(1)$ofile)
  }, error = function(e) NULL)
}

# Method 3: Works in RStudio when running interactively
if (is.null(script_dir) || script_dir == "" || script_dir == ".") {
  script_dir <- tryCatch({
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      dirname(rstudioapi::getActiveDocumentContext()$path)
    } else {
      NULL
    }
  }, error = function(e) NULL)
}

# Method 4: Fallback to getwd() if all else fails
if (is.null(script_dir) || script_dir == "" || script_dir == ".") {
  script_dir <- getwd()
  warning("Could not determine script directory. Using getwd(): ", script_dir)
}

cat("Script directory:", script_dir, "\n")

# Source data generation utilities
source(file.path(script_dir, "utils", "data_generation.R"))

# Load required packages
cat("\nLoading required packages...\n")

required_packages <- c("caret", "pROC")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

library(BART)
cat("Packages loaded successfully!\n\n")

# ==============================================================================
# SECTION 3: Simulation Settings
# ==============================================================================

sample_sizes <- c(1500, 2500)
gene_counts <- c(500, 2000, 3000)

# DART parameters
ntrees <- 200
ndpost <- 1000    # Number of posterior samples to save
nskip <- 2000     # Burn-in iterations
keepevery <- 5    # Thinning: save every 5th iteration (5000 post-burn / 5 = 1000 saved)
sparse <- TRUE
a <- 0.5
b <- 1

# FDR control threshold
alpha_FDR <- 0.05

# Number of CV folds
K_folds <- 5

# True model
true_model <- 1

cat("================================================================================\n")
cat("  Simulation Settings\n")
cat("================================================================================\n")
cat("Replicate ID:", rep_id, "\n")
cat("True Model:", true_model, "(gene-covariate interactions)\n")
cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
cat("Gene counts:", paste(gene_counts, collapse = ", "), "\n")
cat("DART parameters: ntrees =", ntrees, ", ndpost =", ndpost,
    ", nskip =", nskip, ", keepevery =", keepevery, ", sparse =", sparse, "\n")
cat("Sparsity prior: a =", a, ", b =", b, "\n")
cat("FDR threshold:", alpha_FDR, "\n")
cat("CV folds:", K_folds, "\n")
cat("No initial screening - all genes passed to BART\n\n")

# ==============================================================================
# SECTION 4: Run Replicate for Each (n, p) Combination
# ==============================================================================

for (n in sample_sizes) {
  for (p_genes in gene_counts) {

    cat("\n")
    cat("################################################################################\n")
    cat("  Running replicate", rep_id, "for n =", n, ", p =", p_genes, "\n")
    cat("################################################################################\n\n")

    # Results directory (now includes p in path)
    results_dir <- file.path(script_dir, "results", "TrueModel1",
                             paste0("n", n), paste0("p", p_genes), "regular_BART")

    if (!dir.exists(results_dir)) {
      dir.create(results_dir, recursive = TRUE)
    }

    # Check if already done
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
        n = n,
        p_genes = p_genes,
        n_true_genes = 10,
        true_model = true_model,
        seed = 1000 * n + p_genes + rep_id
      )

      X_genes <- sim_data$X_genes
      Z_cov <- sim_data$Z_covariates
      Y <- sim_data$Y
      f_true <- sim_data$f_true  # True latent function values

      gene_names <- colnames(X_genes)
      true_genes <- sim_data$true_gene_names

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
      fold_selected_covariates <- vector("list", K_folds)
      fold_covariate_PIPs <- vector("list", K_folds)

      for (k in 1:K_folds) {
        cat("  Fold", k, "/", K_folds, "...\n")

        val_indices <- fold_assignment[[k]]
        train_indices <- setdiff(1:n, val_indices)

        # Combined design matrix (all genes + covariates)
        X_train_k <- cbind(X_genes[train_indices, ], Z_cov[train_indices, ])
        Y_train_k <- Y[train_indices]

        # Fit probit BART with DART (no test set needed for variable selection)
        fit_k <- pbart(
          x.train = X_train_k,
          y.train = Y_train_k,
          ntree = ntrees,
          ndpost = ndpost,
          nskip = nskip,
          keepevery = keepevery,
          sparse = sparse,
          a = a, b = b,
          rho = ncol(X_train_k)
        )

        # Extract PIPs for genes only (not covariates)
        varcount_k <- fit_k$varcount[, 1:p_genes, drop = FALSE]
        var_used_k <- varcount_k > 0
        PIP_k <- colMeans(var_used_k)
        names(PIP_k) <- gene_names

        fold_PIPs[[k]] <- PIP_k

        # FDR control for genes
        fdr_result_k <- fdr_gene_selection(PIP_k, alpha = alpha_FDR)
        fold_selected_genes[[k]] <- fdr_result_k$selected_genes
        cat("    Selected:", length(fdr_result_k$selected_genes), "genes\n")

        # Extract PIPs for covariates (in parallel with gene selection)
        n_total_vars <- ncol(X_train_k)  # p_genes + 4 covariates
        covariate_indices <- (p_genes + 1):n_total_vars
        varcount_cov_k <- fit_k$varcount[, covariate_indices, drop = FALSE]
        var_used_cov_k <- varcount_cov_k > 0
        PIP_cov_k <- colMeans(var_used_cov_k)
        names(PIP_cov_k) <- colnames(Z_cov)

        fold_covariate_PIPs[[k]] <- PIP_cov_k

        # FDR control for covariates
        fdr_result_cov_k <- fdr_gene_selection(PIP_cov_k, alpha = alpha_FDR)
        fold_selected_covariates[[k]] <- fdr_result_cov_k$selected_genes
        cat("    Selected:", length(fdr_result_cov_k$selected_genes), "covariates\n")
      }

      # ========================================================================
      # Step 3: Union of selected genes and covariates
      # ========================================================================
      cat("Step 3: Computing union of selected genes and covariates...\n")
      union_genes <- unique(unlist(fold_selected_genes))
      cat("  Union genes:", length(union_genes), "\n")

      # Gene frequency across folds
      if (length(union_genes) > 0) {
        gene_fold_frequency <- sapply(union_genes, function(gene) {
          sum(sapply(fold_selected_genes, function(fg) gene %in% fg))
        })
      } else {
        gene_fold_frequency <- integer(0)
      }

      # Union of selected covariates
      union_covariates <- unique(unlist(fold_selected_covariates))
      cat("  Union covariates:", length(union_covariates), "\n")

      # Covariate frequency across folds
      if (length(union_covariates) > 0) {
        covariate_fold_frequency <- sapply(union_covariates, function(cov) {
          sum(sapply(fold_selected_covariates, function(fc) cov %in% fc))
        })
      } else {
        covariate_fold_frequency <- integer(0)
      }

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
        # Create design matrix with union genes + covariates
        X_union <- cbind(X_genes[, union_genes, drop = FALSE], Z_cov)

        # Fit final model on full data
        final_fit <- pbart(
          x.train = X_union,
          y.train = Y,
          ntree = ntrees,
          ndpost = ndpost,
          nskip = nskip,
          keepevery = keepevery,
          sparse = sparse,
          a = a, b = b,
          rho = ncol(X_union)
        )

        # Predicted probabilities (posterior mean)
        f_hat <- colMeans(final_fit$yhat.train)
        prob_hat <- pnorm(f_hat)

        # AUC
        if (length(unique(Y)) == 2) {
          roc_obj <- tryCatch({
            roc(Y, prob_hat, quiet = TRUE, levels = c(0, 1), direction = "<")
          }, error = function(e) NULL)
          final_AUC <- if (!is.null(roc_obj)) as.numeric(auc(roc_obj)) else NA
        } else {
          final_AUC <- NA
        }

        # Brier Score
        final_Brier <- mean((Y - prob_hat)^2)

        # RMSE of latent function (f_hat vs f_true)
        final_RMSE <- sqrt(mean((f_hat - f_true)^2))

        cat("  Final AUC:", round(final_AUC, 3), "\n")
        cat("  Final Brier:", round(final_Brier, 4), "\n")
        cat("  Final RMSE:", round(final_RMSE, 4), "\n")

      } else {
        # No genes selected - use null model (intercept only via covariates)
        cat("  WARNING: No genes selected. Fitting model with covariates only.\n")

        final_fit <- pbart(
          x.train = Z_cov,
          y.train = Y,
          ntree = ntrees,
          ndpost = ndpost,
          nskip = nskip,
          keepevery = keepevery,
          sparse = sparse,
          a = a, b = b,
          rho = ncol(Z_cov)
        )

        f_hat <- colMeans(final_fit$yhat.train)
        prob_hat <- pnorm(f_hat)

        if (length(unique(Y)) == 2) {
          roc_obj <- tryCatch({
            roc(Y, prob_hat, quiet = TRUE, levels = c(0, 1), direction = "<")
          }, error = function(e) NULL)
          final_AUC <- if (!is.null(roc_obj)) as.numeric(auc(roc_obj)) else NA
        } else {
          final_AUC <- NA
        }

        final_Brier <- mean((Y - prob_hat)^2)
        final_RMSE <- sqrt(mean((f_hat - f_true)^2))

        cat("  Final AUC:", round(final_AUC, 3), "\n")
        cat("  Final Brier:", round(final_Brier, 4), "\n")
        cat("  Final RMSE:", round(final_RMSE, 4), "\n")
      }

      # ========================================================================
      # Step 6: Save results
      # ========================================================================
      rep_end <- Sys.time()
      runtime <- difftime(rep_end, rep_start, units = "mins")

      result <- list(
        # Gene selection
        union_genes = union_genes,
        n_union_genes = length(union_genes),
        gene_fold_frequency = gene_fold_frequency,
        fold_selected_genes = fold_selected_genes,
        fold_PIPs = fold_PIPs,

        # NOTE: Covariate selection is computed in-loop for diagnostic purposes
        # but not persisted here. Writeup Table 7.3 reports covariate selection
        # for DGM 2 only (script 03_TrueModel2_regular_BART.R). DGM 1 has no
        # corresponding writeup table.

        # Gene selection metrics
        TP = TP, FP = FP, TN = TN, FN = FN,
        sensitivity = sensitivity,
        specificity = specificity,
        FDR_actual = FDR_actual,
        PPV = PPV,

        # Final model metrics (single values)
        final_AUC = final_AUC,
        final_Brier = final_Brier,
        final_RMSE = final_RMSE,

        # Ground truth
        true_gene_names = true_genes,
        n_true_genes = length(true_genes),
        beta_true = sim_data$beta_true,

        # Metadata
        sample_size = n,
        p_genes = p_genes,
        replicate_id = rep_id,
        true_model = true_model,
        method = "regular_BART",
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
