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
# PREREQUISITE: Run 00_setup_BART_package.R first
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

cat("\nLoading required packages...\n")

required_packages <- c("MASS", "truncnorm", "caret", "pROC")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

library(BART)

# Verify modified BART package
required_functions <- c("wbart_create", "wbart_update_response",
                        "wbart_run_iteration", "wbart_destroy")
missing_functions <- required_functions[!sapply(required_functions, exists)]
if (length(missing_functions) > 0) {
  stop("Missing BART functions: ", paste(missing_functions, collapse = ", "),
       "\nPlease run 00_setup_BART_package.R first.")
}

cat("Packages loaded successfully!\n\n")

# ==============================================================================
# SECTION 3: spBART Gibbs Sampler
# ==============================================================================

#' Semi-parametric Probit BART Gibbs Sampler
#' Consistent with real data analysis (Dec05) version
#'
#' Gibbs sampler:
#'   1. U | Y, f, β  (Albert-Chib latent variable augmentation)
#'   2. β | U, f     (Conjugate Gaussian)
#'   3. f | U, β     (BART with stateful updates)
fit_semiparametric_probit_DART <- function(X, Z, Y,
                                           X_test = NULL, Z_test = NULL,
                                           n_burn = 2000, n_iter = 1000, n_thin = 5,
                                           n_trees = 200, sparse = TRUE,
                                           theta = 0, omega = 1,
                                           a = 0.5, b = 1, rho = NULL,
                                           numcut = 100, usequants = FALSE,
                                           sigma_beta = sqrt(10)) {

  N <- length(Y)
  D <- ncol(X)
  J <- ncol(Z)

  if (is.null(rho)) rho <- D

  X <- as.matrix(X)
  Z <- as.matrix(Z)
  Y <- as.integer(Y)

  # Handle test data
  has_test <- !is.null(X_test) && !is.null(Z_test)
  if (has_test) {
    X_test <- as.matrix(X_test)
    Z_test <- as.matrix(Z_test)
    N_test <- nrow(X_test)
  }

  # Initialize parameters
  beta <- rep(0, J)
  U <- rnorm(N)
  f_train <- rep(0, N)

  # Storage for posterior draws
  beta_draws <- matrix(NA, nrow = n_iter, ncol = J)
  colnames(beta_draws) <- colnames(Z)
  f_draws <- matrix(NA, nrow = n_iter, ncol = N)
  varcount <- matrix(0, nrow = n_iter, ncol = D)
  colnames(varcount) <- colnames(X)

  # Storage for test predictions
  if (has_test) {
    f_test_draws <- matrix(NA, nrow = n_iter, ncol = N_test)
    prob_test_draws <- matrix(NA, nrow = n_iter, ncol = N_test)
  }

  # Initialize BART sampler
  Y_adjusted_init <- U - Z %*% beta

  bart_sampler <- wbart_create(
    x.train = X,
    y.train = Y_adjusted_init,
    x.test = if (has_test) X_test else X[1, , drop = FALSE],
    sparse = sparse,
    theta = theta,
    omega = omega,
    a = a, b = b, rho = rho,
    augment = FALSE,
    ntree = n_trees,
    numcut = numcut,
    usequants = usequants
  )

  n_total <- n_burn + (n_iter * n_thin)

  for (iter in 1:n_total) {

    if (iter %% 500 == 0) {
      message("  Iteration ", iter, "/", n_total)
    }

    # Step 1: Update latent variables U | Y, f, β (Albert-Chib)
    mu_U <- f_train + Z %*% beta
    for (i in 1:N) {
      if (Y[i] == 1) {
        U[i] <- rtruncnorm(1, a = 0, b = Inf, mean = mu_U[i], sd = 1)
      } else {
        U[i] <- rtruncnorm(1, a = -Inf, b = 0, mean = mu_U[i], sd = 1)
      }
    }

    # Step 2: Update β | U, f (Conjugate Gaussian)
    residuals <- U - f_train
    precision_prior <- diag(J) / (sigma_beta^2)
    precision_post <- t(Z) %*% Z + precision_prior
    V_beta <- solve(precision_post)
    mean_beta <- V_beta %*% (t(Z) %*% residuals)
    beta <- MASS::mvrnorm(n = 1, mu = mean_beta, Sigma = V_beta)

    # Step 3: Update DART function f | U, β (STATEFUL UPDATE)
    Y_adjusted <- as.numeric(U - Z %*% beta)
    wbart_update_response(bart_sampler, Y_adjusted)
    bart_result <- wbart_run_iteration(bart_sampler)
    f_train <- bart_result$yhat.train

    if (has_test) {
      f_test <- bart_result$yhat.test
    }

    # Store posterior draws (after burn-in, with thinning)
    if (iter > n_burn && (iter - n_burn) %% n_thin == 0) {
      post_idx <- (iter - n_burn) / n_thin
      beta_draws[post_idx, ] <- beta
      f_draws[post_idx, ] <- f_train
      varcount[post_idx, ] <- bart_result$varcount

      if (has_test) {
        f_test_draws[post_idx, ] <- f_test
        mu_test <- f_test + as.numeric(Z_test %*% beta)
        prob_test_draws[post_idx, ] <- pnorm(mu_test)
      }
    }
  }

  wbart_destroy(bart_sampler)
  gc()  # Force garbage collection after destroying BART sampler

  # Compute posterior mean of latent function
  f_hat <- colMeans(f_draws)
  beta_hat <- colMeans(beta_draws)
  mu_hat <- f_hat + as.numeric(Z %*% beta_hat)

  result <- list(
    beta_draws = beta_draws,
    f_train_draws = f_draws,
    varcount = varcount,
    n_burn = n_burn,
    n_iter = n_iter,
    n_thin = n_thin,
    f_hat = f_hat,
    beta_hat = beta_hat,
    mu_hat = mu_hat
  )

  if (has_test) {
    result$f_test_draws <- f_test_draws
    result$prob.test <- prob_test_draws  # Match pbart() output format
  }

  return(result)
}

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
sigma_beta <- sqrt(10)

alpha_FDR <- 0.05
K_folds <- 5
true_model <- 1

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

    results_dir <- file.path(script_dir, "results", "TrueModel1",
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

        fit_k <- fit_semiparametric_probit_DART(
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
        final_fit <- fit_semiparametric_probit_DART(
          X = X_union, Z = Z_cov, Y = Y,
          n_burn = n_burn, n_iter = n_iter, n_thin = n_thin,
          n_trees = ntrees, sparse = sparse,
          a = a, b = b, rho = length(union_genes),
          sigma_beta = sigma_beta
        )

        # Predicted probabilities
        mu_hat <- final_fit$mu_hat
        prob_hat <- pnorm(mu_hat)

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

        # RMSE of latent function
        final_RMSE <- sqrt(mean((mu_hat - f_true)^2))

        # Beta estimation (descriptive only — DGM 1 has no true beta, so no
        # truth-comparison metrics are computed; see writeup Sec 7.1)
        final_beta_draws <- final_fit$beta_draws
        final_beta_hat <- colMeans(final_beta_draws)
        final_beta_sd <- apply(final_beta_draws, 2, sd)
        final_beta_ci_lower <- apply(final_beta_draws, 2, quantile, 0.025)
        final_beta_ci_upper <- apply(final_beta_draws, 2, quantile, 0.975)

      } else {
        # No genes selected — fill final metrics with NAs and skip the final fit.
        # (We deliberately do not run a covariate-only spBART fallback here.)
        cat("  WARNING: No genes selected. Skipping final-fit; saving NA metrics.\n")

        final_AUC   <- NA_real_
        final_Brier <- NA_real_
        final_RMSE  <- NA_real_

        J <- ncol(Z_cov)
        na_vec <- setNames(rep(NA_real_, J), colnames(Z_cov))
        final_beta_hat      <- na_vec
        final_beta_sd       <- na_vec
        final_beta_ci_lower <- na_vec
        final_beta_ci_upper <- na_vec
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
        final_RMSE = final_RMSE,

        # Beta estimation from final model (descriptive only — DGM 1 has no
        # true beta, so bias/RMSE/coverage/sign-prob are not computed)
        beta_hat = final_beta_hat,
        beta_sd = final_beta_sd,
        beta_ci_lower = final_beta_ci_lower,
        beta_ci_upper = final_beta_ci_upper,

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
