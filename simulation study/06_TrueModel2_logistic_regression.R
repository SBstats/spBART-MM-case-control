################################################################################
# 06_TrueModel2_logistic_regression.R
#
# PURPOSE:
# Fit parametric regression models (logistic AND probit) on data generated
# from True Model 2 to evaluate recovery of true covariate coefficients (beta).
#
# TRUE MODEL 2: Semi-parametric (Additive Covariate Effects)
# f(X, Z) = f(X) + Z'beta, where f(X) includes gene-gene interactions
# but NO gene-covariate interactions.
# Data generated with PROBIT link: Y ~ Bernoulli(Phi(f_trees(X) + Z'beta))
#
# ANALYSIS PIPELINE:
# 1. Generate data from True Model 2 (probit link)
# 2. Fit TWO models:
#    a) Logistic regression: Y ~ X_genes + Z_covariates (MISSPECIFIED link)
#    b) Probit regression: Y ~ X_genes + Z_covariates (CORRECT link, wrong gene structure)
# 3. Extract beta estimates for covariates and compute:
#    - Bias (beta_hat - beta_true)
#    - MSE/RMSE
#    - Sign correctness
#
# SIMULATION DESIGN:
# - Sample sizes: n = 1500, 2500
# - Number of genes: p = 500, 2000, 3000
# - True signal genes: 10 (fixed across all p values)
#
# SLURM ARRAY JOB:
# - Each array task (1-500) runs one replicate for ALL n,p combinations
#
# OUTPUTS:
# - results/TrueModel2/n{1500,2500}/p{500,2000,3000}/logistic_reg/rep_{001-500}.RDS
# - results/TrueModel2/n{1500,2500}/p{500,2000,3000}/probit_reg/rep_{001-500}.RDS
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
cat("  06_TrueModel2_logistic_regression.R: Parametric Regression Models\n")
cat("  Fits: Logistic Regression AND Probit Regression\n")
cat("  True Model 2: Additive covariate effects (check beta recovery)\n")
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

required_packages <- c("MASS")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("Packages loaded successfully!\n\n")

# ==============================================================================
# SECTION 3: Simulation Settings
# ==============================================================================

sample_sizes <- c(1500, 2500)
gene_counts <- c(500, 2000, 3000)
true_model <- 2  # Only True Model 2 (has true linear covariate effects)

cat("================================================================================\n")
cat("  Simulation Settings\n")
cat("================================================================================\n")
cat("Replicate ID:", rep_id, "\n")
cat("True Model:", true_model, "(additive covariate effects)\n")
cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
cat("Gene counts:", paste(gene_counts, collapse = ", "), "\n")
cat("Methods: Logistic Regression + Probit Regression\n\n")

# ==============================================================================
# SECTION 4: Run Replicate for Each (n, p) Combination
# ==============================================================================

for (n in sample_sizes) {
  for (p_genes in gene_counts) {

    cat("\n")
    cat("################################################################################\n")
    cat("  Running replicate", rep_id, "for n =", n, ", p =", p_genes, "\n")
    cat("################################################################################\n\n")

    # Create results directories for both methods
    results_dir_logistic <- file.path(script_dir, "results", "TrueModel2",
                                      paste0("n", n), paste0("p", p_genes), "logistic_reg")
    results_dir_probit <- file.path(script_dir, "results", "TrueModel2",
                                    paste0("n", n), paste0("p", p_genes), "probit_reg")

    if (!dir.exists(results_dir_logistic)) dir.create(results_dir_logistic, recursive = TRUE)
    if (!dir.exists(results_dir_probit)) dir.create(results_dir_probit, recursive = TRUE)

    output_file_logistic <- file.path(results_dir_logistic, sprintf("rep_%03d.RDS", rep_id))
    output_file_probit <- file.path(results_dir_probit, sprintf("rep_%03d.RDS", rep_id))

    # Check if both already completed
    if (file.exists(output_file_logistic) && file.exists(output_file_probit)) {
      cat("Replicate", rep_id, "already completed for n =", n, ", p =", p_genes, ". Skipping.\n")
      next
    }

    # Seed depends on n, p, and rep_id for reproducibility
    # Use same seed as other TrueModel2 scripts for consistency
    set.seed(2000 * n + p_genes + rep_id)
    rep_start <- Sys.time()

    tryCatch({
      # ========================================================================
      # Step 1: Generate data
      # ========================================================================
      cat("Step 1: Generating simulated data...\n")
      sim_data <- generate_simulation_data(
        n = n, p_genes = p_genes, n_true_genes = 10,
        true_model = true_model, seed = 2000 * n + p_genes + rep_id
      )

      X_genes <- sim_data$X_genes
      Z_cov <- sim_data$Z_covariates
      Y <- sim_data$Y
      beta_true <- sim_data$beta_true

      gene_names <- colnames(X_genes)
      cov_names <- colnames(Z_cov)

      cat("  n =", n, ", p =", ncol(X_genes), ", Y=1:", sum(Y), "\n")
      cat("  True beta:", paste(round(beta_true, 3), collapse = ", "), "\n")

      # ========================================================================
      # Step 2: Fit BOTH parametric models (logistic and probit)
      # ========================================================================
      # Define models to fit
      models_to_fit <- list(
        list(name = "logistic", link = "logit", output_file = output_file_logistic),
        list(name = "probit", link = "probit", output_file = output_file_probit)
      )

      for (model_info in models_to_fit) {
        model_name <- model_info$name
        link_function <- model_info$link
        output_file <- model_info$output_file

        # Skip if already completed
        if (file.exists(output_file)) {
          cat("\n", toupper(model_name), "regression already completed. Skipping.\n", sep = "")
          next
        }

        cat("\n")
        cat("========================================================================\n")
        cat("  Fitting", toupper(model_name), "regression\n")
        cat("========================================================================\n")

        # Combined design matrix: covariates FIRST, then genes
        # This ensures covariates are estimated even when p >> n,
        # because GLM aliases variables from the rightmost columns
        # when the design matrix is rank-deficient.
        X_full <- cbind(Z_cov, X_genes)

        # Create data frame for glm
        # IMPORTANT: check.names=FALSE preserves original column names
        df <- data.frame(Y = Y, X_full, check.names = FALSE)

        # Diagnostic: Print covariate column names in dataframe
        cat("  Covariate columns in dataframe:", paste(cov_names, collapse = ", "), "\n")
        cat("  Checking if covariates exist in df:", all(cov_names %in% colnames(df)), "\n")

        # Fit regression with specified link function
        # Use glm.control to increase max iterations for high-dimensional cases
        fit <- tryCatch({
          glm(Y ~ ., data = df, family = binomial(link = link_function),
              control = glm.control(maxit = 100, epsilon = 1e-8))
        }, error = function(e) {
          cat("  ERROR: GLM fitting failed -", conditionMessage(e), "\n")
          return(NULL)
        }, warning = function(w) {
          cat("  WARNING: GLM fitting -", conditionMessage(w), "\n")
          glm(Y ~ ., data = df, family = binomial(link = link_function),
              control = glm.control(maxit = 100, epsilon = 1e-8))
        })

        if (is.null(fit)) {
          cat("  Skipping", model_name, "regression due to GLM failure.\n")
          next  # Skip to next model
        }

        # Check convergence
        converged <- fit$converged
        cat("  GLM converged:", converged, "\n")

        # Check for separation: if fitted probabilities are all 0 or 1,
        # the model has perfect/quasi-perfect separation and coefficients
        # will be astronomically large and unreliable
        fitted_probs <- fitted(fit)
        separation_detected <- all(fitted_probs < 1e-10 | fitted_probs > 1 - 1e-10)
        if (!separation_detected) {
          # Also check if any non-NA covariate coefficient is extremely large
          coef_all <- coef(fit)
          cov_coef_values <- coef_all[cov_names]
          cov_coef_values <- cov_coef_values[!is.na(cov_coef_values)]
          if (length(cov_coef_values) > 0 && any(abs(cov_coef_values) > 100)) {
            separation_detected <- TRUE
          }
        }
        cat("  Separation detected:", separation_detected, "\n")

        # ======================================================================
        # Extract covariate coefficient estimates
        # ======================================================================
        cat("  Extracting covariate coefficient estimates...\n")

        # Get coefficient estimates - use coef() which includes NAs for aliased coefficients
        coef_full <- coef(fit)  # This includes NAs for aliased variables

        # Get summary for standard errors, z-values, p-values (only for non-aliased)
        coef_summary <- summary(fit)$coefficients

        # Diagnostic: Check coefficient vectors
        cat("  Full coefficient vector length:", length(coef_full), "\n")
        cat("  Coefficient summary matrix dimensions:", nrow(coef_summary), "x", ncol(coef_summary), "\n")
        cat("  Expected number of coefficients:", 1 + p_genes + 4, "\n")

        # Extract coefficients BY NAME (robust to aliasing/dropped coefficients)
        # coef_full includes NAs for aliased variables
        # coef_summary only includes non-aliased variables
        coef_full_names <- names(coef_full)
        coef_summary_names <- rownames(coef_summary)

        # Debug: Check if covariates are in coefficient names
        cat("  Checking for covariates in full coefficient vector...\n")
        covariates_in_full <- sapply(cov_names, function(cov) cov %in% coef_full_names)
        cat("  Covariates in coef_full:", paste(names(which(covariates_in_full)), collapse = ", "), "\n")

        # Additional diagnostic: show last few coefficient names from full vector
        if (length(coef_full_names) > 0) {
          cat("  Last 10 variables in coef_full:", paste(tail(coef_full_names, 10), collapse = ", "), "\n")
        }

        # Initialize vectors with NAs
        beta_hat <- rep(NA_real_, length(cov_names))
        beta_se <- rep(NA_real_, length(cov_names))
        beta_z <- rep(NA_real_, length(cov_names))
        beta_pval <- rep(NA_real_, length(cov_names))
        names(beta_hat) <- cov_names
        names(beta_se) <- cov_names
        names(beta_z) <- cov_names
        names(beta_pval) <- cov_names

        # Extract covariate coefficients that are present
        for (cov in cov_names) {
          # Try to get estimate from coef_full (includes NAs for aliased)
          if (cov %in% coef_full_names) {
            beta_hat[cov] <- coef_full[cov]

            # Get SE, z, p-value from summary (only if not aliased)
            if (cov %in% coef_summary_names) {
              beta_se[cov] <- coef_summary[cov, "Std. Error"]
              beta_z[cov] <- coef_summary[cov, "z value"]
              beta_pval[cov] <- coef_summary[cov, "Pr(>|z|)"]
            } else {
              cat("  NOTE: Covariate", cov, "has coefficient but no SE/p-value (aliased)\n")
            }
          } else {
            cat("  WARNING: Covariate", cov, "not found in coefficient vector (dropped)\n")
          }
        }

        # Check if all covariates were successfully extracted
        n_missing_covariates <- sum(is.na(beta_hat))
        all_covariates_present <- n_missing_covariates == 0

        if (!all_covariates_present) {
          cat("  WARNING:", n_missing_covariates, "out of", length(cov_names),
              "covariates were aliased/dropped by GLM\n")
        }

        cat("  Beta estimates:", paste(round(beta_hat, 3), collapse = ", "), "\n")

        # ======================================================================
        # Compute beta estimation metrics
        # ======================================================================
        cat("  Computing beta estimation metrics...\n")

        # Bias (per covariate) - will be NA for missing coefficients
        beta_bias <- beta_hat - beta_true

        # MSE (per covariate, not aggregated)
        beta_mse <- beta_bias^2

        # RMSE (per covariate)
        beta_rmse <- abs(beta_bias)  # For single replicate, RMSE = |bias|

        # Sign agreement (correct sign indicator)
        # For each covariate, check if estimated sign matches true sign
        # Will be NA for missing coefficients
        beta_sign_correct <- sign(beta_hat) == sign(beta_true)

        # Statistical significance (p-value < 0.05)
        # For each covariate, check if p-value is significant
        # Will be NA for missing coefficients
        beta_significant <- beta_pval < 0.05
        n_significant_covariates <- sum(beta_significant, na.rm = TRUE)

        cat("  Beta bias:", paste(round(beta_bias, 4), collapse = ", "), "\n")
        cat("  Beta MSE (per covariate):", paste(round(beta_mse, 4), collapse = ", "), "\n")
        cat("  Beta RMSE (per covariate):", paste(round(beta_rmse, 4), collapse = ", "), "\n")
        cat("  Beta sign correct:", paste(as.integer(beta_sign_correct), collapse = ", "), "\n")
        cat("  Beta significant (p<0.05):", paste(as.integer(beta_significant), collapse = ", "), "\n")
        cat("  Number of significant covariates:", n_significant_covariates, "out of", length(cov_names), "\n")

        # ======================================================================
        # Compute prediction metrics (optional, for comparison)
        # ======================================================================
        cat("  Computing prediction metrics...\n")

        # Predicted probabilities
        prob_hat <- predict(fit, newdata = df, type = "response")

        # Compute basic prediction metrics
        pred_class <- ifelse(prob_hat > 0.5, 1, 0)
        accuracy <- mean(pred_class == Y)

        # Brier Score
        brier_score <- mean((Y - prob_hat)^2)

        cat("  Accuracy:", round(accuracy, 3), "\n")
        cat("  Brier Score:", round(brier_score, 4), "\n")

        # ======================================================================
        # Save results
        # ======================================================================
        model_end <- Sys.time()
        model_runtime <- difftime(model_end, rep_start, units = "mins")

        result <- list(
          # Beta estimation metrics (per covariate)
          beta_hat = beta_hat,
          beta_se = beta_se,
          beta_z = beta_z,
          beta_pval = beta_pval,
          beta_bias = beta_bias,          # Per covariate
          beta_mse = beta_mse,             # Per covariate (not aggregated)
          beta_rmse = beta_rmse,           # Per covariate (not aggregated)
          beta_sign_correct = beta_sign_correct,
          beta_significant = beta_significant,  # Per covariate (p < 0.05)
          n_significant_covariates = n_significant_covariates,  # Count of significant covariates

          # Prediction metrics
          accuracy = accuracy,
          brier_score = brier_score,

          # Model diagnostics
          converged = converged,
          separation_detected = separation_detected,
          all_covariates_present = all_covariates_present,
          n_missing_covariates = n_missing_covariates,

          # Ground truth
          beta_true = beta_true,

          # Metadata
          sample_size = n,
          p_genes = p_genes,
          replicate_id = rep_id,
          true_model = true_model,
          method = paste0(model_name, "_regression"),
          link_function = link_function,
          runtime_mins = as.numeric(model_runtime),
          seed = 2000 * n + p_genes + rep_id
        )

        saveRDS(result, file = output_file)
        cat("  Saved results to:", basename(output_file), "\n")
        cat("  Runtime:", round(as.numeric(model_runtime), 2), "minutes\n")
      }  # End of models loop

      cat("\nCompleted n =", n, ", p =", p_genes, "\n")

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
