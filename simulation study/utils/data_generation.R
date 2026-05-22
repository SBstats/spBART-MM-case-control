################################################################################
# data_generation.R
#
# PURPOSE:
# Generate simulation data for comparing Regular BART (with DART sparsity prior)
# vs Semi-parametric BART (spBART) for high-dimensional binary outcomes.
#
# DATA GENERATING PROCESS:
# - Binary outcome Y via probit link: Y ~ Bernoulli(Phi(f(X, Z)))
# - Gene expression matrix X: p = 10,000 genes (10 true + 9,990 null)
# - Clinical covariates Z: age, sex, race, BMI
#
# TWO TRUE MODELS (Tree-based DGP with balanced gene signals):
#
# 1. True Model 1: 20 trees using BOTH genes and covariates (NO linear Z'beta)
#    - 10 depth-1 trees with pure marginal gene effects (±0.35 each)
#    - 10 depth-2 trees with gene-covariate interactions (±0.40 each, SYMMETRIC)
#    - Y ~ Bernoulli(Phi(g_1 + ... + g_20))  -- PURE SUM OF TREES
#    - spBART is MISSPECIFIED (cannot capture gene-cov interactions within trees)
#    - Each gene appears exactly 2 times (BALANCED signal)
#    - Interaction contribution ~53% of gene signal (strong misspecification test)
#
# 2. True Model 2: 15 trees using ONLY genes + linear covariate effects
#    - 10 depth-1 trees with pure marginal gene effects (±0.40 each)
#    - 5 depth-2 trees with gene-gene interactions (±0.35 each, SYMMETRIC)
#    - Linear covariate effects: Z'beta
#    - Y ~ Bernoulli(Phi(g_1 + ... + g_15 + Z'beta))
#    - spBART is CORRECTLY SPECIFIED (additive gene trees + linear covariates)
#    - Each gene appears exactly 2 times (BALANCED signal)
#
# CUTPOINT STRATEGY:
# - Cutpoints are FIXED population quantiles from N(0,1), NOT sample-dependent
# - This ensures the true function is fixed across simulations
# - Uses varied quantiles (35th-65th percentiles) to avoid oracle-optimal splits
#
# AUTHOR: Simulation Study
# DATE: January 2026
################################################################################

#' Generate Simulation Data with Tree-based DGP
#'
#' @param n Sample size
#' @param p_genes Total number of genes (500, 5000, or 10000)
#' @param n_true_genes Number of true signal genes (default: 10, fixed across all p)
#' @param true_model Which true model to use (1 or 2)
#' @param seed Random seed for reproducibility
#'
#' @return A list containing:
#'   - Y: Binary outcome vector
#'   - prob_Y: True P(Y=1) for each observation
#'   - X_genes: Gene expression matrix (n x p_genes)
#'   - Z_covariates: Clinical covariate matrix (n x 4)
#'   - true_gene_indices: Indices of true signal genes
#'   - true_gene_names: Names of true signal genes
#'   - beta_true: True covariate coefficients (for True Model 2)
#'   - tree_contributions: Matrix of individual tree contributions
#'   - n, p_genes, n_true_genes: Design parameters
#'
#' @note The true signal genes are always genes 1-10 regardless of p_genes value.
#'       This allows fair comparison across different dimensionality settings.
#'
generate_simulation_data <- function(n = 500,
                                     p_genes = 10000,
                                     n_true_genes = 10,
                                     true_model = 1,
                                     seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # ============================================================================
  # STEP 1: Generate Gene Expression Matrix X (Independent Genes)
  # ============================================================================
  # All genes ~ N(0, 1) independently
  X_genes <- matrix(rnorm(n * p_genes), nrow = n, ncol = p_genes)
  colnames(X_genes) <- paste0("gene_", 1:p_genes)

  # True signal genes are the first n_true_genes columns
  true_gene_indices <- 1:n_true_genes
  true_gene_names <- paste0("gene_", true_gene_indices)

  # Extract true genes for convenience
  X_true <- X_genes[, true_gene_indices, drop = FALSE]

  # ============================================================================
  # STEP 2: Generate Clinical Covariates Z
  # ============================================================================
  # z1 (age): Standardized, ~ N(0, 1)
  # z2 (sex): Binary, ~ Bernoulli(0.5)
  # z3 (race): Binary, ~ Bernoulli(0.3)
  # z4 (BMI): Standardized, ~ N(0, 1)

  z1_age <- rnorm(n, mean = 0, sd = 1)
  z2_sex <- rbinom(n, size = 1, prob = 0.5)
  z3_race <- rbinom(n, size = 1, prob = 0.3)
  z4_bmi <- rnorm(n, mean = 0, sd = 1)

  Z_covariates <- cbind(age = z1_age, sex = z2_sex, race = z3_race, bmi = z4_bmi)

  # ============================================================================
  # STEP 3: Define Tree Functions Based on True Model
  # ============================================================================
  # Cut points are defined relative to variable ranges (using quantiles)
  # Trees have depth 2-3 for moderate complexity

  # ==========================================================================
  # FIXED POPULATION QUANTILE CUTPOINTS
  # ==========================================================================
  # Use theoretical quantiles from N(0,1) rather than sample quantiles
  # This ensures the true function is FIXED across simulations
  # Different genes have different cutpoints for realism

  # Gene cutpoints: varied population quantiles (not all at median)
  # Gene 1:  45th percentile -> qnorm(0.45) = -0.126
  # Gene 2:  55th percentile -> qnorm(0.55) =  0.126
  # Gene 3:  40th percentile -> qnorm(0.40) = -0.253
  # Gene 4:  60th percentile -> qnorm(0.60) =  0.253
  # Gene 5:  50th percentile -> qnorm(0.50) =  0.000 (median)
  # Gene 6:  50th percentile -> qnorm(0.50) =  0.000 (median)
  # Gene 7:  35th percentile -> qnorm(0.35) = -0.385
  # Gene 8:  65th percentile -> qnorm(0.65) =  0.385
  # Gene 9:  42nd percentile -> qnorm(0.42) = -0.202
  # Gene 10: 58th percentile -> qnorm(0.58) =  0.202

  c_gene <- c(qnorm(0.45), qnorm(0.55), qnorm(0.40), qnorm(0.60),
              qnorm(0.50), qnorm(0.50), qnorm(0.35), qnorm(0.65),
              qnorm(0.42), qnorm(0.58))

  # Covariate cutpoints: also fixed population quantiles
  c_age <- qnorm(0.55)  # = 0.126 (45% "older")
  c_bmi <- qnorm(0.40)  # = -0.253 (60% "higher BMI")

  if (true_model == 1) {
    # ==========================================================================
    # TRUE MODEL 1: 20 Trees using BOTH genes AND covariates (PURE SUM OF TREES)
    # ==========================================================================
    # This model has trees that split on both genes and clinical covariates,
    # making spBART's assumption (genes-only in trees) MISSPECIFIED.
    # NOTE: NO linear Z'beta component - this is a pure sum-of-trees model.
    #
    # DESIGN (Jan 2026 - v4 with balanced signals):
    # - 10 depth-1 trees: One marginal tree per gene (equal baseline)
    # - 10 depth-2 trees: One gene-covariate interaction per gene
    # - NO linear covariate main effects (pure tree model)
    # - Each gene appears exactly twice: 1 marginal + 1 interaction
    # - Marginal terminal values: ±0.35 per tree
    # - Interaction terminal values: ±0.40 per tree (SYMMETRIC)
    # - Interaction contribution is ~53% of gene signal (stronger test)
    #
    # Gene appearance counts: ALL genes appear exactly 2 times (balanced)

    x1 <- X_true[, 1]; x2 <- X_true[, 2]; x3 <- X_true[, 3]
    x4 <- X_true[, 4]; x5 <- X_true[, 5]; x6 <- X_true[, 6]
    x7 <- X_true[, 7]; x8 <- X_true[, 8]; x9 <- X_true[, 9]
    x10 <- X_true[, 10]

    # =========================================================================
    # DEPTH-1 TREES: Pure Marginal Gene Effects (Trees 1-10)
    # =========================================================================
    # Each gene gets exactly one marginal tree with terminal value ±0.35
    # This provides baseline signal that spBART CAN capture

    g1 <- ifelse(x1 <= c_gene[1], -0.35, 0.35)   # Gene 1
    g2 <- ifelse(x2 <= c_gene[2], -0.35, 0.35)   # Gene 2
    g3 <- ifelse(x3 <= c_gene[3], -0.35, 0.35)   # Gene 3
    g4 <- ifelse(x4 <= c_gene[4], -0.35, 0.35)   # Gene 4
    g5 <- ifelse(x5 <= c_gene[5], -0.35, 0.35)   # Gene 5
    g6 <- ifelse(x6 <= c_gene[6], -0.35, 0.35)   # Gene 6
    g7 <- ifelse(x7 <= c_gene[7], -0.35, 0.35)   # Gene 7
    g8 <- ifelse(x8 <= c_gene[8], -0.35, 0.35)   # Gene 8
    g9 <- ifelse(x9 <= c_gene[9], -0.35, 0.35)   # Gene 9
    g10 <- ifelse(x10 <= c_gene[10], -0.35, 0.35) # Gene 10

    # =========================================================================
    # DEPTH-2 TREES: Gene-Covariate Interactions (Trees 11-20)
    # =========================================================================
    # Each gene gets exactly one interaction tree with terminal value ±0.40
    # This is the signal that spBART CANNOT capture (misspecification)
    # All interactions are SYMMETRIC (same magnitude for both directions)

    # Tree 11: Gene 1 x Age interaction
    g11 <- ifelse(x1 > c_gene[1] & z1_age > c_age, 0.40,
                  ifelse(x1 <= c_gene[1] & z1_age <= c_age, -0.40, 0))

    # Tree 12: Gene 2 x Sex interaction
    g12 <- ifelse(x2 > c_gene[2] & z2_sex == 1, 0.40,
                  ifelse(x2 <= c_gene[2] & z2_sex == 0, -0.40, 0))

    # Tree 13: Gene 3 x BMI interaction
    g13 <- ifelse(x3 > c_gene[3] & z4_bmi > c_bmi, 0.40,
                  ifelse(x3 <= c_gene[3] & z4_bmi <= c_bmi, -0.40, 0))

    # Tree 14: Gene 4 x Race interaction
    g14 <- ifelse(x4 > c_gene[4] & z3_race == 1, 0.40,
                  ifelse(x4 <= c_gene[4] & z3_race == 0, -0.40, 0))

    # Tree 15: Gene 5 x Age interaction
    g15 <- ifelse(x5 > c_gene[5] & z1_age > c_age, 0.40,
                  ifelse(x5 <= c_gene[5] & z1_age <= c_age, -0.40, 0))

    # Tree 16: Gene 6 x Sex interaction
    g16 <- ifelse(x6 > c_gene[6] & z2_sex == 1, 0.40,
                  ifelse(x6 <= c_gene[6] & z2_sex == 0, -0.40, 0))

    # Tree 17: Gene 7 x BMI interaction
    g17 <- ifelse(x7 > c_gene[7] & z4_bmi > c_bmi, 0.40,
                  ifelse(x7 <= c_gene[7] & z4_bmi <= c_bmi, -0.40, 0))

    # Tree 18: Gene 8 x Race interaction
    g18 <- ifelse(x8 > c_gene[8] & z3_race == 1, 0.40,
                  ifelse(x8 <= c_gene[8] & z3_race == 0, -0.40, 0))

    # Tree 19: Gene 9 x Age interaction
    g19 <- ifelse(x9 > c_gene[9] & z1_age > c_age, 0.40,
                  ifelse(x9 <= c_gene[9] & z1_age <= c_age, -0.40, 0))

    # Tree 20: Gene 10 x Sex interaction
    g20 <- ifelse(x10 > c_gene[10] & z2_sex == 1, 0.40,
                  ifelse(x10 <= c_gene[10] & z2_sex == 0, -0.40, 0))

    # Store tree contributions (all 20 trees)
    tree_contributions <- cbind(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
                                g11, g12, g13, g14, g15, g16, g17, g18,
                                g19, g20)

    # =========================================================================
    # PURE SUM OF TREES - NO LINEAR Z'BETA
    # =========================================================================
    # Model 1 is a pure sum-of-trees model with NO linear covariate component.
    # The covariates only appear through tree splits (gene-covariate interactions).
    # This tests whether spBART fails when covariates interact with genes in trees.
    f_true <- rowSums(tree_contributions)

    # Mark beta as NA to indicate no linear covariate effects in Model 1
    beta_true <- c(age = NA, sex = NA, race = NA, bmi = NA)

  } else if (true_model == 2) {
    # ==========================================================================
    # TRUE MODEL 2: 15 Trees using ONLY genes + Linear Covariate Effects
    # ==========================================================================
    # This model has:
    #   - Trees that split ONLY on genes (no covariate splits)
    #   - Additive linear covariate effects Z'beta
    #
    # This is the spBART model, so spBART is CORRECTLY SPECIFIED.
    #
    # f(X, Z) = g_1(X) + ... + g_15(X) + beta_1*Z_1 + ... + beta_4*Z_4
    #
    # DESIGN (Jan 2026 - v4 with balanced signals):
    # - 10 depth-1 trees: One marginal tree per gene (equal baseline)
    # - 5 depth-2 trees: Gene-gene interactions (maintains complexity)
    # - Each gene appears exactly 1-2 times for balanced signal
    # - Marginal terminal values: ±0.40 per tree
    # - Interaction terminal values: ±0.35 per tree (SYMMETRIC)
    #
    # Gene appearance counts: ALL genes appear 1-2 times (balanced)

    x1 <- X_true[, 1]; x2 <- X_true[, 2]; x3 <- X_true[, 3]
    x4 <- X_true[, 4]; x5 <- X_true[, 5]; x6 <- X_true[, 6]
    x7 <- X_true[, 7]; x8 <- X_true[, 8]; x9 <- X_true[, 9]
    x10 <- X_true[, 10]

    # =========================================================================
    # DEPTH-1 TREES: Pure Marginal Gene Effects (Trees 1-10)
    # =========================================================================
    # Each gene gets exactly one marginal tree with terminal value ±0.40

    g1 <- ifelse(x1 <= c_gene[1], -0.40, 0.40)   # Gene 1
    g2 <- ifelse(x2 <= c_gene[2], -0.40, 0.40)   # Gene 2
    g3 <- ifelse(x3 <= c_gene[3], -0.40, 0.40)   # Gene 3
    g4 <- ifelse(x4 <= c_gene[4], -0.40, 0.40)   # Gene 4
    g5 <- ifelse(x5 <= c_gene[5], -0.40, 0.40)   # Gene 5
    g6 <- ifelse(x6 <= c_gene[6], -0.40, 0.40)   # Gene 6
    g7 <- ifelse(x7 <= c_gene[7], -0.40, 0.40)   # Gene 7
    g8 <- ifelse(x8 <= c_gene[8], -0.40, 0.40)   # Gene 8
    g9 <- ifelse(x9 <= c_gene[9], -0.40, 0.40)   # Gene 9
    g10 <- ifelse(x10 <= c_gene[10], -0.40, 0.40) # Gene 10

    # =========================================================================
    # DEPTH-2 TREES: Gene-Gene Interactions (Trees 11-15)
    # =========================================================================
    # 5 interaction trees pairing genes
    # All interactions are SYMMETRIC (same magnitude for both directions)

    # Tree 11: Gene 1 x Gene 2 interaction
    g11 <- ifelse(x1 > c_gene[1] & x2 > c_gene[2], 0.35,
                  ifelse(x1 <= c_gene[1] & x2 <= c_gene[2], -0.35, 0))

    # Tree 12: Gene 3 x Gene 4 interaction
    g12 <- ifelse(x3 > c_gene[3] & x4 > c_gene[4], 0.35,
                  ifelse(x3 <= c_gene[3] & x4 <= c_gene[4], -0.35, 0))

    # Tree 13: Gene 5 x Gene 6 interaction
    g13 <- ifelse(x5 > c_gene[5] & x6 > c_gene[6], 0.35,
                  ifelse(x5 <= c_gene[5] & x6 <= c_gene[6], -0.35, 0))

    # Tree 14: Gene 7 x Gene 8 interaction
    g14 <- ifelse(x7 > c_gene[7] & x8 > c_gene[8], 0.35,
                  ifelse(x7 <= c_gene[7] & x8 <= c_gene[8], -0.35, 0))

    # Tree 15: Gene 9 x Gene 10 interaction
    g15 <- ifelse(x9 > c_gene[9] & x10 > c_gene[10], 0.35,
                  ifelse(x9 <= c_gene[9] & x10 <= c_gene[10], -0.35, 0))

    # Store tree contributions (all 15 trees)
    tree_contributions <- cbind(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
                                g11, g12, g13, g14, g15)

    # Sum of trees (gene component)
    f_genes <- rowSums(tree_contributions)

    # Linear covariate effects
    beta_true <- c(age = 0.50, sex = 0.60, race = -0.40, bmi = 0.45)
    Z_beta <- Z_covariates %*% beta_true

    # Total function
    f_true <- f_genes + as.vector(Z_beta)

  } else {
    stop("true_model must be 1 or 2")
  }

  # ============================================================================
  # STEP 4: Generate Binary Outcome Y
  # ============================================================================
  # Y ~ Bernoulli(Phi(f_true))
  # where Phi is the standard normal CDF (probit link)

  # Compute P(Y = 1 | X, Z) = Phi(f(X, Z))
  prob_Y <- pnorm(f_true)

  # Check for extreme probabilities - warn but DO NOT alter f_true
  # The true function must remain fixed to match the designed DGP
  min_prob <- min(prob_Y)
  max_prob <- max(prob_Y)

  if (min_prob > 0.90 || max_prob < 0.10) {
    warning(sprintf(
      "Extreme probability range detected: [%.3f, %.3f]. Consider adjusting tree coefficients.",
      min_prob, max_prob
    ))
  }

  # Generate binary outcome
  Y <- rbinom(n, size = 1, prob = prob_Y)

  # Safety check: ensure not all 0s or all 1s
  max_attempts <- 10
  attempt <- 1
  while ((sum(Y) == 0 || sum(Y) == n) && attempt <= max_attempts) {
    # Regenerate with adjusted probabilities
    Y <- rbinom(n, size = 1, prob = prob_Y)
    attempt <- attempt + 1
  }

  if (sum(Y) == 0 || sum(Y) == n) {
    warning("Generated Y has no variation. Consider adjusting tree coefficients.")
  }

  # ============================================================================
  # STEP 5: Return Results
  # ============================================================================

  result <- list(
    # Outcome
    Y = Y,
    prob_Y = prob_Y,
    f_true = f_true,

    # Tree contributions
    tree_contributions = tree_contributions,

    # Predictors
    X_genes = X_genes,
    Z_covariates = Z_covariates,

    # True gene information
    true_gene_indices = true_gene_indices,
    true_gene_names = true_gene_names,

    # Model parameters
    beta_true = beta_true,
    true_model = true_model,

    # Design parameters
    n = n,
    p_genes = p_genes,
    n_true_genes = n_true_genes
  )

  return(result)
}


#' Summarize Generated Data
#'
#' @param data Output from generate_simulation_data()
#'
#' @return Prints summary statistics
#'
summarize_simulation_data <- function(data) {

  cat("\n")
  cat("================================================================================\n")
  cat(sprintf("  Simulation Data Summary (True Model %d)\n", data$true_model))
  cat("================================================================================\n\n")

  # Sample size and dimensions
  cat(sprintf("Sample size: n = %d\n", data$n))
  cat(sprintf("Gene matrix: %d x %d\n", nrow(data$X_genes), ncol(data$X_genes)))
  cat(sprintf("True genes: %d (%s)\n", data$n_true_genes,
              paste(data$true_gene_names, collapse = ", ")))

  # Class balance
  n_Y1 <- sum(data$Y == 1)
  n_Y0 <- sum(data$Y == 0)
  prop_Y1 <- mean(data$Y)
  cat(sprintf("\nClass balance:\n"))
  cat(sprintf("  Y = 1: %d (%.1f%%)\n", n_Y1, 100 * prop_Y1))
  cat(sprintf("  Y = 0: %d (%.1f%%)\n", n_Y0, 100 * (1 - prop_Y1)))

  # True function f
  cat(sprintf("\nTrue function f(X, Z):\n"))
  cat(sprintf("  Mean: %.3f\n", mean(data$f_true)))
  cat(sprintf("  SD: %.3f\n", sd(data$f_true)))
  cat(sprintf("  Range: [%.3f, %.3f]\n", min(data$f_true), max(data$f_true)))

  # True probabilities
  cat(sprintf("\nTrue probabilities P(Y=1 | X, Z):\n"))
  cat(sprintf("  Mean: %.3f\n", mean(data$prob_Y)))
  cat(sprintf("  SD: %.3f\n", sd(data$prob_Y)))
  cat(sprintf("  Range: [%.3f, %.3f]\n", min(data$prob_Y), max(data$prob_Y)))

  # Tree contributions
  cat(sprintf("\nTree contributions (mean absolute value):\n"))
  tree_means <- colMeans(abs(data$tree_contributions))
  for (i in 1:length(tree_means)) {
    cat(sprintf("  g_%d: %.3f\n", i, tree_means[i]))
  }

  # Covariates
  cat(sprintf("\nCovariates:\n"))
  cat(sprintf("  age (z1): mean = %.3f, SD = %.3f\n",
              mean(data$Z_covariates[, "age"]), sd(data$Z_covariates[, "age"])))
  cat(sprintf("  sex (z2): prop = %.3f\n", mean(data$Z_covariates[, "sex"])))
  cat(sprintf("  race (z3): prop = %.3f\n", mean(data$Z_covariates[, "race"])))
  cat(sprintf("  bmi (z4): mean = %.3f, SD = %.3f\n",
              mean(data$Z_covariates[, "bmi"]), sd(data$Z_covariates[, "bmi"])))

  # True beta (for model 2 only)
  if (data$true_model == 2) {
    cat(sprintf("\nTrue beta coefficients (linear covariate effects):\n"))
    cat(sprintf("  %s\n", paste(names(data$beta_true), "=",
                                 round(data$beta_true, 3), collapse = ", ")))
  } else {
    cat("\nNote: Model 1 is a PURE SUM OF TREES with no linear Z'beta component.\n")
    cat("      Covariates only enter through gene-covariate interactions in trees.\n")
  }

  cat("\n")
}


#' Check Class Balance Across Multiple Seeds
#'
#' @param n_seeds Number of seeds to test
#' @param n Sample size
#' @param true_model True model (1 or 2)
#' @param base_seed Base seed for generating seed sequence
#'
#' @return Data frame with class balance statistics
#'
check_class_balance <- function(n_seeds = 20, n = 500, true_model = 1,
                                base_seed = 1000) {

  class_props <- numeric(n_seeds)
  all_zeros <- 0
  all_ones <- 0

  for (s in 1:n_seeds) {
    seed <- base_seed * n + s
    data_s <- generate_simulation_data(n = n, true_model = true_model, seed = seed)
    class_props[s] <- mean(data_s$Y)

    if (sum(data_s$Y) == 0) all_zeros <- all_zeros + 1
    if (sum(data_s$Y) == n) all_ones <- all_ones + 1
  }

  result <- data.frame(
    true_model = true_model,
    n = n,
    n_seeds = n_seeds,
    mean_prop_Y1 = mean(class_props),
    sd_prop_Y1 = sd(class_props),
    min_prop_Y1 = min(class_props),
    max_prop_Y1 = max(class_props),
    min_class_count = round(min(pmin(class_props, 1 - class_props)) * n),
    all_zeros = all_zeros,
    all_ones = all_ones
  )

  cat(sprintf("\nClass Balance Check (True Model %d, n = %d, %d seeds):\n",
              true_model, n, n_seeds))
  cat(sprintf("  P(Y=1): mean = %.3f, SD = %.3f, range = [%.3f, %.3f]\n",
              result$mean_prop_Y1, result$sd_prop_Y1,
              result$min_prop_Y1, result$max_prop_Y1))
  cat(sprintf("  Min class count: %d\n", result$min_class_count))
  cat(sprintf("  Datasets with all 0s: %d, all 1s: %d\n", all_zeros, all_ones))

  return(result)
}


# ==============================================================================
# SECTION: Initial Gene Screening Functions (OPTIONAL - Not Used in Main Simulation)
# ==============================================================================
# NOTE: These screening functions are NOT used in the main simulation study.
# The simulation passes all p genes directly to BART without pre-screening.
# These functions are retained for potential alternative analyses.
# ==============================================================================

#' Probit regression p-value for a single gene
#'
#' @param y Binary outcome vector
#' @param x Gene expression vector
#' @param Z Covariate matrix
#' @return p-value for gene coefficient
probit_pvalue <- function(y, x, Z) {
  df <- data.frame(y = y, gene = x, Z)
  tryCatch({
    fit <- suppressWarnings(
      glm(y ~ ., data = df, family = binomial(link = "probit"),
          control = glm.control(maxit = 100))
    )
    summary(fit)$coefficients["gene", "Pr(>|z|)"]
  }, error = function(e) return(1))
}

#' Probit regression coefficient for a single gene
#'
#' @param y Binary outcome vector
#' @param x Gene expression vector
#' @param Z Covariate matrix
#' @return Coefficient for gene
probit_coef <- function(y, x, Z) {
  df <- data.frame(y = y, gene = x, Z)
  tryCatch({
    fit <- suppressWarnings(
      glm(y ~ ., data = df, family = binomial(link = "probit"),
          control = glm.control(maxit = 100))
    )
    gene_coef <- coef(fit)["gene"]
    if (is.na(gene_coef)) return(0)
    gene_coef
  }, error = function(e) return(0))
}

#' Initial gene screening: probit p-value < 0.05
#'
#' Single-step screening procedure:
#' Fit univariate probit model for each gene, keep genes with p-value < 0.05
#'
#' @param X_genes Gene expression matrix (n x p)
#' @param Z Covariate matrix (n x q)
#' @param Y Binary outcome vector
#' @return List containing screened genes, p-values, coefficients, and screening info
initial_gene_screening <- function(X_genes, Z, Y) {

  p_genes <- ncol(X_genes)
  gene_names <- colnames(X_genes)

  # ==========================================================================
  # Fit univariate probit models and filter by p < 0.05
  # ==========================================================================
  cat("    Fitting univariate probit models for", p_genes, "genes...\n")

  pvals <- sapply(1:p_genes, function(j) probit_pvalue(Y, X_genes[, j], Z))
  coefs <- sapply(1:p_genes, function(j) probit_coef(Y, X_genes[, j], Z))
  names(pvals) <- gene_names
  names(coefs) <- gene_names

  # Filter by p < 0.05
  screened_genes <- names(pvals)[pvals < 0.05]
  n_screened <- length(screened_genes)

  cat("    Screening result:", n_screened, "genes with p-value < 0.05\n")

  return(list(
    screened_genes = screened_genes,
    n_screened = n_screened,
    pvalues = pvals[screened_genes],
    coefficients = coefs[screened_genes]
  ))
}


#' FDR-controlled gene selection from PIPs
#'
#' @param PIPs Named vector of posterior inclusion probabilities
#' @param alpha FDR control level (default: 0.05)
#' @return List with selected genes, R_star, tau, and FDR_hat
fdr_gene_selection <- function(PIPs, alpha = 0.05) {
  PIP_sorted <- sort(PIPs, decreasing = TRUE)
  D <- length(PIP_sorted)

  if (D == 0) {
    return(list(selected_genes = character(0), R_star = 0, tau = NA, FDR_hat = NA))
  }

  # Compute FDR at each rank
  FDR <- numeric(D)
  for (r in 1:D) {
    FDR[r] <- sum(1 - PIP_sorted[1:r]) / r
  }

  # Find optimal R*
  valid_R <- which(FDR <= alpha)

  if (length(valid_R) > 0) {
    R_star <- max(valid_R)
    selected_genes <- names(PIP_sorted)[1:R_star]
    tau <- PIP_sorted[R_star]
    FDR_hat <- FDR[R_star]
  } else {
    R_star <- 0
    selected_genes <- character(0)
    tau <- NA
    FDR_hat <- NA
  }

  return(list(
    selected_genes = selected_genes,
    R_star = R_star,
    tau = tau,
    FDR_hat = FDR_hat
  ))
}


#' Visualize Tree Structure (for documentation)
#'
#' @param true_model Which model (1 or 2)
#'
#' @return Prints tree structure description
#'
describe_trees <- function(true_model = 1) {

  cat("\n")
  n_trees <- ifelse(true_model == 1, 20, 15)
  cat("================================================================================\n")
  cat(sprintf("  Tree Structure for True Model %d (%d Trees)\n", true_model, n_trees))
  cat("================================================================================\n\n")

  cat("CUTPOINT STRATEGY: Fixed population quantiles from N(0,1)\n")
  cat("  Gene 1:  45th percentile -> qnorm(0.45) = -0.126\n")
  cat("  Gene 2:  55th percentile -> qnorm(0.55) =  0.126\n")
  cat("  Gene 3:  40th percentile -> qnorm(0.40) = -0.253\n")
  cat("  Gene 4:  60th percentile -> qnorm(0.60) =  0.253\n")
  cat("  Gene 5:  50th percentile -> qnorm(0.50) =  0.000\n")
  cat("  Gene 6:  50th percentile -> qnorm(0.50) =  0.000\n")
  cat("  Gene 7:  35th percentile -> qnorm(0.35) = -0.385\n")
  cat("  Gene 8:  65th percentile -> qnorm(0.65) =  0.385\n")
  cat("  Gene 9:  42nd percentile -> qnorm(0.42) = -0.202\n")
  cat("  Gene 10: 58th percentile -> qnorm(0.58) =  0.202\n")
  cat("  Covariates: age at 55th %ile (0.126), BMI at 40th %ile (-0.253)\n\n")

  if (true_model == 1) {
    cat("Model 1: Trees use BOTH genes AND covariates (PURE SUM OF TREES)\n")
    cat("         (spBART MISSPECIFIED - cannot capture gene-cov interactions)\n")
    cat("Design: 10 marginal trees + 10 gene-covariate interaction trees (NO Z'beta)\n")
    cat("Terminal values: ±0.35 per marginal, ±0.40 per interaction (SYMMETRIC)\n")
    cat("Interaction contribution: ~53% of gene signal (strong misspecification test)\n\n")

    cat("--- DEPTH-1 TREES: Pure Marginal Effects (Trees 1-10) ---\n\n")

    cat("Tree g_1:  gene_1 marginal [±0.35] (cut at -0.126)\n")
    cat("Tree g_2:  gene_2 marginal [±0.35] (cut at  0.126)\n")
    cat("Tree g_3:  gene_3 marginal [±0.35] (cut at -0.253)\n")
    cat("Tree g_4:  gene_4 marginal [±0.35] (cut at  0.253)\n")
    cat("Tree g_5:  gene_5 marginal [±0.35] (cut at  0.000)\n")
    cat("Tree g_6:  gene_6 marginal [±0.35] (cut at  0.000)\n")
    cat("Tree g_7:  gene_7 marginal [±0.35] (cut at -0.385)\n")
    cat("Tree g_8:  gene_8 marginal [±0.35] (cut at  0.385)\n")
    cat("Tree g_9:  gene_9 marginal [±0.35] (cut at -0.202)\n")
    cat("Tree g_10: gene_10 marginal [±0.35] (cut at  0.202)\n\n")

    cat("--- DEPTH-2 TREES: Gene-Covariate Interactions (Trees 11-20) ---\n\n")

    cat("Tree g_11: gene_1 x age  [±0.40] (symmetric)\n")
    cat("Tree g_12: gene_2 x sex  [±0.40] (symmetric)\n")
    cat("Tree g_13: gene_3 x BMI  [±0.40] (symmetric)\n")
    cat("Tree g_14: gene_4 x race [±0.40] (symmetric)\n")
    cat("Tree g_15: gene_5 x age  [±0.40] (symmetric)\n")
    cat("Tree g_16: gene_6 x sex  [±0.40] (symmetric)\n")
    cat("Tree g_17: gene_7 x BMI  [±0.40] (symmetric)\n")
    cat("Tree g_18: gene_8 x race [±0.40] (symmetric)\n")
    cat("Tree g_19: gene_9 x age  [±0.40] (symmetric)\n")
    cat("Tree g_20: gene_10 x sex [±0.40] (symmetric)\n\n")

    cat("Gene appearance summary:\n")
    cat("  ALL genes: exactly 2 trees (1 marginal + 1 interaction) - BALANCED\n\n")

    cat("Covariate effects:\n")
    cat("  NO linear Z'beta component - covariates only appear via tree splits\n\n")

  } else {
    cat("Model 2: Trees use ONLY genes + linear covariates (spBART CORRECT)\n")
    cat("Design: 10 marginal trees + 5 gene-gene interaction trees + Z'beta\n")
    cat("Terminal values: ±0.40 per marginal, ±0.35 per interaction (SYMMETRIC)\n\n")

    cat("--- DEPTH-1 TREES: Pure Marginal Effects (Trees 1-10) ---\n\n")

    cat("Tree g_1:  gene_1 marginal [±0.40] (cut at -0.126)\n")
    cat("Tree g_2:  gene_2 marginal [±0.40] (cut at  0.126)\n")
    cat("Tree g_3:  gene_3 marginal [±0.40] (cut at -0.253)\n")
    cat("Tree g_4:  gene_4 marginal [±0.40] (cut at  0.253)\n")
    cat("Tree g_5:  gene_5 marginal [±0.40] (cut at  0.000)\n")
    cat("Tree g_6:  gene_6 marginal [±0.40] (cut at  0.000)\n")
    cat("Tree g_7:  gene_7 marginal [±0.40] (cut at -0.385)\n")
    cat("Tree g_8:  gene_8 marginal [±0.40] (cut at  0.385)\n")
    cat("Tree g_9:  gene_9 marginal [±0.40] (cut at -0.202)\n")
    cat("Tree g_10: gene_10 marginal [±0.40] (cut at  0.202)\n\n")

    cat("--- DEPTH-2 TREES: Gene-Gene Interactions (Trees 11-15) ---\n\n")

    cat("Tree g_11: gene_1 x gene_2 [±0.35] (symmetric)\n")
    cat("Tree g_12: gene_3 x gene_4 [±0.35] (symmetric)\n")
    cat("Tree g_13: gene_5 x gene_6 [±0.35] (symmetric)\n")
    cat("Tree g_14: gene_7 x gene_8 [±0.35] (symmetric)\n")
    cat("Tree g_15: gene_9 x gene_10 [±0.35] (symmetric)\n\n")

    cat("Gene appearance summary:\n")
    cat("  ALL genes: exactly 2 trees (1 marginal + 1 interaction) - BALANCED\n\n")

    cat("Linear covariate effects:\n")
    cat("  beta = (0.50, 0.60, -0.40, 0.45) for (age, sex, race, BMI)\n\n")
  }

  cat("================================================================================\n")
}
