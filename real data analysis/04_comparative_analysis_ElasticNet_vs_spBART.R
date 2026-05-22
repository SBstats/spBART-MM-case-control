################################################################################
# Comparative Analysis: Elastic Net vs. Semi-Parametric Probit BART (spBART)
################################################################################
#
# PURPOSE:
# Compare variable selection between Elastic Net and the proposed spBART model
# by examining gene selection at different FDR thresholds.
#
# METHODOLOGY:
# 1. Run Elastic Net analysis to get reference gene set
#    - Using pipeline from 02_naive_analysis...Dec05 file
#    - Screening with age, sex_binary ONLY
#    - Final model with all 4 clinical covariates (age, sex, race, bmi)
#    - study_indicator in penalized X (like genes)
#    - BMI cutoff >= 25 (matches writeup Table 8.4)
#
# 2. Run spBART analysis to get PIPs and FDR curve
#    - Using pipeline from 03_real_data_analysis...Dec05 file
#    - 5-fold CV with union selection criterion
#    - study_indicator in BART X (like genes)
#    - Clinical covariates (Z) with age, sex during CV
#    - Final model with all 4 clinical covariates
#
# 3. Compute spBART gene selection at multiple FDR alpha levels
# 4. Find the alpha level where spBART matches Elastic Net gene count
# 5. Visualize comparison and overlap statistics
#
# SEEDS (matching original files for consistency):
# - Data partitioning: seed 789 (test-first sampling)
# - Elastic Net replications: seed 789 + rep (rep = 1:100)
# - Elastic Net final model: seed 456
# - spBART CV fold creation: seed 890
#
# OUTPUTS:
# - comparative_analysis_results.RDS (all results)
# - gene_selection_comparison.pdf (main comparison plot)
# - comparison_table.pdf and comparison_table.tex (detailed table)
# - overlap_statistics.pdf and overlap_statistics.tex
#
################################################################################

# ==============================================================================
# SECTION 1: Setup and Load Dependencies
# ==============================================================================

script_start_time <- Sys.time()
cat("\n")
cat("================================================================================\n")
cat("  Comparative Analysis: Elastic Net vs. spBART\n")
cat("  Dec 2024 Version - Binary BMI (>=25), Study Indicator in BART\n")
cat("================================================================================\n")
cat("Script started at:", format(script_start_time), "\n")

# Determine script location (consistent with spBART analysis files)
script_path <- tryCatch({
  # Method 1: Works when running via Rscript
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    normalizePath(dirname(sub("--file=", "", file_arg)))
  } else {
    NULL
  }
}, error = function(e) NULL)

# Try alternative methods if Method 1 failed
if (is.null(script_path) || length(script_path) == 0) {
  script_path <- tryCatch({
    # Method 2: Works when sourcing the script
    dirname(sys.frame(1)$ofile)
  }, error = function(e) NULL)
}

if (is.null(script_path) || length(script_path) == 0) {
  script_path <- tryCatch({
    # Method 3: Works in RStudio
    if (requireNamespace("rstudioapi", quietly = TRUE)) {
      dirname(rstudioapi::getActiveDocumentContext()$path)
    } else {
      NULL
    }
  }, error = function(e) NULL)
}

# Method 4: Fallback to current working directory
if (is.null(script_path) || length(script_path) == 0 || script_path == "") {
  script_path <- getwd()
  cat("Note: Using current working directory as script location\n")
}

cat("Script location:", script_path, "\n")

# Set working directory to script location (ensures data/ paths work correctly)
setwd(script_path)
cat("Working directory set to:", getwd(), "\n")

# Create results directory
todays_date <- format(Sys.Date(), "%Y_%m_%d")
results_dir <- file.path(script_path, paste0("results_comparative_analysis_", todays_date))

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
  cat("Created results directory:", results_dir, "\n")
} else {
  cat("Using existing results directory:", results_dir, "\n")
}

result_path <- function(filename) {
  file.path(results_dir, filename)
}

# Load required packages
cat("\nLoading required packages...\n")

library(dplyr)
library(pROC)
library(MASS)
library(glmnet)
library(caret)
library(mclust)
library(readxl)
library(data.table)
library(truncnorm)

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
library(DESeq2)

# Install modified BART package
cat("\n=== Installing Modified BART Package ===\n")

if ("BART" %in% rownames(installed.packages())) {
  cat("  Removing existing BART package...\n")
  remove.packages("BART")
}

# The package should be a tarball called "BART_2.9.9.tar.gz" in the same directory as this script
bart_path <- file.path(script_path, "BART_2.9.9.tar.gz")

if (!file.exists(bart_path)) {
  stop(paste0("ERROR: Modified BART package not found at: ", bart_path,
              "\n  Please ensure the BART_2.9.9.tar.gz file is in the same directory as this script."))
}

cat(sprintf("  Installing modified BART package from: %s\n", bart_path))

install.packages(bart_path, repos = NULL, type = "source", INSTALL_opts = "--no-multiarch")
library(BART)

required_functions <- c("wbart_create", "wbart_update_response", "wbart_run_iteration", "wbart_destroy")
missing_functions <- required_functions[!sapply(required_functions, exists)]
if (length(missing_functions) > 0) {
  stop(paste0("ERROR: Missing BART functions: ", paste(missing_functions, collapse = ", ")))
}

cat("All packages loaded successfully!\n\n")


# ==============================================================================
# SECTION 2: Data Preprocessing (Common to Both Methods)
# ==============================================================================

cat("================================================================================\n")
cat("  SECTION 2: Data Preprocessing\n")
cat("================================================================================\n\n")

# Load cohort data
errc_seq_id_cohort <- as.data.frame(read_excel("data/ERRC.ID.xlsx"))
genebody_data_cohort <- as.data.frame(readRDS("data/genebody_377.RDS"))
load("data/patient_metadata_full_n797.RData")
patient_metadata_cohort <- as.data.frame(analysis_data)

patient_metadata_cohort_subset <- patient_metadata_cohort %>%
  filter(errcid %in% errc_seq_id_cohort$ERRC.ID)

cohort_workdf <- data.frame(
  errcid = patient_metadata_cohort_subset$errcid,
  age = patient_metadata_cohort_subset$age_diag,
  race = patient_metadata_cohort_subset$race_composite,
  sex = patient_metadata_cohort_subset$sex_emr,
  bmi = patient_metadata_cohort_subset$bmi_dx_emr,
  mm_types = patient_metadata_cohort_subset$dx_errc
)
cohort_workdf <- na.omit(cohort_workdf)

cohort_workdf <- cohort_workdf[cohort_workdf$mm_types == "Multiple Myeloma", ]
cohort_workdf$mm_status <- rep("CASE", nrow(cohort_workdf))

cohort_workdf <- cohort_workdf[cohort_workdf$race == "White" | cohort_workdf$race == "Black/African-American", ]
cohort_workdf$race_binary <- ifelse(cohort_workdf$race == "White", 1, 0)
# BMI cutoff >= 25 (matches writeup Table 8.4)
cohort_workdf$bmi_binary <- ifelse(cohort_workdf$bmi >= 25, 1, 0)
cohort_workdf$sex_binary <- ifelse(cohort_workdf$sex == "M", 1, 0)
cohort_workdf$PooledData_ID <- 1:nrow(cohort_workdf)

seq_ID_cohort_final <- errc_seq_id_cohort$Sequencing.ID[errc_seq_id_cohort$ERRC.ID %in% cohort_workdf$errcid]
cohort_workdf$sequencing_ID <- seq_ID_cohort_final

genebody_data_cohort_unfiltered_unnormalized <- genebody_data_cohort[,
  colnames(genebody_data_cohort) %in% seq_ID_cohort_final]

seq_to_pooled_map <- setNames(cohort_workdf$PooledData_ID, cohort_workdf$sequencing_ID)
colnames(genebody_data_cohort_unfiltered_unnormalized) <- seq_to_pooled_map[colnames(genebody_data_cohort_unfiltered_unnormalized)]
genebody_data_cohort_unfiltered_unnormalized <- genebody_data_cohort_unfiltered_unnormalized[,
  order(as.numeric(colnames(genebody_data_cohort_unfiltered_unnormalized)))]

message("Cohort data loaded: ", nrow(cohort_workdf), " patients")

# Load case-control data
genebody_data_case_control <- as.data.frame(readRDS("data/Canada_case_control_genebody_count.RDS"))
patient_metadata_case_control <- readRDS("data/MM_Questionnaire_BC_for_BART_projects.RDS")
Canada_case_control_sample_key <- as.data.frame(read.csv("data/Canada_case_control_sample_key.csv"))

Canada_case_control_key <- data.frame(
  Study.ID = Canada_case_control_sample_key$Study.ID,
  Assigned.ID = Canada_case_control_sample_key$Assigned.ID
)

patient_metadata_case_control_subset <- patient_metadata_case_control %>%
  filter(StudyID %in% Canada_case_control_key$Study.ID)

HtCurrIn_clean <- ifelse(is.na(patient_metadata_case_control_subset$HtCurrIn), 0,
                         patient_metadata_case_control_subset$HtCurrIn)
HtCurrentInches <- patient_metadata_case_control_subset$HtCurrFt * 12 + HtCurrIn_clean
BMI_current <- ifelse(is.na(patient_metadata_case_control_subset$HtCurrFt) |
                        is.na(patient_metadata_case_control_subset$WtCurrLb),
                      NA,
                      (patient_metadata_case_control_subset$WtCurrLb / (HtCurrentInches^2)) * 703)
patient_metadata_case_control_subset$BMI_current <- BMI_current

caseControl_workdf <- data.frame(
  StudyID = patient_metadata_case_control_subset$StudyID,
  age = patient_metadata_case_control_subset$Age,
  race = patient_metadata_case_control_subset$EthSelf,
  sex = patient_metadata_case_control_subset$Sex,
  bmi = patient_metadata_case_control_subset$BMI_current,
  mm_status = patient_metadata_case_control_subset$MM_STATUS
)
caseControl_workdf <- na.omit(caseControl_workdf)

caseControl_workdf <- caseControl_workdf[caseControl_workdf$race == 1 | caseControl_workdf$race == 5, ]
caseControl_workdf$race_binary <- ifelse(caseControl_workdf$race == 5, 0, caseControl_workdf$race)
# BMI cutoff >= 25 (matches writeup Table 8.4)
caseControl_workdf$bmi_binary <- ifelse(caseControl_workdf$bmi >= 25, 1, 0)
caseControl_workdf$sex_binary <- ifelse(caseControl_workdf$sex == "M", 1, 0)
caseControl_workdf$PooledData_ID <- (nrow(cohort_workdf) + 1):(nrow(cohort_workdf) + nrow(caseControl_workdf))

seq_ID_caseControl_final <- Canada_case_control_key$Assigned.ID[Canada_case_control_key$Study.ID %in% caseControl_workdf$StudyID]
caseControl_workdf$sequencing_ID <- seq_ID_caseControl_final

genebody_data_caseControl_unfiltered_unnormalized <- genebody_data_case_control[,
  colnames(genebody_data_case_control) %in% seq_ID_caseControl_final]

seq_to_pooled_map_caseControl <- setNames(caseControl_workdf$PooledData_ID, caseControl_workdf$sequencing_ID)
colnames(genebody_data_caseControl_unfiltered_unnormalized) <- seq_to_pooled_map_caseControl[colnames(genebody_data_caseControl_unfiltered_unnormalized)]
genebody_data_caseControl_unfiltered_unnormalized <- genebody_data_caseControl_unfiltered_unnormalized[,
  order(as.numeric(colnames(genebody_data_caseControl_unfiltered_unnormalized)))]

message("Case-control data loaded: ", nrow(caseControl_workdf), " patients")

# Pool metadata
cohort_selected <- cohort_workdf %>%
  dplyr::select(PooledData_ID, race_binary, bmi_binary, sex_binary, age, mm_status)
# study_indicator coding: 0 = UCMM (UChicago cohort), 1 = BC (British Columbia case-control)
cohort_selected$study_indicator <- 0

caseControl_selected <- caseControl_workdf %>%
  dplyr::select(PooledData_ID, race_binary, bmi_binary, sex_binary, age, mm_status)
# study_indicator coding: 0 = UCMM (UChicago cohort), 1 = BC (British Columbia case-control)
caseControl_selected$study_indicator <- 1

pooled_metadata <- rbind(cohort_selected, caseControl_selected) %>%
  arrange(PooledData_ID)
pooled_metadata$outcome <- ifelse(pooled_metadata$mm_status == "CASE", 1, 0)

# Pool genebody data
pooled_genebody_data_unfiltered_unnormalized <- cbind(
  genebody_data_cohort_unfiltered_unnormalized,
  genebody_data_caseControl_unfiltered_unnormalized
)

message("Pooled data: ", nrow(pooled_metadata), " patients, ",
        nrow(pooled_genebody_data_unfiltered_unnormalized), " genes")

# Gene filtering
n_samples <- ncol(pooled_genebody_data_unfiltered_unnormalized)
threshold_pct <- 0.05
threshold_samples <- ceiling(n_samples * threshold_pct)
low_count_per_gene <- rowSums(pooled_genebody_data_unfiltered_unnormalized < 10)
genes_to_keep <- low_count_per_gene <= threshold_samples
pooled_genebody_data_filtered_unnormalized <- pooled_genebody_data_unfiltered_unnormalized[genes_to_keep, ]

message("Genes after filtering: ", nrow(pooled_genebody_data_filtered_unnormalized))

# DESeq2 normalization
colData <- data.frame(
  PooledData_ID = colnames(pooled_genebody_data_filtered_unnormalized),
  row.names = colnames(pooled_genebody_data_filtered_unnormalized)
)

dds <- DESeqDataSetFromMatrix(
  countData = round(pooled_genebody_data_filtered_unnormalized),
  colData = colData,
  design = ~1
)

vsd <- vst(dds, blind = FALSE)
pooled_genebody_data_filtered_normalized <- as.data.frame(assay(vsd))

message("Normalization complete: ", nrow(pooled_genebody_data_filtered_normalized), " genes x ",
        ncol(pooled_genebody_data_filtered_normalized), " samples\n")


# ==============================================================================
# SECTION 3: Elastic Net Analysis
# ==============================================================================

cat("================================================================================\n")
cat("  SECTION 3: Elastic Net Analysis\n")
cat("================================================================================\n\n")

# Use SAME data split as original files (seed 789, sample TEST first)
set.seed(789)

# Data partitioning (CONSISTENT with Dec05 files)
outcome_vec <- pooled_metadata$outcome
n_total <- nrow(pooled_metadata)
N_CV <- 500
N_test <- n_total - N_CV

outcome_1_ids <- pooled_metadata$PooledData_ID[outcome_vec == 1]
outcome_0_ids <- pooled_metadata$PooledData_ID[outcome_vec == 0]

# Stratified sampling: sample TEST set first (same as original files)
prop_cases <- sum(outcome_vec == 1) / n_total
n_test_cases <- round(N_test * prop_cases)
n_test_controls <- N_test - n_test_cases

test_ids_cases <- sample(outcome_1_ids, n_test_cases, replace = FALSE)
test_ids_controls <- sample(outcome_0_ids, n_test_controls, replace = FALSE)
val_ids <- c(test_ids_cases, test_ids_controls)

# Remaining IDs go to development/CV pool
dev_ids <- setdiff(pooled_metadata$PooledData_ID, val_ids)

message("Data split (consistent with Dec05 files, seed=789):")
message("  Development set: ", length(dev_ids), " patients")
message("  Test/Validation set: ", length(val_ids), " patients")

# Prepare data matrices
dev_indices <- which(pooled_metadata$PooledData_ID %in% dev_ids)
X_dev_genes <- t(pooled_genebody_data_filtered_normalized[, dev_indices])
y_dev <- pooled_metadata$outcome[dev_indices]

# Clinical covariates for screening: age, sex_binary ONLY (matching Dec05 file)
Z_dev_screening <- as.matrix(pooled_metadata[dev_indices, c("age", "sex_binary")])
# Clinical covariates for final model: all 4 (age, sex, race, bmi)
Z_dev_full <- as.matrix(pooled_metadata[dev_indices, c("age", "sex_binary", "race_binary", "bmi_binary")])

# Extract study_indicator separately to add to X matrix
study_indicator_dev <- pooled_metadata$study_indicator[dev_indices]

colnames(X_dev_genes) <- rownames(pooled_genebody_data_filtered_normalized)

# Initial screening: logistic regression + GMM (matching Dec05 file)
message("\nElastic Net - Stage 1: Logistic regression screening (age, sex_binary only)...")

logistic_pvalue <- function(y, x, Z) {
  df <- data.frame(y = y, gene = x, Z)
  tryCatch({
    fit <- suppressWarnings(glm(y ~ ., data = df, family = binomial(link = "logit")))
    summary(fit)$coefficients["gene", "Pr(>|z|)"]
  }, error = function(e) return(1))
}

pvals <- apply(X_dev_genes, 2, function(x) logistic_pvalue(y_dev, x, Z_dev_screening))
names(pvals) <- colnames(X_dev_genes)

dev_coefs <- sapply(colnames(X_dev_genes), function(gene) {
  df <- data.frame(y = y_dev, gene = X_dev_genes[, gene], Z_dev_screening)
  tryCatch({
    fit <- suppressWarnings(glm(y ~ ., data = df, family = binomial(link = "logit")))
    coef(fit)["gene"]
  }, error = function(e) return(0))
})
names(dev_coefs) <- colnames(X_dev_genes)

significant_genes_dev <- names(pvals)[pvals < 0.05]
message("  Genes with p < 0.05: ", length(significant_genes_dev))

# GMM clustering
# NOTE: Mclust uses hierarchical clustering for initialization, which is deterministic
# But we set seed here to ensure any internal randomness is controlled
if (length(significant_genes_dev) >= 10) {
  dev_coefs_significant <- dev_coefs[significant_genes_dev]
  dev_coefs_significant <- dev_coefs_significant[!is.na(dev_coefs_significant)]
  abs_coefs_dev <- abs(dev_coefs_significant)

  gmm_dev <- Mclust(abs_coefs_dev, G = 2:4, modelNames = "V", verbose = FALSE)

  if (!is.null(gmm_dev)) {
    cluster_assignments_dev <- gmm_dev$classification
    cluster_means_dev <- sapply(1:gmm_dev$G, function(k) mean(abs_coefs_dev[cluster_assignments_dev == k]))
    strongest_cluster_dev <- which.max(cluster_means_dev)
    candidate_genes <- names(dev_coefs_significant)[cluster_assignments_dev == strongest_cluster_dev]
    message("  GMM clusters: ", gmm_dev$G, ", strongest cluster: ", length(candidate_genes), " genes")
  } else {
    median_coef_dev <- median(abs(dev_coefs_significant))
    candidate_genes <- names(dev_coefs_significant)[abs(dev_coefs_significant) > median_coef_dev]
  }
} else {
  candidate_genes <- significant_genes_dev
}

screened_genes <- candidate_genes
message("  Screened genes for elastic net: ", length(screened_genes))
message("  (Compare with original file: should match 02_naive_analysis...Dec05 screened genes count)")

# Elastic net with 100 replications (matching Dec05 file)
message("\nElastic Net - Stage 2: Stability selection (100 replications)...")

# Add study_indicator to the penalized X matrix (along with genes)
# Using Z_dev_screening (age, sex_binary) as unpenalized covariates
X_dev_screened <- cbind(Z_dev_screening,
                        study_indicator = study_indicator_dev,
                        X_dev_genes[, screened_genes, drop = FALSE])

# Penalty factors: 0 for clinical covariates (unpenalized), 1 for study_indicator and genes (penalized)
penalty_factors_cv <- c(
  rep(0, ncol(Z_dev_screening)),       # age, sex_binary: NOT penalized
  1,                                   # study_indicator: PENALIZED (like genes)
  rep(1, length(screened_genes))       # Genes: PENALIZED
)

message("  Penalized variables: study_indicator + ", length(screened_genes), " genes")
message("  Unpenalized variables: ", paste(colnames(Z_dev_screening), collapse = ", "))

n_iterations <- 100
# Track selection of study_indicator AND genes
gene_selection_matrix <- matrix(0, nrow = n_iterations, ncol = length(screened_genes) + 1)
colnames(gene_selection_matrix) <- c("study_indicator", screened_genes)

for (rep in 1:n_iterations) {
  if (rep %% 20 == 0) message("  Replication ", rep, "/", n_iterations)

  # Seed matching Dec05 file: 789 + rep
  set.seed(789 + rep)

  cv_enet <- cv.glmnet(
    x = X_dev_screened,
    y = y_dev,
    family = "binomial",
    alpha = 0.5,
    penalty.factor = penalty_factors_cv,
    nfolds = 10,
    type.measure = "deviance"
  )

  coefs <- coef(cv_enet, s = "lambda.1se")[-1, ]
  selected_vars <- names(coefs)[coefs != 0]

  # Track study_indicator selection
  if ("study_indicator" %in% selected_vars) {
    gene_selection_matrix[rep, "study_indicator"] <- 1
  }

  # Track gene selection
  selected_genes_rep <- selected_vars[selected_vars %in% screened_genes]
  if (length(selected_genes_rep) > 0) {
    gene_selection_matrix[rep, selected_genes_rep] <- 1
  }
}

gene_selection_freq <- colSums(gene_selection_matrix)
selection_threshold <- 50

# Report study_indicator selection frequency
message("  study_indicator selected in ", gene_selection_freq["study_indicator"], "/", n_iterations, " replications")

# Get genes (excluding study_indicator) selected in >=50 replications
gene_freq_only <- gene_selection_freq[names(gene_selection_freq) != "study_indicator"]
enet_selected_genes_stability <- names(gene_freq_only)[gene_freq_only >= selection_threshold]

message("  Genes selected in >=50 replications: ", length(enet_selected_genes_stability))

# Final prognostic model (matching Dec05 file)
message("\nElastic Net - Stage 3: Final prognostic model (all 4 clinical covariates)...")

if (length(enet_selected_genes_stability) > 0) {
  # Include study_indicator as penalized variable in final model
  # Using Z_dev_full (all 4 clinical covariates) as unpenalized
  X_dev_final <- cbind(Z_dev_full,
                       study_indicator = study_indicator_dev,
                       X_dev_genes[, enet_selected_genes_stability, drop = FALSE])

  penalty_factors_final <- c(
    rep(0, ncol(Z_dev_full)),              # age, sex, race, bmi: NOT penalized
    1,                                      # study_indicator: PENALIZED
    rep(1, length(enet_selected_genes_stability))  # Genes: PENALIZED
  )

  # Seed matching Dec05 file: 456
  set.seed(456)
  prognostic_model_cv <- cv.glmnet(
    x = X_dev_final,
    y = y_dev,
    family = "binomial",
    alpha = 0.5,
    penalty.factor = penalty_factors_final,
    nfolds = 10,
    type.measure = "deviance"
  )

  coef_final <- coef(prognostic_model_cv, s = "lambda.1se")[-1, ]

  # Check if study_indicator was selected in final model
  study_indicator_selected <- coef_final["study_indicator"] != 0
  message("  study_indicator selected in final model: ", study_indicator_selected)

  # Get final genes (excluding study_indicator)
  # This matches original file's significant_genes_only computation
  enet_final_genes <- names(coef_final)[coef_final != 0 & names(coef_final) %in% enet_selected_genes_stability]

  # Diagnostic: Report stability selection count vs final model count
  message("  Genes from stability selection (>=50 reps): ", length(enet_selected_genes_stability))
  message("  Genes with non-zero coef in final model: ", length(enet_final_genes))

  # Check for genes that were shrunk to zero in final model
  genes_shrunk_to_zero <- setdiff(enet_selected_genes_stability, enet_final_genes)
  if (length(genes_shrunk_to_zero) > 0) {
    message("  Genes shrunk to zero in final model: ", length(genes_shrunk_to_zero))
    message("    ", paste(genes_shrunk_to_zero, collapse = ", "))
  }
} else {
  enet_final_genes <- character(0)
  message("  WARNING: No genes selected by Elastic Net")
}

# Store Elastic Net results
n_enet <- length(enet_final_genes)
message("\n=== ELASTIC NET REFERENCE: ", n_enet, " genes selected ===\n")


# ==============================================================================
# SECTION 4: spBART Analysis
# ==============================================================================

cat("================================================================================\n")
cat("  SECTION 4: spBART Analysis\n")
cat("================================================================================\n\n")

# spBART Gibbs sampler function
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

  beta <- rep(0, J)
  U <- rnorm(N)
  f_train <- rep(0, N)

  beta_draws <- matrix(NA, nrow = n_iter, ncol = J)
  f_draws <- matrix(NA, nrow = n_iter, ncol = N)
  varcount <- matrix(0, nrow = n_iter, ncol = D)

  if (!is.null(X_test)) {
    N_test <- nrow(X_test)
    f_test_draws <- matrix(NA, nrow = n_iter, ncol = N_test)
    prob_test_draws <- matrix(NA, nrow = n_iter, ncol = N_test)
  }

  Y_adjusted_init <- U - Z %*% beta

  bart_sampler <- wbart_create(
    x.train = X,
    y.train = Y_adjusted_init,
    x.test = X_test,
    sparse = sparse,
    theta = theta, omega = omega, a = a, b = b, rho = rho,
    augment = FALSE, ntree = n_trees, numcut = numcut, usequants = usequants
  )

  n_total <- n_burn + (n_iter * n_thin)

  for (iter in 1:n_total) {
    if (iter %% 1000 == 0) message("  Iteration ", iter, "/", n_total)

    mu_U <- f_train + Z %*% beta
    for (i in 1:N) {
      if (Y[i] == 1) {
        U[i] <- rtruncnorm(1, a = 0, b = Inf, mean = mu_U[i], sd = 1)
      } else {
        U[i] <- rtruncnorm(1, a = -Inf, b = 0, mean = mu_U[i], sd = 1)
      }
    }

    residuals <- U - f_train
    precision_prior <- diag(J) / (sigma_beta^2)
    precision_post <- t(Z) %*% Z + precision_prior
    V_beta <- solve(precision_post)
    mean_beta <- V_beta %*% (t(Z) %*% residuals)
    beta <- MASS::mvrnorm(n = 1, mu = mean_beta, Sigma = V_beta)

    Y_adjusted <- as.numeric(U - Z %*% beta)
    wbart_update_response(bart_sampler, Y_adjusted)
    bart_result <- wbart_run_iteration(bart_sampler)

    f_train <- bart_result$yhat.train
    if (!is.null(X_test)) f_test <- bart_result$yhat.test

    if (iter > n_burn && (iter - n_burn) %% n_thin == 0) {
      post_idx <- (iter - n_burn) / n_thin
      beta_draws[post_idx, ] <- beta
      f_draws[post_idx, ] <- f_train
      varcount[post_idx, ] <- bart_result$varcount

      if (!is.null(X_test)) {
        f_test_draws[post_idx, ] <- f_test
        mu_test <- f_test + Z_test %*% beta
        prob_test_draws[post_idx, ] <- pnorm(mu_test)
      }
    }
  }

  wbart_destroy(bart_sampler)

  result <- list(
    beta_draws = beta_draws,
    f_train_draws = f_draws,
    varcount = varcount,
    n_burn = n_burn, n_iter = n_iter, n_thin = n_thin
  )

  if (!is.null(X_test)) {
    result$f_test_draws <- f_test_draws
    result$prob.test <- prob_test_draws
  }

  return(result)
}

# Data partitioning for spBART (use SAME split as elastic net for fair comparison)
# CRITICAL: Reset seed to 789 to match the RNG state of the original file
# The Elastic Net section above consumed random numbers, so we must reset
# to ensure screening/GMM produces identical results to the original file
set.seed(789)

# CRITICAL: Consume the same random numbers as original file to match RNG state for Mclust
# Original file has stratified sampling at lines 1104-1105:
#   test_ids_cases <- sample(outcome_1_ids, n_test_cases, replace = FALSE)
#   test_ids_controls <- sample(outcome_0_ids, n_test_controls, replace = FALSE)
# We must replicate these sample() calls to match the RNG state before Mclust
# Use same variables and sizes as original to ensure identical RNG consumption
outcome_vec_temp <- pooled_metadata$outcome
outcome_1_ids_temp <- pooled_metadata$PooledData_ID[outcome_vec_temp == 1]
outcome_0_ids_temp <- pooled_metadata$PooledData_ID[outcome_vec_temp == 0]
N_total_temp <- nrow(pooled_metadata)
N_test_temp <- N_total_temp - 500  # Same as original
prop_cases_temp <- sum(outcome_vec_temp == 1) / N_total_temp
n_test_cases_temp <- round(N_test_temp * prop_cases_temp)
n_test_controls_temp <- N_test_temp - n_test_cases_temp

# These sample() calls consume random numbers to match original file's RNG state
dummy_test_cases <- sample(outcome_1_ids_temp, n_test_cases_temp, replace = FALSE)
dummy_test_controls <- sample(outcome_0_ids_temp, n_test_controls_temp, replace = FALSE)
# The dummy variables are not used - they just synchronize RNG state with original file

cv_ids <- dev_ids
test_ids <- val_ids

cv_indices <- which(pooled_metadata$PooledData_ID %in% cv_ids)
test_indices <- which(pooled_metadata$PooledData_ID %in% test_ids)

# Prepare data with study_indicator in X (for BART)
X_genes_cv <- t(pooled_genebody_data_filtered_normalized[, cv_indices])
X_genes_test <- t(pooled_genebody_data_filtered_normalized[, test_indices])

study_indicator_cv <- pooled_metadata$study_indicator[cv_indices]
study_indicator_test <- pooled_metadata$study_indicator[test_indices]

X_cv <- cbind(X_genes_cv, study_indicator = study_indicator_cv)
X_test <- cbind(X_genes_test, study_indicator = study_indicator_test)

# Clinical covariates for CV: sex_binary, age ONLY (MATCHING ORIGINAL Dec05 file EXACTLY)
# Original file line 1137: clinical_vars_cv <- c("sex_binary", "age")
clinical_vars_cv <- c("sex_binary", "age")
Z_cv <- as.matrix(pooled_metadata[cv_indices, clinical_vars_cv])
Z_test <- as.matrix(pooled_metadata[test_indices, clinical_vars_cv])

# Clinical covariates for final model: all 4 (MATCHING ORIGINAL Dec05 file EXACTLY)
# Original file line 1138: clinical_vars_full <- c("race_binary", "bmi_binary", "sex_binary", "age")
clinical_vars_full <- c("race_binary", "bmi_binary", "sex_binary", "age")
Z_cv_full <- as.matrix(pooled_metadata[cv_indices, clinical_vars_full])
Z_test_full <- as.matrix(pooled_metadata[test_indices, clinical_vars_full])

Y_cv <- pooled_metadata$outcome[cv_indices]
Y_test <- pooled_metadata$outcome[test_indices]

message("spBART - CV set: ", length(Y_cv), " patients")
message("spBART - Test set: ", length(Y_test), " patients")

# Initial screening for spBART (same approach as elastic net)
message("\nspBART - Initial gene screening (age, sex_binary only)...")

probit_pvalue <- function(y, x, Z) {
  df <- data.frame(y = y, gene = x, Z)
  tryCatch({
    fit <- suppressWarnings(glm(y ~ ., data = df, family = binomial(link = "probit"),
                                control = glm.control(maxit = 100)))
    summary(fit)$coefficients["gene", "Pr(>|z|)"]
  }, error = function(e) return(1))
}

# IMPORTANT: Exclude study_indicator from screening (matching Dec05 file)
gene_columns <- setdiff(colnames(X_cv), "study_indicator")
X_cv_genes_only <- X_cv[, gene_columns, drop = FALSE]

dev_pvals_screen <- apply(X_cv_genes_only, 2, function(x) probit_pvalue(Y_cv, x, Z_cv))
names(dev_pvals_screen) <- colnames(X_cv_genes_only)

# Compute coefficients for genes only (MATCHING ORIGINAL with explicit NA handling)
dev_coefs_screen <- sapply(colnames(X_cv_genes_only), function(gene) {
  df <- data.frame(y = Y_cv, gene = X_cv_genes_only[, gene], Z_cv)
  tryCatch({
    fit <- suppressWarnings(glm(y ~ ., data = df, family = binomial(link = "probit"),
                                control = glm.control(maxit = 100)))
    # Return coefficient even if non-converged or large
    # Separated variables often have large coefficients - these are real signals!
    gene_coef <- coef(fit)["gene"]

    # Only return 0 if coefficient is genuinely missing (NA)
    if (is.na(gene_coef)) return(0)

    gene_coef
  }, error = function(e) return(0))
})
names(dev_coefs_screen) <- colnames(X_cv_genes_only)

# Filter 1: p < 0.05
significant_genes_spbart <- names(dev_pvals_screen)[dev_pvals_screen < 0.05]
message("  Genes with p < 0.05: ", length(significant_genes_spbart))

# Stage 2: GMM clustering on |coefficients| (MATCHING ORIGINAL nested structure)
if (length(significant_genes_spbart) > 0) {
  dev_coefs_sig <- dev_coefs_screen[significant_genes_spbart]
  dev_coefs_sig <- dev_coefs_sig[!is.na(dev_coefs_sig)]

  if (length(dev_coefs_sig) >= 10) {
    abs_coefs_spbart <- abs(dev_coefs_sig)

    gmm_screen <- Mclust(abs_coefs_spbart, G = 2:4, modelNames = "V", verbose = FALSE)

    if (!is.null(gmm_screen)) {
      cluster_assignments_screen <- gmm_screen$classification
      cluster_means_screen <- sapply(1:gmm_screen$G, function(k) mean(abs_coefs_spbart[cluster_assignments_screen == k]))
      strongest_cluster_screen <- which.max(cluster_means_screen)
      screened_genes_spbart <- names(dev_coefs_sig)[cluster_assignments_screen == strongest_cluster_screen]
    } else {
      # GMM failed, use median threshold
      median_coef_screen <- median(abs_coefs_spbart)
      screened_genes_spbart <- names(dev_coefs_sig)[abs(dev_coefs_sig) > median_coef_screen]
    }
  } else {
    # Too few genes for GMM, use median threshold
    median_coef_screen <- median(abs(dev_coefs_sig))
    screened_genes_spbart <- names(dev_coefs_sig)[abs(dev_coefs_sig) > median_coef_screen]
  }
} else {
  screened_genes_spbart <- character(0)
}

message("  Screened genes for spBART: ", length(screened_genes_spbart))

# Filter X matrices to screened genes
X_cv_screened <- X_cv[, screened_genes_spbart, drop = FALSE]
X_test_screened <- X_test[, screened_genes_spbart, drop = FALSE]

# Add study_indicator back (matching Dec05 file)
X_cv_screened <- cbind(X_cv_screened, study_indicator = X_cv[, "study_indicator"])
X_test_screened <- cbind(X_test_screened, study_indicator = X_test[, "study_indicator"])

message("  Final X_cv for spBART: ", ncol(X_cv_screened), " variables (genes + study_indicator)")

# MCMC settings (matches writeup: 2000 burn-in, 5000 post-burn-in, thin by 5 -> Q=1000)
n_burn <- 2000
n_iter <- 1000
n_thin <- 5

# DART parameters - MUST MATCH ORIGINAL Dec05 FILE EXACTLY
# (See original file lines 1488-1500)
n_trees <- 200
sparse <- TRUE     # Use sparse (DART) prior
theta <- 0         # theta parameter for sparse prior (0 = automatic selection)
omega <- 1         # omega parameter
a <- 0.5           # Hyper-parameters for the beta distribution
b <- 1             # Hyper-parameters for the beta distribution
numcut <- 100      # Number of cut points for continuous variables
usequants <- FALSE # Use uniform cutpoints (not quantiles)

# Prior for clinical coefficients β
sigma_beta <- sqrt(10)  # β ~ N(0, 10*I)

rho_cv <- ncol(X_cv_screened)  # rho = number of screened genes + study_indicator

message(sprintf("DART: %d trees, sparse=%s", n_trees, sparse))
message(sprintf("Sparsity parameters: theta=%d, omega=%d, a=%.1f, b=%d", theta, omega, a, b))
message(sprintf("Prior: β ~ N(0, %.1f*I)", sigma_beta^2))

message("\n=== spBART - 5-Fold Cross-Validation (EXACTLY MATCHING Dec05 methodology) ===")

# ------------------------------------------------------------------------------
# Step 1: Create 5-Fold CV Partitions (EXACTLY MATCHING ORIGINAL)
# ------------------------------------------------------------------------------
K <- 5
# Seed matching Dec05 file: 890
set.seed(890)

fold_assignment <- createFolds(Y_cv, k = K, list = TRUE, returnTrain = FALSE)

message("Created ", K, "-fold stratified CV partitions:")
for (k in 1:K) {
  fold_k_indices <- fold_assignment[[k]]
  n_cases_k <- sum(Y_cv[fold_k_indices] == 1)
  n_controls_k <- sum(Y_cv[fold_k_indices] == 0)
  message(sprintf("Fold %d: %d patients (%d cases, %d controls)",
                  k, length(fold_k_indices), n_cases_k, n_controls_k))
}

# Define alpha levels for multi-threshold comparison
alpha_levels <- seq(0.05, 0.15, by = 0.025)

# Storage for fold-specific results
cv_PIPs_by_fold <- matrix(0, nrow = K, ncol = ncol(X_cv_screened))
colnames(cv_PIPs_by_fold) <- colnames(X_cv_screened)

# Storage for selected genes at ALL alpha levels for each fold
fold_selected_by_alpha <- vector("list", K)
for (k in 1:K) {
  fold_selected_by_alpha[[k]] <- list()
}

# ------------------------------------------------------------------------------
# Step 2: Run 5-Fold CV Loop
# ------------------------------------------------------------------------------
message("\nRunning 5-Fold CV for spBART variable selection...")
message("  Computing gene selection at ALL alpha levels: ", paste(alpha_levels, collapse = ", "))

# Loop over folds (EXACTLY MATCHING ORIGINAL lines 1530-1710)
for (k in 1:K) {

  message(sprintf("\n--- Fold %d/%d ---", k, K))

  # Step: Split data into training and validation (MATCHING ORIGINAL)
  val_k_indices <- fold_assignment[[k]]
  train_k_indices <- setdiff(1:length(Y_cv), val_k_indices)

  # Training data for fold k
  X_train_k <- X_cv_screened[train_k_indices, ]
  Z_train_k <- Z_cv[train_k_indices, ]
  Y_train_k <- Y_cv[train_k_indices]

  # Validation data for fold k
  X_val_k <- X_cv_screened[val_k_indices, ]
  Z_val_k <- Z_cv[val_k_indices, ]
  Y_val_k <- Y_cv[val_k_indices]

  message(sprintf("  Training: %d patients", length(Y_train_k)))
  message(sprintf("  Validation: %d patients", length(Y_val_k)))

  # Fit semi-parametric probit DART model (MATCHING ORIGINAL)
  message("  Fitting TRUE semi-parametric probit DART model...")
  message("  Using custom Gibbs sampler with explicit f(X) + Z'β decomposition")

  # Fit semi-parametric model using custom Gibbs sampler
  bart_fit_k <- fit_semiparametric_probit_DART(
    X = X_train_k,             # Gene expression only (high-dimensional)
    Z = Z_train_k,             # Clinical covariates only (low-dimensional)
    Y = Y_train_k,
    X_test = X_val_k,
    Z_test = Z_val_k,
    n_burn = n_burn,
    n_iter = n_iter,
    n_thin = n_thin,
    n_trees = n_trees,
    sparse = sparse,
    theta = theta,
    omega = omega,
    a = a,
    b = b,
    rho = rho_cv,              # Use rho = screened genes + 1 (study_indicator in X_cv_screened)
    numcut = numcut,
    usequants = usequants,
    sigma_beta = sigma_beta    # Prior std for clinical coefficients
  )

  message("  Model fitting complete.")

  # Compute fold-specific PIPs (MATCHING ORIGINAL lines 1622-1707)
  if (!is.null(bart_fit_k$varcount)) {
    # Compute PIP for each gene: proportion of draws where gene was used in DART
    var_used <- bart_fit_k$varcount > 0  # n_iter × D binary matrix
    PIP_k <- colMeans(var_used)  # Average over MCMC draws

    # Ensure names are set correctly
    names(PIP_k) <- colnames(X_cv_screened)

    # Store fold-specific PIPs (genes only)
    cv_PIPs_by_fold[k, ] <- PIP_k

    message(sprintf("  Top 5 genes by PIP in fold %d:", k))
    top5_idx <- order(PIP_k, decreasing = TRUE)[1:5]
    for (i in 1:5) {
      gene_idx <- top5_idx[i]
      message(sprintf("    %s: PIP = %.3f", colnames(X_cv_screened)[gene_idx], PIP_k[gene_idx]))
    }

    # Fold-specific FDR control for variable selection (MATCHING ORIGINAL)
    message(sprintf("  Applying fold-specific Bayesian FDR control..."))

    # Rank genes by decreasing PIP for this fold
    PIP_k_sorted <- sort(PIP_k, decreasing = TRUE)
    D_k <- length(PIP_k_sorted)

    # Compute fold-specific FDR for each R
    FDR_k <- numeric(D_k)
    for (r in 1:D_k) {
      FDR_k[r] <- sum(1 - PIP_k_sorted[1:r]) / r
    }

    # Diagnostic: Show FDR for top genes
    if (D_k >= 5) {
      message(sprintf("  FDR diagnostics for top genes:"))
      for (r in c(1, 5, 10, 20, 50)) {
        if (r <= D_k) {
          message(sprintf("    Top %d genes: FDR = %.4f, min PIP = %.4f",
                          r, FDR_k[r], PIP_k_sorted[r]))
        }
      }
    }

    # Apply FDR control at ALL alpha levels and store selected genes
    for (alpha in alpha_levels) {
      alpha_str <- as.character(alpha)
      valid_R_k <- which(FDR_k <= alpha)

      if (length(valid_R_k) > 0) {
        R_star_k <- max(valid_R_k)
        fold_selected_by_alpha[[k]][[alpha_str]] <- names(PIP_k_sorted)[1:R_star_k]
      } else {
        fold_selected_by_alpha[[k]][[alpha_str]] <- character(0)
      }
    }

    # Report summary for this fold
    n_at_005 <- length(fold_selected_by_alpha[[k]][["0.05"]])
    n_at_015 <- length(fold_selected_by_alpha[[k]][["0.15"]])
    message(sprintf("  Fold %d: genes selected at alpha=0.05: %d, alpha=0.15: %d",
                    k, n_at_005, n_at_015))
  } else {
    # No variable counts - set empty selections for all alpha levels
    for (alpha in alpha_levels) {
      fold_selected_by_alpha[[k]][[as.character(alpha)]] <- character(0)
    }
    message(sprintf("  WARNING: No variable counts available for fold %d", k))
  }

  message(sprintf("  Fold %d complete.", k))
}

message("\n=== 5-Fold Cross-Validation Complete ===")

# ------------------------------------------------------------------------------
# Step 3: Compute Union of Selected Genes at Each Alpha Level
# ------------------------------------------------------------------------------
message("\n=== Computing Union of Selected Genes at Each Alpha Level ===")

# For each alpha, compute union across folds
union_genes_by_alpha <- list()

for (alpha in alpha_levels) {
  alpha_str <- as.character(alpha)

  # Collect selected genes from all folds at this alpha
  union_at_alpha <- character(0)
  for (k in 1:K) {
    union_at_alpha <- union(union_at_alpha, fold_selected_by_alpha[[k]][[alpha_str]])
  }

  union_genes_by_alpha[[alpha_str]] <- union_at_alpha
}

# IMPORTANT: Ensure study_indicator is included in union_genes for each alpha level
# (Matching original Dec05 file behavior: study_indicator must be in final model)
message("\n=== Ensuring study_indicator is in union genes (matching Dec05 file) ===")
for (alpha in alpha_levels) {
  alpha_str <- as.character(alpha)
  if (!("study_indicator" %in% union_genes_by_alpha[[alpha_str]])) {
    message(sprintf("  alpha = %.3f: study_indicator NOT in union - adding it", alpha))
    union_genes_by_alpha[[alpha_str]] <- c(union_genes_by_alpha[[alpha_str]], "study_indicator")
  } else {
    message(sprintf("  alpha = %.3f: study_indicator already in union", alpha))
  }
}

# Report union sizes at different alpha levels
message("\nUnion sizes at different alpha levels (after adding study_indicator):")
for (alpha in alpha_levels) {
  alpha_str <- as.character(alpha)
  message(sprintf("  alpha = %.3f: %d genes in union", alpha, length(union_genes_by_alpha[[alpha_str]])))
}

# ------------------------------------------------------------------------------
# Step 4: Fit Final Prognostic Model for Each Alpha Level
# ------------------------------------------------------------------------------
message("\n=== Fitting Final Prognostic Models at Each Alpha Level ===")
message("  (One model per alpha level, using union genes for that alpha)")
message("  Using FULL clinical covariate set for final models: age, sex, race, bmi")

# Storage for final model results at each alpha
final_model_results <- list()

for (alpha in alpha_levels) {
  alpha_str <- as.character(alpha)
  union_genes <- union_genes_by_alpha[[alpha_str]]

  message(sprintf("\n--- Alpha = %.3f: %d union genes ---", alpha, length(union_genes)))

  if (length(union_genes) > 0) {
    # Subset X matrices to union genes for this alpha
    X_cv_union <- X_cv_screened[, union_genes, drop = FALSE]
    X_test_union <- X_test_screened[, union_genes, drop = FALSE]

    rho_alpha <- ncol(X_cv_union)

    # Fit final model on union genes for this alpha
    # Using Z_cv_full (all 4 clinical covariates) for final model
    # MUST PASS ALL PARAMETERS MATCHING ORIGINAL Dec05 FILE
    bart_fit_alpha <- fit_semiparametric_probit_DART(
      X = X_cv_union,
      Z = Z_cv_full,
      Y = Y_cv,
      X_test = X_test_union,
      Z_test = Z_test_full,
      n_burn = n_burn,
      n_iter = n_iter,
      n_thin = n_thin,
      n_trees = n_trees,
      sparse = sparse,
      theta = theta,
      omega = omega,
      a = a,
      b = b,
      rho = rho_alpha,
      numcut = numcut,
      usequants = usequants,
      sigma_beta = sigma_beta
    )

    # Compute PIPs for this model
    var_used_alpha <- bart_fit_alpha$varcount > 0
    PIP_alpha <- colMeans(var_used_alpha)
    names(PIP_alpha) <- union_genes

    # Sort by PIP and compute FDR curve
    PIP_alpha_sorted <- sort(PIP_alpha, decreasing = TRUE)
    D_alpha <- length(PIP_alpha_sorted)

    FDR_alpha <- numeric(D_alpha)
    for (k in 1:D_alpha) {
      FDR_alpha[k] <- sum(1 - PIP_alpha_sorted[1:k]) / k
    }

    # Apply FDR control at THIS alpha level to get final selected genes
    valid_K <- which(FDR_alpha <= alpha)
    if (length(valid_K) > 0) {
      K_star <- max(valid_K)
      selected_genes <- names(PIP_alpha_sorted)[1:K_star]
    } else {
      selected_genes <- character(0)
    }

    message(sprintf("  Final model: %d genes selected at FDR alpha = %.3f", length(selected_genes), alpha))

    # Store results
    final_model_results[[alpha_str]] <- list(
      union_genes = union_genes,
      n_union = length(union_genes),
      PIPs = PIP_alpha_sorted,
      FDR = FDR_alpha,
      selected_genes = selected_genes,
      n_selected = length(selected_genes)
    )
  } else {
    message("  No union genes - skipping model fitting")
    final_model_results[[alpha_str]] <- list(
      union_genes = character(0),
      n_union = 0,
      PIPs = numeric(0),
      FDR = numeric(0),
      selected_genes = character(0),
      n_selected = 0
    )
  }
}

message("\n=== All final models fitted ===")


# ==============================================================================
# SECTION 5: Multi-Threshold Gene Selection Comparison
# ==============================================================================

cat("\n================================================================================\n")
cat("  SECTION 5: Multi-Threshold Gene Selection Comparison\n")
cat("================================================================================\n\n")

# Storage for results
comparison_results <- data.frame(
  alpha = numeric(),
  n_union = integer(),
  n_spbart = integer(),
  n_enet = integer(),
  n_overlap = integer(),
  overlap_pct = numeric(),
  jaccard = numeric(),
  stringsAsFactors = FALSE
)

spbart_genes_by_alpha <- list()

message("Comparing spBART and Elastic Net gene selection at each FDR alpha level...")
message("Target (Elastic Net): ", n_enet, " genes\n")

for (alpha in alpha_levels) {
  alpha_str <- as.character(alpha)

  # Get results from final model at this alpha
  result_alpha <- final_model_results[[alpha_str]]
  selected_genes <- result_alpha$selected_genes

  # IMPORTANT: Exclude study_indicator from gene counts for fair comparison
  # study_indicator is a covariate, not a gene, so should not be counted in "genes selected"
  selected_genes_only <- setdiff(selected_genes, "study_indicator")
  n_union_genes_only <- result_alpha$n_union - ifelse("study_indicator" %in% result_alpha$union_genes, 1, 0)
  n_spbart <- length(selected_genes_only)

  # Compute overlap with Elastic Net (using genes only, not study_indicator)
  overlap_genes <- intersect(selected_genes_only, enet_final_genes)
  n_overlap <- length(overlap_genes)

  # Overlap percentage (relative to spBART selection)
  overlap_pct <- ifelse(n_spbart > 0, 100 * n_overlap / n_spbart, 0)

  # Jaccard similarity (using genes only, not study_indicator)
  union_count <- length(union(selected_genes_only, enet_final_genes))
  jaccard <- ifelse(union_count > 0, n_overlap / union_count, 0)

  # Store results
  comparison_results <- rbind(comparison_results, data.frame(
    alpha = alpha,
    n_union = n_union_genes_only,
    n_spbart = n_spbart,
    n_enet = n_enet,
    n_overlap = n_overlap,
    overlap_pct = round(overlap_pct, 1),
    jaccard = round(jaccard, 3),
    stringsAsFactors = FALSE
  ))

  # Store genes only (excluding study_indicator) for downstream analysis
  spbart_genes_by_alpha[[alpha_str]] <- selected_genes_only

  message(sprintf("  alpha = %.3f: union = %d, spBART = %d genes, overlap = %d (%.1f%%)",
                  alpha, n_union_genes_only, n_spbart, n_overlap, overlap_pct))
}

# Find the matching alpha
matching_idx <- which(comparison_results$n_spbart >= n_enet)[1]
if (!is.na(matching_idx)) {
  matching_alpha <- comparison_results$alpha[matching_idx]
  message("\nMatching alpha level: ", matching_alpha)
} else {
  matching_alpha <- max(comparison_results$alpha)
  message("\nNo exact match found. Maximum alpha tested: ", matching_alpha)
}


# ==============================================================================
# SECTION 6: Visualization
# ==============================================================================

cat("\n================================================================================\n")
cat("  SECTION 6: Visualization\n")
cat("================================================================================\n\n")

# 6a. Main comparison plot
message("Creating gene_selection_comparison.pdf...")

pdf(result_path("gene_selection_comparison.pdf"), width = 12, height = 8)

par(mar = c(6, 5, 4, 5))  # Increased bottom margin for percentages

# Calculate y-axis range with padding below 0
y_max <- max(c(n_enet, comparison_results$n_spbart)) * 1.25
y_min <- -y_max * 0.08  # Padding below 0 for percentages

# Plot spBART gene count with explicit x-axis labels
plot(comparison_results$alpha, comparison_results$n_spbart,
     type = "b", pch = 19, col = "steelblue", lwd = 2, cex = 1.5,
     xlab = "FDR Alpha Level", ylab = "Number of Selected Genes",
     main = "Gene Selection: Elastic Net vs. spBART at Different FDR Thresholds\n(Binary BMI >= 25, Study Indicator in BART)",
     xlim = c(min(comparison_results$alpha) - 0.01, max(comparison_results$alpha) + 0.01),
     ylim = c(y_min, y_max),
     cex.lab = 1.3, cex.main = 1.2,
     xaxt = "n")  # Suppress default x-axis

# Add explicit x-axis labels
axis(1, at = comparison_results$alpha,
     labels = sprintf("%.3f", comparison_results$alpha),
     cex.axis = 0.8, las = 1)

# Add horizontal line at y=0
abline(h = 0, col = "black", lwd = 1)

# Elastic Net reference line
abline(h = n_enet, col = "firebrick", lwd = 3, lty = 2)

# Overlap bars
bar_width <- 0.008
for (i in 1:nrow(comparison_results)) {
  rect(comparison_results$alpha[i] - bar_width,
       0,
       comparison_results$alpha[i] + bar_width,
       comparison_results$n_overlap[i],
       col = adjustcolor("forestgreen", alpha.f = 0.5),
       border = "darkgreen")
}

# Add overlap percentage below y=0 line
for (i in 1:nrow(comparison_results)) {
  overlap_pct <- comparison_results$overlap_pct[i]
  text(comparison_results$alpha[i],
       y_min * 0.5,  # Position below y=0
       labels = sprintf("%.1f%%", overlap_pct),
       cex = 0.7, col = "darkgreen")
}

# Vertical line at matching alpha
if (!is.na(matching_idx)) {
  abline(v = matching_alpha, col = "purple", lwd = 2, lty = 3)
}

# Legend
legend("topleft",
       legend = c(
         sprintf("Elastic Net (n = %d)", n_enet),
         "spBART gene count",
         "Overlap with Elastic Net",
         ifelse(!is.na(matching_idx),
                sprintf("Match at Alpha = %.3f", matching_alpha),
                "No match in range")
       ),
       col = c("firebrick", "steelblue", "forestgreen", "purple"),
       lty = c(2, 1, NA, 3),
       lwd = c(3, 2, NA, 2),
       pch = c(NA, 19, 15, NA),
       pt.cex = c(NA, 1.5, 2, NA),
       bg = "white", cex = 1.0)

# Add grid
grid(col = "gray80", lty = "dotted")

# Add text annotations at key points (spBART gene count above points)
text(comparison_results$alpha, comparison_results$n_spbart + max(comparison_results$n_spbart) * 0.05,
     labels = comparison_results$n_spbart, cex = 0.8, col = "steelblue")

dev.off()
message("  Saved: gene_selection_comparison.pdf")


# 6b. Overlap statistics plot
message("Creating overlap_statistics.pdf...")

pdf(result_path("overlap_statistics.pdf"), width = 12, height = 10)

par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))

# Plot 1: Overlap count
barplot(comparison_results$n_overlap,
        names.arg = sprintf("%.3f", comparison_results$alpha),
        col = "forestgreen", border = "darkgreen",
        main = "Number of Overlapping Genes by FDR Alpha",
        xlab = "FDR Alpha Level", ylab = "Number of Overlapping Genes",
        cex.names = 0.8, cex.lab = 1.2, cex.main = 1.3)
grid(col = "gray80", lty = "dotted")

# Plot 2: Overlap percentage
barplot(comparison_results$overlap_pct,
        names.arg = sprintf("%.3f", comparison_results$alpha),
        col = "dodgerblue", border = "darkblue",
        main = "Overlap Percentage (of spBART selection)",
        xlab = "FDR Alpha Level", ylab = "Overlap %",
        cex.names = 0.8, cex.lab = 1.2, cex.main = 1.3,
        ylim = c(0, 100))
grid(col = "gray80", lty = "dotted")

# Plot 3: Jaccard similarity
barplot(comparison_results$jaccard,
        names.arg = sprintf("%.3f", comparison_results$alpha),
        col = "coral", border = "darkred",
        main = "Jaccard Similarity Index",
        xlab = "FDR Alpha Level", ylab = "Jaccard Index",
        cex.names = 0.8, cex.lab = 1.2, cex.main = 1.3,
        ylim = c(0, 1))
grid(col = "gray80", lty = "dotted")

# Plot 4: Gene count comparison (line plot)
plot(comparison_results$alpha, comparison_results$n_spbart,
     type = "b", pch = 19, col = "steelblue", lwd = 2,
     main = "Gene Count: spBART vs Elastic Net",
     xlab = "FDR Alpha Level", ylab = "Number of Genes",
     cex.lab = 1.2, cex.main = 1.3,
     ylim = c(0, max(c(n_enet, comparison_results$n_spbart)) * 1.1))
abline(h = n_enet, col = "firebrick", lwd = 2, lty = 2)
legend("bottomright",
       legend = c("spBART", sprintf("Elastic Net (n=%d)", n_enet)),
       col = c("steelblue", "firebrick"),
       lty = c(1, 2), pch = c(19, NA), lwd = 2)
grid(col = "gray80", lty = "dotted")

dev.off()
message("  Saved: overlap_statistics.pdf")


# 6c. Comparison table (PDF)
message("Creating comparison_table.pdf...")

pdf(result_path("comparison_table.pdf"), width = 10, height = 8)

par(mar = c(1, 1, 3, 1))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

title(main = "Comparison of Gene Selection: Elastic Net vs. spBART\n(Binary BMI >= 25, Study Indicator in BART)",
      cex.main = 1.5, font.main = 2)

# Table header
header_y <- 0.9
col_x <- c(0.08, 0.22, 0.36, 0.50, 0.64, 0.80)
text(col_x, header_y, c("FDR Alpha", "spBART", "Elastic Net", "Overlap", "Overlap %", "Jaccard"),
     font = 2, cex = 1.0)

# Draw header line
lines(c(0.02, 0.95), c(header_y - 0.03, header_y - 0.03), lwd = 2)

# Table rows
row_height <- 0.06
for (i in 1:nrow(comparison_results)) {
  row_y <- header_y - 0.05 - (i * row_height)

  # Highlight matching row
  if (!is.na(matching_idx) && i == matching_idx) {
    rect(0.02, row_y - 0.025, 0.95, row_y + 0.025,
         col = adjustcolor("yellow", alpha.f = 0.3), border = NA)
  }

  text(col_x[1], row_y, sprintf("%.3f", comparison_results$alpha[i]), cex = 0.9)
  text(col_x[2], row_y, comparison_results$n_spbart[i], cex = 0.9)
  text(col_x[3], row_y, comparison_results$n_enet[i], cex = 0.9)
  text(col_x[4], row_y, comparison_results$n_overlap[i], cex = 0.9)
  text(col_x[5], row_y, sprintf("%.1f%%", comparison_results$overlap_pct[i]), cex = 0.9)
  text(col_x[6], row_y, sprintf("%.3f", comparison_results$jaccard[i]), cex = 0.9)
}

# Footer
text(0.5, 0.05, sprintf("Elastic Net reference: %d genes | Matching alpha = %.3f",
                        n_enet, matching_alpha),
     cex = 1.0, font = 3)

dev.off()
message("  Saved: comparison_table.pdf")


# 6d. Comparison table (TeX)
message("Creating comparison_table.tex...")

tex_file <- result_path("comparison_table.tex")

cat("\\begin{table}[ht]\n", file = tex_file)
cat("\\centering\n", file = tex_file, append = TRUE)
cat("\\caption{Comparison of Gene Selection: Elastic Net vs. spBART at Different FDR Thresholds (Binary BMI $\\geq$ 25)}\n",
    file = tex_file, append = TRUE)
cat("\\label{tab:comparison}\n", file = tex_file, append = TRUE)
cat("\\begin{tabular}{c c c c c c}\n", file = tex_file, append = TRUE)
cat("\\hline\n", file = tex_file, append = TRUE)
cat("\\textbf{FDR $\\alpha$} & \\textbf{spBART} & \\textbf{Elastic Net} & \\textbf{Overlap} & \\textbf{Overlap \\%} & \\textbf{Jaccard} \\\\\n",
    file = tex_file, append = TRUE)
cat("\\hline\n", file = tex_file, append = TRUE)

for (i in 1:nrow(comparison_results)) {
  cat(sprintf("%.3f & %d & %d & %d & %.1f\\%% & %.3f \\\\\n",
              comparison_results$alpha[i],
              comparison_results$n_spbart[i],
              comparison_results$n_enet[i],
              comparison_results$n_overlap[i],
              comparison_results$overlap_pct[i],
              comparison_results$jaccard[i]),
      file = tex_file, append = TRUE)
}

cat("\\hline\n", file = tex_file, append = TRUE)
cat("\\end{tabular}\n", file = tex_file, append = TRUE)
cat(sprintf("\\caption*{Note: Elastic Net selected %d genes. spBART matches at $\\alpha$ = %.3f.}\n",
            n_enet, matching_alpha),
    file = tex_file, append = TRUE)
cat("\\end{table}\n", file = tex_file, append = TRUE)

message("  Saved: comparison_table.tex")


# 6e. Overlap statistics table (TeX)
message("Creating overlap_statistics.tex...")

tex_file2 <- result_path("overlap_statistics.tex")

cat("\\begin{table}[ht]\n", file = tex_file2)
cat("\\centering\n", file = tex_file2, append = TRUE)
cat("\\caption{Overlap Statistics Between Elastic Net and spBART Gene Selection (Binary BMI $\\geq$ 25)}\n",
    file = tex_file2, append = TRUE)
cat("\\label{tab:overlap}\n", file = tex_file2, append = TRUE)
cat("\\begin{tabular}{c c c c}\n", file = tex_file2, append = TRUE)
cat("\\hline\n", file = tex_file2, append = TRUE)
cat("\\textbf{FDR $\\alpha$} & \\textbf{Overlap Count} & \\textbf{Overlap \\%} & \\textbf{Jaccard Index} \\\\\n",
    file = tex_file2, append = TRUE)
cat("\\hline\n", file = tex_file2, append = TRUE)

for (i in 1:nrow(comparison_results)) {
  cat(sprintf("%.3f & %d & %.1f\\%% & %.3f \\\\\n",
              comparison_results$alpha[i],
              comparison_results$n_overlap[i],
              comparison_results$overlap_pct[i],
              comparison_results$jaccard[i]),
      file = tex_file2, append = TRUE)
}

cat("\\hline\n", file = tex_file2, append = TRUE)
cat("\\end{tabular}\n", file = tex_file2, append = TRUE)
cat("\\end{table}\n", file = tex_file2, append = TRUE)

message("  Saved: overlap_statistics.tex")


# ==============================================================================
# SECTION 7: Gene Rank Analysis - spBART Genes in Elastic Net |β| Ranking
# ==============================================================================

cat("\n================================================================================\n")
cat("  SECTION 7: Gene Rank Analysis\n")
cat("  Where Do spBART Genes Rank in Elastic Net |β| Ranking?\n")
cat("================================================================================\n\n")

# ------------------------------------------------------------------------------
# Step 7.1: Extract Gene Sets and Compute |β| for Elastic Net Selected Genes
# ------------------------------------------------------------------------------

message("=== Step 7.1: Extracting Gene Sets ===\n")

# S_enet: Already defined as enet_final_genes
S_enet <- enet_final_genes
message("S_enet (Elastic Net final genes): ", length(S_enet), " genes")

# S_spbart: Genes selected at FDR α=0.05 from final spBART model
S_spbart <- final_model_results[["0.05"]]$selected_genes
# Remove study_indicator if present (we only want genes)
S_spbart <- setdiff(S_spbart, "study_indicator")
message("S_spbart (spBART FDR α=0.05 genes): ", length(S_spbart), " genes")

# Get |β| ONLY for genes in S_enet (final Elastic Net selection)
# coef_final was computed during elastic net analysis
coef_final_vec <- coef(prognostic_model_cv, s = "lambda.1se")[-1, ]

# Extract only the gene coefficients (exclude clinical covariates and study_indicator)
gene_coef_names <- intersect(names(coef_final_vec), S_enet)
beta_magnitude <- abs(coef_final_vec[gene_coef_names])
names(beta_magnitude) <- gene_coef_names

message("Genes with |β| computed: ", length(beta_magnitude))

# ------------------------------------------------------------------------------
# Step 7.2: Rank Genes by |β|
# ------------------------------------------------------------------------------

message("\n=== Step 7.2: Ranking Genes by |β| ===\n")

# Rank by decreasing |β| (rank 1 = largest |β|)
beta_rank <- rank(-beta_magnitude, ties.method = "average")
names(beta_rank) <- names(beta_magnitude)

# Create sorted version for visualization
beta_sorted <- sort(beta_magnitude, decreasing = TRUE)
gene_order <- names(beta_sorted)

n_enet_genes <- length(S_enet)
message("Total Elastic Net genes being ranked: ", n_enet_genes)

# ------------------------------------------------------------------------------
# Step 7.3: Classify Elastic Net Genes by spBART Overlap
# ------------------------------------------------------------------------------

message("\n=== Step 7.3: Classifying Genes by Overlap ===\n")

# Within S_enet, which genes are also in S_spbart?
S_overlap <- intersect(S_enet, S_spbart)
S_enet_only <- setdiff(S_enet, S_spbart)

message("Overlap (in both methods): ", length(S_overlap), " genes")
message("Elastic Net only: ", length(S_enet_only), " genes")

# Also report spBART-only genes (for context, not plotted in ranking)
S_spbart_only <- setdiff(S_spbart, S_enet)
message("spBART only (not in Elastic Net): ", length(S_spbart_only), " genes")

# Create classification vector for S_enet genes
gene_group <- rep("Enet only", length(S_enet))
names(gene_group) <- S_enet
gene_group[S_overlap] <- "Both"

# ------------------------------------------------------------------------------
# Step 7.4: Summary Statistics
# ------------------------------------------------------------------------------

message("\n=== Step 7.4: Summary Statistics ===\n")

# Compute summary statistics for each group
if (length(S_enet_only) > 0 && length(S_overlap) > 0) {
  summary_stats <- data.frame(
    Group = c("Enet only", "Also in spBART"),
    N = c(length(S_enet_only), length(S_overlap)),
    Mean_Beta = c(mean(beta_magnitude[S_enet_only]), mean(beta_magnitude[S_overlap])),
    Median_Beta = c(median(beta_magnitude[S_enet_only]), median(beta_magnitude[S_overlap])),
    Min_Beta = c(min(beta_magnitude[S_enet_only]), min(beta_magnitude[S_overlap])),
    Max_Beta = c(max(beta_magnitude[S_enet_only]), max(beta_magnitude[S_overlap])),
    Mean_Rank = c(mean(beta_rank[S_enet_only]), mean(beta_rank[S_overlap])),
    Median_Rank = c(median(beta_rank[S_enet_only]), median(beta_rank[S_overlap])),
    Pct_in_Top_Quartile = c(
      100 * mean(beta_rank[S_enet_only] <= n_enet_genes / 4),
      100 * mean(beta_rank[S_overlap] <= n_enet_genes / 4)
    )
  )

  message("Summary Statistics:")
  print(summary_stats)
} else {
  message("WARNING: One or both groups are empty. Skipping summary statistics.")
  summary_stats <- NULL
}

# ------------------------------------------------------------------------------
# Step 7.5: Statistical Test - Wilcoxon Rank-Sum Test
# ------------------------------------------------------------------------------

message("\n=== Step 7.5: Statistical Test ===\n")

if (length(S_overlap) >= 2 && length(S_enet_only) >= 2) {
  # Wilcoxon test: Are overlap genes ranked higher than Enet-only genes?
  wilcox_result <- wilcox.test(
    beta_rank[S_overlap],
    beta_rank[S_enet_only],
    alternative = "less"  # H1: overlap genes have lower (better) ranks
  )

  message("Wilcoxon Rank-Sum Test:")
  message("  H0: Overlap genes have same ranks as Enet-only genes")
  message("  H1: Overlap genes have lower (better) ranks")
  message("  W statistic: ", wilcox_result$statistic)
  message("  p-value: ", sprintf("%.4f", wilcox_result$p.value))
  message("  Interpretation: ",
          ifelse(wilcox_result$p.value < 0.05,
                 "Overlap genes are significantly higher ranked (p < 0.05)",
                 "No significant difference in ranks (p >= 0.05)"))
} else {
  message("WARNING: Insufficient samples for Wilcoxon test.")
  wilcox_result <- NULL
}

# ------------------------------------------------------------------------------
# Step 7.6: Visualizations (Manuscript-Ready)
# ------------------------------------------------------------------------------

message("\n=== Step 7.6: Creating Manuscript-Ready Visualizations ===\n")

# Define consistent colors
col_overlap <- "steelblue"
col_enet_only <- "gray60"
col_median_overlap <- "dodgerblue4"
col_median_enet <- "gray30"

# ============================================================================
# Plot A: Ranked Bar Plot with Overlap Highlighting
# ============================================================================

message("Creating beta_ranking_barplot.pdf...")

pdf(result_path("beta_ranking_barplot.pdf"), width = 10, height = 6)

# Set margins: bottom, left, top, right
par(mar = c(4, 5, 4, 1), mgp = c(3, 0.8, 0))

# Assign colors based on group membership
bar_colors <- rep(col_enet_only, length(gene_order))
names(bar_colors) <- gene_order
bar_colors[intersect(gene_order, S_overlap)] <- col_overlap

# Create barplot
bp <- barplot(beta_sorted,
              col = bar_colors[gene_order],
              border = NA,
              xaxt = "n",
              ylab = expression("|" * beta * "| (Coefficient Magnitude)"),
              main = "",
              cex.lab = 1.2,
              cex.axis = 1.0,
              las = 1)

# Add title with proper line spacing
title(main = expression("Elastic Net Gene Ranking by |" * beta * "|"),
      cex.main = 1.4, line = 2.5)
title(main = expression("Highlighting Overlap with spBART (FDR " * alpha * " = 0.05)"),
      cex.main = 1.1, line = 1.2, font.main = 1)

# Add x-axis label
mtext("Genes (ranked by decreasing coefficient magnitude)",
      side = 1, line = 2.5, cex = 1.0)

# Add legend in top right with proper positioning
legend("topright",
       legend = c(
         sprintf("Also in spBART (n = %d)", length(S_overlap)),
         sprintf("Elastic Net only (n = %d)", length(S_enet_only))
       ),
       fill = c(col_overlap, col_enet_only),
       border = c(col_overlap, col_enet_only),
       bty = "n",
       cex = 1.0,
       inset = c(0.02, 0.02))

# Add rotated gene names on top of each bar
text(x = bp,
     y = beta_sorted + 0.02 * max(beta_sorted),  # Slightly above each bar
     labels = gene_order,
     srt = 90,       # Rotate 90 degrees (vertical)
     adj = c(0, 0.5), # Left-justify and vertically center
     cex = 0.6,      # Smaller font size for readability
     xpd = TRUE)     # Allow drawing outside plot region

dev.off()
message("  Saved: beta_ranking_barplot.pdf")

# ============================================================================
# Plot B: Rank Position Strip Plot
# ============================================================================

message("Creating rank_strip_plot.pdf...")

pdf(result_path("rank_strip_plot.pdf"), width = 10, height = 5)

par(mar = c(5, 10, 4, 2), mgp = c(3, 0.8, 0))

plot(NULL,
     xlim = c(0.5, n_enet_genes + 0.5),
     ylim = c(0.3, 2.7),
     xlab = "",
     ylab = "",
     yaxt = "n",
     xaxt = "n",
     main = "",
     bty = "n")

# Add title
title(main = expression("Rank Distribution: Where Do Overlapping Genes Fall in |" * beta * "| Ranking?"),
      cex.main = 1.3, line = 2)

# Add x-axis
axis(1, at = pretty(c(1, n_enet_genes)), cex.axis = 1.0)
mtext(expression("Elastic Net |" * beta * "| Rank (1 = highest)"),
      side = 1, line = 3, cex = 1.1)

# Add horizontal bands
rect(0, 0.6, n_enet_genes + 1, 1.4, col = adjustcolor(col_enet_only, 0.2), border = NA)
rect(0, 1.6, n_enet_genes + 1, 2.4, col = adjustcolor(col_overlap, 0.2), border = NA)

# Plot rank positions as vertical lines
if (length(S_enet_only) > 0) {
  segments(x0 = beta_rank[S_enet_only], y0 = 0.7, y1 = 1.3,
           col = col_enet_only, lwd = 2)
}
if (length(S_overlap) > 0) {
  segments(x0 = beta_rank[S_overlap], y0 = 1.7, y1 = 2.3,
           col = col_overlap, lwd = 2)
}

# Y-axis labels
axis(2, at = c(1, 2),
     labels = c("Elastic Net\nonly", "Also in\nspBART"),
     las = 1, tick = FALSE, cex.axis = 1.0, line = -1)

# Add median rank markers with vertical dashed lines
if (length(S_enet_only) > 0) {
  median_enet <- median(beta_rank[S_enet_only])
  abline(v = median_enet, col = col_median_enet, lty = 2, lwd = 2)
}
if (length(S_overlap) > 0) {
  median_overlap <- median(beta_rank[S_overlap])
  abline(v = median_overlap, col = col_median_overlap, lty = 2, lwd = 2)
}

# Add legend with median ranks
if (length(S_enet_only) > 0 && length(S_overlap) > 0) {
  legend("topright",
         legend = c(
           sprintf("Enet only median: %.1f (n = %d)", median(beta_rank[S_enet_only]), length(S_enet_only)),
           sprintf("Overlap median: %.1f (n = %d)", median(beta_rank[S_overlap]), length(S_overlap))
         ),
         col = c(col_median_enet, col_median_overlap),
         lty = 2, lwd = 2,
         bty = "n",
         cex = 0.95,
         inset = c(0.02, 0.05))
}

dev.off()
message("  Saved: rank_strip_plot.pdf")

# ============================================================================
# Plot C: Boxplot Comparing Rank Distributions
# ============================================================================

message("Creating rank_boxplot.pdf...")

pdf(result_path("rank_boxplot.pdf"), width = 6, height = 7)

par(mar = c(5, 5, 4, 2), mgp = c(3, 0.8, 0))

if (length(S_enet_only) > 0 && length(S_overlap) > 0) {

  boxplot(
    beta_rank[S_enet_only],
    beta_rank[S_overlap],
    names = c("Elastic Net\nonly", "Also in\nspBART"),
    col = c(adjustcolor(col_enet_only, 0.7), adjustcolor(col_overlap, 0.7)),
    border = c(col_median_enet, col_median_overlap),
    ylab = "",
    main = "",
    cex.axis = 1.0,
    cex.lab = 1.2,
    las = 1,
    outline = TRUE,
    boxwex = 0.6,
    lwd = 1.5
  )

  # Add y-axis label with proper expression
  mtext(expression("Elastic Net |" * beta * "| Rank (lower = stronger effect)"),
        side = 2, line = 3.5, cex = 1.1)

  # Add title
  title(main = "Rank Distribution by\nGene Selection Method",
        cex.main = 1.3, line = 1.5)

  # Add individual points with jitter
  set.seed(123)  # For reproducible jitter
  stripchart(
    list(beta_rank[S_enet_only], beta_rank[S_overlap]),
    vertical = TRUE,
    method = "jitter",
    jitter = 0.12,
    pch = 19,
    cex = 0.7,
    col = adjustcolor("black", 0.4),
    add = TRUE
  )

  # Add sample sizes below boxes
  mtext(sprintf("n = %d", length(S_enet_only)), side = 1, at = 1, line = 2.8, cex = 0.95)
  mtext(sprintf("n = %d", length(S_overlap)), side = 1, at = 2, line = 2.8, cex = 0.95)

  # Add p-value annotation if test was performed
  if (!is.null(wilcox_result)) {
    p_text <- ifelse(wilcox_result$p.value < 0.001,
                     "p < 0.001",
                     sprintf("p = %.3f", wilcox_result$p.value))

    # Add bracket and p-value
    y_max <- max(beta_rank) * 1.05
    segments(1, y_max, 2, y_max, lwd = 1.5)
    segments(1, y_max * 0.98, 1, y_max, lwd = 1.5)
    segments(2, y_max * 0.98, 2, y_max, lwd = 1.5)
    text(1.5, y_max * 1.08, p_text, cex = 0.95)
  }

} else {
  plot.new()
  text(0.5, 0.5, "Insufficient data for boxplot", cex = 1.2)
}

dev.off()
message("  Saved: rank_boxplot.pdf")

# ============================================================================
# Plot D: Cumulative Proportion in Top-K
# ============================================================================

message("Creating cumulative_overlap_curve.pdf...")

pdf(result_path("cumulative_overlap_curve.pdf"), width = 8, height = 7)

par(mar = c(5, 5, 4, 2), mgp = c(3, 0.8, 0))

if (length(S_overlap) > 0) {

  # For each rank position K, compute proportion of overlap genes in top K
  gene_order_ranked <- names(sort(beta_rank))  # Genes ordered by rank (1 to n_enet)

  cumulative_overlap <- cumsum(gene_order_ranked %in% S_overlap)
  proportion_overlap <- cumulative_overlap / length(S_overlap)

  # Create plot
  plot(1:n_enet_genes, proportion_overlap,
       type = "l",
       lwd = 3,
       col = col_overlap,
       xlab = "",
       ylab = "",
       main = "",
       cex.axis = 1.0,
       las = 1,
       xlim = c(1, n_enet_genes),
       ylim = c(0, 1.05),
       xaxs = "i",
       yaxs = "i")

  # Add axis labels
  mtext(expression("Top K Genes by Elastic Net |" * beta * "| Rank"),
        side = 1, line = 3, cex = 1.1)
  mtext("Cumulative Proportion of\nspBART Overlap Genes Recovered",
        side = 2, line = 3, cex = 1.1)

  # Add title
  title(main = "Recovery of spBART Genes\nin Elastic Net Ranking",
        cex.main = 1.3, line = 1.5)

  # Reference: diagonal line (expected if overlap genes were randomly distributed)
  abline(a = 0, b = 1/n_enet_genes, lty = 2, col = "gray50", lwd = 2)

  # Add horizontal reference lines at 50% and 100%
  abline(h = 0.5, lty = 3, col = "gray70", lwd = 1)
  abline(h = 1.0, lty = 3, col = "gray70", lwd = 1)

  # Find K where we capture 50% of overlap genes
  K_50pct <- which(proportion_overlap >= 0.5)[1]
  if (!is.na(K_50pct)) {
    points(K_50pct, 0.5, pch = 19, col = "firebrick", cex = 1.5)

    # Add annotation with proper positioning
    text_x <- K_50pct + n_enet_genes * 0.05
    if (text_x > n_enet_genes * 0.8) text_x <- K_50pct - n_enet_genes * 0.05
    text(text_x, 0.58,
         sprintf("50%% at K = %d", K_50pct),
         col = "firebrick", cex = 0.95, adj = 0)
  }

  # Find K where we capture 100% of overlap genes
  K_100pct <- which(proportion_overlap >= 1.0)[1]
  if (!is.na(K_100pct)) {
    points(K_100pct, 1.0, pch = 19, col = "darkgreen", cex = 1.5)

    # Add annotation
    text(K_100pct, 0.92,
         sprintf("100%% at K = %d", K_100pct),
         col = "darkgreen", cex = 0.95, adj = 0.5)
  }

  # Add legend
  legend("bottomright",
         legend = c("Observed recovery", "Random expectation"),
         col = c(col_overlap, "gray50"),
         lty = c(1, 2),
         lwd = c(3, 2),
         bty = "n",
         cex = 1.0,
         inset = c(0.02, 0.02))

  # Add grid
  grid(col = "gray90", lty = 1)

  # Redraw the line on top of grid
  lines(1:n_enet_genes, proportion_overlap, lwd = 3, col = col_overlap)

} else {
  plot.new()
  text(0.5, 0.5, "No overlap genes to plot", cex = 1.2)
}

dev.off()
message("  Saved: cumulative_overlap_curve.pdf")

# ============================================================================
# Create Summary Table (LaTeX)
# ============================================================================

message("Creating rank_analysis_summary.tex...")

tex_file_rank <- result_path("rank_analysis_summary.tex")

cat("\\begin{table}[ht]\n", file = tex_file_rank)
cat("\\centering\n", file = tex_file_rank, append = TRUE)
cat("\\caption{Elastic Net $|\\beta|$ Rank Distribution by Gene Selection Method}\n",
    file = tex_file_rank, append = TRUE)
cat("\\label{tab:rank_analysis}\n", file = tex_file_rank, append = TRUE)
cat("\\begin{tabular}{l c c}\n", file = tex_file_rank, append = TRUE)
cat("\\hline\n", file = tex_file_rank, append = TRUE)
cat("\\textbf{Statistic} & \\textbf{Elastic Net Only} & \\textbf{Also in spBART} \\\\\n",
    file = tex_file_rank, append = TRUE)
cat("\\hline\n", file = tex_file_rank, append = TRUE)

if (!is.null(summary_stats)) {
  cat(sprintf("Number of genes & %d & %d \\\\\n",
              summary_stats$N[1], summary_stats$N[2]),
      file = tex_file_rank, append = TRUE)
  cat(sprintf("Mean $|\\beta|$ & %.4f & %.4f \\\\\n",
              summary_stats$Mean_Beta[1], summary_stats$Mean_Beta[2]),
      file = tex_file_rank, append = TRUE)
  cat(sprintf("Median $|\\beta|$ & %.4f & %.4f \\\\\n",
              summary_stats$Median_Beta[1], summary_stats$Median_Beta[2]),
      file = tex_file_rank, append = TRUE)
  cat(sprintf("Mean rank & %.1f & %.1f \\\\\n",
              summary_stats$Mean_Rank[1], summary_stats$Mean_Rank[2]),
      file = tex_file_rank, append = TRUE)
  cat(sprintf("Median rank & %.1f & %.1f \\\\\n",
              summary_stats$Median_Rank[1], summary_stats$Median_Rank[2]),
      file = tex_file_rank, append = TRUE)
  cat(sprintf("\\%% in top quartile & %.1f\\%% & %.1f\\%% \\\\\n",
              summary_stats$Pct_in_Top_Quartile[1], summary_stats$Pct_in_Top_Quartile[2]),
      file = tex_file_rank, append = TRUE)
}

cat("\\hline\n", file = tex_file_rank, append = TRUE)

if (!is.null(wilcox_result)) {
  p_text <- ifelse(wilcox_result$p.value < 0.001, "$<$ 0.001", sprintf("%.3f", wilcox_result$p.value))
  cat(sprintf("\\multicolumn{3}{l}{Wilcoxon rank-sum test: $p$ = %s} \\\\\n", p_text),
      file = tex_file_rank, append = TRUE)
}

cat("\\hline\n", file = tex_file_rank, append = TRUE)
cat("\\end{tabular}\n", file = tex_file_rank, append = TRUE)
cat("\\caption*{Note: Lower rank indicates stronger effect (larger $|\\beta|$). ",
    file = tex_file_rank, append = TRUE)
cat("spBART selection at FDR $\\alpha$ = 0.05.}\n", file = tex_file_rank, append = TRUE)
cat("\\end{table}\n", file = tex_file_rank, append = TRUE)

message("  Saved: rank_analysis_summary.tex")

# ============================================================================
# Create Detailed Gene Table (CSV)
# ============================================================================

message("Creating gene_rank_table.csv...")

gene_table <- data.frame(
  Gene = gene_order,
  Rank = beta_rank[gene_order],
  Beta_Magnitude = round(beta_magnitude[gene_order], 6),
  In_spBART = gene_order %in% S_spbart,
  Group = ifelse(gene_order %in% S_overlap, "Both", "Enet_only"),
  stringsAsFactors = FALSE
)

# Sort by rank
gene_table <- gene_table[order(gene_table$Rank), ]

write.csv(gene_table, result_path("gene_rank_table.csv"), row.names = FALSE)
message("  Saved: gene_rank_table.csv")

# Store rank analysis results for saving
rank_analysis_results <- list(
  S_enet = S_enet,
  S_spbart = S_spbart,
  S_overlap = S_overlap,
  S_enet_only = S_enet_only,
  S_spbart_only = S_spbart_only,
  beta_magnitude = beta_magnitude,
  beta_rank = beta_rank,
  gene_order = gene_order,
  summary_stats = summary_stats,
  wilcox_result = wilcox_result,
  gene_table = gene_table
)

message("\n=== Gene Rank Analysis Complete ===\n")


# ==============================================================================
# SECTION 8: Save Results
# ==============================================================================

cat("\n================================================================================\n")
cat("  SECTION 8: Save Results\n")
cat("================================================================================\n\n")

# Save comprehensive results
saveRDS(list(
  # Elastic Net results
  enet_final_genes = enet_final_genes,
  n_enet = n_enet,
  enet_stability_genes = enet_selected_genes_stability,
  enet_screened_genes = screened_genes,

  # spBART results (new structure with results per alpha)
  union_genes_by_alpha = union_genes_by_alpha,
  final_model_results = final_model_results,
  fold_selected_by_alpha = fold_selected_by_alpha,
  cv_PIPs_by_fold = cv_PIPs_by_fold,
  spbart_genes_by_alpha = spbart_genes_by_alpha,

  # Comparison results
  comparison_table = comparison_results,
  matching_alpha = matching_alpha,

  # At matching alpha
  spbart_genes_at_match = spbart_genes_by_alpha[[as.character(matching_alpha)]],
  overlap_genes_at_match = intersect(
    spbart_genes_by_alpha[[as.character(matching_alpha)]],
    enet_final_genes
  ),

  # Rank analysis results (Section 7)
  rank_analysis = rank_analysis_results,

  # Metadata
  n_cv = length(cv_ids),
  n_test = length(test_ids),
  alpha_levels = alpha_levels,
  seeds_used = list(
    data_partitioning = 789,
    elastic_net_replications = "789 + rep (rep = 1:100)",
    elastic_net_final_model = 456,
    spbart_cv_folds = 890
  ),
  bmi_cutoff = ">=25 (binary)",
  study_indicator_location = "In BART X (penalized like genes)",
  date = Sys.Date()

), result_path("comparative_analysis_results.RDS"))

message("Results saved to: comparative_analysis_results.RDS")

# Print summary
cat("\n================================================================================\n")
cat("  SUMMARY\n")
cat("================================================================================\n\n")

cat("Configuration:\n")
cat("  - BMI cutoff: >= 25 (binary)\n")
cat("  - Study indicator: In BART X (penalized like genes)\n")
cat("  - Screening covariates: age, sex_binary\n")
cat("  - Final model covariates: age, sex_binary, race_binary, bmi_binary\n\n")

cat("Seeds used (matching Dec05 files):\n")
cat("  - Data partitioning: 789\n")
cat("  - Elastic Net replications: 789 + rep (rep = 1:100)\n")
cat("  - Elastic Net final model: 456\n")
cat("  - spBART CV folds: 890\n\n")

cat("Elastic Net:\n")
cat("  - Final selected genes: ", n_enet, "\n\n")

cat("spBART at matching alpha (", matching_alpha, "):\n")
if (!is.na(matching_idx)) {
  cat("  - Selected genes: ", comparison_results$n_spbart[matching_idx], "\n")
  cat("  - Overlap with Elastic Net: ", comparison_results$n_overlap[matching_idx], "\n")
  cat("  - Overlap percentage: ", comparison_results$overlap_pct[matching_idx], "%\n")
  cat("  - Jaccard similarity: ", comparison_results$jaccard[matching_idx], "\n\n")
} else {
  cat("  - No matching alpha found in tested range\n\n")
}

cat("Output files saved to: ", results_dir, "\n")
cat("  - comparative_analysis_results.RDS\n")
cat("  - gene_selection_comparison.pdf\n")
cat("  - overlap_statistics.pdf\n")
cat("  - comparison_table.pdf\n")
cat("  - comparison_table.tex\n")
cat("  - overlap_statistics.tex\n")
cat("  - beta_ranking_barplot.pdf\n")
cat("  - rank_strip_plot.pdf\n")
cat("  - rank_boxplot.pdf\n")
cat("  - cumulative_overlap_curve.pdf\n")
cat("  - rank_analysis_summary.tex\n")
cat("  - gene_rank_table.csv\n")

script_end_time <- Sys.time()
cat("\nScript completed at:", format(script_end_time), "\n")
cat("Total runtime:", round(difftime(script_end_time, script_start_time, units = "mins"), 1), "minutes\n")
