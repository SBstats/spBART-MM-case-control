################################################################################
# 00_setup_BART_package.R
#
# PURPOSE:
# One-time setup script to install the modified BART package (BART_2.9.9.tar.gz)
# from the local Simulation Study directory.
#
# This script should be run ONCE before running any of the simulation scripts
# (01-04). The modified BART package includes:
#   - All standard BART functions (pbart, wbart, etc.) for Regular BART
#   - Stateful sampler functions (wbart_create, wbart_update_response,
#     wbart_run_iteration, wbart_destroy) for spBART
#
# USAGE:
# Run this script once in your R session before running scripts 01-04 in
# parallel across multiple RStudio sessions.
#
################################################################################

cat("\n")
cat("================================================================================\n")
cat("  00_setup_BART_package.R: Install Modified BART Package\n")
cat("================================================================================\n")
cat("Script started at:", format(Sys.time()), "\n\n")

# ==============================================================================
# SECTION 1: Determine Script Directory
# ==============================================================================

# Get script directory
script_dir <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    getwd()
  }
})

if (is.null(script_dir) || script_dir == "") {
  script_dir <- getwd()
}

cat("Script directory:", script_dir, "\n")

# ==============================================================================
# SECTION 2: Install Modified BART Package
# ==============================================================================

cat("\n=== Installing Modified BART Package ===\n")

# Step 1: Remove existing BART package if it exists
if ("BART" %in% rownames(installed.packages())) {
  cat("  Removing existing BART package from library...\n")
  remove.packages("BART")
  cat("  Existing BART package removed successfully.\n")
} else {
  cat("  No existing BART package found.\n")
}

# Step 2: Locate modified BART package tarball
# Check in current directory and parent directory
bart_path <- file.path(script_dir, "BART_2.9.9.tar.gz")
if (!file.exists(bart_path)) {
  bart_path <- file.path(dirname(script_dir), "Simulation Study Initial Screening New DGP Dec 30", "BART_2.9.9.tar.gz")
}

if (!file.exists(bart_path)) {
  stop(paste0("ERROR: Modified BART package not found.\n",
              "  Please copy BART_2.9.9.tar.gz to: ", script_dir))
}

cat(sprintf("  Found modified BART package at: %s\n", bart_path))

# Step 3: Install modified BART package from local tarball
cat("  Installing modified BART package...\n")

install.packages(bart_path, repos = NULL, type = "source", INSTALL_opts = "--no-multiarch")

cat("  Modified BART package installed successfully!\n\n")

# ==============================================================================
# SECTION 3: Verify Installation
# ==============================================================================

cat("=== Verifying BART Package Installation ===\n")

library(BART)

# Check version
cat("  BART package version:", as.character(packageVersion("BART")), "\n")

# Verify stateful sampler functions (for spBART)
spbart_functions <- c("wbart_create", "wbart_update_response",
                      "wbart_run_iteration", "wbart_destroy")
missing_spbart <- spbart_functions[!sapply(spbart_functions, exists)]

if (length(missing_spbart) > 0) {
  stop(paste0("ERROR: Missing spBART functions: ", paste(missing_spbart, collapse = ", ")))
}
cat("  [OK] spBART functions available:", paste(spbart_functions, collapse = ", "), "\n")

# Verify standard BART functions (for Regular BART)
regular_bart_functions <- c("pbart", "wbart")
missing_regular <- regular_bart_functions[!sapply(regular_bart_functions, function(f) {
  exists(f) || exists(f, where = asNamespace("BART"), mode = "function")
})]

# Check if pbart exists
if (exists("pbart", mode = "function") ||
    tryCatch({get("pbart", envir = asNamespace("BART")); TRUE}, error = function(e) FALSE)) {
  cat("  [OK] Regular BART functions available: pbart, wbart\n")
} else {
  warning("  [WARNING] pbart function not found - Regular BART may not work")
}

# ==============================================================================
# SECTION 4: Summary
# ==============================================================================

cat("\n================================================================================\n")
cat("  BART Package Setup Complete!\n")
cat("================================================================================\n")
cat("\n")
cat("The modified BART package (v2.9.9) has been installed and verified.\n")
cat("This package supports:\n")
cat("  - Regular BART with sparsity (DART) via pbart()\n")
cat("  - Semi-parametric BART (spBART) via wbart_create/update/run/destroy()\n")
cat("\n")
cat("You can now run scripts 01-04 in parallel across 4 separate RStudio sessions.\n")
cat("\n")

# ==============================================================================
# SECTION 5: Test Data Generation
# ==============================================================================

cat("\n================================================================================\n")
cat("  Testing Data Generation (Tree-based DGP)\n")
cat("================================================================================\n\n")

# Source data generation functions
source(file.path(script_dir, "utils", "data_generation.R"))

# Test both true models
for (tm in 1:2) {
  cat(sprintf("=== True Model %d ===\n", tm))

  # Generate example data
  test_data <- generate_simulation_data(
    n = 500,
    p_genes = 10000,
    n_true_genes = 10,
    true_model = tm,
    seed = 12345
  )

  # Use summary function
  summarize_simulation_data(test_data)
}

# Test class balance across multiple seeds
cat("\n=== Class Balance Stability Check ===\n\n")

check_class_balance(n_seeds = 20, n = 500, true_model = 1, base_seed = 1000)
check_class_balance(n_seeds = 20, n = 500, true_model = 2, base_seed = 2000)

# Show tree structure documentation
describe_trees(true_model = 1)
describe_trees(true_model = 2)

cat("\n")
cat("Script completed at:", format(Sys.time()), "\n")
