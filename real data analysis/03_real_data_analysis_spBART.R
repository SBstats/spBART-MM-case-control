################################################################################
# Semi-Parametric Probit DART Analysis for High-Dimensional Gene Expression
################################################################################
# 
#
# Usage on HPC:
#   Rscript 03_real_data_analysis.R
#
# Requirements:
#   - R >= 4.0.0
#   - Packages: BART, dplyr, pROC, MASS, future.apply, matrixStats, data.table,
#               readxl, gtools, BiocManager, DESeq2, caret, truncnorm, glmnet, mclust
#   - Memory: ~16GB RAM recommended
#   - Cores: Will use all available cores via future::availableCores()
#   - Runtime: 15-20 mins depending on hardware
################################################################################

# ==============================================================================
# HPC SETUP: Environment Configuration
# ==============================================================================

# Set options for batch/non-interactive execution
options(
  warn = 1,                    # Print warnings as they occur
  error = function() {         # Enhanced error reporting for debugging
    cat("\n========== ERROR OCCURRED ==========\n")
    cat("Error in script execution at:", date(), "\n")
    traceback(2)
    quit(status = 1)
  }
)

# Suppress interactive prompts
options(
  menu.graphics = FALSE,
  device.ask.default = FALSE
)

# Start timing
script_start_time <- Sys.time()
cat("\n")
cat("================================================================================\n")
cat("  Semi-Parametric Probit DART Analysis\n")
cat("================================================================================\n")
cat("Script started at:", format(script_start_time), "\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("Working directory:", getwd(), "\n")

# Detect available resources
n_cores <- parallel::detectCores()
available_memory_gb <- round(as.numeric(system("free -g | awk '/^Mem:/{print $2}'",
                                               intern = TRUE, ignore.stderr = TRUE)), 1)
if (length(available_memory_gb) == 0 || is.na(available_memory_gb)) {
  available_memory_gb <- "Unknown"
}

cat("Available CPU cores:", n_cores, "\n")
cat("Available memory:", available_memory_gb, "GB\n")
cat("================================================================================\n\n")

# ==============================================================================
# SETUP: Create Results Directory at Script Location
# ==============================================================================
# This script will create ONE folder called "results_spDART" at the same
# location as this script file. All outputs will be saved to this single folder:
#
# Data files (.RDS):
#   - cv_performance_metrics.RDS
#   - cv_PIPs.RDS
#   - selected_genes_FDR_control.RDS
#   - PIP_full_refit.RDS
#   - beta_posterior_draws.RDS
#   - beta_summary.RDS
#   - test_set_performance.RDS
#
# Visualization files (.pdf):
#   - roc_curves.pdf
#   - performance_distributions.pdf
#   - pip_analysis.pdf
#   - beta_coefficients.pdf
#   - calibration.pdf
#   - prediction_distributions.pdf
#
# LaTeX tables (.tex):
#   - performance_summary.tex
#   - beta_coefficients.tex
#   - top_selected_genes.tex
# ==============================================================================

# Get the directory where this script is located
script_path <- tryCatch({
  # Method 1: Works when using Rscript (HPC batch mode)
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

# Create single results directory with date stamp (all outputs in one folder)
# Format: results_spDART_YYYY-MM-DD (e.g., results_spDART_2025-01-15)
todays_date <- format(Sys.Date(), "%Y_%m_%d")
results_folder_name <- paste0("results_spDART_", todays_date)
results_dir <- file.path(script_path, results_folder_name)

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
  cat("Created results directory:", results_dir, "\n")
} else {
  cat("Using existing results directory:", results_dir, "\n")
  cat("Note: Results from a previous run on", todays_date, "will be overwritten\n")
}

# Helper function to construct file paths (all files in same folder)
result_path <- function(filename) {
  file.path(results_dir, filename)
}

cat("\n=== All results (data + plots) will be saved to: ===\n")
cat("  ", results_dir, "\n\n")





################################################################################
# 1) Load Required Libraries
################################################################################
cat("Loading required R packages...\n")

# Setup personal R library if on HPC (for package installations)
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") {
  user_lib <- file.path(Sys.getenv("HOME"), "R", "library", paste0(R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1]))
}
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  cat("Created personal R library:", user_lib, "\n")
}
.libPaths(c(user_lib, .libPaths()))
cat("R library paths:", paste(.libPaths(), collapse = ", "), "\n\n")

# Function to safely load packages with error handling
safe_library <- function(pkg, repos = "https://cloud.r-project.org") {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("  Installing", pkg, "...\n")
    install.packages(pkg, repos = repos, dependencies = TRUE, quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
  cat("  [OK]", pkg, "\n")
}



################################################################################
# CRITICAL: Install Modified BART Package with Stateful Sampler
################################################################################
# We use a MODIFIED version of the BART package that implements stateful
# sampling for use in Gibbs samplers. This is ESSENTIAL for semi-parametric
# probit DART models.
#
# Steps:
# 1. Remove any existing BART package from CRAN/internet
# 2. Install our modified BART package from local directory
#
# This ensures the supercomputer uses the correct version!
################################################################################

cat("\n=== Installing Modified BART Package ===\n")

# Step 1: Remove existing BART package if it exists
if ("BART" %in% rownames(installed.packages())) {
  cat("  Removing existing BART package from library...\n")
  remove.packages("BART")
  cat("  Existing BART package removed successfully.\n")
}

# Step 2: Install our modified BART package from local tarball
# The package should be a tarball called "BART_2.9.9.tar.gz" in the same directory as this script
bart_path <- file.path(script_path, "BART_2.9.9.tar.gz")

if (!file.exists(bart_path)) {
  stop(paste0("ERROR: Modified BART package not found at: ", bart_path,
              "\n  Please ensure the BART_2.9.9.tar.gz file is in the same directory as this script."))
}

cat(sprintf("  Installing modified BART package from: %s\n", bart_path))

install.packages(bart_path, repos = NULL, type = "source", INSTALL_opts = "--no-multiarch")

cat("  Modified BART package installed successfully!\n")
cat("  Package includes stateful sampler functions:\n")
cat("    - wbart_create()\n")
cat("    - wbart_update_response()\n")
cat("    - wbart_run_iteration()\n")
cat("    - wbart_destroy()\n\n")

# Step 3: Load the package and verify stateful functions are available
library(BART)

# Verify that our stateful functions are exported
required_functions <- c("wbart_create", "wbart_update_response",
                        "wbart_run_iteration", "wbart_destroy")

missing_functions <- required_functions[!sapply(required_functions, exists)]

if (length(missing_functions) > 0) {
  stop(paste0("ERROR: Modified BART package is missing required functions: ",
              paste(missing_functions, collapse = ", "),
              "\nPlease check that the BART package modifications were successful."))
}

cat("  [VERIFIED] All stateful sampler functions are available.\n\n")

# Load CRAN packages
safe_library("BART")            # Modified BART with stateful sampler (already loaded, but safe_library will verify)
safe_library("dplyr")           # Data manipulation with pipe operator (%>%)
safe_library("pROC")            # ROC curve analysis
safe_library("MASS")            # mvrnorm for multivariate normal
safe_library("future.apply")    # Parallel processing
safe_library("matrixStats")     # Fast matrix operations
safe_library("data.table")      # Data manipulation
safe_library("gtools")          # Miscellaneous tools
safe_library("readxl")          # Read Excel files
safe_library("glmnet")          # Elastic net (for comparison)
safe_library("caret")           # Cross-validation utilities
safe_library("mclust")          # Gaussian Mixture Models for coefficient clustering

# Load Bioconductor packages
cat("  Installing BiocManager if needed...\n")
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
}
library(BiocManager)
cat("  [OK] BiocManager\n")

if (!require("DESeq2", quietly = TRUE)) {
  cat("  Installing DESeq2 from Bioconductor...\n")
  BiocManager::install("DESeq2", ask = FALSE, update = FALSE, quiet = TRUE)
}
library(DESeq2)
cat("  [OK] DESeq2\n")

# Load truncnorm for semi-parametric model (will be loaded later but check now)
safe_library("truncnorm")

cat("All packages loaded successfully!\n\n")


#install.packages("remotes")
#library(remotes)

# Install the package from GitHub:
#install_github("rsparapa/bnptools", subdir="BART3")
#library(BART3)

################################################################################################################

################################################################################
# 2) Data pre-processing
################################################################################






################################################################################
# 2a) Pre-process UCHICAGO cohort data
################################################################################



################################################################################
# Load cohort data files


# Load ERRC ID data (Excel file): contains mapping of ERRC ID and sequencing ID
errc_seq_id_cohort <- as.data.frame(read_excel("data/ERRC.ID.xlsx"))

# Load gene expression data (RDS file): gene data for each sequencing ID
genebody_data_cohort <- as.data.frame(readRDS("data/genebody_377.RDS"))

# Load patient metadata (RData file): clinical data for each ERRC ID
load("data/patient_metadata_full_n797.RData")  #this data is named "analysis_data" # nolint
patient_metadata_cohort <- as.data.frame(analysis_data)

#head(patient_metadata_cohort)


# Subset patient metadata based on ERRC.ID from errc_seq_ID_cohort
patient_metadata_cohort_subset <- patient_metadata_cohort %>%
  filter(errcid %in% errc_seq_id_cohort$ERRC.ID)


#cohort working dataframe
cohort_workdf = data.frame(errcid = patient_metadata_cohort_subset$errcid,
                           age = patient_metadata_cohort_subset$age_diag,
                           race = patient_metadata_cohort_subset$race_composite,    
                           sex = patient_metadata_cohort_subset$sex_emr,
                           bmi = patient_metadata_cohort_subset$bmi_dx_emr,
                           #max_bmi = patient_metadata_cohort_subset$bmi_adultmax_qx,
                           mm_types = patient_metadata_cohort_subset$dx_errc)
#treatment = patient_metadata_cohort_subset$dtq,
#iss_stage = patient_metadata_cohort_subset$iss_derived
#response_raw = patient_metadata_cohort_subset$response_abstracted)
cohort_workdf = na.omit(cohort_workdf) # remove any rows with NAs



##################Pre-process cohort clinical data in order#########################################

################### 1) code response_binary to binary CR (1) vs non-CR (0)##################
#response_binary = ifelse(cohort_workdf$response_raw== "complete response",1,0)
#response_binary[is.na(response_binary)] = 0. #removing NAs with 0

#cohort_workdf$response_binary = response_binary



# 1) restrict to MM subtype only
#cohort_workdf = cohort_workdf[cohort_workdf$mm_types=="Multiple Myeloma"| cohort_workdf$mm_types=="MGUS",]
cohort_workdf = cohort_workdf[cohort_workdf$mm_types=="Multiple Myeloma",]
cohort_workdf$mm_status = rep("CASE", nrow(cohort_workdf))     #each patient is a case in the cohort study once we restrict to MM only

# 2) restrict to White and AA  only
cohort_workdf = cohort_workdf[cohort_workdf$race=="White" | cohort_workdf$race=="Black/African-American",]
#  code race to binary AA (0) vs W (1)
race_binary = ifelse(cohort_workdf$race== "White",1,0)
cohort_workdf$race_binary = race_binary

# 3) code bmi to binary (>=25 vs. <25) using BMI >= 25 cutoff (matches writeup Table 8.4)
bmi_binary = ifelse(cohort_workdf$bmi>= 25,1,0)
cohort_workdf$bmi_binary = bmi_binary


# 4) code sex to binary F (0) vs M (1)
sex_binary = ifelse(cohort_workdf$sex== "M",1,0)
#sex_binary[is.na(sex_binary)] = 0. #removing NAs with 0
cohort_workdf$sex_binary = sex_binary


# 5) Define a new identifier for pooled analysis
pooled_id_chorot = c(1:nrow(cohort_workdf))
cohort_workdf$PooledData_ID = pooled_id_chorot


#################### 4) restrict to dtq treatment only, combine triplet+quad, make it binary doublet (0) vs triplet+quad (1) ###################
#cohort_workdf = cohort_workdf[!is.na(cohort_workdf$treatment), ] #removes NAs and only keeps dtq trt patients
#cohort_workdf$treatment_binary = ifelse(cohort_workdf$treatment== "doub",0,1)

# cohort_workdf is a N X num_covariates df







##################Pre-process cohort 5-hmC data #########################################

# get sequencing ID of patients from the clinical dataset subset

seq_ID_cohort_final = errc_seq_id_cohort$Sequencing.ID[errc_seq_id_cohort$ERRC.ID %in%  cohort_workdf$errcid] #vec of length N





message("After data-preprocessing, the sample size from the cohort study is: ", length(seq_ID_cohort_final), ".\n")

cohort_workdf$sequencing_ID = seq_ID_cohort_final              # adding sequencing ID as a column in cohort workdf




genebody_data_cohort_unfiltered_unnormalized = genebody_data_cohort[,
                                                                    colnames(genebody_data_cohort) %in%
                                                                      seq_ID_cohort_final] #  num_5mC_seq X N df

# Map column names from sequencing_ID to PooledData_ID
# Create a mapping from sequencing_ID to PooledData_ID
seq_to_pooled_map <- setNames(cohort_workdf$PooledData_ID, cohort_workdf$sequencing_ID)

# Get current column names (sequencing IDs)
current_colnames <- colnames(genebody_data_cohort_unfiltered_unnormalized)

# Map them to PooledData_IDs
new_colnames <- seq_to_pooled_map[current_colnames]

# Replace column names
colnames(genebody_data_cohort_unfiltered_unnormalized) <- new_colnames

# Order columns by ascending column name (PooledData_ID)
genebody_data_cohort_unfiltered_unnormalized <- genebody_data_cohort_unfiltered_unnormalized[, order(as.numeric(colnames(genebody_data_cohort_unfiltered_unnormalized)))]

message("Columns ordered by ascending PooledData_ID")
message("First few column names after ordering: ", paste(head(colnames(genebody_data_cohort_unfiltered_unnormalized)), collapse = ", "))
message("Last few column names after ordering: ", paste(tail(colnames(genebody_data_cohort_unfiltered_unnormalized)), collapse = ", "))




################################################################################################################








################################################################################
# 2b) Pre-process CANADA case-control study
################################################################################



################################################################################
# Load case-control study files

# Load gene expression data (RDS file): gene data for each sequencing ID (19100 x 734 df)
genebody_data_case_control <- as.data.frame(readRDS("data/Canada_case_control_genebody_count.RDS"))


# Load patient metadata (password-protected xlsx file)
#library(xlsx)

# Read password-protected Excel file using xlsx package
# patient_metadata_case_control <- read.xlsx(
#  file = "data/MM_Questionnaire_BC for BART projects.xlsx",
#  sheetIndex = 1,
#  password = "MMVan2019"
# )

# Save as RDS file (R's native format, no password needed for future use)
#saveRDS(patient_metadata_case_control, "data/MM_Questionnaire_BC_for_BART_projects.RDS")

# Load patient metadata from RDS file (no password needed) (861 x 41 df)
patient_metadata_case_control <- readRDS("data/MM_Questionnaire_BC_for_BART_projects.RDS")

message("Table of Cases and Control before pre-processing:")
table(patient_metadata_case_control$MM_STATUS)


#head(patient_metadata_case_control)

# Load key connecting Study.ID (used to identify clinical covariates) and Assigned.ID (used to identify genebody)
Canada_case_control_sample_key <- as.data.frame(read.csv("data/Canada_case_control_sample_key.csv"))
#(734 x 2 df)
Canada_case_control_key = data.frame(Study.ID = Canada_case_control_sample_key$Study.ID,
                                     Assigned.ID = Canada_case_control_sample_key$Assigned.ID)




# Subset patient metadata based on Assigned.ID from errc_seq_ID_cohort (734 x 41 df)
patient_metadata_case_control_subset <- patient_metadata_case_control %>%
  filter(StudyID %in% Canada_case_control_key$Study.ID)


##################Compute current BMI using current height in inches and current weight in Lb#################

# Check HtCurrFt (Height in Feet)
message("Table of HtCurrFt:")
print(table(patient_metadata_case_control_subset$HtCurrFt))
message("Number of NAs in HtCurrFt: ", sum(is.na(patient_metadata_case_control_subset$HtCurrFt)))

# Check HtCurrIn (Height in Inches)
message("\nTable of HtCurrIn:")
print(table(patient_metadata_case_control_subset$HtCurrIn))
message("Number of NAs in HtCurrIn: ", sum(is.na(patient_metadata_case_control_subset$HtCurrIn)))

# Check WtCurrLb (Weight in Pounds)
message("\nTable of WtCurrLb:")
print(table(patient_metadata_case_control_subset$WtCurrLb))
message("Number of NAs in WtCurrLb: ", sum(is.na(patient_metadata_case_control_subset$WtCurrLb)))

# Compute BMI_current
# Step 1: Handle NAs in HtCurrIn (set to 0 if NA)
HtCurrIn_clean <- ifelse(is.na(patient_metadata_case_control_subset$HtCurrIn),
                         0,
                         patient_metadata_case_control_subset$HtCurrIn)

# Step 2: Compute total height in inches
HtCurrentInches <- patient_metadata_case_control_subset$HtCurrFt * 12 + HtCurrIn_clean

# Step 3: Compute BMI using formula: BMI = (weight_lb / height_inches^2) * 703
# If HtCurrFt or WtCurrLb are NA, BMI_current will be NA
BMI_current <- ifelse(is.na(patient_metadata_case_control_subset$HtCurrFt) |
                        is.na(patient_metadata_case_control_subset$WtCurrLb),
                      NA,
                      (patient_metadata_case_control_subset$WtCurrLb / (HtCurrentInches^2)) * 703)


patient_metadata_case_control_subset$BMI_current <- BMI_current
message("\nBMI_current computed:")
message("Number of valid BMI values: ", sum(!is.na(BMI_current)))
message("Number of NA BMI values: ", sum(is.na(BMI_current)))




################################################################################################################


#case-control working dataframe (734 x  6 df)
caseControl_workdf = data.frame(StudyID = patient_metadata_case_control_subset$StudyID,
                                age = patient_metadata_case_control_subset$Age,
                                race = patient_metadata_case_control_subset$EthSelf,    
                                sex = patient_metadata_case_control_subset$Sex,
                                bmi = patient_metadata_case_control_subset$BMI_current,
                                mm_status = patient_metadata_case_control_subset$MM_STATUS)


caseControl_workdf = na.omit(caseControl_workdf) # remove any rows with NAs (645 x  6 df)


# 1) restrict to White and AA  only
caseControl_workdf = caseControl_workdf[caseControl_workdf$race==1 | caseControl_workdf$race==5,]    # 1 = White, 5 = AA
#  code race to binary AA (0) vs W (1)
race_binary_case_control = ifelse(caseControl_workdf$race== 5,0,caseControl_workdf$race)
caseControl_workdf$race_binary = race_binary_case_control

# 2) code bmi to binary (>=25 vs. <25) using BMI >= 25 cutoff (matches writeup Table 8.4)
bmi_binary_case_control = ifelse(caseControl_workdf$bmi>= 25,1,0)
caseControl_workdf$bmi_binary = bmi_binary_case_control


# 3) code sex to binary F (0) vs M (1)
sex_binary_case_control = ifelse(caseControl_workdf$sex== "M",1,0)
caseControl_workdf$sex_binary = sex_binary_case_control


# 4) Define a new identifier for pooled analysis
pooled_id_caseControl = c( (nrow(cohort_workdf)+1): (nrow(cohort_workdf)+nrow(caseControl_workdf)) )
caseControl_workdf$PooledData_ID = pooled_id_caseControl


message("The dimension of caseControl_workdf after pre-processing is: ", nrow(caseControl_workdf), "x", ncol(caseControl_workdf), ". \n")
################################################################################################################





##################Pre-process cohort 5-hmC data #########################################



# get sequencing ID of patients from the clinical dataset subset

seq_ID_caseControl_final = Canada_case_control_key$Assigned.ID[Canada_case_control_key$Study.ID
                                                               %in%  caseControl_workdf$StudyID] #vec of length N



message("After data-preprocessing, the sample size from the Canada case-control study is: ", length(seq_ID_caseControl_final), ".\n")

caseControl_workdf$sequencing_ID = seq_ID_caseControl_final              # adding sequencing ID as a column in case control workdf




genebody_data_caseControl_unfiltered_unnormalized = genebody_data_case_control[,
                                                                               colnames(genebody_data_case_control) %in%
                                                                                 seq_ID_caseControl_final] #  num_5mC_seq X N df

# Map column names from sequencing_ID to PooledData_ID
# Create a mapping from sequencing_ID to PooledData_ID
seq_to_pooled_map_caseControl <- setNames(caseControl_workdf$PooledData_ID, caseControl_workdf$sequencing_ID)

# Get current column names (sequencing IDs)
current_colnames_caseControl <- colnames(genebody_data_caseControl_unfiltered_unnormalized)

# Map them to PooledData_IDs
new_colnames_caseControl <- seq_to_pooled_map_caseControl[current_colnames_caseControl]

# Replace column names
colnames(genebody_data_caseControl_unfiltered_unnormalized) <- new_colnames_caseControl

# Order columns by ascending column name (PooledData_ID)
genebody_data_caseControl_unfiltered_unnormalized <- genebody_data_caseControl_unfiltered_unnormalized[,
                                                                                                       order(as.numeric(colnames(genebody_data_caseControl_unfiltered_unnormalized)))]

message("Columns ordered by ascending PooledData_ID")
message("First few column names after ordering: ", paste(head(colnames(genebody_data_caseControl_unfiltered_unnormalized)), collapse = ", "))
message("Last few column names after ordering: ", paste(tail(colnames(genebody_data_caseControl_unfiltered_unnormalized)), collapse = ", "))
















################################################################################################################
# 2c) Pool metadata and genebody data
################################################################################################################



####################Pool metadata############################################

# Select variables from cohort_workdf
cohort_selected <- cohort_workdf %>%
  dplyr::select(PooledData_ID, race_binary, bmi_binary, sex_binary, age, mm_status)

# study_indicator coding: 0 = UCMM (UChicago cohort), 1 = BC (British Columbia case-control)
cohort_selected$study_indicator = 0


# Select variables from caseControl_workdf
caseControl_selected <- caseControl_workdf %>%
  dplyr::select(PooledData_ID, race_binary, bmi_binary, sex_binary, age, mm_status)

# study_indicator coding: 0 = UCMM (UChicago cohort), 1 = BC (British Columbia case-control)
caseControl_selected$study_indicator = 1

# Combine the two dataframes and arrange by PooledData_ID
pooled_metadata <- rbind(cohort_selected, caseControl_selected) %>%
  arrange(PooledData_ID)

####################Define binary outcome############################################
pooled_metadata_outcome = ifelse(pooled_metadata$mm_status=="CASE", 1, 0)
pooled_metadata$outcome = pooled_metadata_outcome


message("Pooled metadata created:")
head(pooled_metadata)
message("  The number of cases in the pooled dataset is: ", nrow(pooled_metadata[pooled_metadata$mm_status=="CASE",]))
message("  The number of controls in the pooled dataset is: ", nrow(pooled_metadata[pooled_metadata$mm_status=="CONTROL",]))
message("  Cohort rows: ", nrow(cohort_selected))
message("  Case-Control rows: ", nrow(caseControl_selected))
message("  Total pooled rows: ", nrow(pooled_metadata))






#############Pool genebody data##############################################



# Combine the two dataframes and arrange by PooledData_ID
pooled_genebody_data_unfiltered_unnormalized <- cbind(genebody_data_cohort_unfiltered_unnormalized,
                                                      genebody_data_caseControl_unfiltered_unnormalized) 

message("Pooled genebody data created:")
message("  Total pooled rows: ", nrow(pooled_genebody_data_unfiltered_unnormalized))
message("  Total pooled columns: ", ncol(pooled_genebody_data_unfiltered_unnormalized))









################################################################################################################
# 2d) Gene expression data preprocessing: filtering and normalization
################################################################################################################


# Step 1: Filter genes with <10 counts in >5% of samples
# Calculate threshold: 5% of 263 samples = 13.15, so >5% means >14 samples
n_samples <- ncol(pooled_genebody_data_unfiltered_unnormalized)
threshold_pct <- 0.05
threshold_samples <- ceiling(n_samples * threshold_pct)  # 14 samples

cat("Total samples:", n_samples, "\n")
cat("Filtering genes with <10 counts in more than", threshold_samples, "samples\n")

# For each gene (row), count how many samples have <10 counts
low_count_per_gene <- rowSums(pooled_genebody_data_unfiltered_unnormalized < 10)

# Keep genes where <10 counts occur in <=5% of samples 
genes_to_keep <- low_count_per_gene <= threshold_samples

cat("Genes before filtering:", nrow(pooled_genebody_data_unfiltered_unnormalized), "\n")
cat("Genes after filtering:", sum(genes_to_keep), "\n")
cat("Genes removed:", sum(!genes_to_keep), "\n")

# Filter the data
pooled_genebody_data_filtered_unnormalized <- pooled_genebody_data_unfiltered_unnormalized[genes_to_keep, ]
################################################################################################################


# Step 2: DESeq2 Normalization
# DESeq2 requires: counts matrix (genes x samples), column data with sample info


# Create colData (sample metadata) - minimal required structure
# DESeq2 needs this to create the dataset, but for normalization we just need patient IDs
colData <- data.frame(
  PooledData_ID = colnames(pooled_genebody_data_filtered_unnormalized),
  row.names = colnames(pooled_genebody_data_filtered_unnormalized)
)

message("Creating DESeq2 dataset for normalization...")
message("  Genes: ", nrow(pooled_genebody_data_filtered_unnormalized))
message("  Samples: ", ncol(pooled_genebody_data_filtered_unnormalized))

# Create DESeq2 dataset
# Note: DESeq2 expects integer counts
dds <- DESeqDataSetFromMatrix(
  countData = round(pooled_genebody_data_filtered_unnormalized),  # ensure integer counts
  colData = colData,
  design = ~ 1  # design = ~1 means no design formula (intercept-only model for normalization)
)


# Variance stabilizing transformation for downstream analysis
vsd <- vst(dds, blind = FALSE)
vsd_counts <- assay(vsd)

# Save results
pooled_genebody_data_filtered_normalized <- as.data.frame(vsd_counts)

cat("\nNormalization complete!\n")
cat("The range of normalized counts is : [", range(pooled_genebody_data_filtered_normalized) ,"]\n")
cat("Normalized data dimensions:", nrow(pooled_genebody_data_filtered_normalized), "genes x",
    ncol(pooled_genebody_data_filtered_normalized), "samples\n")

# Summary of dataframe
cat("\nSummary of normalized counts:\n")
head(pooled_genebody_data_filtered_normalized[,1:10])




################################################################################################################
################################################################################################################









#########################################################################################################################
# 3) Fit proposed semi-parametric probit DART model with 5-fold cross validation and held-out test strategy
#########################################################################################################################




# ==============================================================================
# MODEL SPECIFICATION
# ==============================================================================
# Semi-parametric probit DART model:
#   Y_i ~ Bernoulli(Φ(U_i))
#   U_i = f(X_i) + Z_i'β + ε_i,  ε_i ~ N(0,1)
#
# Where:
#   - f(X_i) = DART tree ensemble capturing nonlinear gene effects
#   - Z_i'β  = Linear effects of clinical covariates (interpretable)
#   - Φ(·)   = Standard normal CDF (probit link)
#
# Priors:
#   - β ~ N(0, σ_β²I) with σ_β² = 10
#   - (π₁,...,πD) ~ Dirichlet(α/D,...,α/D) for variable splitting probabilities
#   - Standard BART priors on tree structure and leaf parameters
# ==============================================================================

# ------------------------------------------------------------------------------
# Step 3.1: Data Partitioning (Stratified Train/Test Split)
# ------------------------------------------------------------------------------
# Strategy: Hold out N_test patients for final external validation
#           Use remaining N_CV patients for 5-fold cross-validation



################################################################################
# True Semi-Parametric Gibbs Sampler Implementation
################################################################################
# This implements the TRUE semi-parametric model with explicit separation
# of linear (β) and nonparametric (f) components:
#
# Model: Y ~ Bernoulli(Φ(f(X) + Z'β + ε)), ε ~ N(0,1)
#   where f(X) = DART function on genes (X)
#         Z'β = linear effects for clinical covariates (Z)
#
# Three-step Gibbs sampler:
#   1. Update U | Y, f, β via Albert-Chib truncated normal
#   2. Update β | U, f via conjugate Gaussian posterior
#   3. Update f | U, β via DART on residuals (U - Z'β)
################################################################################

# Load required packages
if (!require("truncnorm")) install.packages("truncnorm")
library(truncnorm)

################################################################################
# Semi-Parametric Probit DART Gibbs Sampler 
################################################################################
# This implementation uses the MODIFIED BART package with stateful sampling
# that properly maintains tree state across Gibbs iterations.
#
# KEY IMPROVEMENT over standard wbart():
#  Create sampler ONCE, then call run_iteration() (trees evolve properly) 
#
# This ensures:
#   1. Trees maintain their structure across iterations
#   2. MCMC mixing is efficient
#   3. Linero (2018) Dirichlet sparsity is fully functional
#   4. Proper integration into Gibbs sampler for semi-parametric probit model
#
# Model:
#   Y_i ~ Bernoulli(Φ(f(X_i) + Z_i'β))
#   f(X) ~ BART with Dirichlet sparsity (Linero 2018)
#   β ~ N(0, σ²_β I)
#
# where:
#   - Y: binary outcome
#   - X: high-dimensional gene expression (p >> n)
#   - Z: clinical covariates (low-dimensional)
#   - f: nonparametric function via BART
#   - β: parametric coefficients
#
# Gibbs sampler:
#   1. U | Y, f, β  (Albert-Chib latent variable augmentation)
#   2. β | U, f     (Conjugate Gaussian)
#   3. f | U, β     (BART with stateful updates)
################################################################################

fit_semiparametric_probit_DART <- function(X, Z, Y,
                                           X_test = NULL, Z_test = NULL,
                                           n_burn = 2000, n_iter = 1000, n_thin = 5,
                                           n_trees = 200, sparse = TRUE,
                                           theta = 0, omega = 1,
                                           a = 0.5, b = 1, rho = NULL,
                                           numcut = 100, usequants = FALSE,
                                           sigma_beta = sqrt(10)) {
  
  message("\n================================================================================")
  message("Semi-Parametric Probit DART with Stateful BART Sampler")
  message("================================================================================")
  
  # Data dimensions
  N <- length(Y)
  D <- ncol(X)  # Number of genes
  J <- ncol(Z)  # Number of clinical covariates
  
  message(sprintf("Data dimensions:"))
  message(sprintf("  N = %d observations", N))
  message(sprintf("  D = %d gene predictors", D))
  message(sprintf("  J = %d clinical covariates", J))
  message(sprintf("  Test set: %d observations", ifelse(is.null(X_test), 0, nrow(X_test))))
  
  # Set default rho
  if (is.null(rho)) {
    rho <- D  # Standard choice: concentration = number of predictors
  }
  
  message(sprintf("\nMCMC configuration:"))
  message(sprintf("  Burn-in: %d iterations", n_burn))
  message(sprintf("  Sampling: %d iterations (thinned every %d)", n_iter, n_thin))
  message(sprintf("  Total MCMC: %d iterations", n_burn + n_iter * n_thin))
  message(sprintf("  Trees: %d", n_trees))
  message(sprintf("  Dirichlet sparsity: %s (rho = %.1f)", ifelse(sparse, "ENABLED", "DISABLED"), rho))
  
  # Initialize parameters
  beta <- rep(0, J)  # Clinical coefficients
  U <- rnorm(N)      # Latent variables
  f_train <- rep(0, N)  # DART function values
  
  # Storage for posterior draws
  beta_draws <- matrix(NA, nrow = n_iter, ncol = J)
  f_draws <- matrix(NA, nrow = n_iter, ncol = N)
  varcount <- matrix(0, nrow = n_iter, ncol = D)  # Track variable usage for PIPs
  
  if (!is.null(X_test)) {
    N_test <- nrow(X_test)
    f_test_draws <- matrix(NA, nrow = n_iter, ncol = N_test)
    prob_test_draws <- matrix(NA, nrow = n_iter, ncol = N_test)
  }
  
  # ============================================================================
  # CRITICAL: Create STATEFUL BART sampler (only once!)
  # ============================================================================
  message("\n>>> Creating stateful BART sampler (trees initialized once)...")
  
  # Initial adjusted response
  Y_adjusted_init <- U - Z %*% beta
  
  # Create the sampler using MODIFIED BART package
  bart_sampler <- wbart_create(
    x.train = X,
    y.train = Y_adjusted_init,
    x.test = X_test,
    sparse = sparse,
    theta = theta,
    omega = omega,
    a = a,
    b = b,
    rho = rho,
    augment = FALSE,
    ntree = n_trees,
    numcut = numcut,
    usequants = usequants
  )
  
  message(">>> Sampler created successfully! Trees will evolve across iterations.\n")
  
  # ============================================================================
  # Gibbs sampler iterations
  # ============================================================================
  message("Starting Gibbs sampler...")
  message("================================================================================\n")
  
  n_total <- n_burn + (n_iter * n_thin)  # Total iterations including thinning
  
  for (iter in 1:n_total) {
    
    if (iter %% 500 == 0) {
      message(sprintf("  Iteration %d/%d (%.1f%%)",
                      iter, n_total, 100 * iter / n_total))
    }
    
    # ========================================================================
    # Step 1: Update latent variables U | Y, f, β (Albert-Chib)
    # ========================================================================
    mu_U <- f_train + Z %*% beta
    
    for (i in 1:N) {
      if (Y[i] == 1) {
        # Y=1: sample from truncated normal on (0, ∞)
        U[i] <- rtruncnorm(1, a = 0, b = Inf, mean = mu_U[i], sd = 1)
      } else {
        # Y=0: sample from truncated normal on (-∞, 0]
        U[i] <- rtruncnorm(1, a = -Inf, b = 0, mean = mu_U[i], sd = 1)
      }
    }
    
    # ========================================================================
    # Step 2: Update β | U, f (Conjugate Gaussian)
    # ========================================================================
    # Residual after removing DART function
    residuals <- U - f_train
    
    # Posterior variance: V_β = (Z'Z + σ_β^(-2) I)^(-1)
    precision_prior <- diag(J) / (sigma_beta^2)
    precision_post <- t(Z) %*% Z + precision_prior
    V_beta <- solve(precision_post)
    
    # Posterior mean: V_β Z' (U - f(X))
    mean_beta <- V_beta %*% (t(Z) %*% residuals)
    
    # Sample β
    beta <- MASS::mvrnorm(n = 1, mu = mean_beta, Sigma = V_beta)
    
    # ========================================================================
    # Step 3: Update DART function f | U, β (STATEFUL UPDATE - KEY CHANGE!)
    # ========================================================================
    # Create adjusted response: U - Z'β
    Y_adjusted <- as.numeric(U - Z %*% beta)
    
    # Update response in sampler (DOES NOT reinitialize trees!)
    wbart_update_response(bart_sampler, Y_adjusted)
    
    # Run ONE iteration from current tree state
    bart_result <- wbart_run_iteration(bart_sampler)
    
    # Extract fitted values
    f_train <- bart_result$yhat.train
    if (!is.null(X_test)) {
      f_test <- bart_result$yhat.test
    }
    
    # ========================================================================
    # Store posterior draws (after burn-in, with thinning)
    # ========================================================================
    # Only save every n_thin-th iteration after burn-in
    if (iter > n_burn && (iter - n_burn) %% n_thin == 0) {
      post_idx <- (iter - n_burn) / n_thin
      beta_draws[post_idx, ] <- beta
      f_draws[post_idx, ] <- f_train
      
      # Track which variables were used in this iteration
      varcount[post_idx, ] <- bart_result$varcount
      
      if (!is.null(X_test)) {
        f_test_draws[post_idx, ] <- f_test
        
        # Compute predicted probabilities for test set
        mu_test <- f_test + Z_test %*% beta
        prob_test_draws[post_idx, ] <- pnorm(mu_test)
      }
    }
  }
  
  # ============================================================================
  # Clean up
  # ============================================================================
  message("\n>>> Cleaning up sampler...")
  wbart_destroy(bart_sampler)
  
  message("================================================================================")
  message("Gibbs sampler completed successfully!")
  message("================================================================================\n")
  
  # Return results
  result <- list(
    beta_draws = beta_draws,
    f_train_draws = f_draws,
    varcount = varcount,  # For computing PIPs
    n_burn = n_burn,
    n_iter = n_iter,
    n_thin = n_thin
  )
  
  if (!is.null(X_test)) {
    result$f_test_draws <- f_test_draws
    result$prob.test <- prob_test_draws  # Match pbart() output format
  }
  
  return(result)
}


################################################################################
################################################################################













message("\n=== Step 3.1: Data Partitioning ===")

# Set random seed for reproducibility
set.seed(789)

# Define split sizes
N_total <- nrow(pooled_metadata)
N_CV    <-  500 # Cross-validation training pool
N_test  <- N_total - N_CV  # External held-out test set



message(sprintf("Total patients: %d", N_total))
message(sprintf("  - CV training pool: %d patients", N_CV))
message(sprintf("  - Held-out test set: %d patients", N_test))

# Stratified sampling to preserve outcome proportions
outcome_vec <- pooled_metadata$outcome
outcome_1_ids <- pooled_metadata$PooledData_ID[outcome_vec == 1]
outcome_0_ids <- pooled_metadata$PooledData_ID[outcome_vec == 0]

# Calculate proportion of cases
prop_cases <- sum(outcome_vec == 1) / N_total
n_test_cases <- round(N_test * prop_cases)
n_test_controls <- N_test - n_test_cases

# Sample test set with balanced outcomes
test_ids_cases <- sample(outcome_1_ids, n_test_cases, replace = FALSE)
test_ids_controls <- sample(outcome_0_ids, n_test_controls, replace = FALSE)
test_ids <- c(test_ids_cases, test_ids_controls)

# Remaining IDs go to CV pool
cv_ids <- setdiff(pooled_metadata$PooledData_ID, test_ids)

message(sprintf("\nTest set: %d cases, %d controls (%.1f%% cases)",
                n_test_cases, n_test_controls, 100*prop_cases))
message(sprintf("CV pool: %d cases, %d controls",
                sum(pooled_metadata$outcome[pooled_metadata$PooledData_ID %in% cv_ids] == 1),
                sum(pooled_metadata$outcome[pooled_metadata$PooledData_ID %in% cv_ids] == 0)))

# Extract data matrices for CV pool and test set
cv_indices <- which(pooled_metadata$PooledData_ID %in% cv_ids)
test_indices <- which(pooled_metadata$PooledData_ID %in% test_ids)

# High-dimensional gene expression (transpose: rows=patients, cols=genes)
# PLUS study_indicator as an additional predictor in the BART model
X_genes_cv <- t(pooled_genebody_data_filtered_normalized[, cv_indices])
X_genes_test <- t(pooled_genebody_data_filtered_normalized[, test_indices])

# Add study_indicator as a column to X (goes into BART along with genes)
study_indicator_cv <- pooled_metadata$study_indicator[cv_indices]
study_indicator_test <- pooled_metadata$study_indicator[test_indices]

X_cv <- cbind(X_genes_cv, study_indicator = study_indicator_cv)
X_test <- cbind(X_genes_test, study_indicator = study_indicator_test)

# Low-dimensional clinical covariates (NO study_indicator - it's now in X for BART)
# NOTE: No manual intercept added - BART internally includes an intercept
# For initial screening and model development (5-fold CV): use only age and sex_binary
# For final prognostic model: add bmi_binary and race_binary back
clinical_vars_cv <- c("sex_binary", "age")  # Only age and sex for CV
clinical_vars_full <- c("race_binary", "bmi_binary", "sex_binary", "age")  # All covariates for final model

Z_cv <- as.matrix(pooled_metadata[cv_indices, clinical_vars_cv])
Z_test <- as.matrix(pooled_metadata[test_indices, clinical_vars_cv])

# Also prepare full Z matrices for final prognostic model (will be used later)
Z_cv_full <- as.matrix(pooled_metadata[cv_indices, clinical_vars_full])
Z_test_full <- as.matrix(pooled_metadata[test_indices, clinical_vars_full])

# Binary outcomes
Y_cv <- pooled_metadata$outcome[cv_indices]
Y_test <- pooled_metadata$outcome[test_indices]

message(sprintf("\nData matrices prepared:"))
message(sprintf("  X_cv: %d × %d (patients × [genes + study_indicator])", nrow(X_cv), ncol(X_cv)))
message(sprintf("    - %d genes + 1 study_indicator (included in BART)", ncol(X_genes_cv)))
message(sprintf("  Z_cv: %d × %d (patients × clinical covariates for CV: sex, age)", nrow(Z_cv), ncol(Z_cv)))
message(sprintf("  Z_cv_full: %d × %d (patients × all clinical covariates for final model: race, bmi, sex, age)", nrow(Z_cv_full), ncol(Z_cv_full)))
message(sprintf("  Y_cv: %d binary outcomes", length(Y_cv)))
message(sprintf("\nNOTE: Initial screening and 5-fold CV use Z_cv (sex, age only)"))
message(sprintf("      Final prognostic model will use Z_cv_full (race, bmi, sex, age)"))




# ------------------------------------------------------------------------------
# Step 3.1: Initial Gene Screening using Probit Regression with GMM
# ------------------------------------------------------------------------------
# Screening on development set (X_cv) only:
#   1. Fit probit regression: Y ~ gene + Z for each gene
#   2. Filter 1: Keep genes with p < 0.05
#   3. Filter 2: Fit GMM to |coefficients|, keep strongest cluster (highest mean |β|)
# All downstream DART analysis will use only the screened genes

message("\n=== Step 3.1: Initial Gene Screening (Probit Regression + GMM) ===")

# Helper function: probit regression p-value
# Suppress warnings but keep all coefficients (even non-converged)
# NOTE: Z does not include intercept column - let glm add intercept automatically
probit_pvalue <- function(y, x, Z) {
  df <- data.frame(y = y, gene = x, Z)
  tryCatch({
    # Suppress convergence warnings and limit iterations
    # Use y ~ . to include intercept in the probit screening model
    fit <- suppressWarnings(
      glm(y ~ ., data = df, family = binomial(link = "probit"),
          control = glm.control(maxit = 100))
    )
    
    # Return p-value even if non-converged
    # For separated variables, p-value will be very small (as it should be)
    summary(fit)$coefficients["gene", "Pr(>|z|)"]
  }, error = function(e) return(1))
}




# Stage 1: Compute p-values for all genes (EXCLUDING study_indicator)
message("\nStage 1: Computing p-values for all genes...")

# IMPORTANT: Exclude study_indicator from screening - it will be added back later
gene_columns <- setdiff(colnames(X_cv), "study_indicator")
X_cv_genes_only <- X_cv[, gene_columns, drop = FALSE]

message(sprintf("  Analyzing %d genes in development set (%d patients)", ncol(X_cv_genes_only), nrow(X_cv_genes_only)))
message(sprintf("  NOTE: study_indicator excluded from screening (will be added back after GMM clustering)"))

# Compute p-values for genes only (not study_indicator)
dev_pvals_screen <- apply(X_cv_genes_only, 2, function(x) probit_pvalue(Y_cv, x, Z_cv))
names(dev_pvals_screen) <- colnames(X_cv_genes_only)

# Compute coefficients for genes only (not study_indicator)
dev_coefs_screen <- sapply(colnames(X_cv_genes_only), function(gene) {
  df <- data.frame(y = Y_cv, gene = X_cv_genes_only[, gene], Z_cv)
  tryCatch({
    # Suppress convergence warnings and limit iterations
    # Use y ~ . to include intercept in the probit screening model
    fit <- suppressWarnings(
      glm(y ~ ., data = df, family = binomial(link = "probit"),
          control = glm.control(maxit = 100))
    )

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
significant_genes_dev_screen <- names(dev_pvals_screen)[dev_pvals_screen < 0.05]
message(sprintf("  Genes with p < 0.05: %d out of %d",
                length(significant_genes_dev_screen), ncol(X_cv)))

# Stage 2: GMM clustering on |coefficients|
message("\nStage 2: Gaussian Mixture Model clustering on |coefficients|...")

if (length(significant_genes_dev_screen) > 0) {
  dev_coefs_sig <- dev_coefs_screen[significant_genes_dev_screen]
  dev_coefs_sig <- dev_coefs_sig[!is.na(dev_coefs_sig)]
  
  if (length(dev_coefs_sig) >= 10) {
    message(sprintf("  Fitting GMM to %d significant gene coefficients...",
                    length(dev_coefs_sig)))
    
    # Use absolute values for clustering
    abs_coefs_dev_screen <- abs(dev_coefs_sig)
    
    # Fit GMM with 3-5 components (BIC selects optimal, minimum 3 clusters)
    # 3 clusters ensures: weak, moderate, strong effect separation
    gmm_screen <- Mclust(abs_coefs_dev_screen, G = 2:4, modelNames = "V", verbose = FALSE)
    
    if (!is.null(gmm_screen)) {
      # Get cluster assignments
      cluster_assignments_screen <- gmm_screen$classification
      n_clusters_screen <- gmm_screen$G
      
      message(sprintf("  Optimal number of clusters: %d", n_clusters_screen))
      
      # Compute mean |coefficient| for each cluster
      cluster_means_screen <- sapply(1:n_clusters_screen, function(k) {
        mean(abs_coefs_dev_screen[cluster_assignments_screen == k])
      })
      
      # Identify strongest cluster (highest mean |β|)
      strongest_cluster_screen <- which.max(cluster_means_screen)
      
      # Display cluster statistics
      message("  Cluster interpretation:")
      cluster_labels <- c("Weakly predictive", "Moderately predictive", "Strongly predictive")
      cluster_order <- order(cluster_means_screen)
      
      for (i in 1:n_clusters_screen) {
        k <- cluster_order[i]
        n_genes_k <- sum(cluster_assignments_screen == k)
        mean_k <- cluster_means_screen[k]
        sd_k <- sd(abs_coefs_dev_screen[cluster_assignments_screen == k])
        label <- if(i <= length(cluster_labels)) cluster_labels[i] else paste("Cluster", i)
        is_strongest <- ifelse(k == strongest_cluster_screen, " <- STRONGEST (SELECTED)", "")
        
        message(sprintf("    Cluster %d (%s): %d genes, mean |coef| = %.4f (SD = %.4f)%s",
                        k, label, n_genes_k, mean_k, sd_k, is_strongest))
      }
      
      # Keep genes in strongest cluster
      screened_genes_final <- names(dev_coefs_sig)[cluster_assignments_screen == strongest_cluster_screen]
      
      message(sprintf("\n  Genes in strongest cluster (Cluster %d): %d out of %d",
                      strongest_cluster_screen, length(screened_genes_final),
                      length(dev_coefs_sig)))
      
      # Save GMM results
      gmm_screening_results <- list(
        model = gmm_screen,
        n_clusters = n_clusters_screen,
        cluster_means = cluster_means_screen,
        strongest_cluster = strongest_cluster_screen,
        cluster_assignments = cluster_assignments_screen,
        abs_coefficients = abs_coefs_dev_screen
      )
      
    } else {
      # GMM failed, use median threshold
      message("  Warning: GMM fitting failed, using median |coefficient| threshold")
      median_coef_screen <- median(abs(dev_coefs_sig))
      screened_genes_final <- names(dev_coefs_sig)[abs(dev_coefs_sig) > median_coef_screen]
      gmm_screening_results <- NULL
      
      message(sprintf("  Median |coefficient|: %.4f", median_coef_screen))
      message(sprintf("  Genes with |coef| > median: %d out of %d",
                      length(screened_genes_final), length(dev_coefs_sig)))
    }
    
  } else {
    # Too few genes for GMM
    message(sprintf("  Too few significant genes (%d) for GMM, using median threshold",
                    length(dev_coefs_sig)))
    median_coef_screen <- median(abs(dev_coefs_sig))
    screened_genes_final <- names(dev_coefs_sig)[abs(dev_coefs_sig) > median_coef_screen]
    gmm_screening_results <- NULL
    
    message(sprintf("  Median |coefficient|: %.4f", median_coef_screen))
    message(sprintf("  Genes with |coef| > median: %d out of %d",
                    length(screened_genes_final), length(dev_coefs_sig)))
  }
} else {
  screened_genes_final <- character(0)
  gmm_screening_results <- NULL
}

# Check if we have genes to proceed
if (length(screened_genes_final) == 0) {
  stop("No genes passed screening. Analysis cannot proceed.")
}

# Display summary
message("\n=== Screening Summary ===")
if (!is.null(gmm_screening_results)) {
  message(sprintf("Development set: %d → %d genes (p<0.05 + GMM strongest cluster)",
                  length(significant_genes_dev_screen), length(screened_genes_final)))
  message(sprintf("  GMM: %d clusters identified, keeping Cluster %d (highest mean |coef|)",
                  gmm_screening_results$n_clusters, gmm_screening_results$strongest_cluster))
} else {
  message(sprintf("Development set: %d → %d genes (p<0.05 + |coef|>median)",
                  length(significant_genes_dev_screen), length(screened_genes_final)))
}

# Filter X_cv and X_test to only include screened genes
message("\nFiltering gene expression matrices to screened genes only...")
X_cv_original <- X_cv  # Keep original for reference
X_test_original <- X_test

X_cv <- X_cv[, screened_genes_final, drop = FALSE]
X_test <- X_test[, screened_genes_final, drop = FALSE]

message(sprintf("  Development set (X_cv): %d × %d → %d × %d (patients × genes)",
                nrow(X_cv_original), ncol(X_cv_original), nrow(X_cv), ncol(X_cv)))
message(sprintf("  Test set (X_test): %d × %d → %d × %d (patients × genes)",
                nrow(X_test_original), ncol(X_test_original), nrow(X_test), ncol(X_test)))

# IMPORTANT: Add study_indicator back to X_cv and X_test after gene screening
# Study indicator was excluded from GMM clustering but must be included in BART model
message("\n  Adding study_indicator back to X matrices (excluded from screening, included in BART)...")
X_cv <- cbind(X_cv, study_indicator = X_cv_original[, "study_indicator"])
X_test <- cbind(X_test, study_indicator = X_test_original[, "study_indicator"])

message(sprintf("  Final X_cv: %d × %d (patients × [screened genes + study_indicator])",
                nrow(X_cv), ncol(X_cv)))
message(sprintf("  Final X_test: %d × %d (patients × [screened genes + study_indicator])",
                nrow(X_test), ncol(X_test)))

# Save screening results
if (!dir.exists("dart_results")) {
  dir.create("dart_results", recursive = TRUE)
}

saveRDS(list(
  # Development set results
  dev_significant_genes = significant_genes_dev_screen,
  dev_pvalues = dev_pvals_screen[significant_genes_dev_screen],
  dev_coefficients = dev_coefs_screen[significant_genes_dev_screen],
  
  # GMM results (if used)
  gmm_results = gmm_screening_results,
  
  # Final screened genes
  screened_genes = screened_genes_final,
  screened_genes_coefficients = dev_coefs_screen[screened_genes_final],
  screened_genes_pvalues = dev_pvals_screen[screened_genes_final],
  
  # Data info
  n_dev_patients = nrow(X_cv_original),
  n_test_patients = nrow(X_test_original),
  original_n_genes = ncol(X_cv_original),
  final_n_genes = length(screened_genes_final)
), "dart_results/probit_gmm_screening_results.RDS")

message("\nScreening results saved to: dart_results/probit_gmm_screening_results.RDS")

# Set rho for CV based on number of screened genes entering BART
rho_cv <- ncol(X_cv)
message(sprintf("\nSetting rho for CV: rho = %d (screened genes + 1 for study_indicator in X_cv)", rho_cv))





















# ------------------------------------------------------------------------------
# Step 3.2: 5-Fold Cross-Validation Setup
# ------------------------------------------------------------------------------
# Create stratified K-fold partitions on CV pool

message("\n=== Step 3.2: Creating 5-Fold CV Partitions ===")

K <- 5  # Number of folds
set.seed(890)

# Create stratified folds
library(caret)
fold_assignment <- createFolds(Y_cv, k = K, list = TRUE, returnTrain = FALSE)

# Verify fold balance
for (k in 1:K) {
  fold_k_indices <- fold_assignment[[k]]
  n_cases_k <- sum(Y_cv[fold_k_indices] == 1)
  n_controls_k <- sum(Y_cv[fold_k_indices] == 0)
  message(sprintf("Fold %d: %d patients (%d cases, %d controls)",
                  k, length(fold_k_indices), n_cases_k, n_controls_k))
}


# ------------------------------------------------------------------------------
# Step 3.3: MCMC Settings for DART
# ------------------------------------------------------------------------------

message("\n=== Step 3.3: MCMC Settings ===")

# MCMC parameters
n_burn <- 2000      # Burn-in iterations (matches writeup: 2000 burn-in)
n_iter <- 1000      # Number of saved posterior draws (M in roadmap)
n_thin <- 5         # Thinning: save every 5th iteration
n_total <- n_burn + (n_iter * n_thin)  # Total iterations: 2000 + 5000 = 7000

# ============================================================================
# DART/BART parameters for HIGH-DIMENSIONAL variable selection
# ============================================================================
# TUNING FOR HIGHER PIPs: Adjust these parameters to increase variable selection
#
# 1. n_trees: FEWER trees = HIGHER PIPs (more concentration per variable)
#    - Default BART: 200 trees
#    - For high PIPs: 50-100 trees
#    - Current: 50 (AGGRESSIVE for high PIPs)
#
# 2. rho: DART concentration parameter (α in Dirichlet prior)
#    - Controls sparsity: SMALLER rho = MORE concentrated selection
#    - Default: rho = D (uniform over all variables)
#    - For sparsity: rho = 10, 50, or 100
#    - Current: 50 (MODERATE sparsity)
#
# 3. theta/omega: Additional sparsity controls
#    - theta: probability that a variable is available for splitting
#    - omega: tuning parameter for theta's prior
#
# 4. a, b: Tree structure prior P(split at depth d) = a(1+d)^(-b)
#    - Larger trees can capture more interactions but dilute PIPs
#    - Smaller b = shallower trees = more concentrated PIPs
# ============================================================================

n_trees <- 200
sparse <- TRUE     # Use sparse (DART) prior
# NOTE: rho will be set dynamically based on number of genes entering BART
# rho = ncol(X_cv) for CV, rho = ncol(X_cv_filtered) for full model
theta <- 0         # theta parameter for sparse prior (0 = automatic selection)
omega <- 1         # omega parameter
a <- 0.5           # Hyper-parameters for the beta distribution
b <- 1             # Hyper-parameters for the beta distribution
numcut <- 100      # Number of cut points for continuous variables
usequants <- FALSE # Use uniform cutpoints (not quantiles)

# Prior for clinical coefficients β
sigma_beta <- sqrt(10)  # β ~ N(0, 10*I)

message(sprintf("MCMC: %d burn-in + %d saved draws (thin=%d) = %d total iterations",
                n_burn, n_iter, n_thin, n_total))
message(sprintf("DART: %d trees, sparse=%s", n_trees, sparse))
message(sprintf("Sparsity parameter rho: Will be set dynamically based on number of genes"))
message(sprintf("Tree prior: a=%.1f, b=%.1f", a, b))
message(sprintf("Prior: β ~ N(0, %.1f*I)", sigma_beta^2))


# ------------------------------------------------------------------------------
# Step 3.4: Cross-Validation Loop
# ------------------------------------------------------------------------------

message("\n=== Step 3.4: Running 5-Fold Cross-Validation ===")

# Initialize storage for CV results
cv_predictions <- matrix(NA, nrow = length(Y_cv), ncol = n_iter)
cv_PIPs_by_fold <- matrix(0, nrow = K, ncol = ncol(X_cv))
colnames(cv_PIPs_by_fold) <- colnames(X_cv)

# Storage for fold-specific performance metrics (per MCMC draw)
cv_AUC_by_draw <- matrix(NA, nrow = K, ncol = n_iter)
cv_Brier_by_draw <- matrix(NA, nrow = K, ncol = n_iter)

# NEW: Storage for fold-specific selected genes (for union criterion)
fold_selected_genes <- vector("list", K)  # List to store selected genes per fold
fold_FDR_info <- vector("list", K)        # Store FDR control details per fold

# Loop over folds
for (k in 1:K) {
  
  message(sprintf("\n--- Fold %d/%d ---", k, K))
  
  # Step 3.4.1: Split data into training and validation
  val_k_indices <- fold_assignment[[k]]
  train_k_indices <- setdiff(1:length(Y_cv), val_k_indices)
  
  # Training data for fold k
  X_train_k <- X_cv[train_k_indices, ]
  Z_train_k <- Z_cv[train_k_indices, ]
  Y_train_k <- Y_cv[train_k_indices]
  
  # Validation data for fold k
  X_val_k <- X_cv[val_k_indices, ]
  Z_val_k <- Z_cv[val_k_indices, ]
  Y_val_k <- Y_cv[val_k_indices]
  
  message(sprintf("  Training: %d patients", length(Y_train_k)))
  message(sprintf("  Validation: %d patients", length(Y_val_k)))
  
  # Step 3.4.2: Fit semi-parametric probit DART model
  # Model: Y ~ Φ(f(X) + Z'β + ε), ε ~ N(0,1)
  #
  # IMPLEMENTATION:
  # We use a custom three-step Gibbs sampler that explicitly separates:
  #   - f(X): DART function on gene expression (X) - nonparametric
  #   - Z'β: Linear effects for clinical covariates (Z) - parametric
  #
  # Three-step Gibbs update:
  #   1) Update latent U | Y, f, β via Albert-Chib truncated normal
  #   2) Update β | U, f via conjugate Gaussian posterior
  #   3) Update f | U, β via DART on residuals (U - Z'β)
  
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
    rho = rho_cv,              # Use rho = number of screened genes
    numcut = numcut,
    usequants = usequants,
    sigma_beta = sigma_beta    # Prior std for clinical coefficients
  )
  
  message("  Model fitting complete.")
  
  # Step 3.4.3: Extract posterior predictions on validation fold
  # bart_fit_k$prob.test contains posterior predictive probabilities
  # Dimensions: n_iter × n_validation
  prob_val_k <- bart_fit_k$prob.test  # n_iter × n_val
  
  # Store predictions for later aggregation
  for (m in 1:n_iter) {
    cv_predictions[val_k_indices, m] <- prob_val_k[m, ]
  }
  
  # Step 3.4.4: Compute per-draw AUC and Brier score for fold k
  for (m in 1:n_iter) {
    p_m <- prob_val_k[m, ]
    
    # AUC for draw m
    if (length(unique(Y_val_k)) == 2) {
      roc_m <- roc(Y_val_k, p_m, quiet = TRUE, levels = c(0, 1), direction = "<")
      cv_AUC_by_draw[k, m] <- as.numeric(auc(roc_m))
    } else {
      cv_AUC_by_draw[k, m] <- NA
    }
    
    # Brier score for draw m
    cv_Brier_by_draw[k, m] <- mean((Y_val_k - p_m)^2)
  }
  
  # Step 3.4.5: Compute fold-specific PIPs
  # Extract variable counts from DART fit
  # bart_fit_k$varcount: n_iter × D matrix of variable usage counts (genes only)
  # Note: Clinical covariates (Z) are handled parametrically via β, not tracked in PIPs
  
  if (!is.null(bart_fit_k$varcount)) {
    # Compute PIP for each gene: proportion of draws where gene was used in DART
    var_used <- bart_fit_k$varcount > 0  # n_iter × D binary matrix
    PIP_k <- colMeans(var_used)  # Average over MCMC draws
    
    # Ensure names are set correctly
    names(PIP_k) <- colnames(X_cv)
    
    # Store fold-specific PIPs (genes only)
    cv_PIPs_by_fold[k, ] <- PIP_k
    
    message(sprintf("  Top 5 genes by PIP in fold %d:", k))
    top5_idx <- order(PIP_k, decreasing = TRUE)[1:5]
    for (i in 1:5) {
      gene_idx <- top5_idx[i]
      message(sprintf("    %s: PIP = %.3f", colnames(X_cv)[gene_idx], PIP_k[gene_idx]))
    }
    
    # NEW: Fold-specific FDR control for variable selection
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
    
    # Find optimal R* for this fold at FDR level α = 0.15
    # NOTE: Using same threshold at fold level and final model for consistency
    alpha_FDR <- 0.05
    valid_R_k <- which(FDR_k <= alpha_FDR)
    
    if (length(valid_R_k) > 0) {
      R_star_k <- max(valid_R_k)
      tau_k <- PIP_k_sorted[R_star_k]
      fold_selected_genes[[k]] <- names(PIP_k_sorted)[1:R_star_k]
      
      message(sprintf("  Fold %d FDR control: R* = %d genes, tau = %.4f, FDR = %.4f",
                      k, R_star_k, tau_k, FDR_k[R_star_k]))
    } else {
      # No genes meet FDR threshold in this fold - keep selection empty
      R_star_k <- 0
      tau_k <- NA
      fold_selected_genes[[k]] <- character(0)
      
      message(sprintf("  Fold %d FDR control: WARNING - No genes meet FDR threshold (alpha = %.2f)", k, alpha_FDR))
      message(sprintf("  Keeping fold selection EMPTY (no fallback)"))
    }
    
    # Store FDR control information for this fold
    fold_FDR_info[[k]] <- list(
      R_star = R_star_k,
      tau = tau_k,
      FDR_hat = if (length(valid_R_k) > 0) FDR_k[R_star_k] else NA,
      selected_genes = fold_selected_genes[[k]],
      all_FDR = FDR_k,
      all_PIPs = PIP_k_sorted
    )
  } else {
    # If varcount is NULL, initialize empty selection for this fold
    fold_selected_genes[[k]] <- character(0)
    fold_FDR_info[[k]] <- list(
      R_star = 0,
      tau = NA,
      FDR_hat = NA,
      selected_genes = character(0),
      all_FDR = numeric(0),
      all_PIPs = numeric(0)
    )
    message(sprintf("  WARNING: No variable counts available for fold %d", k))
  }
  
  message(sprintf("  Fold %d complete.", k))
}

message("\n=== 5-Fold Cross-Validation Complete ===")


# ------------------------------------------------------------------------------
# Step 3.5: Aggregate CV Performance Metrics
# ------------------------------------------------------------------------------

message("\n=== Step 3.5: Aggregating CV Performance Metrics ===")

# Posterior-averaged predictions for each patient in CV pool
cv_predictions_mean <- rowMeans(cv_predictions, na.rm = TRUE)

# Point estimates: AUC and Brier on full CV pool
cv_roc <- roc(Y_cv, cv_predictions_mean, quiet = TRUE, levels = c(0, 1), direction = "<")
cv_AUC_point <- as.numeric(auc(cv_roc))
cv_Brier_point <- mean((Y_cv - cv_predictions_mean)^2)

message(sprintf("CV Performance (point estimates):"))
message(sprintf("  AUC: %.4f", cv_AUC_point))
message(sprintf("  Brier Score: %.4f", cv_Brier_point))

# Uncertainty quantification: posterior distributions
# Flatten across folds and draws
cv_AUC_all_draws <- as.vector(cv_AUC_by_draw)
cv_Brier_all_draws <- as.vector(cv_Brier_by_draw)

# Remove NAs
cv_AUC_all_draws <- cv_AUC_all_draws[!is.na(cv_AUC_all_draws)]
cv_Brier_all_draws <- cv_Brier_all_draws[!is.na(cv_Brier_all_draws)]

# Compute credible intervals
cv_AUC_mean <- mean(cv_AUC_all_draws)
cv_AUC_CrI <- quantile(cv_AUC_all_draws, c(0.025, 0.975))

cv_Brier_mean <- mean(cv_Brier_all_draws)
cv_Brier_CrI <- quantile(cv_Brier_all_draws, c(0.025, 0.975))

message(sprintf("\nCV Performance (posterior summaries):"))
message(sprintf("  AUC: %.4f [%.4f, %.4f]", cv_AUC_mean, cv_AUC_CrI[1], cv_AUC_CrI[2]))
message(sprintf("  Brier: %.4f [%.4f, %.4f]", cv_Brier_mean, cv_Brier_CrI[1], cv_Brier_CrI[2]))

# Save CV results
saveRDS(list(
  predictions_mean = cv_predictions_mean,
  predictions_full = cv_predictions,
  AUC_point = cv_AUC_point,
  Brier_point = cv_Brier_point,
  AUC_draws = cv_AUC_all_draws,
  Brier_draws = cv_Brier_all_draws,
  AUC_mean = cv_AUC_mean,
  AUC_CrI = cv_AUC_CrI,
  Brier_mean = cv_Brier_mean,
  Brier_CrI = cv_Brier_CrI
), result_path("cv_performance_metrics.RDS"))


# ------------------------------------------------------------------------------
# Step 3.6: Compute Union of Selected Genes Across Folds
# ------------------------------------------------------------------------------

message("\n=== Step 3.6: Computing Union of Selected Genes ===")

# METHODOLOGY: Retain genes selected in ANY fold
# This ensures comprehensive coverage - a gene is kept if it shows signal in at least one fold

# Display per-fold selections
message("Per-fold selections:")
for (k in 1:K) {
  n_selected <- length(fold_selected_genes[[k]])
  message(sprintf("  Fold %d: %d genes selected", k, n_selected))
}

# Compute union: genes selected in ANY of the K folds
if (any(sapply(fold_selected_genes, length) > 0)) {
  # Start with empty set
  union_genes <- character(0)
  
  # Add genes from each fold (automatically handles duplicates)
  for (k in 1:K) {
    union_genes <- union(union_genes, fold_selected_genes[[k]])
  }
  
  message(sprintf("\nUnion criterion:"))
  message(sprintf("  Genes selected in ANY of the %d folds: %d genes", K, length(union_genes)))
  
  # Compute frequency: how many folds selected each gene
  gene_fold_frequency <- sapply(union_genes, function(gene) {
    sum(sapply(fold_selected_genes, function(fold_genes) gene %in% fold_genes))
  })
  
  # Sort genes by frequency (most frequently selected first)
  union_genes_sorted <- union_genes[order(gene_fold_frequency, decreasing = TRUE)]
  gene_fold_frequency_sorted <- gene_fold_frequency[order(gene_fold_frequency, decreasing = TRUE)]
  
  if (length(union_genes) > 0) {
    message(sprintf("\nGenes selected in union (showing up to 50, sorted by frequency):"))
    n_display <- min(50, length(union_genes_sorted))
    for (i in 1:n_display) {
      gene_name <- union_genes_sorted[i]
      n_folds <- gene_fold_frequency_sorted[i]
      # Compute average PIP across folds for this gene
      avg_PIP <- mean(cv_PIPs_by_fold[, gene_name])
      message(sprintf("  %2d. %s: selected in %d/%d folds, avg PIP = %.4f",
                      i, gene_name, n_folds, K, avg_PIP))
    }
    if (length(union_genes_sorted) > 50) {
      message(sprintf("  ... and %d more genes", length(union_genes_sorted) - 50))
    }
  }
  
  # Use union genes for downstream analysis
  union_genes <- union_genes  # Keep variable name for compatibility with downstream code

} else {
  union_genes <- character(0)
  warning("No folds have selected genes. Union is empty.")
  message("This should not happen if FDR control is working properly.")
}

# IMPORTANT: Ensure study_indicator is included in union_genes for final prognostic model
# Study indicator must be part of the final model regardless of whether it was selected during CV
if ("study_indicator" %in% union_genes) {
  message("\n  study_indicator was selected during CV - already in union_genes")
  study_indicator_selected_in_cv <- TRUE
} else {
  message("\n  study_indicator was NOT selected during CV - adding to union_genes for final model")
  union_genes <- c(union_genes, "study_indicator")
  study_indicator_selected_in_cv <- FALSE
}
message(sprintf("  Final union_genes includes study_indicator: %d total variables", length(union_genes)))

# Save union and fold-specific selections
saveRDS(list(
  union_genes = union_genes,  # Union genes (using union_genes variable name for compatibility)
  fold_selected_genes = fold_selected_genes,
  fold_FDR_info = fold_FDR_info,
  PIP_by_fold = cv_PIPs_by_fold,
  gene_fold_frequency = gene_fold_frequency_sorted,  # How many folds selected each gene
  union_genes_sorted = union_genes_sorted  # Genes sorted by frequency
), result_path("cv_gene_selection.RDS"))


# ------------------------------------------------------------------------------
# Step 3.7: Refit Model on Full Development Set with Union Genes
# ------------------------------------------------------------------------------

message("\n=== Step 3.7: Refitting Model with Union Genes ===")

# METHODOLOGY: Refit model on full development set using union genes (selected in ANY fold)
# This provides full-data PIPs for final FDR control and reporting
# IMPORTANT: Final prognostic model uses ALL clinical covariates (race, bmi, sex, age)
#            whereas CV used only (sex, age)

if (length(union_genes) > 0) {
  message(sprintf("Refitting semi-parametric model with %d union genes...",
                  length(union_genes)))
  message(sprintf("Using FULL clinical covariate set for final model: race_binary, bmi_binary, sex_binary, age"))

  # Subset X to only include union genes
  X_cv_filtered <- X_cv[, union_genes, drop = FALSE]
  X_test_filtered <- X_test[, union_genes, drop = FALSE]

  # Set rho for full model based on number of union genes entering BART
  rho_full <- ncol(X_cv_filtered)
  message(sprintf("Setting rho for full model: rho = %d (union genes + 1 for study_indicator in X_cv_filtered)", rho_full))

  # Fit final model on full development set with filtered genes
  # NOTE: Using Z_cv_full and Z_test_full (all 4 clinical covariates) instead of Z_cv/Z_test
  bart_fit_full <- fit_semiparametric_probit_DART(
    X = X_cv_filtered,         # ONLY union genes
    Z = Z_cv_full,             # ALL clinical covariates (race, bmi, sex, age)
    Y = Y_cv,
    X_test = X_test_filtered,  # ONLY union genes
    Z_test = Z_test_full,      # ALL clinical covariates (race, bmi, sex, age)
    n_burn = n_burn,
    n_iter = n_iter,
    n_thin = n_thin,
    n_trees = n_trees,
    sparse = sparse,
    theta = theta,
    omega = omega,
    a = a,
    b = b,
    rho = rho_full,            # Use rho = number of union genes
    numcut = numcut,
    usequants = usequants,
    sigma_beta = sigma_beta
  )
  
  message("Full model fitting complete.")
  
  # Diagnostic: Check dimensions of fitted model outputs
  message("\n=== Model Output Diagnostics ===")
  if (!is.null(bart_fit_full$prob.test)) {
    message(sprintf("  prob.test dimensions: %d × %d (n_iter × n_test)",
                    nrow(bart_fit_full$prob.test), ncol(bart_fit_full$prob.test)))
    message(sprintf("  Expected n_test: %d (Y_test length)", length(Y_test)))
    if (ncol(bart_fit_full$prob.test) != length(Y_test)) {
      warning("DIMENSION MISMATCH: prob.test columns != Y_test length!")
    }
  } else {
    warning("prob.test is NULL - test predictions were not generated!")
  }
  
  # Extract full-data PIPs for union genes
  if (!is.null(bart_fit_full$varcount)) {
    var_used_full <- bart_fit_full$varcount > 0  # n_iter × length(union_genes)
    PIP_full <- colMeans(var_used_full)
    names(PIP_full) <- union_genes
    
    # Sort by decreasing PIP
    PIP_full_sorted <- sort(PIP_full, decreasing = TRUE)
    
    message("\nFull-data PIPs for union genes:")
    for (i in 1:min(20, length(PIP_full_sorted))) {
      message(sprintf("  %2d. %s: PIP = %.4f", i, names(PIP_full_sorted)[i],
                      PIP_full_sorted[i]))
    }
    
    # Apply Bayesian FDR control to full-data PIPs
    message("\n=== Applying Bayesian FDR Control to Full-Data PIPs ===")
    
    alpha_FDR <- 0.05  # FDR control for final prognostic model 
    D_full <- length(PIP_full_sorted)
    FDR_hat_full <- numeric(D_full)
    
    for (k_genes in 1:D_full) {
      top_K_PIPs <- PIP_full_sorted[1:k_genes]
      FDR_hat_full[k_genes] <- sum(1 - top_K_PIPs) / k_genes
    }
    
    # Diagnostic: Show FDR for top genes
    message(sprintf("  FDR diagnostics for full-data model:"))
    for (r in c(1, 5, 10, 20, 50, 100)) {
      if (r <= D_full) {
        message(sprintf("    Top %d genes: FDR = %.4f, min PIP = %.4f",
                        r, FDR_hat_full[r], PIP_full_sorted[r]))
      }
    }
    
    # Find optimal K
    valid_K_full <- which(FDR_hat_full <= alpha_FDR)
    
    if (length(valid_K_full) > 0) {
      K_star_full <- max(valid_K_full)
      tau_full <- PIP_full_sorted[K_star_full]
      selected_genes_final <- names(PIP_full_sorted)[1:K_star_full]
      
      message(sprintf("FDR control at α = %.2f:", alpha_FDR))
      message(sprintf("  Optimal K: %d genes", K_star_full))
      message(sprintf("  PIP threshold: τ = %.4f", tau_full))
      message(sprintf("  Expected FDR: %.4f", FDR_hat_full[K_star_full]))
    } else {
      # All union genes have too low PIPs - use all of them
      warning("No union genes meet FDR threshold. Using all union genes.")
      K_star_full <- D_full
      tau_full <- PIP_full_sorted[K_star_full]
      selected_genes_final <- names(PIP_full_sorted)[1:K_star_full]
      
      message(sprintf("FDR control at α = %.2f:", alpha_FDR))
      message(sprintf("  WARNING: No genes meet FDR threshold"))
      message(sprintf("  Using all %d union genes", K_star_full))
    }
    
    # Save final selected genes
    saveRDS(list(
      selected_genes = selected_genes_final,
      K_star = K_star_full,
      tau = tau_full,
      FDR_hat = if (length(valid_K_full) > 0) FDR_hat_full[K_star_full] else NA,
      all_FDR = FDR_hat_full,
      all_PIPs_full = PIP_full_sorted,
      union_genes = union_genes  # Union genes used for full refit
    ), result_path("selected_genes_FDR_control.RDS"))
    
    # ========================================================================
    # Visualization: Posterior Distributions of Gene Usage (PIPs)
    # ========================================================================
    message("\n=== Creating PIP Posterior Distribution Visualizations ===")

    if (K_star_full > 0) {
      # Extract variable usage counts for selected genes
      selected_gene_indices <- match(selected_genes_final, union_genes)
      varcount_selected <- bart_fit_full$varcount[, selected_gene_indices, drop = FALSE]
      colnames(varcount_selected) <- selected_genes_final

      # IMPORTANT: Also extract study_indicator for separate plotting
      # Study indicator is always included in union_genes and final model
      study_indicator_idx <- which(union_genes == "study_indicator")
      study_indicator_varcount <- bart_fit_full$varcount[, study_indicator_idx]
      study_indicator_pip <- PIP_full["study_indicator"]

      # Create list of variables to plot: selected genes + study_indicator
      # Ensure study_indicator is included even if not in selected_genes_final
      if ("study_indicator" %in% selected_genes_final) {
        vars_to_plot <- selected_genes_final
        message("  study_indicator is among FDR-selected genes")
      } else {
        vars_to_plot <- c(selected_genes_final, "study_indicator")
        message("  study_indicator added to PIP plots (not FDR-selected, but included in model)")
      }

      # Determine grid layout: Fixed 3 rows × 4 columns for consistent sizing
      n_vars <- length(vars_to_plot)
      n_cols <- 4  # Fixed: 4 columns across all pages
      n_rows <- 3  # Fixed: 3 rows across all pages
      vars_per_page <- n_rows * n_cols  # 12 variables per page
      n_pages <- ceiling(n_vars / vars_per_page)

      message(sprintf("Creating posterior distribution plots for %d variables (%d genes + study_indicator) across %d pages",
                      n_vars, length(selected_genes_final), n_pages))
      message(sprintf("  Layout: %d rows × %d columns per page (fixed for consistent sizing)", n_rows, n_cols))

      pdf(result_path("posterior_pip_distributions.pdf"), width = 16, height = 12)

      for (page in 1:n_pages) {
        start_idx <- (page - 1) * vars_per_page + 1
        end_idx <- min(page * vars_per_page, n_vars)
        vars_this_page <- vars_to_plot[start_idx:end_idx]
        n_vars_this_page <- length(vars_this_page)

        # Always use fixed layout (3 rows × 4 columns) for consistent plot sizes
        par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1), oma = c(2, 2, 4, 1))

        for (var_name in vars_this_page) {
          # Get variable usage counts and PIP
          var_idx <- which(union_genes == var_name)
          var_usage <- bart_fit_full$varcount[, var_idx]
          var_pip <- PIP_full[var_name]

          # Display name: use "STUDY INDICATOR" label to distinguish it
          if (var_name == "study_indicator") {
            display_name <- "STUDY INDICATOR"
          } else {
            display_name <- var_name
          }

          # Create histogram of usage counts (identical styling for all variables)
          hist(var_usage,
               breaks = 30,
               col = "steelblue",
               border = "white",
               main = sprintf("%s\n(PIP = %.3f)", display_name, var_pip),
               xlab = "Variable Usage Count per Iteration",
               ylab = "Frequency",
               cex.main = 1.0)

          # Add vertical line at mean
          abline(v = mean(var_usage), col = "red", lwd = 2, lty = 2)

          # Add text showing mean and PIP
          legend("topright",
                 legend = c(sprintf("Mean = %.1f", mean(var_usage)),
                            sprintf("PIP = %.3f", var_pip)),
                 bty = "n", cex = 0.8)

          grid(col = "gray80", lty = "dotted")
        }

        # Fill remaining slots with empty plots to maintain consistent layout
        n_empty <- vars_per_page - n_vars_this_page
        if (n_empty > 0) {
          for (i in 1:n_empty) {
            plot.new()  # Empty placeholder plot
          }
        }

        # Add overall title
        mtext("Posterior Distributions of Variable Usage (PIPs)",
              outer = TRUE, cex = 1.4, font = 2, line = 1.5)
        mtext(sprintf("(Page %d of %d) - Includes selected genes + study_indicator", page, n_pages),
              outer = TRUE, cex = 1.0, line = 0.2)
      }

      dev.off()

      message(sprintf("PIP posterior distributions saved to: %s",
                      result_path("posterior_pip_distributions.pdf")))
      message(sprintf("  study_indicator PIP = %.4f", study_indicator_pip))
    } else {
      message("No genes selected. Skipping PIP distribution visualization.")
    }
    
    # ========================================================================
    # Convergence Diagnostics: Traceplots for Selected Genes and Covariates
    # ========================================================================
    message("\n=== Step 3.8: MCMC Convergence Diagnostics ===")
    
    # Select 5 significant genes randomly from top genes
    set.seed(123)  # For reproducibility
    if (K_star_full >= 5) {
      sampled_gene_indices <- sample(1:min(20, K_star_full), size = 5, replace = FALSE)
      sampled_genes <- names(PIP_full_sorted)[sampled_gene_indices]
    } else {
      sampled_genes <- names(PIP_full_sorted)[1:K_star_full]
    }
    
    # Select 1 clinical covariate (just select the first one)
    # NOTE: Use Z_cv_full since final model uses all 4 clinical covariates
    beta_names <- colnames(Z_cv_full)
    sampled_covariate_idx <- 1
    sampled_covariate <- beta_names[sampled_covariate_idx]
    
    message(sprintf("\nCreating traceplots for MCMC convergence diagnostics:"))
    message(sprintf("  Sampled genes (%d):", length(sampled_genes)))
    for (gene in sampled_genes) {
      message(sprintf("    - %s (PIP = %.4f)", gene, PIP_full[gene]))
    }
    message(sprintf("  Sampled covariate: %s", sampled_covariate))
    
    # Create traceplots
    pdf(result_path("mcmc_traceplots.pdf"), width = 12, height = 10)
    
    # Layout: 3 rows × 2 columns (5 genes + 1 covariate)
    par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
    
    # Plot 1-5: Gene usage traceplots (binary indicator)
    for (i in 1:length(sampled_genes)) {
      gene <- sampled_genes[i]
      gene_idx <- which(union_genes == gene)
      gene_usage <- bart_fit_full$varcount[, gene_idx]  # n_iter × 1 binary vector
      
      plot(1:n_iter, gene_usage, type = "l",
           xlab = "MCMC Iteration (after burn-in)",
           ylab = "Variable Usage Count",
           main = sprintf("Gene: %s (PIP = %.4f)", gene, PIP_full[gene]),
           col = "darkblue", lwd = 0.8)
      abline(h = mean(gene_usage), col = "red", lty = 2, lwd = 2)
      legend("topright", legend = sprintf("Mean = %.2f", mean(gene_usage)),
             col = "red", lty = 2, lwd = 2, bty = "n")
    }
    
    # Plot 6: Clinical covariate traceplot
    beta_chain <- bart_fit_full$beta_draws[, sampled_covariate_idx]
    
    plot(1:n_iter, beta_chain, type = "l",
         xlab = "MCMC Iteration (after burn-in)",
         ylab = "Coefficient Value",
         main = sprintf("Covariate: %s", sampled_covariate),
         col = "darkgreen", lwd = 0.8)
    abline(h = mean(beta_chain), col = "red", lty = 2, lwd = 2)
    legend("topright", legend = sprintf("Mean = %.4f", mean(beta_chain)),
           col = "red", lty = 2, lwd = 2, bty = "n")
    
    dev.off()
    
    message(sprintf("\nTraceplots saved to: %s", result_path("mcmc_traceplots.pdf")))
    
    # Print summary statistics
    message("\n=== MCMC Summary Statistics ===")
    message("\nGene Usage Statistics:")
    for (gene in sampled_genes) {
      gene_idx <- which(union_genes == gene)
      gene_usage <- bart_fit_full$varcount[, gene_idx]
      message(sprintf("  %s: Mean = %.2f, SD = %.2f, PIP = %.4f",
                      gene, mean(gene_usage), sd(gene_usage), PIP_full[gene]))
    }
    
    message(sprintf("\nClinical Covariate Statistics (%s):", sampled_covariate))
    message(sprintf("  Mean = %.4f", mean(beta_chain)))
    message(sprintf("  SD = %.4f", sd(beta_chain)))
    message(sprintf("  95%% Credible Interval: [%.4f, %.4f]",
                    quantile(beta_chain, 0.025), quantile(beta_chain, 0.975)))
    
  } else {
    stop("bart_fit_full$varcount is NULL - cannot compute PIPs")
  }
  
} else {
  # No union genes - cannot proceed with this methodology
  warning("No union genes available. Skipping full model refit.")
  message("The union is empty. This should not happen if FDR control is working.")
  message("Consider:")
  message("  1. More aggressive tuning (n_trees = 20-30)")
  message("  2. Checking fold-specific FDR control settings")
  message("  3. Verifying screening step is working correctly")
  
  # Create empty placeholder
  bart_fit_full <- NULL
  selected_genes_final <- character(0)
}


# ------------------------------------------------------------------------------
# Step 3.8: Extract Clinical Coefficients from Full Model
# ------------------------------------------------------------------------------

message("\n=== Step 3.8: Extracting Clinical Coefficients ===")

# Extract and save posterior draws of clinical coefficients β
# (Only if full model was fitted successfully)
if (!is.null(bart_fit_full) && !is.null(bart_fit_full$beta_draws)) {
  message("\nPosterior estimates for clinical coefficients (β):")
  beta_posterior_mean <- colMeans(bart_fit_full$beta_draws)
  beta_posterior_sd <- apply(bart_fit_full$beta_draws, 2, sd)
  beta_credible_intervals <- apply(bart_fit_full$beta_draws, 2, function(x) quantile(x, c(0.025, 0.975)))
  
  # Determine beta names based on actual number of columns in beta_draws
  n_beta_cols <- ncol(bart_fit_full$beta_draws)
  n_Z_cols <- ncol(Z_cv_full)  # Use Z_cv_full since final model uses all 4 clinical covariates

  message(sprintf("  beta_draws has %d columns, Z_cv_full has %d columns", n_beta_cols, n_Z_cols))

  # Check if intercept is included (beta_draws has one more column than Z_cv_full)
  if (n_beta_cols == n_Z_cols + 1) {
    # Intercept is first column
    beta_names <- c("intercept", colnames(Z_cv_full))
    message("  Intercept detected in beta_draws")
  } else if (n_beta_cols == n_Z_cols) {
    # No intercept, just the covariates
    beta_names <- colnames(Z_cv_full)
    message("  No intercept in beta_draws")
  } else {
    # Fallback: generic names
    beta_names <- paste0("beta_", 1:n_beta_cols)
    warning(sprintf("Unexpected beta_draws dimensions: %d cols vs %d Z_cv_full cols. Using generic names.",
                    n_beta_cols, n_Z_cols))
  }
  
  beta_summary <- data.frame(
    Covariate = beta_names,
    Mean = beta_posterior_mean,
    SD = beta_posterior_sd,
    CI_Lower = beta_credible_intervals[1, ],
    CI_Upper = beta_credible_intervals[2, ]
  )
  
  print(beta_summary)
  saveRDS(bart_fit_full$beta_draws, result_path("beta_posterior_draws.RDS"))
  saveRDS(beta_summary, result_path("beta_summary.RDS"))
} else {
  warning("Full model was not fitted. Skipping clinical coefficient extraction.")
}


# ------------------------------------------------------------------------------
# Step 3.9: External Test Set Evaluation
# ------------------------------------------------------------------------------

message("\n=== Step 3.9: External Test Set Evaluation ===")

# Extract posterior predictions on test set
# (Only if full model was fitted successfully)
if (!is.null(bart_fit_full) && !is.null(bart_fit_full$prob.test)) {
  # bart_fit_full$prob.test: n_iter × n_test
  prob_test_full <- bart_fit_full$prob.test
  
  # Compute per-draw AUC and Brier on test set
  test_AUC_by_draw <- numeric(n_iter)
  test_Brier_by_draw <- numeric(n_iter)
  
  for (m in 1:n_iter) {
    p_m <- prob_test_full[m, ]
    
    # AUC
    if (length(unique(Y_test)) == 2) {
      roc_m <- roc(Y_test, p_m, quiet = TRUE, levels = c(0, 1), direction = "<")
      test_AUC_by_draw[m] <- as.numeric(auc(roc_m))
    } else {
      test_AUC_by_draw[m] <- NA
    }
    
    # Brier
    test_Brier_by_draw[m] <- mean((Y_test - p_m)^2)
  }
  
  # Posterior summaries
  test_AUC_mean <- mean(test_AUC_by_draw, na.rm = TRUE)
  test_AUC_CrI <- quantile(test_AUC_by_draw, c(0.025, 0.975), na.rm = TRUE)
  
  test_Brier_mean <- mean(test_Brier_by_draw)
  test_Brier_CrI <- quantile(test_Brier_by_draw, c(0.025, 0.975))
  
  message(sprintf("\nTest Set Performance (N=%d):", length(Y_test)))
  message(sprintf("  AUC: %.4f [%.4f, %.4f]", test_AUC_mean, test_AUC_CrI[1], test_AUC_CrI[2]))
  message(sprintf("  Brier: %.4f [%.4f, %.4f]", test_Brier_mean, test_Brier_CrI[1], test_Brier_CrI[2]))
  
  # Posterior-averaged predictions
  # prob_test_full is n_iter × n_test, so colMeans gives one prediction per test patient
  test_predictions_mean <- colMeans(prob_test_full)
  
  # Save test results
  saveRDS(list(
    predictions_mean = test_predictions_mean,
    predictions_full = prob_test_full,
    AUC_draws = test_AUC_by_draw,
    Brier_draws = test_Brier_by_draw,
    AUC_mean = test_AUC_mean,
    AUC_CrI = test_AUC_CrI,
    Brier_mean = test_Brier_mean,
    Brier_CrI = test_Brier_CrI,
    Y_test = Y_test
  ), result_path("test_set_performance.RDS"))
} else {
  warning("Full model was not fitted. Skipping test set evaluation.")
  message("Test set evaluation requires successful model fit with union genes.")
}

message("\n=== Semi-Parametric Probit DART Analysis Complete! ===")
message("Results saved to: ", results_dir)


################################################################################
# 4) Visualization of Results
################################################################################

message("\n=== Generating Visualizations ===")

# Only generate visualizations if full model was successfully fitted AND required variables exist
viz_ready <- (!is.null(bart_fit_full) &&
                exists("test_predictions_mean") && exists("cv_predictions_mean") &&
                exists("test_AUC_mean") && exists("cv_AUC_mean") &&
                exists("test_AUC_CrI") && exists("cv_AUC_CrI"))

if (viz_ready) {
  
  # Diagnostic checks
  message(sprintf("Dimension checks before ROC:"))
  message(sprintf("  Y_test length: %d", length(Y_test)))
  message(sprintf("  test_predictions_mean length: %d", length(test_predictions_mean)))
  message(sprintf("  Y_cv length: %d", length(Y_cv)))
  message(sprintf("  cv_predictions_mean length: %d", length(cv_predictions_mean)))
  
  # Verify dimensions match
  if (length(Y_test) != length(test_predictions_mean)) {
    warning(sprintf("Dimension mismatch: Y_test (%d) != test_predictions_mean (%d). Skipping test ROC.",
                    length(Y_test), length(test_predictions_mean)))
    skip_test_roc <- TRUE
  } else {
    skip_test_roc <- FALSE
  }
  
  if (length(Y_cv) != length(cv_predictions_mean)) {
    warning(sprintf("Dimension mismatch: Y_cv (%d) != cv_predictions_mean (%d). Skipping CV ROC.",
                    length(Y_cv), length(cv_predictions_mean)))
    skip_cv_roc <- TRUE
  } else {
    skip_cv_roc <- FALSE
  }
  
  # ------------------------------------------------------------------------------
  # 1. ROC Curves for CV and Test Sets
  # ------------------------------------------------------------------------------
  
  if (!skip_cv_roc || !skip_test_roc) {
    message("Creating ROC curves...")
    
    pdf(result_path("roc_curves.pdf"), width = 12, height = 6)
    
    if (!skip_cv_roc && !skip_test_roc) {
      par(mfrow = c(1, 2))
    } else {
      par(mfrow = c(1, 1))
    }
    
    # CV ROC curve
    if (!skip_cv_roc) {
      cv_roc <- roc(Y_cv, cv_predictions_mean, quiet = TRUE, levels = c(0, 1), direction = "<")
      plot(cv_roc, main = "ROC Curve: Cross-Validation",
           col = "blue", lwd = 2,
           xlab = "False Positive Rate (1 - Specificity)",
           ylab = "True Positive Rate (Sensitivity)")
      text(0.6, 0.2, sprintf("AUC = %.3f [%.3f, %.3f]",
                             cv_AUC_mean, cv_AUC_CrI[1], cv_AUC_CrI[2]),
           col = "blue", cex = 1.2)
      abline(0, 1, lty = 2, col = "gray")
      grid()
    }
    
    # Test ROC curve
    if (!skip_test_roc) {
      test_roc <- roc(Y_test, test_predictions_mean, quiet = TRUE, levels = c(0, 1), direction = "<")
      plot(test_roc, main = "ROC Curve: Test Set",
           col = "red", lwd = 2,
           xlab = "False Positive Rate (1 - Specificity)",
           ylab = "True Positive Rate (Sensitivity)")
      text(0.6, 0.2, sprintf("AUC = %.3f [%.3f, %.3f]",
                             test_AUC_mean, test_AUC_CrI[1], test_AUC_CrI[2]),
           col = "red", cex = 1.2)
      abline(0, 1, lty = 2, col = "gray")
      grid()
    }
    
    dev.off()
  } else {
    message("Skipping ROC curves due to dimension mismatches.")
  }


  # ------------------------------------------------------------------------------
  # 1b. Performance Metrics Across Folds (Manuscript-Ready Grid Plot)
  # ------------------------------------------------------------------------------
  # 1x2 grid: AUC (left) and Brier Score (right)
  # X-axis: Fold 1-5 + Validation Set
  # Y-axis: Posterior mean with 95% credible interval band
  # Colors: Black, white, gray only

  message("Creating manuscript-ready performance metrics grid plot...")

  # Compute fold-specific posterior summaries
  fold_AUC_means <- numeric(K)
  fold_AUC_lower <- numeric(K)
  fold_AUC_upper <- numeric(K)
  fold_Brier_means <- numeric(K)
  fold_Brier_lower <- numeric(K)
  fold_Brier_upper <- numeric(K)

  for (k in 1:K) {
    # AUC summaries
    fold_AUC_draws <- cv_AUC_by_draw[k, ]
    fold_AUC_draws <- fold_AUC_draws[!is.na(fold_AUC_draws)]
    fold_AUC_means[k] <- mean(fold_AUC_draws)
    fold_AUC_CrI <- quantile(fold_AUC_draws, c(0.025, 0.975))
    fold_AUC_lower[k] <- fold_AUC_CrI[1]
    fold_AUC_upper[k] <- fold_AUC_CrI[2]

    # Brier summaries
    fold_Brier_draws <- cv_Brier_by_draw[k, ]
    fold_Brier_draws <- fold_Brier_draws[!is.na(fold_Brier_draws)]
    fold_Brier_means[k] <- mean(fold_Brier_draws)
    fold_Brier_CrI <- quantile(fold_Brier_draws, c(0.025, 0.975))
    fold_Brier_lower[k] <- fold_Brier_CrI[1]
    fold_Brier_upper[k] <- fold_Brier_CrI[2]
  }

  # Combine with validation set (test set) metrics
  all_AUC_means <- c(fold_AUC_means, test_AUC_mean)
  all_AUC_lower <- c(fold_AUC_lower, test_AUC_CrI[1])
  all_AUC_upper <- c(fold_AUC_upper, test_AUC_CrI[2])

  all_Brier_means <- c(fold_Brier_means, test_Brier_mean)
  all_Brier_lower <- c(fold_Brier_lower, test_Brier_CrI[1])
  all_Brier_upper <- c(fold_Brier_upper, test_Brier_CrI[2])

  # X-axis positions and labels
  x_positions <- 1:6
  x_labels <- c("Fold 1", "Fold 2", "Fold 3", "Fold 4", "Fold 5", "Validation")

  # Create PDF
  pdf(result_path("performance_metrics_grid.pdf"), width = 10, height = 5)

  # Set up 1x2 grid layout with proper margins for manuscript
  par(mfrow = c(1, 2),
      mar = c(5, 5, 3, 1.5),    # bottom, left, top, right
      oma = c(0, 0, 2, 0),       # outer margins for overall title
      family = "serif")          # serif font for manuscript

  # ============================================================================
  # Panel A: AUC
  # ============================================================================

  # Determine y-axis range with padding
  auc_ylim <- c(min(all_AUC_lower) - 0.02, max(all_AUC_upper) + 0.02)

  # Create empty plot
  plot(x_positions, all_AUC_means,
       type = "n",
       xlim = c(0.5, 6.5),
       ylim = auc_ylim,
       xlab = "",
       ylab = "AUC",
       xaxt = "n",
       las = 1,
       cex.lab = 1.3,
       cex.axis = 1.1)

  # Add x-axis labels
  axis(1, at = x_positions, labels = x_labels, cex.axis = 1.0, las = 2)
  mtext("Dataset", side = 1, line = 4, cex = 1.1)

  # Add credible interval band (polygon)
  polygon(c(x_positions, rev(x_positions)),
          c(all_AUC_lower, rev(all_AUC_upper)),
          col = "gray85", border = NA)

  # Add connecting line
  lines(x_positions, all_AUC_means, col = "gray40", lwd = 2)

  # Add points for posterior means
  points(x_positions, all_AUC_means, pch = 19, col = "black", cex = 1.5)

  # Add credible interval error bars
  for (i in 1:6) {
    segments(x_positions[i], all_AUC_lower[i],
             x_positions[i], all_AUC_upper[i],
             col = "black", lwd = 1.5)
    # Add caps to error bars
    segments(x_positions[i] - 0.1, all_AUC_lower[i],
             x_positions[i] + 0.1, all_AUC_lower[i],
             col = "black", lwd = 1.5)
    segments(x_positions[i] - 0.1, all_AUC_upper[i],
             x_positions[i] + 0.1, all_AUC_upper[i],
             col = "black", lwd = 1.5)
  }

  # Add subtle grid
  abline(h = seq(floor(auc_ylim[1] * 10) / 10, ceiling(auc_ylim[2] * 10) / 10, by = 0.05),
         col = "gray90", lty = 1)

  # Add panel label
  mtext("A", side = 3, line = 0.5, at = 0.5, cex = 1.5, font = 2)

  # Add box
  box(lwd = 1)

  # ============================================================================
  # Panel B: Brier Score
  # ============================================================================

  # Determine y-axis range with padding
  brier_ylim <- c(min(all_Brier_lower) - 0.01, max(all_Brier_upper) + 0.01)

  # Create empty plot
  plot(x_positions, all_Brier_means,
       type = "n",
       xlim = c(0.5, 6.5),
       ylim = brier_ylim,
       xlab = "",
       ylab = "Brier Score",
       xaxt = "n",
       las = 1,
       cex.lab = 1.3,
       cex.axis = 1.1)

  # Add x-axis labels
  axis(1, at = x_positions, labels = x_labels, cex.axis = 1.0, las = 2)
  mtext("Dataset", side = 1, line = 4, cex = 1.1)

  # Add credible interval band (polygon)
  polygon(c(x_positions, rev(x_positions)),
          c(all_Brier_lower, rev(all_Brier_upper)),
          col = "gray85", border = NA)

  # Add connecting line
  lines(x_positions, all_Brier_means, col = "gray40", lwd = 2)

  # Add points for posterior means
  points(x_positions, all_Brier_means, pch = 19, col = "black", cex = 1.5)

  # Add credible interval error bars
  for (i in 1:6) {
    segments(x_positions[i], all_Brier_lower[i],
             x_positions[i], all_Brier_upper[i],
             col = "black", lwd = 1.5)
    # Add caps to error bars
    segments(x_positions[i] - 0.1, all_Brier_lower[i],
             x_positions[i] + 0.1, all_Brier_lower[i],
             col = "black", lwd = 1.5)
    segments(x_positions[i] - 0.1, all_Brier_upper[i],
             x_positions[i] + 0.1, all_Brier_upper[i],
             col = "black", lwd = 1.5)
  }

  # Add subtle grid
  abline(h = seq(floor(brier_ylim[1] * 100) / 100, ceiling(brier_ylim[2] * 100) / 100, by = 0.02),
         col = "gray90", lty = 1)

  # Add panel label
  mtext("B", side = 3, line = 0.5, at = 0.5, cex = 1.5, font = 2)

  # Add box
  box(lwd = 1)

  # Add overall title
  mtext("Model Performance Across Cross-Validation Folds and Validation Set",
        outer = TRUE, cex = 1.3, font = 2, line = 0.8)

  # Add subtitle with sample size information
  mtext("5-Fold CV: 400 training + 100 validation per fold | Held-out validation set: n = 369",
        outer = TRUE, cex = 0.95, font = 1, line = -0.3)

  dev.off()

  message("  Saved: performance_metrics_grid.pdf")
  message(sprintf("  AUC range: %.3f - %.3f", min(all_AUC_means), max(all_AUC_means)))
  message(sprintf("  Brier range: %.4f - %.4f", min(all_Brier_means), max(all_Brier_means)))


  # ------------------------------------------------------------------------------
  # 2. Posterior Distributions of Performance Metrics
  # ------------------------------------------------------------------------------

  message("Creating performance metric distributions...")
  
  # Manuscript-ready: Compact layout with readable text
  pdf(result_path("performance_distributions.pdf"), width = 10, height = 7)
  
  # ============================================================================
  # Page 1: Fold-Specific AUC Distributions (5 folds only)
  # ============================================================================

  # Layout: 2 rows × 3 columns (5 folds, last slot empty)
  # Compact margins suitable for manuscript
  par(mfrow = c(2, 3),
      mar = c(4, 4, 2.5, 0.5),  # Tighter margins: bottom, left, top, right
      oma = c(0, 0, 3, 0))       # Outer margin for title only

  # Define common x-axis limits for AUC (for comparability)
  auc_xlim <- c(0.5, 1.0)

  # Plot each fold's AUC distribution
  for (k in 1:K) {
    fold_AUC_draws <- cv_AUC_by_draw[k, ]
    fold_AUC_draws <- fold_AUC_draws[!is.na(fold_AUC_draws)]

    fold_AUC_mean <- mean(fold_AUC_draws)
    fold_AUC_CrI <- quantile(fold_AUC_draws, c(0.025, 0.975))

    hist(fold_AUC_draws, breaks = 30,
         col = "lightblue", border = "darkblue",
         main = sprintf("Fold %d", k),
         xlab = "AUC", ylab = "Frequency",
         xlim = auc_xlim,
         cex.main = 1.1, cex.lab = 1.0, cex.axis = 0.9)
    abline(v = fold_AUC_mean, col = "red", lwd = 2, lty = 1)
    abline(v = fold_AUC_CrI, col = "red", lwd = 2, lty = 2)

    # Add summary text (left side for AUC)
    text(0.55, max(hist(fold_AUC_draws, breaks = 30, plot = FALSE)$counts) * 0.85,
         sprintf("%.3f\n[%.3f, %.3f]",
                 fold_AUC_mean, fold_AUC_CrI[1], fold_AUC_CrI[2]),
         cex = 0.85, adj = 0)
    grid(col = "gray85", lty = "dotted")
  }

  # Leave 6th slot empty (no Overall CV plot)
  plot.new()

  # Overall title for page 1
  mtext("Cross-Validation: Posterior Distributions of AUC by Fold",
        outer = TRUE, cex = 1.3, font = 2, line = 1)
  
  
  # ============================================================================
  # Page 2: Fold-Specific Brier Score Distributions (5 folds only)
  # ============================================================================

  # New page with same compact layout
  par(mfrow = c(2, 3),
      mar = c(4, 4, 2.5, 0.5),
      oma = c(0, 0, 3, 0))

  # Define common x-axis limits for Brier (for comparability)
  brier_xlim <- range(cv_Brier_all_draws)

  # Plot each fold's Brier distribution
  for (k in 1:K) {
    fold_Brier_draws <- cv_Brier_by_draw[k, ]
    fold_Brier_draws <- fold_Brier_draws[!is.na(fold_Brier_draws)]

    fold_Brier_mean <- mean(fold_Brier_draws)
    fold_Brier_CrI <- quantile(fold_Brier_draws, c(0.025, 0.975))

    hist(fold_Brier_draws, breaks = 30,
         col = "lightgreen", border = "darkgreen",
         main = sprintf("Fold %d", k),
         xlab = "Brier Score", ylab = "Frequency",
         xlim = brier_xlim,
         cex.main = 1.1, cex.lab = 1.0, cex.axis = 0.9)
    abline(v = fold_Brier_mean, col = "red", lwd = 2, lty = 1)
    abline(v = fold_Brier_CrI, col = "red", lwd = 2, lty = 2)

    # Add summary text (to the right of upper credible interval)
    text(fold_Brier_CrI[2] + diff(brier_xlim) * 0.02,
         max(hist(fold_Brier_draws, breaks = 30, plot = FALSE)$counts) * 0.85,
         sprintf("%.4f\n[%.4f, %.4f]",
                 fold_Brier_mean, fold_Brier_CrI[1], fold_Brier_CrI[2]),
         cex = 0.85, adj = 0)
    grid(col = "gray85", lty = "dotted")
  }

  # Leave 6th slot empty (no Overall CV plot)
  plot.new()

  # Overall title for page 2
  mtext("Cross-Validation: Posterior Distributions of Brier Score by Fold",
        outer = TRUE, cex = 1.3, font = 2, line = 1)
  
  
  # ============================================================================
  # Page 3: Final Validation Set Distributions (for comparison with CV)
  # ============================================================================

  # Compact layout for 2 plots
  par(mfrow = c(1, 2),
      mar = c(4.5, 4.5, 3, 1),
      oma = c(0, 0, 3, 0))

  # Define x-axis limits for Final Validation AUC (0.7 to 1.0 for better visualization)
  test_auc_xlim <- c(0.7, 1.0)

  # Final validation AUC distribution
  hist(test_AUC_by_draw, breaks = 50, col = "lightcoral", border = "darkred",
       main = "AUC",
       xlab = "AUC", ylab = "Frequency", xlim = test_auc_xlim,
       cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.0)
  abline(v = test_AUC_mean, col = "blue", lwd = 2.5, lty = 1)
  abline(v = test_AUC_CrI, col = "blue", lwd = 2, lty = 2)
  legend("topleft",
         legend = c("Mean", "95% CrI"),
         col = c("blue", "blue"), lty = c(1, 2), lwd = c(2.5, 2),
         cex = 0.95, bg = "white")
  text(0.72, max(hist(test_AUC_by_draw, breaks = 50, plot = FALSE)$counts) * 0.85,
       sprintf("%.3f\n[%.3f, %.3f]",
               test_AUC_mean, test_AUC_CrI[1], test_AUC_CrI[2]),
       cex = 0.95, adj = 0)
  grid(col = "gray85", lty = "dotted")

  # Final validation Brier distribution
  hist(test_Brier_by_draw, breaks = 50, col = "lightyellow", border = "orange",
       main = "Brier Score",
       xlab = "Brier Score", ylab = "Frequency",
       cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.0)
  abline(v = test_Brier_mean, col = "blue", lwd = 2.5, lty = 1)
  abline(v = test_Brier_CrI, col = "blue", lwd = 2, lty = 2)
  legend("topright",
         legend = c("Mean", "95% CrI"),
         col = c("blue", "blue"), lty = c(1, 2), lwd = c(2.5, 2),
         cex = 0.95, bg = "white")
  # Text to the left of upper credible interval (right-aligned to fit within margins)
  text(test_Brier_CrI[2] - diff(range(test_Brier_by_draw)) * 0.02,
       max(hist(test_Brier_by_draw, breaks = 50, plot = FALSE)$counts) * 0.85,
       sprintf("%.4f\n[%.4f, %.4f]",
               test_Brier_mean, test_Brier_CrI[1], test_Brier_CrI[2]),
       cex = 0.95, adj = 1)
  grid(col = "gray85", lty = "dotted")

  # Overall title for page 3
  mtext("Final Validation Set: Posterior Distributions of Performance Metrics",
        outer = TRUE, cex = 1.3, font = 2, line = 1)

  dev.off()
  
  
  # ------------------------------------------------------------------------------
  # 3. Posterior Inclusion Probabilities (PIPs)
  # ------------------------------------------------------------------------------
  
  message("Creating PIP visualizations...")
  
  pdf(result_path("pip_analysis.pdf"), width = 14, height = 10)
  par(mfrow = c(2, 2))
  
  # 3a. Top 30 genes by full-data PIP
  top_30_genes <- head(PIP_full_sorted, 30)
  barplot(top_30_genes,
          main = "Top 30 Genes by Posterior Inclusion Probability (Full-Data Model)",
          ylab = "PIP", xlab = "Gene",
          col = "steelblue", border = "navy",
          las = 2, cex.names = 0.6)
  abline(h = tau_full, col = "red", lwd = 2, lty = 2)
  text(5, tau_full + 0.05, sprintf("FDR threshold (τ = %.3f)", tau_full),
       col = "red", cex = 0.8)
  grid()
  
  # 3b. PIP distribution across all union genes
  hist(PIP_full, breaks = 100, col = "lightblue", border = "blue",
       main = "Distribution of PIPs Across Union Genes (Full-Data Model)",
       xlab = "Posterior Inclusion Probability", ylab = "Frequency")
  abline(v = tau_full, col = "red", lwd = 2, lty = 2)
  text(tau_full + 0.05, max(hist(PIP_full, breaks = 100, plot = FALSE)$counts) * 0.9,
       sprintf("FDR threshold\n(τ = %.3f)", tau_full), col = "red", cex = 0.8)
  grid()
  
  # 3c. FDR curve (full-data model)
  plot(1:min(1000, length(FDR_hat_full)), FDR_hat_full[1:min(1000, length(FDR_hat_full))],
       type = "l", lwd = 2, col = "darkblue",
       main = "Expected FDR vs. Number of Selected Genes (Full-Data Model)",
       xlab = "Number of Selected Genes (K)", ylab = "Expected FDR")
  abline(h = 0.05, col = "red", lwd = 2, lty = 2)  # alpha_FDR = 0.05 for full model
  abline(v = K_star_full, col = "green", lwd = 2, lty = 2)
  legend("topleft",
         legend = c("FDR = 0.05",
                    sprintf("Optimal K = %d", K_star_full)),
         col = c("red", "green"), lty = 2, lwd = 2)
  grid()
  
  # 3d. Fold consistency of PIPs (top 50 genes from union)
  # Show how PIPs varied across folds for top genes from full model
  top_50_genes_names <- names(head(PIP_full_sorted, min(50, length(PIP_full_sorted))))
  pip_top50_by_fold <- cv_PIPs_by_fold[, top_50_genes_names, drop = FALSE]
  
  boxplot(t(pip_top50_by_fold),
          main = "PIP Consistency Across Folds (Top 50 Genes)",
          xlab = "Gene", ylab = "PIP",
          las = 2, cex.axis = 0.5, col = "lightgreen", border = "darkgreen",
          outline = FALSE)
  grid()
  
  dev.off()
  
  
  # ------------------------------------------------------------------------------
  # 4. Clinical Coefficients (β) Posterior Distributions
  # ------------------------------------------------------------------------------
  
  message("Creating beta coefficient visualizations...")
  
  # Use actual number of betas from beta_summary (already computed above)
  n_betas <- nrow(beta_summary)
  message(sprintf("  Plotting posterior distributions for %d clinical covariates", n_betas))
  
  # Dynamically determine layout based on number of coefficients
  if (n_betas <= 5) {
    n_rows <- 1
    n_cols <- n_betas
    pdf_width <- 3.5 * n_betas
    pdf_height <- 5
  } else if (n_betas <= 10) {
    n_rows <- 2
    n_cols <- ceiling(n_betas / 2)
    pdf_width <- 3.5 * n_cols
    pdf_height <- 8
  } else {
    n_rows <- 3
    n_cols <- ceiling(n_betas / 3)
    pdf_width <- 3.5 * n_cols
    pdf_height <- 11
  }
  
  # Single-page PDF suitable for manuscript inclusion
  pdf(result_path("beta_coefficients.pdf"), width = pdf_width, height = pdf_height)
  par(mfrow = c(n_rows, n_cols),
      mar = c(4.5, 4.5, 3, 1),   # Margins: bottom, left, top, right
      oma = c(0, 0, 2, 0))        # Outer margins for overall title
  
  for (j in 1:n_betas) {
    covariate_name <- beta_summary$Covariate[j]
    beta_draws_j <- bart_fit_full$beta_draws[, j]
    
    hist(beta_draws_j, breaks = 50, col = "lavender", border = "purple",
         main = sprintf("%s", covariate_name),
         xlab = expression(beta), ylab = "Frequency",
         cex.main = 1.1, cex.lab = 1.0, cex.axis = 0.9)
    
    # Add posterior mean and credible interval
    abline(v = beta_posterior_mean[j], col = "red", lwd = 2.5, lty = 1)
    abline(v = beta_credible_intervals[, j], col = "red", lwd = 2, lty = 2)
    abline(v = 0, col = "black", lwd = 1.5, lty = 3)
    
    legend("topright",
           legend = c("Mean", "95% CrI", "Zero"),
           col = c("red", "red", "black"),
           lty = c(1, 2, 3), lwd = c(2.5, 2, 1.5),
           cex = 0.7, bg = "white")
    grid(col = "gray80", lty = "dotted")
  }
  
  # Add overall title
  mtext("Posterior Distributions of Clinical Coefficients",
        outer = TRUE, cex = 1.5, font = 2, line = 0.5)
  
  dev.off()
  message(sprintf("  All %d coefficient plots on single page (manuscript-ready)", n_betas))


  # ------------------------------------------------------------------------------
  # 4b. Posterior Probability of Positive Effect (P(β > 0)) for Clinical Covariates
  # ------------------------------------------------------------------------------

  message("Computing posterior probability of positive effect for clinical covariates...")

  # Compute P(β > 0) for each clinical covariate
  posterior_prob_positive <- numeric(n_betas)
  for (j in 1:n_betas) {
    beta_draws_j <- bart_fit_full$beta_draws[, j]
    posterior_prob_positive[j] <- mean(beta_draws_j > 0)
  }

  # Create summary data frame
  beta_direction_summary <- data.frame(
    Covariate = beta_summary$Covariate,
    Mean = beta_summary$Mean,
    SD = beta_summary$SD,
    CI_Lower = beta_summary$CI_Lower,
    CI_Upper = beta_summary$CI_Upper,
    Prob_Positive = posterior_prob_positive,
    Prob_Negative = 1 - posterior_prob_positive
  )

  message("Posterior probabilities of positive effect:")
  for (j in 1:n_betas) {
    message(sprintf("  %s: P(β > 0) = %.4f", beta_direction_summary$Covariate[j],
                    beta_direction_summary$Prob_Positive[j]))
  }

  # --- PDF Output: Posterior Probability Table ---
  message("Creating posterior probability PDF...")

  pdf(result_path("beta_posterior_probability.pdf"), width = 10, height = 6)

  # Set up blank plot for table
  par(mar = c(1, 1, 3, 1))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  # Title
  title(main = "Posterior Probability of Positive Effect for Clinical Covariates",
        cex.main = 1.4, font.main = 2)

  # Table header
  header_y <- 0.85
  col_positions <- c(0.05, 0.25, 0.40, 0.55, 0.75, 0.90)
  headers <- c("Covariate", "Mean", "95% CrI", "P(β > 0)", "P(β < 0)", "Direction")

  text(col_positions, rep(header_y, 6), headers, font = 2, adj = 0, cex = 1.0)
  segments(0.02, header_y - 0.03, 0.98, header_y - 0.03, lwd = 1.5)

  # Table rows
  row_height <- 0.10
  for (j in 1:n_betas) {
    row_y <- header_y - 0.08 - (j - 1) * row_height

    # Covariate name
    text(col_positions[1], row_y, beta_direction_summary$Covariate[j], adj = 0, cex = 0.95)

    # Mean
    text(col_positions[2], row_y, sprintf("%.4f", beta_direction_summary$Mean[j]), adj = 0, cex = 0.95)

    # 95% CrI
    text(col_positions[3], row_y,
         sprintf("[%.3f, %.3f]", beta_direction_summary$CI_Lower[j], beta_direction_summary$CI_Upper[j]),
         adj = 0, cex = 0.95)

    # P(β > 0)
    prob_pos <- beta_direction_summary$Prob_Positive[j]
    text(col_positions[4], row_y, sprintf("%.4f", prob_pos), adj = 0, cex = 0.95)

    # P(β < 0)
    text(col_positions[5], row_y, sprintf("%.4f", beta_direction_summary$Prob_Negative[j]), adj = 0, cex = 0.95)

    # Direction indicator
    if (prob_pos > 0.95) {
      direction <- "Strong +"
      dir_col <- "darkgreen"
    } else if (prob_pos > 0.80) {
      direction <- "Moderate +"
      dir_col <- "green3"
    } else if (prob_pos < 0.05) {
      direction <- "Strong -"
      dir_col <- "darkred"
    } else if (prob_pos < 0.20) {
      direction <- "Moderate -"
      dir_col <- "red3"
    } else {
      direction <- "Uncertain"
      dir_col <- "gray50"
    }
    text(col_positions[6], row_y, direction, adj = 0, cex = 0.95, col = dir_col, font = 2)
  }

  # Add footnote
  text(0.5, 0.08,
       "Note: P(β > 0) is the posterior probability that the coefficient is positive.",
       adj = 0.5, cex = 0.85, font = 3)
  text(0.5, 0.03,
       "Direction: Strong (+/-) if P > 0.95 or P < 0.05; Moderate if P > 0.80 or P < 0.20; Uncertain otherwise.",
       adj = 0.5, cex = 0.85, font = 3)

  dev.off()
  message("  Saved: beta_posterior_probability.pdf")

  # --- TeX Output: Posterior Probability Table ---
  message("Creating posterior probability TeX table...")

  tex_prob_file <- result_path("beta_posterior_probability.tex")

  cat("\\begin{table}[ht]\n", file = tex_prob_file)
  cat("\\centering\n", file = tex_prob_file, append = TRUE)
  cat("\\caption{Posterior Probability of Positive Effect for Clinical Covariates}\n",
      file = tex_prob_file, append = TRUE)
  cat("\\label{tab:beta_posterior_probability}\n", file = tex_prob_file, append = TRUE)
  cat("\\begin{tabular}{l c c c c}\n", file = tex_prob_file, append = TRUE)
  cat("\\hline\n", file = tex_prob_file, append = TRUE)
  cat("\\textbf{Covariate} & \\textbf{Mean} & \\textbf{95\\% CrI} & \\textbf{P($\\beta > 0$)} & \\textbf{P($\\beta < 0$)} \\\\\n",
      file = tex_prob_file, append = TRUE)
  cat("\\hline\n", file = tex_prob_file, append = TRUE)

  for (j in 1:n_betas) {
    covariate_name <- gsub("_", "\\\\_", beta_direction_summary$Covariate[j])
    cat(sprintf("%s & %.4f & [%.4f, %.4f] & %.4f & %.4f \\\\\n",
                covariate_name,
                beta_direction_summary$Mean[j],
                beta_direction_summary$CI_Lower[j],
                beta_direction_summary$CI_Upper[j],
                beta_direction_summary$Prob_Positive[j],
                beta_direction_summary$Prob_Negative[j]),
        file = tex_prob_file, append = TRUE)
  }

  cat("\\hline\n", file = tex_prob_file, append = TRUE)
  cat("\\end{tabular}\n", file = tex_prob_file, append = TRUE)
  cat("\\vspace{2mm}\n", file = tex_prob_file, append = TRUE)
  cat("\\parbox{0.9\\textwidth}{\\footnotesize \\textit{Note:} P($\\beta > 0$) is the posterior probability that the coefficient is positive, computed as the proportion of MCMC draws where $\\beta > 0$.}\n",
      file = tex_prob_file, append = TRUE)
  cat("\\end{table}\n", file = tex_prob_file, append = TRUE)

  message("  Saved: beta_posterior_probability.tex")

  # Save the summary data frame as RDS for future use
  saveRDS(beta_direction_summary, result_path("beta_posterior_probability.RDS"))
  message("  Saved: beta_posterior_probability.RDS")


  # ------------------------------------------------------------------------------
  # 5. Calibration Plots
  # ------------------------------------------------------------------------------
  
  message("Creating calibration plots...")
  
  pdf(result_path("calibration.pdf"), width = 12, height = 6)
  par(mfrow = c(1, 2))
  
  # CV calibration
  cv_bins <- cut(cv_predictions_mean, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
  cv_calibration <- aggregate(Y_cv, by = list(cv_bins), FUN = mean)
  cv_predicted <- aggregate(cv_predictions_mean, by = list(cv_bins), FUN = mean)
  
  plot(cv_predicted$x, cv_calibration$x,
       xlim = c(0, 1), ylim = c(0, 1),
       pch = 19, cex = 2, col = "blue",
       main = "Calibration Plot: Cross-Validation",
       xlab = "Predicted Probability", ylab = "Observed Frequency")
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  grid()
  text(0.2, 0.8, "Perfect calibration", col = "red", srt = 45, cex = 0.9)
  
  # Test calibration
  test_bins <- cut(test_predictions_mean, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
  test_calibration <- aggregate(Y_test, by = list(test_bins), FUN = mean)
  test_predicted <- aggregate(test_predictions_mean, by = list(test_bins), FUN = mean)
  
  plot(test_predicted$x, test_calibration$x,
       xlim = c(0, 1), ylim = c(0, 1),
       pch = 19, cex = 2, col = "red",
       main = "Calibration Plot: Test Set",
       xlab = "Predicted Probability", ylab = "Observed Frequency")
  abline(0, 1, col = "blue", lwd = 2, lty = 2)
  grid()
  text(0.2, 0.8, "Perfect calibration", col = "blue", srt = 45, cex = 0.9)
  
  dev.off()
  
  
  # ------------------------------------------------------------------------------
  # 6. Prediction Distribution by Outcome
  # ------------------------------------------------------------------------------
  
  message("Creating prediction distribution plots...")
  
  pdf(result_path("prediction_distributions.pdf"), width = 12, height = 6)
  par(mfrow = c(1, 2))
  
  # CV predictions by outcome
  boxplot(cv_predictions_mean ~ Y_cv,
          col = c("lightblue", "lightcoral"),
          main = "CV: Predicted Probabilities by Outcome",
          xlab = "True Outcome", ylab = "Predicted Probability",
          names = c("Control (Y=0)", "Case (Y=1)"))
  grid()
  
  # Test predictions by outcome
  boxplot(test_predictions_mean ~ Y_test,
          col = c("lightblue", "lightcoral"),
          main = "Test: Predicted Probabilities by Outcome",
          xlab = "True Outcome", ylab = "Predicted Probability",
          names = c("Control (Y=0)", "Case (Y=1)"))
  grid()
  
  dev.off()
  
  
  # ------------------------------------------------------------------------------
  # 7. Summary Table (LaTeX format)
  # ------------------------------------------------------------------------------
  
  message("Creating LaTeX summary tables...")
  
  # Performance summary table
  perf_table_file <- result_path("performance_summary.tex")
  cat("\\begin{table}[ht]\n", file = perf_table_file)
  cat("\\centering\n", file = perf_table_file, append = TRUE)
  cat("\\caption{Semi-Parametric Probit DART Model Performance}\n",
      file = perf_table_file, append = TRUE)
  cat("\\label{tab:performance_summary}\n",
      file = perf_table_file, append = TRUE)
  cat("\\begin{tabular}{l c c c}\n",
      file = perf_table_file, append = TRUE)
  cat("\\hline\n", file = perf_table_file, append = TRUE)
  cat("\\textbf{Dataset} & \\textbf{Metric} & \\textbf{Mean} & \\textbf{95\\% CrI} \\\\\n",
      file = perf_table_file, append = TRUE)
  cat("\\hline\n", file = perf_table_file, append = TRUE)
  
  cat(sprintf("Cross-Validation & AUC & %.4f & [%.4f, %.4f] \\\\\n",
              cv_AUC_mean, cv_AUC_CrI[1], cv_AUC_CrI[2]),
      file = perf_table_file, append = TRUE)
  cat(sprintf("Cross-Validation & Brier & %.4f & [%.4f, %.4f] \\\\\n",
              cv_Brier_mean, cv_Brier_CrI[1], cv_Brier_CrI[2]),
      file = perf_table_file, append = TRUE)
  cat("\\hline\n", file = perf_table_file, append = TRUE)
  cat(sprintf("Test Set & AUC & %.4f & [%.4f, %.4f] \\\\\n",
              test_AUC_mean, test_AUC_CrI[1], test_AUC_CrI[2]),
      file = perf_table_file, append = TRUE)
  cat(sprintf("Test Set & Brier & %.4f & [%.4f, %.4f] \\\\\n",
              test_Brier_mean, test_Brier_CrI[1], test_Brier_CrI[2]),
      file = perf_table_file, append = TRUE)
  
  cat("\\hline\n", file = perf_table_file, append = TRUE)
  cat("\\end{tabular}\n", file = perf_table_file, append = TRUE)
  cat("\\end{table}\n", file = perf_table_file, append = TRUE)
  
  
  # Clinical coefficients table
  beta_table_file <- result_path("beta_coefficients.tex")
  cat("\\begin{table}[ht]\n", file = beta_table_file)
  cat("\\centering\n", file = beta_table_file, append = TRUE)
  cat("\\caption{Posterior Estimates of Clinical Coefficients}\n",
      file = beta_table_file, append = TRUE)
  cat("\\label{tab:beta_coefficients}\n",
      file = beta_table_file, append = TRUE)
  cat("\\begin{tabular}{l c c c}\n",
      file = beta_table_file, append = TRUE)
  cat("\\hline\n", file = beta_table_file, append = TRUE)
  cat("\\textbf{Covariate} & \\textbf{Mean} & \\textbf{SD} & \\textbf{95\\% CrI} \\\\\n",
      file = beta_table_file, append = TRUE)
  cat("\\hline\n", file = beta_table_file, append = TRUE)
  
  for (j in 1:nrow(beta_summary)) {
    covariate_name <- gsub("_", "\\\\_", beta_summary$Covariate[j])
    cat(sprintf("%s & %.4f & %.4f & [%.4f, %.4f] \\\\\n",
                covariate_name,
                beta_summary$Mean[j],
                beta_summary$SD[j],
                beta_summary$CI_Lower[j],
                beta_summary$CI_Upper[j]),
        file = beta_table_file, append = TRUE)
  }
  
  cat("\\hline\n", file = beta_table_file, append = TRUE)
  cat("\\end{tabular}\n", file = beta_table_file, append = TRUE)
  cat("\\end{table}\n", file = beta_table_file, append = TRUE)
  
  
  # Top selected genes table
  if (exists("selected_genes_final") && length(selected_genes_final) > 0) {
    top_20_selected <- head(selected_genes_final, 20)
    # Get PIPs from full model
    top_20_pips <- PIP_full_sorted[top_20_selected]
    
    genes_table_file <- result_path("top_selected_genes.tex")
    cat("\\begin{table}[ht]\n", file = genes_table_file)
    cat("\\centering\n", file = genes_table_file, append = TRUE)
    cat(sprintf("\\caption{Top 20 Selected Genes (FDR = %.2f, K = %d)}\n",
                0.05, length(selected_genes_final)),
        file = genes_table_file, append = TRUE)
    cat("\\label{tab:top_selected_genes}\n",
        file = genes_table_file, append = TRUE)
    cat("\\begin{tabular}{r l c}\n",
        file = genes_table_file, append = TRUE)
    cat("\\hline\n", file = genes_table_file, append = TRUE)
    cat("\\textbf{Rank} & \\textbf{Gene} & \\textbf{PIP} \\\\\n",
        file = genes_table_file, append = TRUE)
    cat("\\hline\n", file = genes_table_file, append = TRUE)
    
    for (i in 1:length(top_20_selected)) {
      gene_name <- gsub("_", "\\\\_", top_20_selected[i])
      cat(sprintf("%d & %s & %.4f \\\\\n", i, gene_name, top_20_pips[i]),
          file = genes_table_file, append = TRUE)
    }
    
    cat("\\hline\n", file = genes_table_file, append = TRUE)
    cat("\\end{tabular}\n", file = genes_table_file, append = TRUE)
    cat("\\end{table}\n", file = genes_table_file, append = TRUE)
  }
  
  # ==============================================================================
  # Final Prognostic Model: Complete List of Selected Genes (PDF and TeX)
  # ==============================================================================
  if (exists("selected_genes_final") && length(selected_genes_final) > 0) {
    
    message("\n=== Generating Final Prognostic Model Gene List ===")
    
    n_selected <- length(selected_genes_final)
    message(sprintf("  Total genes selected: %d", n_selected))
    
    # --- PDF Output: Simple list of gene names ---
    pdf(result_path("final_prognostic_genes.pdf"), width = 8.5, height = 11)
    
    # Calculate layout: genes per page
    genes_per_col <- 50
    cols_per_page <- 3
    genes_per_page <- genes_per_col * cols_per_page
    n_pages <- ceiling(n_selected / genes_per_page)
    
    for (page in 1:n_pages) {
      # Set up margins for clean text layout
      par(mar = c(2, 2, 3, 2))
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1))
      
      # Title
      title(main = sprintf("Final Prognostic Model: Selected Genes (Page %d of %d)", page, n_pages),
            cex.main = 1.2, font.main = 2)
      mtext(sprintf("Total: %d genes | FDR α = 0.05", n_selected), side = 3, line = 0, cex = 0.9)
      
      # Get genes for this page
      start_idx <- (page - 1) * genes_per_page + 1
      end_idx <- min(page * genes_per_page, n_selected)
      genes_this_page <- selected_genes_final[start_idx:end_idx]
      
      # Print genes in columns
      col_width <- 1 / cols_per_page
      line_height <- 0.95 / genes_per_col
      
      for (i in seq_along(genes_this_page)) {
        col <- ((i - 1) %% cols_per_page) + 1
        row <- ((i - 1) %/% cols_per_page) + 1
        
        if (row <= genes_per_col) {
          x_pos <- (col - 1) * col_width + 0.02
          y_pos <- 0.95 - (row - 1) * line_height
          
          # Gene number and name
          gene_label <- sprintf("%d. %s", start_idx + i - 1, genes_this_page[i])
          text(x_pos, y_pos, gene_label, adj = c(0, 0.5), cex = 0.6, family = "mono")
        }
      }
    }
    
    dev.off()
    message("  Saved: final_prognostic_genes.pdf")
    
    # --- TeX Output: Simple list of gene names ---
    tex_file <- result_path("final_prognostic_genes.tex")
    
    cat("\\begin{table}[ht]\n", file = tex_file)
    cat("\\centering\n", file = tex_file, append = TRUE)
    cat(sprintf("\\caption{Complete List of Selected Genes from Final Prognostic Model (n = %d, FDR $\\alpha$ = 0.05)}\n", n_selected),
        file = tex_file, append = TRUE)
    cat("\\label{tab:final_prognostic_genes}\n", file = tex_file, append = TRUE)
    cat("\\begin{tabular}{r l}\n", file = tex_file, append = TRUE)
    cat("\\hline\n", file = tex_file, append = TRUE)
    cat("\\textbf{Rank} & \\textbf{Gene Name} \\\\\n", file = tex_file, append = TRUE)
    cat("\\hline\n", file = tex_file, append = TRUE)
    
    for (i in seq_along(selected_genes_final)) {
      gene_name <- gsub("_", "\\\\_", selected_genes_final[i])
      cat(sprintf("%d & %s \\\\\n", i, gene_name), file = tex_file, append = TRUE)
    }
    
    cat("\\hline\n", file = tex_file, append = TRUE)
    cat("\\end{tabular}\n", file = tex_file, append = TRUE)
    cat("\\end{table}\n", file = tex_file, append = TRUE)
    
    message("  Saved: final_prognostic_genes.tex")

    # --- RDS Output: Simple vector of gene names only ---
    # Save just the gene names as a character vector (no ranks, no PIPs)
    saveRDS(selected_genes_final, result_path("selected_gene_names.RDS"))
    message("  Saved: selected_gene_names.RDS")

    # --- PDF Output: Simple list of gene names only (no ranks) ---
    message("Creating simple gene names list PDF...")

    pdf(result_path("selected_gene_names.pdf"), width = 8.5, height = 11)

    # Calculate layout: genes per page
    genes_per_col <- 45
    cols_per_page <- 3
    genes_per_page <- genes_per_col * cols_per_page
    n_pages_simple <- ceiling(n_selected / genes_per_page)

    for (page in 1:n_pages_simple) {
      # Set up margins for clean text layout
      par(mar = c(2, 2, 4, 2))
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1))

      # Title
      title(main = "Selected Genes from Final Prognostic Model",
            cex.main = 1.4, font.main = 2)
      mtext(sprintf("Total: %d genes | FDR α = 0.05", n_selected), side = 3, line = 0.5, cex = 1.0)
      if (n_pages_simple > 1) {
        mtext(sprintf("(Page %d of %d)", page, n_pages_simple), side = 3, line = -0.5, cex = 0.9)
      }

      # Get genes for this page
      start_idx <- (page - 1) * genes_per_page + 1
      end_idx <- min(page * genes_per_page, n_selected)
      genes_this_page <- selected_genes_final[start_idx:end_idx]

      # Print genes in columns (names only, no numbers)
      col_width <- 1 / cols_per_page
      line_height <- 0.92 / genes_per_col

      for (i in seq_along(genes_this_page)) {
        col <- ((i - 1) %% cols_per_page) + 1
        row <- ((i - 1) %/% cols_per_page) + 1

        if (row <= genes_per_col) {
          x_pos <- (col - 1) * col_width + 0.02
          y_pos <- 0.92 - (row - 1) * line_height

          # Gene name only (no number)
          text(x_pos, y_pos, genes_this_page[i], adj = c(0, 0.5), cex = 0.65, family = "mono")
        }
      }
    }

    dev.off()
    message("  Saved: selected_gene_names.pdf")

  } else {
    message("  Skipping final prognostic gene list - no genes selected")
  }

  message("\n=== Visualization Complete! ===")
  message("All files saved to: ", results_dir)
  message("  - roc_curves.pdf")
  message("  - performance_metrics_grid.pdf")
  message("  - performance_distributions.pdf")
  message("  - pip_analysis.pdf")
  message("  - beta_coefficients.pdf")
  message("  - beta_posterior_probability.pdf")
  message("  - beta_posterior_probability.tex")
  message("  - calibration.pdf")
  message("  - prediction_distributions.pdf")
  message("  - performance_summary.tex")
  message("  - beta_coefficients.tex")
  message("  - top_selected_genes.tex")
  message("  - final_prognostic_genes.pdf")
  message("  - final_prognostic_genes.tex")
  message("  - selected_gene_names.pdf")
  message("  - selected_gene_names.RDS")
  
} else {
  warning("Visualizations skipped - missing required variables.")
  message("This can happen if:")
  message("  1. The full model was not fitted successfully")
  message("  2. Test set evaluation was skipped")
  message("  3. CV predictions were not generated")
  message("Check earlier sections for errors or warnings.")
}

























