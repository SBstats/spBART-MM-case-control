# ------------------------------------------------------------------------------
# Naive Analysis: Elastic Net with Separate Penalties for Clinical vs. Genes
# ------------------------------------------------------------------------------
#
# IMPORTANT MODIFICATION (Added):
# This script now uses penalty.factor in glmnet to handle clinical covariates
# and gene expression data differently:
#
# TWO-STAGE COVARIATE HANDLING:
#
# STAGE 1 - Screening & Elastic Net Variable Selection:
# - Clinical covariates (age, sex_binary ONLY): penalty.factor = 0 (UNPENALIZED)
#   * Minimal adjustment set for gene screening
#   * Prevents over-adjustment during variable selection
#
# STAGE 2 - Final Prognostic Model:
# - Clinical covariates (age, sex_binary, race_binary, bmi_binary): penalty.factor = 0 (UNPENALIZED)
#   * Full covariate adjustment in final model
#   * Ensures stable, interpretable coefficients
#   * Maintains biologically expected directions (e.g., age/BMI positive)
#
# - Study indicator: penalty.factor = 1 (STANDARD PENALTY, treated like genes)
#   * Moved to X vector (penalized feature set)
#   * Subject to elastic net regularization for variable selection
#   * Can be shrunk to zero if not predictive after accounting for other features
#
# - Gene expression variables: penalty.factor = 1 (STANDARD PENALTY)
#   * Subject to elastic net regularization for variable selection
#   * Shrinkage applied to control overfitting
#
# This approach addresses the issue where clinical coefficients had unexpected
# negative signs in standard elastic net, likely due to high-dimensional
# confounding from correlated gene expression data.
#
# NOTE: study_indicator is now in the X (penalized) vector, NOT the Z (unpenalized)
# clinical covariates vector. This allows the model to determine if the study
# indicator provides additional predictive information beyond the genes and
# clinical covariates.
#
# ------------------------------------------------------------------------------
# DATA PREPROCESSING
#





########################################################
# 1) Load required libraries
########################################################
# ------------------------------------------------------------------------------
# Requirements:
#   install.packages(c("BART","pROC","future.apply","matrixStats","MASS","tidyverse"))
# ------------------------------------------------------------------------------
library(tidyverse)
library(BART)            # pbart/wbart 
library(pROC)
library(MASS)            # mvrnorm
library(future.apply)    # parallel
library(matrixStats)
library(tidyverse)
library(data.table)
library(gtools)
#library(openxlsx)
library(dplyr)
library(readxl)

# if (!require("xlsx", character.only = TRUE)) {
#   install.packages("xlsx", dependencies = TRUE)
#   library(xlsx)
# } else {
#   library(xlsx)
# }


if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!require("DESeq2", quietly = TRUE)){
  BiocManager::install("DESeq2")
}

library(DESeq2)

if (!require("glmnet")) install.packages("glmnet")
library(glmnet)

if (!require("caret")) install.packages("caret")
library(caret)

# For Gaussian Mixture Model clustering of coefficients
if (!require("mclust")) install.packages("mclust")
library(mclust)

# Note: MASS package not needed since we're using logistic regression via glm()


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

# 3) code bmi to binary non-obese (0) vs obese (1)
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

# 2) code bmi to binary non-obese (0) vs obese (1)
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

cohort_selected$study_indicator = 0


# Select variables from caseControl_workdf
caseControl_selected <- caseControl_workdf %>%
  dplyr::select(PooledData_ID, race_binary, bmi_binary, sex_binary, age, mm_status)

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














#########################################################################################################################
# 3) Elastic Net Regularization with Two-Stage Logistic Regression Screening
#########################################################################################################################
# Methodology:
# - Split data: 500 model development set + 369 final validation set
# - Two-stage screening:
#   * Stage 1: Logistic regression on development set (p < 0.05) adjusting for covariates
#   * Stage 2: One dim clustering of selected candidate genes' covariates and select the more significant ones
# - Fit elastic net 100 times on development set with screened genes
#   * Each fit: 10-fold CV within glmnet to select optimal lambda
#   * Use alpha=0.95 (nearly LASSO) for more aggressive L1 penalty
#   * Select variables at lambda.1se (more conservative) for each replication
#   * Stability selection: Keep genes selected in ≥50 replications
# - Evaluate final model on held-out validation set
#########################################################################################################################

# Set seed for reproducibility
# NOTE: Using seed 789 and TEST-FIRST sampling to match comparative analysis file
#       (04_comparative_analysis_ElasticNet_vs_spBART_continuousBMI_studyIndicator_in_BART_Nov28.R)
set.seed(789)

# Step 1: Create stratified split (500 development, 369 validation)
# IMPORTANT: Sample TEST/VALIDATION set FIRST to match comparative analysis methodology
message("\n=== Step 1: Data Partitioning ===")
message("Using seed=789 and TEST-FIRST sampling (consistent with comparative analysis)")
outcome_vec <- pooled_metadata$outcome
n_total <- nrow(pooled_metadata)
N_CV <- 500  # Development/CV set size
N_test <- n_total - N_CV  # Test/validation set size (369)

# Stratified sampling: sample TEST set FIRST (same as comparative analysis)
outcome_1_ids <- pooled_metadata$PooledData_ID[outcome_vec == 1]
outcome_0_ids <- pooled_metadata$PooledData_ID[outcome_vec == 0]

prop_cases <- sum(outcome_vec == 1) / n_total
n_test_cases <- round(N_test * prop_cases)
n_test_controls <- N_test - n_test_cases

# Sample test/validation set first
test_ids_cases <- sample(outcome_1_ids, n_test_cases, replace = FALSE)
test_ids_controls <- sample(outcome_0_ids, n_test_controls, replace = FALSE)
val_ids <- c(test_ids_cases, test_ids_controls)

# Development set: remaining IDs (after test set is sampled)
dev_ids <- setdiff(pooled_metadata$PooledData_ID, val_ids)

message(sprintf("Model development set: %d patients (%d cases, %d controls)",
                length(dev_ids),
                sum(pooled_metadata$outcome[pooled_metadata$PooledData_ID %in% dev_ids] == 1),
                sum(pooled_metadata$outcome[pooled_metadata$PooledData_ID %in% dev_ids] == 0)))
message(sprintf("Final validation set: %d patients (%d cases, %d controls)",
                length(val_ids),
                sum(pooled_metadata$outcome[pooled_metadata$PooledData_ID %in% val_ids] == 1),
                sum(pooled_metadata$outcome[pooled_metadata$PooledData_ID %in% val_ids] == 0)))

# Prepare development set data
dev_indices <- which(pooled_metadata$PooledData_ID %in% dev_ids)
X_dev_genes <- t(pooled_genebody_data_filtered_normalized[, dev_indices])  # 500 x 14147
y_dev <- pooled_metadata$outcome[dev_indices]

# Clinical covariates setup:
# - Z_dev_screening: age, sex_binary ONLY (used for initial screening and elastic net)
# - Z_dev_full: age, sex_binary, race_binary, bmi_binary (used ONLY for final prognostic model)
# NOTE: study_indicator is treated as part of the gene/feature set (penalized) rather than clinical covariates (unpenalized)
Z_dev_full <- as.matrix(pooled_metadata[dev_indices, c("age", "sex_binary", "race_binary", "bmi_binary")])
Z_dev_screening <- as.matrix(pooled_metadata[dev_indices, c("age", "sex_binary")])  # Only age and sex for screening/elastic net

# Extract study_indicator separately to add to X matrix
study_indicator_dev <- pooled_metadata$study_indicator[dev_indices]

# Prepare final validation set data
val_indices <- which(pooled_metadata$PooledData_ID %in% val_ids)
X_val_genes <- t(pooled_genebody_data_filtered_normalized[, val_indices])  # 369 x 14147
y_val <- pooled_metadata$outcome[val_indices]
Z_val_full <- as.matrix(pooled_metadata[val_indices, c("age", "sex_binary", "race_binary", "bmi_binary")])
Z_val_screening <- as.matrix(pooled_metadata[val_indices, c("age", "sex_binary")])  # Only age and sex for screening/elastic net

# Extract study_indicator separately to add to X matrix
study_indicator_val <- pooled_metadata$study_indicator[val_indices]

# Ensure both matrices have the same column names (gene IDs)
colnames(X_dev_genes) <- rownames(pooled_genebody_data_filtered_normalized)
colnames(X_val_genes) <- rownames(pooled_genebody_data_filtered_normalized)

# Verify column names match
message(sprintf("X_dev_genes columns: %d", ncol(X_dev_genes)))
message(sprintf("X_val_genes columns: %d", ncol(X_val_genes)))
message(sprintf("Column names match: %s", identical(colnames(X_dev_genes), colnames(X_val_genes))))











# Step 2: Initial Logistic Regression Screening on Development Set
message("\n=== Step 2: Initial Logistic Regression Screening on Development Set ===")

# Helper function: logistic regression p-value
logistic_pvalue <- function(y, x, Z) {
  # Fit logistic model: y ~ x + Z (gene + covariates)
  df <- data.frame(y = y, gene = x, Z)
  tryCatch({
    fit <- glm(y ~ ., data = df, family = binomial(link = "logit"))
    summary(fit)$coefficients["gene", "Pr(>|z|)"]
  }, error = function(e) return(1))  # Return non-significant if error
}

# Stage 1: Logistic regression screening (p < 0.05 AND clustering) on development set
message("Stage 1: Development set screening (p < 0.05 AND clustering)...")

# Compute p-values for all genes (using all clinical covariates: age, sex_binary, race_binary, bmi_binary, study_indicator)
pvals <- apply(X_dev_genes, 2, function(x) logistic_pvalue(y_dev, x, Z_dev_screening))
names(pvals) <- colnames(X_dev_genes)  # Ensure names are set

# Compute coefficients for all genes
dev_coefs <- sapply(colnames(X_dev_genes), function(gene) {
  df <- data.frame(y = y_dev, gene = X_dev_genes[, gene], Z_dev_screening)
  tryCatch({
    fit <- glm(y ~ ., data = df, family = binomial(link = "logit"))
    coef(fit)["gene"]
  }, error = function(e) return(0))
})
names(dev_coefs) <- colnames(X_dev_genes)  # Explicitly set names

# Filter 1: p < 0.05
significant_genes_dev <- names(pvals)[pvals < 0.05]
message(sprintf("  Genes with p < 0.05: %d out of %d",
                length(significant_genes_dev), ncol(X_dev_genes)))

# Filter 2: Gaussian Mixture Model (GMM) clustering on |coefficients|
# Assumption: Coefficients arise from mixture of effect distributions (null, weak, strong)
if (length(significant_genes_dev) > 0) {
  dev_coefs_significant <- dev_coefs[significant_genes_dev]
  # Remove NAs if any
  dev_coefs_significant <- dev_coefs_significant[!is.na(dev_coefs_significant)]

  if (length(dev_coefs_significant) >= 10) {  # Need sufficient genes for GMM
    message(sprintf("  Fitting Gaussian Mixture Model to %d coefficients...",
                    length(dev_coefs_significant)))

    # Use absolute values of coefficients for clustering
    abs_coefs_dev <- abs(dev_coefs_significant)

    # Fit GMM with 2-4 components (let mclust select optimal number via BIC)
    gmm_dev <- Mclust(abs_coefs_dev, G = 2:4, modelNames = "V", verbose = FALSE)

    if (!is.null(gmm_dev)) {
      # Get cluster assignments
      cluster_assignments_dev <- gmm_dev$classification
      n_clusters_dev <- gmm_dev$G

      message(sprintf("  Optimal number of clusters: %d", n_clusters_dev))

      # Compute mean |coefficient| for each cluster
      cluster_means_dev <- sapply(1:n_clusters_dev, function(k) {
        mean(abs_coefs_dev[cluster_assignments_dev == k])
      })

      # Identify strongest cluster (highest mean |β|)
      strongest_cluster_dev <- which.max(cluster_means_dev)

      # Display cluster statistics
      message("  Cluster statistics:")
      for (k in 1:n_clusters_dev) {
        n_genes_k <- sum(cluster_assignments_dev == k)
        mean_k <- cluster_means_dev[k]
        sd_k <- sd(abs_coefs_dev[cluster_assignments_dev == k])
        message(sprintf("    Cluster %d: %d genes, mean |coef| = %.4f (SD = %.4f)%s",
                        k, n_genes_k, mean_k, sd_k,
                        ifelse(k == strongest_cluster_dev, " <- STRONGEST", "")))
      }

      # Keep genes in strongest cluster
      candidate_genes <- names(dev_coefs_significant)[cluster_assignments_dev == strongest_cluster_dev]

      message(sprintf("  Genes in strongest cluster (Cluster %d): %d out of %d",
                      strongest_cluster_dev, length(candidate_genes),
                      length(dev_coefs_significant)))

      # Save GMM results for development set
      gmm_dev_results <- list(
        model = gmm_dev,
        n_clusters = n_clusters_dev,
        cluster_means = cluster_means_dev,
        strongest_cluster = strongest_cluster_dev,
        cluster_assignments = cluster_assignments_dev
      )

    } else {
      # GMM failed, fall back to median approach
      message("  Warning: GMM fitting failed, using median |coefficient| threshold")
      median_coef_dev <- median(abs(dev_coefs_significant))
      candidate_genes <- names(dev_coefs_significant)[abs(dev_coefs_significant) > median_coef_dev]
      gmm_dev_results <- NULL

      message(sprintf("  Median |coefficient|: %.4f", median_coef_dev))
      message(sprintf("  Genes with |coef| > median: %d out of %d",
                      length(candidate_genes), length(dev_coefs_significant)))
    }

  } else {
    # Too few genes for GMM, use median approach
    message(sprintf("  Too few genes (%d) for GMM, using median threshold",
                    length(dev_coefs_significant)))
    median_coef_dev <- median(abs(dev_coefs_significant))
    candidate_genes <- names(dev_coefs_significant)[abs(dev_coefs_significant) > median_coef_dev]
    gmm_dev_results <- NULL

    message(sprintf("  Median |coefficient|: %.4f", median_coef_dev))
    message(sprintf("  Genes with |coef| > median: %d out of %d",
                    length(candidate_genes), length(dev_coefs_significant)))
  }
} else {
  candidate_genes <- character(0)
  gmm_dev_results <- NULL
}

if (length(candidate_genes) == 0) {
  stop("No genes passed development set screening. Analysis cannot proceed.")
}



screened_genes <- candidate_genes

message(sprintf("  Final screened genes: %d", length(screened_genes)))

# Display filtering summary
message("\n=== Screening Summary ===")
if (!is.null(gmm_dev_results)) {
  message(sprintf("Development set: %d → %d genes (p<0.05 + GMM strongest cluster)",
                  length(significant_genes_dev), length(screened_genes)))
  message(sprintf("  GMM: %d clusters identified, keeping Cluster %d (highest mean |coef|)",
                  gmm_dev_results$n_clusters, gmm_dev_results$strongest_cluster))
} else {
  message(sprintf("Development set: %d → %d genes (p<0.05 + |coef|>median)",
                  length(significant_genes_dev), length(screened_genes)))
}














if (length(screened_genes) == 0) {
  stop("No genes were significant in BOTH development and validation sets. Analysis cannot proceed.")
}

# Save comprehensive screening results
if (!dir.exists("exploratory_analysis_elastic_net_results")) {
  dir.create("exploratory_analysis_elastic_net_results", recursive = TRUE)
}
saveRDS(list(
  # Development set results
  dev_significant_genes = significant_genes_dev,
  dev_candidate_genes = candidate_genes,
  dev_pvalues = pvals[candidate_genes],
  dev_coefficients = dev_coefs[candidate_genes],
  dev_median_coef = if(exists("median_coef_dev")) median_coef_dev else NA,

  # GMM results (if used)
  gmm_results = gmm_dev_results,

  # Final screened genes
  screened_genes = screened_genes
), "exploratory_analysis_elastic_net_results/screening_results.RDS")













# Step 3: Fit Elastic Net 100 Times on Development Set with Screened Genes
message("\n=== Step 3: Elastic Net Selection (100 replications on screened genes) ===")
message("Using clinical covariates: age, sex_binary (unpenalized) - NOTE: race_binary, bmi_binary added only in final prognostic model")
message("Using study_indicator as part of penalized feature set (like genes)")

# Combine screened genes with all clinical covariates for elastic net
# NOTE: study_indicator is added to the gene/feature matrix (penalized), not clinical covariates (unpenalized)
X_dev_screened <- cbind(Z_dev_screening, study_indicator = study_indicator_dev, X_dev_genes[, screened_genes, drop = FALSE])

message(sprintf("Running 100 elastic net fits with %d screened genes + study_indicator (penalized) + %d clinical covariates (age, sex_binary - unpenalized)...",
                length(screened_genes), ncol(Z_dev_screening)))

# IMPORTANT: Create penalty factors to NOT penalize clinical covariates
# This ensures clinical coefficients retain correct signs
# and are not artificially shrunk or flipped due to gene collinearity
# NOTE: study_indicator IS penalized (like genes) to allow for potential shrinkage/selection
penalty_factors_cv <- c(
  rep(0, ncol(Z_dev_screening)),           # 0 = NO penalty for clinical covariates (age, sex_binary only)
  1,                                        # 1 = standard penalty for study_indicator
  rep(1, length(screened_genes))           # 1 = standard penalty for genes (subject to selection)
)

message(sprintf("Penalty structure: %d unpenalized clinical covariates (age, sex_binary), 1 penalized study_indicator, %d penalized genes",
                ncol(Z_dev_screening), length(screened_genes)))

# Initialize selection tracking matrices
n_iterations <- 100
# Include study_indicator in gene selection matrix since it's now penalized
gene_selection_matrix <- matrix(0, nrow = n_iterations, ncol = length(screened_genes) + 1)
colnames(gene_selection_matrix) <- c("study_indicator", screened_genes)
covariate_selection_matrix <- matrix(0, nrow = n_iterations, ncol = ncol(Z_dev_screening))
colnames(covariate_selection_matrix) <- colnames(Z_dev_screening)

# Run elastic net 100 times with different random seeds
for (rep in 1:n_iterations) {
  if (rep %% 10 == 0) message(sprintf("  Replication %d/%d", rep, n_iterations))

  # Set seed for reproducibility of this replication
  set.seed(789 + rep)

  # Fit elastic net using logistic regression with separate penalties
  # Use cv.glmnet to select optimal lambda via internal cross-validation
  cv_enet <- cv.glmnet(
    x = X_dev_screened,
    y = y_dev,
    family = "binomial",     # Logistic regression (logit link)
    alpha = 0.5,             # Elastic net with equal L1/L2 balance
    penalty.factor = penalty_factors_cv,  # KEY: No penalty on clinical covariates
    nfolds = 10,             # 10-fold CV for lambda selection
    type.measure = "deviance"
  )
  
  # Extract selected variables at lambda.1se
  coefs <- coef(cv_enet, s = "lambda.1se")[-1, ]  # Remove intercept
  selected_vars <- names(coefs)[coefs != 0]

  # Separate genes (including study_indicator) and clinical covariates
  # study_indicator is now tracked with genes since it's penalized
  selected_genes_rep <- selected_vars[selected_vars %in% c("study_indicator", screened_genes)]
  selected_covars_rep <- selected_vars[selected_vars %in% colnames(Z_dev_screening)]  # Clinical covariates only (age, sex, race, bmi)

  # Record gene selections (including study_indicator)
  if (length(selected_genes_rep) > 0) {
    gene_selection_matrix[rep, selected_genes_rep] <- 1
  }

  # Record covariate selections (age, sex, race, bmi - NOT study_indicator)
  if (length(selected_covars_rep) > 0) {
    covariate_selection_matrix[rep, selected_covars_rep] <- 1
  }
}

# Aggregate results: keep genes/covariates selected in ≥50 replications (more stringent)
message("\n=== Aggregating Selection Results ===")
gene_selection_freq <- colSums(gene_selection_matrix)
covariate_selection_freq <- colSums(covariate_selection_matrix)

# Use higher threshold for more stringent selection
selection_threshold <- 50

final_selected_genes <- names(gene_selection_freq)[gene_selection_freq >= selection_threshold]
final_selected_covariates <- names(covariate_selection_freq)[covariate_selection_freq >= selection_threshold]

# Check if study_indicator was selected
study_indicator_selected <- "study_indicator" %in% final_selected_genes
# Separate study_indicator from actual genes for reporting
final_selected_genes_only <- final_selected_genes[final_selected_genes != "study_indicator"]

message(sprintf("Genes selected in ≥%d replications: %d out of %d screened",
                selection_threshold, length(final_selected_genes_only), length(screened_genes)))
message(sprintf("Study indicator selected: %s (frequency: %d/100)",
                ifelse(study_indicator_selected, "YES", "NO"), gene_selection_freq["study_indicator"]))
message(sprintf("Covariates selected in ≥%d replications: %d out of %d (age, sex_binary)",
                selection_threshold, length(final_selected_covariates), ncol(Z_dev_screening)))

if (length(final_selected_covariates) > 0) {
  message(sprintf("  Covariate names: %s", paste(final_selected_covariates, collapse = ", ")))
}

# Save results
if (!dir.exists("exploratory_analysis_elastic_net_results")) {
  dir.create("exploratory_analysis_elastic_net_results", recursive = TRUE)
}

saveRDS(final_selected_genes, "exploratory_analysis_elastic_net_results/final_selected_genes.RDS")
saveRDS(final_selected_covariates, "exploratory_analysis_elastic_net_results/final_selected_covariates.RDS")
saveRDS(gene_selection_freq, "exploratory_analysis_elastic_net_results/gene_selection_frequency.RDS")
saveRDS(covariate_selection_freq, "exploratory_analysis_elastic_net_results/covariate_selection_frequency.RDS")

# Save selection matrix for detailed analysis
saveRDS(list(
  gene_selection_matrix = gene_selection_matrix,
  covariate_selection_matrix = covariate_selection_matrix,
  gene_selection_freq = gene_selection_freq,
  covariate_selection_freq = covariate_selection_freq,
  final_selected_genes = final_selected_genes,
  final_selected_covariates = final_selected_covariates
), "exploratory_analysis_elastic_net_results/elastic_net_selection_details.RDS")





################################################################################################################
# 4a) Marginal Analysis: Y vs Each Clinical Covariate (Univariate)
################################################################################################################

message("\n=== Step 4a: Marginal Analysis of Clinical Covariates ===")
message("Fitting univariate logistic regression: Y ~ covariate (one at a time)")
message("Purpose: Investigate if coefficient signs are correct BEFORE adding genes")
message("Note: Analyzing all 4 clinical covariates (age, sex_binary, race_binary, bmi_binary)")
message("      These will all be used as unpenalized covariates in the final prognostic model\n")

marginal_results <- data.frame(
  covariate = character(),
  coefficient = numeric(),
  std_error = numeric(),
  z_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (covar in colnames(Z_dev_full)) {
  # Fit marginal model: Y ~ covariate
  marginal_data <- data.frame(
    outcome = y_dev,
    covariate = Z_dev_full[, covar]
  )

  marginal_fit <- glm(
    outcome ~ covariate,
    data = marginal_data,
    family = binomial(link = "logit")
  )

  # Extract coefficient for the covariate (not intercept)
  coef_summary <- summary(marginal_fit)$coefficients
  coef_row <- coef_summary["covariate", ]

  marginal_results <- rbind(marginal_results, data.frame(
    covariate = covar,
    coefficient = coef_row["Estimate"],
    std_error = coef_row["Std. Error"],
    z_value = coef_row["z value"],
    p_value = coef_row["Pr(>|z|)"],
    stringsAsFactors = FALSE
  ))
}

# Display marginal results
message("=== Marginal (Univariate) Results ===")
for (i in 1:nrow(marginal_results)) {
  covar <- marginal_results$covariate[i]
  coef_val <- marginal_results$coefficient[i]
  se <- marginal_results$std_error[i]
  p_val <- marginal_results$p_value[i]
  sig_star <- ifelse(p_val < 0.001, "***",
                     ifelse(p_val < 0.01, "**",
                            ifelse(p_val < 0.05, "*", "")))

  # Flag if sign is negative for age or BMI
  warning_flag <- ""
  if ((covar %in% c("age", "bmi_binary")) && coef_val < 0) {
    warning_flag <- " <- WARNING: Negative (unexpected!)"
  }

  message(sprintf("  %s: coef = %7.4f (SE = %.4f, p = %.4f) %s%s",
                  covar, coef_val, se, p_val, sig_star, warning_flag))
}

# Save marginal results
saveRDS(marginal_results, "exploratory_analysis_elastic_net_results/marginal_clinical_results.RDS")
message("\nMarginal results saved to: marginal_clinical_results.RDS\n")


################################################################################################################
# 4b) Run a final prognostic logistic regression model on the development set using the selected genes and covariates
################################################################################################################

message("\n=== Step 4b: Final Prognostic Elastic Net Model ===")
message("NOTE: Now using ALL 4 clinical covariates (age, sex_binary, race_binary, bmi_binary) as unpenalized")
message("      This differs from screening/elastic net which only used age and sex_binary")

if (length(final_selected_genes) > 0) {
  # Prepare development set with selected genes + study_indicator (penalized) + clinical covariates (unpenalized)
  # NOTE: final_selected_genes may include "study_indicator" if it was selected
  # We need to handle study_indicator separately since it's not in X_dev_genes
  genes_only_selected <- final_selected_genes[final_selected_genes != "study_indicator"]
  X_dev_selected_genes <- X_dev_genes[, genes_only_selected, drop = FALSE]

  # Always include study_indicator in the penalized part (regardless of whether it passed selection threshold)
  # This allows it to be shrunk/selected in the final model
  # NOTE: Using Z_dev_full (all 4 covariates) for final prognostic model
  X_dev_final <- cbind(Z_dev_full, study_indicator = study_indicator_dev, X_dev_selected_genes)

  message(sprintf("Fitting elastic net with %d genes + study_indicator (penalized) + %d clinical covariates (unpenalized)...",
                  length(genes_only_selected), ncol(Z_dev_full)))
  message("Using elastic net (alpha=0.5) with SEPARATE penalties:")
  message("  - Clinical covariates (age, sex_binary, race_binary, bmi_binary): NO penalty (forced in model)")
  message("  - Study indicator + Genes: Standard penalty (subject to selection)")

  # IMPORTANT: Create penalty factors for final prognostic model
  # Clinical covariates are NOT penalized to ensure stable, interpretable coefficients
  # study_indicator IS penalized (like genes)
  penalty_factors_final <- c(
    rep(0, ncol(Z_dev_full)),                  # 0 = NO penalty for clinical covariates (age, sex, race, bmi_binary)
    1,                                          # 1 = standard penalty for study_indicator
    rep(1, length(genes_only_selected))        # 1 = standard penalty for genes
  )

  message(sprintf("Penalty structure: %d unpenalized clinical covariates, 1 penalized study_indicator, %d penalized genes",
                  ncol(Z_dev_full), length(genes_only_selected)))

  # Fit elastic net logistic regression with cv for lambda selection
  # Use cv.glmnet with alpha=0.5 to get coefficients and predictions
  set.seed(456)
  prognostic_model_cv <- cv.glmnet(
    x = X_dev_final,
    y = y_dev,
    family = "binomial",
    alpha = 0.5,  # Elastic net mixing parameter (equal L1/L2 balance)
    penalty.factor = penalty_factors_final,  # KEY: No penalty on clinical covariates
    nfolds = 10,
    type.measure = "deviance"
  )

  # Extract coefficients at lambda.1se
  coef_final_prognostic <- coef(prognostic_model_cv, s = "lambda.1se")
  coef_values <- as.vector(coef_final_prognostic)[-1]  # Remove intercept
  coef_names <- rownames(coef_final_prognostic)[-1]

  # Create coefficient table
  # Identify important predictors by coefficient magnitude
  coef_table <- data.frame(
    variable = coef_names,
    coefficient = coef_values,
    abs_coefficient = abs(coef_values),
    stringsAsFactors = FALSE
  )
  coef_table <- coef_table[order(coef_table$abs_coefficient, decreasing = TRUE), ]

  # For prognostic model, keep predictors with non-zero coefficients
  significant_predictors <- coef_table$variable[coef_table$abs_coefficient != 0]

  # Separate significant genes (including study_indicator if selected) and clinical covariates
  # study_indicator is now in the "gene" group (penalized)
  significant_genes <- significant_predictors[significant_predictors %in% c("study_indicator", genes_only_selected)]
  significant_covariates <- significant_predictors[significant_predictors %in% colnames(Z_dev_full)]

  # Further separate study_indicator from actual genes for reporting
  significant_genes_only <- significant_genes[significant_genes != "study_indicator"]
  study_indicator_significant <- "study_indicator" %in% significant_genes

  # Store the glmnet model as prognostic_model
  prognostic_model <- prognostic_model_cv

  message(sprintf("\nSignificant predictors (|coef| > 0):"))
  message(sprintf("  Genes: %d out of %d", length(significant_genes_only), length(genes_only_selected)))
  message(sprintf("  Study indicator: %s", ifelse(study_indicator_significant, "YES (non-zero coefficient)", "NO (shrunk to zero)")))
  message(sprintf("  Clinical covariates: %d out of %d (age, sex_binary, race_binary, bmi_binary - all unpenalized in final model)", length(significant_covariates), ncol(Z_dev_full)))

  # Display clinical covariates from elastic net
  message("\n=== Clinical Covariates (Elastic Net with Genes) ===")
  covar_coefs <- coef_table[coef_table$variable %in% colnames(Z_dev_full), ]
  for (i in 1:nrow(covar_coefs)) {
    covar <- covar_coefs$variable[i]
    coef_val <- covar_coefs$coefficient[i]

    # Flag if sign is negative for age or BMI_binary
    warning_flag <- ""
    if ((covar %in% c("age", "bmi_binary")) && coef_val < 0) {
      warning_flag <- " <- WARNING: Negative (compare with marginal!)"
    }

    message(sprintf("  %s: coef = %7.4f%s",
                    covar, coef_val, warning_flag))
  }

  # Display study_indicator coefficient
  message("\n=== Study Indicator (Penalized) ===")
  study_ind_coef <- coef_table[coef_table$variable == "study_indicator", ]
  if (nrow(study_ind_coef) > 0) {
    message(sprintf("  study_indicator: coef = %7.4f %s",
                    study_ind_coef$coefficient,
                    ifelse(study_ind_coef$coefficient == 0, "(shrunk to zero)", "")))
  }

  if (length(significant_genes_only) > 0) {
    message("\n=== Significant Genes (|coef| > 0) ===")
    gene_coefs <- coef_table[coef_table$variable %in% significant_genes_only, ]
    for (i in 1:min(20, nrow(gene_coefs))) {
      gene <- gene_coefs$variable[i]
      coef_val <- gene_coefs$coefficient[i]
      message(sprintf("  %2d. %s: coef = %7.4f", i, gene, coef_val))
    }
    if (nrow(gene_coefs) > 20) {
      message(sprintf("  ... and %d more genes with non-zero coefficients", nrow(gene_coefs) - 20))
    }
  } else {
    message("\n=== Significant Genes ===")
    message("  No genes with non-zero coefficients")
  }
  
  # Save prognostic model results
  saveRDS(list(
    model = prognostic_model,
    model_type = "glmnet",  # Elastic net model
    lambda_min = prognostic_model$lambda.min,
    lambda_1se = prognostic_model$lambda.1se,
    coefficients = coef_table,
    significant_predictors = significant_predictors,
    significant_genes = significant_genes,  # Includes study_indicator if selected
    significant_genes_only = significant_genes_only,  # Actual genes only (no study_indicator)
    significant_covariates = significant_covariates,
    study_indicator_significant = study_indicator_significant,
    selected_genes_input = final_selected_genes,  # Includes study_indicator if selected
    genes_only_selected = genes_only_selected,  # Actual genes only (no study_indicator)
    all_covariates = colnames(Z_dev_full),  # 4 clinical covariates (age, sex, race, bmi_binary)
    screening_covariates = colnames(Z_dev_screening),  # Same as all_covariates
    study_indicator_in_X = TRUE,  # Flag indicating study_indicator is in X (penalized)
    marginal_results = marginal_results  # Include marginal analysis
  ), "exploratory_analysis_elastic_net_results/prognostic_model_results.RDS")
  
  message("\nPrognostic model results saved to: prognostic_model_results.RDS")
  
} else {
  warning("No genes selected from elastic net. Cannot fit prognostic model.")
  significant_genes <- character(0)
  significant_covariates <- character(0)
  prognostic_model <- NULL
}


################################################################################################################
# 5) Evaluate Model on Final Validation Set
################################################################################################################
# Using retained genes, evaluate predictive performance on held-out validation set
# Metrics: AUC, Sensitivity, Specificity, Brier Score
################################################################################################################

if (!is.null(prognostic_model)) {
  message(sprintf("\n=== Step 5: Final Validation Set Evaluation ==="))
  message(sprintf("Evaluating prognostic model with %d genes + study_indicator (penalized) + %d clinical covariates (unpenalized)",
                  length(genes_only_selected), ncol(Z_dev_full)))

  # Prepare validation set with selected genes + study_indicator (penalized) + clinical covariates (unpenalized)
  # Must match the structure used in training: Z_val_full + study_indicator + genes
  X_val_selected_genes <- X_val_genes[, genes_only_selected, drop = FALSE]
  X_val_final <- cbind(Z_val_full, study_indicator = study_indicator_val, X_val_selected_genes)

  # Predict on validation set using elastic net prognostic model
  pred_prob_val <- as.vector(predict(prognostic_model, newx = X_val_final,
                                     s = "lambda.1se", type = "response"))
  pred_class_val <- ifelse(pred_prob_val > 0.5, 1, 0)
  
  # Performance metrics
  message("\n=== Validation Set Performance ===")
  
  # 1. AUC
  roc_val <- roc(y_val, pred_prob_val, quiet = TRUE, levels = c(0, 1), direction = "<")
  auc_val <- as.numeric(auc(roc_val))
  message(sprintf("AUC: %.4f (95%% CI: %.4f - %.4f)",
                  auc_val, ci.auc(roc_val)[1], ci.auc(roc_val)[3]))
  
  # 2. Sensitivity and Specificity
  confusion_matrix <- table(Predicted = pred_class_val, Actual = y_val)
  message("\nConfusion Matrix:")
  print(confusion_matrix)
  
  if (all(c(0, 1) %in% pred_class_val) && all(c(0, 1) %in% y_val)) {
    sensitivity <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
    specificity <- confusion_matrix[1, 1] / sum(confusion_matrix[, 1])
    ppv <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
    npv <- confusion_matrix[1, 1] / sum(confusion_matrix[1, ])
    
    message(sprintf("\nSensitivity: %.4f", sensitivity))
    message(sprintf("Specificity: %.4f", specificity))
    message(sprintf("PPV: %.4f", ppv))
    message(sprintf("NPV: %.4f", npv))
  } else {
    sensitivity <- NA
    specificity <- NA
    ppv <- NA
    npv <- NA
    message("\nWarning: Cannot compute sensitivity/specificity (missing prediction class)")
  }
  
  # 3. Brier Score
  brier_score <- mean((y_val - pred_prob_val)^2)
  message(sprintf("\nBrier Score: %.4f", brier_score))
  
  # 4. Accuracy
  accuracy_val <- mean(pred_class_val == y_val)
  message(sprintf("Accuracy: %.4f", accuracy_val))
  
  # Use coefficient table from Step 4 
  coef_df_final <- coef_table
  
  # Separate genes (excluding study_indicator), study_indicator, and clinical covariates
  gene_coefs_final <- coef_df_final[coef_df_final$variable %in% genes_only_selected, ]
  study_ind_coef_final <- coef_df_final[coef_df_final$variable == "study_indicator", ]
  clinical_coefs_final <- coef_df_final[coef_df_final$variable %in% colnames(Z_dev_full), ]

  gene_coefs_final_for_paper_and_sim <- data.frame(gene = gene_coefs_final$variable, coefficient = gene_coefs_final$coefficient)
  gene_coefs_final_for_paper_and_sim <- gene_coefs_final_for_paper_and_sim[gene_coefs_final_for_paper_and_sim$coefficient != 0,]

  clinical_coefs_final_for_paper_and_sim <- data.frame(var = clinical_coefs_final$variable, coefficient = clinical_coefs_final$coefficient)

  # Also save study_indicator coefficient separately
  study_ind_coef_for_paper_and_sim <- data.frame(var = "study_indicator",
                                                  coefficient = ifelse(nrow(study_ind_coef_final) > 0,
                                                                       study_ind_coef_final$coefficient, 0))

  message(sprintf("\n=== Prognostic Model Summary  ==="))
  message(sprintf("Total predictors: %d genes + study_indicator (penalized) + %d clinical covariates (unpenalized)",
                  length(genes_only_selected), ncol(Z_dev_full)))
  message(sprintf("  Important genes (|coef| > 0): %d", length(significant_genes_only)))
  message(sprintf("  Study indicator: %s", ifelse(study_indicator_significant, "selected (non-zero)", "not selected (zero)")))
  message(sprintf("  Important covariates (|coef| > 0): %d", length(significant_covariates)))
  
  if (nrow(gene_coefs_final) > 0) {
    message("\nTop 10 genes by coefficient magnitude:")
    print(head(gene_coefs_final[, c("variable", "coefficient", "abs_coefficient")], 10))
  }
  
  if (nrow(clinical_coefs_final) > 0) {
    message("\nClinical covariate coefficients:")
    print(clinical_coefs_final[, c("variable", "coefficient", "abs_coefficient")])
  }
  
  # Save all results
  saveRDS(prognostic_model, "exploratory_analysis_elastic_net_results/final_prognostic_model.RDS")
  saveRDS(coef_df_final, "exploratory_analysis_elastic_net_results/final_model_coefficients.RDS")
  saveRDS(gene_coefs_final, "exploratory_analysis_elastic_net_results/final_model_gene_coefficients.RDS")
  saveRDS(study_ind_coef_final, "exploratory_analysis_elastic_net_results/final_model_study_indicator_coefficient.RDS")
  saveRDS(clinical_coefs_final, "exploratory_analysis_elastic_net_results/final_model_clinical_coefficients.RDS")

  saveRDS(gene_coefs_final_for_paper_and_sim, "exploratory_analysis_elastic_net_results/gene_coefs_final_for_paper_and_sim.RDS")
  saveRDS(study_ind_coef_for_paper_and_sim, "exploratory_analysis_elastic_net_results/study_ind_coef_for_paper_and_sim.RDS")
  saveRDS(clinical_coefs_final_for_paper_and_sim, "exploratory_analysis_elastic_net_results/clinical_coefs_final_for_paper_and_sim.RDS")
  
  
  # Save validation predictions
  validation_results <- data.frame(
    patient_id = val_ids,
    true_outcome = y_val,
    predicted_prob = pred_prob_val,
    predicted_class = pred_class_val
  )
  saveRDS(validation_results, "exploratory_analysis_elastic_net_results/validation_predictions.RDS")
  
  # Save performance metrics
  validation_performance <- data.frame(
    metric = c("AUC", "Sensitivity", "Specificity", "PPV", "NPV", "Brier_Score", "Accuracy"),
    value = c(auc_val, sensitivity, specificity, ppv, npv, brier_score, accuracy_val),
    stringsAsFactors = FALSE
  )
  saveRDS(validation_performance, "exploratory_analysis_elastic_net_results/validation_performance.RDS")
  
  message("\n=== Results saved to: exploratory_analysis_elastic_net_results/ ===")
  message("\n===============================================================")
  message("=== ANALYSIS COMPLETE ===")
  message("===============================================================")
  message(sprintf("Screening: %d genes passed two-stage screening", length(screened_genes)))
  message(sprintf("Elastic Net: %d genes + study_indicator selected (≥%d/100 replications)",
                  length(final_selected_genes_only), selection_threshold))
  message(sprintf("  Study indicator selected in stability selection: %s",
                  ifelse(study_indicator_selected, "YES", "NO")))
  message(sprintf("Prognostic Model: %d significant genes (|coef| > 0)", length(significant_genes_only)))
  message(sprintf("  Study indicator in final model: %s",
                  ifelse(study_indicator_significant, "YES (non-zero)", "NO (shrunk to zero)")))
  message(sprintf("Validation AUC: %.4f", auc_val))
  message("===============================================================")
  message("NOTE: study_indicator is penalized (in X vector) like genes.")
  message("      Screening/Elastic Net: Only age, sex_binary as unpenalized covariates")
  message("      Final Prognostic Model: All 4 covariates (age, sex_binary, race_binary, bmi_binary) unpenalized")
  message("===============================================================")

  # ========================================================================
  # Step 6: Visualize Gene Expression vs Case/Control Ratio
  # NON-OVERLAPPING WINDOW APPROACH (statistically independent bins)
  # 1) Line up all 500 patients from lowest gene expression -> highest gene expression
  # 2) Look at the first 25 patients (lowest expression). Count how many are cases
  # 3) Move over by 25 patients. Look at patients 26-50. Count how many are cases
  # 4) Keep sliding: patients 51-75, 76-100, 101-125, etc. (no overlap)
  # 5) For each window, compute:
  #      case_ratio = number_of_cases / total_patients_in_window
  #      X-axis = median gene expression in that window
  #      Y-axis = case ratio (proportion of cases)
  # What we see:
  # i) If line goes UP -> higher expression = more cases = gene predicts disease
  # ii) If line goes DOWN -> higher expression = fewer cases = gene protects against disease
  # iii) If line is FLAT -> expression doesn't matter = gene not predictive
  # ========================================================================
  message("\n=== Step 6: Gene Expression Trend Visualization ===")

  # Visualize ONLY genes with non-zero coefficients in final prognostic model
  # Exclude genes with coefficient = 0
  gene_coefs_nonzero <- gene_coefs_final[gene_coefs_final$coefficient != 0, ]

  # Sort genes by absolute coefficient
  gene_coefs_sorted <- gene_coefs_nonzero[order(-gene_coefs_nonzero$abs_coefficient), ]
  genes_to_plot <- gene_coefs_sorted$variable
  n_genes_to_plot <- length(genes_to_plot)

  if (n_genes_to_plot > 0) {
    message(sprintf("\nCreating non-overlapping window visualization for %d genes with non-zero coefficients...", n_genes_to_plot))
    message(sprintf("(Excluding %d genes with zero coefficients from elastic net)",
                    nrow(gene_coefs_final) - n_genes_to_plot))
    message(sprintf("(Using 25-patient windows with no overlap for statistical independence)"))
    message(sprintf("(This will create a multi-page PDF with 12 genes per page)"))

    # Function: compute case/control ratio with NON-OVERLAPPING windows
    compute_case_control_ratio <- function(gene_expr, outcome, window_size = 25) {
      # Sort by gene expression
      order_idx <- order(gene_expr)
      gene_sorted <- gene_expr[order_idx]
      outcome_sorted <- outcome[order_idx]

      n <- length(gene_expr)
      # Non-overlapping windows: step size = window size
      n_windows <- floor(n / window_size)

      window_midpoints <- numeric(n_windows)
      case_ratios <- numeric(n_windows)
      case_counts <- numeric(n_windows)
      control_counts <- numeric(n_windows)

      for (i in 1:n_windows) {
        # Non-overlapping: each window starts where previous ended
        start_idx <- (i - 1) * window_size + 1
        end_idx <- i * window_size

        window_outcomes <- outcome_sorted[start_idx:end_idx]
        window_exprs <- gene_sorted[start_idx:end_idx]

        # Midpoint of expression window
        window_midpoints[i] <- median(window_exprs)

        # Count cases and controls
        n_cases <- sum(window_outcomes == 1)
        n_controls <- sum(window_outcomes == 0)

        case_counts[i] <- n_cases
        control_counts[i] <- n_controls

        # Case ratio (proportion of cases in window)
        case_ratios[i] <- n_cases / (n_cases + n_controls)
      }

      return(list(
        midpoints = window_midpoints,
        case_ratio = case_ratios,
        n_cases = case_counts,
        n_controls = control_counts
      ))
    }

    # Create visualization
    pdf("exploratory_analysis_elastic_net_results/gene_expression_case_control_trends.pdf",
        width = 16, height = 12)

    # Plot 12 genes per page (3 rows × 4 columns)
    genes_per_page <- 12
    n_pages <- ceiling(n_genes_to_plot / genes_per_page)

    overall_case_ratio <- mean(y_dev)

    for (page in 1:n_pages) {
      # Determine genes for this page
      start_idx <- (page - 1) * genes_per_page + 1
      end_idx <- min(page * genes_per_page, n_genes_to_plot)
      genes_this_page <- genes_to_plot[start_idx:end_idx]

      # Set up layout for this page
      par(mfrow = c(3, 4), mar = c(4, 4, 3, 1), oma = c(2, 2, 3, 1))

      for (gene in genes_this_page) {
        # Get gene expression from development set
        gene_expr <- X_dev_genes[, gene]

        # Get coefficient for this gene
        gene_coef <- gene_coefs_sorted[gene_coefs_sorted$variable == gene, "coefficient"]

        # Compute non-overlapping window statistics
        window_stats <- compute_case_control_ratio(
          gene_expr = gene_expr,
          outcome = y_dev,
          window_size = 25  # 25 patients per window, no overlap
        )

        # Determine color based on coefficient sign
        plot_col <- ifelse(gene_coef > 0, "darkblue", "darkred")
        coef_direction <- ifelse(gene_coef > 0, "↑ Cases", "↓ Cases")

        # Plot case ratio vs gene expression
        plot(window_stats$midpoints, window_stats$case_ratio,
             type = "o", pch = 19, col = plot_col, lwd = 2,
             xlab = "Normalized Gene Expression (25-pt non-overlapping bins, low -> high)",
             ylab = "Case Ratio (proportion)",
             main = sprintf("%s (Coeff = %.3f, %s)", gene, gene_coef, coef_direction),
             ylim = c(0, 1), cex.main = 0.9)

        # Add horizontal line at overall case ratio
        abline(h = overall_case_ratio, col = "gray50", lty = 2, lwd = 1.5)

        # Add smoothed trend line
        if (length(window_stats$midpoints) > 3) {
          smooth_fit <- loess(window_stats$case_ratio ~ window_stats$midpoints, span = 0.5)
          lines(window_stats$midpoints, predict(smooth_fit), col = "orange", lwd = 3)
        }

        # Add legend on first plot of first page only
        if (page == 1 && gene == genes_this_page[1]) {
          legend("topleft",
                 legend = c("Non-overlapping bins", "Smooth trend", "Overall case ratio"),
                 col = c(plot_col, "orange", "gray50"),
                 lty = c(1, 1, 2), lwd = c(2, 3, 1.5),
                 pch = c(19, NA, NA), bty = "n", cex = 0.7)
        }

        # Add grid
        grid(col = "gray80", lty = "dotted")
      }

      # Overall title for this page
      mtext(sprintf("Gene Expression vs Case/Control Ratio (Page %d of %d)", page, n_pages),
            outer = TRUE, cex = 1.5, font = 2, line = 0.5)
    }

    dev.off()

    message(sprintf("\nTrend visualization saved to: exploratory_analysis_elastic_net_results/gene_expression_case_control_trends.pdf"))
    message(sprintf("  Total: %d genes across %d pages (12 genes per page)", n_genes_to_plot, n_pages))

    # Print interpretation guide
    message("\n=== Interpretation Guide ===")
    message("Windowing approach:")
    message("  - NON-OVERLAPPING bins of 25 patients each")
    message("  - Windows: 1-25, 26-50, 51-75, ..., 476-500")
    message("  - Total: ~20 independent bins per gene")
    message("  - Each patient appears in exactly ONE bin (statistical independence)")
    message("\nColor coding:")
    message("  - BLUE: Positive coefficient (gene ↑ → case probability ↑)")
    message("  - RED: Negative coefficient (gene ↑ → case probability ↓)")
    message("\nFor predictive genes, you should see:")
    message("  - Blue genes: UPWARD trend (gene predicts cases)")
    message("  - Red genes: DOWNWARD trend (gene predicts controls)")
    message("  - Flat trend suggests gene may not be strongly predictive")
    message("\nPlot elements:")
    message("  - Gray dashed line = overall case ratio in development set")
    message("  - Colored points = case ratio in each 25-patient bin (independent)")
    message("  - Orange line = smoothed trend (LOESS)")
    message(sprintf("\nGenes are sorted by absolute coefficient (strongest first)"))

  } else {
    message("\nNo genes available for visualization.")
  }

} else {
  message("\n=== ERROR: No genes were selected ===")
  message("Cannot proceed with validation. Check earlier steps.")
}

message("\n=== Multi-Stage Variable Selection Pipeline Completed! ===")





################################################################################################################
# 5) Generate LaTeX Tables for Results
################################################################################################################

message("\n=== Step 5: Generating LaTeX Tables ===")

# Function to generate LaTeX table from coefficient dataframe
generate_coef_latex_table <- function(coef_data, caption, label, filename) {
  if (nrow(coef_data) == 0) {
    message(sprintf("Skipping %s: no data", filename))
    return(NULL)
  }
  
  # Select only variable and coefficient columns
  table_data <- coef_data[, c("variable", "coefficient"), drop = FALSE]
  
  # Format coefficients to 4 decimal places
  table_data$coefficient <- sprintf("%.4f", table_data$coefficient)
  
  # Open file connection
  sink(filename)
  
  cat("\\begin{table}[ht]\n")
  cat("\\centering\n")
  cat(sprintf("\\caption{%s}\n", caption))
  cat(sprintf("\\label{%s}\n", label))
  cat("\\begin{tabular}{l r}\n")
  cat("\\hline\n")
  cat("\\textbf{Variable} & \\textbf{Coefficient} \\\\\n")
  cat("\\hline\n")
  
  # Add rows
  for (i in 1:nrow(table_data)) {
    # Escape underscores in variable names for LaTeX
    var_name <- gsub("_", "\\\\_", table_data$variable[i])
    cat(sprintf("%s & %s \\\\\n", var_name, table_data$coefficient[i]))
  }
  
  cat("\\hline\n")
  cat("\\end{tabular}\n")
  cat("\\end{table}\n")
  
  # Close file connection
  sink()
  
  message(sprintf("  LaTeX table saved: %s", filename))
}

# Function to generate performance metrics table
generate_performance_table <- function(perf_data, caption, label, filename) {
  sink(filename)
  
  cat("\\begin{table}[ht]\n")
  cat("\\centering\n")
  cat(sprintf("\\caption{%s}\n", caption))
  cat(sprintf("\\label{%s}\n", label))
  cat("\\begin{tabular}{l r}\n")
  cat("\\hline\n")
  cat("\\textbf{Metric} & \\textbf{Value} \\\\\n")
  cat("\\hline\n")
  
  for (i in 1:nrow(perf_data)) {
    metric_name <- gsub("_", "\\\\_", perf_data$metric[i])
    metric_value <- ifelse(is.na(perf_data$value[i]), "NA",
                           sprintf("%.4f", perf_data$value[i]))
    cat(sprintf("%s & %s \\\\\n", metric_name, metric_value))
  }
  
  cat("\\hline\n")
  cat("\\end{tabular}\n")
  cat("\\end{table}\n")
  
  sink()
  
  message(sprintf("  Performance table saved: %s", filename))
}

# Table 1: Final model coefficients (all)
if (exists("coef_df_final") && nrow(coef_df_final) > 0) {
  generate_coef_latex_table(
    coef_data = coef_df_final,
    caption = "Non-zero coefficients from final elastic net model on validation set",
    label = "tab:final_model_coefficients",
    filename = "exploratory_analysis_elastic_net_results/table_final_coefficients.tex"
  )
}

# Table 2: Gene coefficients only
if (exists("gene_coefs_final") && nrow(gene_coefs_final) > 0) {
  generate_coef_latex_table(
    coef_data = gene_coefs_final,
    caption = "Gene coefficients from final elastic net model (genes only)",
    label = "tab:final_gene_coefficients",
    filename = "exploratory_analysis_elastic_net_results/table_final_gene_coefficients.tex"
  )
}

# Table 3: Clinical covariate coefficients only
if (exists("clinical_coefs_final") && nrow(clinical_coefs_final) > 0) {
  generate_coef_latex_table(
    coef_data = clinical_coefs_final,
    caption = "Clinical covariate coefficients from final elastic net model",
    label = "tab:final_clinical_coefficients",
    filename = "exploratory_analysis_elastic_net_results/table_final_clinical_coefficients.tex"
  )
}

# Table 4: Validation performance metrics
if (exists("validation_performance")) {
  generate_performance_table(
    perf_data = validation_performance,
    caption = "Predictive performance on final validation set (N=369)",
    label = "tab:validation_performance",
    filename = "exploratory_analysis_elastic_net_results/table_validation_performance.tex"
  )
}








true_genes_old <- c("ZBED3", "VSX1", "THRAP3", "TCP10L2", "ST5", "SPTLC1",
                    "SPEN", "SLC5A1", "RDH11", "PPP2R2D", "PPARD", "ORC1",
                    "NRIP3", "NR2E1", "MPHOSPH8", "GTF2I", "FAM53C", "ARHGAP40")




final_gene_set = readRDS("exploratory_analysis_elastic_net_results/gene_coefs_final_for_paper_and_sim.RDS")

true_genes_old_in_final_gene_set = true_genes_old[true_genes_old %in% final_gene_set$gene]








# ========================================================================
# Visualization: True Genes from Brian's Cox PH Model Paper
# ========================================================================
message("\n=== Creating Visualization for True Genes (Brian's Cox PH Model) ===")

# Check if true genes exist in the data
true_genes_in_data <- true_genes_old[true_genes_old %in% colnames(X_dev_genes)]
n_true_genes <- length(true_genes_in_data)

message(sprintf("Found %d out of %d true genes in the dataset", n_true_genes, length(true_genes_old)))

if (n_true_genes > 0) {
  # Reuse the same windowing function
  compute_case_control_ratio_true <- function(gene_expr, outcome, window_size = 25) {
    order_idx <- order(gene_expr)
    gene_sorted <- gene_expr[order_idx]
    outcome_sorted <- outcome[order_idx]
    
    n <- length(gene_expr)
    n_windows <- floor(n / window_size)
    
    window_midpoints <- numeric(n_windows)
    case_ratios <- numeric(n_windows)
    
    for (i in 1:n_windows) {
      start_idx <- (i - 1) * window_size + 1
      end_idx <- i * window_size
      
      window_outcomes <- outcome_sorted[start_idx:end_idx]
      window_exprs <- gene_sorted[start_idx:end_idx]
      
      window_midpoints[i] <- median(window_exprs)
      n_cases <- sum(window_outcomes == 1)
      n_controls <- sum(window_outcomes == 0)
      case_ratios[i] <- n_cases / (n_cases + n_controls)
    }
    
    return(list(midpoints = window_midpoints, case_ratio = case_ratios))
  }
  
  # Create PDF
  pdf("exploratory_analysis_elastic_net_results/true_genes_brian_cox_ph_trends.pdf",
      width = 16, height = 12)
  
  # Calculate number of pages needed
  genes_per_page <- 12
  n_pages <- ceiling(n_true_genes / genes_per_page)
  overall_case_ratio <- mean(y_dev)
  
  for (page in 1:n_pages) {
    start_idx <- (page - 1) * genes_per_page + 1
    end_idx <- min(page * genes_per_page, n_true_genes)
    genes_this_page <- true_genes_in_data[start_idx:end_idx]
    
    par(mfrow = c(3, 4), mar = c(4, 4, 3, 1), oma = c(2, 2, 4, 1))
    
    for (gene in genes_this_page) {
      # Get gene expression
      gene_expr <- X_dev_genes[, gene]
      
      # Compute window statistics
      window_stats <- compute_case_control_ratio_true(
        gene_expr = gene_expr,
        outcome = y_dev,
        window_size = 25
      )
      
      # Plot
      plot(window_stats$midpoints, window_stats$case_ratio,
           type = "o", pch = 19, col = "darkgreen", lwd = 2,
           xlab = "Gene Expression (25-pt non-overlapping bins, sorted low to high)",
           ylab = "Case Ratio (proportion)",
           main = gene,
           ylim = c(0, 1), cex.main = 1.0)
      
      # Add baseline
      abline(h = overall_case_ratio, col = "gray50", lty = 2, lwd = 1.5)
      
      # Add smoothed trend
      if (length(window_stats$midpoints) > 3) {
        smooth_fit <- loess(window_stats$case_ratio ~ window_stats$midpoints, span = 0.5)
        lines(window_stats$midpoints, predict(smooth_fit), col = "orange", lwd = 3)
      }
      
      # Legend on first plot only
      if (page == 1 && gene == genes_this_page[1]) {
        legend("topleft",
               legend = c("Non-overlapping bins", "Smooth trend", "Overall case ratio"),
               col = c("darkgreen", "orange", "gray50"),
               lty = c(1, 1, 2), lwd = c(2, 3, 1.5),
               pch = c(19, NA, NA), bty = "n", cex = 0.7)
      }
      
      grid(col = "gray80", lty = "dotted")
    }
    
    # Overall title
    mtext("Gene Expression vs Case/Control Ratio Using Genes from Brian's Cox PH Model Paper",
          outer = TRUE, cex = 1.3, font = 2, line = 1.5)
    mtext(sprintf("(Page %d of %d)", page, n_pages),
          outer = TRUE, cex = 1.0, line = 0.2)
  }
  
  dev.off()
  
  message(sprintf("\nTrue genes visualization saved to: exploratory_analysis_elastic_net_results/true_genes_brian_cox_ph_trends.pdf"))
  message(sprintf("  Total: %d true genes across %d pages", n_true_genes, n_pages))
  
} else {
  message("No true genes found in the dataset. Skipping visualization.")
}

