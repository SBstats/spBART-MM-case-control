################################################################################
# 00_prepare_analysis_base.R
#
# Builds the pooled analysis data ONCE for the analyses that need the
# real 5hmC matrix, without changing 03_real_data_analysis.R:
#   - real_data_analysis/analysis_base/analysis_base.RDS
#       development (n = 500) / validation (n = 369) split, VST genes, study
#       indicator, covariates, and the probit + GMM screened genes
#       (used by 05_BART_application.R)
#   - simulation/data_semisynthetic/semisynthetic_base.RDS
#       standardized real X (115 screened genes; 2,000 random genes) and real Z
#       for the semi-synthetic simulations (studies 07, 09 and 11)
#
# Cohort, BMI (at diagnosis; cutoff 25), gene filtering, VST, the split
# (seed 789) and the screening (age + sex probit, GMM 2-4 components) are
# exactly those of 03_real_data_analysis.R.
################################################################################

script_path <- tryCatch({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE)))  # Rscript encodes spaces as ~+~
  } else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    getwd()
  }
}, error = function(e) getwd())
setwd(script_path)
spbart_root <- dirname(script_path)
source(file.path(spbart_root, "R", "spBART_sampler.R"))   # spbart_user_lib()
.libPaths(unique(c(spbart_user_lib(), .libPaths())))

suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(DESeq2); library(mclust)
})

#-------------------------------------------------------------------------------
# Cohort construction, filtering and VST (identical to 03 / 04)
#-------------------------------------------------------------------------------
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
# BMI cutoff >= 25 (binaryBMInewCutoff)
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
# BMI cutoff >= 25 (binaryBMInewCutoff)
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
cohort_selected$study_indicator <- 1        # S = 1: UCMM Study (manuscript Eq. 1)

caseControl_selected <- caseControl_workdf %>%
  dplyr::select(PooledData_ID, race_binary, bmi_binary, sex_binary, age, mm_status)
caseControl_selected$study_indicator <- 0   # S = 0: BC Case-control Study

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

#-------------------------------------------------------------------------------
# Development / validation split (identical to 03: seed 789, test sampled first)
#-------------------------------------------------------------------------------
set.seed(789)
N_total <- nrow(pooled_metadata)
N_CV    <- 500
N_test  <- N_total - N_CV
outcome_1_ids <- pooled_metadata$PooledData_ID[pooled_metadata$outcome == 1]
outcome_0_ids <- pooled_metadata$PooledData_ID[pooled_metadata$outcome == 0]
prop_cases <- mean(pooled_metadata$outcome == 1)
n_test_cases <- round(N_test * prop_cases)
n_test_controls <- N_test - n_test_cases
test_ids <- c(sample(outcome_1_ids, n_test_cases, replace = FALSE),
              sample(outcome_0_ids, n_test_controls, replace = FALSE))
cv_ids <- setdiff(pooled_metadata$PooledData_ID, test_ids)
cv_indices <- which(pooled_metadata$PooledData_ID %in% cv_ids)
test_indices <- which(pooled_metadata$PooledData_ID %in% test_ids)

X_genes_cv   <- t(pooled_genebody_data_filtered_normalized[, cv_indices])
X_genes_test <- t(pooled_genebody_data_filtered_normalized[, test_indices])
Z_cv       <- as.matrix(pooled_metadata[cv_indices, c("sex_binary", "age")])
Z_cv_full  <- as.matrix(pooled_metadata[cv_indices, c("race_binary", "bmi_binary", "sex_binary", "age")])
Z_test_full <- as.matrix(pooled_metadata[test_indices, c("race_binary", "bmi_binary", "sex_binary", "age")])
Y_cv   <- pooled_metadata$outcome[cv_indices]
Y_test <- pooled_metadata$outcome[test_indices]
cat("Split: development", length(Y_cv), "(", sum(Y_cv), "cases ), validation", length(Y_test),
    "(", sum(Y_test), "cases )\n")

#-------------------------------------------------------------------------------
# Screening (identical rule to 03): probit per gene adjusted for sex and age,
# p < 0.05, then the strongest GMM cluster of |coefficients|
#-------------------------------------------------------------------------------
cat("Screening", ncol(X_genes_cv), "genes ...\n")
screen_one <- function(x) {
  df <- data.frame(y = Y_cv, gene = x, Z_cv)
  fit <- tryCatch(suppressWarnings(glm(y ~ ., data = df, family = binomial(link = "probit"),
                                       control = glm.control(maxit = 100))), error = function(e) NULL)
  if (is.null(fit)) return(c(coef = 0, p = 1))
  cf <- summary(fit)$coefficients
  if (!"gene" %in% rownames(cf)) return(c(coef = 0, p = 1))
  c(coef = cf["gene", "Estimate"], p = cf["gene", "Pr(>|z|)"])
}
scr <- t(apply(X_genes_cv, 2, screen_one))
sig <- rownames(scr)[scr[, "p"] < 0.05]
abs_coef <- abs(scr[sig, "coef"])
gmm <- Mclust(abs_coef, G = 2:4, modelNames = "V", verbose = FALSE)
cl_means <- sapply(seq_len(gmm$G), function(k) mean(abs_coef[gmm$classification == k]))
screened_genes <- sig[gmm$classification == which.max(cl_means)]
cat(sprintf("Screening: %d genes with p < 0.05; BIC chose %d clusters; %d genes kept\n",
            length(sig), gmm$G, length(screened_genes)))

#-------------------------------------------------------------------------------
# Save the analysis base
#-------------------------------------------------------------------------------
dir.create("analysis_base", showWarnings = FALSE)
saveRDS(list(X_genes_cv = X_genes_cv, X_genes_test = X_genes_test,
             study_cv = pooled_metadata$study_indicator[cv_indices],
             study_test = pooled_metadata$study_indicator[test_indices],
             Z_cv = Z_cv, Z_cv_full = Z_cv_full, Z_test_full = Z_test_full,
             Y_cv = Y_cv, Y_test = Y_test, screened_genes = screened_genes,
             screening = scr, cv_ids = cv_ids, test_ids = test_ids),
        file.path("analysis_base", "analysis_base.RDS"))

#-------------------------------------------------------------------------------
# Semi-synthetic base: standardized real X and real Z (age standardized)
#-------------------------------------------------------------------------------
set.seed(2026)
random_genes <- sample(colnames(X_genes_cv), 2000)
stdz <- function(A, ref) sweep(sweep(A, 2, colMeans(ref)), 2, apply(ref, 2, sd), "/")
mk_Z <- function(Zf) cbind(age = (Zf[, "age"] - mean(Z_cv_full[, "age"])) / sd(Z_cv_full[, "age"]),
                           sex = Zf[, "sex_binary"], race = Zf[, "race_binary"], bmi = Zf[, "bmi_binary"])
semisyn <- list(
  X_dev = list(screened   = stdz(X_genes_cv[, screened_genes], X_genes_cv[, screened_genes]),
               random2000 = stdz(X_genes_cv[, random_genes],   X_genes_cv[, random_genes])),
  X_val = list(screened   = stdz(X_genes_test[, screened_genes], X_genes_cv[, screened_genes]),
               random2000 = stdz(X_genes_test[, random_genes],   X_genes_cv[, random_genes])),
  Z_dev = mk_Z(Z_cv_full), Z_val = mk_Z(Z_test_full),
  note = "Standardized with development-set means/s.d.; Z = (age std, sex, race, BMI >= 25)")
out_dir <- file.path(spbart_root, "simulation", "data_semisynthetic")
dir.create(out_dir, showWarnings = FALSE)
saveRDS(semisyn, file.path(out_dir, "semisynthetic_base.RDS"))
cat("Saved analysis_base/analysis_base.RDS and simulation/data_semisynthetic/semisynthetic_base.RDS\n")
