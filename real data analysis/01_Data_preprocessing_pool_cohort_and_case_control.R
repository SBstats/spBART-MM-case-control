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



# Load required package for probit regression
if (!require("MASS")) install.packages("MASS")
library(MASS)


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
                           mm_types = patient_metadata_cohort_subset$dx_errc
                           )
#iss_stage= patient_metadata_cohort_subset$iss_derived,
#education = patient_metadata_cohort_subset$education_qx,
#ldh_level = as.numeric(patient_metadata_cohort_subset$`LACTIC DEHYDROGENASE_value_(U/L)`),
#estimated_gfr = patient_metadata_cohort_subset$egfr_ckdepi,
#asct_abstracted = patient_metadata_cohort_subset$asct_abstracted
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



# 3) code bmi to binary non-overweight (0) vs overweight/obese (1)
bmi_binary = ifelse(cohort_workdf$bmi>= 25,1,0)
cohort_workdf$bmi_binary = bmi_binary


# 4) code sex to binary F (0) vs M (1)
sex_binary = ifelse(cohort_workdf$sex== "M",1,0)
#sex_binary[is.na(sex_binary)] = 0. #removing NAs with 0
cohort_workdf$sex_binary = sex_binary


cohort_workdf$iss_stage <- patient_metadata_cohort_subset$iss_derived[match(cohort_workdf$errcid, patient_metadata_cohort_subset$errcid)]
cohort_workdf$education <- patient_metadata_cohort_subset$education_qx[match(cohort_workdf$errcid, patient_metadata_cohort_subset$errcid)]
cohort_workdf$ldh_level <- as.numeric(patient_metadata_cohort_subset$`LACTIC DEHYDROGENASE_value_(U/L)`[match(cohort_workdf$errcid, patient_metadata_cohort_subset$errcid)])
cohort_workdf$estimated_gfr <- patient_metadata_cohort_subset$egfr_ckdepi[match(cohort_workdf$errcid, patient_metadata_cohort_subset$errcid)]
cohort_workdf$asct_abstracted <- patient_metadata_cohort_subset$asct_abstracted[match(cohort_workdf$errcid, patient_metadata_cohort_subset$errcid)]
cohort_workdf$kappa_lambda_ratio <- patient_metadata_cohort_subset$`KAPPA/LAMBDA RATIO, S_value_(NA)`[match(cohort_workdf$errcid, patient_metadata_cohort_subset$errcid)]
cohort_workdf$smoking <- patient_metadata_cohort_subset$tobacco_status_qx[match(cohort_workdf$errcid, patient_metadata_cohort_subset$errcid)]



# 5) code education 
cohort_workdf$education_3cat <- ifelse(
  cohort_workdf$education %in% c("6 years (grade school)", 
                                 "10-12 years (high schol)",
                                 "Technical training (beyond high school)"),
  "Below college",
  ifelse(cohort_workdf$education %in% c("13-15 years (some college)",
                                        "16 years (completed college)"),
         "Some college or completed college",
         ifelse(cohort_workdf$education == "More than 16 years (graduate or professional degree)",
                "Graduate or professional degree",
                NA_character_))
)




# 6) code Smoking to two categories
table(patient_metadata_cohort_subset$Q11_CigaretteUse)


table(patient_metadata_cohort_subset$tobacco_status_qx)

smoking_2cat_cohort <- case_when(
  cohort_workdf$smoking %in% c("Current Smoker", 
                               "Former Smoker") ~ "Yes",
  cohort_workdf$smoking == "Never Tobacco User" ~ "No",
  TRUE ~ NA_character_
)

cohort_workdf$smoking_2cat = smoking_2cat_cohort



# 7) code LDH level
ldh_binary =  ifelse(is.na(cohort_workdf$ldh_level), NA,
                                   ifelse(cohort_workdf$ldh_level >= 240, 1, 0)) #>= 240 (elevated) coded 1 and < 240 (normal) coded 0
cohort_workdf$ldh_binary = ldh_binary


# 8) code estimated GFR (Estimated GFR (mL/min/1.73 m2), median) GFR = glomerular filtration rate
estimated_gfr_binary =  ifelse(is.na(cohort_workdf$estimated_gfr), NA,
                     ifelse(cohort_workdf$estimated_gfr >= 60, 1, 0)) #>= 60 (elevated) coded 1 and < 240 (normal) coded 0
cohort_workdf$estimated_gfr_binary = estimated_gfr_binary






# 9) Code autologous stem cell transplantation (ASCT)
asct_3cat <- ifelse(
  patient_metadata_cohort_subset$asct_abstracted == "yes", 1,
  ifelse(patient_metadata_cohort_subset$asct_abstracted == "no", 0, NA_character_)
)

patient_metadata_cohort_subset$asct_3cat = as.numeric(asct_3cat)


# 10) Code Serum free light chain (kappa/lambda) ratio

kappa_lambda_ratio_raw <- cohort_workdf$kappa_lambda_ratio

# Extract the first number from each string
kappa_lambda_ratio_parsed <- as.numeric(str_extract(kappa_lambda_ratio_raw, "[0-9]+\\.?[0-9]*"))

# Now categorize
cohort_workdf$kappa_lambda_ratio_3cat <- case_when(
  is.na(kappa_lambda_ratio_parsed) ~ NA_character_,
  kappa_lambda_ratio_parsed < 0.26 ~ "Low (<0.26)",
  kappa_lambda_ratio_parsed >= 0.26 & kappa_lambda_ratio_parsed <= 1.65 ~ "Normal (0.26-1.65)",
  kappa_lambda_ratio_parsed > 1.65 ~ "High (>1.65)"
)



# 11) Define a new identifier for pooled analysis
pooled_id_cohort = c(1:nrow(cohort_workdf))
cohort_workdf$PooledData_ID = pooled_id_cohort


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
# patient_metadata_case_control_v2 <- read.xlsx(
#  file = "data/MM_Questionnaire_BC v2 add smoking for BART project.xlsx",        # v2 version of the data contains Smoking variable
#  sheetIndex = 1,
#  password = "MMVan2019"
# )

# Save as RDS file (R's native format, no password needed for future use)
#saveRDS(patient_metadata_case_control_v2, "data/MM_Questionnaire_BC_for_BART_projectsV2.RDS")

# Load patient metadata from RDS file (no password needed) (861 x 41 df)
patient_metadata_case_control <- readRDS("data/MM_Questionnaire_BC_for_BART_projectsV2.RDS")

message("Table of Cases and Control before pre-processing:")
table(patient_metadata_case_control$MM_STATUS)


#head(patient_metadata_case_control)

# Load key connecting Study.ID (used to identify clinical covariates) and Assigned.ID (used to identify genebody)
Canada_case_control_sample_key <- as.data.frame(read.csv("data/Canada_case_control_sample_key.csv"))
#(734 x 2 df)
Canada_case_control_key = data.frame(Study.ID = Canada_case_control_sample_key$Study.ID,
                                     Assigned.ID = Canada_case_control_sample_key$Assigned.ID)




# Subset patient metadata based on Study.ID
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







###########################################################################################################
##################CLinical variables were only collected for cases (MM patients)###########################
# Load patient clinical data (password-protected xlsx file)
#library(xlsx)

# Read password-protected Excel file using xlsx package
# patient_clinicaldata_case_control <- read.xlsx(
#  file = "data/MM_ClinicalData_BC_eligible",
#  sheetIndex = 1,
#  password = "MMVan2019"
# )

# Save as RDS file (R's native format, no password needed for future use)
#saveRDS(patient_metadata_case_control, "data/MM_ClinicalData_BC_eligible.RDS")

# Load patient metadata from RDS file (no password needed) (861 x 41 df)
patient_clinicaldata_case_control <- readRDS("data/MM_ClinicalData_BC_eligible.RDS")     # clinical data version contains CASES only
###########################################################################################################


message("Table of Cases and Control before pre-processing:")
table(patient_clinicaldata_case_control$MM_STATUS)


#head(patient_clinicaldata_case_control)





# Subset patient clinical based on Study.ID
patient_clinicaldata_case_control_subset <- patient_clinicaldata_case_control %>%
  filter(StudyID %in% Canada_case_control_key$Study.ID)














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

# 2) code bmi to binary non-overweight (0) vs overweight/obese (1)
bmi_binary_case_control = ifelse(caseControl_workdf$bmi>= 25,1,0)
caseControl_workdf$bmi_binary = bmi_binary_case_control


# 3) code sex to binary F (0) vs M (1)
sex_binary_case_control = ifelse(caseControl_workdf$sex== "M",1,0)
caseControl_workdf$sex_binary = sex_binary_case_control


# 4) Define a new identifier for pooled analysis
pooled_id_caseControl = c( (nrow(cohort_workdf)+1): (nrow(cohort_workdf)+nrow(caseControl_workdf)) )
caseControl_workdf$PooledData_ID = pooled_id_caseControl


# 5) code education to three categories
caseControl_workdf$education <- patient_metadata_case_control_subset$Educ[match(caseControl_workdf$StudyID, patient_metadata_case_control_subset$StudyID)]
education_3cat_case_control <- case_when(
  caseControl_workdf$education %in% c(1, 2, 3, 4, 5) ~ "Below college",
  caseControl_workdf$education %in% c(6, 7) ~ "Some college or completed college",
  caseControl_workdf$education == 8 ~ "Graduate or professional degree",
  TRUE ~ NA_character_
)

caseControl_workdf$education_3cat = education_3cat_case_control


message("The dimension of caseControl_workdf after pre-processing is: ", nrow(caseControl_workdf), "x", ncol(caseControl_workdf), ". \n")





# 6) code Smoking to two categories
table(patient_metadata_case_control_subset$Smk)
caseControl_workdf$smoking <- patient_metadata_case_control_subset$Smk[match(caseControl_workdf$StudyID, patient_metadata_case_control_subset$StudyID)]
smoking_2cat_case_control <- case_when(
  caseControl_workdf$smoking == 1 ~ "Yes",
  caseControl_workdf$smoking == 2 ~ "No",
  TRUE ~ NA_character_
)

caseControl_workdf$smoking_2cat = smoking_2cat_case_control


message("The dimension of caseControl_workdf after pre-processing is: ", nrow(caseControl_workdf), "x", ncol(caseControl_workdf), ". \n")



##############################################################################################
##############Clinical variables for CASES ONLY in the case-control analysis##################

caseControl_workdf_clinical_vars = data.frame(StudyID = patient_clinicaldata_case_control_subset$StudyID,
                                iss_stage = patient_clinicaldata_case_control_subset$ISS_Stage,
                                kappa_lambda_ratio = patient_clinicaldata_case_control_subset$dFLC..ratio.)


caseControl_workdf_clinical_vars <- caseControl_workdf_clinical_vars[
  caseControl_workdf_clinical_vars$StudyID %in% caseControl_workdf$StudyID[caseControl_workdf$mm_status == "CASE"], ]




# 1) code iss_stage to match cohort categories {1, 2, 3, Could not be calculated}
caseControl_workdf_clinical_vars$iss_stage <- case_when(
  caseControl_workdf_clinical_vars$iss_stage == "I" ~ "1",
  caseControl_workdf_clinical_vars$iss_stage == "II" ~ "2",
  caseControl_workdf_clinical_vars$iss_stage == "III" ~ "3",
  caseControl_workdf_clinical_vars$iss_stage %in% c("UNK", "") ~ "Could not be calculated",
  TRUE ~ NA_character_
)



# 2) Code Serum free light chain (kappa/lambda) ratio

kappa_lambda_ratio_raw_case_control <- caseControl_workdf_clinical_vars$kappa_lambda_ratio

# Extract the first number from each string
kappa_lambda_ratio_case_control_parsed <- as.numeric(str_extract(kappa_lambda_ratio_raw_case_control, "[0-9]+\\.?[0-9]*"))

# Now categorize
caseControl_workdf_clinical_vars$kappa_lambda_ratio_3cat <- case_when(
  is.na(kappa_lambda_ratio_case_control_parsed) ~ NA_character_,
  kappa_lambda_ratio_case_control_parsed < 0.26 ~ "Low (<0.26)",
  kappa_lambda_ratio_case_control_parsed >= 0.26 & kappa_lambda_ratio_case_control_parsed <= 1.65 ~ "Normal (0.26-1.65)",
  kappa_lambda_ratio_case_control_parsed > 1.65 ~ "High (>1.65)"
)


##############################################################################################



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
  dplyr::select(PooledData_ID, race_binary, bmi_binary, sex_binary, age, mm_status, smoking_2cat, education_3cat)

# Select variables from caseControl_workdf
caseControl_selected <- caseControl_workdf %>%
  dplyr::select(PooledData_ID, race_binary, bmi_binary, sex_binary, age, mm_status, smoking_2cat, education_3cat)

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
# For normalization ONLY, we don't need real experimental conditions - just placeholder data

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
# 3) Generate Table 1 for Manuscript: Demographic and clinical characteristics (MM Cases vs Controls)
#    This table matches the format in the manuscript (sn-article.tex)
################################################################################################################

message("\n")
message("================================================================================")
message("GENERATING TABLE 1 FOR MANUSCRIPT: Demographic and clinical characteristics")
message("================================================================================\n")

# Add study indicator to pooled_metadata
pooled_metadata$study <- ifelse(pooled_metadata$PooledData_ID <= nrow(cohort_selected),
                                "UChicago", "Canada")

# Separate cases and controls
cases <- pooled_metadata[pooled_metadata$outcome == 1, ]
controls <- pooled_metadata[pooled_metadata$outcome == 0, ]

n_cases <- nrow(cases)
n_controls <- nrow(controls)

message("Sample sizes:")
message("  MM Cases: ", n_cases)
message("  Controls: ", n_controls)
message("  Total: ", n_cases + n_controls)

# ------------------------------------------------------------------------------
# Helper functions for Table 1 (Manuscript version)
# ------------------------------------------------------------------------------

# Format count and percentage
format_n_pct_manuscript <- function(count, total) {
  pct <- round(100 * count / total, 1)
  sprintf("%d (%.1f%%)", count, pct)
}

# Format mean and SD
format_mean_sd <- function(x) {
  sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

# Format p-value
format_pval <- function(p) {
  if (p < 0.001) {
    return("<0.001")
  } else {
    return(sprintf("%.3f", p))
  }
}

# ------------------------------------------------------------------------------
# Compute statistics for each characteristic
# ------------------------------------------------------------------------------

# 1. Sample size (N)
n_case_str <- as.character(n_cases)
n_control_str <- as.character(n_controls)
pval_n <- "--"

# 2. Age: mean (SD) with t-test
age_case <- format_mean_sd(cases$age)
age_control <- format_mean_sd(controls$age)
age_ttest <- t.test(cases$age, controls$age)
pval_age <- format_pval(age_ttest$p.value)

# 3. Sex: n (%) with chi-square test
n_male_case <- sum(cases$sex_binary == 1)
n_male_control <- sum(controls$sex_binary == 1)
n_female_case <- sum(cases$sex_binary == 0)
n_female_control <- sum(controls$sex_binary == 0)

sex_table <- matrix(c(n_male_case, n_male_control, n_female_case, n_female_control),
                    nrow = 2, byrow = TRUE)
sex_chisq <- chisq.test(sex_table)
pval_sex <- format_pval(sex_chisq$p.value)

male_case_str <- format_n_pct_manuscript(n_male_case, n_cases)
male_control_str <- format_n_pct_manuscript(n_male_control, n_controls)
female_case_str <- format_n_pct_manuscript(n_female_case, n_cases)
female_control_str <- format_n_pct_manuscript(n_female_control, n_controls)

# 4. Race: n (%) with chi-square test
n_white_case <- sum(cases$race_binary == 1)
n_white_control <- sum(controls$race_binary == 1)
n_aa_case <- sum(cases$race_binary == 0)
n_aa_control <- sum(controls$race_binary == 0)

race_table <- matrix(c(n_white_case, n_white_control, n_aa_case, n_aa_control),
                     nrow = 2, byrow = TRUE)
race_chisq <- chisq.test(race_table)
pval_race <- format_pval(race_chisq$p.value)

white_case_str <- format_n_pct_manuscript(n_white_case, n_cases)
white_control_str <- format_n_pct_manuscript(n_white_control, n_controls)
aa_case_str <- format_n_pct_manuscript(n_aa_case, n_cases)
aa_control_str <- format_n_pct_manuscript(n_aa_control, n_controls)

# 5. BMI category: n (%) with chi-square test
n_bmi_ge30_case <- sum(cases$bmi_binary == 1)
n_bmi_ge30_control <- sum(controls$bmi_binary == 1)
n_bmi_lt30_case <- sum(cases$bmi_binary == 0)
n_bmi_lt30_control <- sum(controls$bmi_binary == 0)

bmi_table <- matrix(c(n_bmi_lt30_case, n_bmi_lt30_control, n_bmi_ge30_case, n_bmi_ge30_control),
                    nrow = 2, byrow = TRUE)
bmi_chisq <- chisq.test(bmi_table)
pval_bmi <- format_pval(bmi_chisq$p.value)

bmi_lt30_case_str <- format_n_pct_manuscript(n_bmi_lt30_case, n_cases)
bmi_lt30_control_str <- format_n_pct_manuscript(n_bmi_lt30_control, n_controls)
bmi_ge30_case_str <- format_n_pct_manuscript(n_bmi_ge30_case, n_cases)
bmi_ge30_control_str <- format_n_pct_manuscript(n_bmi_ge30_control, n_controls)

# 6. Study: n (%) - no p-value (design variable)
n_uchicago_case <- sum(cases$study == "UChicago")
n_uchicago_control <- sum(controls$study == "UChicago")
n_canada_case <- sum(cases$study == "Canada")
n_canada_control <- sum(controls$study == "Canada")

uchicago_case_str <- format_n_pct_manuscript(n_uchicago_case, n_cases)
uchicago_control_str <- format_n_pct_manuscript(n_uchicago_control, n_controls)
canada_case_str <- format_n_pct_manuscript(n_canada_case, n_cases)
canada_control_str <- format_n_pct_manuscript(n_canada_control, n_controls)

# 7. Smoking: n (%) with chi-square test (excluding NAs from test)
n_smoking_yes_case <- sum(cases$smoking_2cat == "Yes", na.rm = TRUE)
n_smoking_yes_control <- sum(controls$smoking_2cat == "Yes", na.rm = TRUE)
n_smoking_no_case <- sum(cases$smoking_2cat == "No", na.rm = TRUE)
n_smoking_no_control <- sum(controls$smoking_2cat == "No", na.rm = TRUE)
n_smoking_na_case <- sum(is.na(cases$smoking_2cat))
n_smoking_na_control <- sum(is.na(controls$smoking_2cat))

smoking_table <- matrix(c(n_smoking_yes_case, n_smoking_yes_control,
                           n_smoking_no_case, n_smoking_no_control),
                         nrow = 2, byrow = TRUE)
smoking_chisq <- chisq.test(smoking_table)
pval_smoking <- format_pval(smoking_chisq$p.value)

smoking_yes_case_str <- format_n_pct_manuscript(n_smoking_yes_case, n_cases)
smoking_yes_control_str <- format_n_pct_manuscript(n_smoking_yes_control, n_controls)
smoking_no_case_str <- format_n_pct_manuscript(n_smoking_no_case, n_cases)
smoking_no_control_str <- format_n_pct_manuscript(n_smoking_no_control, n_controls)
smoking_na_case_str <- format_n_pct_manuscript(n_smoking_na_case, n_cases)
smoking_na_control_str <- format_n_pct_manuscript(n_smoking_na_control, n_controls)

# 8. Education: n (%) with chi-square test (excluding NAs from test)
n_edu_below_case <- sum(cases$education_3cat == "Below college", na.rm = TRUE)
n_edu_below_control <- sum(controls$education_3cat == "Below college", na.rm = TRUE)
n_edu_some_case <- sum(cases$education_3cat == "Some college or completed college", na.rm = TRUE)
n_edu_some_control <- sum(controls$education_3cat == "Some college or completed college", na.rm = TRUE)
n_edu_grad_case <- sum(cases$education_3cat == "Graduate or professional degree", na.rm = TRUE)
n_edu_grad_control <- sum(controls$education_3cat == "Graduate or professional degree", na.rm = TRUE)
n_edu_na_case <- sum(is.na(cases$education_3cat))
n_edu_na_control <- sum(is.na(controls$education_3cat))

edu_table <- matrix(c(n_edu_below_case, n_edu_below_control,
                       n_edu_some_case, n_edu_some_control,
                       n_edu_grad_case, n_edu_grad_control),
                     nrow = 3, byrow = TRUE)
edu_chisq <- chisq.test(edu_table)
pval_edu <- format_pval(edu_chisq$p.value)

edu_below_case_str <- format_n_pct_manuscript(n_edu_below_case, n_cases)
edu_below_control_str <- format_n_pct_manuscript(n_edu_below_control, n_controls)
edu_some_case_str <- format_n_pct_manuscript(n_edu_some_case, n_cases)
edu_some_control_str <- format_n_pct_manuscript(n_edu_some_control, n_controls)
edu_grad_case_str <- format_n_pct_manuscript(n_edu_grad_case, n_cases)
edu_grad_control_str <- format_n_pct_manuscript(n_edu_grad_control, n_controls)
edu_na_case_str <- format_n_pct_manuscript(n_edu_na_case, n_cases)
edu_na_control_str <- format_n_pct_manuscript(n_edu_na_control, n_controls)

# ------------------------------------------------------------------------------
# Create Table 1 data frame (Manuscript version)
# ------------------------------------------------------------------------------

table1_manuscript <- data.frame(
  Characteristic = c(
    "N",
    "Age, mean (SD)",
    "Sex, n (%)",
    "  Male",
    "  Female",
    "Race, n (%)",
    "  White",
    "  African American",
    "BMI category, n (%)",
    "  <25 kg/m2",
    "  >=25 kg/m2",
    "Smoking history, n (%)",
    "  Ever smoker",
    "  Never smoker",
    "  Missing",
    "Education, n (%)",
    "  Below college",
    "  Some college or completed college",
    "  Graduate or professional degree",
    "  Missing",
    "Study, n (%)",
    "  UChicago",
    "  Canada"
  ),
  MM_Cases = c(
    n_case_str,
    age_case,
    "",
    male_case_str,
    female_case_str,
    "",
    white_case_str,
    aa_case_str,
    "",
    bmi_lt30_case_str,
    bmi_ge30_case_str,
    "",
    smoking_yes_case_str,
    smoking_no_case_str,
    smoking_na_case_str,
    "",
    edu_below_case_str,
    edu_some_case_str,
    edu_grad_case_str,
    edu_na_case_str,
    "",
    uchicago_case_str,
    canada_case_str
  ),
  Controls = c(
    n_control_str,
    age_control,
    "",
    male_control_str,
    female_control_str,
    "",
    white_control_str,
    aa_control_str,
    "",
    bmi_lt30_control_str,
    bmi_ge30_control_str,
    "",
    smoking_yes_control_str,
    smoking_no_control_str,
    smoking_na_control_str,
    "",
    edu_below_control_str,
    edu_some_control_str,
    edu_grad_control_str,
    edu_na_control_str,
    "",
    uchicago_control_str,
    canada_control_str
  ),
  P_value = c(
    pval_n,
    pval_age,
    pval_sex,
    "",
    "",
    pval_race,
    "",
    "",
    pval_bmi,
    "",
    "",
    pval_smoking,
    "",
    "",
    "",
    pval_edu,
    "",
    "",
    "",
    "",
    "--",
    "",
    ""
  ),
  stringsAsFactors = FALSE
)

# Print the table
message("\n=== Table 1 (Manuscript): Demographic and clinical characteristics ===\n")
print(table1_manuscript, row.names = FALSE)

# ------------------------------------------------------------------------------
# Generate LaTeX code for the manuscript
# ------------------------------------------------------------------------------

message("\n")
message("================================================================================")
message("LaTeX code for Table 1 (copy to sn-article.tex):")
message("================================================================================\n")

latex_table <- paste0(
  "\\begin{table}[h]\n",
  "\\caption{Demographic and clinical characteristics of the study population}\\label{tab:demographics}\n",
  "\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}lccc}\n",
  "\\toprule\n",
  "\\textbf{Characteristic} & \\textbf{MM Cases} & \\textbf{Controls} & \\textbf{P-value} \\\\\n",
  "\\midrule\n",
  "N & ", n_case_str, " & ", n_control_str, " & -- \\\\\n",
  "Age, mean (SD) & ", age_case, " & ", age_control, " & ", pval_age, " \\\\\n",
  "Sex, n (\\%) & & & ", pval_sex, " \\\\\n",
  "\\quad Male & ", male_case_str, " & ", male_control_str, " & \\\\\n",
  "\\quad Female & ", female_case_str, " & ", female_control_str, " & \\\\\n",
  "Race, n (\\%) & & & ", pval_race, " \\\\\n",
  "\\quad White & ", white_case_str, " & ", white_control_str, " & \\\\\n",
  "\\quad African American & ", aa_case_str, " & ", aa_control_str, " & \\\\\n",
  "BMI category, n (\\%) & & & ", pval_bmi, " \\\\\n",
  "\\quad $<$25 kg/m$^2$ & ", bmi_lt30_case_str, " & ", bmi_lt30_control_str, " & \\\\\n",
  "\\quad $\\geq$25 kg/m$^2$ & ", bmi_ge30_case_str, " & ", bmi_ge30_control_str, " & \\\\\n",
  "Smoking history, n (\\%)$^{\\dagger}$ & & & ", pval_smoking, " \\\\\n",
  "\\quad Ever smoker & ", smoking_yes_case_str, " & ", smoking_yes_control_str, " & \\\\\n",
  "\\quad Never smoker & ", smoking_no_case_str, " & ", smoking_no_control_str, " & \\\\\n",
  "\\quad Missing & ", smoking_na_case_str, " & ", smoking_na_control_str, " & \\\\\n",
  "Education, n (\\%)$^{\\dagger}$ & & & ", pval_edu, " \\\\\n",
  "\\quad Below college & ", edu_below_case_str, " & ", edu_below_control_str, " & \\\\\n",
  "\\quad Some college or completed college & ", edu_some_case_str, " & ", edu_some_control_str, " & \\\\\n",
  "\\quad Graduate or professional degree & ", edu_grad_case_str, " & ", edu_grad_control_str, " & \\\\\n",
  "\\quad Missing & ", edu_na_case_str, " & ", edu_na_control_str, " & \\\\\n",
  "Study, n (\\%) & & & -- \\\\\n",
  "\\quad UChicago & ", uchicago_case_str, " & ", uchicago_control_str, " & \\\\\n",
  "\\quad Canada & ", canada_case_str, " & ", canada_control_str, " & \\\\\n",
  "\\botrule\n",
  "\\end{tabular*}\n",
  "\\footnotetext{SD, standard deviation; BMI, body mass index. P-values calculated using t-test for continuous variables and chi-square test for categorical variables. $^{\\dagger}$Percentages are calculated using the total N as the denominator and may not sum to 100\\% due to missing data. Chi-square tests for smoking and education were performed on non-missing observations only.}\n",
  "\\end{table}\n"
)

cat(latex_table)

# Create output directory if it doesn't exist
if (!dir.exists("output")) {
  dir.create("output")
}

# Save LaTeX table to file
writeLines(latex_table, "output/table1_demographics.tex")
message("\nLaTeX table saved to: output/table1_demographics.tex")

# ------------------------------------------------------------------------------
# Additional summary statistics for the manuscript
# ------------------------------------------------------------------------------

message("\n")
message("================================================================================")
message("Additional statistics for manuscript text:")
message("================================================================================\n")

message("Total pooled sample size: ", nrow(pooled_metadata))
message("  MM Cases: ", n_cases, " (", round(100*n_cases/nrow(pooled_metadata), 1), "%)")
message("  Controls: ", n_controls, " (", round(100*n_controls/nrow(pooled_metadata), 1), "%)")
message("")
message("By study:")
message("  UChicago: ", sum(pooled_metadata$study == "UChicago"),
        " (Cases: ", n_uchicago_case, ", Controls: ", n_uchicago_control, ")")
message("  Canada: ", sum(pooled_metadata$study == "Canada"),
        " (Cases: ", n_canada_case, ", Controls: ", n_canada_control, ")")
message("")
message("Development set size: ", round(0.576 * nrow(pooled_metadata)))
message("Validation set size: ", nrow(pooled_metadata) - round(0.576 * nrow(pooled_metadata)))




################################################################################################################
# 3b) Generate Table 1 (Original version): Summary of clinical covariates by study
################################################################################################################

# Helper function to create formatted count and percentage
format_n_pct <- function(count, total) {
  pct <- round(100 * count / total, 1)
  sprintf("%d (%.1f\\%%)", count, pct)
}

# Function to generate summary statistics for a dataset
generate_summary <- function(data, dataset_name) {
  n_total <- nrow(data)
  
  # Participant type
  n_case <- sum(data$outcome == 1)
  n_control <- sum(data$outcome == 0)
  
  # Sex
  n_male <- sum(data$sex_binary == 1)
  n_female <- sum(data$sex_binary == 0)
  
  # Age groups
  n_age_lt60 <- sum(data$age < 60)
  n_age_60_70 <- sum(data$age >= 60 & data$age <= 70)
  n_age_gt70 <- sum(data$age > 70)
  
  # BMI
  n_bmi_ge30 <- sum(data$bmi_binary == 1)
  n_bmi_lt30 <- sum(data$bmi_binary == 0)
  
  # Race
  n_white <- sum(data$race_binary == 1)
  n_black <- sum(data$race_binary == 0)
  
  # Create summary data frame
  summary_df <- data.frame(
    Characteristic = c(
      "Participant type", "",
      "Sex", "",
      "Age (years)", "", "",
      "BMI", "",
      "Race", "",
      "Total"
    ),
    Category = c(
      "Case ($Y_i = 1$)", "Control ($Y_i = 0$)",
      "Male", "Female",
      "$<60$", "60--70", "$>70$",
      "$\\geq 25$", "$<25$",
      "White", "Black/African American",
      ""
    ),
    Value = c(
      format_n_pct(n_case, n_total),
      format_n_pct(n_control, n_total),
      format_n_pct(n_male, n_total),
      format_n_pct(n_female, n_total),
      format_n_pct(n_age_lt60, n_total),
      format_n_pct(n_age_60_70, n_total),
      format_n_pct(n_age_gt70, n_total),
      format_n_pct(n_bmi_ge30, n_total),
      format_n_pct(n_bmi_lt30, n_total),
      format_n_pct(n_white, n_total),
      format_n_pct(n_black, n_total),
      sprintf("%d (100.0\\%%)", n_total)
    ),
    stringsAsFactors = FALSE
  )
  
  return(summary_df)
}





# Generate summaries for each dataset
cohort_summary <- generate_summary(cohort_selected %>%
                                     dplyr::select(outcome = mm_status, sex_binary, age,
                                                   bmi_binary, race_binary) %>%
                                     mutate(outcome = ifelse(cohort_selected$mm_status == "CASE", 1, 0)),
                                   "UChicago")

# For cohort, all are cases, so we need to handle this specially
cohort_summary_special <- data.frame(
  Characteristic = c(
    "Participant type", "",
    "Sex", "",
    "Age (years)", "", "",
    "BMI", "",
    "Race", "",
    "Total"
  ),
  Category = c(
    "Case ($Y_i = 1$)", "Control ($Y_i = 0$)",
    "Male", "Female",
    "$<60$", "60--70", "$>70$",
    "$\\geq 25$", "$<25$",
    "White", "Black/African American",
    ""
  ),
  UChicago = c(
    format_n_pct(nrow(cohort_selected), nrow(cohort_selected)),
    format_n_pct(0, nrow(cohort_selected)),
    format_n_pct(sum(cohort_selected$sex_binary == 1), nrow(cohort_selected)),
    format_n_pct(sum(cohort_selected$sex_binary == 0), nrow(cohort_selected)),
    format_n_pct(sum(cohort_selected$age < 60), nrow(cohort_selected)),
    format_n_pct(sum(cohort_selected$age >= 60 & cohort_selected$age <= 70), nrow(cohort_selected)),
    format_n_pct(sum(cohort_selected$age > 70), nrow(cohort_selected)),
    format_n_pct(sum(cohort_selected$bmi_binary == 1), nrow(cohort_selected)),
    format_n_pct(sum(cohort_selected$bmi_binary == 0), nrow(cohort_selected)),
    format_n_pct(sum(cohort_selected$race_binary == 1), nrow(cohort_selected)),
    format_n_pct(sum(cohort_selected$race_binary == 0), nrow(cohort_selected)),
    sprintf("%d (100.0\\%%)", nrow(cohort_selected))
  ),
  stringsAsFactors = FALSE
)






# Case-control summary (mm_status is "CASE" or "CONTROL")
caseControl_summary_special <- data.frame(
  Canada = c(
    format_n_pct(sum(caseControl_selected$mm_status == "CASE"), nrow(caseControl_selected)),
    format_n_pct(sum(caseControl_selected$mm_status == "CONTROL"), nrow(caseControl_selected)),
    format_n_pct(sum(caseControl_selected$sex_binary == 1), nrow(caseControl_selected)),
    format_n_pct(sum(caseControl_selected$sex_binary == 0), nrow(caseControl_selected)),
    format_n_pct(sum(caseControl_selected$age < 60), nrow(caseControl_selected)),
    format_n_pct(sum(caseControl_selected$age >= 60 & caseControl_selected$age <= 70), nrow(caseControl_selected)),
    format_n_pct(sum(caseControl_selected$age > 70), nrow(caseControl_selected)),
    format_n_pct(sum(caseControl_selected$bmi_binary == 1), nrow(caseControl_selected)),
    format_n_pct(sum(caseControl_selected$bmi_binary == 0), nrow(caseControl_selected)),
    format_n_pct(sum(caseControl_selected$race_binary == 1), nrow(caseControl_selected)),
    format_n_pct(sum(caseControl_selected$race_binary == 0), nrow(caseControl_selected)),
    sprintf("%d (100.0\\%%)", nrow(caseControl_selected))
  ),
  stringsAsFactors = FALSE
)





# Pooled summary
pooled_summary_special <- data.frame(
  Pooled = c(
    format_n_pct(sum(pooled_metadata$outcome == 1), nrow(pooled_metadata)),
    format_n_pct(sum(pooled_metadata$outcome == 0), nrow(pooled_metadata)),
    format_n_pct(sum(pooled_metadata$sex_binary == 1), nrow(pooled_metadata)),
    format_n_pct(sum(pooled_metadata$sex_binary == 0), nrow(pooled_metadata)),
    format_n_pct(sum(pooled_metadata$age < 60), nrow(pooled_metadata)),
    format_n_pct(sum(pooled_metadata$age >= 60 & pooled_metadata$age <= 70), nrow(pooled_metadata)),
    format_n_pct(sum(pooled_metadata$age > 70), nrow(pooled_metadata)),
    format_n_pct(sum(pooled_metadata$bmi_binary == 1), nrow(pooled_metadata)),
    format_n_pct(sum(pooled_metadata$bmi_binary == 0), nrow(pooled_metadata)),
    format_n_pct(sum(pooled_metadata$race_binary == 1), nrow(pooled_metadata)),
    format_n_pct(sum(pooled_metadata$race_binary == 0), nrow(pooled_metadata)),
    sprintf("%d (100.0\\%%)", nrow(pooled_metadata))
  ),
  stringsAsFactors = FALSE
)

# Combine all summaries
table1_data <- cbind(cohort_summary_special, caseControl_summary_special, pooled_summary_special)

# Print the table
message("\n=== Table 1: Summary of Clinical Covariates ===")
print(table1_data)



################################################################################################################
# 4) Generate Table: UCMM Cases vs BC Cases (Canada MM Cases only, excluding controls)
#    Comparison of demographic and clinical characteristics between study sites
################################################################################################################

message("\n")
message("================================================================================")
message("GENERATING TABLE: UCMM Cases vs BC Cases (excluding controls)")
message("================================================================================\n")

# Separate UCMM cases (all from cohort) and BC cases (cases only from case-control)
ucmm_cases <- pooled_metadata[pooled_metadata$study == "UChicago" & pooled_metadata$outcome == 1, ]
bc_cases <- pooled_metadata[pooled_metadata$study == "Canada" & pooled_metadata$outcome == 1, ]

n_ucmm <- nrow(ucmm_cases)
n_bc <- nrow(bc_cases)

message("Sample sizes (Cases only):")
message("  UCMM Cases: ", n_ucmm)
message("  BC Cases: ", n_bc)
message("  Total Cases: ", n_ucmm + n_bc)

# ------------------------------------------------------------------------------
# Compute statistics for each characteristic (UCMM vs BC cases)
# ------------------------------------------------------------------------------

# 1. Sample size (N)
n_ucmm_str <- as.character(n_ucmm)
n_bc_str <- as.character(n_bc)
pval_n_cases <- "--"

# 2. Age: mean (SD) with t-test
age_ucmm <- format_mean_sd(ucmm_cases$age)
age_bc <- format_mean_sd(bc_cases$age)
age_ttest_cases <- t.test(ucmm_cases$age, bc_cases$age)
pval_age_cases <- format_pval(age_ttest_cases$p.value)

# 3. Sex: n (%) with chi-square test
n_male_ucmm <- sum(ucmm_cases$sex_binary == 1)
n_male_bc <- sum(bc_cases$sex_binary == 1)
n_female_ucmm <- sum(ucmm_cases$sex_binary == 0)
n_female_bc <- sum(bc_cases$sex_binary == 0)

sex_table_cases <- matrix(c(n_male_ucmm, n_male_bc, n_female_ucmm, n_female_bc),
                          nrow = 2, byrow = TRUE)
sex_chisq_cases <- chisq.test(sex_table_cases)
pval_sex_cases <- format_pval(sex_chisq_cases$p.value)

male_ucmm_str <- format_n_pct_manuscript(n_male_ucmm, n_ucmm)
male_bc_str <- format_n_pct_manuscript(n_male_bc, n_bc)
female_ucmm_str <- format_n_pct_manuscript(n_female_ucmm, n_ucmm)
female_bc_str <- format_n_pct_manuscript(n_female_bc, n_bc)

# 4. Race: n (%) with chi-square test
n_white_ucmm <- sum(ucmm_cases$race_binary == 1)
n_white_bc <- sum(bc_cases$race_binary == 1)
n_aa_ucmm <- sum(ucmm_cases$race_binary == 0)
n_aa_bc <- sum(bc_cases$race_binary == 0)

race_table_cases <- matrix(c(n_white_ucmm, n_white_bc, n_aa_ucmm, n_aa_bc),
                           nrow = 2, byrow = TRUE)
race_chisq_cases <- chisq.test(race_table_cases)
pval_race_cases <- format_pval(race_chisq_cases$p.value)

white_ucmm_str <- format_n_pct_manuscript(n_white_ucmm, n_ucmm)
white_bc_str <- format_n_pct_manuscript(n_white_bc, n_bc)
aa_ucmm_str <- format_n_pct_manuscript(n_aa_ucmm, n_ucmm)
aa_bc_str <- format_n_pct_manuscript(n_aa_bc, n_bc)

# 5. BMI category: n (%) with chi-square test
n_bmi_ge30_ucmm <- sum(ucmm_cases$bmi_binary == 1)
n_bmi_ge30_bc <- sum(bc_cases$bmi_binary == 1)
n_bmi_lt30_ucmm <- sum(ucmm_cases$bmi_binary == 0)
n_bmi_lt30_bc <- sum(bc_cases$bmi_binary == 0)

bmi_table_cases <- matrix(c(n_bmi_lt30_ucmm, n_bmi_lt30_bc, n_bmi_ge30_ucmm, n_bmi_ge30_bc),
                          nrow = 2, byrow = TRUE)
bmi_chisq_cases <- chisq.test(bmi_table_cases)
pval_bmi_cases <- format_pval(bmi_chisq_cases$p.value)

bmi_lt30_ucmm_str <- format_n_pct_manuscript(n_bmi_lt30_ucmm, n_ucmm)
bmi_lt30_bc_str <- format_n_pct_manuscript(n_bmi_lt30_bc, n_bc)
bmi_ge30_ucmm_str <- format_n_pct_manuscript(n_bmi_ge30_ucmm, n_ucmm)
bmi_ge30_bc_str <- format_n_pct_manuscript(n_bmi_ge30_bc, n_bc)

# 6. Smoking: n (%) with chi-square test (excluding NAs from test)
n_smoking_yes_ucmm <- sum(ucmm_cases$smoking_2cat == "Yes", na.rm = TRUE)
n_smoking_yes_bc <- sum(bc_cases$smoking_2cat == "Yes", na.rm = TRUE)
n_smoking_no_ucmm <- sum(ucmm_cases$smoking_2cat == "No", na.rm = TRUE)
n_smoking_no_bc <- sum(bc_cases$smoking_2cat == "No", na.rm = TRUE)
n_smoking_na_ucmm <- sum(is.na(ucmm_cases$smoking_2cat))
n_smoking_na_bc <- sum(is.na(bc_cases$smoking_2cat))

smoking_table_cases <- matrix(c(n_smoking_yes_ucmm, n_smoking_yes_bc,
                                 n_smoking_no_ucmm, n_smoking_no_bc),
                               nrow = 2, byrow = TRUE)
smoking_chisq_cases <- chisq.test(smoking_table_cases)
pval_smoking_cases <- format_pval(smoking_chisq_cases$p.value)

smoking_yes_ucmm_str <- format_n_pct_manuscript(n_smoking_yes_ucmm, n_ucmm)
smoking_yes_bc_str <- format_n_pct_manuscript(n_smoking_yes_bc, n_bc)
smoking_no_ucmm_str <- format_n_pct_manuscript(n_smoking_no_ucmm, n_ucmm)
smoking_no_bc_str <- format_n_pct_manuscript(n_smoking_no_bc, n_bc)
smoking_na_ucmm_str <- format_n_pct_manuscript(n_smoking_na_ucmm, n_ucmm)
smoking_na_bc_str <- format_n_pct_manuscript(n_smoking_na_bc, n_bc)

# 7. Education: n (%) with chi-square test (excluding NAs from test)
n_edu_below_ucmm <- sum(ucmm_cases$education_3cat == "Below college", na.rm = TRUE)
n_edu_below_bc <- sum(bc_cases$education_3cat == "Below college", na.rm = TRUE)
n_edu_some_ucmm <- sum(ucmm_cases$education_3cat == "Some college or completed college", na.rm = TRUE)
n_edu_some_bc <- sum(bc_cases$education_3cat == "Some college or completed college", na.rm = TRUE)
n_edu_grad_ucmm <- sum(ucmm_cases$education_3cat == "Graduate or professional degree", na.rm = TRUE)
n_edu_grad_bc <- sum(bc_cases$education_3cat == "Graduate or professional degree", na.rm = TRUE)
n_edu_na_ucmm <- sum(is.na(ucmm_cases$education_3cat))
n_edu_na_bc <- sum(is.na(bc_cases$education_3cat))

edu_table_cases <- matrix(c(n_edu_below_ucmm, n_edu_below_bc,
                              n_edu_some_ucmm, n_edu_some_bc,
                              n_edu_grad_ucmm, n_edu_grad_bc),
                            nrow = 3, byrow = TRUE)
edu_chisq_cases <- chisq.test(edu_table_cases)
pval_edu_cases <- format_pval(edu_chisq_cases$p.value)

edu_below_ucmm_str <- format_n_pct_manuscript(n_edu_below_ucmm, n_ucmm)
edu_below_bc_str <- format_n_pct_manuscript(n_edu_below_bc, n_bc)
edu_some_ucmm_str <- format_n_pct_manuscript(n_edu_some_ucmm, n_ucmm)
edu_some_bc_str <- format_n_pct_manuscript(n_edu_some_bc, n_bc)
edu_grad_ucmm_str <- format_n_pct_manuscript(n_edu_grad_ucmm, n_ucmm)
edu_grad_bc_str <- format_n_pct_manuscript(n_edu_grad_bc, n_bc)
edu_na_ucmm_str <- format_n_pct_manuscript(n_edu_na_ucmm, n_ucmm)
edu_na_bc_str <- format_n_pct_manuscript(n_edu_na_bc, n_bc)

# 8. ISS Stage: n (%) with chi-square test (excluding NAs from test)
# UCMM: from cohort_workdf$iss_stage
n_iss1_ucmm <- sum(cohort_workdf$iss_stage == "1", na.rm = TRUE)
n_iss2_ucmm <- sum(cohort_workdf$iss_stage == "2", na.rm = TRUE)
n_iss3_ucmm <- sum(cohort_workdf$iss_stage == "3", na.rm = TRUE)
n_iss_cnc_ucmm <- sum(cohort_workdf$iss_stage == "Could not be calculated", na.rm = TRUE)
n_iss_na_ucmm <- sum(is.na(cohort_workdf$iss_stage))

# BC: from caseControl_workdf_clinical_vars$iss_stage (already harmonized)
n_iss1_bc <- sum(caseControl_workdf_clinical_vars$iss_stage == "1", na.rm = TRUE)
n_iss2_bc <- sum(caseControl_workdf_clinical_vars$iss_stage == "2", na.rm = TRUE)
n_iss3_bc <- sum(caseControl_workdf_clinical_vars$iss_stage == "3", na.rm = TRUE)
n_iss_cnc_bc <- sum(caseControl_workdf_clinical_vars$iss_stage == "Could not be calculated", na.rm = TRUE)
n_iss_na_bc <- sum(is.na(caseControl_workdf_clinical_vars$iss_stage))

iss_table_cases <- matrix(c(n_iss1_ucmm, n_iss1_bc,
                             n_iss2_ucmm, n_iss2_bc,
                             n_iss3_ucmm, n_iss3_bc),
                           nrow = 3, byrow = TRUE)
iss_chisq_cases <- chisq.test(iss_table_cases)
pval_iss_cases <- format_pval(iss_chisq_cases$p.value)

iss1_ucmm_str <- format_n_pct_manuscript(n_iss1_ucmm, n_ucmm)
iss1_bc_str <- format_n_pct_manuscript(n_iss1_bc, n_bc)
iss2_ucmm_str <- format_n_pct_manuscript(n_iss2_ucmm, n_ucmm)
iss2_bc_str <- format_n_pct_manuscript(n_iss2_bc, n_bc)
iss3_ucmm_str <- format_n_pct_manuscript(n_iss3_ucmm, n_ucmm)
iss3_bc_str <- format_n_pct_manuscript(n_iss3_bc, n_bc)
iss_cnc_ucmm_str <- format_n_pct_manuscript(n_iss_cnc_ucmm, n_ucmm)
iss_cnc_bc_str <- format_n_pct_manuscript(n_iss_cnc_bc, n_bc)
iss_na_ucmm_str <- format_n_pct_manuscript(n_iss_na_ucmm, n_ucmm)
iss_na_bc_str <- format_n_pct_manuscript(n_iss_na_bc, n_bc)

# 7. Kappa/Lambda Ratio: n (%) with chi-square test (excluding NAs from test)
# UCMM: from cohort_workdf$kappa_lambda_ratio_3cat
n_kl_low_ucmm <- sum(cohort_workdf$kappa_lambda_ratio_3cat == "Low (<0.26)", na.rm = TRUE)
n_kl_normal_ucmm <- sum(cohort_workdf$kappa_lambda_ratio_3cat == "Normal (0.26-1.65)", na.rm = TRUE)
n_kl_high_ucmm <- sum(cohort_workdf$kappa_lambda_ratio_3cat == "High (>1.65)", na.rm = TRUE)
n_kl_na_ucmm <- sum(is.na(cohort_workdf$kappa_lambda_ratio_3cat))

# BC: from caseControl_workdf_clinical_vars$kappa_lambda_ratio_3cat
n_kl_low_bc <- sum(caseControl_workdf_clinical_vars$kappa_lambda_ratio_3cat == "Low (<0.26)", na.rm = TRUE)
n_kl_normal_bc <- sum(caseControl_workdf_clinical_vars$kappa_lambda_ratio_3cat == "Normal (0.26-1.65)", na.rm = TRUE)
n_kl_high_bc <- sum(caseControl_workdf_clinical_vars$kappa_lambda_ratio_3cat == "High (>1.65)", na.rm = TRUE)
n_kl_na_bc <- sum(is.na(caseControl_workdf_clinical_vars$kappa_lambda_ratio_3cat))

kl_table_cases <- matrix(c(n_kl_low_ucmm, n_kl_low_bc,
                            n_kl_normal_ucmm, n_kl_normal_bc,
                            n_kl_high_ucmm, n_kl_high_bc),
                          nrow = 3, byrow = TRUE)
kl_chisq_cases <- chisq.test(kl_table_cases)
pval_kl_cases <- format_pval(kl_chisq_cases$p.value)

kl_low_ucmm_str <- format_n_pct_manuscript(n_kl_low_ucmm, n_ucmm)
kl_low_bc_str <- format_n_pct_manuscript(n_kl_low_bc, n_bc)
kl_normal_ucmm_str <- format_n_pct_manuscript(n_kl_normal_ucmm, n_ucmm)
kl_normal_bc_str <- format_n_pct_manuscript(n_kl_normal_bc, n_bc)
kl_high_ucmm_str <- format_n_pct_manuscript(n_kl_high_ucmm, n_ucmm)
kl_high_bc_str <- format_n_pct_manuscript(n_kl_high_bc, n_bc)
kl_na_ucmm_str <- format_n_pct_manuscript(n_kl_na_ucmm, n_ucmm)
kl_na_bc_str <- format_n_pct_manuscript(n_kl_na_bc, n_bc)

# ------------------------------------------------------------------------------
# Create Table data frame (UCMM vs BC Cases)
# ------------------------------------------------------------------------------

table_ucmm_vs_bc <- data.frame(
  Characteristic = c(
    "N",
    "Age, mean (SD)",
    "Sex, n (%)",
    "  Male",
    "  Female",
    "Race, n (%)",
    "  White",
    "  African American",
    "BMI category, n (%)",
    "  <25 kg/m2",
    "  >=25 kg/m2",
    "Smoking history, n (%)",
    "  Ever smoker",
    "  Never smoker",
    "  Missing",
    "Education, n (%)",
    "  Below college",
    "  Some college or completed college",
    "  Graduate or professional degree",
    "  Missing",
    "ISS Stage, n (%)",
    "  Stage I",
    "  Stage II",
    "  Stage III",
    "  Could not be calculated",
    "  Missing",
    "Serum free light chain (kappa/lambda) ratio, n (%)",
    "  Normal (0.26-1.65)",
    "  High (>1.65)",
    "  Low (<0.26)",
    "  Missing"
  ),
  UCMM_Cases = c(
    n_ucmm_str,
    age_ucmm,
    "",
    male_ucmm_str,
    female_ucmm_str,
    "",
    white_ucmm_str,
    aa_ucmm_str,
    "",
    bmi_lt30_ucmm_str,
    bmi_ge30_ucmm_str,
    "",
    smoking_yes_ucmm_str,
    smoking_no_ucmm_str,
    smoking_na_ucmm_str,
    "",
    edu_below_ucmm_str,
    edu_some_ucmm_str,
    edu_grad_ucmm_str,
    edu_na_ucmm_str,
    "",
    iss1_ucmm_str,
    iss2_ucmm_str,
    iss3_ucmm_str,
    iss_cnc_ucmm_str,
    iss_na_ucmm_str,
    "",
    kl_normal_ucmm_str,
    kl_high_ucmm_str,
    kl_low_ucmm_str,
    kl_na_ucmm_str
  ),
  BC_Cases = c(
    n_bc_str,
    age_bc,
    "",
    male_bc_str,
    female_bc_str,
    "",
    white_bc_str,
    aa_bc_str,
    "",
    bmi_lt30_bc_str,
    bmi_ge30_bc_str,
    "",
    smoking_yes_bc_str,
    smoking_no_bc_str,
    smoking_na_bc_str,
    "",
    edu_below_bc_str,
    edu_some_bc_str,
    edu_grad_bc_str,
    edu_na_bc_str,
    "",
    iss1_bc_str,
    iss2_bc_str,
    iss3_bc_str,
    iss_cnc_bc_str,
    iss_na_bc_str,
    "",
    kl_normal_bc_str,
    kl_high_bc_str,
    kl_low_bc_str,
    kl_na_bc_str
  ),
  P_value = c(
    pval_n_cases,
    pval_age_cases,
    pval_sex_cases,
    "",
    "",
    pval_race_cases,
    "",
    "",
    pval_bmi_cases,
    "",
    "",
    pval_smoking_cases,
    "",
    "",
    "",
    pval_edu_cases,
    "",
    "",
    "",
    "",
    pval_iss_cases,
    "",
    "",
    "",
    "",
    "",
    pval_kl_cases,
    "",
    "",
    "",
    ""
  ),
  stringsAsFactors = FALSE
)

# Print the table
message("\n=== Table: UCMM Cases vs BC Cases (Demographic and clinical characteristics) ===\n")
print(table_ucmm_vs_bc, row.names = FALSE)

# ------------------------------------------------------------------------------
# Generate LaTeX code for UCMM vs BC Cases table
# ------------------------------------------------------------------------------

message("\n")
message("================================================================================")
message("LaTeX code for UCMM vs BC Cases Table:")
message("================================================================================\n")

latex_table_ucmm_bc <- paste0(
  "\\begin{table}[h]\n",
  "\\caption{Demographic and clinical characteristics: UCMM Cases vs BC Cases}\\label{tab:ucmm_vs_bc}\n",
  "\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}lccc}\n",
  "\\toprule\n",
  "\\textbf{Characteristic} & \\textbf{UCMM Cases} & \\textbf{BC Cases} & \\textbf{P-value} \\\\\n",
  "\\midrule\n",
  "N & ", n_ucmm_str, " & ", n_bc_str, " & -- \\\\\n",
  "Age, mean (SD) & ", age_ucmm, " & ", age_bc, " & ", pval_age_cases, " \\\\\n",
  "Sex, n (\\%) & & & ", pval_sex_cases, " \\\\\n",
  "\\quad Male & ", male_ucmm_str, " & ", male_bc_str, " & \\\\\n",
  "\\quad Female & ", female_ucmm_str, " & ", female_bc_str, " & \\\\\n",
  "Race, n (\\%) & & & ", pval_race_cases, " \\\\\n",
  "\\quad White & ", white_ucmm_str, " & ", white_bc_str, " & \\\\\n",
  "\\quad African American & ", aa_ucmm_str, " & ", aa_bc_str, " & \\\\\n",
  "BMI category, n (\\%) & & & ", pval_bmi_cases, " \\\\\n",
  "\\quad $<$25 kg/m$^2$ & ", bmi_lt30_ucmm_str, " & ", bmi_lt30_bc_str, " & \\\\\n",
  "\\quad $\\geq$25 kg/m$^2$ & ", bmi_ge30_ucmm_str, " & ", bmi_ge30_bc_str, " & \\\\\n",
  "Smoking history, n (\\%)$^{\\dagger}$ & & & ", pval_smoking_cases, " \\\\\n",
  "\\quad Ever smoker & ", smoking_yes_ucmm_str, " & ", smoking_yes_bc_str, " & \\\\\n",
  "\\quad Never smoker & ", smoking_no_ucmm_str, " & ", smoking_no_bc_str, " & \\\\\n",
  "\\quad Missing & ", smoking_na_ucmm_str, " & ", smoking_na_bc_str, " & \\\\\n",
  "Education, n (\\%)$^{\\dagger}$ & & & ", pval_edu_cases, " \\\\\n",
  "\\quad Below college & ", edu_below_ucmm_str, " & ", edu_below_bc_str, " & \\\\\n",
  "\\quad Some college or completed college & ", edu_some_ucmm_str, " & ", edu_some_bc_str, " & \\\\\n",
  "\\quad Graduate or professional degree & ", edu_grad_ucmm_str, " & ", edu_grad_bc_str, " & \\\\\n",
  "\\quad Missing & ", edu_na_ucmm_str, " & ", edu_na_bc_str, " & \\\\\n",
  "ISS Stage, n (\\%)$^{\\dagger}$ & & & ", pval_iss_cases, " \\\\\n",
  "\\quad Stage I & ", iss1_ucmm_str, " & ", iss1_bc_str, " & \\\\\n",
  "\\quad Stage II & ", iss2_ucmm_str, " & ", iss2_bc_str, " & \\\\\n",
  "\\quad Stage III & ", iss3_ucmm_str, " & ", iss3_bc_str, " & \\\\\n",
  "\\quad Could not be calculated & ", iss_cnc_ucmm_str, " & ", iss_cnc_bc_str, " & \\\\\n",
  "\\quad Missing & ", iss_na_ucmm_str, " & ", iss_na_bc_str, " & \\\\\n",
  "Serum free light chain ($\\kappa$/$\\lambda$) ratio, n (\\%)$^{\\dagger}$ & & & ", pval_kl_cases, " \\\\\n",
  "\\quad Normal (0.26--1.65) & ", kl_normal_ucmm_str, " & ", kl_normal_bc_str, " & \\\\\n",
  "\\quad High ($>$1.65) & ", kl_high_ucmm_str, " & ", kl_high_bc_str, " & \\\\\n",
  "\\quad Low ($<$0.26) & ", kl_low_ucmm_str, " & ", kl_low_bc_str, " & \\\\\n",
  "\\quad Missing & ", kl_na_ucmm_str, " & ", kl_na_bc_str, " & \\\\\n",
  "\\botrule\n",
  "\\end{tabular*}\n",
  "\\footnotetext{UCMM, University of Chicago Multiple Myeloma; BC, British Columbia (Canada). SD, standard deviation; BMI, body mass index; ISS, International Staging System. P-values calculated using t-test for continuous variables and chi-square test for categorical variables. $^{\\dagger}$Percentages are calculated using the total N as the denominator and may not sum to 100\\% due to missing data. Chi-square tests for smoking, education, ISS stage, and kappa/lambda ratio were performed on non-missing observations only (excluding ``Could not be calculated'' and ``Missing'' categories).}\n",
  "\\end{table}\n"
)

cat(latex_table_ucmm_bc)

# Save LaTeX table to file
writeLines(latex_table_ucmm_bc, "output/table_ucmm_vs_bc_cases.tex")
message("\nLaTeX table saved to: output/table_ucmm_vs_bc_cases.tex")

# ------------------------------------------------------------------------------
# Additional summary statistics for UCMM vs BC comparison
# ------------------------------------------------------------------------------

message("\n")
message("================================================================================")
message("Additional statistics for UCMM vs BC Cases comparison:")
message("================================================================================\n")

message("Total MM cases: ", n_ucmm + n_bc)
message("  UCMM Cases: ", n_ucmm, " (", round(100*n_ucmm/(n_ucmm + n_bc), 1), "%)")
message("  BC Cases: ", n_bc, " (", round(100*n_bc/(n_ucmm + n_bc), 1), "%)")
message("")
message("Mean age difference: ", round(mean(ucmm_cases$age) - mean(bc_cases$age), 1), " years")
message("  UCMM mean age: ", round(mean(ucmm_cases$age), 1))
message("  BC mean age: ", round(mean(bc_cases$age), 1))



################################################################################################################
# 5) Generate Table: Model Development Set vs Validation Set
#    Using same data splitting mechanism as in 03_real_data_analysis_Dec09 file
################################################################################################################

message("\n")
message("================================================================================")
message("GENERATING TABLE: Model Development Set vs Validation Set")
message("================================================================================\n")

# Set random seed for reproducibility (same as in 03_real_data_analysis_Dec09)
set.seed(789)

# Define split sizes (same as in 03_real_data_analysis_Dec09)
N_total_split <- nrow(pooled_metadata)
N_CV <- 500  # Cross-validation training pool (model development set)
N_test <- N_total_split - N_CV  # External held-out test set (validation set)

message(sprintf("Total patients: %d", N_total_split))
message(sprintf("  - Model development set: %d patients", N_CV))
message(sprintf("  - Validation set: %d patients", N_test))

# Stratified sampling to preserve outcome proportions
outcome_vec_split <- pooled_metadata$outcome
outcome_1_ids_split <- pooled_metadata$PooledData_ID[outcome_vec_split == 1]
outcome_0_ids_split <- pooled_metadata$PooledData_ID[outcome_vec_split == 0]

# Calculate proportion of cases
prop_cases_split <- sum(outcome_vec_split == 1) / N_total_split
n_test_cases_split <- round(N_test * prop_cases_split)
n_test_controls_split <- N_test - n_test_cases_split

# Sample test set with balanced outcomes (same mechanism as 03_real_data_analysis_Dec09)
test_ids_cases_split <- sample(outcome_1_ids_split, n_test_cases_split, replace = FALSE)
test_ids_controls_split <- sample(outcome_0_ids_split, n_test_controls_split, replace = FALSE)
test_ids_split <- c(test_ids_cases_split, test_ids_controls_split)

# Remaining IDs go to CV pool (model development set)
cv_ids_split <- setdiff(pooled_metadata$PooledData_ID, test_ids_split)

# Create subsets for model development and validation sets
dev_set <- pooled_metadata[pooled_metadata$PooledData_ID %in% cv_ids_split, ]
val_set <- pooled_metadata[pooled_metadata$PooledData_ID %in% test_ids_split, ]

n_dev <- nrow(dev_set)
n_val <- nrow(val_set)

message(sprintf("\nModel Development Set: %d patients", n_dev))
message(sprintf("  - Cases: %d, Controls: %d", sum(dev_set$outcome == 1), sum(dev_set$outcome == 0)))
message(sprintf("Validation Set: %d patients", n_val))
message(sprintf("  - Cases: %d, Controls: %d", sum(val_set$outcome == 1), sum(val_set$outcome == 0)))

# ------------------------------------------------------------------------------
# Compute statistics for each characteristic (Development vs Validation)
# ------------------------------------------------------------------------------

# 1. Sample size (N)
n_dev_str <- as.character(n_dev)
n_val_str <- as.character(n_val)
pval_n_split <- "--"

# 2. Age: mean (SD) with t-test
age_dev <- format_mean_sd(dev_set$age)
age_val <- format_mean_sd(val_set$age)
age_ttest_split <- t.test(dev_set$age, val_set$age)
pval_age_split <- format_pval(age_ttest_split$p.value)

# 3. Sex: n (%) with chi-square test
n_male_dev <- sum(dev_set$sex_binary == 1)
n_male_val <- sum(val_set$sex_binary == 1)
n_female_dev <- sum(dev_set$sex_binary == 0)
n_female_val <- sum(val_set$sex_binary == 0)

sex_table_split <- matrix(c(n_male_dev, n_male_val, n_female_dev, n_female_val),
                          nrow = 2, byrow = TRUE)
sex_chisq_split <- chisq.test(sex_table_split)
pval_sex_split <- format_pval(sex_chisq_split$p.value)

male_dev_str <- format_n_pct_manuscript(n_male_dev, n_dev)
male_val_str <- format_n_pct_manuscript(n_male_val, n_val)
female_dev_str <- format_n_pct_manuscript(n_female_dev, n_dev)
female_val_str <- format_n_pct_manuscript(n_female_val, n_val)

# 4. Race: n (%) with chi-square test
n_white_dev <- sum(dev_set$race_binary == 1)
n_white_val <- sum(val_set$race_binary == 1)
n_aa_dev <- sum(dev_set$race_binary == 0)
n_aa_val <- sum(val_set$race_binary == 0)

race_table_split <- matrix(c(n_white_dev, n_white_val, n_aa_dev, n_aa_val),
                           nrow = 2, byrow = TRUE)
race_chisq_split <- chisq.test(race_table_split)
pval_race_split <- format_pval(race_chisq_split$p.value)

white_dev_str <- format_n_pct_manuscript(n_white_dev, n_dev)
white_val_str <- format_n_pct_manuscript(n_white_val, n_val)
aa_dev_str <- format_n_pct_manuscript(n_aa_dev, n_dev)
aa_val_str <- format_n_pct_manuscript(n_aa_val, n_val)

# 5. BMI category: n (%) with chi-square test
n_bmi_ge30_dev <- sum(dev_set$bmi_binary == 1)
n_bmi_ge30_val <- sum(val_set$bmi_binary == 1)
n_bmi_lt30_dev <- sum(dev_set$bmi_binary == 0)
n_bmi_lt30_val <- sum(val_set$bmi_binary == 0)

bmi_table_split <- matrix(c(n_bmi_lt30_dev, n_bmi_lt30_val, n_bmi_ge30_dev, n_bmi_ge30_val),
                          nrow = 2, byrow = TRUE)
bmi_chisq_split <- chisq.test(bmi_table_split)
pval_bmi_split <- format_pval(bmi_chisq_split$p.value)

bmi_lt30_dev_str <- format_n_pct_manuscript(n_bmi_lt30_dev, n_dev)
bmi_lt30_val_str <- format_n_pct_manuscript(n_bmi_lt30_val, n_val)
bmi_ge30_dev_str <- format_n_pct_manuscript(n_bmi_ge30_dev, n_dev)
bmi_ge30_val_str <- format_n_pct_manuscript(n_bmi_ge30_val, n_val)

# 6. MM Status (Case/Control): n (%) - stratified by design, so no p-value needed
n_case_dev <- sum(dev_set$outcome == 1)
n_case_val <- sum(val_set$outcome == 1)
n_control_dev <- sum(dev_set$outcome == 0)
n_control_val <- sum(val_set$outcome == 0)

case_dev_str <- format_n_pct_manuscript(n_case_dev, n_dev)
case_val_str <- format_n_pct_manuscript(n_case_val, n_val)
control_dev_str <- format_n_pct_manuscript(n_control_dev, n_dev)
control_val_str <- format_n_pct_manuscript(n_control_val, n_val)

# 7. Study (UChicago/Canada): n (%) with chi-square test
n_uchicago_dev <- sum(dev_set$study == "UChicago")
n_uchicago_val <- sum(val_set$study == "UChicago")
n_canada_dev <- sum(dev_set$study == "Canada")
n_canada_val <- sum(val_set$study == "Canada")

study_table_split <- matrix(c(n_uchicago_dev, n_uchicago_val, n_canada_dev, n_canada_val),
                            nrow = 2, byrow = TRUE)
study_chisq_split <- chisq.test(study_table_split)
pval_study_split <- format_pval(study_chisq_split$p.value)

uchicago_dev_str <- format_n_pct_manuscript(n_uchicago_dev, n_dev)
uchicago_val_str <- format_n_pct_manuscript(n_uchicago_val, n_val)
canada_dev_str <- format_n_pct_manuscript(n_canada_dev, n_dev)
canada_val_str <- format_n_pct_manuscript(n_canada_val, n_val)

# ------------------------------------------------------------------------------
# Create Table data frame (Development vs Validation)
# ------------------------------------------------------------------------------

table_dev_vs_val <- data.frame(
  Characteristic = c(
    "N",
    "Age, mean (SD)",
    "Sex, n (%)",
    "  Male",
    "  Female",
    "Race, n (%)",
    "  White",
    "  African American",
    "BMI category, n (%)",
    "  <25 kg/m2",
    "  >=25 kg/m2",
    "MM Status, n (%)",
    "  Case",
    "  Control",
    "Study, n (%)",
    "  UChicago",
    "  Canada"
  ),
  Development_Set = c(
    n_dev_str,
    age_dev,
    "",
    male_dev_str,
    female_dev_str,
    "",
    white_dev_str,
    aa_dev_str,
    "",
    bmi_lt30_dev_str,
    bmi_ge30_dev_str,
    "",
    case_dev_str,
    control_dev_str,
    "",
    uchicago_dev_str,
    canada_dev_str
  ),
  Validation_Set = c(
    n_val_str,
    age_val,
    "",
    male_val_str,
    female_val_str,
    "",
    white_val_str,
    aa_val_str,
    "",
    bmi_lt30_val_str,
    bmi_ge30_val_str,
    "",
    case_val_str,
    control_val_str,
    "",
    uchicago_val_str,
    canada_val_str
  ),
  P_value = c(
    pval_n_split,
    pval_age_split,
    pval_sex_split,
    "",
    "",
    pval_race_split,
    "",
    "",
    pval_bmi_split,
    "",
    "",
    "--",
    "",
    "",
    pval_study_split,
    "",
    ""
  ),
  stringsAsFactors = FALSE
)

# Print the table
message("\n=== Table: Model Development Set vs Validation Set ===\n")
print(table_dev_vs_val, row.names = FALSE)

# ------------------------------------------------------------------------------
# Generate LaTeX code for Development vs Validation table
# ------------------------------------------------------------------------------

message("\n")
message("================================================================================")
message("LaTeX code for Development vs Validation Set Table:")
message("================================================================================\n")

latex_table_dev_val <- paste0(
  "\\begin{table}[h]\n",
  "\\caption{Demographic and clinical characteristics: Model Development Set vs Validation Set}\\label{tab:dev_vs_val}\n",
  "\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}lccc}\n",
  "\\toprule\n",
  "\\textbf{Characteristic} & \\textbf{Development Set} & \\textbf{Validation Set} & \\textbf{P-value} \\\\\n",
  " & \\textbf{(n = ", n_dev_str, ")} & \\textbf{(n = ", n_val_str, ")} & \\\\\n",
  "\\midrule\n",
  "N & ", n_dev_str, " & ", n_val_str, " & -- \\\\\n",
  "Age, mean (SD) & ", age_dev, " & ", age_val, " & ", pval_age_split, " \\\\\n",
  "Sex, n (\\%) & & & ", pval_sex_split, " \\\\\n",
  "\\quad Male & ", male_dev_str, " & ", male_val_str, " & \\\\\n",
  "\\quad Female & ", female_dev_str, " & ", female_val_str, " & \\\\\n",
  "Race, n (\\%) & & & ", pval_race_split, " \\\\\n",
  "\\quad White & ", white_dev_str, " & ", white_val_str, " & \\\\\n",
  "\\quad African American & ", aa_dev_str, " & ", aa_val_str, " & \\\\\n",
  "BMI category, n (\\%) & & & ", pval_bmi_split, " \\\\\n",
  "\\quad $<$25 kg/m$^2$ & ", bmi_lt30_dev_str, " & ", bmi_lt30_val_str, " & \\\\\n",
  "\\quad $\\geq$25 kg/m$^2$ & ", bmi_ge30_dev_str, " & ", bmi_ge30_val_str, " & \\\\\n",
  "MM Status, n (\\%) & & & -- \\\\\n",
  "\\quad Case & ", case_dev_str, " & ", case_val_str, " & \\\\\n",
  "\\quad Control & ", control_dev_str, " & ", control_val_str, " & \\\\\n",
  "Study, n (\\%) & & & ", pval_study_split, " \\\\\n",
  "\\quad UChicago & ", uchicago_dev_str, " & ", uchicago_val_str, " & \\\\\n",
  "\\quad Canada & ", canada_dev_str, " & ", canada_val_str, " & \\\\\n",
  "\\botrule\n",
  "\\end{tabular*}\n",
  "\\footnotetext{Development set used for 5-fold cross-validation model training; Validation set used for external held-out testing. SD, standard deviation; BMI, body mass index; MM, Multiple Myeloma. P-values calculated using t-test for continuous variables and chi-square test for categorical variables. Stratified sampling was used to preserve outcome proportions across sets.}\n",
  "\\end{table}\n"
)

cat(latex_table_dev_val)

# Save LaTeX table to file
writeLines(latex_table_dev_val, "output/table_dev_vs_val.tex")
message("\nLaTeX table saved to: output/table_dev_vs_val.tex")

# ------------------------------------------------------------------------------
# Additional summary statistics for Development vs Validation comparison
# ------------------------------------------------------------------------------

message("\n")
message("================================================================================")
message("Additional statistics for Development vs Validation comparison:")
message("================================================================================\n")

message("Data split summary:")
message("  Total patients: ", N_total_split)
message("  Development set: ", n_dev, " (", round(100*n_dev/N_total_split, 1), "%)")
message("  Validation set: ", n_val, " (", round(100*n_val/N_total_split, 1), "%)")
message("")
message("Case proportions (stratified sampling preserved):")
message("  Development set: ", round(100*n_case_dev/n_dev, 1), "% cases")
message("  Validation set: ", round(100*n_case_val/n_val, 1), "% cases")
message("  Overall: ", round(100*sum(pooled_metadata$outcome == 1)/N_total_split, 1), "% cases")


################################################################################################################
# 6) Generate Table: Summary of demographic and lifestyle variables across studies
#    UCMM Study, BC Case-Control Study, and Pooled Data
################################################################################################################

message("\n")
message("================================================================================")
message("GENERATING TABLE: Summary of demographic and lifestyle variables across studies")
message("================================================================================\n")

# Get sample sizes for each study
n_ucmm_total <- nrow(cohort_selected)
n_bc_total <- nrow(caseControl_selected)
n_pooled_total <- nrow(pooled_metadata)

# ------------------------------------------------------------------------------
# Compute statistics for UCMM Study
# ------------------------------------------------------------------------------

# Participant type
n_case_ucmm <- sum(cohort_selected$mm_status == "CASE")
n_control_ucmm <- sum(cohort_selected$mm_status == "CONTROL")

# Age
age_mean_ucmm <- mean(cohort_selected$age, na.rm = TRUE)
age_sd_ucmm <- sd(cohort_selected$age, na.rm = TRUE)

# Sex
n_male_ucmm_total <- sum(cohort_selected$sex_binary == 1)
n_female_ucmm_total <- sum(cohort_selected$sex_binary == 0)

# BMI
n_bmi_ge30_ucmm_total <- sum(cohort_selected$bmi_binary == 1)
n_bmi_lt30_ucmm_total <- sum(cohort_selected$bmi_binary == 0)

# Race
n_white_ucmm_total <- sum(cohort_selected$race_binary == 1)
n_aa_ucmm_total <- sum(cohort_selected$race_binary == 0)

# ------------------------------------------------------------------------------
# Compute statistics for BC Case-Control Study
# ------------------------------------------------------------------------------

# Participant type
n_case_bc <- sum(caseControl_selected$mm_status == "CASE")
n_control_bc <- sum(caseControl_selected$mm_status == "CONTROL")

# Age
age_mean_bc <- mean(caseControl_selected$age, na.rm = TRUE)
age_sd_bc <- sd(caseControl_selected$age, na.rm = TRUE)

# Sex
n_male_bc_total <- sum(caseControl_selected$sex_binary == 1)
n_female_bc_total <- sum(caseControl_selected$sex_binary == 0)

# BMI
n_bmi_ge30_bc_total <- sum(caseControl_selected$bmi_binary == 1)
n_bmi_lt30_bc_total <- sum(caseControl_selected$bmi_binary == 0)

# Race
n_white_bc_total <- sum(caseControl_selected$race_binary == 1)
n_aa_bc_total <- sum(caseControl_selected$race_binary == 0)

# ------------------------------------------------------------------------------
# Compute statistics for Pooled Data
# ------------------------------------------------------------------------------

# Participant type
n_case_pooled <- sum(pooled_metadata$outcome == 1)
n_control_pooled <- sum(pooled_metadata$outcome == 0)

# Age
age_mean_pooled <- mean(pooled_metadata$age, na.rm = TRUE)
age_sd_pooled <- sd(pooled_metadata$age, na.rm = TRUE)

# Sex
n_male_pooled <- sum(pooled_metadata$sex_binary == 1)
n_female_pooled <- sum(pooled_metadata$sex_binary == 0)

# BMI
n_bmi_ge30_pooled <- sum(pooled_metadata$bmi_binary == 1)
n_bmi_lt30_pooled <- sum(pooled_metadata$bmi_binary == 0)

# Race
n_white_pooled <- sum(pooled_metadata$race_binary == 1)
n_aa_pooled <- sum(pooled_metadata$race_binary == 0)

# ------------------------------------------------------------------------------
# Create formatted strings for LaTeX table
# ------------------------------------------------------------------------------

# Helper function for count and percentage (LaTeX escaped)
format_n_pct_latex <- function(count, total) {
  pct <- round(100 * count / total, 1)
  sprintf("%d (%.1f\\%%)", count, pct)
}

# Helper function for mean and SD
format_mean_sd_latex <- function(mean_val, sd_val) {
  sprintf("%.1f (%.1f)", mean_val, sd_val)
}

# ------------------------------------------------------------------------------
# Generate LaTeX code for the summary table
# ------------------------------------------------------------------------------

latex_table_summary <- paste0(
  "\\begin{table}[h]\n",
  "\\centering\n",
  "\\caption{Summary of demographic and lifestyle variables across the UCMM Study, BC Case-control Study, and the pooled dataset. Continuous variables are summarized as mean (standard deviation), and categorical variables as count (percentage).}\n",
  "\\label{Table1}\n",
  "\\renewcommand{\\arraystretch}{1.2}\n",
  "\\resizebox{1\\textwidth}{!}{%\n",
  "\\begin{tabular}{l l c c c}\n",
  "\\hline\n",
  "\\textbf{Characteristic} &  & \\textbf{UCMM Study} & \\textbf{BC Case-Control Study} & \\textbf{Pooled Data} \\\\\n",
  "\\hline\n",
  "\\multirow{2}{*}{Participant type} \n",
  "  & Case ($Y_i =1$)   & ", format_n_pct_latex(n_case_ucmm, n_ucmm_total), " & ", format_n_pct_latex(n_case_bc, n_bc_total), " & ", format_n_pct_latex(n_case_pooled, n_pooled_total), " \\\\\n",
  "  & Control ($Y_i =0$) & ", format_n_pct_latex(n_control_ucmm, n_ucmm_total), " & ", format_n_pct_latex(n_control_bc, n_bc_total), " & ", format_n_pct_latex(n_control_pooled, n_pooled_total), " \\\\\n",
  "\\hline\n",
  "Age (years)\n",
  "  & mean (SD)   & ", format_mean_sd_latex(age_mean_ucmm, age_sd_ucmm), " & ", format_mean_sd_latex(age_mean_bc, age_sd_bc), " & ", format_mean_sd_latex(age_mean_pooled, age_sd_pooled), " \\\\\n",
  "\\hline\n",
  "\\multirow{2}{*}{Sex} \n",
  "  & Male   & ", format_n_pct_latex(n_male_ucmm_total, n_ucmm_total), " & ", format_n_pct_latex(n_male_bc_total, n_bc_total), " & ", format_n_pct_latex(n_male_pooled, n_pooled_total), " \\\\\n",
  "  & Female & ", format_n_pct_latex(n_female_ucmm_total, n_ucmm_total), " & ", format_n_pct_latex(n_female_bc_total, n_bc_total), " & ", format_n_pct_latex(n_female_pooled, n_pooled_total), " \\\\\n",
  "\\hline\n",
  "\\multirow{2}{*}{BMI} \n",
  "  & $\\geq 25$    & ", format_n_pct_latex(n_bmi_ge30_ucmm_total, n_ucmm_total), " & ", format_n_pct_latex(n_bmi_ge30_bc_total, n_bc_total), " & ", format_n_pct_latex(n_bmi_ge30_pooled, n_pooled_total), " \\\\\n",
  "  & $<25$   & ", format_n_pct_latex(n_bmi_lt30_ucmm_total, n_ucmm_total), " & ", format_n_pct_latex(n_bmi_lt30_bc_total, n_bc_total), " & ", format_n_pct_latex(n_bmi_lt30_pooled, n_pooled_total), " \\\\\n",
  "\\hline\n",
  "\\multirow{2}{*}{Race} \n",
  "  & White     & ", format_n_pct_latex(n_white_ucmm_total, n_ucmm_total), " & ", format_n_pct_latex(n_white_bc_total, n_bc_total), " & ", format_n_pct_latex(n_white_pooled, n_pooled_total), " \\\\\n",
  "  & African American      & ", format_n_pct_latex(n_aa_ucmm_total, n_ucmm_total), " & ", format_n_pct_latex(n_aa_bc_total, n_bc_total), " & ", format_n_pct_latex(n_aa_pooled, n_pooled_total), " \\\\\n",
  "\\hline\n",
  "\\multirow{1}{*}{Total}\n",
  "  &      & ", format_n_pct_latex(n_ucmm_total, n_ucmm_total), " & ", format_n_pct_latex(n_bc_total, n_bc_total), " & ", format_n_pct_latex(n_pooled_total, n_pooled_total), " \\\\\n",
  "\\hline\n",
  "\\end{tabular}}\n",
  "\\end{table}\n"
)

message("\n")
message("================================================================================")
message("LaTeX code for Summary Table across Studies:")
message("================================================================================\n")

cat(latex_table_summary)

# Save LaTeX table to file
writeLines(latex_table_summary, "output/table_summary_across_studies.tex")
message("\nLaTeX table saved to: output/table_summary_across_studies.tex")

# ------------------------------------------------------------------------------
# Print summary statistics for verification
# ------------------------------------------------------------------------------

message("\n")
message("================================================================================")
message("Summary statistics for verification:")
message("================================================================================\n")

message("UCMM Study (n = ", n_ucmm_total, "):")
message("  Cases: ", n_case_ucmm, ", Controls: ", n_control_ucmm)
message("  Age: ", round(age_mean_ucmm, 1), " (", round(age_sd_ucmm, 1), ")")
message("  Male: ", n_male_ucmm_total, ", Female: ", n_female_ucmm_total)
message("  BMI >=25: ", n_bmi_ge30_ucmm_total, ", BMI <25: ", n_bmi_lt30_ucmm_total)
message("  White: ", n_white_ucmm_total, ", African American: ", n_aa_ucmm_total)

message("\nBC Case-Control Study (n = ", n_bc_total, "):")
message("  Cases: ", n_case_bc, ", Controls: ", n_control_bc)
message("  Age: ", round(age_mean_bc, 1), " (", round(age_sd_bc, 1), ")")
message("  Male: ", n_male_bc_total, ", Female: ", n_female_bc_total)
message("  BMI >=25: ", n_bmi_ge30_bc_total, ", BMI <25: ", n_bmi_lt30_bc_total)
message("  White: ", n_white_bc_total, ", African American: ", n_aa_bc_total)

message("\nPooled Data (n = ", n_pooled_total, "):")
message("  Cases: ", n_case_pooled, ", Controls: ", n_control_pooled)
message("  Age: ", round(age_mean_pooled, 1), " (", round(age_sd_pooled, 1), ")")
message("  Male: ", n_male_pooled, ", Female: ", n_female_pooled)
message("  BMI >=25: ", n_bmi_ge30_pooled, ", BMI <25: ", n_bmi_lt30_pooled)
message("  White: ", n_white_pooled, ", African American: ", n_aa_pooled)


################################################################################################################
# 7) t-SNE visualization: Cases vs Controls using 8 selected genes + 4 covariates
################################################################################################################

if (!require("Rtsne")) install.packages("Rtsne")
library(Rtsne)

message("\n")
message("================================================================================")
message("GENERATING t-SNE PLOT: Cases vs Controls")
message("================================================================================\n")

# Define the 8 selected genes
selected_genes <- c("IL1RAP", "ST5", "HERC6", "KL", "MYO1E", "ELK3", "CAPN2", "UBR4")

# Extract gene expression data for selected genes (genes x samples -> transpose to samples x genes)
gene_idx <- rownames(pooled_genebody_data_filtered_normalized) %in% selected_genes
tsne_genes <- t(pooled_genebody_data_filtered_normalized[gene_idx, ])

message("Selected genes found in data: ", paste(colnames(tsne_genes), collapse = ", "))
message("Gene matrix dimensions: ", nrow(tsne_genes), " samples x ", ncol(tsne_genes), " genes")

# Extract 4 covariates from pooled_metadata
tsne_covariates <- data.frame(
  age = pooled_metadata$age,
  sex_binary = pooled_metadata$sex_binary,
  bmi_binary = pooled_metadata$bmi_binary,
  race_binary = pooled_metadata$race_binary
)

# Combine genes and covariates into a single feature matrix
tsne_input <- cbind(tsne_genes, tsne_covariates)

message("Combined feature matrix dimensions: ", nrow(tsne_input), " samples x ", ncol(tsne_input), " features")
message("Features: ", paste(colnames(tsne_input), collapse = ", "))

# Scale features before t-SNE
tsne_input_scaled <- scale(tsne_input)

# Run t-SNE
set.seed(42)
tsne_result <- Rtsne(tsne_input_scaled, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 1000)

# Create data frame for plotting
tsne_df <- data.frame(
  tSNE1 = tsne_result$Y[, 1],
  tSNE2 = tsne_result$Y[, 2],
  Status = ifelse(pooled_metadata$outcome == 1, "MM Case", "Control")
)

# Plot (high-quality base R)
# Define colors: deep blue for Controls, deep red for MM Cases
col_case <- adjustcolor("#D73027", alpha.f = 0.65)
col_control <- adjustcolor("#4575B4", alpha.f = 0.65)
pt_colors <- ifelse(tsne_df$Status == "MM Case", col_case, col_control)
pt_shapes <- ifelse(tsne_df$Status == "MM Case", 16, 17)

# Plot controls first, then cases on top
plot_order <- order(tsne_df$Status == "MM Case")

png("output/tsne_cases_vs_controls.png", width = 7, height = 6, units = "in", res = 300)
par(mar = c(5, 5, 1.5, 1.5), family = "sans", cex.lab = 1.3, cex.axis = 1.1)
plot(tsne_df$tSNE1[plot_order], tsne_df$tSNE2[plot_order],
     col = pt_colors[plot_order],
     pch = pt_shapes[plot_order],
     cex = 1.2,
     xlab = "t-SNE 1",
     ylab = "t-SNE 2",
     main = "",
     las = 1,
     bty = "l",
     tcl = -0.3,
     mgp = c(3, 0.6, 0))
legend("topright",
       legend = c("MM Case", "Control"),
       col = c("#D73027", "#4575B4"),
       pch = c(16, 17),
       pt.cex = 1.5,
       cex = 1.1,
       bty = "n")
dev.off()

pdf("output/tsne_cases_vs_controls.pdf", width = 7, height = 6)
par(mar = c(5, 5, 1.5, 1.5), family = "sans", cex.lab = 1.3, cex.axis = 1.1)
plot(tsne_df$tSNE1[plot_order], tsne_df$tSNE2[plot_order],
     col = pt_colors[plot_order],
     pch = pt_shapes[plot_order],
     cex = 1.2,
     xlab = "t-SNE 1",
     ylab = "t-SNE 2",
     main = "",
     las = 1,
     bty = "l",
     tcl = -0.3,
     mgp = c(3, 0.6, 0))
legend("topright",
       legend = c("MM Case", "Control"),
       col = c("#D73027", "#4575B4"),
       pch = c(16, 17),
       pt.cex = 1.5,
       cex = 1.1,
       bty = "n")
dev.off()

message("\nt-SNE plot saved to: output/tsne_cases_vs_controls.png")
message("  Total samples plotted: ", nrow(tsne_df))
message("  MM Cases: ", sum(tsne_df$Status == "MM Case"))
message("  Controls: ", sum(tsne_df$Status == "Control"))


# ------------------------------------------------------------------------------
# Multi-panel t-SNE: color by each covariate to identify cluster drivers
# ------------------------------------------------------------------------------

# Add covariate labels to tsne_df
tsne_df$Study <- pooled_metadata$study
tsne_df$Race <- ifelse(pooled_metadata$race_binary == 1, "White", "African American")
tsne_df$Sex <- ifelse(pooled_metadata$sex_binary == 1, "Male", "Female")
tsne_df$BMI <- ifelse(pooled_metadata$bmi_binary == 1, ">=25", "<25")
tsne_df$Age_group <- cut(pooled_metadata$age, breaks = c(0, 55, 65, 75, Inf),
                          labels = c("<55", "55-65", "65-75", ">75"))

# Helper function for one t-SNE panel
plot_tsne_panel <- function(tsne_df, color_var, colors, labels, title, pch_val = 16) {
  color_vec <- colors[as.character(tsne_df[[color_var]])]
  plot(tsne_df$tSNE1, tsne_df$tSNE2,
       col = color_vec,
       pch = pch_val, cex = 0.9,
       xlab = "t-SNE 1", ylab = "t-SNE 2",
       main = title, cex.main = 1.2,
       las = 1, bty = "l", tcl = -0.3, mgp = c(2.5, 0.5, 0),
       cex.lab = 1.1, cex.axis = 0.9)
  legend("topright", legend = labels, col = colors[labels],
         pch = pch_val, pt.cex = 1.3, cex = 0.85, bty = "n")
}

# --- PNG version ---
png("output/tsne_multipanel.png", width = 14, height = 10, units = "in", res = 300)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), family = "sans")

# Panel A: MM Status
plot_tsne_panel(tsne_df, "Status",
  colors = c("MM Case" = adjustcolor("#D73027", 0.65), "Control" = adjustcolor("#4575B4", 0.65)),
  labels = c("MM Case", "Control"), title = "A) MM Status")

# Panel B: Study
plot_tsne_panel(tsne_df, "Study",
  colors = c("UChicago" = adjustcolor("#E66101", 0.65), "Canada" = adjustcolor("#5E3C99", 0.65)),
  labels = c("UChicago", "Canada"), title = "B) Study")

# Panel C: Race
plot_tsne_panel(tsne_df, "Race",
  colors = c("White" = adjustcolor("#1B9E77", 0.65), "African American" = adjustcolor("#D95F02", 0.65)),
  labels = c("White", "African American"), title = "C) Race")

# Panel D: Sex
plot_tsne_panel(tsne_df, "Sex",
  colors = c("Male" = adjustcolor("#7570B3", 0.65), "Female" = adjustcolor("#E7298A", 0.65)),
  labels = c("Male", "Female"), title = "D) Sex")

# Panel E: BMI
plot_tsne_panel(tsne_df, "BMI",
  colors = c(">=25" = adjustcolor("#A6761D", 0.65), "<25" = adjustcolor("#66A61E", 0.65)),
  labels = c(">=25", "<25"), title = "E) BMI (kg/m2)")

# Panel F: Age group
plot_tsne_panel(tsne_df, "Age_group",
  colors = c("<55" = adjustcolor("#FEE08B", 0.8), "55-65" = adjustcolor("#FDAE61", 0.8),
             "65-75" = adjustcolor("#F46D43", 0.8), ">75" = adjustcolor("#A50026", 0.8)),
  labels = c("<55", "55-65", "65-75", ">75"), title = "F) Age group (years)")

dev.off()

# --- PDF version ---
pdf("output/tsne_multipanel.pdf", width = 14, height = 10)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), family = "sans")

plot_tsne_panel(tsne_df, "Status",
  colors = c("MM Case" = adjustcolor("#D73027", 0.65), "Control" = adjustcolor("#4575B4", 0.65)),
  labels = c("MM Case", "Control"), title = "A) MM Status")

plot_tsne_panel(tsne_df, "Study",
  colors = c("UChicago" = adjustcolor("#E66101", 0.65), "Canada" = adjustcolor("#5E3C99", 0.65)),
  labels = c("UChicago", "Canada"), title = "B) Study")

plot_tsne_panel(tsne_df, "Race",
  colors = c("White" = adjustcolor("#1B9E77", 0.65), "African American" = adjustcolor("#D95F02", 0.65)),
  labels = c("White", "African American"), title = "C) Race")

plot_tsne_panel(tsne_df, "Sex",
  colors = c("Male" = adjustcolor("#7570B3", 0.65), "Female" = adjustcolor("#E7298A", 0.65)),
  labels = c("Male", "Female"), title = "D) Sex")

plot_tsne_panel(tsne_df, "BMI",
  colors = c(">=25" = adjustcolor("#A6761D", 0.65), "<25" = adjustcolor("#66A61E", 0.65)),
  labels = c(">=25", "<25"), title = "E) BMI (kg/m2)")

plot_tsne_panel(tsne_df, "Age_group",
  colors = c("<55" = adjustcolor("#FEE08B", 0.8), "55-65" = adjustcolor("#FDAE61", 0.8),
             "65-75" = adjustcolor("#F46D43", 0.8), ">75" = adjustcolor("#A50026", 0.8)),
  labels = c("<55", "55-65", "65-75", ">75"), title = "F) Age group (years)")

dev.off()

message("\nMulti-panel t-SNE plot saved to: output/tsne_multipanel.png and output/tsne_multipanel.pdf")


# ------------------------------------------------------------------------------
# Multi-panel t-SNE using ONLY 4 demographic covariates (no genes)
# ------------------------------------------------------------------------------

message("\n")
message("================================================================================")
message("GENERATING t-SNE (covariates only): 4 demographic variables, no genes")
message("================================================================================\n")

# Build feature matrix with only 4 covariates
tsne_input_cov_only <- data.frame(
  age = pooled_metadata$age,
  sex_binary = pooled_metadata$sex_binary,
  bmi_binary = pooled_metadata$bmi_binary,
  race_binary = pooled_metadata$race_binary
)

tsne_input_cov_only_scaled <- scale(tsne_input_cov_only)

# Run t-SNE
set.seed(42)
tsne_result_cov <- Rtsne(tsne_input_cov_only_scaled, dims = 2, perplexity = 30,
                          verbose = TRUE, max_iter = 1000, check_duplicates = FALSE)

# Create data frame for plotting
tsne_df_cov <- data.frame(
  tSNE1 = tsne_result_cov$Y[, 1],
  tSNE2 = tsne_result_cov$Y[, 2],
  Status = ifelse(pooled_metadata$outcome == 1, "MM Case", "Control"),
  Race = ifelse(pooled_metadata$race_binary == 1, "White", "African American"),
  Sex = ifelse(pooled_metadata$sex_binary == 1, "Male", "Female"),
  BMI = ifelse(pooled_metadata$bmi_binary == 1, ">=25", "<25"),
  Age_group = cut(pooled_metadata$age, breaks = c(0, 55, 65, 75, Inf),
                   labels = c("<55", "55-65", "65-75", ">75"))
)

# --- PNG version ---
png("output/tsne_multipanel_covariates_only.png", width = 14, height = 5, units = "in", res = 300)
par(mfrow = c(1, 5), mar = c(4, 4, 3, 1), family = "sans")

plot_tsne_panel(tsne_df_cov, "Status",
  colors = c("MM Case" = adjustcolor("#D73027", 0.65), "Control" = adjustcolor("#4575B4", 0.65)),
  labels = c("MM Case", "Control"), title = "A) MM Status")

plot_tsne_panel(tsne_df_cov, "Race",
  colors = c("White" = adjustcolor("#1B9E77", 0.65), "African American" = adjustcolor("#D95F02", 0.65)),
  labels = c("White", "African American"), title = "B) Race")

plot_tsne_panel(tsne_df_cov, "Sex",
  colors = c("Male" = adjustcolor("#7570B3", 0.65), "Female" = adjustcolor("#E7298A", 0.65)),
  labels = c("Male", "Female"), title = "C) Sex")

plot_tsne_panel(tsne_df_cov, "BMI",
  colors = c(">=25" = adjustcolor("#A6761D", 0.65), "<25" = adjustcolor("#66A61E", 0.65)),
  labels = c(">=25", "<25"), title = "D) BMI (kg/m2)")

plot_tsne_panel(tsne_df_cov, "Age_group",
  colors = c("<55" = adjustcolor("#FEE08B", 0.8), "55-65" = adjustcolor("#FDAE61", 0.8),
             "65-75" = adjustcolor("#F46D43", 0.8), ">75" = adjustcolor("#A50026", 0.8)),
  labels = c("<55", "55-65", "65-75", ">75"), title = "E) Age group (years)")

dev.off()

# --- PDF version ---
pdf("output/tsne_multipanel_covariates_only.pdf", width = 14, height = 5)
par(mfrow = c(1, 5), mar = c(4, 4, 3, 1), family = "sans")

plot_tsne_panel(tsne_df_cov, "Status",
  colors = c("MM Case" = adjustcolor("#D73027", 0.65), "Control" = adjustcolor("#4575B4", 0.65)),
  labels = c("MM Case", "Control"), title = "A) MM Status")

plot_tsne_panel(tsne_df_cov, "Race",
  colors = c("White" = adjustcolor("#1B9E77", 0.65), "African American" = adjustcolor("#D95F02", 0.65)),
  labels = c("White", "African American"), title = "B) Race")

plot_tsne_panel(tsne_df_cov, "Sex",
  colors = c("Male" = adjustcolor("#7570B3", 0.65), "Female" = adjustcolor("#E7298A", 0.65)),
  labels = c("Male", "Female"), title = "C) Sex")

plot_tsne_panel(tsne_df_cov, "BMI",
  colors = c(">=25" = adjustcolor("#A6761D", 0.65), "<25" = adjustcolor("#66A61E", 0.65)),
  labels = c(">=25", "<25"), title = "D) BMI (kg/m2)")

plot_tsne_panel(tsne_df_cov, "Age_group",
  colors = c("<55" = adjustcolor("#FEE08B", 0.8), "55-65" = adjustcolor("#FDAE61", 0.8),
             "65-75" = adjustcolor("#F46D43", 0.8), ">75" = adjustcolor("#A50026", 0.8)),
  labels = c("<55", "55-65", "65-75", ">75"), title = "E) Age group (years)")

dev.off()

message("\nCovariates-only multi-panel t-SNE saved to: output/tsne_multipanel_covariates_only.png and output/tsne_multipanel_covariates_only.pdf")


# ------------------------------------------------------------------------------
# 1x3 comparison panel: Cases vs Controls under 3 feature sets
# A) 8 genes + 4 covariates   B) 4 covariates only   C) 8 genes only
# ------------------------------------------------------------------------------

message("\n")
message("================================================================================")
message("GENERATING t-SNE comparison panel: 3 feature sets")
message("================================================================================\n")

# t-SNE with 8 genes only
tsne_genes_only_scaled <- scale(tsne_genes)

set.seed(42)
tsne_result_genes <- Rtsne(tsne_genes_only_scaled, dims = 2, perplexity = 30,
                            verbose = TRUE, max_iter = 1000)

tsne_df_genes <- data.frame(
  tSNE1 = tsne_result_genes$Y[, 1],
  tSNE2 = tsne_result_genes$Y[, 2],
  Status = ifelse(pooled_metadata$outcome == 1, "MM Case", "Control")
)

# Helper for the comparison panels (cases vs controls only)
plot_tsne_comparison <- function(df, title) {
  col_vec <- ifelse(df$Status == "MM Case",
                    adjustcolor("#D73027", 0.65),
                    adjustcolor("#4575B4", 0.65))
  pch_vec <- ifelse(df$Status == "MM Case", 16, 17)
  plot_ord <- order(df$Status == "MM Case")
  plot(df$tSNE1[plot_ord], df$tSNE2[plot_ord],
       col = col_vec[plot_ord], pch = pch_vec[plot_ord],
       cex = 1.0,
       xlab = "t-SNE 1", ylab = "t-SNE 2",
       main = title, cex.main = 1.2,
       las = 1, bty = "l", tcl = -0.3, mgp = c(2.5, 0.5, 0),
       cex.lab = 1.1, cex.axis = 0.9)
  legend("topleft", legend = c("MM Case", "Control"),
         col = c("#D73027", "#4575B4"), pch = c(16, 17),
         pt.cex = 1.3, cex = 0.9, bty = "n")
}

# --- PNG version ---
png("output/tsne_comparison_3panels.png", width = 14, height = 5, units = "in", res = 300)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1), family = "sans")

plot_tsne_comparison(tsne_df, "A) 8 Genes + 4 Covariates")
plot_tsne_comparison(tsne_df_cov, "B) 4 Covariates Only")
plot_tsne_comparison(tsne_df_genes, "C) 8 Genes Only")

dev.off()

# --- PDF version ---
pdf("output/tsne_comparison_3panels.pdf", width = 14, height = 5)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1), family = "sans")

plot_tsne_comparison(tsne_df, "A) 8 Genes + 4 Covariates")
plot_tsne_comparison(tsne_df_cov, "B) 4 Covariates Only")
plot_tsne_comparison(tsne_df_genes, "C) 8 Genes Only")

dev.off()

message("\n3-panel t-SNE comparison saved to: output/tsne_comparison_3panels.png and output/tsne_comparison_3panels.pdf")

