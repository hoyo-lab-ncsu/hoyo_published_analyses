# THIS SCRIPT WILL OUTLINE THE CpG- AND ICR-LEVEL ANALYSES TO IDENTIFY SIGNIFICANT ASSOCIATIONS WITH SUSTAINED
# OBESITY IN UMBILICAL CORD BLOOD

# Load libraries
library(dplyr)
library(tdhia)
library(qvalue)
library(tibble)

# 1. First complete 0_download_GEO_data.R. This script relies on the downloaded data
# ============================================================================================================
# If starting with pre-processed CpG beta dataframe:
predictor_cpg_df <- processed_cpg_df_umbilical %>%
  rename(sample = Sample.name)
row.names(predictor_cpg_df) <- predictor_cpg_df$sample # rename the rownames to the sample name

# 2. Load ICR map
# ============================================================================================================
icr_mapping = tdhia::mapping_cpg_icr_ids # This file is required to connect CpG sites to their ICR ID

icr_mapping_unique <- icr_mapping %>% # Some CpGs have >1 probe so they are repeated in the icr_mapping table
  distinct(CpG_id, .keep_all = TRUE) # This will remove the duplicates

# 3. Prepare study data
# ============================================================================================================
# Rename columns to easier variable names to work with
study_data <- study_data %>%
  rename(child_obesity = Childhood.Obesity.Status,
         mat_smoking = Maternal.Smoking,
         mat_obesity = Maternal.Obesity,
         mat_education = Maternal.Education,
         mat_race_eth = Maternal.Race.Ethnicity,
         breastfeeding = Ever.Breastfed,
         child_sex = Baby.Sex,
         child_age = Child.age.at.sample.collection..months.,
         sample = X.Sample.name..from.processed.data.table.,
         tissue = X..tissue
  )

# Ensure all categorical variables are factors
study_data <- study_data %>%
  mutate(child_obesity = as.factor(child_obesity),
         mat_smoking = as.factor(mat_smoking),
         mat_obesity = as.factor(mat_obesity),
         mat_education = as.factor(mat_education),
         mat_race_eth  = as.factor(mat_race_eth),
         breastfeeding = as.factor(breastfeeding),
         child_sex = as.factor(child_sex)
  )

row.names(study_data) <- study_data$sample # rename the rownames to the sample name

# Filter for only umbilical cord blood (collected at birth/age 0)
study_data_birth <- study_data %>%
  filter(tissue == "Umbilical cord blood")

# Filter for only peripheral blood (collected from children aged 8-16 years old)
study_data_late_childhood <- study_data %>%
  filter(tissue == "Peripheral blood")

# ============================================================================================================
# STEPS 4-7 WILL OUTLINE THE CpG- AND ICR-LEVEL ANALYSES TO IDENTIFY SIGNIFICANT ASSOCIATIONS WITH SUSTAINED
# OBESITY IN UMBILICAL CORD BLOOD
# ============================================================================================================

# 4. CpG Level logistic regression to identify CpGs at birth that are associated with sustained childhood obesity
# ============================================================================================================
# Filter data for only sustained obesity and non-obese (no intermittent obesity children)
study_data_birth_sustained_vs_non <- study_data_birth %>%
  filter(child_obesity == "Sustained obesity" | child_obesity == "Non-obese")

# Filter the data and cpg_df so that they have the same samples
common_rows <- intersect(rownames(study_data_birth_sustained_vs_non), rownames(predictor_cpg_df)) #find common rows in metadata and cpg matrix
predictor_cpg_birth_sustained_vs_non  <- predictor_cpg_df[common_rows, , drop = FALSE] #only keeps common rows in CpG matrix
study_data_birth_sustained_vs_non <- study_data_birth_sustained_vs_non[common_rows, , drop = FALSE] #only keeps common rows in study data

# Remove the sample name column so it's not treated as a variable
predictor_cpg_birth_sustained_vs_non <- predictor_cpg_birth_sustained_vs_non %>% dplyr::select(-sample)

#Check that the study data and CpG beta dataframe have the same order of patients
identical(rownames(study_data_birth_sustained_vs_non), rownames(predictor_cpg_birth_sustained_vs_non))

# Run CpG-level models using the analyze_association package in TDHIA
cpg_model_sustained_obesity <- tdhia::analyze_association(R= dplyr::select(study_data_birth_sustained_vs_non, child_obesity),
                                                          P= predictor_cpg_birth_sustained_vs_non,
                                                          Pe=dplyr::select(study_data_birth_sustained_vs_non, mat_smoking, mat_obesity,
                                                                           mat_education, mat_race_eth, breastfeeding, child_sex),
                                                          rm.na.Pe=TRUE,
                                                          family='binomial',
                                                          max_p_val = 0.05,
                                                          impute_na = FALSE, # the way I processed the CpG beta matrix does not have any NAs so no data imputation would occur even if TRUE
                                                          n_p_adj = 6711) # the input CpG beta dataframe has exactly 6711 CpG sites
# Save results and calculate q-values
results_cpg_sustained_obesity <- list()
results_cpg_sustained_obesity <- cpg_model_sustained_obesity$imp_site
qobj <- qvalue(p = results_cpg_sustained_obesity$P_VAL, fdr.level = 0.05)
results_cpg_sustained_obesity$CpG_qvalue <- qobj$qvalues

# Calculate estimates, odds ratios, and 95% confidence intervals (CIs) per 1% increase in methylation
# 95% CIs for odds ratios were derived using the Wald method based on the
# approximate normal distribution of logistic regression coefficient estimates.
results_cpg_sustained_obesity <- results_cpg_sustained_obesity %>%
  mutate(StdError_1pct = StdError/100,
         Estimate_1pct = Estimate/100,
         NormalApprox_Estimate_Lower95 = Estimate_1pct - 1.96*StdError_1pct,
         NormalApprox_Estimate_Upper95 = Estimate_1pct + 1.96*StdError_1pct,
         Odds_ratio_1pct = exp(Estimate_1pct),
         NormalApprox_OR_Lower95 = exp(Estimate_1pct - 1.96*StdError_1pct),
         NormalApprox_OR_Upper95 = exp(Estimate_1pct + 1.96*StdError_1pct)
  )

# Annotate CpGs with their corresponding ICR information
results_cpg_sustained_obesity <- results_cpg_sustained_obesity %>%
  left_join(y = dplyr::select(icr_mapping_unique, c(CpG_id, CpG_chr, CpG_start, CpG_stop, ICR_id, ICR_chr, ICR_start, ICR_stop)),
            by = dplyr::join_by("Variable" == "CpG_id"))


# 5. ICR Level PC Regression to identify ICRs at birth that are associated with sustained childhood obesity
# ============================================================================================================
pc_regression_sustained_obesity <- pc_regression_test(
  cpg_beta = predictor_cpg_birth_sustained_vs_non,
  m_value_transform = FALSE, # I used beta values throughout my analyses
  data_norm_type = "n1",
  pct_variance = 0.8,
  df_study = study_data_birth_sustained_vs_non,
  outcome = "child_obesity",
  covariates = c("mat_smoking", "mat_obesity", "mat_education",
                 "mat_race_eth", "breastfeeding", "child_sex"),
  Patient_ID = "sample",
  family = "binomial",
  icr_ids = NULL,
  min_cpg = 1,
  verbose = TRUE,
  n.cores = 1
)


# 6. ICR Level KM Regression to identify ICRs at birth that are associated with sustained childhood obesity
# ============================================================================================================
# The SKAT function requires cpg_id to be rownames so the dataframe must be transposed:
transposed_predictor_cpg_birth_sustained_vs_non <- as.data.frame(t(predictor_cpg_birth_sustained_vs_non))

# For SKAT, coerce factors into numerical
skat_data <- study_data_birth_sustained_vs_non %>%
  mutate(mat_race_eth = case_when(
    mat_race_eth == "Hispanic" ~ 0,
    mat_race_eth == "NHB" ~ 1,
    mat_race_eth == "NHW" ~ 2)) %>%
  mutate(mat_education = case_when(
    mat_education == "Not a college graduate" ~ 0,
    mat_education == "College graduate" ~ 1)) %>%
  mutate(child_obesity = case_when(
    child_obesity == "Sustained obesity" ~ 0,
    child_obesity == "Non-obese" ~ 1)) %>%
  mutate(across(c("mat_smoking", "mat_obesity","breastfeeding", "child_sex",
                  "mat_race_eth", "mat_education"), ~ as.numeric(as.character(.x))))

# Note that this SKAT function typically takes around 25 minutes to run
skat_sustained_obesity <- skat_icr_test(
  cpg_betas = transposed_predictor_cpg_birth_sustained_vs_non,
  df_study = skat_data,
  response = "child_obesity",
  predictors = c("mat_smoking", "mat_obesity", "mat_education",
                 "mat_race_eth", "breastfeeding", "child_sex"),
  method = "optimal.adj",
  out_type = "D", # "D": dichotomous.
  icr_ids = NULL,
  min_cpg = 1,
  db_flag = TRUE,
  m_value_transform = FALSE, # I used beta values throughout my analyses
  scaling = TRUE,
  verbose = TRUE,
  n.cores = 1
)

# 7. Merge all results files (CpG analysis, PC regression analysis, and KM regression analysis) into a single dataframe
# ============================================================================================================
# rename some columns so it's clear where the data is from
pc_regression_sustained_obesity <- pc_regression_sustained_obesity %>%
  rename(PC_Regression_p_value = p_value,
         PC_Regression_adj_p = adj_p_value,
         PC_Regression_q_value = q_value)

intermediate_results_df <- merge(results_cpg_sustained_obesity,
                                 pc_regression_sustained_obesity, by="ICR_id", all.x=TRUE)

final_results_sustained_obesity_birth <- merge(intermediate_results_df,
                                               skat_sustained_obesity, by.x="ICR_id", by.y="icr_id")

  
