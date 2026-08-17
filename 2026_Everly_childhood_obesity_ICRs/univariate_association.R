# Association of maternal and child characteristics with childhood obesity
# One example is given for a categorical variable (maternal smoking)

# First complete 0_download_GEO_data.R. This script relies on the downloaded data

# Prepare study data (this is identical to Step #3 in "2_umbilical_ICRs_associated_w_sustained_obesity.R", so if you 
# have already completed this step you do not need to repeat)
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


# Calculate marginal association of maternal smoking with sustained and intermittent childhood obesity
# ============================================================================================================
# Categorical example: maternal smoking
# Model for sustained obesity vs. non-obese
# ----------------------------------------------
study_data_birth_sustained_vs_non <- study_data_birth %>%
  filter(child_obesity == "Sustained obesity" | child_obesity == "Non-obese")
table(study_data_birth_sustained_vs_non$child_obesity, study_data_birth_sustained_vs_non$mat_smoking) #sample sizes
univariate_model_maternal_smoking_sustained <- glm(child_obesity ~ mat_smoking, data = study_data_birth_sustained_vs_non, family=binomial)
summary(univariate_model_maternal_smoking_sustained) # gives estimates and p-values
confint(univariate_model_maternal_smoking_sustained) # gives 95% CIs

# Model for intermittent obesity vs. non-obese
# ----------------------------------------------
study_data_birth_intermittent_vs_non <- study_data_birth %>%
  filter(child_obesity == "Intermittent obesity" | child_obesity == "Non-obese")
table(study_data_birth_intermittent_vs_non$child_obesity, study_data_birth_intermittent_vs_non$mat_smoking) #sample sizes
univariate_model_maternal_smoking_intermittent <- glm(child_obesity ~ mat_smoking, data = study_data_birth_intermittent_vs_non, family=binomial)
summary(univariate_model_maternal_smoking_intermittent) # gives estimates and p-values
confint(univariate_model_maternal_smoking_intermittent) # gives 95% CIs
