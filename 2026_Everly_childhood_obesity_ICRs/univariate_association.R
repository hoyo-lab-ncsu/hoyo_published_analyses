# Association of maternal and child characteristics with childhood obesity


# one example is given for a categorical variable
# -----------------------------------------------------------------------------------------------
# Categorical example: maternal smoking
# model for sustained obesity vs. non-obese
table(study_data_birth_sustained_vs_non$child_obesity, study_data_birth_sustained_vs_non$mat_smoking) #sample sizes
univariate_model_maternal_smoking_sustained <- glm(child_obesity ~ mat_smoking, data = study_data_birth_sustained_vs_non, family=binomial)
summary(univariate_model_maternal_smoking_sustained) # gives estimates and p-values
confint(univariate_model_maternal_smoking_sustained) # gives 95% CIs

# model for intermittent obesity vs. non-obese
study_data_birth_intermittent_vs_non <- study_data_birth %>%
  filter(child_obesity == "Intermittent obesity" | child_obesity == "Non-obese")
table(study_data_birth_intermittent_vs_non$child_obesity, study_data_birth_intermittent_vs_non$mat_smoking) #sample sizes
univariate_model_maternal_smoking_intermittent <- glm(child_obesity ~ mat_smoking, data = study_data_birth_intermittent_vs_non, family=binomial)
summary(univariate_model_maternal_smoking_intermittent) # gives estimates and p-values
confint(univariate_model_maternal_smoking_intermittent) # gives 95% CIs
