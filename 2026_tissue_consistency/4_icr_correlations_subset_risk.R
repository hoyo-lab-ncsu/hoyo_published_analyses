

# Imprintome data: icr, cpg
imp_str = "ICR"

# Marginal direciton of correlation: patient, probe
marg_dir = "probe"


# Set imprintome data for function calls below
if (imp_str == "icr"){
  imp_long = icr_long
} else {
  imp_long = cpg_long
}


# Pearson's rho split by tissue type, icr confidence                      ######
#_______________________________________________________________________________

df_rho <- marginal_grouped_agreement(imp_long = imp_long, marg_dir = marg_dir, 
                                     grouping_var = "icr_conf", agreement = "pearson",
                                     db_flag = TRUE, grouping_add_all = TRUE)
df_stat <- pairwise_stats(df = df_rho, db_flag = TRUE, partition_vars = "icr_conf", pairwise_var = "tissue")
sum_df <- subset(df_rho) %>% 
  group_by(icr_conf, tissue) %>%  summarize(mean = mean(rho), sd = sd(rho), sem = sd(rho)/sqrt(length(rho)), n = n())

plot_agreement_by_confidence(df_rho, sum_df, df_stat, x_var = "tissue",
                             subset_level = NA, fill_var = NA, facet_var = "icr_conf",
                             xlab_str = "Tissue", 
                             ylab_str = paste0("Rho (",toupper(imp_str)," Beta)"),
                             output_dir_path = output_dir_path,
                             export_filename = paste0("4_",marg_dir,"-wise_", toupper(imp_str),"_Rho.jpg"),
                             fig_dim = c(6,3), db_flag = TRUE, pretty_x = T)




# Pearson's rho split by tissue type, icr confidence, and risk factor: BMI ######
#_______________________________________________________________________________


df_rho <- marginal_grouped_agreement(
  imp_long = imp_long, marg_dir = marg_dir, grouping_var = 
    ifelse(marg_dir == "patient", "icr_conf", "mat_bmi_hi"), 
  agreement = "pearson", db_flag = TRUE, grouping_add_all = TRUE)
# df_stat <- pairwise_stats(df = df_rho,db_flag = TRUE, partition_vars = "icr_conf", pairwise_var = "tissue")
sum_df <- subset(df_rho) %>% 
  group_by(icr_conf, tissue, mat_bmi_hi) %>%  
  summarize(mean = mean(rho), sd = sd(rho), sem = sd(rho)/sqrt(length(rho)), n = n())

plot_agreement_by_confidence(df_rho = df_rho, sum_df = sum_df, df_stat = NULL, x_var = "tissue",
                             subset_level = NA, fill_var = "mat_bmi_hi", facet_var = "icr_conf",
                             xlab_str = "Tissue", 
                             ylab_str = paste0("Rho (",toupper(imp_str)," Beta)"),
                             output_dir_path = output_dir_path,
                             export_filename = paste0("4_",marg_dir,"-wise_", toupper(imp_str),"_MAT_BMI_Rho.jpg"),
                             fig_dim = c(5,2.75), db_flag = TRUE, pretty_x = T)


mod <- glm(data = filter(df_rho, icr_conf != "All"), formula = "rho ~ mat_bmi_hi + tissue + icr_conf")
res<-summary(mod)
View(res$coefficients)
res

# Pearson's rho split by tissue type, icr confidence, and risk factor: aces ######
#_______________________________________________________________________________

df_rho <- marginal_grouped_agreement(
  imp_long = imp_long, marg_dir = marg_dir, grouping_var = 
    ifelse(marg_dir == "patient", "icr_conf", "mat_aces_hi"), 
  agreement = "pearson", db_flag = TRUE, grouping_add_all = TRUE)
df_rho <- df_rho %>% filter(!is.na(mat_aces_hi))
# df_stat <- pairwise_stats(df = df_rho,db_flag = TRUE, partition_vars = "icr_conf", pairwise_var = "tissue")
sum_df <- subset(df_rho) %>% 
  group_by(icr_conf, tissue, mat_aces_hi) %>%  summarize(mean = mean(rho), sd = sd(rho), sem = sd(rho)/sqrt(length(rho)), n = n())

plot_agreement_by_confidence(df_rho = df_rho, sum_df = sum_df, df_stat = NULL, x_var = "tissue",
                             subset_level = NA, fill_var = "mat_aces_hi", facet_var = "icr_conf",
                             xlab_str = "Tissue", 
                             ylab_str = paste0("Rho (",toupper(imp_str)," Beta)"),
                             output_dir_path = output_dir_path,
                             export_filename = paste0("4_",marg_dir,"-wise_", toupper(imp_str),"_MAT_ACES_Rho.jpg"),
                             fig_dim = c(5,2.75), db_flag = TRUE, pretty_x = T)

mod <- glm(data = filter(df_rho, icr_conf != "All"), formula = "rho ~ mat_aces_hi + tissue + icr_conf")
res<-summary(mod)
View(res$coefficients)
res
