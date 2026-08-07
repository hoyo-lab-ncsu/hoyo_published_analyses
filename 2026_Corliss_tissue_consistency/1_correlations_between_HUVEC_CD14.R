
# User Settings
# Marginal direction of correlation: patient, probe
marg_dir = "patient"


# Figure 2A) Agreement of ICRs between HUVECS and CD14            ####
#_______________________________________________________________________________

df_rho <- marginal_grouped_agreement(imp_long = icr_long, marg_dir = marg_dir, 
                                      grouping_var = "icr_conf", agreement = "pearson",
                                      db_flag = TRUE, grouping_add_all = TRUE)
df_stat <- pairwise_stats(df = df_rho, db_flag = TRUE, partition_vars = "tissue", pairwise_var = "icr_conf")  
sum_df <- df_rho %>%  group_by(tissue, icr_conf) %>%  
  summarize(mean = mean(rho), sd = sd(rho), sem = sd(rho)/sqrt(length(rho)), n = n())
plot_agreement_by_confidence(df_rho, sum_df, df_stat, x_var = "icr_conf",
                                         subset_level = "huvec_cd14", 
                                         xlab_str = "ICR Confidence", 
                                         ylab_str = "Rho (ICR Betas)",
                                         output_dir_path = output_dir_path,
                                         export_filename = paste(marg_dir, "_ICR_Rho_HUVEC-CD14.jpg"),
                             fig_dim = c(4,2), db_flag = TRUE)
  
  
# Figure 2A_2) Agreement of ICRs with and without zinc finger  ####
#_______________________________________________________________________________

# Calculate agreement between beta values of ICRs with and without zinc finger sites
df_rho <- marginal_grouped_agreement(imp_long = icr_long, marg_dir = marg_dir, 
                                     grouping_var = "is_icr_zinc", agreement = "pearson",
                                     db_flag = TRUE, grouping_add_all = FALSE)
# Pairwise and 1 sample t-tests
df_stat <- pairwise_stats(df = df_rho, db_flag = TRUE, partition_vars = "tissue", 
                          pairwise_var = "is_icr_zinc")  
# Summary states (mean, sd, sem)
sum_df <- df_rho %>% group_by(tissue, is_icr_zinc) %>% 
  summarize( mean = mean(rho), sd = sd(rho), sem = sd(rho)/sqrt(length(rho)),
             n=length(rho))
# Plot
plot_agreement_by_confidence(df_rho, sum_df, df_stat,  x_var = "is_icr_zinc",
                             subset_level = "huvec_cd14", 
                             xlab_str = ".", 
                             ylab_str = "Rho (ICR Betas)",
                             output_dir_path = output_dir_path,
                             export_filename = paste(marg_dir, "_ICR_Rho_HUVEC-CD14_zinc_finger.jpg"), 
                             fig_dim = c(2,2), db_flag = FALSE)





# Figure 2B) Agreement of cpgs between HUVECS and CD14            ####
#_______________________________________________________________________________

df_rho <- marginal_grouped_agreement(imp_long = cpg_long, marg_dir = marg_dir, 
                                     grouping_var = "icr_conf", agreement = "pearson",
                                     db_flag = TRUE, grouping_add_all = TRUE)
df_stat <- pairwise_stats(df = df_rho, db_flag = TRUE, partition_vars = "tissue", pairwise_var = "icr_conf")  
sum_df <- df_rho %>% group_by(tissue, icr_conf) %>%  
  summarize(mean = mean(rho), sd = sd(rho), sem = sd(rho)/sqrt(length(rho)), n = n())
plot_agreement_by_confidence(df_rho, sum_df, df_stat, x_var = "icr_conf",
                             subset_level = "huvec_cd14", 
                             xlab_str = "ICR Confidence", 
                             ylab_str = "Rho (CpG Betas)",
                             output_dir_path = output_dir_path,
                             export_filename = paste(marg_dir, "_CPG_Rho_HUVEC-CD14.jpg"),
                             fig_dim = c(4,2))


# Figure 2B_2) Agreement of ICRs with and without zinc finger  ####
#_______________________________________________________________________________

# Calculate agreement between beta values of ICRs with and without zinc finger sites
df_rho <- marginal_grouped_agreement(imp_long = cpg_long, marg_dir = marg_dir, 
                                     grouping_var = "is_icr_zinc", agreement = "pearson",
                                     db_flag = TRUE, grouping_add_all = FALSE)
# Pairwise and 1 sample t-tests
df_stat <- pairwise_stats(df = df_rho, partition_vars = "tissue", pairwise_var = "is_icr_zinc")  
# Summary states (mean, sd, sem)
sum_df <- df_rho %>% group_by(tissue, is_icr_zinc) %>%  
  summarize(mean = mean(rho), sd = sd(rho),sem = sd(rho)/sqrt(length(rho)), n=length(rho))
# plot
plot_agreement_by_confidence(df_rho, sum_df, df_stat,  x_var = "is_icr_zinc",
                             subset_level = "huvec_cd14", 
                             xlab_str = ".", 
                             ylab_str = "Rho of ICR Betas",
                             output_dir_path = output_dir_path,
                             export_filename = paste(marg_dir, "_CPG_Rho_HUVEC-CD14_zinc_finger.jpg"), 
                             fig_dim = c(2,2), db_flag = FALSE)



# Figure 2C) How consistent are the Cpgs shared with 850K between HUVECS and CD14   ####
#_______________________________________________________________________________
# Calculate rho on continuous beta values
unq_cpg_850k <- readRDS(file = paste0(support_data_path, "/unq_cpg_850k.rds"))

df_rho <- marginal_grouped_agreement(imp_long = cpg_long %>% filter(cpg_id %in% unq_cpg_850k),
                                     marg_dir = marg_dir, 
                                     grouping_var = "icr_conf", agreement = "pearson",
                                     db_flag = TRUE, grouping_add_all = TRUE)
df_stat <- pairwise_stats(df = df_rho,db_flag = TRUE, partition_vars = "tissue", pairwise_var = "icr_conf")  
sum_df <- df_rho %>% group_by(tissue, icr_conf) %>% 
  summarize(mean = mean(rho), sd = sd(rho), sem = sd(rho)/sqrt(length(rho)), n = n())
plot_agreement_by_confidence(df_rho, sum_df, df_stat, x_var = "icr_conf",
                             subset_level = "huvec_cd14", 
                             xlab_str = "ICR Confidence", 
                             ylab_str = "Rho of 850K CpG Betas",
                             output_dir_path = output_dir_path,
                             export_filename = paste(marg_dir, "_CPG_850k_Rho_HUVEC-CD14.jpg"),
                             fig_dim = c(4,2))


# Calculate agreement between beta values of ICRs with and without zinc finger sites
df_rho <- marginal_grouped_agreement(imp_long = cpg_long %>% filter(cpg_id %in% unq_cpg_850k),
                                     marg_dir = marg_dir, 
                                     grouping_var = "is_icr_zinc", agreement = "pearson",
                                     db_flag = TRUE, grouping_add_all = FALSE)
# Pairwise and 1 sample t-tests
df_stat <- pairwise_stats(df = df_rho, partition_vars = "tissue", pairwise_var = "is_icr_zinc")
# Summary states (mean, sd, sem)
sum_df <- df_rho %>% group_by(tissue, is_icr_zinc) %>%  
  summarize(mean = mean(rho), sd = sd(rho),sem = sd(rho)/sqrt(length(rho)), n=length(rho))
# plot
plot_agreement_by_confidence(df_rho, sum_df, df_stat,  x_var = "is_icr_zinc",
                             subset_level = "huvec_cd14", 
                             xlab_str = ".", 
                             ylab_str = "Rho of 850K CpG Betas",
                             output_dir_path = output_dir_path,
                             export_filename = paste(marg_dir, "_CPG_850K_Rho_HUVEC-CD14_zinc_finger.jpg"), 
                             fig_dim = c(2,2), db_flag = FALSE)









# Figure 2D) How consistent are the Cpgs shared with CORSIVS between HUVECS and CD14   ####
#_______________________________________________________________________________
# Calculate rho on continuous beta values

df_rho <- marginal_grouped_agreement(imp_long = cpg_long %>% filter(is_corsiv),
                                     marg_dir = marg_dir, 
                                     grouping_var = "icr_conf", agreement = "pearson",
                                     db_flag = TRUE, grouping_add_all = TRUE)
df_stat <- pairwise_stats(df = df_rho,db_flag = TRUE, partition_vars = "tissue", pairwise_var = "icr_conf")  
sum_df <- df_rho %>% group_by(tissue, icr_conf) %>% 
  summarize(mean = mean(rho), sd = sd(rho), sem = sd(rho)/sqrt(length(rho)), n = n())
plot_agreement_by_confidence(df_rho, sum_df, df_stat, x_var = "icr_conf",
                             subset_level = "huvec_cd14", 
                             xlab_str = "ICR Confidence", 
                             ylab_str = "Rho of CorSIV CpG Betas",
                             output_dir_path = output_dir_path,
                             export_filename = paste(marg_dir, "_CPG_CORSIV_Rho_HUVEC-CD14.jpg"),
                             fig_dim = c(4,2))


# Calculate agreement between beta values of ICRs with and without zinc finger sites
df_rho <- marginal_grouped_agreement(imp_long = cpg_long %>% filter(is_corsiv),
                                     marg_dir = marg_dir, 
                                     grouping_var = "is_icr_zinc", agreement = "pearson",
                                     db_flag = TRUE, grouping_add_all = FALSE)
# Pairwise and 1 sample t-tests
df_stat <- pairwise_stats(df = df_rho, partition_vars = "tissue", pairwise_var = "is_icr_zinc")
# Summary states (mean, sd, sem)
sum_df <- df_rho %>% group_by(tissue, is_icr_zinc) %>%  
  summarize(mean = mean(rho), sd = sd(rho),sem = sd(rho)/sqrt(length(rho)), n=length(rho))
# plot
plot_agreement_by_confidence(df_rho, sum_df, df_stat,  x_var = "is_icr_zinc",
                             subset_level = "huvec_cd14", 
                             xlab_str = ".", 
                             ylab_str = "Rho of CorSIV ICR Betas",
                             output_dir_path = output_dir_path,
                             export_filename = paste(marg_dir, "_CPG_CORSIV_Rho_HUVEC-CD14_zinc_finger.jpg"), 
                             fig_dim = c(2,2), db_flag = FALSE,
                             include_neg_y_axis = F)













# Agreement between ICRs between CD14 and HUVECs, with design score threshold  ####
#_______________________________________________________________________________

# Calculate rho with cpg beta values with all probes versus high design score probes
marg_dir = "patient"
# All probes
df_rho <- marginal_grouped_agreement(imp_long = cpg_long,# %>% filter(cpg_id %in% unq_cpg_850k),
                                     marg_dir = marg_dir, 
                                     grouping_var = "icr_conf", agreement = "pearson",
                                     db_flag = TRUE, grouping_add_all = TRUE)

# Reprocess beta values
# 2) Filter probes that are not mapped and discard poor signal
# Output: rows: probe_id     x     columns: patient
data_tri$filt_probe_beta2 <-
  filter_probes(probe_beta = data_tri$probe_beta, discard_unmapped_probes = TRUE ,
                max_sig_pval = 0.2, set_failed_betas_na = FALSE, 
                max_probe_fail_rate = 0.25, min_design_score = .5, 
                discard_failed_probes = TRUE, db_flag = FALSE, max_patient_fail_rate = 1,
                discard_failed_patients = F)
#3  Convert probe beta matrix to a cpg beta matrix
# Output: rows: cpg_id     x     columns: patient
data_tri$cpg_beta2 <- convert_probes_to_cpgs(data_tri$filt_probe_beta2, quantile_norm = FALSE, 
                                   db_flag = TRUE, smooth_adj_cpgs = FALSE)

cpg_long2 <- beta_mat_to_long(data_tri$cpg_beta2$cpg_beta_df, filt_study_data, db_flag = T)

# Just high design score probes
df_rho2 <- marginal_grouped_agreement(imp_long = cpg_long2,
                                      marg_dir = marg_dir, 
                                     grouping_var = "icr_conf", agreement = "pearson",
                                     db_flag = TRUE, grouping_add_all = TRUE)


df_rho_both <- rbind(df_rho %>% mutate(design_score = FALSE),
                     df_rho2 %>% mutate(design_score = TRUE))
# Plot of agreement partitioned by icr confidence
df_stat <- rbind(pairwise_stats(df = df_rho, partition_vars = "tissue", pairwise_var = "icr_conf") %>% 
                   mutate(design_score = FALSE), 
                 pairwise_stats(df = df_rho2, partition_vars = "tissue", pairwise_var = "icr_conf") %>%
                   mutate(design_score = TRUE))
sum_df <- df_rho_both %>% group_by(tissue, icr_conf, design_score) %>%  
  summarize(mean = mean(rho), sd = sd(rho), sem = sd(rho)/sqrt(length(rho)))
# plot
plot_agreement_by_confidence(df_rho_both, sum_df, df_stat,  x_var = "icr_conf",
                             fill_var = "design_score",
                             subset_level = "huvec_cd14", 
                             xlab_str = "ICR Confidence", 
                             ylab_str = "Rho of CpG Betas",
                             output_dir_path = output_dir_path,
                             export_filename = paste(marg_dir, "_ICR_Rho_Design_SCore_HUVEC-CD14.jpg"), 
                             fig_dim = c(4,2), db_flag = TRUE)

# Fit model to see interactions with
# rho ~ icr_conf + design_score
df_model <- df_rho_both %>% subset(tissue == "huvec_cd14")
design_model <- glm(formula = "rho ~ icr_conf + design_score", family = "gaussian", data = df_model %>% filter(icr_conf != "All"))
design_summary <- summary(design_model)
coeffs<-design_summary$coefficients







