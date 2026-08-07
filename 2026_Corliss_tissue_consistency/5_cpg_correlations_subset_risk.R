






# Pearson's rho split by tissue type, icr confidence                      ######
#_______________________________________________________________________________
df_cpg_pcor_cont <- pcor_beta(df_beta = cpg_beta$cpg_beta_df,
                              conf_group = df_cpg_conf$icr_conf,
                              df_idat_sources = df_idat_sources)
df_stat <- pcorr_stats(df = df_pcor_cont, pairwise= TRUE, group_by_var = "icr_conf",
                       pairwise_var = "tissue", comp_levels = levels(df_pcor_cont$tissue))
sum_df <- df_pcor_cont %>% group_by(tissue, icr_conf) %>%  
  summarize(mean = mean(rho), sd = sd(rho))

# Plotting
g1 <- ggplot(data = df_cpg_pcor_cont, aes(x = tissue, y = rho)) +
  geom_violin(fill = "grey90", color = NA) +
  geom_errorbar(data = sum_df, width = 0.2, aes(
    x = tissue, y = mean, ymin = mean-sd, ymax = mean+sd)) +
  geom_point(data = sum_df, shape = 3, size = 5, aes(x = tissue, y = mean)) +
  geom_text(data = filter(df_stat, p_val_type == "pvalue_0"),
            aes(x = tissue, y = 0.90, label = pstr),
            color = "black") +
  
  geom_text(data = filter(df_stat, p_val_type == "pvalue_comp"),
            aes(x = group1, label = sig_strs, y = 1.1),
            color = "black", size = 3) +
  facet_wrap(~icr_conf) +
  xlab("Tissue Comparison") + ylab("Rho (CpG Beta)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(0, 1.15)) + scale_y_continuous(breaks=seq(0,1,.25))
g1
save_plot(paste0(output_dir_path,"/4_Cpg_Rho_Tissues_All.jpg"), g1,
          base_width = 5, base_height = 2.5)

# Pearson's rho split by tissue type, icr confidence, and risk factor: BMI ######
#_______________________________________________________________________________
sum_df <- df_pcor_cont %>% group_by(tissue, icr_conf, mat_bmi_hi) %>%  
  summarize(mean = mean(rho), sd = sd(rho))
g2 <- ggplot(data = df_cpg_pcor_cont, aes(x = tissue, y = rho, fill = mat_bmi_hi)) +
  geom_violin(color = NA, alpha = 0.75) +
  geom_errorbar(data = sum_df, position = position_dodge(width = 0.9), aes(
    x = tissue, y = mean, ymin = mean-sd, ymax = mean+sd, group = mat_bmi_hi),
    width = 0.2) +
  geom_point(data = sum_df, position = position_dodge(width = 0.9), shape = 3, 
             size = 5, aes(x = tissue, y = mean, group = mat_bmi_hi)) +
  facet_wrap(~icr_conf) +
  scale_fill_manual(values = c("#41D3D7", "#FA9892")) +
  xlab("Tissue Comparison") + ylab("Rho (CpG Beta)") +
  guides(fill=guide_legend(title="Maternal BMI")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1))
g2
save_plot(paste0(output_dir_path,"/4_CPG_Rho_Tissues_Maternal BMI.jpg"), g2,
          base_width = 4.6, base_height = 2.5)
mod <- glm(data =  filter(df_pcor_cont, icr_conf != "All"), 
           formula = "rho ~ mat_aces_hi + tissue + icr_conf")
res<-summary(mod)
View(res$coefficients)
res

# Pearson's rho split by tissue type, icr confidence, and risk factor: aces ######
#_______________________________________________________________________________
#_______________________________________________________________________________
sum_df <- na.omit(df_pcor_cont) %>% group_by(tissue, icr_conf, mat_aces_hi) %>%  
  summarize(mean = mean(rho), sd = sd(rho))
g3 <- ggplot(data = na.omit(df_cpg_pcor_cont), aes(x = tissue, y = rho, fill = mat_aces_hi)) +
  geom_violin(color = NA, alpha = 0.75) +
  geom_errorbar(data = sum_df, position = position_dodge(width = 0.9), aes(
    x = tissue, y = mean, ymin = mean-sd, ymax = mean+sd, group = mat_aces_hi),
    width = 0.2) +
  geom_point(data = sum_df, position = position_dodge(width = 0.9), shape = 3, 
             size = 5, aes(x = tissue, y = mean, group = mat_aces_hi)) +
  facet_wrap(~icr_conf) +
  scale_fill_manual(values = c("#41D3D7", "#FA9892")) +
  xlab("Tissue Comparison") + ylab("Rho (CpG Beta)") +
  guides(fill=guide_legend(title="Maternal ACEs")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1))
g3
save_plot(paste0(output_dir_path,"/4_CPG_Rho_Tissues_Maternal ACEs.jpg"), g3,
          base_width = 4.6,base_height = 2.5)
mod <- glm(data =  filter(df_pcor_cont, icr_conf != "All"), formula = "rho ~ mat_aces_hi + tissue + icr_conf")
res<-summary(mod)
View(res$coefficients)
res

