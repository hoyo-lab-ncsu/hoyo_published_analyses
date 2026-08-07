

# Use original CPG beta from original


# # Get all ICRs from CD14 results
# df_imp_sig <- df_imp_all %>% filter()
# # Significant ICRs
# df_imp_sig <- filter(select(df_imp_all, -Formula), ADJ_P_VAL < 0.05)

# Extract all cd14 cpg sites with significant results
df_imp_sig <- select(df_imp_all, -Formula) %>% filter((tissue == "cd14") & (ADJ_P_VAL < 0.05))

# Add gene and zinf finger metadata
df_imp_sig <- cbind(df_imp_sig,add_metadata_to_imp_sites(df_imp_sig$Variable, imp_type = "cpg"))
df_imp_sig$icr_name = str_replace(df_imp_sig$icr_id, "_"," ")
df_imp_sig$group_label <- 
  str_to_title(apply(X = cbind(str_split(df_imp_sig$group,"_", simplify = TRUE),
                               df_imp_sig$tissue)[,c(1,3,2)], MARGIN = 1, FUN = paste,
                     collapse=" "))


# Get list of unique associated ICRs with CD14
icr_compare_list <- unique(df_imp_sig$icr_id)
icr_whitelist <- rownames(data_tri$icr_beta$icr_beta_df) %in% icr_compare_list

cpg_compare_list <- unique(df_imp_sig$cpg_id)
cpg_whitelist <- rownames(data_tri$cpg_beta$cpg_beta_df) %in% cpg_compare_list

# Export the ICRs and CPGs from the association analysis for use in other parts of paper
save(icr_compare_list, cpg_compare_list, file = paste0(output_dir_path, "/cd14_ICR_correlation_hits.rda"))

# Repeat correlation figures with this subset of ICRs
#_______________________________________________________________________________


# Figure 2A) How consistent are the ICRs between HUVECS and CD14            ####
#_______________________________________________________________________________

# Calculate rho on continuous beta values
df_rho <- pcor_beta(df_beta = data_tri$icr_beta$icr_beta_df[icr_whitelist,], 
                    conf_group = left_join(x = data.frame(icr_id = rownames(data_tri$icr_beta$icr_beta_df[icr_whitelist,])), 
                                           y = df_icr_conf %>% select("icr_id", "icr_conf"),
                                           by = join_by("icr_id"), keep = F) %>% pull("icr_conf"), 
                    df_idat_sources = df_idat_sources, db_flag = TRUE)
df_stat <- pcorr_stats(df = df_rho)
sum_df <- subset(df_rho, tissue == "HUVEC Vs.\nCD14") %>% 
  group_by(icr_conf) %>%  summarize(mean = mean(rho), sd = sd(rho), sem = sd(rho)/sqrt(length(rho)))
# Plot of agreement partitioned by icr confidence
g0 <- ggplot(data = subset(df_rho, tissue == "HUVEC Vs.\nCD14"), 
             aes(x = icr_conf, y = rho)) +
  geom_violin(fill = "grey90", color = NA) +
  geom_errorbar(data = sum_df, width = 0.2, aes(
    x = icr_conf, y = mean, ymin = mean-sd, ymax = mean+sd)) +
  geom_point(data = sum_df, shape = 3, size = 5, aes(x = icr_conf, y = mean)) +
  geom_text(data = df_stat %>% subset((tissue == "HUVEC Vs.\nCD14") & (p_val_type == "pvalue_0")),
            aes(x = group1, label = pstr, y = rep(0.90,4)),
            color = "black") +
  geom_text(data = df_stat %>% subset((tissue == "HUVEC Vs.\nCD14") & (p_val_type == "pvalue_comp")),
            aes(x = group1, label = sig_strs, y = rep(1.0,4)),
            color = "black", size = 3.5) +
  xlab("ICR Confidence") + ylab("Rho of ICR Betas") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1.03)) 
# scale_fill_manual(values=c(icr_conf_colors,"grey"))
g0
save_plot(paste0(output_dir_path,"/subsetICR_sig_CD14_HUVEC_v_CD14_confidence.png"), g0,
          base_width = 4, base_height = 2)




# Figure 2A_2) What is agreement between ICRs with and without zinc finger    ####
df_rho <- pcor_beta(df_beta = data_tri$icr_beta$icr_beta_df[icr_whitelist,], 
                    conf_group = left_join(x = data.frame(icr_id = rownames(data_tri$icr_beta$icr_beta_df[icr_whitelist,])), 
                                           y = df_icr_conf %>% select("icr_id", "zinc_finger"),
                                           by = join_by("icr_id"), keep = F) %>% pull("zinc_finger"), 
                    df_idat_sources = df_idat_sources)
df_rho <- filter(df_rho, n !=0) 
df_rho$icr_conf <- factor(as.character(df_rho$icr_conf), 
                          labels =  c("-Zinc Finger", "+Zinc Finger") )
df_stat <- pcorr_stats(df = df_rho, comp_levels = c("-Zinc Finger", "+Zinc Finger"))
sum_df <- subset(df_rho, tissue == "HUVEC Vs.\nCD14") %>% 
  group_by(icr_conf) %>%  summarize(mean = mean(rho), sd = sd(rho),sem = sd(rho)/sqrt(length(rho)), n=length(rho))
# What is agree ment like for all 
g0 <- ggplot(data = subset(df_rho, tissue == "HUVEC Vs.\nCD14"), 
             aes(x = icr_conf, y = rho)) +
  geom_violin(fill = "grey90", color = NA) +
  geom_errorbar(data = sum_df, width = 0.2, aes(
    x = icr_conf, y = mean, ymin = mean-sd, ymax = mean+sd)) +
  geom_point(data = sum_df, shape = 3, size = 5, aes(x = icr_conf, y = mean)) +
  geom_text(data = df_stat %>% subset((tissue == "HUVEC Vs.\nCD14") & (p_val_type == "pvalue_0")),
            aes(x = group1, label = pstr, y = rep(0.90,2)),
            color = "black") +
  geom_text(data = df_stat %>% subset((tissue == "HUVEC Vs.\nCD14") & (p_val_type == "pvalue_comp")),
            aes(x = group1, label = sig_strs, y = rep(1.0,2)),
            color = "black", size = 3.5) +
  xlab(".") + ylab("Rho of ICR Betas") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1)) 
# scale_fill_manual(values=c(icr_conf_colors,"grey"))
g0
save_plot(paste0(output_dir_path,"/subsetICR_sig_CD14_HUVEC_v_CD14_zinc_finger.png"), g0,
          base_width = 2, base_height = 2)









# Figure 2B) How consistent are the Cpgs between HUVECS and CD14            ####
#_______________________________________________________________________________

# Calculate rho on continuous beta values
df_rho <- pcor_beta(df_beta = data_tri$cpg_beta$cpg_beta_df[cpg_whitelist,], 
                    conf_group = left_join(x = data.frame(cpg_id = rownames(data_tri$cpg_beta$cpg_beta_df[cpg_whitelist,])), 
                                           y = df_cpg_conf %>% select("cpg_id", "icr_conf"),
                                           by = join_by("cpg_id"), keep = F) %>% pull("icr_conf"), 
                    df_idat_sources = df_idat_sources, db_flag = TRUE)
df_stat <- pcorr_stats(df = df_rho)
sum_df <- subset(df_rho, tissue == "HUVEC Vs.\nCD14") %>% 
  group_by(icr_conf) %>%  summarize(mean = mean(rho), sd = sd(rho))
# What is agree ment like for all 
g0 <- ggplot(data = subset(df_rho, tissue == "HUVEC Vs.\nCD14"), 
             aes(x = icr_conf, y = rho)) +
  geom_violin(fill = "grey90", color = NA) +
  geom_errorbar(data = sum_df, width = 0.2, aes(
    x = icr_conf, y = mean, ymin = mean-sd, ymax = mean+sd)) +
  geom_point(data = sum_df, shape = 3, size = 5, aes(x = icr_conf, y = mean)) +
  geom_text(data = df_stat %>% subset((tissue == "HUVEC Vs.\nCD14") & (p_val_type == "pvalue_0")),
            aes(x = group1, label = pstr, y = rep(0.90,4)),
            color = "black") +
  geom_text(data = df_stat %>% subset((tissue == "HUVEC Vs.\nCD14") & (p_val_type == "pvalue_comp")),
            aes(x = group1, label = sig_strs, y = rep(1.00,4)),
            color = "black", size = 3.5) +
  xlab("ICR Confidence") + ylab("Rho of CpG Betas") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1)) 
# scale_fill_manual(values=c(icr_conf_colors,"grey"))
g0
save_plot(paste0(output_dir_path,"/subsetCPG_sig_CD14_HUVEC_v_CD14_confidence.png"), g0,
          base_width = 4, base_height = 2)


# Figure 2B_2) What is agreement between ICRs with and without zinc finger    ####
df_rho <- pcor_beta(df_beta = data_tri$cpg_beta$cpg_beta_df, 
                    conf_group = left_join(x = data.frame(cpg_id = rownames(data_tri$cpg_beta$cpg_beta_df[cpg_whitelist,])), 
                                           y = df_cpg_conf %>% select("cpg_id", "zinc_finger"),
                                           by = join_by("cpg_id"), keep = F) %>% pull("zinc_finger"), 
                    df_idat_sources = df_idat_sources)
df_rho <- filter(df_rho, n !=0) 
df_rho$icr_conf <- factor(as.character(df_rho$icr_conf), 
                          labels =  c("-Zinc Finger", "+Zinc Finger") )
df_stat <- pcorr_stats(df = df_rho, comp_levels = c("-Zinc Finger", "+Zinc Finger"))
sum_df <- subset(df_rho, tissue == "HUVEC Vs.\nCD14") %>% 
  group_by(icr_conf) %>%  summarize(mean = mean(rho), sd = sd(rho))
# What is agree ment like for all 
g0 <- ggplot(data = subset(df_rho, tissue == "HUVEC Vs.\nCD14"), 
             aes(x = icr_conf, y = rho)) +
  geom_violin(fill = "grey90", color = NA) +
  geom_errorbar(data = sum_df, width = 0.2, aes(
    x = icr_conf, y = mean, ymin = mean-sd, ymax = mean+sd)) +
  geom_point(data = sum_df, shape = 3, size = 5, aes(x = icr_conf, y = mean)) +
  geom_text(data = df_stat %>% subset((tissue == "HUVEC Vs.\nCD14") & (p_val_type == "pvalue_0")),
            aes(x = group1, label = pstr, y = rep(0.90,2)),
            color = "black") +
  geom_text(data = df_stat %>% subset((tissue == "HUVEC Vs.\nCD14") & (p_val_type == "pvalue_comp")),
            aes(x = group1, label = sig_strs, y = rep(1.0,2)),
            color = "black", size = 3.5) +
  xlab(".") + ylab("Rho of CpG Betas") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1)) 
# scale_fill_manual(values=c(icr_conf_colors,"grey"))
g0
save_plot(paste0(output_dir_path,"/subsetCPG_sig_CD14_HUVEC_v_CD14_zinc_finger.png"), g0,
          base_width = 2, base_height = 2)
