


# For each ICR and Response combination, record info, and:
#  1) Minimum adjusted p-value across cpg sites (group_by, summarize)
#  2) Max magnitude of model estimate
df_icr_sig <- df_imp_sig_out %>%
  group_by(icr_name, group_label) %>%
  summarise(ADJ_P_VAL = min(ADJ_P_VAL), meanADJ_P_VAL = mean(ADJ_P_VAL),
            Estimate = max(abs(Estimate)) * sign(Estimate[which.max(abs(Estimate))])) %>% 
  merge(y = select(df_imp_sig_out,-c(Response,group_label, Estimate, ADJ_P_VAL)) %>% distinct(), 
        by = "icr_name") %>% arrange(desc(group_label))

df_icr_sig$Freq <- sapply( df_icr_sig$icr_name, function(x) sum(x==df_icr_sig$icr_name))
df_icr_sig$Responses.Shared <- sapply( df_icr_sig$icr_name, function(x) 
  paste(df_icr_sig$group_label[x==df_icr_sig$icr_name], collapse = ", "))


# # For each ICR and Response combination, record info, and:
# #  1) Minimum adjusted p-value across cpg sites
# #  2) Max magnitude of model estimate
# df_icr_sig <- df_imp_sig_out %>%
#   group_by(model_group, icr_id, Response) %>%
#   summarise(ADJ_P_VAL = min(ADJ_P_VAL), 
#             meanADJ_P_VAL = mean(ADJ_P_VAL),
#             Estimate = max(abs(Estimate)) * sign(Estimate[which.max(abs(Estimate))])) %>% 
#   merge(y = select(df_imp_sig_out,-c(model_group, Response, Estimate, ADJ_P_VAL)) %>% distinct(), 
#         by = "icr_id") %>% arrange(desc(Response))
# df_icr_sig$model_response <- paste0( df_icr_sig$Response, " (",  
#                                      substr(df_icr_sig$model_group,1 ,1), ")")
# # Only use data from Male and Female models
# df_icr_sig <- df_icr_sig %>% filter(model_group != "All")
# df_icr_sig$Freq <- sapply( df_icr_sig$icr_id, function(x) sum(x==df_icr_sig$icr_id))

# Order ICRs by their frequency
df_icr_sig <- df_icr_sig %>% arrange(desc(Freq), meanADJ_P_VAL)

df_icr_sig_wide <- df_icr_sig %>% select(icr_name, Estimate, group_label) %>% 
  pivot_wider(names_from = group_label, values_from = Estimate) %>%
  column_to_rownames("icr_name") %>% 
  select(sort(names(.)))

# COnvert to matrix, set NaNs to 0 since heatmap requires all numeric
mat_min_p <- as.matrix(df_icr_sig_wide) %>% sign()
# Sort by rows by number of non-NA values
mat_min_p[is.na(mat_min_p)] <- 0
rownames(mat_min_p) <- str_replace(rownames(mat_min_p), "_", "-")
colnames(mat_min_p) <- sub("(\\s.*?)\\s", "\\1\n", colnames(mat_min_p))


gh <- pheatmap::pheatmap(mat_min_p[0:100,],  col = c("#e41a1c", "white","#377eb8"),
                         legend_breaks = 0:2, legend_labels = c("", "", ""), 
                         show_rownames = T,  treeheight_row = 0, treeheight_col = 0,
                         legend = FALSE,border_color = "black",cluster_cols = F, cluster_rows = F,
                         fontsize_row = 7, fontsize_col = 7)
save_plot(filename = paste0(output_dir_path, "/sign_heatmap_1.png"),
          plot = gh, base_height = 9, base_width = 3.25)


gh <- pheatmap::pheatmap(mat_min_p[101:200,],  col = c("#e41a1c", "white","#377eb8"),
                         legend_breaks = 0:2, legend_labels = c("", "", ""), 
                         show_rownames = T,  treeheight_row = 0, treeheight_col = 0,
                         legend = FALSE,border_color = "black",cluster_cols = F, cluster_rows = F,
                         fontsize_row = 7, fontsize_col = 7)
save_plot(filename = paste0(output_dir_path, "/sign_heatmap_2.png"),
          plot = gh, base_height = 9, base_width = 3.25)


gh <- pheatmap::pheatmap(mat_min_p[201:300,],  col = c("#e41a1c", "white","#377eb8"),
                         legend_breaks = 0:2, legend_labels = c("", "", ""), 
                         show_rownames = T,  treeheight_row = 0, treeheight_col = 0,
                         legend = FALSE,border_color = "black",cluster_cols = F, cluster_rows = F,
                         fontsize_row = 7, fontsize_col = 7)
save_plot(filename = paste0(output_dir_path, "/sign_heatmap_3.png"),
          plot = gh, base_height = 9, base_width = 3.25)


gh <- pheatmap::pheatmap(mat_min_p[301:400,],  col = c("#e41a1c", "white","#377eb8"),
                         legend_breaks = 0:2, legend_labels = c("", "", ""), 
                         show_rownames = T,  treeheight_row = 0, treeheight_col = 0,
                         legend = FALSE, border_color = "black",cluster_cols = F, cluster_rows = F,
                         fontsize_row = 7, fontsize_col = 7)
save_plot(filename = paste0(output_dir_path, "/sign_heatmap_4.png"),
          plot = gh, base_height = 9, base_width = 3.25)

