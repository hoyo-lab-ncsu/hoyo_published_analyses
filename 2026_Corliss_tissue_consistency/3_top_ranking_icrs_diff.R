



# df_icr_conf


# Figure 3A) Top top ranking ICRS                                   ############
#_______________________________________________________________________________


icr_mean_diff <- pairwise_dist_beta(df_beta = data_tri$icr_beta$icr_beta_df,
                                    conf_group = df_icr_conf$icr_conf,
                                    df_idat_sources = df_idat_sources,
                                    dist_func = \(x,y) y-x)

xlabel = "Mean Beta Diff."
# xlabel = "Hemi-methylation Beta Dist."


icr_conf_colors = c("#fbb4ae", "#ccebc5", "#b3cde3")


# Matrix of distance between HUVEC vs CD14 betas for each patient
df_pdiff_all = as.data.frame(icr_mean_diff$HUVEC_CD14)


# Summarize the range of distances across all patients for each ICR 
df_pdiff_summary = data.frame( mean = rowMeans(df_pdiff_all), abs_mean = abs(rowMeans(df_pdiff_all)))

df_pdiff_summary$p25 <- sapply(1:nrow(df_pdiff_all),FUN = function(x) 
  quantile(as.matrix(df_pdiff_all[x,]),0.25))
df_pdiff_summary$p50 <- sapply(1:nrow(df_pdiff_all),FUN = function(x) 
  quantile(as.matrix(df_pdiff_all[x,]),0.50))
df_pdiff_summary$p75 <- sapply(1:nrow(df_pdiff_all),FUN = function(x) 
  quantile(as.matrix(df_pdiff_all[x,]),0.75))
df_pdiff_summary$sd <- sapply(1:nrow(df_pdiff_all),FUN = function(x) 
  sd(as.matrix(df_pdiff_all[x,])) )
df_pdiff_summary$sem <- df_pdiff_summary$sd/sqrt(nrow(df_pdiff_all))
# df_pdiff_summary$icr_conf = factor(add_metadata_to_imp_sites(
#   imp_ids = rownames(df_pdiff_summary),imp_type = "icr")$icr_conf,ordered = TRUE)


df_pdiff_summary <- cbind(df_pdiff_summary, add_metadata_to_imp_sites(imp_ids = rownames(df_pdiff_summary),imp_type = "icr"))
# Add column that denotes zinc status of each ICR (close to Zinf finger site)
# df_pdiff_summary$zinc =  add_metadata_to_imp_sites(imp_ids = rownames(df_pdiff_summary),imp_type = "icr")$is_icr_zinc 

  
# Sort entries by mean distance
sort_ind <- order(df_pdiff_summary$abs_mean, decreasing = FALSE)
df_pdiff_summary <- df_pdiff_summary[sort_ind,]
df_pdiff_all <- df_pdiff_all[sort_ind,]


# Add icr_id as a column from row names
df_pdiff_summary$icr_id <- str_replace(rownames(df_pdiff_summary), "_", " ")
df_pdiff_summary$icr_id <- factor(x = df_pdiff_summary$icr_id, 
                                  levels = df_pdiff_summary$icr_id, ordered = TRUE)
df_pdiff_summary$icr_name <- df_pdiff_summary$icr_id
# df_pdiff_summary$icr_confidence = factor(df_pdiff_summary$icr_conf, levels = c("1", "2", "3"), labels = c("High", "Medium", "Low"))


# Add closest gene to each ICR
# df_pdiff_summary <- merge(x = df_pdiff_summary, y = imp_whole %>%
#                             select(icr_name, ),
#                           by = "icr_name", all.x = TRUE, all.y = FALSE, sort = FALSE)

# # Data table exports                                                  ##########
# #_______________________________________________________________________________
# Making pretty version of table
export_table_diff <- select(df_pdiff_summary, c(icr_name, icr_conf, is_icr_zinc, Nearest.Transcript,
                                                     Distance.to.Nearest.Transcript, mean, abs_mean))
rownames(export_table_diff) = NULL
write_rds(x = export_table_diff, file = paste0(output_dir_path, "/ICR_ranks_diff_CD14_HUVEC.rds"))
colnames(export_table_diff) <- c("ICR ID",	"ICR Confidence",	"Zinc Finger Proximity",
                                      "Nearest Transcript",	"Transcript Distance",	"Mean Diff", "|Mean Diff|")


write_csv(x = export_table_diff, file = paste0(output_dir_path, "/ICR_ranks_diff_CD14_HUVEC.csv"))

# # If outptu data from CD14 stats analysis exist, then find top 10 and export to disk
# if (file.exists(paste0(output_dir_path, "/cd14_ICR_correlation_hits.rda"))) {
#   load(file = paste0(output_dir_path, "/cd14_ICR_correlation_hits.rda"))
#   # icr_compare_list, cpg_compare_list
#   
#   write.csv(x = export_table_diff[export_table_diff$`ICR ID` %in%
#                                          str_replace(icr_compare_list, "_", " "), ][1:10,], 
#             file = paste0(output_dir_path, "/top10_icrs_diff_from_cd14_study_subset.csv"))
#   
# }
# # save(icr_compare_list, cpg_compare_list, file = paste0(output_dir_path, "/cd14_ICR_correlation_hits.rda"))
# 
# 
# 
# 
# 
# # Output all ICRs sorted by both metrics               #########################
# #_______________________________________________________________________________
# # Find top ICRs by hemi-distance and diff
# df_shared <- select(df_pdist_summary, -c("p25", "p50", "p75", "sd", "sem"))
# df_shared <- df_shared %>% merge(y=select(df_pdiff_summary, c("icr_name", "mean")), by = "icr_name", all.x = TRUE)
# 
# 
# df_top_shared <- select(df_pdist_summary, -c("p25", "p50", "p75", "sd", "sem"))
# df_top_shared <- df_top_shared %>% merge(y = select(df_pdiff_summary, 
#                                                     c("icr_name", "mean")), by = "icr_name", 
#                                          all.x = TRUE)
# 
# df_top_shared$mean <- rowMeans(preprocessCore::normalize.quantiles(as.matrix(
#   select(df_top_shared, c("mean.x", "mean.y"))), copy=TRUE))
# 
# df_top_shared <- df_top_shared[order(df_top_shared$mean),]
# 
# rownames(df_top_shared) <- 1: nrow(df_top_shared)
# # Export spreadsheet of top icrs by methylation distance
# export_table_hemi_both <-select(df_top_shared, c(icr_name, icr_confidence, zinc, Nearest.Transcript, 
#                                                  Distance.to.Nearest.Transcript, mean))
# colnames(export_table_hemi_both) <- c("ICR ID",	"ICR Confidence",	"Zinc Finger Proximity",
#                                       "Nearest Transcript",	"Transcript Distance",	"Mean Hemi-Dist")
# 
# # Write single xlsx for all of the individual ICR rankings
# WriteXLS::WriteXLS(list(export_table_hemi_both, export_table_hemi_dist, export_table_diff), 
#                    ExcelFileName = paste0(output_dir_path,"/ICR_Rankings_HUVEC-CD14.xlsx"), 
#                    SheetNames = list("Hemi_Both", "Hemi_Dist", "Hemi_Diff"))
# 
# 
# 
# # Export just the ICRS shared in the top 100 of both rankings   ################
# #_______________________________________________________________________________
# # Output shared ICRs in both tables
# shared_top_icrs <- as.character(df_pdist_summary$icr_name[1:100][
#   is.element(df_pdist_summary$icr_name[1:100], df_pdiff_summary$icr_name[1:100])])
# 
# write.csv(x = select(df_top_shared, c(icr_name, icr_confidence, zinc, Nearest.Transcript, 
#                                       Distance.to.Nearest.Transcript, mean)) %>% 
#             filter(icr_name %in% shared_top_icrs),
#           file = paste0(output_dir_path,"/top_icrs_both_top100_HUVEC-CD14.csv"), row.names = TRUE)


df_pdiff_summary$icr_conf <- as.factor(df_pdiff_summary$icr_conf)
# ICR Rank 1-50                                                    #############
#_______________________________________________________________________________
# Plot
g2a <- ggplot(data = df_pdiff_summary[1:50,], aes(y = icr_id, x = mean) ) + 
  geom_tile(aes(width = Inf, height = 1, fill = icr_conf)) + 
  geom_tile(data = filter(df_pdiff_summary[1:50,], is_icr_zinc  == TRUE), 
            aes(y = icr_id, x = abs(mean), width = Inf, height = 1), 
            fill = "transparent", color = "black", linewidth = .5) +
  geom_linerange(aes(xmin = mean-2.6*sem, xmax = mean+2.6*sem),
                 linewidth = 1) + geom_point() +
  scale_y_discrete(limits = rev) + ylab("") + xlab(xlabel)  +
  scale_fill_manual(values = icr_conf_colors[1:3]) +
  labs(fill="ICR Conf") + #coord_cartesian(xlim = c(0, .15)) +
  coord_cartesian(xlim = c(-0.05, 0.05)) +
  geom_vline(xintercept = 0) + 
  theme_classic() + theme(legend.position = "none")
g2a
save_plot(paste0(output_dir_path,"/ICR_diff_rho_HUVEC-CD14_1-50.jpg"), g2a,
          base_width = 3, base_height = 7)


# ICR Rank 51-100                                                    #############
#_______________________________________________________________________________
# Plot 51:100,    df_pdiff_summary[1000:nrow(df_pdiff_summary),]
g2b <- ggplot(data = df_pdiff_summary[51:100,], aes(y = icr_id, x = mean) ) + 
  geom_tile(aes(width = Inf, height = 1, fill = icr_conf)) + 
  geom_tile(data = filter(df_pdiff_summary[51:100,], is_icr_zinc  == TRUE), 
            aes(y = icr_id, x = abs(mean), width = Inf, height = 1), 
            fill = "transparent", color = "black", linewidth = .5) +
  geom_linerange(aes(xmin = mean-2.6*sem, xmax = mean+2.6*sem),
                 linewidth = 1) + geom_point() +
  scale_y_discrete(limits = rev) + ylab("") + xlab(xlabel)  +
  scale_fill_manual(values = icr_conf_colors[1:3]) +
  labs(fill="ICR Conf") + #coord_cartesian(xlim = c(0, .15)) +
  coord_cartesian(xlim = c(-0.05, 0.05)) +
  geom_vline(xintercept = 0) + 
  theme_classic() + theme(legend.position = "none")
g2b
save_plot(paste0(output_dir_path,"/ICR_diff_rho_HUVEC-CD14_51-100.jpg"), g2b,
          base_width = 3, base_height = 7)


# Make an order tiled image of all icrs                            #############
#_______________________________________________________________________________
# df_pdiff_all
# make r and c indices
df_tile <- select(df_pdiff_summary, icr_id, mean, icr_conf, is_icr_zinc)
df_tile$row <- rep(1:50, ceiling(nrow(df_tile)/50))[1:nrow(df_tile)]
df_tile$column <- rep(1:ceiling(nrow(df_tile)/50), each=50)[1:nrow(df_tile)]

g2c <- ggplot(df_tile, aes(x=column, y = row, fill= icr_conf)) + 
  scale_fill_manual(values=icr_conf_colors) + labs(fill="ICR Conf") +
  geom_tile() +  theme_classic() + theme(legend.position = "none") +
  geom_tile(data = filter(df_tile, is_icr_zinc  == TRUE), aes(x=column, y = row), fill = "transparent", 
            color = "black", linewidth = .25) +
  ylab("ICR Rank") + xlab("Column Starting Rank") + 
  scale_x_continuous(breaks = seq(1,21,by=2), expand=c(0,0), position = "top",
                     labels =  seq(0, 1000, by = 100))  + 
  scale_y_continuous(expand=c(0,0), trans = "reverse", breaks = c(1,10,20,30,40,50))
g2c
save_plot(paste0(output_dir_path,"/ICR_diff_rho_HUVEC-CD14_all.jpg"), g2c,
          base_width = 4, base_height = 2.5)






