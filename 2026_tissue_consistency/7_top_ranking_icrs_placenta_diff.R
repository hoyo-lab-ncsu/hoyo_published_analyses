



# df_icr_conf


# Figure 3A) Top top ranking ICRS                                   ############
#_______________________________________________________________________________


icr_mean_diff <- pairwise_dist_beta(df_beta = data_tri$icr_beta$icr_beta_df,
                                    conf_group = df_icr_conf$icr_conf,
                                    df_idat_sources = df_idat_sources,
                                    dist_func = \(x,y) abs(x-y))

xlabel = "Grand Mean Beta Diff."
# xlabel = "Hemi-methylation Beta Dist."


icr_conf_colors = c("#fbb4ae", "#ccebc5", "#b3cde3")


# Do this for each pair of tissues
# df_pdist_HUVEC_CD14 = as.data.frame(icr_mean_diff$HUVEC_CD14)
df_pdiff_HUVEC_PLAC = as.data.frame(icr_mean_diff$HUVEC_Placenta)
df_pdiff_CD14_PLAC = as.data.frame(icr_mean_diff$CD14_Placenta)

# df_pdiff_all <- (df_pdiff_HUVEC_PLAC + df_pdiff_CD14_PLAC)/2

# Summarize the range of distances across all patients for each ICR 
df_pdiff_summary = data.frame( mean = (rowMeans(df_pdiff_HUVEC_PLAC) + 
                                         rowMeans(df_pdiff_CD14_PLAC))/2)

# Standard deviations combined with sqrt of sum of squares
df_pdiff_summary$sd <- sapply(1:nrow(df_pdiff_HUVEC_PLAC), FUN = function(x) 
  sqrt( sd(as.matrix(df_pdiff_HUVEC_PLAC[x,]))^2 +
          sd(as.matrix(df_pdiff_CD14_PLAC [x,]))^2) )

df_pdiff_summary$sem <- df_pdiff_summary$sd/sqrt(nrow(df_pdiff_HUVEC_PLAC))
df_pdiff_summary$icr_conf = factor(df_icr_conf$icr_conf, levels = c(1,2,3), ordered = TRUE)

# Add column that denotes zinc status of each ICR (close to Zinf finger site)
df_pdiff_summary$zinc = as.numeric(str_replace(rownames(df_pdiff_summary),"ICR_", "")) %in% 
  zinc_finger_icrs


# Sort entries by mean distance
sort_ind <- order(df_pdiff_summary$mean, decreasing = FALSE)
df_pdiff_summary <- df_pdiff_summary[sort_ind,]
# df_pdiff_all <- df_pdiff_all[sort_ind,]


# Add icr_id as a column from row names
df_pdiff_summary$icr_id <- str_replace(rownames(df_pdiff_summary), "_", " ")
df_pdiff_summary$icr_id <- factor(x = df_pdiff_summary$icr_id, 
                                  levels = df_pdiff_summary$icr_id, ordered = TRUE)
df_pdiff_summary$icr_name <- df_pdiff_summary$icr_id
df_pdiff_summary$icr_confidence = factor(df_pdiff_summary$icr_conf, levels = c("1", "2", "3"),
                                         labels = c("High", "Medium", "Low"))

# Add closest gene to each ICR
df_pdiff_summary <- merge(x = df_pdiff_summary, y = imp_whole %>%
                            select(icr_name, Nearest.Transcript, Distance.to.Nearest.Transcript),
                          by = "icr_name", all.x = TRUE, all.y = FALSE, sort = FALSE)

# Data table exports                                                  ##########
#_______________________________________________________________________________
# Making pretty version of table 
export_table_hemi_diff <- select(df_pdiff_summary, c(icr_name, icr_confidence, zinc, Nearest.Transcript,
                                                     Distance.to.Nearest.Transcript, mean))
colnames(export_table_hemi_diff) <- c("ICR ID",	"ICR Confidence",	"Zinc Finger Proximity",
                                      "Nearest Transcript",	"Transcript Distance",	"Mean Hemi-Dist")




# Output all ICRs sorted by both metrics               #########################
#_______________________________________________________________________________
# Find top ICRs by hemi-distance and diff
df_shared <- select(df_pdist_summary, -c("sd", "sem"))
df_shared <- df_shared %>% merge(y=select(df_pdiff_summary, c("icr_name", "mean")), by = "icr_name", all.x = TRUE)


df_top_shared <- select(df_pdist_summary, -c("sd", "sem"))
df_top_shared <- df_top_shared %>% merge(y = select(df_pdiff_summary, c(
  "icr_name", "mean")), by = "icr_name",
  all.x = TRUE)


df_top_shared$mean <- rowMeans(preprocessCore::normalize.quantiles(as.matrix(
  select(df_top_shared, c("mean.x", "mean.y"))), copy=TRUE))

df_top_shared <- df_top_shared[order(df_top_shared$mean),]

rownames(df_top_shared) <- 1: nrow(df_top_shared)
# Export spreadsheet of top icrs by methylation distance
# Export spreadsheet of top icrs by methylation distance
export_table_hemi_both <-select(df_top_shared, c(icr_name, icr_confidence, zinc, Nearest.Transcript, 
                                                 Distance.to.Nearest.Transcript, mean))
colnames(export_table_hemi_both) <- c("ICR ID",	"ICR Confidence",	"Zinc Finger Proximity",
                                      "Nearest Transcript",	"Transcript Distance",	"Mean Hemi-Dist")

# Write single xlsx for all of the individual ICR rankings
WriteXLS::WriteXLS(list(export_table_hemi_both, export_table_hemi_dist, export_table_hemi_diff), 
                   ExcelFileName = paste0(output_dir_path,"/ICR_Rankings_Placenta.xlsx"), 
                   SheetNames = list("Hemi_Both", "Hemi_Dist", "Hemi_Diff"))




# Export just the ICRS shared in the top 100 of both rankings   ################
#_______________________________________________________________________________

# Output shared ICRs in both tables
shared_top_icrs <- as.character(df_pdist_summary$icr_name[1:100][
  is.element(df_pdist_summary$icr_name[1:100], df_pdiff_summary$icr_name[1:100])])


write.csv(x = select(df_top_shared, c(icr_name, icr_confidence, zinc, Nearest.Transcript, 
                                      Distance.to.Nearest.Transcript, mean)) %>% 
            filter(icr_name %in% shared_top_icrs),
          file = paste0(output_dir_path,"/top_icrs_both_top100_placenta.csv"), row.names = TRUE)






# ICR Rank 1-50
#_______________________________________________________________________________
# Plot
g2a <- ggplot(data = df_pdiff_summary[1:50,], aes(y = icr_id, x = abs(mean)) ) + 
  geom_tile(aes(width = Inf, height = 1, fill = icr_conf)) + 
  geom_tile(data = filter(df_pdiff_summary[1:50,], zinc == TRUE), 
            aes(y = icr_id, x = abs(mean), width = Inf, height = 1), 
            fill = "transparent", color = "black", linewidth = .5) +
  geom_linerange(aes(xmin = abs(mean)-2.6*sem, xmax = abs(mean)+2.6*sem),
                 linewidth = 1) + geom_point() +
  scale_y_discrete(limits = rev) + ylab("") + xlab(xlabel)  +
  scale_fill_manual(values = icr_conf_colors[1:3]) +
  labs(fill="ICR Conf") + coord_cartesian(xlim = c(0, .15)) +
  theme_classic() + theme(legend.position = "none")
g2a
save_plot(paste0(output_dir_path,"/top_icrs_hemi_diff_placenta_1-50.jpg"), g2a,
          base_width = 3, base_height = 7)


# ICR Rank 51-100
#_______________________________________________________________________________
# Plot
g2b <- ggplot(data = df_pdiff_summary[51:100,], aes(y = icr_id, x = abs(mean)) ) + 
  geom_tile(aes(width = Inf, height = 1, fill = icr_conf)) + 
  geom_tile(data = filter(df_pdiff_summary[51:100,], zinc == TRUE), 
            aes(y = icr_id, x = abs(mean), width = Inf, height = 1), 
            fill = "transparent", color = "black", linewidth = .5) +
  geom_linerange(aes(xmin = abs(mean)-2.6*sem, xmax = abs(mean)+2.6*sem),
                 linewidth = 1) + geom_point() +
  scale_y_discrete(limits = rev) + ylab("") + xlab(xlabel)  +
  scale_fill_manual(values=icr_conf_colors) +
  labs(fill="ICR Conf") + coord_cartesian(xlim = c(0, .15)) +
  theme_classic() + theme(legend.position = "none")
g2b
save_plot(paste0(output_dir_path,"/top_icrs_hemi_diff_placenta_51-100.jpg"), g2b,
          base_width = 3, base_height = 7)


# Make an order tiled image of all icrs
#_______________________________________________________________________________
# df_pdiff_all
# make r and c indices
df_tile <- select(df_pdiff_summary, icr_id, mean, icr_conf, zinc)
df_tile$row <- rep(1:50, ceiling(nrow(df_tile)/50))[1:nrow(df_tile)]
df_tile$column <- rep(1:ceiling(nrow(df_tile)/50), each=50)[1:nrow(df_tile)]

g2c <- ggplot(df_tile, aes(x=column, y = row, fill= icr_conf)) + 
  scale_fill_manual(values=icr_conf_colors) + labs(fill="ICR Conf") +
  geom_tile() +  theme_classic() + theme(legend.position = "none") +
  geom_tile(data = filter(df_tile, zinc == TRUE), aes(x=column, y = row), fill = "transparent", 
            color = "black", linewidth = .25) +
  ylab("ICR Rank") + xlab("Column Starting Rank") + 
  scale_x_continuous(breaks = seq(1,21,by=2), expand=c(0,0), position = "top",
                     labels =  seq(0, 1000, by = 100))  + 
  scale_y_continuous(expand=c(0,0), trans = "reverse", breaks = c(1,10,20,30,40,50))
g2c
save_plot(paste0(output_dir_path,"/top_icrs_hemi_diff_placenta.jpg"), g2c,
          base_width = 4, base_height = 2.5)







  
