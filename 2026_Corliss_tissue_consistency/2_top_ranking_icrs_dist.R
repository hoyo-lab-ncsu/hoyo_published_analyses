

# Figure 2A) Top top ranking ICRS                                   ############
#_______________________________________________________________________________
icr_mean_diff <- pairwise_dist_beta(df_beta = data_tri$icr_beta$icr_beta_df,
                            conf_group = df_icr_conf$icr_conf,
                            df_idat_sources = df_idat_sources)
# xlabel = "Mean Beta Diff."
xlabel = "Hemi-methylation Beta Dist.    "
icr_conf_colors = c("#fbb4ae", "#ccebc5", "#b3cde3")


# Matrix of distance between HUVEC vs CD14 betas for each patient
df_pdist_all = as.data.frame(icr_mean_diff$HUVEC_CD14)


# Summarize the range of distances across all patients for each ICR 
df_pdist_summary = data.frame( mean = rowMeans(df_pdist_all))

df_pdist_summary$p25 <- sapply(1:nrow(df_pdist_all),FUN = function(x) 
  quantile(as.matrix(df_pdist_all[x,]),0.25))
df_pdist_summary$p50 <- sapply(1:nrow(df_pdist_all),FUN = function(x) 
  quantile(as.matrix(df_pdist_all[x,]),0.50))
df_pdist_summary$p75 <- sapply(1:nrow(df_pdist_all),FUN = function(x) 
  quantile(as.matrix(df_pdist_all[x,]),0.75))
df_pdist_summary$sd <- sapply(1:nrow(df_pdist_all),FUN = function(x) 
  sd(as.matrix(df_pdist_all[x,])) )
df_pdist_summary$sem <- df_pdist_summary$sd/sqrt(nrow(df_pdist_all))
df_pdist_summary$icr_conf = factor(df_icr_conf$icr_conf, levels = c(1,2,3), ordered = TRUE)

# Add column that denotes zinc status of each ICR (close to Zinf finger site)
df_pdist_summary$zinc = as.numeric(str_replace(rownames(df_pdist_summary),"ICR_", "")) %in% 
  zinc_finger_icrs


# Sort entries by mean distance
sort_ind <- order(df_pdist_summary$mean, decreasing = FALSE)
df_pdist_summary <- df_pdist_summary[sort_ind,]
df_pdist_all <- df_pdist_all[sort_ind,]



# Add icr_id as a column from row names
df_pdist_summary$icr_id <- str_replace(rownames(df_pdist_summary), "_", " ")
df_pdist_summary$icr_id <- factor(x = df_pdist_summary$icr_id, 
                                  levels = df_pdist_summary$icr_id, ordered = TRUE)
df_pdist_summary$icr_name <- df_pdist_summary$icr_id
df_pdist_summary$icr_confidence = factor(df_pdist_summary$icr_conf, levels = c("1", "2", "3"), labels = c("High", "Medium", "Low"))
df_pdist_summary$rank <- 1:nrow(df_pdist_summary)

# Quick summaries for abstract submission
df_pdist_summary %>% group_by(icr_confidence) %>% summarize(mean_rank = mean(rank), mean_dist = mean(mean), 
                                                            sd_dist = sd(mean),
                                                            sem_dist = sd(mean)/sqrt(length(mean)))
nrow(df_pdist_summary %>% filter(mean < 0.191 & (icr_confidence != "High")))
nrow(df_pdist_summary %>% filter(mean < 0.191 & (icr_confidence == "Medium")))
nrow(df_pdist_summary %>% filter(mean < 0.191 & (icr_confidence == "Low")))


# Add closest gene to each ICR
df_pdist_summary <- merge(x = df_pdist_summary, y = imp_whole %>%
        select(icr_name, Nearest.Transcript, Distance.to.Nearest.Transcript),
      by = "icr_name", all.x = TRUE, all.y = FALSE, sort = FALSE)





# Data table exports
#_______________________________________________________________________________
# Making pretty version of table 
export_table_hemi_dist <- select(df_pdist_summary, c(icr_name, icr_confidence, zinc, Nearest.Transcript,
                                                     Distance.to.Nearest.Transcript, mean))
colnames(export_table_hemi_dist) <- c("ICR ID",	"ICR Confidence",	"Zinc Finger Proximity",
                                      "Nearest Transcript",	"Transcript Distance",	"Mean Hemi-Dist")



# ICR Rank 1-50
#_______________________________________________________________________________
# Plot
g2a <- ggplot(data = df_pdist_summary[1:50,], aes(y = icr_id, x = abs(mean)) ) + 
  geom_tile(aes(width = Inf, height = 1, fill = icr_conf)) + 
  geom_tile(data = filter(df_pdist_summary[1:50,], zinc == TRUE), 
            aes(y = icr_id, x = abs(mean), width = Inf, height = 1), 
            fill = "transparent", color = "black", linewidth = .5) +
  geom_linerange(aes(xmin = abs(mean)-2.6*sem, xmax = abs(mean)+2.6*sem),
                 linewidth = 1) + geom_point() +
  scale_y_discrete(limits = rev) + ylab("") + xlab(xlabel)  +
  scale_fill_manual(values = icr_conf_colors[sort(unique(df_pdist_summary[1:50,]$icr_conf))]) +
  labs(fill="ICR Conf") + coord_cartesian(xlim = c(0, .15)) +
  theme_classic() + theme(legend.position = "none")
g2a
save_plot(paste0(output_dir_path,"/ICR_hemi_rho_HUVEC-CD14_1-50.jpg"), g2a,
          base_width = 3, base_height = 7)



# ICR Rank 51-100
#_______________________________________________________________________________
# Plot
g2b <- ggplot(data = df_pdist_summary[51:100,], aes(y = icr_id, x = abs(mean)) ) + 
  geom_tile(aes(width = Inf, height = 1, fill = icr_conf)) + 
  geom_tile(data = filter(df_pdist_summary[51:100,], zinc == TRUE), 
            aes(y = icr_id, x = abs(mean), width = Inf, height = 1), 
            fill = "transparent", color = "black", linewidth = .5) +
  geom_linerange(aes(xmin = abs(mean)-2.6*sem, xmax = abs(mean)+2.6*sem),
                 linewidth = 1) + geom_point() +
  scale_y_discrete(limits = rev) + ylab("") + xlab(xlabel)  +
  scale_fill_manual(values=icr_conf_colors[sort(unique(df_pdist_summary[51:100,]$icr_conf))]) +
  labs(fill="ICR Conf") + coord_cartesian(xlim = c(0, .15)) +
  theme_classic() + theme(legend.position = "none")
g2b
save_plot(paste0(output_dir_path,"/ICR_hemi_rho_HUVEC-CD14_51-100.jpg"), g2b,
          base_width = 3, base_height = 7)


# Make an order tiled image of all icrs
#_______________________________________________________________________________
# df_pdist_all
# make r and c indices
df_tile <- select(df_pdist_summary, icr_id, mean, icr_conf, zinc)
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
save_plot(paste0(output_dir_path,"/ICR_hemi_rho_HUVEC-CD14_all.jpg"), g2c,
          base_width = 4, base_height = 2.5)




