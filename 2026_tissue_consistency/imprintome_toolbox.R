
library(broom)


#' Calculate agreement in beta values (either probe-wise or patient-wise). 
#' ICRs can be placed into subgroups (usually ICR confidence). Option to add an 
#' additional level to grouping_var to encompass all data without subgrouping.
#' 
#' @marg_dir marginal direction to measure correlations ("patient","probe"). Default: "patient".
#' @param grouping_var column name in imp_long for grouping beta measurements 
#' together before calculating correlations. Example includes:
#'   patient_wise: icr_cong, is_icr_zinc
#'   probe-wise: mat_aces_hi, mat_bmi_high
#' @param agreement type of correlation performed (default: pearson)
#' @param grouping_add_all for the grouping_var, option to add another group with
#'  all beta measurements included and add "All" as an additional factor level.
#' 
#' @return dataframe, where each row contains the correlation coefficient for a 
#' particular cpg site,. particular patient, and pair of tissue types.
marginal_grouped_agreement <- function(
    imp_long, marg_dir = "patient", grouping_var = "icr_conf", 
    agreement = "pearson", db_flag = TRUE, grouping_add_all = TRUE, 
    obs_handling = "complete.obs",  patientwise_cols = c("mat_bmi_hi", "mat_aces_hi"),
    probewise_cols = c("icr_conf", "is_icr_zinc")) 
  {
  
  if (db_flag) {save(list = ls(all.names = TRUE), file = "marginal_grouped_agreement.RData")}
  # load(file = "marginal_grouped_agreement.RData")
  
  if (marg_dir == "patient") {
    # Probe-wise
    if (any(grouping_var %in% patientwise_cols)) grouping_var = NA
    
    # Convert beta values to wide format for tissues to comptue correlations
    imp_wide <- imp_long %>% pivot_wider(names_from = tissue, values_from = beta)
    
    # Calculate sample-correlation: seprate by patient_id, then whatever column is specified by grouping_var
    imp_corr_each <- imp_wide %>% group_by(patient_id, !!sym(grouping_var)) %>%
      summarize(huvec_cd14 = cor(x = huvec, y = cd14, use = obs_handling, method = "pearson"), 
                huvec_plac = cor(x = huvec, y = plac, use = obs_handling, method = "pearson"), 
                cd14_plac =  cor(x = cd14, y = plac, use = obs_handling, method = "pearson"), n = n())
    
    if (grouping_add_all) {
      # add all group
      imp_corr_all <- imp_wide %>% group_by(patient_id) %>%
        summarize(huvec_cd14 = cor(x = huvec, y = cd14, use = obs_handling, method = "pearson"), 
                  huvec_plac = cor(x = huvec, y = plac, use = obs_handling, method = "pearson"), 
                  cd14_plac =  cor(x = cd14, y = plac, use = obs_handling, method = "pearson"), n = n()) %>% 
        mutate(!!grouping_var := factor("All"))
      imp_corr <- rbind(imp_corr_each,imp_corr_all)
    } else {imp_corr = imp_corr_each}
    
    # Convert correlation results to long forms
    imp_corr <- imp_corr %>% pivot_longer(cols = c(
      "huvec_cd14", "huvec_plac", "cd14_plac"), names_to = "tissue", values_to = "rho")

    # Recover all of the original extra columns in the input dataframe
    # Get list of repeat columns  (and probe-wise columns)
    rm_cols <- intersect(x = c("tissue", "beta", "icr_id", "cpg_id", probewise_cols), 
                         y = colnames(imp_long))
    imp_corr <- imp_corr %>% left_join(y = dplyr::select(imp_long %>% ungroup(), -!!rm_cols), 
                                      by = join_by(patient_id), multiple = "first", keep = FALSE)
    # Convert tissue to factor
    imp_corr$tissue <- factor(imp_corr$tissue)
    
    
  } else {
    # Probe-wise correlations                                              #####
    #__________________________________________________________________________#
    
    # If patient-wise variables specified for grouping_var, then set to NA. This
    # is for convenience for the analysis scripts, using one function call for both
    # probe-wise and patient-wise correlations with only changing marg_dir and no
    # other arguments
    
    
    # Detect whether imprintome site variable is icr_id or cpg_id
    imp_id_name = colnames(imp_long)[grepl(pattern = "(cpg_id)", colnames(imp_long))]
    if (length(imp_id_name)==0) imp_id_name = colnames(imp_long)[grepl(
      pattern = "(icr_id)", colnames(imp_long))]
  
    
    # probe
    imp_wide <- imp_long %>% pivot_wider(names_from = tissue, values_from = beta)
    # If grouping_var is NA, then make dummy column for the grouping variable to 
    # select for the correlations
    if (any(grouping_var %in% probewise_cols) || is.na (grouping_var))  { grouping_var = "dummy_col"
     imp_wide[[grouping_var]] = T
    }
    
    imp_corr_each <- imp_wide %>% group_by(!!sym(imp_id_name), !!sym(grouping_var)) %>%
      summarize(huvec_cd14 = cor(x = huvec, y = cd14, use = obs_handling, method = "pearson"), 
                huvec_plac = cor(x = huvec, y = plac, use = obs_handling, method = "pearson"), 
                cd14_plac =  cor(x = cd14, y = plac, use = obs_handling, method = "pearson"),
                n = n())
    
    if (any("dummy_col" %in% imp_corr_each)) imp_corr_each <- imp_corr_each %>% 
      dplyr::select(-!!sym(grouping_var))
    
    # Debugging Single test result
    # single_test<-imp_wide %>% subset( !!sym(imp_id_name)==imp_wide[[imp_id_name]][3])
    # cor(x = single_test$huvec, y = single_test$cd14, use = obs_handling, method = "pearson")
    # cor(x = single_test$huvec, y = single_test$plac, use = obs_handling, method = "pearson")
    # cor(x = single_test$cd14, y = single_test$plac, use = obs_handling, method = "pearson")
    # plot(single_test$huvec, single_test$cd14)
    # cg08627972
    

    
    # imp_corr_each = imp_corr_each
    # Convert correlation results to long forms
    imp_corr_each <- imp_corr_each %>% pivot_longer(cols = c(
      "huvec_cd14", "huvec_plac", "cd14_plac"), names_to = "tissue", values_to = "rho")
    
    # Recover all of the original extra columns in the input dataframe
    # Get list of repeat columns (and patient_wise columns)
    rm_cols <- intersect(x = c("tissue", "beta", patientwise_cols),  
                         y = colnames(imp_long))
    imp_corr_each <- imp_corr_each %>% left_join(y = dplyr::select(imp_long %>% ungroup(), -!!rm_cols), 
                                       by = join_by(!!sym(imp_id_name)), multiple = "first", keep = FALSE)
    
    
    if (grouping_add_all) {
      # add all group
      imp_corr_all <- imp_corr_each %>% mutate(icr_conf  = "All")
      imp_corr = rbind(imp_corr_each, imp_corr_all)
    } else {imp_corr = imp_corr_each}
    
    
    # Convert tissue to factor
    imp_corr$tissue <- factor(imp_corr$tissue)
    
  }
  
  return(imp_corr)
}




#' Perform 1s and 2s t.tests of beta values.
#' 
#' 1-sample analysis: data gets partitioned by both the partition_vars and
#'  the pairwise_var, and data is t.tested against 0.
#' 2-sample analysis: partition data by partition vars, and then perform 
#' pairwise t-test with the pairwise_var within each partition.
#' 
#' @param df dataframe assumed to have all columns specified.
#' @param comp_column
#' @param pairwise boolean, when true calculate pairwise analysis
#' @param partition_vars the column names to partition the variables into 
#' distinct subgroups that are not compared.
#' @param pairwise_var the column name for the grouping variable for 
#' the pairwise t.test.
#' @param p.adjust.1s_n manual override of the n used to adjust for multiple 
#' comparisons for 1s t.tests.
pairwise_stats <- function(df, comp_column = "group2", 
                           pairwise = T,
                           partition_vars = NA,
                           pairwise_var = "tissue",
                           p.adjust.1s_n = NA, db_flag = T) {
  if (db_flag) {save(list = ls(all.names = TRUE), file = "pairwise_stats.RData")}
  # load(file = "pairwise_stats.RData")
  
  
  # Calculate 1 tailed p-value against zero for each group
  df_1s <- df %>% 
    group_by(!!!syms(c(partition_vars, pairwise_var))) %>% 
    summarize(pval = t.test(x=rho, y=NULL)$p.value) %>% 
    rename(group = !!pairwise_var)  %>%
    mutate(grouping_var = !!pairwise_var)
    
  df_1s$pval <-custom_p.adjust(df_1s$pval, method = "holm", 
                                   n=length(unique(df[[pairwise_var]])))
  df_1s$pstr <- pval_str(df_1s$pval)
  
  
  if (pairwise) {
  
    # Calculate pairwise p-values within each partition independently
    df_pairs <-  df %>% group_by(!!!syms(partition_vars)) %>%
      summarize(out = list(stats::pairwise.t.test(x = rho, g = !!!syms(pairwise_var), 
                                           p.adjust.method = "holm") %>% tidy)) %>% unnest(c(out))
    df_pairs$group1 <- factor(df_pairs$group1)
    df_pairs$group2 <- factor(df_pairs$group2)
    df_pairs$n <- length( levels(df[[pairwise_var]]) )
    
    # Make a copy of df_pairs and reverse group order, merge into one dataframe
    df_pairs2 <- df_pairs
    df_pairs2 <- df_pairs2 %>% 
      rename(temp2 = group2, temp1 = group1) %>% 
      rename(group1=temp2, group2=temp1)
    df_mirrored <- rbind(df_pairs, df_pairs2) 
    df_mirrored$p_val_sig = df_mirrored$p.value < 0.05
    
    # Assign a letter for each group for comparisons
    #   based on order of group2 (defined by level order in df_1s)
    df_mirrored <- df_mirrored %>%  group_by(!!!syms(partition_vars)) %>%
      mutate(sig_letter = 
               sapply(!!!syms(comp_column), function(x) 
                 letters[ which(as.character(x) == as.character(df_1s$group))[1] ]  )  )  #[1:n[1]
    # Remove any letter with insignificant p-values
    df_mirrored$sig_letter[df_mirrored$p.value >= 0.05] = ""
    df_mirrored$sig_letter[is.na(df_mirrored$p.value)]  = ""
    
    # Within each partition,
    df_2s <- df_mirrored %>%  group_by(!!!syms(partition_vars), group1)  %>% 
      summarize(sig_strs = paste0(sig_letter, collapse = "")) %>%
      rename(group = group1) %>%
      mutate(pairwise_var = pairwise_var)
  
  } else {
    df_2s = NULL
  }
  

  # Include 1s and 2s results  
  out <- rbind(df_1s %>% mutate(p_val_type = "1s"),
               df_2s %>% mutate(p_val_type = "2s"))
  
  # Rename "group" var to variable name, this is needed downstream for plotting 
  # faceting for ggplot when using multiple input data frame, facet variable 
  # must be the same.
  out <- out %>% rename(!!sym(pairwise_var) := group) 
  
  return(out)
}





#'
plot_agreement_by_confidence <- function(df_rho, sum_df, df_stat = NULL, x_var = NA,
                                         fill_var = NA, facet_var = NA,
                                         subset_level = NA, 
                                         xlab_str = "ICR Confidence", 
                                         ylab_str = "Rho of ICR Betas",
                                         output_dir_path,
                                         include_neg_y_axis = NA,
                                         export_filename = "ICR_Rho_HUVEC-CD14.jpg",
                                         fig_dim = c(4,2), db_flag = F, pretty_x = F,
                                         ylim = c(0, 1.03)) {
  if (db_flag) {save(list = ls(all.names = TRUE), file = "plot_agreement_by_confidence.RData")}
  # load(file = "plot_agreement_by_confidence.RData")
  
  if (is.na(include_neg_y_axis)) {include_neg_y_axis = (quantile(df_rho$rho, 0.05)<0)}
  
  
  # IF there is a fill variable, we have to dodge positions on x axis, otherwise don't
  if (!is.na(fill_var)){
    ggpos = position_dodge(width=0.9)
  } else { ggpos = position_identity()}
  
  # If design score is used, only keep some of the pairwise stat entries
  if (!is.na(fill_var) && fill_var == "design_score") {
    df_stat <- df_stat %>% filter(design_score == F)
    df_stat$sig_strs=""
  }
  
  # Subset the dataframe to a specific tissue if desired
  if (!is.na(subset_level)) {
    df_rho <-  df_rho %>%  subset(tissue == subset_level)
    sum_df <-  sum_df %>%  subset(tissue == subset_level)
    df_stat <- df_stat %>% subset(tissue == subset_level)
  }
  
  
  # Pretty print x labels with function
  if (pretty_x) {
    pretty_fun = function(df, x) toupper(str_replace(string = df[[x]],pattern = "_", replacement = " \n"))
    df_rho[[x_var]] <- pretty_fun(df_rho, x_var)
    sum_df[[x_var]] <- pretty_fun(sum_df, x_var)
    if (!is.null(df_stat))  df_stat[[x_var]] <- pretty_fun(df_stat, x_var)
  }
  

  # Fill statically if there is no fill var, otherwise place within AES with colors specified
  if (!is.na(fill_var)) {
    g0 <- ggplot(data = df_rho, aes(x = .data[[x_var]], y = rho, fill = .data[[fill_var]])) +
      geom_violin(color = NA, position = ggpos) +
      scale_fill_manual(values = c("#41D3D7", "#FA9892"))
  } else {
    g0 <- ggplot(data = df_rho, aes(x = .data[[x_var]], y = rho)) + 
      geom_violin(fill = "#E5E5E5", color = NA, position = ggpos)
  }
  

  # Plot of agreement partitioned by icr confidence
  g0 <- g0 +
    geom_errorbar(data = sum_df, width = 0.2, position = ggpos, aes(
      x = .data[[x_var]], y = mean, ymin = mean-sd, ymax = mean+sd)) +
    geom_point(data = sum_df, shape = "\U2014", size = 3, position = ggpos, aes(x = !!sym(x_var), y = mean))
  
  if (!is.null(df_stat)) {
    g0 <- g0 + geom_text(data = df_stat %>% subset( (p_val_type == "1s")),
                         aes(x = .data[[x_var]], label = pstr, y = 1-0.1-0.1*include_neg_y_axis-0.1*!is.na(facet_var)),
                         color = "black") +
      geom_text(data = df_stat %>% subset(p_val_type == "2s"),
                aes(x = .data[[x_var]], label = sig_strs, y = 1.0),
                color = "black", size = 3.5)
  }
  
  # faceting if requested
  if (!is.na(facet_var)) g0 <- g0 + facet_wrap(~get(facet_var))
  
  # Add axis labels, theme, and extra space for significance labels
  if (include_neg_y_axis) ylim = c(-1,1.03)
   g0 <- g0 +
    xlab(xlab_str) + ylab(ylab_str) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = ylim) 
   if (include_neg_y_axis)  g0 <- g0 + geom_hline(yintercept = 0, color = "grey20")
   
  print(g0)
  
  save_plot(paste0(output_dir_path,"/",export_filename), g0,
            base_width = fig_dim[1], base_height = fig_dim[2])
  return(g0)
}


beta_mat_to_long <- function(cpg_beta, filt_study_data, db_flag = FALSE) {
  if (db_flag) {save(list = ls(all.names = TRUE), file = "beta_mat_to_long.RData")}
  # load(file = "beta_mat_to_long.RData")
  
  temp_list = list()
  cpg_meta_data <- add_metadata_to_imp_sites(rownames(cpg_beta),imp_type = "cpg") %>% 
    dplyr::select(icr_id, icr_conf, is_icr_zinc)
  cpg_meta_data$icr_conf = factor(cpg_meta_data$icr_conf, labels = c("High", "Medium", "Low"))
  cpg_meta_data$is_icr_zinc  = factor(cpg_meta_data$is_icr_zinc, labels = c("-Zinc Finger", "+Zinc Finger"))
  
  # unq_icr_conf <- cpg_meta_data %>% distinct()
  # df_icr_conf
  # test <- left_join(x = df_icr_conf, y = unq_icr_conf, by = join_by("icr_id"),keep = FALSE,multiple = "first")
  
  
  # For each patient, create a table of ICRs that have beta values from all 3 tissues
  for (n in 1:nrow(filt_study_data)) {
    temp_list[[n]] <- 
      cpg_meta_data %>% cbind(
        data.frame(patient_id = filt_study_data$ship_ID[n], 
                   mat_bmi_hi = filt_study_data$mat_bmi_hi[n],
                   mat_aces_hi = filt_study_data$mat_aces_hi[n],
                   cpg_id = factor(rownames(cpg_beta)),
                   huvec = cpg_beta[ , filt_study_data$HUVEC_Patient_ID[n]],
                   cd14 = cpg_beta[ ,  filt_study_data$CD14_Patient_ID[n]],
                   plac = cpg_beta[ ,  filt_study_data$Placental_Patient_ID[n]]
        )  )
    
  }
  cpg_wide <- do.call(rbind, temp_list)
  
  # Convert to long format
  cpg_long <- cpg_wide %>% pivot_longer(cols = c("huvec", "cd14", "plac"), names_to = "tissue", values_to = "beta")
  cpg_long$tissue = factor(cpg_long$tissue)
  cpg_long$mat_bmi_hi  <- factor(cpg_long$mat_bmi_hi)
  cpg_long$mat_aces_hi <-factor(cpg_long$mat_aces_hi)
  
  
 return(cpg_long)
  
}





plot_time_consistency_dual <- function(df, meth = "cpg", plot_dim = c(1.9, 2)) {
  
  g2 <-  ggplot(df, aes(x = diff12)) + 
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.05) +
    geom_vline(xintercept = quantile(df$diff12, probs = c(0.05, 0.95)), color = "black") +
    xlab(paste0("Beta diff")) + ylab("Rel. Frequency") +
    coord_cartesian(expand = FALSE) +
    theme_classic() 
  save_plot(filename = paste0("fig/", meth, "_beta_diff_12.jpg"),plot = g2, base_height = plot_dim[1], base_width = plot_dim[2])  
  print(sprintf("Q95 T1 v T2: %.3f, %.3f", quantile(df$diff12, probs = 0.05), quantile(df$diff12, probs = 0.95)))
  
  
  g3 <-  ggplot(df, aes(x = diff12, y= -log(adj.p.value12),color = adj.p.value12<0.05) ) + 
    geom_point(size = 0.5) +
    geom_vline(xintercept = quantile(df$diff12, probs = c(0.05, 0.95)), color = "black") +
    xlab("Beta diff") + ylab("-log(Beta p-value)") + theme_classic() + theme(legend.position = "none")
  g3
  save_plot(filename = paste0("fig/", meth, "_volcano_12.jpg"),plot = g3, base_height = plot_dim[1], base_width = plot_dim[2])   
  
  
  
  
  g2 <-  ggplot(df, aes(x = diff23)) + 
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.05) +
    geom_vline(xintercept = quantile(df$diff23, probs = c(0.05, 0.95)), color = "black") +
    xlab("Beta diff") + ylab("Rel. Frequency") +
    coord_cartesian(expand = FALSE) +
    theme_classic() 
  save_plot(filename = paste0("fig/", meth, "_beta_diff_23.jpg"),plot = g2, base_height = plot_dim[1], base_width = plot_dim[2])  
  print(sprintf("Q95 T1 v T2: %.3f, %.3f", quantile(df$diff23, probs = 0.05), quantile(df$diff23, probs = 0.95)))
  
  g3 <-  ggplot(df, aes(x = diff23, y= -log(adj.p.value23),color = adj.p.value23<0.05) ) + 
    geom_point(size = 0.5) +
    geom_vline(xintercept = quantile(df$diff23, probs = c(0.05, 0.95)), color = "black") +
    xlab("Beta diff") + ylab("-log(Beta p-value)") + theme_classic() + theme(legend.position = "none")
  g3
  save_plot(filename = paste0("fig/", meth, "_volcano_23.jpg"),plot = g3, base_height = plot_dim[1], base_width = plot_dim[2])    
 
  
  
  g2 <-  ggplot(df, aes(x = diff13)) + 
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.05) +
    geom_vline(xintercept = quantile(df$diff13, probs = c(0.05, 0.95)), color = "black") +
    xlab("Beta diff") + ylab("Rel. Frequency") +
    coord_cartesian(expand = FALSE) +
    theme_classic() 
  save_plot(filename = paste0("fig/", meth, "_beta_diff_13.jpg"),plot = g2, base_height = plot_dim[1], base_width = plot_dim[2])  
  print(sprintf("Q95 T1 v T3: %.3f, %.3f", quantile(df$diff13, probs = 0.05), quantile(df$diff13, probs = 0.95)))
  
  g3 <-  ggplot(df, aes(x = diff13, y= -log(adj.p.value13),color = adj.p.value13<0.05) ) + 
    geom_point(size = 0.5) +
    geom_vline(xintercept = quantile(df$diff13, probs = c(0.05, 0.95)), color = "black") +
    xlab("Beta diff") + ylab("-log(Beta p-value)") + theme_classic() + theme(legend.position = "none")
  g3
  save_plot(filename = paste0("fig/", meth, "_volcano_13.jpg"),plot = g3, base_height = plot_dim[1], base_width = plot_dim[2])  
}



get_value_in_range <- function(values, ages, age_target, age_range, Chip) {
  min.ind <- which.min(abs(age_target - ages))
  if ((ages[min.ind] - age_target) > age_range) {out <- list (value = NA, age_months = 0, chip = NA)
  } else { out <- list (value = values[min.ind], age_months = ages[min.ind],
                        chip = ifelse(length(Chip)==0, NA, Chip[min.ind]) )}
  return(out)
}
