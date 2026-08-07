
library(broom)

#' pcor_beta 
#' @description Define function to calculate binary and continuous agreement (can be used at ICR
#' level or CPG)
#' 1. Continuous: pearson rho
#' 2. Continuous: R2
#' Binary: Positive agreement.
pcor_beta <- function(df_beta, conf_group, df_idat_sources, agreement = "pearson",
                      conf_group_levels = c("High", "Medium", "Low", "All"),
                      db_flag = FALSE) {
  if (db_flag) {save(list = ls(all.names = TRUE), file = "pcor_beta.RData")}
  # load(file = "pcor_beta.RData")
  
  # Table: ship_id, group, tissue 1, tissue2, rho
  df_template = data.frame(ship_id = 0, icr_conf = 0, tissue = "none", n = 0, rho = NA)
  
  # Pairwise correlation for each icr group
  pcor_icr_groups = function(df_beta, df_idat_source, df_template, agreement) {
    # Loop for each icr confidence group
    for (g in 1:4) {
      
      # Group 4 is all icrs, group 1-3 matches confidence level
      if (g != 4) {conf_bv = conf_group == g  
      } else { conf_bv = rep(TRUE, length(conf_group)) }
      
      # Initialize data.frame
      df_pcor_conts = df_template
      
      # Agreement is calculated base on variable type
      if(agreement == "positive") {
        # Binary agreement measured with positive agreement
        sub_mat <- df_beta[conf_bv, as.character(df_idat_source[2:4])]
        cormat <- matrix(rep(0,9), ncol = 3, dimnames =
                           list(as.character(df_idat_source[2:4]),
                                as.character(df_idat_source[2:4])))
        cormat[2,1] = sum(sub_mat[,1] == sub_mat[,2])/nrow(sub_mat)
        cormat[3,1] = sum(sub_mat[,1] == sub_mat[,3])/nrow(sub_mat)
        cormat[3,2] = sum(sub_mat[,2] == sub_mat[,3])/nrow(sub_mat)
        diag(cormat) = 1
        cormat = Matrix::forceSymmetric(cormat,uplo="L")
        
      } else if (agreement == "pearson") {
        # Continuous Agreement with Pearson (cols: huvec, cd14, placenta)
        cormat <- cor(df_beta[conf_bv, as.character(df_idat_source[2:4])],
                      use = "complete.obs", method = "pearson")
        
     
        # x$estimate
        # x$p.value
      } else if (agreement == "rsquared") {
        sub_mat <- df_beta[conf_bv, as.character(df_idat_source[2:4])]
        cormat <- matrix(rep(0,9), ncol = 3, dimnames =
                           list(as.character(df_idat_source[2:4]),
                                as.character(df_idat_source[2:4])))
        
        cormat[2,1] = summary(lm(sub_mat[,1]~sub_mat[,2]))$r.squared
        cormat[3,1] = summary(lm(sub_mat[,1]~sub_mat[,3]))$r.squared
        cormat[3,2] = summary(lm(sub_mat[,2]~sub_mat[,3]))$r.squared
        diag(cormat) = 1
        cormat = Matrix::forceSymmetric(cormat,uplo="L")
      }
      
      
      
      # Fill in agreement for each of the 3 pair-wise combinations
      # Rho: some measure of correlation, difference, or R^2
      df_pcor_conts[1,] <- data.frame(ship_id = df_idat_source$ship_ID, icr_conf = g,
                                      tissue = "HUVEC Vs.\nCD14", n = sum(conf_bv), 
                                      rho = cormat[2,1])
      df_pcor_conts[2,] <- data.frame(ship_id = df_idat_source$ship_ID, icr_conf = g,
                                      tissue = "HUVEC Vs.\nPlacenta", n = sum(conf_bv), 
                                      rho = cormat[3,1])
      df_pcor_conts[3,] <- data.frame(ship_id = df_idat_source$ship_ID, icr_conf = g,
                                      tissue = "CD14 Vs.\nPlacenta", n = sum(conf_bv), 
                                      rho = cormat[3,2])
      list_df_pcor_cont[[g]] <- df_pcor_conts
    }
    # Bind data across icr confidence groups
    df_pcor_group <- do.call(rbind, list_df_pcor_cont) 
    df_pcor_group$tissue <- factor(df_pcor_group$tissue)
    return(df_pcor_group)
  }
  # Initialize list
  list_ship_id_df_pcor_cont <- list()
  
  # Loop through each of the patients
  for (n in 1:nrow(df_idat_sources)) {
    list_df_pcor_cont = list()
    list_ship_id_df_pcor_cont[[n]] = pcor_icr_groups(df_beta, df_idat_source = df_idat_sources[n,], 
                                                     df_template,agreement)
  }
  # Bind data across patient
  df_pcor_cont <- do.call(rbind, list_ship_id_df_pcor_cont)
  df_pcor_cont$icr_conf <- factor(df_pcor_cont$icr_conf)
  levels(df_pcor_cont$icr_conf) <- conf_group_levels
  
  # Join extra metadata found in df_idat_source
  df_pcor_cont <- left_join(df_pcor_cont, select(
    df_idat_sources, -c("HUVEC_Patient_ID", "CD14_Patient_ID", "Placental_Patient_ID")), 
    by = c("ship_id" = "ship_ID"))
  
  return(df_pcor_cont)
}


#' pair_dist_beta 
#' @description calculate average distance between betas for each icr/cpg 
#' 
#' @param df_beta cpg/probes (rows) by patients (columns)
#' @param conf_group confidence groups for each cpg/icr
#' @param df_idat_sources list of patient_ids for each tissue type and patient 
#' @param dist_func a function used to compute distance between beta values from
#'  different tissue from same person. Default is euclidean distance from methylated.
#'  Possible options:
#'  - hemi_methylation distance: \(x,y) sqrt((0.5-x)^2 + (0.5-y)^2)
#'  - mean absolute hemi_methylation difference: \(x,y) (abs(0.5-x) + abs(0.5-y))/2
#'  - absolute diff: \(x,y) abs(x-y)
#' @return
#' 
#'  
pairwise_dist_beta <-
  function(df_beta, conf_group, df_idat_sources, tissues = c("HUVEC", "CD14", "Placenta"), 
           dist_func = \(x,y) sqrt((0.5-x)^2 + (0.5-y)^2), db_flag = FALSE) {
    if (db_flag) {save(list = ls(all.names = TRUE), file = "pairwise_dist_beta.RData")}
    # load(file = "pairwise_dist_beta.RData")
  
  # Table: ship_id, group, tissue 1, tissue2, rho
  df_template = data.frame(cpg_id = 0, icr_conf = 0, tissue1 = "", 
                           tissue2 = "", diff = 0)
  
  # Matrix cpg_ids and mean diff across all patients, for each tissue
  df_dist = matrix(rep(0, nrow(df_beta) * 3),  ncol = 3)
  colnames(df_dist) <- c(paste0(tissues[1], '_', paste0(tissues[2])),
                         paste0(tissues[1], '_', paste0(tissues[3])),
                         paste0(tissues[2], '_', paste0(tissues[3])))
  rownames(df_dist) <- rownames(df_beta)
  
  # CpG_id/icr_id x patient_id, val = diff
  df_diffs_1 <- df_diffs_2 <- df_diffs_3 <- 
    matrix(rep(0, nrow(df_beta) * nrow(df_idat_sources)), ncol =
             nrow(df_idat_sources), dimnames =  list(rownames(df_beta), 
                                                     df_idat_sources$ship_ID))
  for (n in 1:nrow(df_idat_sources)) {
    df_diffs_1[,n] = dist_func(df_beta[,df_idat_sources[n,2]], 
                           df_beta[,df_idat_sources[n,3]])
    df_diffs_2[,n] = dist_func(df_beta[,df_idat_sources[n,2]], 
                           df_beta[,df_idat_sources[n,4]])
    df_diffs_3[,n] = dist_func(df_beta[,df_idat_sources[n,3]], 
                           df_beta[,df_idat_sources[n,4]])
  }
  df_dist[,1] = rowMeans(abs(df_diffs_1))
  df_dist[,2] = rowMeans(abs(df_diffs_2))
  df_dist[,3] = rowMeans(abs(df_diffs_3))
  
  out <- list(df_dist, df_diffs_1, df_diffs_2, df_diffs_3)
  names(out) <- c("mean_abs_diffs", colnames(df_dist))
  
  return(out)
}


pval_str = function(x) {
  p1 <- ifelse(x <0.05, "*","" )
  p2 <- ifelse(x <0.01, "*","" )
  p3 <- ifelse(x <0.001, "*","" )
  pstr = paste0(p1,p2,p3)
  return(pstr)
}



pairwise_stats_old <- function(df, comp_column="group2", comp_levels = c("High", "Medium", "Low", "All"), 
                           pairwise = TRUE, group_by_var = "tissue", pairwise_var = "icr_conf", db_flag = TRUE) {
  if (db_flag) {save(list = ls(all.names = TRUE), file = "pcorr_stats.RData")}
  # load(file = "pcorr_stats.RData")
  
  
  # Calculate 1 tailed p-value
  df_pval_0 <- df %>% group_by(tissue, icr_conf) %>% 
    summarize(pval = t.test(x=rho, y=NULL)$p.value)
  
  df_pval_0$pval <-custom_p.adjust(df_pval_0$pval, method = "holm", 
                                   n=length(unique(df$tissue)))
  
  df_0 <- df_pval_0 %>% mutate(group1 = icr_conf) %>% mutate(group2=factor(NA)) %>% 
    relocate(group2, .before = pval)
  df_0$pstr <- pval_str(df_0$pval)
  
  
  if (pairwise) {
    
    # Calculate pairwise p-value
    df_pairs <-  df %>% group_by(!!!syms(group_by_var)) %>%
      summarize(out = list(pairwise.t.test(x = rho, g = !!!syms(pairwise_var), 
                                           p.adjust.method = "holm") %>% 
                             tidy)) %>%  unnest(c(out))
    
    df_pairs2 <- df_pairs
    df_pairs2 <- df_pairs2 %>% 
      rename(temp2 = group2, temp1 = group1) %>% 
      rename(group1=temp2, group2=temp1)
    
    df_mirrored <- rbind(df_pairs, df_pairs2)
    
    
    # Assign letters characters of pairwise p-value comparisons based on th groups specified by comp_levels
    df_mirrored$sig_letter = unname(sapply(df_mirrored[[comp_column]], function(x) 
      letters[1:length(comp_levels)][which(x ==comp_levels)]))
    df_mirrored$sig_letter[df_mirrored$p.value >= 0.05] = ""
    
    
    df_compare <- df_mirrored %>% group_by(!!!syms(group_by_var), group1) %>% 
      summarize(sig_strs = paste0(sig_letter, collapse = ""))
  } else {
    df_compare = NULL
  }
  
  # out <- list(pvalue_0 = df_0, pvalue_comp = df_compare)
  
  out <- rbind(df_0 %>% mutate(p_val_type = "pvalue_0"),
               df_compare %>% mutate(p_val_type = "pvalue_comp"))
  
  return(out)
}



