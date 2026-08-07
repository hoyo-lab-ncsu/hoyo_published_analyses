


# User: First set working directory to base of tdhia package!

# Load libraries                                                ################
#_______________________________________________________________________________

library(tidyverse)
library(devtools)
library(foreach)
library(doParallel)
library(gplots)
library(cowplot)
library(pheatmap)
library(GGally)
library(haven)
library(preprocessCore)
library(ggVennDiagram)
devtools::load_all(paste0(dirname(dirname(getwd())), "/tdhia"))


# Change these file paths for folder with IDAT files and output folder
# Define local variables
github_path <- dirname(dirname(getwd()))
OVERWRITE_TEMP_DATA = F
# idat_dir_paths = paste0(github_path, "/data/IDAT_files 15Nov23")
# 
# C:\Users\bruce\Documents\GitHub\data\
dir.create("fig", showWarnings = F)

idat_dir_paths <- paste0(github_path, "/data/idats_2024_ship-next_imprintome") #list.dirs(paste0(github_path, "/data/idats_2024_ship-next_imprintome/IDAT"))[-1]

output_dir_path = paste0(github_path, "/data/study_tissue_consistency")
study_dir_path <- paste0(github_path, "/data/study_tissue_consistency/study_data")
ref_path = paste0(github_path, "/tdhia_scripts/study_tissue_consistency")
source(paste0(ref_path, "/imprintome_toolbox.R"))

# Perform association study with
# 1) HUVECs
# 2) CD13
# And compare hits, direction of hits



# load study data
# df_study <- haven::read_sas(paste0(output_dir_path, "/ship_bp_hr_06dec23.sas7bdat"))

study_data <- haven::read_sas(paste0(output_dir_path, "/study_data/ship_trudiag_ids_followup.sas7bdat"))
# write.csv(study_data, file = "imprintome_time_tissue_consistency.csv")

study_data$tp_num<- study_data %>%
  group_by(ship_id) %>% # Group by the columns that define a duplicate record
  reframe(tp_num= 1:n()) %>% pull(tp_num)

study_data$total_tp<- study_data %>%
  group_by(ship_id) %>% # Group by the columns that define a duplicate record
  reframe(total_tp= rep(n(),n())) %>% pull(total_tp)



study_data_reps <- study_data %>%
  group_by(ship_id) %>% # Group by the columns that define a duplicate record
  summarize(total_tp= n())

# Remove records with only 2 replicates
study_data <- study_data %>% filter(total_tp> 1) 


# Study Metadata Pre-processing
#-------------------------------------------------------------------------------

# add unique ID for each row
study_data$id <- paste0(study_data$ship_id, "_", study_data$tp_num)

# Remove patients without associated IDAT files
filt_study_metadata <- study_data %>%
  filter(study_data$Chip != "" & !is.na(child_age_months))
rownames(filt_study_metadata) <- filt_study_metadata$id



# Report how many patients removed because of missing idat data
cat(sprintf("%.0f%% of patients discarded b/c no IDAT file (%.0f/ %.0f).\n",
            100 * (nrow(study_data) - nrow(filt_study_metadata))/nrow(study_data),
            (nrow(study_data) - nrow(filt_study_metadata)), nrow(study_data)))



attributes(filt_study_metadata$Chip) <- NULL



# IDAT processing
#-------------------------------------------------------------------------------

# 1) Load IDATS and convert to probe beta matrix
# Output: rows: probe_id     x     columns: patient

out_path <- paste0(output_dir_path, "/study_temp_data/time_consistency_betas.rda")
dir.create(dirname(out_path), recursive = T, showWarnings = F)
if (!file.exists(out_path) || OVERWRITE_TEMP_DATA) {
  probe_beta <-
    load_idata_to_probes(idat_dir_paths = idat_dir_paths, multicore = TRUE,
                         idat_basenames = filt_study_metadata$Chip,
                         quantile_norm = FALSE, db_flag = T, merge_replicates = "pre_beta")
  save(probe_beta, file = out_path)
} else {load(out_path)}


# 2) Filter probes that are not mapped and discard poor signal
# Output: rows: probe_id     x     columns: patient
filt_probe_beta <-
  filter_probes(probe_beta = probe_beta, discard_unmapped_probes = TRUE ,
                max_sig_pval = 0.2, set_failed_betas_na = TRUE, 
                max_probe_fail_rate = 0.2, min_design_score = NA,
                discard_failed_probes = TRUE, db_flag = TRUE)


#3  Convert probe beta matrix to a cpg beta matrix
# Output: rows: cpg_id     x     columns: patient
cpg_beta <- convert_probes_to_cpgs(filt_probe_beta, quantile_norm = FALSE, 
                                   db_flag = FALSE, smooth_adj_cpgs = FALSE)


#4  Convert probe beta matrix to an icr beta matrix
# Output: rows: icr_id     x     columns: patient
icr_beta <- convert_cpgs_to_icrs(cpg_beta, max_icr_fail_rate = 0.2)


# Read in ICR ranks of mean beta difference, CD14 versus HUVEC
df_rank_diff <- read.csv(file = paste0(output_dir_path, "/ICR_ranks_diff_CD14_HUVEC.csv"))
# Get list of top icrs
top100_icrs <- df_rank_diff$ICR.ID[1:100] %>% str_replace(pattern = "\\s", "_")
top100_cpgs  <-  tdhia::manifest_v1A2_design_scores  %>% filter(icr_id %in% top100_icrs) %>% pull("cpg_id")
# not100_cpgs <-   tdhia::manifest_v1A2_design_scores %>% filter(!(icr_id %in% top100_icrs)) %>% pull("cpg_id")





# Create sample (rows) x cpg site (col) matrix
df_cpgs  <- as.data.frame(t(cpg_beta$cpg_beta_df)) %>% rownames_to_column("Chip")
# Add timepoint and ship_id data
metada_cols <-  c("Chip", "ship_id", "child_age_months", "tp_num", "total_tp")
df_cpgs2 <- df_cpgs %>% 
  left_join(filt_study_metadata %>% select(all_of(metada_cols)), 
           by = join_by("Chip"), keep = F ) 

# Subsets of cpg ids of top100 and those outside of top100
df_cpg_long <- df_cpgs2 %>% #select(c(any_of(metada_cols), any_of(top100_cpgs))) %>%
  pivot_longer(cols = starts_with("cg"), names_to = "cpg_id", values_to = "beta")
df_icr_long <- df_cpg_long %>% left_join(tdhia::manifest_v1A2_design_scores %>% select(cpg_id, icr_id), 
                          by = join_by(cpg_id), keep = F,multiple = "first") %>%
  group_by(ship_id, icr_id, child_age_months) %>% summarize( beta = mean(beta,na.rm = T))


# 
# 
# df_cpg_not100 <- df_cpgs2 %>% select(c(any_of(metada_cols), any_of(not100_cpgs)))


my_cor.test <- \(x, y) {
  out <- try(cor.test(x, y), silent = T)
  if (is(out, "try-error")) { out = list(estimate = NA, p.value = NA)}
return(out)
}

my_regress <- \(x, y) {
  out <- try({
    d <- data.frame(x=x, y=y)
    model <- lm(y ~ x, data = d)
    
    # Summary of the model

    estimate = model$coefficients[2]
    p.value = summary(model)$coefficients[2,4]
    out <- list(estimate = estimate, p.value = p.value) }, 
    silent = T)
  if (is(out, "try-error")) { out = list(estimate = NA, p.value = NA)}
  return(out)
}




# CPG Analysis                                                        ##########
#_______________________________________________________________________________
tps = c(6, 24, 48)


# Calculate beta slope for each patient, for each cpg_id, across time.
# Only include patients with 3+ timepoints
df_cpg_slopes <- df_cpg_long %>% 
  filter( Chip   %in%   (filt_study_metadata %>% filter(total_tp > 2) %>% pull(Chip)) ) %>% 
  filter(cpg_id %in% top100_cpgs) %>%
  group_by(ship_id, cpg_id) %>%
  summarize(mean_age = mean(child_age_months, na.rm = TRUE),
            # cor_coef = my_cor.test(x = child_age_months/12, y = beta)$estimate,
            # cor_p_value = my_cor.test(x = child_age_months/12, y = beta)$p.value,
            # regr_coef = my_regress(x = child_age_months/12, y = beta)$estimate,
            # regr_p_value = my_regress(x = child_age_months/12, y = beta)$p.value,
            # Actual list of datapoints included for debugging
            x = paste0(paste(child_age_months, collapse = ",")),
            y =  paste0(paste(beta, collapse = ",")),
            beta1 = get_value_in_range(beta, child_age_months, age_target = tps[1], age_range = 6, Chip)$value,
            time1 = get_value_in_range(beta, child_age_months, age_target = tps[1], age_range = 6, Chip)$age_months,
            chip1 = get_value_in_range(beta, child_age_months, age_target = tps[1], age_range = 6, Chip)$chip,
            
            beta2 = get_value_in_range(beta, child_age_months, age_target = tps[2], age_range = 6, Chip)$value,
            time2 = get_value_in_range(beta, child_age_months, age_target = tps[2], age_range = 6, Chip)$age_months,
            chip2 = get_value_in_range(beta, child_age_months, age_target = tps[2], age_range = 6, Chip)$chip,
            
            beta3 = get_value_in_range(beta, child_age_months, age_target = tps[3], age_range = 6, Chip)$value,
            time3 = get_value_in_range(beta, child_age_months, age_target = tps[3], age_range = 6, Chip)$age_months,
            chip3 = get_value_in_range(beta, child_age_months, age_target = tps[3], age_range = 6, Chip)$chip,
            
            n = n())
# 243 distinct cpgs
# 48 patients

# df_cpg_slopes_sum <- df_cpg_slopes %>% group_by(ship_id) %>% summarize(diffs = age_diffs[1])

# df_cpg_slope_top100 <- df_cpg_slopes %>% filter(cpg_id %in% top100_cpgs)
# df_cpg_slope_top100$adj_cor_p_value <- tdhia::custom_p.adjust(
#   df_cpg_slope_top100$cor_p_value, "fdr", n = 100 * length(unique(df_cpg_slopes$ship_id)))
# df_cpg_slope_top100$adj_regr_p_value <-  tdhia::custom_p.adjust(
#   df_cpg_slope_top100$regr_p_value, "fdr", n = 100 * length(unique(df_cpg_slopes$ship_id)))
# df_cpg_slope_top100 <- df_cpg_slope_top100 %>%  left_join(y = tdhia::manifest_v1A2_design_scores %>% 
#                                                             select(cpg_id, icr_id), by = join_by(cpg_id), keep = F, multiple = "first")
# df_cpg_slope_top100$icr_num = as.numeric(stringr::str_replace(df_cpg_slope_top100$icr_id, "^ICR_", ""))



# 2 Sample test at cpg level
# df_cpg_2s_top100 <- df_cpg_slopes 
df_cpg_2s_test <- df_cpg_slopes %>% filter(cpg_id %in% top100_cpgs) %>% 
  filter(!is.na(beta1) & !is.na(beta2)) %>% group_by(cpg_id) %>% 
  summarize(
  p.value12 = t.test(beta1, beta2, paired = T)$p.value, 
  diff12 = diff(t.test(beta1, beta2)$estimate), 
  n12 = sum(!is.na(beta1) & !is.na(beta2)),
  p.value23 = t.test(beta2, beta3, paired = T)$p.value, 
  diff23 = diff(t.test(beta2, beta3)$estimate), 
  n23 = sum(!is.na(beta2) & !is.na(beta3)),
  p.value13 = t.test(beta1, beta3, paired = T)$p.value, 
  diff13 = diff(t.test(beta1, beta3)$estimate), 
  n13 = sum(!is.na(beta1) & !is.na(beta3))
  )

df_cpg_2s_test$adj.p.value12 = tdhia::custom_p.adjust(p = df_cpg_2s_test$p.value12, method = "fdr", n = 76)
df_cpg_2s_test$adj.p.value23 = tdhia::custom_p.adjust(p = df_cpg_2s_test$p.value23, method = "fdr", n = 76)
df_cpg_2s_test$adj.p.value13 = tdhia::custom_p.adjust(p = df_cpg_2s_test$p.value13, method = "fdr", n = 76)

plot_time_consistency_dual(df_cpg_2s_test, meth = "cpg")
  





# ICR Analysis                                                       ##########
#_______________________________________________________________________________

# Calculate beta slope for each patient, for each cpg_id, across time.
# Also grab timepoints 1, 2, and 3 for discrete tests
df_icr_slopes <- df_icr_long %>% 
  filter( ship_id %in% (filt_study_metadata %>% filter(total_tp> 2) %>% pull(ship_id))) %>% 
  filter(icr_id %in% top100_icrs) %>%
  group_by(ship_id, icr_id) %>%
  summarize(mean_age = mean(child_age_months, na.rm = TRUE),
            # cor_coef = my_cor.test(x = child_age_months/12, y = beta)$estimate,
            # cor_p_value = my_cor.test(x = child_age_months/12, y = beta)$p.value,
            # regr_coef = my_regress(x = child_age_months/12, y = beta)$estimate,
            # regr_p_value = my_regress(x = child_age_months/12, y = beta)$p.value,
            x = paste0(paste(child_age_months, collapse = ",")),
            y =  paste0(paste(beta, collapse = ",")),
            beta1 = get_value_in_range(beta, child_age_months, age_target = tps[1], age_range = 6, NULL)$value,
            time1 = get_value_in_range(beta, child_age_months, age_target = tps[1], age_range = 6, NULL)$age_months,
            beta2 = get_value_in_range(beta, child_age_months, age_target = tps[2], age_range = 6, NULL)$value,
            time2 = get_value_in_range(beta, child_age_months, age_target = tps[2], age_range = 6, NULL)$age_months,
            beta3 = get_value_in_range(beta, child_age_months, age_target = tps[3], age_range = 6, NULL)$value,
            time3 = get_value_in_range(beta, child_age_months, age_target = tps[3], age_range = 6, NULL)$age_months,
            n = n())


df_icr_2s_test <-  df_icr_slopes %>% filter(icr_id %in% top100_icrs) %>%
  filter(!is.na(beta1) & !is.na(beta2)) %>%
  group_by(icr_id) %>% 
  summarize(
    p.value12 = t.test(beta1, beta2, paired = T)$p.value, 
    diff12 = diff(t.test(beta1, beta2)$estimate), 
    n12 = sum(!is.na(beta1) & !is.na(beta2)),  
    p.value23 = t.test(beta2, beta3, paired = T)$p.value, 
    diff23 = diff(t.test(beta2, beta3)$estimate), 
    n23 = sum(!is.na(beta2) & !is.na(beta3)),
    p.value13 = t.test(beta1, beta3, paired = T)$p.value, 
    diff13 = diff(t.test(beta1, beta3)$estimate), 
    n13 = sum(!is.na(beta1) & !is.na(beta3))  
    )
df_icr_2s_test$adj.p.value12 = tdhia::custom_p.adjust(p = df_icr_2s_test$p.value12, method = "fdr", n = 76)
df_icr_2s_test$adj.p.value23 = tdhia::custom_p.adjust(p = df_icr_2s_test$p.value23, method = "fdr", n = 76)
df_icr_2s_test$adj.p.value13 = tdhia::custom_p.adjust(p = df_icr_2s_test$p.value13, method = "fdr", n = 76)
plot_time_consistency_dual(df = df_icr_2s_test,meth = "icr")


# GEO Export
#_______________________________________________________________________________
if (FALSE) {
  
  geo_path = paste0(output_dir_path, "/geo_time_consistency")
  dir.create(geo_path, showWarnings = F,recursive = T)
  
  
  # Prep the metadata table
  geo_metadata <- filt_study_metadata %>%  rename(title = id, tissue = sample_type, sample_name = Chip)
  geo_metadata$molecule = "genomic DNA"
  geo_metadata$label= "biotin"
  geo_metadata$idat1 = paste0(geo_metadata$sample_name, "_Grn.idat")
  geo_metadata$idat2 = paste0(geo_metadata$sample_name, "_Red.idat")
  geo_metadata$tissue = str_replace(geo_metadata$tissue, "_Patient_ID", "")
  
  write.csv(x = geo_metadata, file = 
              paste0(geo_path, "/geo_metadata_time_consist.csv"), row.names = F, col.names = F)

  # Export probe beta values
  idat_out_path <- paste0(output_dir_path, "/data/time_consistency_study_probe.rda")
  no_merge <- tdhia_pipeline(
    idat_dir_paths = idat_dir_paths, probe_data_cache_path = idat_out_path,
    idat_basenames = filt_study_metadata$Chip, merge_replicates = NULL, 
    OVERWRITE_TEMP_DATA = F,db_flag = T)
  write.csv(x = no_merge$probe_beta$probe_beta_df %>% t() %>% as.data.frame() %>%
              rownames_to_column("ID_REF") %>% t(), file = 
              paste0(geo_path, "/processed_beta_values_probe.csv"),row.names = T)
  
  # Export cpg beta matrix
  write.csv(x = cpg_beta$cpg_beta_df %>% t() %>% as.data.frame() %>%
              rownames_to_column("ID_REF") %>% t(), file = 
              paste0(geo_path, "/processed_beta_values_cpg.csv"), row.names = F)
  
  
  # Convert patient_id list to two files
  used_idat_paths <- probe_beta$idat_filepaths[basename(probe_beta$idat_filepaths) %in%  geo_metadata$sample_name]
 
  
  # Copy files to output folder
  dir.create(paste0(geo_path, "/idats_time_consist/"), showWarnings = F, recursive = T)
  for (n in seq_along(probe_beta$idat_filepaths)) {
    file.copy(from = paste0(probe_beta$idat_filepaths[n], "_Grn.idat"), 
              to = paste0(geo_path, "/", 
                          basename(probe_beta$idat_filepaths[n]), "_Grn.idat"),overwrite = T)
    file.copy(from = paste0(probe_beta$idat_filepaths[n], "_Red.idat"), 
              to = paste0(geo_path, "/", 
                          basename(probe_beta$idat_filepaths[n]), "_Red.idat"), overwrite = T)
  }
  
}


