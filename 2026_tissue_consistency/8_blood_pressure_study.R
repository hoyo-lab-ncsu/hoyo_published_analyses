






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

idat_dir_paths = paste0(github_path, "/data/IDAT_files 15Nov23")
output_dir_path = paste0(github_path, "/data/study_tissue_consistency")
study_dir_path <- paste0(github_path, "/data/study_tissue_consistency/study_data")
ref_path = paste0(github_path, "/tdhia_scripts/study_tissue_consistency")
source(paste0(ref_path, "/imprintome_toolbox.R"))

# Perform association study with
# 1) HUVECs
# 2) CD13
# And compare hits, direction of hits

df_icr_ranks <- read_rds(file = paste0(output_dir_path, "/ICR_ranks_diff_CD14_HUVEC.rds"))


# load study data
# df_study <- haven::read_sas(paste0(output_dir_path, "/ship_bp_hr_06dec23.sas7bdat"))

study_data <- haven::read_sas(paste0(output_dir_path, "/study_data/ship_imprint_tissue_consistency.sas7bdat"))
study_data <- study_data %>% filter(!is.na(syspct) & !is.na(diaspct))
attr(study_data, "format.sas") <- NULL
study_data$use_cig_preg[(is.na(study_data$use_cig_preg)) & (study_data$smoke_ever_100 == 0)] <- 0


hist(study_data$diaspct,xlab = "Diastolic Pct",main = "")
hist(study_data$syspct,xlab = "Systolic Pct",main = "")

# Extract rowss for HUVEC
# study_huvec <- study_data %>% filter(HUVEC_Patient_ID != "")
# # Extract rows for CD14
# study_cd14 <- study_data %>% filter(CD14_Patient_ID != "")


# HUVEC idat processing                                          ###############
#_______________________________________________________________________________
data_bps = list()
# 1) Load IDATS and convert to probe beta matrix
# Output: rows: probe_id     x     columns: patient
out_file <- paste0(output_dir_path, "/data/HUVEC-CD14_blood_pressure.rda")
if (!file.exists(out_file)) {
  data_bps$probe_beta <-
    load_idata_to_probes(idat_dir_paths = idat_dir_paths, multicore = FALSE,
                         idat_basenames = c(study_data$HUVEC_Patient_ID[1:nrow(study_data)],
                                            study_data$CD14_Patient_ID[1:nrow(study_data)]),
                         merge_replicates = "pre_beta",
                         quantile_norm = FALSE, db_flag = FALSE, enforce_idat_names = FALSE)
  save(data_bps, file = out_file)
} else {load(out_file)}
# 2) Filter probes that are not mapped and discard poor signal
data_bps$filt_probe_beta <-
  filter_probes(probe_beta = data_bps$probe_beta, discard_unmapped_probes = TRUE ,
                max_sig_pval = 0.2, set_failed_betas_na = TRUE, 
                max_probe_fail_rate = 0.25, min_design_score = NA, 
                discard_failed_probes = TRUE, db_flag = FALSE)
#3  Convert probe beta matrix to a cpg beta matrix
data_bps$cpg_beta <- convert_probes_to_cpgs(data_bps$filt_probe_beta, quantile_norm = FALSE, 
                                   db_flag = TRUE, smooth_adj_cpgs = FALSE)
#4  Convert probe beta matrix to an icr beta matrix
data_bps$icr_beta <- convert_cpgs_to_icrs(data_bps$cpg_beta, max_icr_fail_rate = 0.2)




# Association Studies Vs BP Percentile                              ############
#_______________________________________________________________________________
# Adjust p-=values to number of ICRs (not # of cpgs)
n_p_adj = nrow(data_bps$icr_beta$icr_beta_df)

n.cores = max(c(parallel::detectCores()-2,1))

# 1) HUVEC

# By Sex
# By Races
# systolic_bp_percentile ~ cpg +  maternal_bmi  +  maternal_preg_smoke  + maternal_education_level



custom_model <- function(df, response_var, patient_id_colname , filt_var = NA, filt_val = NA) {
  df$filter_bv = rep(T, nrow(df))
  if (!is.na(filt_var)) { df$filter_bv = df[[filt_var]] == filt_val } 
  
  
  covariate_colnames <- c(patient_id_colname, "mat_bmi_lmp", "use_cig_preg", "mat_bt_your_edu_lev")
  # if (filt_var == "base_sex") covariate_colnames = c(covariate_colnames, "race_final")
  # if (filt_var == "race_final") covariate_colnames = c(covariate_colnames, "ba_sex")
  
  out <- analyze_association(
    R = select(df %>% filter(filter_bv), patient_id_colname, response_var) %>% 
      column_to_rownames(patient_id_colname),
    P = as.data.frame(t(data_bps$cpg_beta$cpg_beta_df %>% select(
      df %>% filter(filter_bv) %>% pull(all_of(patient_id_colname)) ))),
    Pe = select(df %>% filter(filter_bv), all_of(covariate_colnames)) %>% 
      column_to_rownames(patient_id_colname),
    family = "gaussian", n_p_adj = n_p_adj, db_flag = FALSE,
    verbose = TRUE, impute_na = TRUE, rm.na.R = TRUE,
    max_p_val = 0.05, n.cores = n.cores) 
  return(out)
}



out_path <- paste0(output_dir_path, "/data/stats_syspct_HUVEC.rda")
if (!file.exists(out_path)) {
  huvec_stat = list()

  huvec_stat$sys_all <- custom_model(
    df = study_data %>% filter(HUVEC_Patient_ID != ""),
    response_var = "syspct", patient_id_colname = "HUVEC_Patient_ID")
  huvec_stat$sys_female <- custom_model(
    df = study_data %>% filter(HUVEC_Patient_ID != ""),
    response_var = "syspct", patient_id_colname = "HUVEC_Patient_ID",
    filt_var = "ba_sex", filt_val = 1)
  huvec_stat$sys_male <- custom_model(
    df = study_data %>% filter(HUVEC_Patient_ID != ""),
    response_var = "syspct", patient_id_colname = "HUVEC_Patient_ID",
    filt_var = "ba_sex", filt_val = 2)
  huvec_stat$sys_black <- custom_model(
    df = study_data %>% filter(HUVEC_Patient_ID != ""),
    response_var = "syspct", patient_id_colname = "HUVEC_Patient_ID",
    filt_var = "race_final", filt_val = "Black")
  huvec_stat$sys_hispanic <- custom_model(
    df = study_data %>% filter(HUVEC_Patient_ID != ""),
    response_var = "syspct", patient_id_colname = "HUVEC_Patient_ID",
    filt_var = "race_final", filt_val = "Hispanic")
  huvec_stat$sys_white <- custom_model(
    df = study_data %>% filter(HUVEC_Patient_ID != ""),
    response_var = "syspct", patient_id_colname = "HUVEC_Patient_ID",
    filt_var = "race_final", filt_val = "White")
  huvec_stat$sys_other <- custom_model(
    df = study_data %>% filter(HUVEC_Patient_ID != ""),
    response_var = "syspct", patient_id_colname = "HUVEC_Patient_ID",
    filt_var = "race_final", filt_val = "Other")

  huvec_stat$dia_all <- custom_model(
    df = study_data %>% filter(HUVEC_Patient_ID != ""),
    response_var = "diaspct", patient_id_colname = "HUVEC_Patient_ID")
  huvec_stat$dia_female <- custom_model(
    df = study_data %>% filter(HUVEC_Patient_ID != ""),
    response_var = "diaspct", patient_id_colname = "HUVEC_Patient_ID",
    filt_var = "ba_sex", filt_val = 1)
  huvec_stat$dia_male <- custom_model(
    df = study_data %>% filter(HUVEC_Patient_ID != ""),
    response_var = "diaspct", patient_id_colname = "HUVEC_Patient_ID",
    filt_var = "ba_sex", filt_val = 2)
  huvec_stat$dia_black <- custom_model(
    df = study_data %>% filter(HUVEC_Patient_ID != ""),
    response_var = "diaspct", patient_id_colname = "HUVEC_Patient_ID",
    filt_var = "race_final", filt_val = "Black")
  huvec_stat$dia_hispanic <- custom_model(
    df = study_data %>% filter(HUVEC_Patient_ID != ""),
    response_var = "diaspct", patient_id_colname = "HUVEC_Patient_ID",
    filt_var = "race_final", filt_val = "Hispanic")
  huvec_stat$dia_white <- custom_model(
    df = study_data %>% filter(HUVEC_Patient_ID != ""),
    response_var = "diaspct", patient_id_colname = "HUVEC_Patient_ID",
    filt_var = "race_final", filt_val = "White")
  huvec_stat$dia_other <- custom_model(
    df = study_data %>% filter(HUVEC_Patient_ID != ""),
    response_var = "diaspct", patient_id_colname = "HUVEC_Patient_ID",
    filt_var = "race_final", filt_val = "Other")

  save(huvec_stat, file = out_path)
} else {
  load(out_path)
}




out_path <- paste0(output_dir_path, "/data/stats_syspct_CD14.rda")
if (!file.exists(out_path)) {
  cd14_stat = list()
  
  cd14_stat$sys_all <- custom_model(
    df = study_data %>% filter(CD14_Patient_ID != ""), 
    response_var = "syspct", patient_id_colname = "CD14_Patient_ID")
  cd14_stat$sys_male <- custom_model(
    df = study_data %>% filter(CD14_Patient_ID != ""), 
    response_var = "syspct", patient_id_colname = "CD14_Patient_ID", 
    filt_var = "ba_sex", filt_val = 1)
  cd14_stat$sys_female <- custom_model(
    df = study_data %>% filter(CD14_Patient_ID != ""), 
    response_var = "syspct", patient_id_colname = "CD14_Patient_ID", 
    filt_var = "ba_sex", filt_val = 2)
  cd14_stat$sys_black <- custom_model(
    df = study_data %>% filter(CD14_Patient_ID != ""),
    response_var = "syspct", patient_id_colname = "CD14_Patient_ID", 
    filt_var = "race_final", filt_val = "Black")
  cd14_stat$sys_hispanic <- custom_model(
    df = study_data %>% filter(CD14_Patient_ID != ""),
    response_var = "syspct", patient_id_colname = "CD14_Patient_ID", 
    filt_var = "race_final", filt_val = "Hispanic")
  cd14_stat$sys_white <- custom_model(
    df = study_data %>% filter(CD14_Patient_ID != ""),
    response_var = "syspct", patient_id_colname = "CD14_Patient_ID", 
    filt_var = "race_final", filt_val = "White")
  cd14_stat$sys_other <- custom_model(
    df = study_data %>% filter(CD14_Patient_ID != ""),
    response_var = "syspct", patient_id_colname = "CD14_Patient_ID", 
    filt_var = "race_final", filt_val = "Other")
  
  cd14_stat$dia_all <- custom_model(
    df = study_data %>% filter(CD14_Patient_ID != ""), 
    response_var = "diaspct", patient_id_colname = "CD14_Patient_ID")
  cd14_stat$dia_male <- custom_model(
    df = study_data %>% filter(CD14_Patient_ID != ""),
    response_var = "diaspct", patient_id_colname = "CD14_Patient_ID", 
    filt_var = "ba_sex", filt_val = 1)
  cd14_stat$dia_female <- custom_model(
    df = study_data %>% filter(CD14_Patient_ID != ""), 
    response_var = "diaspct", patient_id_colname = "CD14_Patient_ID", 
    filt_var = "ba_sex", filt_val = 2)
  cd14_stat$dia_black <- custom_model(
    df = study_data %>% filter(CD14_Patient_ID != ""),
    response_var = "diaspct", patient_id_colname = "CD14_Patient_ID", 
    filt_var = "race_final", filt_val = "Black")
  cd14_stat$dia_hispanic <- custom_model(
    df = study_data %>% filter(CD14_Patient_ID != ""),
    response_var = "diaspct", patient_id_colname = "CD14_Patient_ID", 
    filt_var = "race_final", filt_val = "Hispanic")
  cd14_stat$dia_white <- custom_model(
    df = study_data %>% filter(CD14_Patient_ID != ""),
    response_var = "diaspct", patient_id_colname = "CD14_Patient_ID", 
    filt_var = "race_final", filt_val = "White")
  cd14_stat$dia_other <- custom_model(
    df = study_data %>% filter(CD14_Patient_ID != ""),
    response_var = "diaspct", patient_id_colname = "CD14_Patient_ID", 
    filt_var = "race_final", filt_val = "Other")
  
  
  save(cd14_stat, file = out_path)
} else {
  load(out_path)
}

df_imp_all <- rbind(
  do.call(rbind,sapply(names(huvec_stat), FUN = function(x) huvec_stat[[x]]$imp_site %>% mutate(group = x, tissue = "huvec",cpg_id = Variable),
                       simplify = FALSE, USE.NAMES = TRUE)),
  do.call(rbind,sapply(names(cd14_stat), FUN = function(x) cd14_stat[[x]]$imp_site %>% mutate(group = x, tissue = "cd14", cpg_id = Variable),
                       simplify = FALSE, USE.NAMES = TRUE)) )

# a<-do.call(rbind,sapply(names(huvec_stat), FUN = function(x) huvec_stat[[x]]$imp_site %>% mutate(group = x, tissue = "huvec"),
#                         simplify = FALSE, USE.NAMES = TRUE))
# a <- filter(select(a, -Formula), ADJ_P_VAL < 0.05)
# b<-do.call(rbind,sapply(names(cd14_stat), FUN = function(x) cd14_stat[[x]]$imp_site %>% mutate(group = x, tissue = "huvec"),
#                         simplify = FALSE, USE.NAMES = TRUE))
# b <- filter(select(b, -Formula), ADJ_P_VAL < 0.05)



# Select cpgs from topN ICRs based on diff results, then correct p-values based on ICRs, return 
nranks = 10
icr_ranks <- read_rds(file = paste0(output_dir_path, "/", "ICR_ranks_diff_CD14_HUVEC.rds"))
top_icr_ids <- str_replace(string = icr_ranks$icr_name[1:nranks], " ", "_")

# # Filter results from particular model
# df_sig_sub <- df_imp_all %>% filter(group == "dia_all" & tissue=="cd14") %>%
#   mutate(icr_id = add_metadata_to_imp_sites(cpg_id, imp_type = "cpg")$icr_id) %>%
#   filter(icr_id %in% top_icr_ids) %>%
#   mutate(ADJ_P_VAL = custom_p.adjust(p = P_VAL, n = nranks)) %>% filter(ADJ_P_VAL <0.05)
# 
# df_sig_sub$ADJ_P_VAL


df_imp_all %>% filter(group == "dia_all" & tissue=="cd14") %>%
  mutate(icr_id = add_metadata_to_imp_sites(cpg_id, imp_type = "cpg")$icr_id) %>%
  filter(icr_id %in% top_icr_ids) %>%
  mutate(ADJ_P_VAL = custom_p.adjust(p = P_VAL, n = nranks)) %>% filter(ADJ_P_VAL <0.05)

df_imp_all %>% filter(group == "dia_all" & tissue=="huvec") %>%
  mutate(icr_id = add_metadata_to_imp_sites(cpg_id, imp_type = "cpg")$icr_id) %>%
  filter(icr_id %in% top_icr_ids) %>%
  mutate(ADJ_P_VAL = custom_p.adjust(p = P_VAL, n = nranks)) %>% filter(ADJ_P_VAL <0.05)

# sum(df_sig_sub$ADJ_P_VAL < 0.05)


# Filter Significance ICRs and add nearest gene
#_______________________________________________________________________________

# Significant ICRs
df_imp_sig <- filter(dplyr::select(df_imp_all, -Formula), ADJ_P_VAL < 0.05)
df_imp_sig<-cbind(df_imp_sig,add_metadata_to_imp_sites(df_imp_sig$Variable, imp_type = "cpg"))
df_imp_sig$icr_name = str_replace(df_imp_sig$icr_id, "_"," ")

df_imp_sig$group_label <- 
  str_to_title(apply(X = cbind(str_split(df_imp_sig$group,"_", simplify = TRUE),
                               df_imp_sig$tissue)[,c(1,3,2)], MARGIN = 1, FUN = paste,
                     collapse=" "))
# paste0(df_imp_sig$group," ", toupper(df_imp_sig$tissue))

# a<-cbind(str_split(df_imp_sig$group,"_", simplify = TRUE),df_imp_sig$tissue)




# do.call("paste",a)
# df_imp_sig$icr_id <- sapply(df_imp_sig$Variable, function(x) 
#   tdhia::mapping_cpg_icr_ids$ICR_id[x==tdhia::mapping_cpg_icr_ids$CpG_id][1])
# df_imp_sig$icr_name <- str_replace(df_imp_sig$icr_id, "_", " ")
# df_imp_sig <- cbind(df_imp_sig,add_metadata_from_cpg(df_imp_sig$Variable))

df_imp_sig_out <- df_imp_sig %>%  dplyr::select(c("group_label", "icr_name", "Response", "Estimate", "ADJ_P_VAL", 
                                           "Genomic.Coordinates", "Nearest.Transcript",
                                           "Distance.to.Nearest.Transcript")) 

# For each ICR and Response combination, record info, and:
#  1) Minimum adjusted p-value across cpg sites (group_by, summarize)
#  2) Max magnitude of model estimate
df_icr_sig <- df_imp_sig_out %>%
  group_by(icr_name, group_label) %>%
  summarise(ADJ_P_VAL = min(ADJ_P_VAL), 
            Estimate = max(abs(Estimate)) * sign(Estimate[which.max(abs(Estimate))])) %>% 
  merge(y = dplyr::select(df_imp_sig_out,-c(Response,group_label, Estimate, ADJ_P_VAL)) %>% distinct(), 
        by = "icr_name") %>% arrange(desc(group_label))

df_icr_sig$Freq <- sapply( df_icr_sig$icr_name, function(x) sum(x==df_icr_sig$icr_name))
df_icr_sig$Responses.Shared <- sapply( df_icr_sig$icr_name, function(x) 
  paste(df_icr_sig$group_label[x==df_icr_sig$icr_name], collapse = ", "))

# # add min-pvalue, max(estimate), genomic position
write.csv(x =  df_icr_sig %>% arrange(group_label, ADJ_P_VAL),
          file = paste0(output_dir_path, paste0(
            "/significant_icr-genes-shared.csv")))

# Count number of significance ICRs in various groupings
df_icr_summary <- df_imp_sig %>% group_by(tissue, Response, group) %>% 
  summarize(
    n_cpg = length(unique(Variable)),
    n_icr = length(unique(icr_id)), 
    n_zinc = length(unique(icr_id[is_icr_zinc])),
    n_highconf = length(unique(icr_id[icr_conf == 1])),
    n_medconf = length(unique(icr_id[icr_conf == 2])),
    nlowconf = length(unique(icr_id[icr_conf == 3])))
df_icr_summary




plot_venn_shared_icrs(df_imp_sig, colname = "group_label", group_names = NULL,
                      output_dir_path = output_dir_path, name_suffix = "SysDia_Groups") 

table(study_data$ba_sex)
table(study_data$race_final)



