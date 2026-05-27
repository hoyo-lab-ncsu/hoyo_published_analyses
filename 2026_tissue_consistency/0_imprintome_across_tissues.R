



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
library(preprocessCore)
library(WriteXLS)
devtools::load_all(paste0(dirname(dirname(getwd())), "/tdhia"))


# Change these file paths for folder with IDAT files and output folder
# Define local variables
github_path <- dirname(dirname(getwd()))

idat_dir_paths = paste0(github_path, "/data/IDAT_files 15Nov23")
output_dir_path = paste0(github_path, "/data/study_tissue_consistency")
study_dir_path <- paste0(github_path, "/data/study_tissue_consistency/study_data")
ref_path = paste0(github_path, "/hoyo_published_analyses/2026_tissue_consistency")
support_data_path <- paste0(ref_path, "/supporting_data")
source(paste0(ref_path, "/imprintome_toolbox.R"))
source(paste0(ref_path, "/imprintome_toolbox_old.R"))

# df_850k <- read.csv(file = paste0(ref_path, "/EPIC-8v2-0_A2.csv"))
# unq_cpg_850k <- unique(str_replace(df_850k$Illumina, "_.*$",""))
# save(unq_cpg_850k, file = paste0(ref_path, "/unq_cpg_850k.RData"))
# write.csv(x = select(df_850k, Illumina), file = "850k_probes.csv")





# Metadata for study (excluding patient ID for true diagnostic imprintome array)
study_data <- read.csv(paste0(study_dir_path, "/20240507_SHIP_Imprint_Tissue_Consistency.csv"))
test_study <- haven::read_sas(paste0(study_dir_path, "/nest_metals_bmi_bp_bruce.sas7bdat"))

# Convert ship_id column to study data.
rownames(study_data) <- study_data$ship_id

# Coerce into factors
study_data$race_final <- as.factor(study_data$race_final)
study_data$ba_sex <- as.factor(study_data$ba_sex)
study_data$smoke_ever_100 <- as.factor(study_data$smoke_ever_100)
# Fill in 
study_data$use_cig_preg[study_data$smoke_ever_100==0] = 0



cat(sprintf("Patients in study file: %.0f\n", nrow(study_data)))
cat(sprintf("HUVEC :: CD14: %.0f\n", 
            sum(study_data$HUVEC_Patient_ID != "" & study_data$CD14_Patient_ID!="")))

cat(sprintf("HUVEC :: Placenta: %.0f\n", 
            sum(study_data$HUVEC_Patient_ID != "" & study_data$Placental_Patient_ID!="")))
cat(sprintf("CD14 :: Placenta: %.0f\n", 
            sum(study_data$CD14_Patient_ID!="" & study_data$Placental_Patient_ID!="")))
cat(sprintf("All three: %.0f\n", sum(study_data$HUVEC_Patient_ID!="" &
                                     study_data$CD14_Patient_ID!="" &
                                     study_data$Placental_Patient_ID!="")))


# FIltering metadata for correlation study                      ##########
#_______________________________________________________________________________
filt_study_data <- filter(study_data, HUVEC_Patient_ID!="" & CD14_Patient_ID!="" &
                        Placental_Patient_ID!="")


filt_study_data$mat_bmi_hi <- factor(filt_study_data$mat_bmi_lmp > 40, 
                                     levels = c(FALSE, TRUE), 
                                     labels = c("Low", "High"),
                                     ordered = TRUE)
filt_study_data$mat_aces_hi <- factor(filt_study_data$Impute_aces_mom_score > 2, 
                                      levels = c(FALSE, TRUE), 
                                      labels = c("Low", "High"),
                                      ordered = TRUE)
  
  
# Grab required metadata from study
df_idat_sources <- filt_study_data %>% dplyr::select(ship_ID, HUVEC_Patient_ID, 
                          CD14_Patient_ID, Placental_Patient_ID, mat_bmi_hi, mat_aces_hi)

# list all idat_basenames 
df_all_idats <- pivot_longer(df_idat_sources, col = ends_with("_Patient_ID"),
                             names_to = "source", values_to = "patient_id")


# Read list of ICRs that overlap with nearby zing finger
df_zinc_finger <-read.csv(paste0(support_data_path, "/ICR_ZNF_1000bp.csv"))
df_zinc_finger$icr_id <- as.numeric(gsub("ICR_([0-9]+).*","\\1",df_zinc_finger$icr))
df_zinc_finger$near_zinf_finger <- rowSums(!df_zinc_finger[, 3:5] == "", na.rm = TRUE) > 0
zinc_finger_icrs <- df_zinc_finger$icr_id[df_zinc_finger$near_zinf_finger]





# Idat processing                                          #####################
#_______________________________________________________________________________
data_tri = list()
# 1) Load IDATS and convert to probe beta matrix
# Output: rows: probe_id     x     columns: patient
idat_out_path <- paste0(output_dir_path, "/data/consistency_study.rda")
if (!file.exists(idat_out_path)) {
  data_tri$probe_beta <-
    load_idata_to_probes(idat_dir_paths = idat_dir_paths, multicore = FALSE,
                         idat_basenames = df_all_idats$patient_id, merge_replicates = "pre_beta",
                         quantile_norm = FALSE, db_flag = F, enforce_idat_names = FALSE)
  save(data_tri, file = idat_out_path)
} else {load(idat_out_path)}


# 2) Filter probes that are not mapped and discard poor signal
# Output: rows: probe_id     x     columns: patient
data_tri$filt_probe_beta <-
  filter_probes(probe_beta = data_tri$probe_beta, discard_unmapped_probes = TRUE ,
                max_sig_pval = 0.2, set_failed_betas_na = FALSE, max_patient_fail_rate = 1,
                max_probe_fail_rate = 0.25, min_design_score = NA, 
                discard_failed_probes = TRUE, db_flag = F)


#3  Convert probe beta matrix to a cpg beta matrix
# Output: rows: cpg_id     x     columns: patient
data_tri$cpg_beta <- convert_probes_to_cpgs(data_tri$filt_probe_beta, quantile_norm = FALSE, 
                                   db_flag = F, smooth_adj_cpgs = FALSE)


#4  Convert probe beta matrix to an icr beta matrix
# Output: rows: icr_id     x     columns: patient
data_tri$icr_beta <- convert_cpgs_to_icrs(data_tri$cpg_beta, max_icr_fail_rate = 0.2)
# icr_beta$icr_hbeta_df <- abs(0.5 - icr_beta$icr_beta_df)


# filt_study_data <- data_tri$cpg_beta$cpg_beta_df %>% filter()


# Statistical analysis setup                                              ######
#_______________________________________________________________________________


# Load annotated list of whole imprintome
# ICR_IDs that end with "#" are high confidence
imp_whole <- read.csv(paste0(support_data_path, "/whole_imprintome_table.csv")) 
imp_whole$icr_id <- as.numeric(gsub("ICR_([0-9]+).*","\\1",imp_whole$ID))
imp_whole$icr_name <- paste0("ICR ", imp_whole$icr_id)
  
# Get list of ICRs that have evidence of gametic origin for methylation
# Includes "high confidence" (lit validated ICRs)
imp_gamete <- read.csv(paste0(support_data_path, "/table1_icrs_gametic_methyl_origin.csv"))
imp_gamete$icr_id <- as.numeric(gsub("ICR_([0-9]+).*","\\1",imp_gamete$ID))
imp_gamete$icr_name <- paste0("ICR ", imp_gamete$icr_id)


icr_ids_gametic <- as.numeric(gsub("ICR_([0-9]+).*","\\1",tdhia::imprintome_gametic_icrs$icr_id))
# 
# # Scan for previously published icrs
high_conf_icrs <- imp_whole$icr_id[grep('#', imp_whole$ID)]
# # Medium confidence include high confidence
med_conf_icrs <- icr_ids_gametic #setdiff(imp_gamete$icr_id, high_conf_icrs)
low_conf_icrs <- setdiff(imp_whole$icr_id, union(med_conf_icrs, high_conf_icrs))


# # Get numeric ICR IDs
# icr_id <- as.numeric(stringr::str_replace(rownames(icr_beta$icr_beta_df), "ICR_", ""))
# 
# # Create confidence groups for icrs
# icr_conf <- as.numeric(icr_id %in% low_conf_icrs)*3
# icr_conf[icr_id %in% med_conf_icrs] <- 2
# icr_conf[icr_id %in% high_conf_icrs] <- 1
# names(icr_conf) <- icr_id



# List of ICRs and their ICR confidence level
#_______________________________________________________________________________
df_icr_conf <- data.frame(icr_id = rownames(data_tri$icr_beta$icr_beta_df))
df_icr_conf$icr_num <- as.numeric(stringr::str_replace(
  df_icr_conf$icr_id, "ICR_", ""))
df_icr_conf$icr_conf <- as.numeric(df_icr_conf$icr_num %in% low_conf_icrs)*3
df_icr_conf$icr_conf[df_icr_conf$icr_num %in% med_conf_icrs] <- 2
df_icr_conf$icr_conf[df_icr_conf$icr_num %in% high_conf_icrs] <- 1
# Add factor for icrs close to zinc finger proteins
df_icr_conf$zinc_finger <- df_icr_conf$icr_num %in% zinc_finger_icrs

table(df_icr_conf$icr_conf, df_icr_conf$zinc_finger)
colSums(table(df_icr_conf$icr_conf, df_icr_conf$zinc_finger))
rowSums(table(df_icr_conf$icr_conf, df_icr_conf$zinc_finger))
nrow(df_icr_conf)

# List of CpG and their ICR confidence level
#_______________________________________________________________________________
# Create vector of cpg conf groups for measured cpg sites (in order)
df_cpg_conf <- data.frame(CpG_id = rownames(data_tri$cpg_beta$cpg_beta_df))
df_cpg_conf <- left_join(x = df_cpg_conf,
          y = dplyr::select(data_tri$icr_beta$cpg_icr_mapping, c(ICR_id, CpG_id)), by = join_by(CpG_id), 
          keep = FALSE, multiple= "first") %>%
  rename(cpg_id = CpG_id, icr_id = ICR_id)
df_cpg_conf$icr_num <- as.numeric(str_replace(df_cpg_conf$icr_id, "ICR_","" ))
df_cpg_conf = left_join(df_cpg_conf, dplyr::select(df_icr_conf, c(icr_num, icr_conf)), by = join_by(icr_num))
# Add factor for cpgs close to zinc finger proteins
df_cpg_conf$zinc_finger <- df_cpg_conf$icr_num %in% zinc_finger_icrs

table(df_cpg_conf$icr_conf, df_cpg_conf$zinc_finger)
colSums(table(df_cpg_conf$icr_conf, df_cpg_conf$zinc_finger))
rowSums(table(df_cpg_conf$icr_conf, df_cpg_conf$zinc_finger))
nrow(df_cpg_conf)

# Go through each patient, calculate correlation between tissues
fthresh <- function (x) x<0.65 & x>0.35


# Assemble beta matrix for each tissue type, where patients are ordered the sample
cpg_long <- beta_mat_to_long(data_tri$cpg_beta$cpg_beta_df, filt_study_data, db_flag = T)


## CORSIV sites                                                              ##### 
df_corsiv <- read.csv(file = paste0(support_data_path, "/CoRSIV_incl_imprintomeoverlap.csv"))
imp_cpg_corsiv <- df_corsiv %>% filter(!(is.na(In.Imprintome))) %>% pull(In.Imprintome)

## 850K cpg sites                                                           #####
# readRDS(file = paste0(support_data_path, "/unq_cpg_850k.rds"))
# saveRDS(object = unq_cpg_850k, file = paste0(ref_path, "/unq_cpg_850k.rds"))
unq_cpg_850k <- readRDS(file = paste0(support_data_path, "/unq_cpg_850k.rds"))
imp_cpg_850k <- unq_cpg_850k[unq_cpg_850k %in% rownames(data_tri$cpg_beta$cpg_beta_df)]


# Add membership for corsiv and 850k
cpg_long <- cpg_long %>% mutate(is_850k = cpg_long$cpg_id %in% (unq_cpg_850k))
cpg_long <- cpg_long %>% mutate(is_corsiv = cpg_long$cpg_id %in% (imp_cpg_corsiv))


# Export table of cpg metadata
cpg_info_export <-   cpg_long %>% select(c("cpg_id", "is_icr_zinc", "is_corsiv", "is_850k", "icr_id", "icr_conf")) %>%
  distinct() %>% arrange(str_replace(icr_id,"ICR_","") %>% as.numeric) %>% mutate(is_icr_zinc = paste0("'", as.character(is_icr_zinc)))
write.csv(x = cpg_info_export, file = paste0(output_dir_path, "/CpG_Metadata.csv"))


# ICR level mfor plotting
icr_long <- cpg_long %>% dplyr::group_by(icr_id, patient_id, tissue) %>%
  summarize(beta = mean(beta), icr_conf = icr_conf[1], 
            is_icr_zinc = is_icr_zinc[1], mat_bmi_hi = mat_bmi_hi[1], 
            mat_aces_hi  = mat_aces_hi[1], is_corsiv = any(is_corsiv))


# GEO File Export
#_______________________________________________________________________________
if (FALSE) {
  
  geo_path = paste0(output_dir_path, "/geo_tissue_compare")
  dir.create(geo_path, showWarnings = F,recursive = T)
  

  
  # mom_ship_id ship_ID 
  geo_metadata <- filt_study_data %>% pivot_longer(
    cols = c("HUVEC_Patient_ID", "CD14_Patient_ID", "Placental_Patient_ID"), names_to = "tissue", values_to = "sample_name")
  geo_metadata <- geo_metadata %>% dplyr::select("ship_ID", "race_final", "mat_bmi_hi", "mat_aces_hi", 
                                                 "tissue", "sample_name") %>%  rename(title = ship_ID)
  geo_metadata$molecule = "genomic DNA"
  geo_metadata$label= "biotin"
  geo_metadata$idat1 = paste0(geo_metadata$sample_name, "_Grn.idat")
  geo_metadata$idat2 = paste0(geo_metadata$sample_name, "_Red.idat")
  geo_metadata$tissue = str_replace(geo_metadata$tissue, "_Patient_ID", "")
  
  write.csv(x = geo_metadata, file = 
              paste0(geo_path, "/geo_metadata_consistency.csv"), row.names = F, col.names = F)
 
  
  # Export probe beta values
  idat_out_path <- paste0(output_dir_path, "/data/consistency_study_probe.rda")
  no_merge <- tdhia_pipeline(
    idat_dir_paths = idat_dir_paths, probe_data_cache_path = idat_out_path,
    idat_basenames = df_all_idats$patient_id, merge_replicates = NULL, 
    OVERWRITE_TEMP_DATA = F,db_flag = T)
  write.csv(x = no_merge$probe_beta$probe_beta_df %>% t() %>% as.data.frame() %>%
              rownames_to_column("ID_REF") %>% t(), file = 
              paste0(geo_path, "/processed_beta_values_probe.csv"),row.names = T)
  
  
  
  # Export cpg beta values
  write.csv(x = data_tri$cpg_beta$cpg_beta_df %>% t() %>% as.data.frame() %>%
             rownames_to_column("ID_REF") %>% t(), file = 
              paste0(geo_path, "/processed_beta_values_cpg.csv"),row.names = T)#,sep = ",")
  
  
  # Convert patient_id list to two files
  df_all_idats$patient_id
  
  # Copy files to output folder
  dir.create(paste0(geo_path, "/idats_tissue_consist/"), showWarnings = F, recursive = T)
  for (n in seq_along(data_tri$probe_beta$idat_filepaths)) {
    file.copy(from = paste0(data_tri$probe_beta$idat_filepaths[n], "_Grn.idat"), 
              to = paste0(geo_path, "/idats_tissue_consist/", 
                          basename(data_tri$probe_beta$idat_filepaths[n]), "_Grn.idat"),overwrite = T)
    file.copy(from = paste0(data_tri$probe_beta$idat_filepaths[n], "_Red.idat"), 
              to = paste0(geo_path, "/idats_tissue_consist/", 
                          basename(data_tri$probe_beta$idat_filepaths[n]), "_Red.idat"), overwrite = T)
  }
  
}






