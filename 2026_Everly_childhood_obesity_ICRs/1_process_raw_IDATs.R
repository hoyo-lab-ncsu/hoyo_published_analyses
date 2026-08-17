# User: First complete 0_download_GEO_data.R. This script relies on the downloaded data

# Load libraries
library(tdhia)

# CRITICAL NOTE:
# The original preprocessing was performed on the complete cohort (approximately 1,900 samples).
# Probe filtering (removal of probes with detection p-values failing in >10% of samples) and sample
# filtering (removal of samples with >5% failed probes) were therefore based on the full dataset.

# To protect participant privacy and minimize data sharing, the repository includes only the IDAT files 
# for the 586 participants included in the present analysis. Consequently, rerunning the preprocessing 
# pipeline from the shared IDAT files will result in slightly different probe and sample filtering because 
# the percentage-based thresholds are calculated using a different number of samples.

# The processed CpG beta matrix used in the manuscript is provided so that all downstream analyses reproduce 
# the published results exactly.


# Process IDAT Files (please read CRITICAL NOTE first)
# ============================================================================================================
# a)  Load IDATS and convert to probe beta matrix
# Must provide path where the IDAT files are located on your computer.
probe_beta <- load_idata_to_probes(idat_dir_paths = idat_dir_path,
                                   multicore = TRUE,
                                   quantile_norm = FALSE,
                                   enforce_req_idats = FALSE,
                                   mask = FALSE,
                                   db_flag = FALSE)


# b) Filter probes that are not mapped and discard poor signal
filt_probe_beta <- filter_probes(probe_beta = probe_beta,
                                 discard_unmapped_probes =TRUE ,
                                 max_sig_pval = 0.2,
                                 max_probe_fail_rate = 0.1,
                                 max_patient_fail_rate = 0.05,
                                 db_flag = TRUE,
                                 set_failed_betas_na = FALSE,
                                 discard_failed_probes = TRUE,
                                 discard_failed_patients = TRUE)

# c)  Convert probe beta matrix to a cpg beta matrix
cpg_beta_test <- convert_probes_to_cpgs(probe_beta= filt_probe_beta,
                                        quantile_norm = FALSE,
                                        discard_unmapped_cpgs = TRUE,
                                        discard_non_icr_cpgs = TRUE,
                                        smooth_adj_cpgs = FALSE,
                                        sort_cpgs = TRUE,
                                        db_flag = FALSE)


predictor_cpg_df <- as.data.frame(t(cpg_beta_test[["cpg_beta_df"]]))
