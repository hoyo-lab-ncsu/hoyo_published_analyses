if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)
library(readr)

# Create a temporary directory for GEO files
geo_id <- "GSE334366"
tmp <- file.path(tempdir(), geo_id)
dir.create(tmp, showWarnings = FALSE, recursive = TRUE)

# Download GEO supplementary files to temporary directory
getGEOSuppFiles(geo_id, baseDir = tmp)

# Check downloaded files
list.files(tmp, recursive = TRUE)

# Load processed CpG beta-value dataframes and study metadata
processed_cpg_df_umbilical <- read_csv(list.files(tmp, pattern = "umbilical_cord_beta_matrix.csv(\\.gz)?$", recursive = TRUE, full.names = TRUE)[1])
processed_cpg_df_peripheral <- read_csv(list.files(tmp, pattern = "peripheral_blood_beta_matrix.csv(\\.gz)?$", recursive = TRUE, full.names = TRUE)[1])
study_data <- read_csv(list.files(tmp, pattern = "\\Supplementary_data_obesity.csv(\\.gz)?$", recursive = TRUE, full.names = TRUE)[1])

# Load raw IDATs
raw_tar <- list.files(tmp, pattern = "GSE334366_RAW\\.tar$", recursive = TRUE, full.names = TRUE)[1]
idat_dir_path <- file.path(tmp, "IDATs")
dir.create(idat_dir_path, showWarnings = FALSE, recursive = TRUE)
untar(raw_tar, exdir = idat_dir_path)
idat_files <- list.files(idat_dir_path, pattern = "\\.idat(\\.gz)?$", recursive = TRUE, full.names = TRUE)
