# List all files with any extension (so no folders)
matrix_files <- list.files(MATRIX_FOLDER, pattern = "\\.")

# Read in covariates data
patient_images <- read_csv(file.path(DATA_FOLDER, "merged_covariates.csv"))
