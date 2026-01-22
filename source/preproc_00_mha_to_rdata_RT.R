library(here)
library(tidyverse)

data_path <- here::here("../../data/SphericalMaps_database/SphericalMaps_database")

outdir <- file.path(here("processed"))

# Insight Toolkit
# https://itk.org
# documentation: https://docs.itk.org/en/latest/
library(SimpleITK)

# Get list of files/folders in database
folder_list <- list.files(data_path)

# Number of folders
n_folder <- length(folder_list) - 5 #(subtract 4 csv files and mask)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Covariate Handling ----
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# we do this first so that we can check matches as we go

covariate_file_inds <- which(dir.exists(file.path(data_path, folder_list)) == FALSE)
covariate_files <- folder_list[covariate_file_inds]

file1 <- file.path(data_path, covariate_files[1])
cov1 <- read_csv(file1) %>%
  mutate(Study = "CHCO", source_file = covariate_files[1]) %>%
  rename(ID = "Image ID") %>%
  mutate(original_file_suffix = paste0(ID)) %>%
  mutate(Group = "Craniosynostosis")

file2 <- file.path(data_path, covariate_files[2])
cov2 <- read_csv(file2) %>%
  dplyr::select("Subject ID", "Age", "Sex") %>%
  mutate(Study = "CHCO", source_file = covariate_files[2]) %>%
  rename(ID = `Subject ID`) %>%
  mutate(original_file_suffix = paste0(ID)) %>%
  mutate(Group = "Normative")

file3 <- file.path(data_path, covariate_files[3])
cov3 <- read_csv(file3) %>%
  mutate(Study = "CNH", source_file = covariate_files[3]) %>%
  rename(ID = "Image ID") %>%
  mutate(original_file_suffix = paste0(as.numeric(ID))) %>%
  mutate(ID = paste0("CNH_Cranio_", sprintf("%05d", ID) ))  %>%
  mutate(Group = "Craniosynostosis")

file4 <- file.path(data_path, covariate_files[4])
cov4 <- read_csv(file4)  %>%
  mutate(Study = "CNH", source_file = covariate_files[4]) %>%
  rename(ID = "Subject ID")   %>%
  mutate(original_file_suffix = paste0(as.numeric(ID))) %>%
  mutate(ID = paste0("CNH_", ID))  %>%
  mutate(Group = "Normative")


covariates_temp <- bind_rows(cov1, cov2, cov3, cov4) %>%
  dplyr::rename_with(~ gsub(" ", "_", .x)) %>%
  dplyr::select(-c("...13")) %>%
  relocate(Group, .after = ID)

# Now, we need to match the covariates to the corresponding folder
covariates <- covariates_temp %>%
  mutate(
    ID_base = case_when(
      str_starts(ID, "X") ~ str_extract(ID, "^X\\d+"),
      TRUE ~ ID
    )
  ) %>%
  mutate(
    ID_fmt = paste0("sub-", ID_base, "_dol-", Age)
  ) %>%
  mutate(
    ID_extra_character = case_when(
      str_starts(ID, "X") ~ str_replace(ID, "^X\\d+", ""),
      TRUE ~ NA_character_
    )
  ) %>%
  relocate(ID_fmt, ID_base, ID_extra_character, .after = ID)

# Checking on duplicates
covariates %>% group_by(ID_fmt) %>% mutate(n = n()) %>% arrange(desc(n))

# Checking that every folder has its match:
for (ifolder in 1:length(folder_list)){
  folder_name <- folder_list[ifolder]
  folder_path <- file.path(data_path, folder_name)
  # Check if folder or file
  if (!dir.exists(folder_path)){
    next
  }
  # Lookup in covariates table
  match_exists <- folder_name %in% covariates$ID_base
  if (!match_exists){
    stop(paste("Missing match for folder:", folder_name, "; index:", ifolder))
  }
}

# X61 is a good example of a "problem" subject with multiple visits

# Finally, add a column with Rdata_fname
covariates <- covariates %>%
  mutate(RData_file = paste0(ID_fmt, ".RData")) %>%
  dplyr::select(-c(ID_extra_character)) %>%
  relocate(ID_base, ID, RData_file, ID_fmt, Study, source_file) %>%
  rename("ID_raw" = "ID")

write_csv(covariates, file = file.path(outdir, "merged_covariates.csv"))

# MHA Handling ----
# Note - this requires special handling depending on whether there is
# more than one time point per subject
ifolder <- 10
ifolder <- 2379 # testing repeated measurements 5 vs 3?

for (ifolder in 1:length(folder_list)){

  folder_name <- folder_list[ifolder]
  folder_path <- file.path(data_path, folder_name)

  # Check if folder or file
  if (!dir.exists(folder_path)){
    next
  }

  # Lookup in covariates table
  match_exists <- folder_name %in% covariates$ID_base
  if (!match_exists){
    stop(paste("Missing match for folder:", folder_name, "; index:", ifolder))
  }

  # This is potentially multiple rows in the case of repeated visits
  matching_rows <- covariates %>% filter(ID_base == folder_name)

  # List sub-files
  file_list <- list.files(folder_path)

  # Verify that we have one data point per row
  closest_normal_indices  <- which(grepl('Closest', file_list))
  malformation_indices    <- which(grepl('Malformations', file_list))
  spherical_coord_indices <- which(grepl('SphericalCoord', file_list))
  if (nrow(matching_rows) != length(closest_normal_indices)){
    stop(paste("Folder index", ifolder, "has mismatch between number of CN files and matching rows"))
  }
  if (nrow(matching_rows) != length(malformation_indices)){
    stop(paste("Folder index", ifolder, "has mismatch between number of Malform files and matching rows"))
  }
  if (nrow(matching_rows) != length(spherical_coord_indices)){
    stop(paste("Folder index", ifolder, "has mismatch between number of SC files and matching rows"))
  }

  # Loop over matching rows
  for (iVisit in 1:nrow(matching_rows)){

    # Current row we are viewing
    current_row <- matching_rows[iVisit, ]

    # Extract components of "formatted" filename
    FNAME_ID <- current_row$ID_fmt

    # Get the current ID from the covariates file
    image_id <- matching_rows$ID[iVisit]
    image_suffix <- paste0("_", matching_rows$original_file_suffix[iVisit])

    # Find the matching version of each type of file

    ### Closest Normal ----
    closest_normal_index <- which(grepl(image_suffix, file_list[closest_normal_indices]))
    if (length(closest_normal_index) != 1){
      stop(paste0("Issue with folder:", ifolder))
    }
    closest_normal_file  <- file.path(folder_path, file_list[closest_normal_index])
    raw_image            <- SimpleITK::ReadImage(closest_normal_file)
    image_data           <- as.array(raw_image)
    image_data_long <- image_data %>%
      reshape2::melt() %>%
      as_tibble() %>%
      rename(X = Var1,
             Y = Var2,
             Z = Var3,
             value = value)
    closest_normal <- image_data_long

    ### Malformation ----
    malformation_index <- which(grepl(image_suffix, file_list[malformation_indices]))
    if (length(malformation_index) != 1){
      stop(paste0("Issue with folder:", ifolder))
    }
    malformation_file  <- file.path(folder_path, file_list[malformation_index])
    raw_image            <- SimpleITK::ReadImage(malformation_file)
    image_data           <- as.array(raw_image)
    image_data_long <- image_data %>%
      reshape2::melt() %>%
      as_tibble() %>%
      rename(X = Var1,
             Y = Var2,
             value = value)
    malformation_map <- image_data_long

    ### SphericalCoord ----
    spherical_coord_index <- which(grepl(image_suffix, file_list[spherical_coord_indices]))
    if (length(spherical_coord_index) != 1){
      stop(paste0("Issue with folder:", ifolder))
    }
    spherical_coord_file  <- file.path(folder_path, file_list[spherical_coord_index])
    raw_image            <- SimpleITK::ReadImage(spherical_coord_file)
    image_data           <- as.array(raw_image)
    image_data_long <- image_data %>%
      reshape2::melt() %>%
      as_tibble() %>%
      rename(X = Var1,
             Y = Var2,
             Z = Var3,
             value = value)
    spherical_coordinate_map <- image_data_long

    # HERE
    # Final information
    ID = folder_name

    subject_data <- list(
      ID = ID,
      FNAME_ID = FNAME_ID,
      closest_normal = closest_normal,
      malformation_map = malformation_map,
      spherical_coordinate_map = spherical_coordinate_map
    )

    # Save the file
    outdir <- file.path(here("processed/individual_subject_files"))
    output_file <- file.path(outdir, paste0(FNAME_ID, ".RData"))

    save(subject_data, file = output_file)

  } # end of loop over visits within child
} # end of loop over child folders

# Region Information ----

## Region Labels ----
region_label_file  <- file.path(data_path, "RegionLabels.mha")
raw_image          <- SimpleITK::ReadImage(region_label_file)
image_data         <- as.array(raw_image)
image_data_long <- image_data %>%
  reshape2::melt() %>%
  as_tibble() %>%
  rename(X = Var1,
         Y = Var2,
         value = value)
region_labels <- image_data_long

## Region Labels ----
spherical_mask_file  <- file.path(data_path, "SphericalMapMask.mha")
raw_image          <- SimpleITK::ReadImage(spherical_mask_file)
image_data         <- as.array(raw_image)
image_data_long <- image_data %>%
  reshape2::melt() %>%
  as_tibble() %>%
  rename(X = Var1,
         Y = Var2,
         value = value)
spherical_map_mask <- image_data_long

save(region_labels, spherical_map_mask, file = file.path(outdir, "region_labels_and_mask.RData"))



