################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Clean growth matrix data for all images in cleaned image-level file
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load image-level complete file ("cranio")
load(file = here("data", "intermediate", "obs_data_complete.rda"))

# Load image-level subset file ("cranio_sub")
load(file = here("data", "intermediate", "obs_data_subset.rda"))

# Load regions and mask ("region_labels")
load(file = paste0(DATA_FOLDER, "region_labels_and_mask.RData"))


# Load growth matrix files ------------------------------------------------

# Create minimal matrix data
cranio_matrix <- cranio %>%
  select(fname, age, sex, fusion_type) %>%
  rowwise() %>%
  mutate(
    raw_data = load_matrix_data(fname, folder = MATRIX_FOLDER)
  ) %>%
  ungroup()

# Rearrange so each row is a pixel
cranio_matrix %<>%
  unnest(raw_data) %>%
  group_by(row, col) %>%
  nest() %>%
  ungroup()

cranio_matrix %<>%
  left_join(region_labels,
            by = c("row" = "X", "col" = "Y")) %>%
  filter(value != 0) %>%
  rename(region = value)

# Flip axes to match overhead view
cranio_matrix %<>%
  mutate(row = -row,
         col = -col)

save(cranio_matrix,
     file = here("data", "intermediate", "point_data_complete.rda"))

# Isolate only the coordinates we need to build smoother
matrix_coords <- cranio_matrix %>%
  select(row, col)

save(matrix_coords,
     file = here("data", "intermediate", "point_coordinates.rda"))

# Remove duplicate data ---------------------------------------------------

# Find instances where average of all observations is the same
dup_avgs <- cranio_matrix %>%
  unnest(data) %>%
  select(fname, fusion_type, age, row, col, diff) %>%
  group_by(fname, fusion_type, age) %>%
  nest() %>%
  mutate(diff_avg = map_dbl(data, ~mean(.$diff))) %>%
  select(-data) %>%
  group_by(diff_avg) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > 1)

cranio_dup_fnames <- unique(dup_avgs$fname)

save(cranio_dup_fnames,
     file = here("data", "intermediate", "duplicate_fnames.rda"))

# Reshape these to 1 per set of duplicates
dups_collected <- dup_avgs %>%
  group_by(diff_avg, n) %>%
  summarize(fnames = paste(na.omit(unique(fname)), collapse = "/"),
            fusions = paste(na.omit(unique(fusion_type)), collapse = "/"),
            ages = paste(na.omit(unique(age)), collapse = "/")) %>%
  ungroup() %>%
  arrange(desc(n), diff_avg)

# Add information about these duplicate pairs
dups_collected %<>%
  # Check study
  mutate(is_cnh = if_else(str_detect(fnames, "CNH"), 1, 0),
         is_chco = if_else(str_detect(fnames, "CHCO"), 1, 0))

dup_by_study <- dups_collected %>%
  group_by(is_cnh, is_chco) %>%
  summarize(n = n()) %>%
  ungroup()

# Filter data to final subset ---------------------------------------------

# Remove duplicates from subject-level data
cranio_clean <- cranio_sub %>% filter(!fname %in% cranio_dup_fnames)

# Filter pointwise data to same subset
cranio_matrix_clean <- cranio_matrix %>%
  unnest(data) %>%
  filter(fname %in% cranio_clean$fname) %>%
  group_by(row, col, region) %>%
  nest() %>%
  ungroup()

# Save final subsets
save(cranio_matrix_clean,
     file = here("data", "cleaned", "point_data_clean.rda"))

save(cranio_clean,
     file = here("data", "cleaned", "obs_data_clean.rda"))

