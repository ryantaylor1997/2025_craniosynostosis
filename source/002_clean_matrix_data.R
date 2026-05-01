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
load(file = file.path(DATA_FOLDER, "region_labels_and_mask.RData"))

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

cranio_point <- cranio_matrix

rm(cranio_matrix)


# Save complete output ----------------------------------------------------

save(cranio_point,
     file = here("data", "intermediate", "point_data_complete.rda"))

# Remove duplicate data ---------------------------------------------------

# Find instances where average of all observations is the same
dup_avgs <- cranio_point %>%
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
cranio_point_clean <- cranio_point %>%
  unnest(data) %>%
  filter(fname %in% cranio_clean$fname) %>%
  group_by(row, col, region) %>%
  nest() %>%
  ungroup()

# Save final subsets ------------------------------------------------------

save(cranio_point_clean,
     file = here("data", "cleaned", "point_data_clean.rda"))

save(cranio_clean,
     file = here("data", "cleaned", "obs_data_clean.rda"))

# Save other useful matrix data -------------------------------------------

# Isolate coordinates and corresponding region to build smoother
point_coords <- cranio_point_clean %>%
  select(-data)

save(point_coords,
     file = here("data", "intermediate", "point_coordinates.rda"))

# Create spatial data of regions with centroids for plotting region-level work
region_shape <- point_coords %>%
  st_as_sf(coords = c("row", "col")) %>%
  group_by(region) %>%
  summarize(geometry = st_union(geometry)) %>%
  ungroup() %>%
  mutate(geometry = st_convex_hull(geometry)) %>%
  mutate(centroid = st_centroid(geometry))

save(region_shape,
     file = here("data", "intermediate", "region_shape_object.rda"))

## Reshape data to have each image in rows and all pixels nested
cranio_obs_points <- cranio_point_clean %>%
  unnest(data) %>%
  select(fname, fusion_type, row, col, diff) %>%
  group_by(fname, fusion_type) %>%
  nest() %>%
  ungroup()

save(cranio_obs_points,
     file = here("data", "intermediate", "point_data_observation.rda"))

