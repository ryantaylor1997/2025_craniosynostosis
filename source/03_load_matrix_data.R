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

# Load regions and mask
load(paste0(DATA_FOLDER, "region_labels_and_mask.RData"))

cranio_matrix %<>%
  left_join(region_labels,
            by = c("row" = "X", "col" = "Y")) %>%
  filter(value != 0) %>%
  rename(region = value)

# Flip axes to match overhead view
cranio_matrix %<>%
  mutate(row = -row,
         col = -col)

saveRDS(cranio_matrix,
        file = here::here("analysis", "intermediate",
                          "covars_and_matrices.rds"))
