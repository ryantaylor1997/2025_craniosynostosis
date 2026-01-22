# Filter pointwise data to same subset
cranio_matrix_sub <- cranio_matrix %>%
  unnest(data) %>%
  filter(fname %in% cranio_sub$fname) %>%
  group_by(row, col, region) %>%
  nest() %>%
  ungroup()

saveRDS(cranio_matrix_sub,
        file = here::here("analysis", "intermediate",
                          "covars_and_matrices_sub.rds"))
