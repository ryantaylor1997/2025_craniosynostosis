
# Check duplicates --------------------------------------------------------

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

saveRDS(cranio_dup_fnames, here::here("analysis", "intermediate",
                                      "duplicate_fnames.rds"))

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
