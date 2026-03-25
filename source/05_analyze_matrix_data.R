
# Visualize Regions -------------------------------------------------------

# Create dataset of only region assignment for each pixel
matrix_region <- cranio_matrix %>%
  select(row, col, region) %>%
  distinct()

# Identify regions' centroids for labels
matrix_region_center <- matrix_region %>%
  group_by(region) %>%
  summarize(row_center = mean(row),
            col_center = mean(col)) %>%
  ungroup()

# Plot regions
plot_regions <- ggplot() +
  geom_raster(data = matrix_region,
              aes(x = row, y = col, fill = factor(region))) +
  geom_text(data = matrix_region_center,
            aes(x = row_center, y = col_center, label = factor(region))) +
  coord_fixed() +
  theme_void() +
  scale_fill_brewer(palette = "Set3", guide = "none") +
  labs(title = "Region Map", fill = "Region")

# Summary by age ----------------------------------------------------------

# Create dataset of average cell for each fusion
matrix_fusion <- cranio_matrix %>%
  # Extract only fusion type and outcome from nested data
  hoist(data, "fusion_type", "diff") %>%
  select(-c(data, region)) %>%
  unnest_longer(col = c(fusion_type, diff)) %>%
  # Summarize average by pixel and fusion type
  group_by(row, col, fusion_type) %>%
  summarize(diff_avg = mean(diff)) %>%
  ungroup()

# Plot brain map by fusion type
plot_fusion <- ggplot(matrix_fusion) +
  geom_raster(aes(x = row, y = col, fill = diff_avg)) +
  facet_wrap(~fusion_type) +
  coord_fixed() +
  theme_void() +
  scale_fill_viridis_c(option = "turbo") +
  labs(title = "Fusion Average Growth",
       fill = "Avg. Growth") +
  theme(legend.position = "bottom")


# Summary by age ----------------------------------------------------------

### Investigate age

# Create dataset with fusion and age info
matrix_age_fusion <- cranio_matrix %>%
  hoist(data, "fusion_type", "age", "diff") %>%
  select(-c(data, region)) %>%
  unnest_longer(col = c(fusion_type, age, diff)) %>%
  # Add age categories
  mutate(age_bin = cut(age,
                       breaks = cranio_age_breaks,
                       labels = cranio_age_break_labels,
                       right = FALSE, include.lowest = T)) %>%
  # Summarize average in these bins
  group_by(row, col, fusion_type, age_bin) %>%
  summarize(diff_avg = mean(diff)) %>%
  ungroup()

# Plot brain map by fusion type and age
plot_fusion_age <- ggplot(matrix_age_fusion) +
  geom_raster(aes(x = row, y = col, fill = diff_avg)) +
  facet_grid(fusion_type ~ age_bin,
             labeller = label_wrap_gen(width = 10)) +
  coord_fixed() +
  scale_fill_viridis_c(option = "turbo") +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(title = "Fusion - Age Average Growth",
       fill = "Avg. Growth") +
  theme(legend.position = "bottom")

## Determine sample size for these averages

# Count images per fusion and age bin
size_fusion_age <- cranio_sub %>%
  mutate(age_bin = cut(age,
                       breaks = cranio_age_breaks,
                       labels = cranio_age_break_labels,
                       right = FALSE, include.lowest = T)) %>%
  group_by(fusion_type, age_bin) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = age_bin, values_from = n, id_cols = fusion_type)

# Age effect --------------------------------------------------------------

# Check shape of age effect by region and fusion

# Summarize by age, region, and fusion type; have option to round to nearest month
matrix_age_region <- cranio_matrix %>%
  hoist(data, "fusion_type", "age", "diff") %>%
  select(-data) %>%
  unnest_longer(col = c(fusion_type, age, diff)) %>%
  mutate(age_round_month = round(age / 30)) %>%
  group_by(fusion_type, age, region) %>%
  summarize(diff_avg = mean(diff)) %>%
  ungroup()

plot_age <- ggplot(matrix_age_region) +
  # Plot 1 line over age for each region-fusion combo
  geom_smooth(aes(x = age, y = diff_avg,
                  color = fusion_type, fill = fusion_type),
              alpha = 0.2) +
  facet_wrap(~region) +
  scale_color_discrete(palette = fusion_color_dict$color) +
  scale_fill_discrete(palette = fusion_color_dict$color) +
  labs(title = "Fusion - Region Average Growth by Age",
       x = "Age (Days)",
       y = "Avg. Growth",
       color = "Fusion Type",
       fill = "Fusion Type") +
  theme_minimal() +
  theme(legend.position = "bottom")


# Sex effect --------------------------------------------------------------

# Take average pointwise differences in cell by fusion and sex
matrix_sex_fusion <- cranio_matrix %>%
  hoist(data, "fusion_type", "sex", "diff") %>%
  select(-c(data, region)) %>%
  unnest_longer(col = c(fusion_type, sex, diff)) %>%
  # Summarize average in these bins
  group_by(row, col, fusion_type, sex) %>%
  summarize(diff_avg = mean(diff)) %>%
  ungroup()

# Plot brain map by fusion type and age
plot_fusion_sex <- ggplot(matrix_sex_fusion) +
  geom_raster(aes(x = row, y = col, fill = diff_avg)) +
  facet_grid(fusion_type ~ sex,
             labeller = label_wrap_gen(width = 9)) +
  coord_fixed() +
  scale_fill_viridis_c(option = "turbo") +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(title = "Fusion-Sex Average Growth",
       fill = "Avg. Growth") +
  theme(legend.position = "bottom")
