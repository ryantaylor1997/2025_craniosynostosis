
# Fusion-level 2D GAMs ----------------------------------------------------

if(DO_PRELIM_MODEL_FUSION){

  # Run preliminary spatial models
  fusion_models <- cranio_matrix %>%
    unnest(data) %>%
    # Split data by fusion type
    group_by(fusion_type) %>%
    nest() %>%
    ungroup() %>%
    # Fit model with linear sex effect, smoothed age term, and row-col tensor surface
    mutate(model = map(data,
                       ~bam(diff ~ sex + s(age, bs = "cr") +
                              te(row, col, k=c(8,8)),
                            data= .x)))

  saveRDS(
    fusion_models,
    file = here::here("analysis", "intermediate", "prelim_models_fusion_te.rds"))
} else {
  fusion_models <- readRDS(
    file = here::here("analysis", "intermediate", "prelim_models_fusion_te.rds")) }

# Create boundary shape around pixels -------------------------------------

# Create polygon around pixels in cranium shape
cranio_bound_shape <- cranio_matrix %>%
  distinct(row, col) %>%
  st_as_sf(coords = c("row", "col")) %>%
  concaveman(concavity = 1)

# Create boundary as points
cranio_bound_matrix <- cranio_bound_shape %>%
  # Add buffer so boundary is beyond data
  st_buffer(1) %>%
  # Convert to line
  st_boundary() %>%
  # Take points along this line at density = 1 unit
  st_line_sample(density = 1) %>%
  st_coordinates() %>%
  as_tibble() %>%
  rename(row = X, col = Y) %>%
  # Add order of appearance in data
  mutate(appearance_order = 1:n())

# Predict from GAM --------------------------------------------------------

# Set new data for GAM
newdata_fusiongam <- tibble(
  fused_Sagittal = 0,
  fused_Metopic = 0,
  fused_RCoronal = 0,
  fused_LCoronal = 0,
  sex = 0,
  age = 300,
  row = list(seq(from = (min(cranio_bound_matrix$row) - 10),
                 to = (max(cranio_bound_matrix$row) + 10),
                 length.out = 50)),
  col = list(seq(from = (min(cranio_bound_matrix$col) - 10),
                 to = (max(cranio_bound_matrix$col) + 10),
                 length.out = 50))
) %>%
  unnest_longer(row) %>% unnest_longer(col)

fusion_models  %<>%
  # Predict from this model for spatial pattern (1 year old male)
  mutate(pred_data = map(
    model,
    ~obtain_predictions(.x, newdata = newdata_fusiongam)))

# Flatten predictions into their own data
fusion_model_summ <- fusion_models %>%
  select(fusion_type, pred_data) %>%
  unnest(pred_data)

# Plot GAM predictions ----------------------------------------------------

# Determine true data scale range
matrix_diff_range <- cranio_matrix %>%
  unnest(data) %>%
  summarize(low_pct = quantile(diff, 0.1, na.rm = T),
            max = max(diff, na.rm = T))

# Plot these predictions (in a rectangle shape with brain mask)
plot_model_fusion <- ggplot() +
  geom_raster(data = fusion_model_summ,
              aes(x = row, y = col, fill = fit)) +
  geom_sf(data = cranio_bound_shape, fill = NA, linewidth = 1) +
  facet_wrap(~fusion_type, nrow = 2) +
  scale_fill_viridis_c(option = "turbo",
                       limits = quantile(fusion_model_summ$fit,
                                         c(0.01, 0.99))) +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(fill = "Pred. Growth\n(10 mo. F)") +
  theme(legend.position = "bottom")

ggsave(here::here("results", "tensor_smooth_fusion.png"),
       plot_model_fusion,
       height = 4, width = 8, units = "in")

# Predict smoothed age effect ---------------------------------------------

# Add row and column to age smooth new-data (to be cancelled out)
newdata_age_gam <- newdata_cranio_age %>%
  mutate(row = 0, col = 0)

## Extract age smoothed predictions
model_fusion_smooth <- fusion_models %>%
  select(fusion_type, model) %>%
  # Add plot of smoothed term
  mutate(smooth_pred = map(
    model,
    ~obtain_predictions(.x,
                        newdata = newdata_age_gam,
                        rm_effects = c("row", "col")))) %>%
  select(-model) %>%
  unnest(smooth_pred)

# Plot 1 smoothed effect per fusion type
plot_fusion_pred <- ggplot(model_fusion_smooth) +
  geom_line(aes(x = age, y = fit, color = fusion_type)) +
  geom_ribbon(aes(x = age, ymin = ci_ll, ymax = ci_ul,
                  group = fusion_type, fill = fusion_type),
              alpha = 0.5) +
  labs(x = "Age (Days)", y = "Predicted Growth (Male)",
       color = "Fusion Type", fill = "Fusion Type") +
  scale_fill_discrete(palette = fusion_color_dict$color) +
  theme_minimal() +
  theme(legend.position = "bottom")


# Print parametric effects ------------------------------------------------

## Compile parametric coefficients (Intercept and sex)
model_fusion_params <- fusion_models %>%
  select(fusion_type, model) %>%
  mutate(params = map(model, ~tidy(., parametric = TRUE))) %>%
  select(-model) %>%
  unnest(params)
