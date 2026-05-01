################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Summarize 2D GAM models by fusion type
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load fusion-specific 2-D GAM models ("fusion_models")
load(file = here("results", "models_te_fusion.rda"))

# Load shape of surface ("cranio_bound_shape")
load(file = here("data", "intermediate", "boundary_polygon.rda"))

# Load boundary points ("cranio_bound_points") [comes with "bound_list"]
load(file = here("data", "intermediate", "boundary_objects.rda"))

# Load image-level data ("cranio_clean") to help create prediction data
load(file = here("data", "cleaned", "obs_data_clean.rda"))

# Load color mapping for fusion types ("fusion_color_dict")
load(file = here::here("data", "intermediate", "fusion_color_mapping.rda"))

# Create new GAM data -----------------------------------------------------

# Define new data for predictions
newdata_cranio_age <- tibble(
  fused_Sagittal = 0,
  fused_Metopic = 0,
  fused_RCoronal = 0,
  fused_LCoronal = 0,
  sex = 0,
  age = seq(min(cranio_clean$age),
            max(cranio_clean$age),
            length.out = 30)
)

# Set new data for GAM
newdata_fusiongam <- tibble(
  fused_Sagittal = 0,
  fused_Metopic = 0,
  fused_RCoronal = 0,
  fused_LCoronal = 0,
  sex = 0,
  age = 300,
  row = list(seq(from = (min(cranio_bound_points$row) - 10),
                 to = (max(cranio_bound_points$row) + 10),
                 length.out = 50)),
  col = list(seq(from = (min(cranio_bound_points$col) - 10),
                 to = (max(cranio_bound_points$col) + 10),
                 length.out = 50))
) %>%
  unnest_longer(row) %>% unnest_longer(col)

# Add row and column to age smooth new data (to be cancelled out)
newdata_age_gam <- newdata_cranio_age %>%
  mutate(row = 0, col = 0)

# Predict on new data -----------------------------------------------------

# Predict from this model for spatial pattern
fusion_models  %<>%
  mutate(pred_data = map(
    model,
    ~obtain_predictions(.x, newdata = newdata_fusiongam)))

# Flatten predictions into their own data
fusion_model_summ <- fusion_models %>%
  select(fusion_type, pred_data) %>%
  unnest(pred_data)

# Plot GAM predictions ----------------------------------------------------

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


