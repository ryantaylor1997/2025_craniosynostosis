################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Summarize pixel-level regressions
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load all pixel-level models ("pixel_models")
load(file = here("results", "models_pixel_regressions.rda"))

# Load image-level data ("cranio_clean") to help create prediction data
load(file = here("data", "cleaned", "obs_data_clean.rda"))

# Create new data to plot smooth ------------------------------------------

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

# Summarize models --------------------------------------------------------

# Combine pixel-level models and add region column
model_pixel_summ <- pixel_models %>%
  select(row, col, coeffs) %>%
  unnest(coeffs) %>%
  # Remove raw bicoronal interaction term (replaced with Bicoronal combo)
  filter(term != "fused_RCoronal:fused_LCoronal") %>%
  arrange(row, col, -str_detect(term, "fused")) %>%
  mutate(term = fct_inorder(term))

# Extract smoother predictions on existing data
model_pixel_smooth <- pixel_models %>%
  select(row, col, region, model) %>%
  # Add plot of smoothed term
  mutate(smooth_pred = map(
    model,
    ~obtain_predictions(
      model = .x, newdata = newdata_cranio_age))) %>%
  select(-model) %>%
  unnest(smooth_pred)

# Save summaries
save(model_pixel_summ, model_pixel_smooth,
     file = here("data", "intermediate", "models_pixel_summary.rda"))

# Plot pixel-level model results ------------------------------------------

plot_pixel <- ggplot(model_pixel_summ) +
  geom_raster(aes(x = row, y = col, fill = estimate)) +
  facet_wrap(~term) +
  coord_fixed() +
  scale_fill_viridis_c(option = "plasma") +
  theme_void() +
  labs(title = "Pixel-Level Model",
       fill = "Coefficient") +
  theme(legend.position = "bottom")

# Plot 1 age smoothed effect per region
plot_pixel_pred <- ggplot(model_pixel_smooth) +
  geom_line(aes(x = age, y = fit,
                group = interaction(row, col),
                color = factor(region)),
            alpha = 0.05, linewidth = 0.05) +
  labs(x = "Age (Days)", y = "Predicted Growth (No Fusion Female)",
       color = "Region") +
  theme_minimal() +
  scale_color_discrete(palette = "Set3") +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.5))) +
  theme(legend.position = "bottom")
