################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Run 2D GAM models by fusion type
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load file of pixel-level data for all images ("cranio_point_clean")
load(file = here("data", "cleaned", "point_data_clean.rda"))

# Run 2D GAMs -------------------------------------------------------------

# Run spatial tensor product models
fusion_models <- cranio_point_clean %>%
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

save(fusion_models,
     file = here("results", "models_te_fusion.rda"))
