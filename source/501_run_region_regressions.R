################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Run exploratory region-level growth models
################################################################################

source(here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load full growth data ("cranio_point_clean")
load(file = here("data", "cleaned", "point_data_clean.rda"))

# Region-level preliminary analysis ---------------------------------------

# Run (time-intensive) region-level models --------------------------------

## Run models for average pointwise differences in each region
region_models <- cranio_point_clean %>%
  # Rearrange to incorporate subject data
  unnest(data) %>%
  # Take average in each region for each image
  group_by(fname, age, sex, fusion_type, region) %>%
  summarize(diff_avg = mean(diff)) %>%
  ungroup() %>%
  # Add fusion indicators
  left_join(fusion_dict, by = "fusion_type") %>%
  # Nest region-specific data
  group_by(region) %>%
  nest() %>%
  ungroup() %>%
  # Fit model with age as a spline (cubic regression basis for fast estimation)
  mutate(model = map(data,
                     ~gam(diff_avg ~
                            fused_Sagittal + fused_Metopic +
                            fused_RCoronal * fused_LCoronal +
                            sex + s(age, bs = "cr"),
                          data = .x, method = "REML")),
         # Also try linear age
         model_linear = map(data,
                            ~lm(diff_avg ~
                                  fused_Sagittal + fused_Metopic +
                                  fused_RCoronal * fused_LCoronal +
                                  sex + age,
                                data = .x)))

save(region_models,
     file = here::here("results", "models_region_regressions.rda"))
