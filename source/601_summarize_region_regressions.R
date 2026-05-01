################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Summarize region-level regression models
################################################################################

source(here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load regression models in nested data set ("region_models")
load(file = here::here("results", "models_region_regressions.rda"))

# Load region spatial data ("region_shape")
load(file = here("data", "intermediate", "region_shape_object.rda"))

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

# Clean output from these region-level models -----------------------------

region_models %<>%
  ## Add bicoronal variance to final coefficients; output coefficient summary
  mutate(across(matches("model"),
                ~map(.,
                     function(m){

                       # Make constrast vector identifying all coronal variables
                       is_coronal <- as.numeric(str_detect(names(coef(m)), "Coronal"))

                       # Compute variance of bicoronal combination
                       bicoronal_v <- crossprod(is_coronal,  vcov(m)) %*% is_coronal
                       bicoronal_e <- crossprod(is_coronal, coef(m))

                       # Add bicoronal as row in tidy summary
                       bicoronal_row <- tibble(
                         term = "fused_Bicoronal",
                         estimate = bicoronal_e[1,1],
                         std.error = sqrt(bicoronal_v[1,1])) %>%
                         mutate(p.value = 2 * pnorm(abs(estimate / std.error),
                                                    lower.tail = FALSE))

                       # Output coefficient table (smaller than model)
                       return(tidy(m, parametric = TRUE) %>%
                                bind_rows(bicoronal_row))
                     }),
                .names = "{.col}_coeffs")) %>%
  # Rename individual model terms before merging
  mutate(model_coeffs = map(model_coeffs,
                            ~rename_with(.,
                                         .cols = -term,
                                         .fn = ~paste0(., "_smooth"))),
         model_linear_coeffs = map(model_linear_coeffs,
                                   ~rename_with(.,
                                                .cols = -term,
                                                .fn = ~paste0(., "_linear")))) %>%
  # Merge model summaries
  mutate(coeffs_combo = map2(model_coeffs, model_linear_coeffs,
                             ~full_join(.x, .y, by = "term"))) %>%
  # Clean model summary
  mutate(coeffs_combo_cln = map(coeffs_combo,
                                ~filter(., term != "fused_RCoronal:fused_LCoronal") %>%
                                  arrange(-str_detect(term, "fused")) %>%
                                  mutate(term = fct_inorder(term))))

### Merge model summaries with points on a shape for plotting
region_models_shape <- region_models %>%
  select(region, coeffs_combo_cln) %>%
  unnest(coeffs_combo_cln) %>%
  full_join(region_shape, by = "region")

# Plot Region model output ------------------------------------------------

# Plot parametric coefficients
plot_region_smooth <- ggplot(region_models_shape %>% filter(term != "age")) +
  geom_sf(aes(geometry = geometry, fill = estimate_smooth), linewidth = 0) +
  geom_sf_text(aes(geometry = centroid, label = region)) +
  facet_wrap(~term) +
  scale_fill_viridis_c(option = "plasma") +
  theme_void() +
  labs(title = "Region-Level Smooth Age Model",
       fill = "Coefficient") +
  theme(legend.position = "bottom")

# Plot parameters from alternate version of the model
plot_region_linear <- ggplot(region_models_shape) +
  geom_sf(aes(geometry = geometry, fill = estimate_linear), linewidth = 0) +
  geom_sf_text(aes(geometry = centroid, label = region)) +
  facet_wrap(~term) +
  scale_fill_viridis_c(option = "plasma") +
  theme_void() +
  labs(title = "Region-Level Linear Age Model",
       fill = "Coefficient") +
  theme(legend.position = "bottom")

# Plot age smooth ---------------------------------------------------------

### Plot smoothed effect

# Calculate predicted values and Conf Ints
model_region_smooth <- region_models %>%
  select(region, model) %>%
  # Add predictions from smooth term
  mutate(smooth_pred = map(
    model,
    ~obtain_predictions(
      model = .x, newdata = newdata_cranio_age))) %>%
  select(-model) %>%
  unnest(smooth_pred)

# Add linear age effect to compare
model_region_linear <- region_models %>%
  select(region, model_linear) %>%
  mutate(pred_data = list(newdata_cranio_age),
         lin_pred = map(
           model_linear,
           ~predict(.x, newdata = newdata_cranio_age, se.fit = TRUE) %>%
             as.data.frame())) %>%
  select(-model_linear) %>%
  unnest(c(pred_data, lin_pred)) %>%
  mutate(ci_ll = fit - 1.96 * se.fit,
         ci_ul = fit + 1.96 * se.fit)

# Combine smooth and linear predictions
model_region_combo <- model_region_smooth %>%
  select(region, age, fit, ci_ll, ci_ul) %>%
  mutate(type = "Smooth") %>%
  bind_rows(
    model_region_linear %>%
      select(region, age, fit, ci_ll, ci_ul) %>%
      mutate(type = "Linear")
  )

# Plot 1 smoothed effect per region
plot_region_smooth_pred <- ggplot(model_region_combo) +
  geom_line(aes(x = age, y = fit,
                color = factor(region), linetype = rev(type))) +
  geom_ribbon(aes(x = age,
                  ymin = ci_ll, ymax = ci_ul,
                  fill = factor(region),
                  linetype = type),
              alpha = 0.3) +
  facet_wrap(~type, nrow = 1) +
  scale_linetype(guide = "none") +
  scale_color_discrete(palette = "Set3") +
  labs(x = "Age (Days)", y = "Predicted Growth (No Fusion Female)",
       color = "Region", fill = "Region") +
  theme_minimal() +
  theme(legend.position = "bottom")

