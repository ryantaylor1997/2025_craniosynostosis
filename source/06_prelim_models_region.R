
# Region-level preliminary analysis ---------------------------------------

# Recreate pixel map with region labels
matrix_region <- cranio_matrix %>%
  select(row, col, region) %>%
  distinct() %>%
  st_as_sf(coords = c("row", "col")) %>%
  group_by(region) %>%
  summarize(geometry = st_union(geometry)) %>%
  ungroup() %>%
  mutate(geometry = st_convex_hull(geometry)) %>%
  mutate(centroid = st_centroid(geometry))

save(matrix_region, file = here::here("intermediate", "region_shape_object.rda"))


# Run (time-intensive) region-level models --------------------------------

if(DO_PRELIM_MODEL_REGION){

  ## Run models for average pointwise differences in each region
  region_models <- cranio_matrix %>%
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

  saveRDS(
    region_models,
    file = here::here("intermediate",
                      "prelim_models_region.rds"))

} else {
  region_models <- readRDS(
    file = here::here("intermediate",
                      "prelim_models_region.rds")) }


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

### Merge model summaries with matrix for plotting
model_region_matrix <- region_models %>%
  select(region, coeffs_combo_cln) %>%
  unnest(coeffs_combo_cln) %>%
  full_join(matrix_region, by = "region")


# Plot Region model output ------------------------------------------------

# Plot parametric coefficients
plot_region_smooth <- ggplot(model_region_matrix %>% filter(term != "age")) +
  geom_sf(aes(geometry = geometry, fill = estimate_smooth), linewidth = 0) +
  geom_sf_text(aes(geometry = centroid, label = region)) +
  facet_wrap(~term) +
  scale_fill_viridis_c(option = "plasma") +
  theme_void() +
  labs(title = "Region-Level Smooth Age Model",
       fill = "Coefficient") +
  theme(legend.position = "bottom")

# Plot parameters from alternate version of the model
plot_region_linear <- ggplot(model_region_matrix) +
  geom_sf(aes(geometry = geometry, fill = estimate_linear), linewidth = 0) +
  geom_sf_text(aes(geometry = centroid, label = region)) +
  facet_wrap(~term) +
  scale_fill_viridis_c(option = "plasma") +
  theme_void() +
  labs(title = "Region-Level Linear Age Model",
       fill = "Coefficient") +
  theme(legend.position = "bottom")

### Plot smoothed effect

# Define new data for predictions
newdata_cranio_age <- tibble(
  fused_Sagittal = 0,
  fused_Metopic = 0,
  fused_RCoronal = 0,
  fused_LCoronal = 0,
  sex = 0,
  age = seq(min(cranio_sub$age),
            max(cranio_sub$age),
            length.out = 30)
)

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
