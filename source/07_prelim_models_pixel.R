
# Run Pixel-level preliminary models --------------------------------------


if(DO_PRELIM_MODEL_PIXEL){

  ## Run models for pointwise differences in each pixel
  pixel_models <- cranio_matrix %>%
    mutate(
      model = map(data,
                  function(d){
                    d_reg <- d %>% left_join(fusion_dict, by = "fusion_type")

                    mdl <- bam(diff ~
                                 fused_Sagittal + fused_Metopic +
                                 fused_RCoronal * fused_LCoronal +
                                 sex + s(age, bs = "cr"),
                               data = d_reg)

                    return(mdl)
                  })) %>%
    mutate(coeffs = map(
      model,
      function(m){

        # Make constrast vector identifying all coronal variables
        is_coronal <- as.numeric(str_detect(names(coef(m)), "Coronal"))

        # Compute variance of bicoronal combination
        bicoronal_v <- crossprod(is_coronal, vcov(m)) %*% is_coronal
        bicoronal_e <- crossprod(is_coronal, coef(m))

        # Add bicoronal as row in tidy summary
        bicoronal_row <- tibble(
          term = "fused_Bicoronal",
          estimate = bicoronal_e[1,1],
          std.error = sqrt(bicoronal_v[1,1])) %>%
          mutate(statistic = estimate / std.error,
                 p.value = 2 * pnorm(abs(estimate / std.error),
                                     lower.tail = FALSE))

        # Output coefficient table (smaller than model)
        return(tidy(m, parametric = TRUE) %>%
                 bind_rows(bicoronal_row))
      }))

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

  save(model_pixel_summ, model_pixel_smooth,
       file = here::here("intermediate",
                         "prelim_model_pixels_summ.rda"))
} else{ load(file = here::here("intermediate",
                               "prelim_model_pixels_summ.rda")) }

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

