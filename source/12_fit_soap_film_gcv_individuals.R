
# Fit frequentist grid search models on all individuals -------------------

if(DO_INDIV_SOAPS){

  # Add to existing models data set
  cranio_models %<>%
    # Map over every image
    mutate(
      model = map(
        data,
        function(d){

          # Fit estimates for every combination of penalty constants
          fits <- so_constants %>%
            select(-visor) %>%
            # Estimate coefficients
            mutate(coeffs = map(constant_row,
                                ~so_constants$visor[[.x]] %*% d$diff)) %>%
            # Fitted values
            mutate(fit = map(coeffs, ~ so_X %*% .x)) %>%
            # Sum of square residuals
            mutate(ssr = map_dbl(fit, ~sum((d$diff - .x)^2))) %>%
            # GCV Score
            mutate(gcv = (so_n * ssr) / (so_n - hat_trace)^2)

          # Keep the GCV scores for all combinations of constants
          fits_gcv <- fits %>%
            select(sp_bound, sp_film, gcv) %>%
            rename_with(.cols = everything(), ~paste0(.x, "_all"))

          # Keep the optimal model (lowest GCV -> smoothes interior -> smoothest boundary)
          image_model <- fits %>%
            arrange(gcv, desc(sp_film), desc(sp_bound)) %>%
            slice(1) %>%
            select(-c(hat_trace, constant_row)) %>%
            # Append all GCV scores
            mutate(all_scores = list(fits_gcv)) %>%
            tibble_row()

          return(image_model)
        }
      )
    )

  saveRDS(cranio_models,
          file = here::here("intermediate", "individual_soap_film_df.rds"))

} else {
  cranio_models <- readRDS(
    file = here::here("intermediate", "individual_soap_film_df.rds"))
}


# Clean dataset of models and export --------------------------------------


# Extract all model components
cranio_models %<>%
  unnest_wider(model) %>%
  unnest(c(coeffs, fit))

# Extract coefficient estimates
cranio_models_beta <- cranio_models %>%
  select(fname, coeffs) %>%
  mutate(coeffs = map(coeffs, as.vector)) %>%
  rename(beta = coeffs) %>%
  left_join(cranio_sub, by = "fname") %>%
  left_join(fusion_dict, by = "fusion_type") %>%
  select(-beta, everything(), beta) %>%
  unnest_wider(beta, names_sep = "_", simplify = T)

# Save coefficients as csv
write_csv(cranio_models_beta,
          here::here("intermediate", "individual_soap_coeffs.csv"))

# Save estimator requirements as RData
save(cranio_models_beta, so_X,
     file = here::here("intermediate", "individual_soap_film_fits.rda"))
