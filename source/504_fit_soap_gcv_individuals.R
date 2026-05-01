################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Fit individual soap film model to images using GCV grid search
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load growth data nested by observation ("cranio_obs_points")
load(file = here("data", "intermediate", "point_data_observation.rda"))

# Load soap film object ("cranio_soap")
load(file = here("data", "cleaned", "soap_object.rda"))

# Load observation-level data ("cranio_clean")
load(file = here("data", "cleaned", "obs_data_clean.rda"))

# Calculate lambda-based constants ----------------------------------------

# Define grid for penalty term optimization
lambdas_film <- 10^(-7:-1)
lambdas_bound <- lambdas_film

## Calculate other necessary objects that only depend on X and penalty constant:
# 1. "visor" matrix for each penalty constant (hat without the first X)
# 2. Trace of hat matrix (with X cycled to end)

# Next, boundary and interior
so_constants <- tibble(
  expand.grid("sp_film" = lambdas_film, "sp_bound" = lambdas_bound)
) %>%
  mutate(visor = map2(sp_film, sp_bound,
                      ~tcrossprod(solve(crossprod(cranio_soap$X) +
                                          (.x * cranio_soap$S[[1]] +
                                             .y * cranio_soap$S[[2]])),
                                  cranio_soap$X))) %>%
  mutate(hat_trace = map_dbl(visor, ~sum(diag(.x %*% cranio_soap$X)))) %>%
  mutate(constant_row = 1:n())

# Fit frequentist grid search models on all individuals -------------------

cranio_obs_points %<>%
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
            mutate(fit = map(coeffs, ~ cranio_soap$X %*% .x)) %>%
            # Sum of square residuals
            mutate(ssr = map_dbl(fit, ~sum((d$diff - .x)^2))) %>%
            # GCV Score
            mutate(gcv = (nrow(cranio_soap$X) * ssr) / (nrow(cranio_soap$X) - hat_trace)^2)

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

# Clean dataset of models and export --------------------------------------

# Extract all model components
cranio_obs_points %<>%
  unnest_wider(model) %>%
  unnest(c(coeffs, fit))

# Extract coefficient estimates
cranio_models_beta <- cranio_obs_points %>%
  select(fname, coeffs) %>%
  mutate(coeffs = map(coeffs, as.vector)) %>%
  rename(beta = coeffs)

# Merge in observation-level characteristics; order by cranio_clean
cranio_models_beta <- left_join(cranio_clean,
                                cranio_models_beta,
                                by = "fname") %>%
  left_join(fusion_dict, by = "fusion_type") %>%
  select(-beta, everything(), beta) %>%
  unnest_wider(beta, names_sep = "_", simplify = T) %>%
  mutate(fusion_type = factor(fusion_type,
                              levels = levels(cranio_clean$fusion_type)))

# Save coefficients as csv
write_csv(cranio_models_beta,
          here::here("results", "individual_soap_coeffs.csv"))

# Save estimator requirements as RData
save(cranio_models_beta,
     file = here::here("results", "individual_soap_film_fits.rda"))
