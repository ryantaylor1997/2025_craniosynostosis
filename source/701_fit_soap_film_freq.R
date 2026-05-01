################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Fit test soap film models using frequentist methods
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load growth data nested by observation ("cranio_obs_points")
load(file = here("data", "intermediate", "point_data_observation.rda"))

# Load design matrix from soap film object ("so_X")
load(file = here("data", "intermediate", "soap_design_full.rda"))

# Load soap design matrix for only interior basis fns ("so_X_film")
load(file = here("data", "intermediate", "soap_design_interior.rda"))

# Load column indices identifying interior vs boundary ("cols_bound", "cols_film")
load(file = here("data", "intermediate", "soap_basis_indices.rda"))

# Load penalty matrix for only the interior ("S_film")
load(file = here("data", "intermediate", "soap_penalty_interior.rda"))

# Load separate penalty matrices for boundary and interior ("S_bound0", "S_film0")
load(file = here("data", "intermediate", "soap_penalty_separates.rda"))

# Calculate lambda-based constants ----------------------------------------

# Define grid for penalty term optimization
lambdas_film <- 10^(-7:-1)
lambdas_bound <- lambdas_film

## Calculate other necessary objects that only depend on X and penalty constant:
# 1. "visor" matrix for each penalty constant (hat without the first X)
# 2. Trace of hat matrix (with X cycled to end)

# First, only interior film
so_constants_film <- tibble("sp" = lambdas_film) %>%
  mutate(visor = map(sp,
                     ~tcrossprod(solve(crossprod(so_X_film) +
                                         .x * S_film),
                                 so_X_film))) %>%
  mutate(hat_trace = map_dbl(visor, ~sum(diag(.x %*% so_X_film))))

# Next, boundary and interior
so_constants <- tibble(
  expand.grid("sp_film" = lambdas_film, "sp_bound" = lambdas_bound)
) %>%
  mutate(visor = map2(sp_film, sp_bound,
                      ~tcrossprod(solve(crossprod(so_X) +
                                          (.x * S_bound0 + .y * S_film0)),
                                  so_X))) %>%
  mutate(hat_trace = map_dbl(visor, ~sum(diag(.x %*% so_X)))) %>%
  mutate(constant_row = 1:n())

# Create subset for testing -----------------------------------------------

# Test with small subset of data
so_test_df <- cranio_obs_points %>%
  group_by(fusion_type) %>%
  slice(1) %>%
  ungroup()

save(so_test_df,
     file = here("data", "intermediate", "soap_test_sample_data.rda"))

# Set range of colors for consistency
so_diff_range <- (
  so_test_df %>%
    select(data) %>%
    mutate(max_diff = map_dbl(data, ~max(.$diff))) %>%
    summarize(max_d = quantile(max_diff, 0.999)) %>%
    pull(max_d)
) * c(-0.1, 1)

save(so_diff_range,
     file = here("data", "intermediate", "soap_test_sample_range.rda"))

# Fit soap film with a for loop, without boundary -------------------------

# Run loop over individuals and potential penalty constants
film_loop_df <- NULL

for(i in 1:nrow(so_test_df)){
  for(j in 1:nrow(so_constants_film)){

    # Data to fit model on
    d <- so_test_df$data[[i]]

    # File name to identify output
    f <- so_test_df$fname[i]

    # Matrix helpful for influence and coefficients
    reg_mx_film <- so_constants_film$visor[[j]]

    trace_H_film <- so_constants_film$hat_trace[[j]]

    # Outcome vector
    soap_y <- d$diff

    # Coefficient estimates
    soap_b_film <- reg_mx_film %*% soap_y

    # Fitted values and square residuals
    soap_fit <- so_X_film %*% soap_b_film
    soap_sqresid <- (soap_y - soap_fit)^2

    # GCV score
    soap_gcv <- (nrow(so_X) * sum(soap_sqresid)) / (nrow(so_X) - trace_H_film)^2

    out_dl <- tibble_row(
      "fname" = f,
      "sp" = so_constants_film$sp[[j]],
      "beta" = list(soap_b_film),
      "fit" = list(soap_fit),
      "gcv" = soap_gcv)

    film_loop_df %<>% bind_rows(out_dl)
  }}

# Add fusion type to this output
film_loop_df %<>%
  left_join(so_test_df %>% select(-data))

save(film_loop_df,
     file = here("data", "intermediate", "soap_test_interior_loop.rda"))

# Fit model without boundary in nested data -------------------------------

# Add model fits to data frame
film_map_df <- so_test_df %>%
  # Expand data to one row per patient and penalty constant
  mutate(constants = list(so_constants_film)) %>%
  unnest(constants) %>%
  # Estimate coefficients
  mutate(coeffs = map2(visor, data,
                       ~.x %*% .y$diff)) %>%
  # Fitted values
  mutate(fit = map(coeffs, ~ so_X_film %*% .x)) %>%
  # Square residuals
  mutate(ssr = map2_dbl(data, fit,
                        ~ sum((.x$diff - .y)^2))) %>%
  # GCV Score
  mutate(gcv = (nrow(so_X) * ssr) / (nrow(so_X) - hat_trace)^2)

save(film_map_df,
     file = here("data", "intermediate", "soap_test_interior_map.rda"))

# Fit soap film with boundary in loop -------------------------------------

# Run loop over individuals and potential penalty constants
soap_loop_df <- NULL

for(i in 1:nrow(so_test_df)){

  # Data to fit model on
  d <- so_test_df$data[[i]]

  # File name to identify output
  f <- so_test_df$fname[i]

  # Outcome vector
  soap_y <- d$diff

  df_i <- NULL

  for(j in 1:nrow(so_constants)){

    # Matrix helpful for influence and coefficients
    reg_mx <- so_constants$visor[[j]]

    trace_H <- so_constants$hat_trace[j]

    # Coefficient estimates
    soap_b <- reg_mx %*% soap_y

    # Fitted values and square residuals
    soap_fit <- so_X %*% soap_b
    soap_sqresid <- (soap_y - soap_fit)^2

    # GCV score
    soap_gcv <- (nrow(so_X) * sum(soap_sqresid)) / (nrow(so_X) - trace_H)^2

    # Add to image dataset
    df_row_j <- tibble_row(
      "sp_bound" = so_constants$sp_bound[[j]],
      "sp_film" = so_constants$sp_film[[j]],
      "beta" = list(soap_b),
      "fit" = list(soap_fit),
      "gcv" = soap_gcv)

    df_i %<>% bind_rows(df_row_j)
  }

  # Keep all GCV values
  df_i_gcv <- df_i %>%
    select(-c(beta, fit)) %>%
    rename_with(.cols = everything(), ~paste0(., "_all"))

  # Sort so best GCV / most smooth is on top
  df_i %<>%
    arrange(gcv, desc(sp_film), desc(sp_bound))

  # Extract best model; add all gcv scores
  out_df <- df_i %>%
    slice(1) %>%
    mutate("fname" = f, .before = everything()) %>%
    mutate(all_gcv = list(df_i_gcv))

  # Append to final data
  soap_loop_df %<>% bind_rows(out_df)
}

# Add patient info to output
soap_loop_df %<>% left_join(so_test_df, by = c("fname"))

save(soap_loop_df,
     file = here("data", "intermediate", "soap_test_full_loop.rda"))

# Fit model in nested data with boundary ----------------------------------

# Add model fits to data frame
soap_map_df <- so_test_df %>%
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
          mutate(fit = map(coeffs, ~ so_X %*% .x),
                 fit_bound = map(coeffs,
                                 ~ so_X[, cols_bound] %*% .x[cols_bound]),
                 fit_int = map(coeffs,
                               ~ so_X[, cols_film] %*% .x[cols_film])) %>%
          # Sum of square residuals
          mutate(ssr = map_dbl(fit, ~sum((d$diff - .x)^2))) %>%
          # GCV Score
          mutate(gcv = (nrow(so_X) * ssr) / (nrow(so_X) - hat_trace)^2)

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

### Extract all model components
soap_map_df %<>%
  unnest_wider(model) %>%
  unnest(c(coeffs, fit, fit_bound, fit_int))

save(soap_map_df,
     file = here("data", "intermediate", "soap_test_full_map.rda"))
