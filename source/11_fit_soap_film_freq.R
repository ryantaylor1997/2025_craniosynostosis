
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


# Create nested matrix data -----------------------------------------------

## Reshape data to have pointwise differences for each image
cranio_models <- cranio_matrix %>%
  unnest(data) %>%
  select(fname, fusion_type, row, col, diff) %>%
  group_by(fname, fusion_type) %>%
  nest() %>%
  ungroup()

save(cranio_models, file = here::here("intermediate", "growth_maps.rda"))


# Create subset for testing -----------------------------------------------

# Test with small subset of data
so_test_df <- cranio_models %>%
  filter(!fname %in% cranio_dup_fnames) %>%
  group_by(fusion_type) %>%
  slice(1) %>%
  ungroup()

# Identify range of outcome values for consistent plotting
so_diff_range <- (
  cranio_models %>%
    select(data) %>%
    mutate(max_diff = map_dbl(data, ~max(.$diff))) %>%
    summarize(max_d = quantile(max_diff, 0.999)) %>%
    pull(max_d)
) * c(-0.1, 1)


# Fit soap film with a loop, without boundary -----------------------------

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
    soap_gcv <- (so_n * sum(soap_sqresid)) / (so_n - trace_H_film)^2

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


### Plot output from loop without boundary

# Plot GCV from test models
film_loop_gcv_plot <- ggplot(film_loop_df) +
  geom_line(aes(x = sp, y = gcv)) +
  facet_wrap(~fusion_type, scales = "free_y") +
  scale_x_continuous(transform = "log10")

# Take model with lowest GCV for each patient
film_loop_optim <- film_loop_df %>%
  arrange(fname, gcv, desc(sp)) %>%
  group_by(fname) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(so_test_df, by = c("fname", "fusion_type"))

# Extract data and fitted values for plotting
film_loop_optim_fit <- film_loop_optim %>%
  select(fname, fusion_type, fit, data) %>%
  unnest(c(fit, data))

# Plot model with lowest GCV for each patient
plot_loop_film <- ggplot(film_loop_optim_fit) +
  geom_raster(aes(x = row, y = col, fill = fit)) +
  scale_fill_viridis_c(option = "turbo", limits = so_diff_range) +
  labs(x = NULL, y = NULL, fill = "Pred. Growth",
       title = "Interior using Loop") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom") +
  facet_wrap(~fusion_type, nrow = 1)

# Plot actual observations
plot_so_test_data <- ggplot(so_test_df %>% unnest(data)) +
  geom_raster(aes(x = row, y = col, fill = diff)) +
  scale_fill_viridis_c(option = "turbo", limits = so_diff_range) +
  labs(x = NULL, y = NULL, fill = "Growth",
       title = "True Data") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom") +
  facet_wrap(~fusion_type, nrow = 1)

# Combine and print
loop_film_both <- ggarrange(plot_loop_film,
                            plot_so_test_data,
                            ncol = 1)


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
  mutate(gcv = (so_n * ssr) / (so_n - hat_trace)^2)

### Plot GCV from test models
film_map_gcv_plot <- ggplot(film_map_df) +
  geom_line(aes(x = sp, y = gcv)) +
  facet_wrap(~fusion_type, scales = "free_y") +
  scale_x_continuous(transform = "log10")

# Extract best model by GCV
film_map_optim <- film_map_df %>%
  arrange(fname, gcv, desc(sp)) %>%
  group_by(fname) %>%
  slice(1) %>%
  ungroup()

# Extract data and fitted values for plotting
film_map_optim_fit <- film_map_optim %>%
  select(fname, fusion_type, fit, data) %>%
  unnest(c(fit, data))

# Plot model with lowest GCV for each patient
plot_map_nobound <- ggplot(film_map_optim_fit) +
  geom_raster(aes(x = row, y = col, fill = fit)) +
  scale_fill_viridis_c(option = "turbo", limits = so_diff_range) +
  labs(x = NULL, y = NULL,
       fill = "Pred. Growth",
       title = "Predictions: Interior Only") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom") +
  facet_wrap(~fusion_type, nrow = 1)

# Combine and print
map_film_both <- ggarrange(plot_map_nobound,
                           plot_so_test_data,
                           ncol = 1)


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
    soap_gcv <- (so_n * sum(soap_sqresid)) / (so_n - trace_H)^2

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

# Plot GCV from test models
soap_loop_gcv_plot <- ggplot(soap_loop_df %>%
                               unnest(all_gcv)) +
  geom_raster(aes(x = sp_bound_all, y = sp_film_all, fill = gcv_all)) +
  facet_wrap(~fusion_type, scales = "free_y") +
  scale_x_continuous(transform = "log10") +
  scale_y_continuous(transform = "log10") +
  scale_fill_continuous(transform = "log")

# Extract data and fitted values for plotting
soap_loop_fit <- soap_loop_df %>%
  select(fname, fusion_type, fit, data) %>%
  unnest(c(fit, data))

# Plot model with lowest GCV for each patient
plot_loop_soap <- ggplot(soap_loop_fit) +
  geom_raster(aes(x = row, y = col, fill = fit)) +
  scale_fill_viridis_c(option = "turbo", limits = so_diff_range) +
  labs(x = NULL, y = NULL,
       fill = "Pred. Growth",
       title = "Predictions Using Loop") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom") +
  facet_wrap(~fusion_type, nrow = 1)

# Combine and print
loop_soap_both <- ggarrange(plot_loop_soap,
                            plot_so_test_data,
                            ncol = 1)


# Fit model in nested data with boundary ----------------------------------

# Define ranges of penalty block matrices
bound_cols <- which(apply(cranio_soap$S[[1]], 2, sum) != 0)
film_cols <- which(apply(cranio_soap$S[[2]], 2, sum) != 0)

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
                                 ~ so_X[, bound_cols] %*% .x[bound_cols]),
                 fit_int = map(coeffs,
                               ~ so_X[, film_cols] %*% .x[film_cols])) %>%
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

### Extract all model components
soap_map_df %<>%
  unnest_wider(model) %>%
  unnest(c(coeffs, fit, fit_bound, fit_int))

# Extract data and fitted values for plotting
so_test_fit <- soap_map_df %>%
  select(fname, fusion_type, fit, fit_bound, fit_int, data) %>%
  unnest(c(fit, fit_bound, fit_int, data)) %>%
  mutate(across(matches("fit"), as.numeric)) %>%
  # Compare separate components
  mutate(fit_sep = fit_bound + fit_int,
         fit_comp = (fit - fit_sep)^2)


### Plot model with lowest GCV for each patient
plot_map_soap <- ggplot(so_test_fit) +
  geom_raster(aes(x = row, y = col, fill = fit)) +
  scale_fill_viridis_c(option = "turbo", limits = so_diff_range) +
  labs(x = NULL, y = NULL,
       fill = "Pred. Growth",
       title = "Sample Predictions") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom") +
  facet_wrap(~fusion_type, nrow = 1)

# Combine and print
map_soap_both <- ggarrange(plot_map_soap,
                           plot_so_test_data,
                           ncol = 1)

### Compare predictions with only bound, only interior, both splines, and truth
plot_map_build <- so_test_fit %>%
  pivot_longer(c(fit_int, fit_bound, fit, diff)) %>%
  mutate(name_cln = fct_recode(fct_inorder(name),
                               "Data" = "diff",
                               "Interior" = "fit_int",
                               "Boundary" = "fit_bound",
                               "Soap Film" = "fit"),
  ) %>%
  ggplot() +
  geom_raster(aes(x = row, y = col, fill = value)) +
  scale_fill_continuous_divergingx(palette = "RdYlBu", rev = TRUE,
                                   l3 = 0, p3 = 2) +
  labs(x = NULL, y = NULL, fill = "Pred. Growth",
       title = "Predictions") +
  facet_grid(name_cln ~ fusion_type) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom")
