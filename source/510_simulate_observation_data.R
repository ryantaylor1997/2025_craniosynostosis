################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Run simulation to test model fit
################################################################################

source(here::here("source", "000_definitions.R"))

# Simulate observation-level data -----------------------------------------

set.seed(978)

# Generate simulated data of observations
obs_df_sim <- make_sim_data_obs_level(
  n_obs = 80,
  age_max = cranio_max_age,
  sex_prop_1 = 0.5
)

save(obs_df_sim,
     file = here("data", "simulations", "sim_obs_data.rda"))

# Convert observation-level data to design matrix -------------------------

# Specify smooth object (thin plate with shrinkage)
obs_smooth_obj <- s(age, bs = age_bs, k = cranio_knots_age_sim)

# Use function to reference-code smooth by fusion type
obs_smooth_sim_fusion <- construct_reference_smooth(
  sm = obs_smooth_obj, dat = obs_df_sim,
  by_var = "fusion_type", param_formula = ~sex + fusion_type
)

# Determine which penalty columns are non-zero
obs_params_cols <- map(obs_smooth_sim_fusion$S,
                       ~which(
                         apply(., 2,
                               function(x){ sum(x != 0)}) != 0))

save(obs_params_cols,
     file = here("data", "simulations", "sim_obs_smooth_column_list.rda"))

# Scale penalty matrices --------------------------------------------------

# Define scaling term to make variances comparable
eigen_obs_pen_sim <- eigen(
  obs_smooth_sim_fusion$S[[2]][obs_params_cols[[2]], obs_params_cols[[2]]],
  symmetric = T, only.values = T
  )$values

# Set so that trace of inverse is equal to trace of identity mx
scale_obs_pen <- sum(
  1 / eigen_obs_pen_sim[which(abs(eigen_obs_pen_sim) > 1e-10)]
) /
  sum(eigen_obs_pen_sim != 0)

# Multiply penalized terms by this scale factor
obs_smooth_sim_fusion$S[-1] <- map(obs_smooth_sim_fusion$S[-1],
                                   ~.x * scale_obs_pen)

# Set up design matrix for our estimation ---------------------------------

# Combine suture fusion penalties so we estimate the same lambda
obs_pen_list_sim <- list(obs_smooth_sim_fusion$S[[1]],
                         Reduce("+", obs_smooth_sim_fusion$S[-1]))

# Collect design and penalty back into one list
obs_smooth_list_sim <- obs_smooth_sim_fusion
obs_smooth_list_sim[["S"]] <- obs_pen_list_sim

# Save results
save(obs_smooth_list_sim,
     file = here("data", "simulations", "sim_obs_smooth_list.rda"))


# Generate corresponding new data -----------------------------------------

# Generate dense grid over same values that we can compare estimates to
obs_newdata <- obs_df_sim %>%
  expand(fusion_type, sex) %>%
  filter(sex == 0) %>%
  rowwise() %>%
  mutate(age = list(seq(0, cranio_max_age, 10))) %>%
  ungroup() %>%
  unnest(age)

### Generate design matrix for these data based on smooth above
# Include rows that will show deviations from the normative curve

# Set number of non-penalized terms
Q_theta <- length(obs_params_cols[[1]])

# Create predictions by fusion type
obs_new_X_list <- map(obs_smooth_list_sim$smooth,
                      ~PredictMat(.x, data = obs_newdata))

# Combine to create Normative and add unpenalized terms
obs_new_design <- Reduce(
  "cbind", c(
    list(
      model.matrix(~sex + fusion_type, data = obs_newdata),
      Reduce("+", obs_new_X_list)),
    obs_new_X_list[-1]
  ))

# Remove row numbers so they don't get mixed up in the next step
dimnames(obs_new_design) <- NULL

## Add rows for effect of fusion deviation curves

# Identify which rows are not Normative and separate out
obs_newdata_dev <- obs_newdata %>%
  mutate(row_num = row_number()) %>%
  filter(fusion_type != "Normative") %>%
  mutate(fusion_type = paste0(fusion_type, "_Only"))

# Zero out non-smooth terms for these new rows
obs_new_design_dev <- cbind(
  matrix(0,
         nrow(obs_newdata_dev),
         (Q_theta + cranio_knots_age_sim - 1)
  ),
  obs_new_design[
    obs_newdata_dev$row_num,
    (Q_theta + cranio_knots_age_sim):ncol(obs_new_design)]
)

# Add these "deviation-only" rows to other data
obs_df_newdata <- bind_rows(obs_newdata, obs_newdata_dev) %>%
  mutate(fusion_type = fct_inorder(fusion_type))

# Add "deviation-only" rows to other design matrix
obs_new_design <- rbind(obs_new_design, obs_new_design_dev)

# Save results
save(obs_df_newdata, obs_new_design,
     file = here("data", "simulations", "sim_obs_newdata.rda"))
