################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Summarize results from hierarchical model simulation
################################################################################

source(here::here("source", "000_definitions.R"))

# Identify sample of gammas to plot
gammas_to_check <- c(1:2, 6:7, 10:12, 46:48, 64:65)

# Load files --------------------------------------------------------------

# Load model results ("hierarchy_test")
LOAD_GIBBS_DATE <- "2026-06-26"

load(file = here::here("data", "simulations",
                       paste0("hierarchy_test_",
                              LOAD_GIBBS_DATE,
                              ".rda")))

# Load simulated data ("sim_df_model_data")
load(file = here("data", "simulations", "sim_data_model.rda"))

# Load new data to plot model effects on ("obs_df_newdata", "obs_new_design")
load(file = here("data", "simulations", "sim_obs_newdata.rda"))

# Load observation-level data ("obs_df_sim")
load(file = here("data", "simulations", "sim_obs_data.rda"))

# Load dataset of gamma values to plot on new data ("gamma_newdata_df")
load(file = here("data", "simulations", "gamma_df_hierarchical.rda"))

# Load data set of points with row-column coordinates ("point_coords")
load(file = here("data", "intermediate", "point_coordinates.rda"))

# Load region object ("region_shape")
load(file = here("data", "intermediate",
                 "region_shape_object.rda"))

# Load consistent mapping from fusion to colors ("fusion_color_dict")
load(file = here("data", "intermediate",
                 "fusion_color_mapping.rda"))

# Load categories for observation-level parameters ("obs_params_dict")
load(file = here("data", "simulations", "obs_param_ids.rda"))

# Make trace plots --------------------------------------------------------

# Turn parameters into long table for trace plotting
hier_gibbs_tr <- imap(
  # Remove timing from output
  hierarchy_test %>%
    list_modify("timing" = zap(),
                "start-stop" = zap()),
  function(x, xn){

    x_dim <- length(dim(x))

    # Turn array into long format
    x <- reshape2::melt(x) %>% mutate(param_cat = xn)

    # Rename to identify iterations and parameter numbers
    if(x_dim == 2){

      x <- rename(x, param_num = Var1, iter = Var2)

    } else if (x_dim == 3){

      x <- rename(x, row_num = Var1, col_num = Var2, iter = Var3)

    }

    x
  }
) %>%
  # Combine
  list_rbind() %>%
  # Create ID for parameters to be identified
  group_by(param_cat) %>%
  mutate(row_id = str_pad(row_num,
                          ceiling(log10(max(row_num))),
                          "left", "0"),
         col_id = str_pad(col_num,
                          ceiling(log10(max(col_num))),
                          "left", "0"),
         param1_id = str_pad(param_num,
                             ceiling(log10(max(param_num))),
                             "left", "0"),
         is_scalar = if_else(!is.na(max(param_num)) & max(param_num) == 1,
                             1, 0)) %>%
  ungroup() %>%
  # Identify each unique parameter by row and column or by index in matrix
  mutate(param_id = if_else(is.na(param_num),
                            paste0(row_id, "_", col_id),
                            param1_id)) %>%
  select(-c(row_id, col_id, param1_id)) %>%
  # Identify burn-in vs. posterior distribution
  mutate(draw_cat = if_else(str_detect(param_cat, "burn"),
                            "Burn-In", "Post-Burn"),
         param_clean = str_remove(param_cat, "burn_")) %>%
  mutate(parameter = if_else(is_scalar == 1, param_clean,
                             paste0(param_clean, "_", param_id))) %>%
  select(-is_scalar) %>%
  # Add true values
  rowwise() %>%
  mutate(
    value_true = ifelse(
      param_clean == "sigma_sq",  sim_df_model_data$call$sigmasq,
      ifelse(
        param_clean == "tau_sq", sim_df_model_data$call$tausq,
        ifelse(
          param_clean == "lambda_basis", sim_df_model_data$call$lambdas_soap[param_num],
          ifelse(
            param_clean == "lambda_demo", sim_df_model_data$call$lambdas_obs[param_num+1],
            ifelse(
              param_clean == "beta",
              as.matrix(sim_df_model_data$beta)[row_num, col_num],
              ifelse(
                param_clean == "gamma", as.matrix(sim_df_model_data$gamma_true)[row_num, col_num]
              ))))))
  ) %>%
  ungroup()

### Make trace plots

# Trace plot of sigma
basis_trace_hier <- ggplot(hier_gibbs_tr %>%
                             filter(str_detect(parameter, "sigma_sq|lambda_basis"))) +
  geom_line(aes(x = iter, y = value)) +
  geom_hline(aes(yintercept = value_true), color = "blue") +
  facet_wrap(~interaction(draw_cat, parameter),
             ncol = 2, scales = "free") +
  labs(title = "Soap Film Variance Traces")

ggsave(here("results", "sim_plot_trace_scalar_soap.png"),
       basis_trace_hier,
       height = 6, width = 6, units = "in")

# Trace plot of tau
obs_trace_hier <- ggplot(hier_gibbs_tr %>%
                            filter(str_detect(parameter, "tau_sq|lambda_demo"))) +
  geom_line(aes(x = iter, y = value)) +
  geom_hline(aes(yintercept = value_true), color = "blue") +
  facet_wrap(~interaction(draw_cat, parameter),
             ncol = 2, scales = "free") +
  labs(title = "Demographic Variance Traces")

ggsave(here("results", "sim_plot_hier_trace_scalar_obs.png"),
       obs_trace_hier,
       height = 6, width = 6, units = "in")


# Make scatter plots comparing estimates to true --------------------------

# Identify gamma parameters and true values
gamma_vals <- hier_gibbs_tr %>%
  filter(param_clean == "gamma" & draw_cat != "Burn-In") %>%
  group_by(parameter, row_num, col_num, value_true) %>%
  summarize(post_mean = mean(value)) %>%
  ungroup() %>%
  left_join(obs_df_sim %>% mutate(row_num = 1:n())) %>%
  mutate(bias = post_mean - value_true,
         soap_cat = if_else(col_num <= 9, "Boundary", "Interior"))

# Plot gamma estimates against true values
gamma_scatter <- ggplot(gamma_vals,
                        aes(x = value_true,
                            y = post_mean,
                            color = fusion_type,
                            shape = soap_cat)) +
  geom_point(alpha = 0.25) +
  ggpubr::stat_cor() +
  scale_shape_manual(values = c("Boundary" = 17,
                                "Interior" = 16)) +
  facet_grid(fusion_type~soap_cat, scales = "free") +
  theme(legend.position = "bottom") +
  labs(x = "True Value", y = "Estimate",
       title = "Gamma (Soap Film) Coefficients")

ggsave(here("results", "sim_plot_hier_scatter_gammas.png"),
       gamma_scatter,
       height = 6, width = 6, units = "in")

# Identify beta parameters and true values
beta_vals <- hier_gibbs_tr %>%
  filter(param_clean == "beta" & draw_cat != "Burn-In") %>%
  group_by(parameter, row_num, col_num, value_true) %>%
  summarize(post_mean = mean(value)) %>%
  ungroup() %>%
  mutate(bias = post_mean - value_true,
         soap_cat = if_else(col_num <= 9, "Boundary", "Interior")) %>%
  left_join(obs_params_dict, by = c("row_num" = "param_num")) %>%
  mutate(effect_cat = fct_inorder(effect_cat))

beta_scatter <- ggplot(beta_vals,
                       aes(x = value_true,
                           y = post_mean,
                           color = effect_cat,
                           shape = soap_cat)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, alpha = 0.5) +
  geom_point(alpha = 0.25) +
  ggpubr::stat_cor() +
  scale_shape_manual(values = c("Boundary" = 17,
                                "Interior" = 16)) +
  facet_wrap(~soap_cat, scales = "free") +
  theme(legend.position = "bottom") +
  labs(x = "True Value", y = "Estimate",
       title = "Beta (Demographic) Coefficients")

ggsave(here("results", "sim_plot_hier_scatter_betas.png"),
       beta_scatter,
       height = 6, width = 6, units = "in")

# Compare gamma estimates and truth ---------------------------------------

# Extract posterior estimates of gamma coefficients
post_gammas_full <- simplify2array(
  apply(hierarchy_test$beta, 3,
        function(x){ obs_new_design %*% x }, simplify = F))

post_gammas_est <- apply(post_gammas_full, c(1,2), mean)

# Plot as curves with respect to age
gamma_curve_est_df <- post_gammas_est %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(subject = 1:n()) %>%
  pivot_longer(matches("V"), names_to = "Gamma", values_to = "estimate") %>%
  mutate(Gamma = as.numeric(str_remove(Gamma, "V"))) %>%
  left_join(gamma_newdata_df %>% rename(coefficient = value)) %>%
  pivot_longer(c(estimate, coefficient), names_to = "value_type")

gamma_curve_est_plot <- ggplot(gamma_curve_est_df %>%
                                filter(Gamma %in% gammas_to_check &
                                         str_detect(fusion_type, "Only|Norm")),
                              aes(x = age, y = value,
                                  color = fusion_type, linetype = value_type)) +
  geom_line() +
  facet_wrap(~Gamma, labeller = label_both, scales = "free_y") +
  scale_linetype_manual(values = c("coefficient" = 1,
                                   "estimate" = 2)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Age", y = "Gamma Value",
       title = "Simulated Age Effects")

ggsave(here("results", "sim_plot_hier_curves_gammas.png"),
       gamma_curve_est_plot,
       height = 6, width = 6, units = "in")

# Compare predicted growth to true values ---------------------------------

# Find nearest point in image data to each region's centroid
mx_pt_shape <- sf::st_as_sf(point_coords, coords = c("row", "col"))

sample_points <- sf::st_nearest_feature(region_shape$centroid, mx_pt_shape)

# Calculate would-be "true" outcomes for new data
outcomes_newdata <- obs_new_design %*%
  tcrossprod(sim_df_model_data$beta,
             sim_df_model_data$call$soap_design_mx)

# Convert would-be true outcomes to data frame
outcomes_newdata_df <- outcomes_newdata %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(row_num = 1:n()) %>%
  bind_cols(obs_df_newdata) %>%
  pivot_longer(matches("V"), names_to = "cell_num") %>%
  mutate(cell_num = as.numeric(str_remove(cell_num, "V"))) %>%
  # Add sample point identification to true outcomes
  left_join(enframe(sample_points) %>%
              rename(region_num = name,
                     cell_num = value))

# Identify rows in soap design matrix for these sample points
soap_design_sample <- sim_df_model_data$call$soap_design_mx[sample_points,]

# Multiply these by estimated gammas
post_outcomes_full <- simplify2array(
  apply(post_gammas_full, 3,
        function(x){ tcrossprod(x, soap_design_sample) },
        simplify = F))

# Take mean and 95% credible interval(?) from these estimates
post_outcomes_est <- apply(post_outcomes_full, c(1, 2),
                           FUN = function(x) matrix(
                             c(mean(x),
                               quantile(x, 0.05),
                               quantile(x, 0.95)),
                             nrow = 1
                           ))

# Convert this to a dataframe
post_outcomes_df <- post_outcomes_est %>%
  reshape2::melt() %>%
  rename(metric = Var1,
         row_num = Var2,
         region_num = Var3,
         estimate = value) %>%
  mutate(metric = case_match(metric,
                             1 ~ "mean",
                             2 ~ "pct05",
                             3 ~ "pct95"))

# Merge to interpretable data
post_newdata_est <- obs_df_newdata %>%
  mutate(row_num = 1:n()) %>%
  left_join(post_outcomes_df) %>%
  left_join(outcomes_newdata_df %>%
              rename(simulated = value)) %>%
  pivot_longer(c(estimate, simulated), names_to = "value_type")

# Plot estimated vs simulated outcomes ------------------------------------

# Create plot of mean predictions by age
post_estimates_plot <- ggplot(post_newdata_est %>%
                                filter(metric == "mean" &
                                         !str_detect(fusion_type, "Only")),
                              aes(x = age, y = value,
                                  color = fusion_type,
                                  linetype = value_type)) +
  geom_line() +
  facet_wrap(~region_num, nrow = 3, scales = "free_y") +
  scale_linetype_manual(values = c("simulated" = 1,
                                   "estimate" = 2,
                                   "data" = 3)) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(here("results", "sim_plot_hier_estimate_age.png"),
       post_estimates_plot,
       height = 6, width = 6, units = "in")

# Plot estimated vs average simulated on surface --------------------------

# Filter to a subset of sample individuals at interpretable ages
post_surface_data <- obs_df_newdata %>%
  mutate(row_num = row_number()) %>%
  filter(age %in% c(round(seq(0, cranio_max_age, length.out = 7) / 10) * 10))

# Extract gammas specific to these individuals
post_surface_gammas <- post_gammas_full[post_surface_data$row_num, ,]

# Get mean estimates of these gammas
post_gammas_surface_mean <- apply(post_surface_gammas, c(1,2), mean)

# Calculate estimated outcomes for each set of estimated gammas
post_outcomes_surface <- tcrossprod(post_gammas_surface_mean,
                                    sim_df_model_data$call$soap_design_mx)

# Convert to DF and merge in relevant info
post_outcomes_surface_df <- reshape2::melt(post_outcomes_surface) %>%
  rename(row_id = Var1, cell_num = Var2) %>%
  left_join(post_surface_data %>% mutate(row_id = 1:n())) %>%
  left_join(point_coords %>% mutate(cell_num = row_number())) %>%
  mutate(age_header = "Age (Days)")

# Plot these estimates outcomes on surface
outcomes_plot <- ggplot(post_outcomes_surface_df %>%
                          filter(!str_detect(fusion_type, "Only"))) +
  geom_raster(aes(x = row, y = col, fill = value)) +
  facet_nested(age ~ fusion_type) +
  coord_fixed() +
  scale_fill_viridis_c() +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(fill = "Est. Growth") +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 10,
                                    margin = margin(t = 2, b = 3)))

ggsave(here("results", "sim_plot_hier_estimate_surface.png"),
       outcomes_plot,
       height = 6, width = 6, units = "in")
