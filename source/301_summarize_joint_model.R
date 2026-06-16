################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Summarize joint Gibbs sampler function on data
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load most recent joint model fit ("cranio_joint")
FIT_DATE <- "2026-06-11"
load(file = here("results",
                 paste0("joint_model_fit_",
                        format(Sys.time(), "%Y-%m-%d"),
                        ".rda")))

# Load data set of points with row-column coordinates ("point_coords")
load(file = here("data", "intermediate", "point_coordinates.rda"))

# Load region object ("matrix_region")
load(file = here("analysis", "intermediate",
                 "region_shape_object.rda"))

# Load cleaned / filtered subject-level data ("cranio_clean")
load(file = here("data", "cleaned", "obs_data_clean.rda"))

# Load observation-level smooth design matrix and penalty ("obs_smooth_list")
load(file = here("data", "cleaned", "obs_level_smooth.rda"))

# Load soap film object ("cranio_soap")
load(file = here("data", "cleaned", "soap_object.rda"))

# Load consistent mapping from fusion to colors ("fusion_color_dict")
load(file = here("data", "intermediate",
                 "fusion_color_mapping.rda"))

# Check MCMC output -------------------------------------------------------

# Check timing
joint_time <- cranio_joint$timing %>%
  mutate(pct = as.numeric(secs) / last(secs)) %>%
  arrange(desc(pct))

# Extract values by iteration for trace plots
# Turn parameters into long table for trace plotting
joint_gibbs_tr <- imap(
  # Remove timing from output
  cranio_joint %>%
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
  select(-is_scalar)

# Trace plot of sigma
soap_trace_joint <- ggplot(joint_gibbs_tr %>%
                             filter(str_detect(parameter, "sigma_sq|lambda_basis"))) +
  geom_line(aes(x = iter, y = value)) +
  facet_wrap(~interaction(draw_cat, parameter),
             ncol = 2, scales = "free") +
  labs(title = "Soap Film Variance Traces")

ggsave(here("results", paste0(FIT_DATE, "_trace_soap.png")),
       soap_trace_joint,
       height = 6, width = 6, units = "in")

# Trace plot of tau
obs_trace_joint <- ggplot(joint_gibbs_tr %>%
                            filter(str_detect(parameter, "tau_sq|lambda_demo"))) +
  geom_line(aes(x = iter, y = value)) +
  facet_wrap(~interaction(draw_cat, parameter),
             ncol = 2, scales = "free") +
  labs(title = "Observation-Level Variance Traces")

ggsave(here("results", paste0(FIT_DATE, "_trace_obs.png")),
       obs_trace_joint,
       height = 6, width = 6, units = "in")

# Visualize growth by age and suture fusion type --------------------------

# Find nearest point in image data to each region's centroid
mx_pt_shape <- st_as_sf(point_coords, coords = c("row", "col"))

sample_points <- st_nearest_feature(matrix_region$centroid, mx_pt_shape)

## Generate sample data
joint_newdata <- cranio_clean %>%
  mutate(min_age = min(age),
         max_age = max(age)) %>%
  distinct(fusion_type, min_age, max_age) %>%
  rowwise() %>%
  mutate(age = list(seq(min_age, max_age, by = 10))) %>%
  ungroup() %>%
  select(-c(min_age, max_age)) %>%
  mutate(sex = 0) %>%
  unnest(age) %>%
  mutate(sex = factor(sex, levels = levels(cranio_clean$sex)),
         fusion_type = factor(fusion_type, levels = levels(cranio_clean$fusion_type)))

# Convert sample data to design matrix
# Create predictions by fusion type
joint_new_X_list <- map(obs_smooth_list$smooth,
                        ~PredictMat(.x, data = joint_newdata))

# Combine to create Normative and add unpenalized terms
joint_new_design <- Reduce(
  "cbind", c(
    list(
      model.matrix(~sex + fusion_type, data = joint_newdata),
      Reduce("+", joint_new_X_list)),
    joint_new_X_list[-1]
  ))

## Calculate estimated gammas for this sample data from posterior betas
joint_new_gammas <- simplify2array(
  apply(cranio_joint$beta, 3,
        function(x){ joint_new_design %*% x }, simplify = F))

## Calculate estimated outcomes from these estimated gammas

# Extract the design matrix for the sample points we've identified
basis_sample <- cranio_soap$X[sample_points,]

# Calculate estimated outcomes for each set of estimated gammas
joint_outcomes_full <- simplify2array(
  apply(joint_new_gammas, 3,
        function(x){ tcrossprod(x, basis_sample) },
        simplify = F))

# Take mean and 95% credible interval(?) from these estimates
joint_outcomes_est <- apply(joint_outcomes_full, c(1, 2),
                            FUN = function(x) matrix(
                              c(median(x),
                                quantile(x, 0.05),
                                quantile(x, 0.95)),
                              nrow = 1
                            ))

# Reshape estimates
joint_outcomes_df <- joint_outcomes_est %>%
  reshape2::melt() %>%
  rename(metric = Var1,
         row_num = Var2,
         region_num = Var3,
         estimate = value) %>%
  mutate(metric = case_match(metric,
                             1 ~ "median",
                             2 ~ "pct05",
                             3 ~ "pct95"))

# Merge to interpretable data
joint_newdata_est <- joint_newdata %>%
  mutate(row_num = 1:n()) %>%
  left_join(joint_outcomes_df)

# Add scale for alpha
fusion_color_dict %<>%
  mutate(alpha = if_else(fusion_type == "Normative", 1, 0.5))

# Create plots of mean predictions
joint_estimates_plot <- ggplot(joint_newdata_est %>%
                                 filter(metric == "median"),
                               aes(x = age, y = estimate,
                                   color = fusion_type,
                                   alpha = fusion_type)) +
  geom_line() +
  facet_wrap(~region_num, nrow = 2, scales = "free_y") +
  scale_color_manual(values = fusion_color_dict %>%
                       select(fusion_type, color) %>%
                       deframe()) +
  scale_alpha_manual(values = fusion_color_dict %>%
                       select(fusion_type, alpha) %>%
                       deframe()) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Age", y = "Median Estimate",
       color = "Fusion", alpha = "Fusion")

ggsave(here("results", paste0(FIT_DATE, "_pred_paths.png")),
       joint_estimates_plot,
       height = 6, width = 10, units = "in")

# Plot sample predictions from model --------------------------------------

# Filter to a subset of sample individuals at interpretable ages
joint_surface_data <- joint_newdata %>%
  mutate(row_num = row_number()) %>%
  filter(age %in% seq(30, 330, 60))

# Extract gammas specific to these individuals
joint_surface_gammas <- joint_new_gammas[joint_surface_data$row_num,,]

joint_gammas_surface_median <- apply(joint_surface_gammas, c(1,2), median)

# Calculate estimated outcomes for each set of estimated gammas
joint_outcomes_surface <- tcrossprod(joint_gammas_surface_median, cranio_soap$X)

# Convert to DF and merge in relevant info
joint_outcomes_surface_df <- reshape2::melt(joint_outcomes_surface) %>%
  rename(row_num = Var1, cell_num = Var2) %>%
  left_join(joint_surface_data) %>%
  left_join(point_coords %>% mutate(cell_num = row_number())) %>%
  mutate(age_header = "Age (Days)")

# Plot outcomes by suture fusion and age
outcomes_plot_sag <- ggplot(joint_outcomes_surface_df %>%
                              filter(fusion_type %in% c("Normative", "Sagittal"))) +
  geom_raster(aes(x = row, y = col, fill = value)) +
  facet_nested(fusion_type ~ age_header + age,
               nest_line = element_line(linetype = 3)) +
  coord_fixed() +
  scale_fill_viridis_c(option = "turbo") +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(fill = "Est. Growth") +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 10,
                                    margin = margin(t = 2, b = 3)))

ggsave(here("results", paste0(FIT_DATE, "_pred_shapes_sag.png")),
       outcomes_plot_sag,
       height = 2, width = 6, units = "in")

outcomes_plot <- ggplot(joint_outcomes_surface_df) +
  geom_raster(aes(x = row, y = col, fill = value)) +
  facet_nested(fusion_type ~ age_header + age,
               nest_line = element_line(linetype = 3)) +
  coord_fixed() +
  scale_fill_viridis_c(option = "turbo") +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(fill = "Est. Growth") +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 10,
                                    margin = margin(t = 2, b = 3)))

ggsave(here("results", paste0(FIT_DATE, "_pred_shapes_all.png")),
       outcomes_plot,
       height = 6, width = 6, units = "in")
