################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Plot results of simulation to test model fit
################################################################################

source(here::here("source", "000_definitions.R"))

# Identify sample of gammas to plot
gammas_to_check <- c(1:2, 6:7, 10:12, 46:48, 64:65)

# Load data ---------------------------------------------------------------

# Load observation-level simulated data ("obs_df_sim")
load(file = here("data", "simulations", "sim_obs_data.rda"))

# Load simulated data ("sim_df_positive")
load(file = here("data", "simulations", "sim_data_positive.rda"))

# Load new data to plot model effects on ("obs_df_newdata", "obs_new_design")
load(file = here("data", "simulations", "sim_obs_newdata.rda"))

# Load data set of points with row-column coordinates ("point_coords")
load(file = here("data", "intermediate", "point_coordinates.rda"))

# Plot true values for simulated parameters -------------------------------

# Calculate gammas as new design matrix multiplied by simulated betas
gamma_newdata <- obs_new_design %*% sim_df_positive$beta

# Match these gammas to age values
gamma_newdata_df <- obs_df_newdata %>%
  bind_cols(
    gamma_newdata %>%
  as.matrix() %>%
  as.data.frame()
  )%>%
  mutate(subject = 1:n()) %>%
  pivot_longer(matches("V"), names_to = "Gamma") %>%
  mutate(Gamma = as.numeric(str_remove(Gamma, "V")))

save(gamma_newdata_df,
     file = here("data", "simulations", "gamma_df_positive.rda"))

# Plot curves, excluding "deviation-only" values
gamma_curve_plot <- ggplot(gamma_newdata_df %>%
                             filter(Gamma %in% gammas_to_check &
                                      !str_detect(fusion_type, "Only")),
                           aes(x = age, y = value, color = fusion_type)) +
  geom_line() +
  facet_wrap(~Gamma, labeller = label_both) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Age", y = "Gamma Value",
       title = "Simulated Gamma by Age")

ggsave(here("results", "sim_plot_pos_gamma_age.png"),
       gamma_curve_plot,
       height = 6, width = 6, units = "in")

# Plot true outcomes ------------------------------------------------------

# Make data frame out of outcome with no error
# (true beta, true gamma, no epsilon)
outcome_long_df <- sim_df_positive$outcome_noerror %>%
  as.matrix() %>%
  reshape2::melt() %>%
  rename(subj_num = Var1, loc_num = Var2) %>%
  left_join(obs_df_sim %>% mutate(subj_num = 1:n())) %>%
  left_join(point_coords %>% mutate(loc_num = 1:n())) %>%
  mutate(age_cut = quantcut(age, 5))

save(outcome_long_df,
     file = here("data", "simulations", "outcome_df_positive.rda"))

# Summarize this by age bucket
outcome_summ_age_fusion <- outcome_long_df %>%
  group_by(fusion_type, age_cut, loc_num, row, col) %>%
  summarize(mean_val = mean(value)) %>%
  ungroup()

# Plot summary by age and fusion
outcome_fusion_age_plot <- ggplot(outcome_summ_age_fusion,
                                  aes(x = row, y = col, fill = mean_val)) +
  geom_raster() +
  facet_grid(age_cut~fusion_type) +
  coord_fixed() +
  scale_fill_viridis_c() +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text = element_text(angle = 45, hjust = 1)) +
  labs(title = "Age-Fusion Summaries of True Outcomes",
       fill = "Avg. Growth")

ggsave(here("results", "sim_plot_pos_outcome_age.png"),
       outcome_fusion_age_plot,
       height = 6, width = 6, units = "in")
