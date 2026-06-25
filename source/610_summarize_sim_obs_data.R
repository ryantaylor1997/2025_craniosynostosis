################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Plot results of simulation to test model fit
################################################################################

source(here::here("source", "000_definitions.R"))


# Load data ---------------------------------------------------------------

# Load observation-level simulated data ("obs_df_sim")
load(file = here("data", "simulations", "sim_obs_data.rda"))

# Load smooth object for age curves ("obs_smooth_list_sim")
load(file = here("data", "simulations", "sim_obs_smooth_list.rda"))

# Load indices of columns that penalty matrix can be divided into ("obs_params_cols")
load(file = here("data", "simulations", "sim_obs_smooth_column_list.rda"))

# Load dense grid of ages to plot model effects on ("obs_df_newdata")
load(file = here("data", "simulations", "sim_obs_newdata.rda"))

# Visualize demographic data ----------------------------------------------

obs_age_plot <- ggplot(obs_df_sim) +
  geom_area(stat = "density", aes(x = age, fill = fusion_type)) +
  facet_grid(fusion_type ~ sex,
             labeller = labeller(.rows = label_value, .cols = label_both)) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(here("results", "sim_plot_obs_data.png"),
       obs_age_plot,
       height = 6, width = 6, units = "in")

# Visualize age smooth bases ----------------------------------------------

# Reshape design matrix for 1 chunk set of knots and merge to age values
obs_expansion_df <- obs_smooth_list_sim$smooth[[1]]$X0 %>%
  reshape2::melt() %>%
  rename(row = Var1, fn = Var2) %>%
  left_join(obs_df_sim %>% mutate(row = 1:n()))

# Plot basis functions
obs_basis_plot <- ggplot(obs_expansion_df) +
  geom_line(aes(x = age, y = value, color = as.factor(fn))) +
  scale_color_manual(values = unname(pals::alphabet2())) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Age", y = "Function Value", color = "Basis Function")

ggsave(here("results", "sim_plot_obs_bases.png"),
       obs_basis_plot,
       height = 6, width = 6, units = "in")

# Visualize age smooth penalty --------------------------------------------

# Set number of non-penalized terms
Q_theta <- length(obs_params_cols[[1]])

# Map each parameter to a penalty number and classify
obs_params_dict <- tibble(
  block_num = 1:length(obs_params_cols),
  param_num = obs_params_cols
) %>%
  unnest(param_num) %>%
  mutate(effect_cat = case_match(block_num,
                                 1 ~ "Theta",
                                 2 ~ "Normative",
                                 3 ~ "Fusion 1",
                                 4 ~ "Fusion 2")) %>%
  select(-block_num)

save(obs_params_dict,
     file = here("data", "simulations", "obs_param_ids.rda"))

# Reshape penalty matrices into long data
obs_pen_expansion_df <- map(
  1:length(obs_smooth_list_sim$S),
  ~ obs_smooth_list_sim$S[[.x]] %>%
    reshape2::melt() %>%
    rename(row = Var1, col = Var2) %>%
    mutate(index = .x)
) %>%
  bind_rows()

# Extract penalty matrix for one set of knots to make it easier to see pattern
obs_pen_single_df <- obs_pen_expansion_df %>%
  filter(index == 2 &
           row %in% obs_params_cols[[2]] &
           col %in% obs_params_cols[[2]])

# Plot penalty matrix for a single set of knots
obs_pen_single_plot <- obs_pen_single_df %>%
  ggplot(aes(x = col, y = -row, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(transform = scales::transform_pseudo_log(sigma = 1e-6),
                       breaks = c(-10^c(-4:0), 0, 10^c(-4:0))) +
  scale_x_continuous(breaks = scales::breaks_width(1)) +
  scale_y_continuous(breaks = scales::breaks_width(1)) +
  theme(legend.position = "bottom",
        legend.text = element_text(angle = 45, hjust = 1)) +
  labs(title = "Age Curve Penalty Matrix")

ggsave(here("results", "sim_plot_obs_penalty_single.png"),
       obs_pen_single_plot,
       height = 6, width = 6, units = "in")

# Plot penalty matrices for entire age curve matrix
obs_pen_full_plot <- obs_pen_expansion_df %>%
  group_split(index) %>%
  map(
    ~ggplot(., aes(x = col, y = -row, fill = value)) +
      geom_tile() +
      facet_wrap(~index, nrow = 1) +
      scale_fill_gradient2(transform = scales::transform_pseudo_log(sigma = 0.001),
                           breaks = c(-1, -0.1, 0, 0.1, 1)) +
      scale_x_continuous(breaks = scales::breaks_width(4)) +
      scale_y_continuous(breaks = scales::breaks_width(4)) +
      theme_minimal() +
      theme(legend.position = "bottom",
            legend.text = element_text(angle = 45, hjust = 1)) +
      coord_fixed()
  ) %>%
  ggarrange(plotlist = ., nrow = 1)

ggsave(here("results", "sim_plot_obs_penalty_full.png"),
       obs_pen_full_plot,
       height = 6, width = 6, units = "in")

