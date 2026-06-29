################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Plot results from test frequentist soap film
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load fits only on interior using for-loop code ("film_loop_df")
load(file = here("data", "intermediate", "soap_test_interior_loop.rda"))

# Load fits only on interior using nest/map code ("film_map_df")
load(file = here("data", "intermediate", "soap_test_interior_map.rda"))

# Load fits including boundary using for-loop code ("soap_loop_df")
load(file = here("data", "intermediate", "soap_test_full_loop.rda"))

# Load fits includind boundary using nest/map code ("soap_map_df")
load(file = here("data", "intermediate", "soap_test_full_map.rda"))

# Load sample data we used to fit these models ("so_test_df")
load(file = here("data", "intermediate", "soap_test_sample_data.rda"))

# Load range of values we will use for consistency ("so_diff_range")
load(file = here("data", "intermediate", "soap_test_sample_range.rda"))

# Summarize model of only interior, using for loop ------------------------

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


# Summarize model of only interior, using nest/map ------------------------

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


# Summarize model fit on interior and boundary using for-loop -------------

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


# Summarize model with interior and boundary, using nest/map --------------

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
       title = "Sample Predictions using Nested Data") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom") +
  facet_wrap(~fusion_type, nrow = 1)

# Combine and print
map_soap_both <- ggarrange(plot_map_soap,
                           plot_so_test_data,
                           ncol = 1)

# Show building blocks of soap film model ---------------------------------

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

ggsave(here("results", "soap_film_fit_building_blocks.png"),
       plot_map_build,
       height = 6, width = 6, units = "in")
