################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Plot soap film basis functions
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load soap film object ("cranio_soap")
load(file = here("data", "cleaned", "soap_object.rda"))

# Load data set of points with row-column coordinates ("point_coords")
load(file = here("data", "intermediate", "point_coordinates.rda"))

# Load shape of surface ("cranio_bound_shape")
load(file = here("data", "intermediate", "boundary_polygon.rda"))

# Load boundary points and list object ("cranio_bound_points", "bound_list")
load(file = here("data", "intermediate", "boundary_objects.rda"))

# Load grid of knots ("so_grid_k")
load(file = here("data", "intermediate", "soap_knot_grid.rda"))

# Load indices identifying types of basis fns ("cols_bound", "cols_film")
load(file = here("data", "intermediate", "soap_basis_indices.rda"))

# Set up soap film knots so they can be plotted ---------------------------

# Link knots to interior basis functions
so_grid_k %<>% mutate(basis_fn = cols_film)

# Identify sample of interior knots to plot (3 groups of 5)
plot_knots <- cols_film[c(1:5,
                          (round(length(cols_film)/2) + c(-2:2)),
                          ((length(cols_film) - 4):length(cols_film)))]

# Replicate distance along boundary ---------------------------------------

### Link boundary points with nearest matrix point
# Boundary is buffered 1 unit away from matrix; need to match it with real points

# Convert boundary to sf for nearest neighbor matching
# Add in distances determined in smoother object
bound_dists <- cranio_bound_points %>%
  mutate(bound_r = cranio_soap$sd$bnd[[1]]$d[-1]) %>%
  st_as_sf(coords = c("row", "col"), remove = F)

# Convert points within matrix to sf for nearest neighbor
points_sf <- point_coords %>%
  st_as_sf(coords = c("row", "col"), remove = F)

# Identify nearest matrix neighbor to each boundary point
bound_dists %<>%
  mutate(nearest_cell = st_nearest_feature(bound_dists, points_sf)) %>%
  st_drop_geometry() %>%
  select(nearest_cell, bound_r)

# Match matrix points to their order on the boundary
points_dist_df <- point_coords %>%
  mutate(cell_id = 1:n()) %>%
  left_join(bound_dists,
            by = c("cell_id" = "nearest_cell")) %>%
  # When we get multiple matches, take first point on boundary
  arrange(cell_id, bound_r) %>%
  group_by(cell_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(-cell_id)

# Create basis function data set for plotting -----------------------------

# Create data of basis functions and pixel locations for plotting
basis_df <- cranio_soap$X %>%
  data.frame() %>%
  bind_cols(points_dist_df) %>%
  pivot_longer(-c(row, col, region, bound_r),
               names_to = "basis_fn") %>%
  mutate(basis_fn = as.numeric(str_remove(basis_fn, "X")),
         fn_cat = case_when(basis_fn %in% cols_bound ~ "Boundary",
                            basis_fn %in% cols_film ~ "Interior"))

# Plot boundary cyclic bases ----------------------------------------------

# Plot boundary cyclic 1-dimensional splines
soap_basis_cyclic <- ggplot(data = basis_df %>%
                              filter(!is.na(bound_r) & fn_cat == "Boundary")) +
  geom_line(aes(x = bound_r, y = value, color = factor(basis_fn))) +
  scale_color_brewer(palette = "Paired") +
  labs(x = "Distance Along Boundary",
       y = "Function Value",
       color = "Basis Function",
       title = "Boundary Cyclic Basis Functions") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = T))

# Plot boundary-induced functions -----------------------------------------

# Plot boundary-induced 2-dimensional basis functions
soap_basis_bd <- ggplot() +
  geom_raster(data = basis_df %>% filter(fn_cat == "Boundary"),
              aes(x = row, y = col, fill = value)) +
  scale_fill_viridis_c() +
  labs(x = NULL, y = NULL, fill = "Fn. Value",
       title = "Boundary Basis Functions") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom") +
  facet_wrap(~basis_fn)

# Plot some interior functions --------------------------------------------

# Plot interior basis functions with corresponding knots
soap_basis_int <- ggplot() +
  geom_raster(data = basis_df %>% filter(basis_fn %in% plot_knots),
              aes(x = row, y = col, fill = value)) +
  geom_point(data = so_grid_k %>% select(-basis_fn),
             aes(x = row, y = col),
             color ="black", size = 0.01) +
  geom_point(data = so_grid_k %>% filter(basis_fn %in% plot_knots),
             aes(x = row, y = col),
             color ="red", size = 0.01) +
  scale_fill_viridis_c() +
  labs(x = NULL, y = NULL, fill = "Fn. Value",
       title = "Interior Basis Functions") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom") +
  facet_wrap(~basis_fn, ncol = 5)
