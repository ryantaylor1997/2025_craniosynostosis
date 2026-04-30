################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Create soap film smoother object
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load cleaned / filtered data of nested growth matrices ("cranio_matrix_clean")
load(file = here("data", "intermediate", "point_coordinates.rda"))

# Set up boundary object --------------------------------------------------

# Create polygon (concave hull) around pixels in cranium shape
cranio_bound_shape <- matrix_coords %>%
  st_as_sf(coords = c("row", "col")) %>%
  concaveman(concavity = 1)

# Create boundary as points
cranio_bound_matrix <- cranio_bound_shape %>%
  # Add buffer so boundary of the shape is slightly beyond data domain
  st_buffer(1) %>%
  # Convert to line
  st_boundary() %>%
  # Take points along this line at density = 1 unit
  st_line_sample(density = 1) %>%
  st_coordinates() %>%
  as_tibble() %>%
  rename(row = X, col = Y) %>%
  # Add order of appearance on boundary
  mutate(appearance_order = 1:n())

# Visualize how boundary compares to all points
bound_viz <- ggplot() +
  geom_point(data = matrix_coords,
             aes(x = row, y = col), alpha = 0.5, size = 0.0001) +
  geom_point(data = cranio_bound_matrix,
             aes(x = row, y = col, color = appearance_order)) +
  scale_color_viridis_c(option = "mako") +
  coord_fixed() +
  labs(title = "Boundary Sanity Check",
       color = "Order of Appearance") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Set up soap film knots ---------------------------------------------------

## Create boundary as list
bound_list <- list(
  list(row = cranio_bound_matrix$row,
       col = cranio_bound_matrix$col)
)

# Identify row knot points
grid_row <- seq(min(cranio_bound_matrix$row), max(cranio_bound_matrix$row),
                length.out = so_GRID_N)

# Identify column knot points
grid_col <- seq(min(cranio_bound_matrix$col), max(cranio_bound_matrix$col),
                length.out = so_GRID_N)

# Expand to grid of knots
grid_all <- expand.grid("row" = grid_row, "col" = grid_col)

# Crop to cranial boundary
so_grid_k <- grid_all %>%
  st_as_sf(coords = c("row", "col")) %>%
  st_filter(cranio_bound_shape) %>%
  st_coordinates() %>%
  as_tibble() %>%
  rename(row = X, col = Y)

# Visualize grid
grid_viz <- ggplot() +
  geom_point(data = cranio_bound_matrix, aes(x = row, y = col), color = "blue") +
  geom_point(data = grid_all, aes(x = row, y = col)) +
  geom_point(data = so_grid_k, aes(x = row, y = col), color ="red") +
  coord_fixed() +
  labs(title = "Grid Sanity Check")

# Create soap film object -------------------------------------------------

# Construct soap film smoother object
cranio_soap <- smoothCon(
  s(
    row, col, # axes
    bs = "so", # basis
    xt = list(
      bnd = bound_list) # boundary
  ),
  data = matrix_coords, # coordinates
  knots = so_grid_k # knots
)

# Keep the actual smooth list
cranio_soap <- cranio_soap[[1]]

cranio_soap_unscaled <- cranio_soap

save(cranio_soap_unscaled,
     file = here::here("data", "intermediate", "soap_object_unscaled.rda"))

# Scale penalty matrices so trace of inverse = trace of identity matrix
soap_pen_eigen <- map(cranio_soap$S,
                      ~eigen(., symmetric = T, only.values = T)$values)

soap_pen_scale <- map_dbl(soap_pen_eigen,
                          ~sum(1 / .[which(abs(.) > 1e-10)]) / sum(. != 0))

cranio_soap$S <- map2(cranio_soap$S, soap_pen_scale,  ~ .x * .y)


# Save output -------------------------------------------------------------

# Save soap film object
save(cranio_soap,
     file = here::here("data", "cleaned", "soap_object.rda"))

# Save penalty scale terms
save(soap_pen_scale,
     file = here::here("data", "intermediate", "soap_penalty_scales.rda"))

# Save shape of surface
save(cranio_bound_shape,
     file = here("data", "intermediate", "boundary_polygon.rda"))

# Save boundary points and list object
save(cranio_bound_matrix, bound_list,
     file = here("data", "intermediate", "boundary_objects.rda"))

# Save grid of knots
save(so_grid_k,
     file = here("data", "intermediate", "soap_knot_grid.rda"))
