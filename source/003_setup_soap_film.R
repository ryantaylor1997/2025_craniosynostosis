################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Create soap film smoother object
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load data set of points with row-column coordinates ("point_coords")
load(file = here("data", "intermediate", "point_coordinates.rda"))

# Set up boundary object --------------------------------------------------

# Create polygon (concave hull) around pixels in cranium shape
cranio_bound_shape <- point_coords %>%
  select(row, col) %>%
  st_as_sf(coords = c("row", "col")) %>%
  concaveman(concavity = 1)

# Create boundary as points
cranio_bound_points <- cranio_bound_shape %>%
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
  geom_point(data = point_coords,
             aes(x = row, y = col), alpha = 0.5, size = 0.0001) +
  geom_point(data = cranio_bound_points,
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
  list(row = cranio_bound_points$row,
       col = cranio_bound_points$col)
)

# Identify row knot points
grid_row <- seq(min(cranio_bound_points$row), max(cranio_bound_points$row),
                length.out = so_GRID_N)

# Identify column knot points
grid_col <- seq(min(cranio_bound_points$col), max(cranio_bound_points$col),
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
  geom_point(data = cranio_bound_points, aes(x = row, y = col), color = "blue") +
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
  data = point_coords, # coordinates
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

# Save full soap film output -------------------------------------------

# Save soap film object
save(cranio_soap,
     file = here("data", "cleaned", "soap_object.rda"))

# Save penalty scale terms
save(soap_pen_scale,
     file = here("data", "intermediate", "soap_penalty_scales.rda"))

# Save shape of surface
save(cranio_bound_shape,
     file = here("data", "intermediate", "boundary_polygon.rda"))

# Save boundary points and list object
save(cranio_bound_points, bound_list,
     file = here("data", "intermediate", "boundary_objects.rda"))

# Save grid of knots
save(so_grid_k,
     file = here("data", "intermediate", "soap_knot_grid.rda"))

# Create separated soap film matrices -------------------------------------

# Define model matrix
so_X <- cranio_soap$X

# Identify boundary and film penalty matrices
S_bound0 <- cranio_soap$S[[1]]
S_film0 <- cranio_soap$S[[2]]

# Identify different parts of penalty matrix
cols_bound <- which(apply(S_bound0, 2, function(x) sum(x != 0)) != 0)
cols_film <- which(apply(S_film0, 2, function(x) sum(x != 0)) != 0)

# Get interior "film" design and penalty matrices
so_X_film <- so_X[, cols_film]
S_film <- S_film0[cols_film, cols_film]

# Save separated soap components ------------------------------------------

save(so_X,
     file = here("data", "intermediate", "soap_design_full.rda"))

save(so_X_film,
     file = here("data", "intermediate", "soap_design_interior.rda"))

save(cols_bound, cols_film,
     file = here("data", "intermediate", "soap_basis_indices.rda"))

save(S_film,
     file = here("data", "intermediate", "soap_penalty_interior.rda"))

save(S_bound0, S_film0,
     file = here("data", "intermediate", "soap_penalty_separates.rda"))

# Understand derivative matrix --------------------------------------------

# Reconstruct PDE derivative matrix D
# - where penalty S = D'D
soap_D <- t(cranio_soap$sd$P) %*%
  cranio_soap$sd$L %*% cranio_soap$sd$U %*%
  cranio_soap$sd$Q

save(soap_D,
     file = here("data", "intermediate", "soap_derivative_mx.rda"))
