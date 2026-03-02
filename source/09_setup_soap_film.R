
# Set up boundary object --------------------------------------------------

## Create boundary around cranium shape

# NB: Created polygon around pixels and points along boundary before this


# Visualize how boundary compares to all points
bound_viz <- ggplot() +
  geom_point(data = cranio_matrix %>% select(row, col),
             aes(x = row, y = col), alpha = 0.5, size = 0.0001) +
  geom_point(data = cranio_bound_matrix,
             aes(x = row, y = col, color = appearance_order)) +
  scale_color_viridis_c(option = "mako") +
  coord_fixed() +
  labs(title = "Boundary Sanity Check",
       color = "Order of Appearance") +
  theme_minimal() +
  theme(legend.position = "bottom")


# Set up soap film infrastructure -----------------------------------------

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

# Keep only what we need to build smoother
so_matrix_pts <- cranio_matrix %>%
  select(row, col)

save(so_matrix_pts, file = here::here("analysis", "intermediate", "matrix_indices.rda"))


# Create soap film object -------------------------------------------------

# Construct soap film smoother object
cranio_soap <- smoothCon(
  s(
    row, col,
    bs = "so",
    xt = list(
      bnd = bound_list)
  ),
  data = so_matrix_pts, knots = so_grid_k
)

# Keep the actual smooth list
cranio_soap <- cranio_soap[[1]]

# Save soap film object
save(cranio_soap,
     file = here::here("analysis", "intermediate", "soap_object.rda"))

# Create separated soap film matrices -------------------------------------

# Define model matrix
so_X <- cranio_soap$X

# Identify sample size
so_n <- nrow(so_X)

# Identify boundary and film penalty matrices
S_bound0 <- cranio_soap$S[[1]]
S_film0 <- cranio_soap$S[[2]]

# Identify different parts of penalty matrix
cols_bound <- which(apply(S_bound0, 2, function(x) sum(x != 0)) != 0)
cols_film <- which(apply(S_film0, 2, function(x) sum(x != 0)) != 0)

# Get interior "film" design and penalty matrices
so_X_film <- so_X[, cols_film]
S_film <- S_film0[cols_film, cols_film]

save(so_X_film,
     file = here::here("analysis", "intermediate", "design_matrix_interior.rData"))
