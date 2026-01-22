library(here)
library(tidyverse)


data <- read_csv(here("processed/merged_covariates.csv"))

# This first subject, single day old, can be used to obtain the PCA predicted "normative" at birth
healthy_data <- data %>%
  filter(Group == "Normative") %>%
  arrange(Age) %>%
  filter(Sex == 1)

reference_study <- (healthy_data %>% pull(Study))[1]
reference_file  <- (healthy_data %>% pull(ID_fmt))[1]

reference_data <- load(paste0(here("processed/individual_subject_files",
                                   paste0(reference_file, ".RData")
)))

# Reference volume
basic_reference_volume = subject_data$closest_normal

# Convert reference volume to an array
basic_reference_volume_arr <- reshape2::acast(basic_reference_volume,
                                              X ~ Y ~ Z, value.var = "value")



# Sanity check
library(plotly)
x_coords <- basic_reference_volume_arr[, , 1]
y_coords <- basic_reference_volume_arr[, , 2]
z_coords <- basic_reference_volume_arr[, , 3]

# Melt into long format
df <- tibble::tibble(
  X = as.vector(x_coords),
  Y = as.vector(y_coords),
  Z = as.vector(z_coords)
)

# Keep only rows where none of the coordinates are -1
df <- df %>%
  filter(X != -1, Y != -1, Z != -1)

plot_ly(
  data = df,
  x = ~X, y = ~Y, z = ~Z,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 2)
)

reference_point <- c(0, 0, 0)
reference_point[1] <- mean(df$X)
reference_point[2] <- mean(df$Y)
reference_point[3] <- min(df$Z)


############

N <- nrow(data)
for (i in 1:N){


  # Extract filename for this child's spherical coordinates
  filename  <- data$RData_file[i]
  data_file <- here("processed/individual_subject_files", filename)

  # Load this child's data (spherical coordinates)
  load(data_file) # stored as subject_data

  # Cast as 3d array to match original format
  data_3d <- reshape2::acast(subject_data$spherical_coordinate_map,
                             X ~ Y ~ Z, value.var = "value")

  # Calculate difference from newborn reference
  # basic_reference_volume_arr is the reference volume
  coords1 <- matrix(data_3d, ncol = 3)
  coords2 <- matrix(basic_reference_volume_arr, ncol = 3)

  # Mask -1 rows
  valid_rows <- apply(coords1, 1, function(x) all(x != -1)) &
    apply(coords2, 1, function(x) all(x != -1))

  # Compute distances only at valid points
  distances <- sqrt(rowSums((coords1[valid_rows, ] - coords2[valid_rows, ])^2))

  # Get the signs as well
  # note - center point is not 0, so just doing the sign doesnt make sense here
  # instead need to see which is "further" away from center point
  # also using the center of Z isn't quite right, but it is close enough to solve the
  # issues with the center point, reference_point
  # center_matrix <- cbind(
  #   rep(reference_point[1], nrow(coords1[valid_rows, ])),
  #   rep(reference_point[2], nrow(coords1[valid_rows, ])),
  #   rep(reference_point[3], nrow(coords1[valid_rows, ]))
  # )

  #dist1 <- rowSums( (coords1[valid_rows, ] - center_matrix )^2)
  #dist2 <- rowSums( (coords2[valid_rows, ] - center_matrix )^2)
  #signs <- sign(dist1 - dist2)

  # Evaluate signed differences
  #signed_distances <- distances * signs

  # Reshape back to matrix
  pointwise_diff             <- matrix(NA, nrow = 500, ncol = 500)
  pointwise_diff[valid_rows] <- distances

  # Save the pointwise differences
  filename  <- paste0(data$ID_fmt[i], "_desc-unsigneddist.RData")
  output_file <- here("processed/unsigned_distances", filename)
  save(file=output_file, pointwise_diff)

  print(output_file)
}
