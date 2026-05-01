################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Create large outcome data for hierarchical Craniosynostosis model
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load cleaned data of growth, nested by pixel ("cranio_points_clean")
load(file = here("data", "cleaned", "point_data_clean.rda"))

# Extract nested growth outcomes ------------------------------------------

# Extract column of pointwise differences for each pixel
# Column-bind together to get 1 row for each image / 1 column for each pixel
growth_mx <- do.call(
  cbind,
  hoist(cranio_point_clean,
        "data", "diff")$diff
  )


# Save output -------------------------------------------------------------

save(growth_mx,
     file = here("data", "cleaned", "point_growth_matrix.rda"))
