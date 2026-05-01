# Define function to read each file
load_matrix_data <- function(x, folder){

  # Path to data file
  full_path <- file.path(folder, x)

  # Load RData file (pointwise_diff)
  load(full_path)

  # Convert matrix to long format and remove values outside mask
  reduced_data <- pointwise_diff %>%
    reshape2::melt(varnames = c("row", "col"), value.name = "diff") %>%
    as_tibble() %>%
    drop_na(diff)

  return(list(reduced_data))
}
