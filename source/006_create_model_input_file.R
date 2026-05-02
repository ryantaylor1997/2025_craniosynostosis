################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Combine all model inputs into one file
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load soap film object ("cranio_soap")
load(file = here("data", "cleaned", "soap_object.rda"))

# Load observation-level smooth design matrix and penalty ("obs_smooth_list")
load(file = here::here("data", "cleaned", "obs_level_smooth.rda"))

# Load outcome matrix ("growth_mx")
load(file = here("data", "cleaned", "point_growth_matrix.rda"))

save(cranio_soap, obs_smooth_list, growth_mx,
     file = here("data", "cleaned", "joint_model_inputs.rda"))
