################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Run simulation to test model fit
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load soap film object ("cranio_soap")
load(file = here("data", "cleaned", "soap_object.rda"))

# Load smooth object for age curves ("obs_smooth_list_sim")
load(file = here("data", "simulations", "sim_obs_smooth_list.rda"))

# Set scalar parameters ---------------------------------------------------

# Scalar parameters for outcome data simulation
sigmasq_sim <- 0.3

tausq_sim <- 2.5e-3

lambda_soap_sim <- c(8e-5, 2e-4)

lambda_obs_sim <- c(1, 1e-7)

# Generate data -----------------------------------------------------------

set.seed(413)

# Use inputs above to generate simulated data based on model
sim_df_model_data <- make_sim_data_hier_soap(
  soap_design_mx = cranio_soap$X,
  obs_design_mx = obs_smooth_list_sim$X,
  penalty_soap_list = cranio_soap$S,
  penalty_obs_list = obs_smooth_list_sim$S,
  sigmasq = sigmasq_sim,
  tausq = tausq_sim,
  lambdas_soap = lambda_soap_sim,
  lambdas_obs = lambda_obs_sim
)

save(sim_df_model_data,
     file = here("data", "simulations", "sim_model_data.rda"))

