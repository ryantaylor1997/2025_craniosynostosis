################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Run Hierarchical Gibbs model once to test on simulated data
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load simulated data ("sim_df_positive")
load(file = here("data", "simulations", "sim_data_positive.rda"))

# Set initial values ------------------------------------------------------

# Set true values as initial gamma values
gamma_init <- as(sim_df_positive$gamma_true, "sparseVector")

set.seed(617)

gamma_init <- gamma_init +
  rnorm(length(gamma_init), 0, sd(gamma_init))

gamma_init <- as.vector(gamma_init)

# Run model ---------------------------------------------------------------

set.seed(603)

positive_test <- hierarchical_penalized_gibbs(
  outcome_mx = sim_df_positive$outcome,
  basis_mx = sim_df_positive$call$soap_design_mx,
  demo_mx = sim_df_positive$call$obs_design_mx,
  pen_basis = sim_df_positive$call$penalty_soap_list,
  pen_demo = sim_df_positive$call$penalty_obs_list,
  basis_has_unpenalized = FALSE,
  demo_has_unpenalized = TRUE,
  gamma0 = 0,
  beta0 = 0,
  sigmasq0 = 1,
  tausq0 = 1,
  lambda_basis0 = 0.1,
  lambda_demo0 = 0.01,
  a_sigma_pri = 1,
  b_sigma_pri = 1,
  a_gamma_pri = 1,
  b_gamma_pri = 1,
  a_tau_pri = 1,
  b_tau_pri = 1,
  a_alpha_pri = 1,
  b_alpha_pri = 1,
  iters = 200,
  burn_pct = 0.5,
  update_pct = 0.1
)

save(positive_test,
     file = here::here("data", "simulations",
                       paste0("positive_test_",
                              format(Sys.time(), "%Y-%m-%d"),
                              ".rda")))
