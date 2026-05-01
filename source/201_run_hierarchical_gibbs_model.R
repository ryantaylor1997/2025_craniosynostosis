################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Run joint Gibbs sampler function on data
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load outcome matrix ("growth_mx")
load(file = here("data", "cleaned", "point_growth_matrix.rda"))

# Load soap film object ("cranio_soap")
load(file = here("data", "cleaned", "soap_object.rda"))

# Load observation-level smooth design matrix and penalty ("obs_smooth_list")
load(file = here::here("data", "cleaned", "obs_level_smooth.rda"))

cranio_growth_mx <- Matrix(cranio_growth_mx)

if(DO_JOINT_FIT){

# Run hierarchical model
cranio_joint <- hierarchical_penalized_gibbs(
  outcome_mx = growth_mx,
  basis_mx = cranio_soap$X,
  demo_mx = obs_smooth_list$X,
  pen_basis = cranio_soap$S,
  pen_demo = obs_smooth_list$S,
  basis_has_unpenalized = FALSE,
  demo_has_unpenalized = TRUE,
  gamma0 = 0,
  beta0 = 0,
  sigmasq0 = 1,
  tausq0 = 1,
  lambda_basis0 = 0.1,
  lambda_demo0 = 1e-4,
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
  update_pct = 0.05
)

save(cranio_joint,
     file = here("results",
                 paste0("joint_model_fit_",
                        format(Sys.time(), "%Y-%m-%d"),
                        ".rda")))

