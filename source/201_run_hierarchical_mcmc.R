################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Run joint Gibbs sampler function on data
################################################################################

suppressMessages(
  source(here::here("source", "000_definitions.R"))
  )

wd = getwd()

if(substring(wd, 2, 6) == "Users"){
  doLocal = TRUE
}else{
  doLocal = FALSE
}

# Load files --------------------------------------------------------------

path_start <- if(doLocal){ "data" } else { "../data_clean" }

# Load all inputs for joint model
# ("cranio_soap", "obs_smooth_list", "growth_mx")
load(file = here(path_start, "cleaned", "joint_model_inputs.rda"))

# Convert large outcome matrix to sparser format
growth_mx <- Matrix(growth_mx)

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

