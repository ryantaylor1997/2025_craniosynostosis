################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Fit Bayesian GAM of a soap film coefficient on age and fusion type
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load all coefficients fit with GCV, with demographics ("cranio_models_beta")
load(file = here::here("results", "individual_soap_film_fits.rda"))

# Load observation-level design matrix ("obs_smooth_list")
load(file = here::here("data", "cleaned", "obs_level_smooth.rda"))

# Run Gibbs sampler for coefficient models --------------------------------

# Set 1 covariate as outcome
# - (plotted above, should be lower for sagittal and higher for metopic)
coeff_y <- cranio_models_beta$beta_38

# Set default parameters
def_pars_coeff <- c(
  "beta" = 0.0,
  "sigma_sq" = var(coeff_y) / length(coeff_y),
  "lambda" = 1e-3,
  "a" = 1.0,
  "b" = 1.0,
  "c" = 1.0,
  "d" = 1.0
)

# Set iterations
bayes_coeff_iters <- 5e3

# Run function
coeff_bayes1 <- lm_penalized_gibbs(
  y = coeff_y, X = obs_smooth_list$X,
  S = obs_smooth_list$S, has_unpenalized = TRUE,
  beta0 = def_pars_coeff["beta"],
  sigmasq0 = def_pars_coeff["sigma_sq"],
  lambda0 = def_pars_coeff["lambda"],
  a_pri = def_pars_coeff["a"],
  b_pri = def_pars_coeff["b"],
  c_pri = def_pars_coeff["c"],
  d_pri = def_pars_coeff["d"],
  iters = bayes_coeff_iters,
  burn_pct = 0.1
)

save(coeff_bayes1, coeff_y, def_pars_coeff,
     file = here("results", "obs_level_test_coefficient_model.rda"))
