################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Test single-individual Bayesian soap film model
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load sample data we used to fit frequentist soap models ("so_test_df")
load(file = here("data", "intermediate", "soap_test_sample_data.rda"))

# Load soap film object ("cranio_soap")
load(file = here("data", "cleaned", "soap_object.rda"))

# Run Bayesian soap film on 1 individual -----------------------------

### Run Gibbs sampler for soap film

# Set data to train on, especially outcome
data_eg <- so_test_df %>%
  filter(fusion_type == "Sagittal") %>%
  filter(fname == first(fname)) %>%
  unnest(data)

save(data_eg,
     file = here("data", "intermediate", "soap_test_single_obs_data.rda"))

y_mx <- data_eg$diff

# Set default parameters
def_pars <- c(
  "beta" = 0.0,
  "sigma_sq" = var(y_mx),
  "lambda" = 0.1,
  "a" = 1.0,
  "b" = 1.0,
  "c" = 1.0,
  "d" = 1.0
)

# Set iterations
bayes_soap_iters <- 1e3

# Run function
soap_bayes1 <- lm_penalized_gibbs(
  y = y_mx, X = cranio_soap$X, S = cranio_soap$S,
  beta0 = def_pars["beta"],
  sigmasq0 = def_pars["sigma_sq"],
  lambda0 = def_pars["lambda"],
  a_pri = def_pars["a"],
  b_pri = def_pars["b"],
  c_pri = def_pars["c"],
  d_pri = def_pars["d"],
  iters = bayes_soap_iters,
  burn_pct = 0.1
)

save(soap_bayes1,
     file = here("data", "intermediate", "soap_test_single_obs_bayes.rda"))
