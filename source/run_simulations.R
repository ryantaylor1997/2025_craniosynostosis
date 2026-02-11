################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Run Simulations for Hierarchical Bayes Model
################################################################################

rm(list = ls()); gc()

# Initialize --------------------------------------------------------------

# Load functions for design matrix and simulations
source(here::here("source", "construct_smooth_ref_code_fns.R"))
source(here::here("source", "simulate_data_fn.R"))
source(here::here("source", "mv_normal_matrix_fn.R"))

## Load soap film object for basis functions
load(file = here::here("analysis", "intermediate", "soap_object.rda"))

### Set parameters

## Parameters for generating demographic data
# Number of patients to use
n_pats <- 80

# Max age of these patients
age_max <- 365

# Proportion of sex = 1
prop_sex <- 0.5

# Knots to use for each age expansion
demo_age_knots <- 10

# Scalar parameters for outcome data simulation
sigma2_test <- 2

tau2_test <- 3

lambda_basis_vec <- c(0.001, 0.01)

lambda_demo <- 0.005

# Take subset of data for demographics ------------------------------------

# Generate a set of patients of a few fusion types

# Set random seed
set.seed(978)

cranio_demo_test_df <- tibble(
  sex = rbinom(n_pats, 1, prop_sex),
  age = runif(n_pats, max = age_max),
  fusion_type = c(rep("Normative", n_pats / 2),
                  rep("Fusion_1", n_pats / 4),
                  rep("Fusion_2", n_pats / 4))
) %>%
  mutate(fusion_type = factor(fusion_type,
                              levels = c("Normative", "Fusion_1", "Fusion_2")))


# Convert demographic data to model matrix --------------------------------

cranio_demo_smooth_test <- smooth.construct2(
  s(age, bs = "cr", k = demo_age_knots),
  data = cranio_demo_test_df, knots = NULL
)

# Attach name of the dataset to this smooth
attr(cranio_demo_smooth_test, "data") <- "cranio_demo_test_df"

# Use function to reference-code smooth by fusion type
demo_ingredients_test <- setup_reference_model(
  sm = cranio_demo_smooth_test,
  by_var = "fusion_type", param_formula = ~sex
)

# Combine suture fusion penalties so we estimate the same lambda
demo_penalty_list_test <- list(demo_ingredients_test$S[[1]],
                               Reduce("+", demo_ingredients_test$S[-1]))

# Simulate corresponding data ---------------------------------------------

### Generate simulated data from these assumptions
sim_data_soap <- make_sim_data_bayes_soap(
  basis_mx = cranio_soap$X,
  demo_mx = demo_ingredients_test$X,
  penalty_basis_list = cranio_soap$S,
  penalty_demo_list = demo_penalty_list_test,
  sigmasq = sigma2_test,
  tausq = tau2_test,
  lambdas_basis = lambda_basis_vec,
  lambdas_demo = lambda_demo,
  demo_has_unpenalized = TRUE
)

# Use Simulated Data to Run Function --------------------------------------

# Set initial parameters
outcome_mx = sim_data_soap$outcome;
basis_mx = cranio_soap$X;
demo_mx = demo_ingredients_test$X;
pen_basis = cranio_soap$S;
pen_demo = demo_penalty_list_test;
basis_has_unpenalized = FALSE;
demo_has_unpenalized = TRUE;
gamma0 = 0;
beta0 = 0;
sigmasq0 = 1;
tausq0 = 5;
lambda_basis0 = 0.1;
lambda_demo0 = 0.01;
a_sigma_pri = 1;
b_sigma_pri = 1;
a_gamma_pri = 1;
b_gamma_pri = 1;
a_tau_pri = 1;
b_tau_pri = 1;
a_alpha_pri = 1;
b_alpha_pri = 1;
iters = 1e4;
burn_pct = 0.5

