################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Run Simulations for Hierarchical Bayes Model
################################################################################


# Take subset of data for demographics ------------------------------------

# Identify patients of each fusion type
set.seed(978)

cranio_demo_test_df <- cranio_sub %>%
  filter(!fname %in% cranio_dup_fnames) %>%
  group_by(fusion_type) %>%
  slice_sample(n = 10)


# Convert demographic data to model matrix --------------------------------

cranio_demo_smooth_test <- smooth.construct2(
  s(age, bs = "cr", k = cranio_knots_age),
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

### Set scalar parameters
sigma2_test <- 2

tau2_test <- 3

lambda_basis_vec <- c(0.001, 0.01)

lambda_demo <- 0.005

### Generate simulated data from these assumptions
sim_data_soap <- make_sim_data_bayes_soap(
  basis_mx = so_X,
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
basis_mx = so_X;
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

