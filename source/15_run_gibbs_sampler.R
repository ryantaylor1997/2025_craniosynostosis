################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Run joint Gibbs sampler function on data
################################################################################

cranio_growth_mx <- do.call(rbind,
                            hoist(cranio_models,
                                  "data", "diff")$diff)

cranio_growth_mx <- Matrix(cranio_growth_mx)

if(DO_JOINT_FIT){

# Run hierarchical model
cranio_joint <- hierarchical_penalized_gibbs(
  outcome_mx = cranio_growth_mx,
  basis_mx = cranio_soap$X,
  demo_mx = cranio_coeff_ingredients$X,
  pen_basis = cranio_soap$S,
  pen_demo = cranio_coeff_S_combo,
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
  iters = 1000,
  burn_pct = 0.5,
  update_pct = 0.1
)

save(cranio_joint,
     file = here::here("analysis", "intermediate", "joint_mdoel_fit.rda"))

} else {
  load(file = here::here("analysis", "intermediate", "joint_mdoel_fit.rda"))
}
