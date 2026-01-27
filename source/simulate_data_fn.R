################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Write function to simulate data for hierarchical model
################################################################################

make_sim_data_bayes_soap <- function(basis_mx, demo_mx,
                                     penalty_basis_list, penalty_demo_list,
                                     sigmasq, tausq,
                                     lambdas_basis, lambdas_demo,
                                     demo_has_unpenalized = T){

  ### Save dimensions
  N <- nrow(demo_mx)
  M <- nrow(basis_mx)
  K <- ncol(basis_mx)
  Q <- ncol(demo_mx)

  ### Combine parameters with given matrices

  # Add 1 to multiply instead of lambda for unpenalized demographics
  lambda_mult_demo <- if(demo_has_unpenalized){
    c(1, lambdas_demo)
  } else {
    lambdas_demo
  }

  ### Total Penalty Matrices

  # Multiply lambdas by demographic penalty matrix
  penalty_total_demo_list <- map2(lambda_mult_demo, penalty_demo_list,
                                  ~.x * .y)

  penalty_total_basis_list <- map2(lambdas_basis, penalty_basis_list,
                                   ~.x * .y)

  # Sum block diagonal matrices together
  penalty_total_demo_mx <- Reduce("+", penalty_total_demo_list)

  penalty_total_basis_mx <- Reduce("+", penalty_total_basis_list)

  ### Moments for generated parameters

  # Take generalized inverse of total penalty matrices
  penalty_var_demo <- ginv(as.matrix(penalty_total_demo_mx))

  penalty_var_basis <- ginv(as.matrix(penalty_total_basis_mx))

  # Convert to efficient Matrix format
  penalty_var_demo <- Matrix(penalty_var_demo)

  penalty_var_basis <- Matrix(penalty_var_basis)

  # Scale inverse penalty matrices by scalar variance terms
  var_demo_tau <- tausq * penalty_var_demo

  var_basis_sigma <- sigmasq * penalty_var_basis

  ### Draw first-level parameters

  # Draw betas as K iid Q-length vectors
  demo_effect_mx <- Matrix(t(
    rmvnorm(K, mean = rep(0, Q), sigma = as.matrix(var_demo_tau),
            method = "svd")
  ))

  # Draw epsilon-gammas as N iid K-length vectors
  epsilon_gamma <- rmvnorm(N, mean = rep(0, K), sigma = as.matrix(var_basis_sigma),
                           method = "svd")

  # Draw epsilon-Gs as N*M iid scalars
  epsilon_G <- rnorm(N*M, 0, sigmasq)

  epsilon_G_mx <- Matrix(epsilon_G, N, M)

  ### Calculate 2nd-level parameters

  # Calculate soap film coefficients
  gamma_mx <- demo_mx %*% demo_effect_mx + epsilon_gamma

  # Calculate outcome data
  outcome_mx_sim <- tcrossprod(gamma_mx, basis_mx) + epsilon_G_mx

  ### Export
  out_list <- list(
    "outcome" = outcome_mx_sim,
    "gamma" = gamma_mx,
    "beta" = demo_effect_mx
  )

  return(out_list)
}
