################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Write functions to simulate data for hierarchical model
################################################################################

### Simulate demographic data
# Note: fusion types are currently hard-coded
make_sim_data_obs_level <- function(n_obs, age_max, sex_prop_1){
  tibble(
    sex = rbinom(n_obs, 1, sex_prop_1),
    age = runif(n_obs, max = age_max),
    fusion_type = c(rep("Normative", n_obs / 2),
                    rep("Fusion_1", n_obs / 4),
                    rep("Fusion_2", n_obs / 4))
  ) %>%
    mutate(fusion_type = fct_inorder(fusion_type),
           sex = as.factor(sex))
}

### Simulate data according to simple hierarchical data-generating model
# Note:
make_sim_data_hier_soap <- function(soap_design_mx, obs_design_mx,
                                    penalty_soap_list, penalty_obs_list,
                                    sigmasq, tausq,
                                    lambdas_soap, lambdas_obs){

  ### Save dimensions
  N <- nrow(obs_design_mx)
  Q <- ncol(obs_design_mx)
  M <- nrow(soap_design_mx)
  K <- ncol(soap_design_mx)

  ### Combine parameters with given matrices for Total Penalty Matrices

  # Multiply lambdas by demographic penalty matrix
  penalty_obs_total_list <- map2(lambdas_obs, penalty_obs_list,
                                  ~.x * .y)

  penalty_soap_total_list <- map2(lambdas_soap, penalty_soap_list,
                                   ~.x * .y)

  # Sum block diagonal matrices together
  penalty_obs_total_mx <- Reduce("+", penalty_obs_total_list)

  penalty_soap_total_mx <- Reduce("+", penalty_soap_total_list)

  ### Moments for generated parameters

  # Take generalized inverse of total penalty matrices
  penalty_obs_var <- ginv(as.matrix(penalty_obs_total_mx))

  penalty_soap_var <- ginv(as.matrix(penalty_soap_total_mx))

  # Convert to efficient Matrix format
  penalty_obs_var <- Matrix(penalty_obs_var)

  penalty_soap_var <- Matrix(penalty_soap_var)

  # Scale inverse penalty matrices by scalar variance terms
  var_obs_tau <- tausq * penalty_obs_var

  var_soap_sigma <- sigmasq * penalty_soap_var

  ### Draw first-level parameters

  # Draw betas as K iid Q-length vectors
  beta_mx <- Matrix(t(
    rmvnorm(K,
            mean = rep(0, Q),
            sigma = as.matrix(var_obs_tau),
            method = "svd")
  ))

  # Draw epsilon-gammas as N iid K-length vectors
  epsilon_gamma <- rmvnorm(N,
                           mean = rep(0, K),
                           sigma = as.matrix(var_soap_sigma),
                           method = "svd")

  # Draw epsilon-Gs as N*M iid scalars
  epsilon_G <- rnorm(N*M, 0, sigmasq)

  epsilon_G_mx <- Matrix(epsilon_G, N, M)

  ### Calculate 2nd-level parameters

  # Calculate soap film coefficients
  gamma_mx_true <- obs_design_mx %*% beta_mx

  gamma_mx <- gamma_mx_true + epsilon_gamma

  # Calculate outcome data
  outcome_mx_noerror <- tcrossprod(gamma_mx_true, soap_design_mx)
  outcome_mx_nonrandom <- tcrossprod(gamma_mx, soap_design_mx)
  outcome_mx_sim <- outcome_mx_nonrandom + epsilon_G_mx

  ### Export
  out_list <- list(
    "outcome" = outcome_mx_sim,
    "outcome_nonrandom" = outcome_mx_nonrandom,
    "outcome_noerror" = outcome_mx_noerror,
    "gamma" = gamma_mx,
    "gamma_true" = gamma_mx_true,
    "beta" = beta_mx,
    "beta_var" = var_obs_tau,
    "gamma_var" = var_soap_sigma,
    "call" = list(
      "soap_design_mx" = soap_design_mx,
      "obs_design_mx" = obs_design_mx,
      "penalty_soap_list" = penalty_soap_list,
      "penalty_obs_list" = penalty_obs_list,
      "sigmasq" = sigmasq,
      "tausq" = tausq,
      "lambdas_soap" = lambdas_soap,
      "lambdas_obs" = lambdas_obs
    )
  )

  return(out_list)
}


### Simulate data with strictly positive outcomes
# Note:
make_sim_data_pos_soap <- function(soap_design_mx,
                                   obs_data_df,
                                   penalty_soap_list,
                                   sigmasq,
                                   lambdas_soap){

  ### Save dimensions
  M <- nrow(soap_design_mx)
  K <- ncol(soap_design_mx)

  ### Combine parameters with given matrices for Total Penalty Matrices

  # Multiply lambdas by demographic penalty matrix
  penalty_soap_total_list <- map2(lambdas_soap, penalty_soap_list,
                                  ~.x * .y)

  # Sum block diagonal matrices together
  penalty_soap_total_mx <- Reduce("+", penalty_soap_total_list)

  ### Moments for generated parameters

  # Take generalized inverse of total penalty matrices
  penalty_soap_var <- ginv(as.matrix(penalty_soap_total_mx))

  # Convert to efficient Matrix format
  penalty_soap_var <- Matrix(penalty_soap_var)

  # Scale inverse penalty matrices by scalar variance terms
  var_soap_sigma <- sigmasq * penalty_soap_var

  ### Draw first-level parameters

  # Draw betas as K iid Q-length vectors
  beta_mx <- Matrix(t(
    rmvnorm(K,
            mean = rep(0, Q),
            sigma = as.matrix(var_obs_tau),
            method = "svd")
  ))

  # Draw epsilon-gammas as N iid K-length vectors
  epsilon_gamma <- rmvnorm(N,
                           mean = rep(0, K),
                           sigma = as.matrix(var_soap_sigma),
                           method = "svd")

  # Draw epsilon-Gs as N*M iid scalars
  epsilon_G <- rnorm(N*M, 0, sigmasq)

  epsilon_G_mx <- Matrix(epsilon_G, N, M)

  ### Calculate 2nd-level parameters

  # Calculate soap film coefficients
  gamma_mx_true <- #*UPDATE*

  gamma_mx <- gamma_mx_true + epsilon_gamma

  # Calculate outcome data
  outcome_mx_noerror <- tcrossprod(gamma_mx_true, soap_design_mx)
  outcome_mx_nonrandom <- tcrossprod(gamma_mx, soap_design_mx)
  outcome_mx_sim <- outcome_mx_nonrandom + epsilon_G_mx

  ### Export
  out_list <- list(
    "outcome" = outcome_mx_sim,
    "outcome_nonrandom" = outcome_mx_nonrandom,
    "outcome_noerror" = outcome_mx_noerror,
    "gamma" = gamma_mx,
    "gamma_true" = gamma_mx_true,
    "gamma_var" = var_soap_sigma,
    "call" = list(
      "soap_design_mx" = soap_design_mx,
      "penalty_soap_list" = penalty_soap_list,
      "sigmasq" = sigmasq,
      "lambdas_soap" = lambdas_soap
    )
  )

  return(out_list)
}

