### Write Gibbs sampler for Bayesian linear regression

# Notes:
# - Assumes that each penalty matrix is the full K x K but is block diagonal
# - Assumes at least one set of penalized terms
# - Assumes one distinct smoothing parameter lambda for each S penalty matrix
# - If there are unpenalized terms, assumes exactly one matrix is unpenalized
# - Assumes that the one unpenalized matrix is the first in the list

hierarchical_penalized_gibbs <- function(
    outcome_mx, # G
    basis_mx, # B
    demo_mx, # W
    pen_basis, # P_gamma
    pen_demo, # P_beta
    basis_has_unpenalized = FALSE, # Controls 1st lambda in basis penalty list
    demo_has_unpenalized = FALSE, # Controls 1st lambda in demo penalty list
    gamma0, beta0, sigmasq0, lambda_basis0, tausq0, lambda_demo0, # Default parameters
    a_sigma_pri, b_sigma_pri, a_gamma_pri, b_gamma_pri,
    a_tau_pri, b_tau_pri, a_alpha_pri, b_alpha_pri, # Default hyperparameters
    iters = 1e4, burn_pct = 0.1 # MCMC parameters
){

  # Set up time tracking
  time_check <- Sys.time()
  time_pts <- c("start" = time_check)
  timer <- c("setup" = 0,
             "penalty_total" = 0,
             "gamma_update_c" = 0,
             "gamma_draw" = 0,
             "beta_update_c" = 0,
             "beta_draw" = 0,
             "sigma2_update" = 0,
             "lambdas_basis_update" = 0,
             "tau2_update" = 0,
             "lambdas_demo_update" = 0)

  ### Reformat penalty matrices
  if(is.matrix(pen_basis)){
    pen_basis <- list(pen_basis) }

  if(is.matrix(pen_list_demo)){
    pen_demo <- list(pen_demo) }

  ### Calculate dimensions

  ## Data matrices: outcome, basis, demographics
  N <- nrow(outcome_mx)
  M <- ncol(outcome_mx)
  K <- ncol(basis_mx)
  Q <- ncol(demographic_mx)

  stopifnot(nrow(basis_mx) == M, nrow(demographic_mx) == N)

  ### Determine what lambdas we need from penalty matrices

  # number of penalties on spatial basis coefficients
  pen_blocks_basis <- length(pen_basis)
  num_penals_basis <- penalty_blocks_basis - as.numeric(basis_has_unpenalized)

  # number of penalties on demographic coefficients
  pen_blocks_demo <- length(pen_demo)
  num_penals_demo <- penalty_blocks_demo - as.numeric(demo_has_unpenalized)

  # non-zero columns in each penalty matrix
  pen_block_cols_basis <- map(pen_basis, ~which(apply(., 2, sum) != 0))
  pen_block_cols_demo <- map(pen_demo, ~which(apply(., 2, sum) != 0))

  # number of parameters for each penalty
  pen_block_dims_basis <- map_dbl(pen_block_cols_basis, length)
  pen_block_dims_demo <- map_dbl(pen_block_cols_demo, length)

  # First number of parameters might not impact any lambda; if so, remove it
  if(basis_has_unpenalized){ pen_block_dims_basis <- pen_block_dims_basis[-1] }
  if(demo_has_unpenalized){ pen_block_dims_demo <- pen_block_dims_demo[-1] }

  ### Data manipulations

  # Convert outcome matrix to vector
  outcome_vec <- outcome_mx
  dim(outcome_vec) <- c(prod(dim(outcome_vec)), 1)

  # Calculate sum of outcome squared (G transpose * G)
  outcome_vec_sq <- sum(outcome_vec^2)

  # Calculate basis matrix squared


  # Calculate Demographic matrix squared


  # Calculate basis kron I_N and basis transpose kron I_N


  # Calculate I_K kron demographic


  ### Initialize hyperparameters

  # Assign initial values
  gamma <- rep_len(gamma0, N*K)
  beta <- rep_len(beta0, Q*K)

  sigmasq <- sigmasq0
  lambda_basis <- rep_len(lambda_basis0, num_penals_basis)

  tausq <- tausq0
  lambda_demo <- rep_len(lambda_demo0, num_penals_demo)

  # Initialize gamma and beta moments before update
  gamma_expect <- gamma
  gamma_var <- diag(N*K)
  basis_sq_plus_pen <- diag(N*K)

  beta_expect <- beta
  beta_var <- diag(Q*K)

  # Compute hyperparameters that don't update (gamma / Igamma shapes)
  a_sigma_post <- a_sigma_pri + N*M/2 + N*K/2
  a_gamma_post <- a_gamma_pri + N*pen_block_dims_basis/2
  a_tau_post <- a_tau_pri + Q*K/2
  a_alpha_post <- a_alpha_pri + K*pen_block_dims_demo/2

  ### Create objects to save output

  ## Set iterations of burn-in
  iters_burn <- ceiling(iters * burn_pct)

  # Posterior samples
  gamma_out <- array(0, dim = c(iters - iters_burn, N, K))
  beta_out <- array(0, dim = c(iters - iters_burn, Q, K))
  sigmasq_out <- array(0, dim = c(iters - iters_burn, 1))
  lambda_basis_out <- array(0, dim = c(iters - iters_burn, num_penals_basis))
  tausq_out <- array(0, dim = c(iters - iters_burn, 1))
  lambda_demo_out <- array(0, dim = c(iters - iters_burn, num_penals_demo))

  # Burn-in period
  gamma_burn <- array(0, dim = c(iters_burn, N, K))
  beta_burn <- array(0, dim = c(iters_burn, Q, K))
  sigmasq_burn <- array(0, dim = c(iters_burn, 1))
  lambda_basis_burn <- array(0, dim = c(iters_burn, num_penals_basis))
  tausq_burn <- array(0, dim = c(iters_burn, 1))
  lambda_demo_burn <- array(0, dim = c(iters_burn, num_penals_demo))

  timer["setup"] <- timer["setup"] +
    difftime(Sys.time(), time_check, units = "secs")
  time_check <- Sys.time()

  ### Run MCMC loop

  for(iter in 1:iters){

    # Output on checkpoints
    if(iter %% iters_burn == 0) cat("at turn", iter, "\n")

    ### Compile "Total" penalty matrices (including lambdas)

    ## First, total penalty matrix for gamma coefficients

    # If there are unpenalized terms, do not multiply 1st matrix by any lambda
    lambda_basis_mult <- if(basis_has_unpenalized){
      c(1, lambda_basis)
    } else { lambda_basis }

    lambda_demo_mult <- if(demo_has_unpenalized){
      c(1, lambda_demo)
    } else { lambda_demo }

    # Multiply block penalty matrices by lambdas
    pen_basis_lambda <- map2(lambda_basis_mult, pen_basis, ~ .x * .y)

    pen_demo_lambda <- map2(lambda_demo_mult, pen_demo, ~ .x * .y)

    # Sum together these individual (lambda * S) block diagonal matrices
    pen_basis_total <- Reduce("+", pen_basis_lambda)

    pen_demo_total <- Reduce("+", pen_demo_lambda)

    timer["penalty_total"] <- timer["penalty_total"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    ## (1) Update gamma (spatial coefficients)

    # ** [[Write Gamma moments update function]] **

    timer["gamma_update_c"] <- timer["gamma_update_c"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    # ** [[Write Gamma draw function]] **

    timer["gamma_draw"] <- timer["gamma_draw"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    ## (2) Update beta (demographic coefficients)

    # Call to Cpp function for beta update

    ## ** [[Write new C++ beta moments update function]]

    timer["beta_update_c"] <- timer["beta_update_c"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    # Generate new beta sample
    beta <- t(rmvnorm(n = 1, beta_expect, beta_var))

    timer["beta_draw"] <- timer["beta_draw"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    ## (2) Update sigma^2 (residual variance)

    # 2nd parameter of inverse gamma
    b_sigma_post <-
      b_pri - crossprod(y, X) %*% beta +
      0.5 * (sum(y^2) + crossprod(beta, xtx_plus_pen) %*% beta)

    # Generate new sigma squared sample
    precision <- rgamma(1, a_post, b_post)
    sigmasq <- 1 / precision

    timer["sigma2_update"] <- timer["sigma2_update"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    ## (3) Update lambda (beta function smoothness factor)

    # Calculate posterior gamma shape
    d_post <- map(
      # Exclude any unpenalized terms from posterior
      if(has_unpenalized){ S[-1] } else { S },
      function(S_j){
        (1 / (2 * sigmasq)) * crossprod(beta, S_j) %*% beta + d_pri
      })

    # Update lambdas
    lambda <- map2_dbl(c_post, d_post,
                       function(c_new, d_new){
                         rgamma(1, c_new, d_new)
                       })

    timer["lambdas_update"] <- timer["lambdas_update"] +
      difftime(Sys.time(), time_check, units = "mins")
    time_check <- Sys.time()

    ## Save results
    iter_post <- iter - iters_burn

    if(iter_post > 0){
      beta_out[iter_post,] <- beta
      sigmasq_out[iter_post,] <- sigmasq
      lambda_out[iter_post,] <- lambda
    } else{
      beta_burn[iter,] <- beta
      sigmasq_burn[iter,] <- sigmasq
      lambda_burn[iter,] <- lambda
    }
  }

  time_pts <- c(time_pts, "end" = Sys.time())
  timer["total"] <- difftime(time_pts["end"], time_pts["start"], units = "secs")
  timer <- enframe(timer, name = "step", value = "secs")

  return(list("beta" = beta_out,
              "sigma_sq" = sigmasq_out,
              "lambda" = lambda_out,
              "burn_beta" = beta_burn,
              "burn_sigma_sq" = sigmasq_burn,
              "burn_lambda" = lambda_burn,
              "start-stop" = time_pts,
              "timing" = timer))
}
