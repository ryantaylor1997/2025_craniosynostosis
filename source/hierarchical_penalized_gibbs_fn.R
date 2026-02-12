
# Write Gibbs sampler for Bayesian linear regression ----------------------

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
    gamma0, beta0, sigmasq0, tausq0, lambda_basis0, lambda_demo0, # Default parameters
    a_sigma_pri, b_sigma_pri, a_gamma_pri, b_gamma_pri,
    a_tau_pri, b_tau_pri, a_alpha_pri, b_alpha_pri, # Default hyperparameters
    iters = 1e4, burn_pct = 0.1 # MCMC parameters
){

  ### Set up time tracking
  time_check <- Sys.time()
  time_pts <- c("start" = time_check)
  timer <- c("setup" = 0,
             "penalty_total" = 0,
             "gamma_update" = 0,
             "gamma_draw" = 0,
             "beta_update" = 0,
             "beta_draw" = 0,
             "sigma2_update" = 0,
             "tau2_update" = 0,
             "lambdas_basis_update" = 0,
             "lambdas_demo_update" = 0)

  ### Reformat inputs

  ## Reformat penalty matrices as lists
  if(is.matrix(pen_basis)){
    pen_basis <- list(pen_basis) }

  if(is.matrix(pen_demo)){
    pen_demo <- list(pen_demo) }

  ## Reformat data matrices
  outcome_mx <- Matrix(outcome_mx)
  basis_mx <- Matrix(basis_mx)
  demo_mx <- Matrix(demo_mx)

  ### Calculate dimensions

  ## Dimensions of data matrices: outcome, basis, demographics
  N <- nrow(outcome_mx)
  M <- ncol(outcome_mx)
  K <- ncol(basis_mx)
  Q <- ncol(demo_mx)

  stopifnot(nrow(basis_mx) == M, nrow(demo_mx) == N)

  ### Determine what we need from penalty matrices

  # number of penalties on spatial basis coefficients
  num_pen_blocks_basis <- length(pen_basis)
  num_penals_basis <- num_pen_blocks_basis - as.numeric(basis_has_unpenalized)

  # number of penalties on demographic coefficients
  num_pen_blocks_demo <- length(pen_demo)
  num_penals_demo <- num_pen_blocks_demo - as.numeric(demo_has_unpenalized)

  # non-zero columns in each penalty matrix
  pen_block_cols_basis <- map(pen_basis, ~which(apply(., 2, sum) != 0))
  pen_block_cols_demo <- map(pen_demo, ~which(apply(., 2, sum) != 0))

  # number of parameters for each penalty
  pen_block_dims_basis <- map_dbl(pen_block_cols_basis, length)
  pen_block_dims_demo <- map_dbl(pen_block_cols_demo, length)

  # Separate block-diagonal penalty matrices into their specific sections
  pen_block_list_basis <- map2(pen_basis, pen_block_cols_basis,
                               ~ .x[.y, .y])

  pen_block_list_demo <- map2(pen_demo, pen_block_cols_demo,
                              ~ .x[.y, .y])

  # If first block of parameters is not included in lambda update,
  # remove it from dimensions, column indices, and one list of penalty matrices
  if(basis_has_unpenalized){
    pen_block_dims_basis <- pen_block_dims_basis[-1]
    pen_block_cols_basis <- pen_block_cols_basis[-1]
    pen_block_list_basis <- pen_block_list_basis[-1]
  }

  if(demo_has_unpenalized){
    pen_block_dims_demo <- pen_block_dims_demo[-1]
    pen_block_cols_demo <- pen_block_cols_demo[-1]
    pen_block_list_demo <- pen_block_list_demo[-1]
  }

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
  gamma_var <- Diagonal(N*K)
  basis_sq_plus_pen <- Diagonal(N*K)

  beta_expect <- beta
  beta_var <- Diagonal(Q*K)

  # Compute hyperparameters that don't update (gamma / Igamma shapes)
  a_sigma_post <- a_sigma_pri + N*M/2 + N*K/2
  a_gamma_post <- a_gamma_pri + N*pen_block_dims_basis/2
  a_tau_post <- a_tau_pri + Q*K/2
  a_alpha_post <- a_alpha_pri + K*pen_block_dims_demo/2

  ### Data manipulations

  # Allocate storage for versions of gamma and beta in matrix form
  gamma_mx <- Matrix(gamma, N, K)
  beta_mx <- Matrix(beta, Q, K)

  # Allocate storage for subsets of these coefficients in lambda updates
  gamma_mx_list <- map(pen_block_cols_basis, ~gamma_mx[, .x])
  beta_mx_basis_list <- map(pen_block_cols_basis, ~beta_mx[, .x])
  beta_mx_demo_list <- map(pen_block_cols_demo, ~beta_mx[.x, ])

  # Convert outcome matrix to vector
  outcome_vec <- as(outcome_mx, "sparseVector")

  # Calculate sum of outcome squared (G transpose * G)
  outcome_vec_sq <- sum(outcome_vec^2)

  # Calculate basis matrix squared
  basis_sq <- crossprod(basis_mx)

  # Calculate Demographic matrix squared
  demo_sq <- crossprod(demo_mx)

  # Calculate basis kron I_N and basis transpose kron I_N
  basis_kron <- kronecker(basis_mx, Diagonal(N))
  basis_t_kron <- kronecker(t(basis_mx), Diagonal(N))

  # Calculate I_K kron demographic
  demo_kron <- kronecker(Diagonal(K), demo_mx)

  ### Create objects to save output

  ## Set iterations of burn-in
  iters_burn <- ceiling(iters * burn_pct)

  # Posterior samples
  gamma_out <- array(0, dim = c(N, K, iters - iters_burn))
  beta_out <- array(0, dim = c(Q, K, iters - iters_burn))
  sigmasq_out <- array(0, dim = c(1, iters - iters_burn))
  lambda_basis_out <- array(0, dim = c(num_penals_basis, iters - iters_burn))
  tausq_out <- array(0, dim = c(1, iters - iters_burn))
  lambda_demo_out <- array(0, dim = c(num_penals_demo, iters - iters_burn))

  # Burn-in period
  gamma_burn <- array(0, dim = c(N, K, iters_burn))
  beta_burn <- array(0, dim = c(Q, K, iters_burn))
  sigmasq_burn <- array(0, dim = c(1, iters_burn))
  lambda_basis_burn <- array(0, dim = c(num_penals_basis, iters_burn))
  tausq_burn <- array(0, dim = c(1, iters_burn))
  lambda_demo_burn <- array(0, dim = c(num_penals_demo, iters_burn))

  timer["setup"] <- timer["setup"] +
    difftime(Sys.time(), time_check, units = "secs")
  time_check <- Sys.time()

  ### Run MCMC loop

  cat("MCMC Status: \n")

  for(iter in 1:iters){

    # Output on checkpoints
    if(iter %% iters_burn == 0) cat("at turn", iter, "\n")

    ### Compile "Total" penalty matrices (including lambdas)

    ## If there are unpenalized terms, do not multiply 1st matrix by any lambda
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

    pen_basis_total <- Matrix(pen_basis_total)

    pen_demo_total <- Matrix(pen_demo_total)

    timer["penalty_total"] <- timer["penalty_total"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    ### (1) Update gamma (spatial coefficients)

    ## Compute matrices that we use in multiple places
    basis_sq_plus_pen <- basis_sq + pen_basis_total

    basis_sq_plus_pen_inv <- solve(basis_sq_plus_pen)

    basis_sq_plus_pen_inv_kron <- kronecker(basis_sq_plus_pen_inv, Diagonal(N))

    pen_basis_kron_demo <- kronecker(pen_basis_total, demo_mx)

    ## Compute expectation and variance
    gamma_expect <- basis_sq_plus_pen_inv_kron %*%
      (basis_t_kron %*% outcome_vec + pen_basis_kron_demo %*% beta)

    gamma_var <- sigmasq * basis_sq_plus_pen_inv_kron

    timer["gamma_update"] <- timer["gamma_update"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    ## Gamma draw

    # Use custom function to perform Cholesky decomposition on sparse Matrix
    gamma <- t(rmvnorm_Matrix(n = 1, gamma_expect, gamma_var))

    gamma_mx <- Matrix(as(gamma, "sparseVector"), N, K, byrow = F)

    timer["gamma_draw"] <- timer["gamma_draw"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    ### (2) Update beta (demographic coefficients)

    ## Compute matrices based on penalty matrices that we use in multiple places
    pen_basis_kron_demo_sq <- kronecker(pen_basis_total, demo_sq)

    pen_demo_kron <- kronecker(Diagonal(K), pen_demo_total)

    penalties_kron <- pen_basis_kron_demo_sq +
      (sigmasq / tausq) * pen_demo_kron

    pen_basis_kron_demo_t <- kronecker(pen_basis_total, t(demo_mx))

    # Take inverse that comprises most of the variance
    # -> tol = 0 avoids testing for whether it is near-singular
    penalties_kron_inv <- solve(penalties_kron)

    ## Compute expectation and variance
    beta_expect <- penalties_kron_inv %*%
      pen_basis_kron_demo_t %*% gamma

    beta_var <- sigmasq * penalties_kron_inv

    timer["beta_update"] <- timer["beta_update_c"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    ## Generate new beta sample
    beta <- t(rmvnorm_Matrix(n = 1, beta_expect, beta_var))

    beta_mx <- matrix(beta, Q, K, byrow = F)

    timer["beta_draw"] <- timer["beta_draw"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    ### (3) Update sigma^2 (residual variance)

    # 2nd parameter of inverse gamma
    b_sigma_post <-
      b_sigma_pri + 0.5 * (
        outcome_vec_sq +
          crossprod(gamma, kronecker(basis_sq_plus_pen, Diagonal(N)) ) %*% gamma +
          crossprod(beta, pen_basis_kron_demo_sq) %*% beta
      ) -
      crossprod(gamma, basis_t_kron) %*% outcome_vec -
      crossprod(beta, pen_basis_kron_demo_t) %*% gamma

    b_sigma_post <- b_sigma_post[1,1]

    # Generate new sigma squared sample
    precision_data <- rgamma(1, a_sigma_post, rate = b_sigma_post)
    sigmasq <- 1 / precision_data

    timer["sigma2_update"] <- timer["sigma2_update"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    ### (4) Update tau^2 (demographic effect scale)

    # Second parameter of inverse gamma
    b_tau_post <- b_tau_pri +
      0.5 * crossprod(beta, kronecker(Diagonal(K), pen_demo_total) ) %*% beta

    b_tau_post <- b_tau_post[1,1]

    # Draw new tau squared sample
    precision_param <- rgamma(1, a_tau_post, rate = b_tau_post)
    tausq <- 1 / precision_param

    timer["tau2_update"] <- timer["tau2_update"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    ### (5) Update lambdas for soap film (spatial smoothness factors)

    # Update lists of parameter matrices separated into soap film-based blocks
    gamma_mx_list <- map(pen_block_cols_basis, ~gamma_mx[, .x])
    beta_mx_basis_list <- map(pen_block_cols_basis, ~beta_mx[, .x])

    # Calculate posterior gamma shape
    b_gamma_post <- pmap(
      list(
        b_gamma_pri,
        gamma_mx_list,
        pen_block_dims_basis,
        beta_mx_basis_list,
        pen_block_list_basis
      ),
      function(b_pri, gamma_block, pen_block_dim, beta_block, pen_block){

        # Vectorize gamma and beta block matrices
        dim(gamma_block) <- c(prod(dim(gamma_block)), 1)
        dim(beta_block) <- c(prod(dim(beta_block)), 1)

        # Isolate demographic Kronecker product the size of the block
        demo_kron_block <- demo_kron[1:N*pen_block_dim, 1:Q*pen_block_dim]

        # Subtract chunk of demographic effects from chunk of gammas
        block_vec <- gamma_block - demo_kron_block %*% beta_block

        b_post <- b_pri +
        (1 / (2 * sigmasq)) *
          crossprod(block_vec, kronecker(pen_block, Diagonal(N)) ) %*% block_vec

        return(b_post)
      })

    # Draw new basis lambdas
    lambda_basis <- map2_dbl(
      a_gamma_post, b_gamma_post,
      function(a_new, b_new){
        rgamma(1, a_new, rate = b_new)
      })

    timer["lambdas_basis_update"] <- timer["lambdas_basis_update"] +
      difftime(Sys.time(), time_check, units = "mins")
    time_check <- Sys.time()

    ### (6) Update lambda(s) for demographic effects (demographic smoothness factors)

    # Update list of demographic effects in penalized and unpenalized blocks
    beta_mx_demo_list <- map(pen_block_cols_demo, ~beta_mx[.x, ])

    # Calculate posterior gamma shape
    b_alpha_post <- pmap(
      list(
        b_alpha_pri,
        beta_mx_demo_list,
        pen_block_dims_demo,
        pen_block_list_demo
      ),
      function(b_pri, beta_block, pen_block_dim, pen_block){

        # Vectorize beta block tranpose matrices
        beta_block_t <- t(beta_block)
        dim(beta_block_t) <- c(prod(dim(beta_block_t)), 1)

        b_post <- b_pri +
          (1 / (2 * tausq)) *
          crossprod(beta_block_t, kronecker(pen_block, Diagonal(K)) ) %*% beta_block_t

        return(b_post)
      })

    # Draw new basis lambdas
    lambda_demo <- map2_dbl(
      a_alpha_post, b_alpha_post,
      function(a_new, b_new){
        rgamma(1, a_new, rate = b_new)
      })

    timer["lambdas_demo_update"] <- timer["lambdas_demo_update"] +
      difftime(Sys.time(), time_check, units = "mins")
    time_check <- Sys.time()

    ## Save results
    iter_post <- iter - iters_burn

    if(iter_post > 0){
      gamma_out[, , iter_post] <- gamma_mx
      beta_out[, , iter_post] <- beta_mx
      sigmasq_out[, iter_post] <- sigmasq
      tausq_out[, iter_post]
      lambda_basis_out[, iter_post] <- lambda_basis
      lambda_demo_out[, iter_post] <- lambda_demo
    } else{
      gamma_burn[, , iter] <- gamma_mx
      beta_burn[, , iter] <- beta_mx
      sigmasq_burn[, iter] <- sigmasq
      tausq_burn[, iter] <- tausq
      lambda_basis_burn[, iter] <- lambda_basis
      lambda_demo_burn[, iter] <- lambda_demo
    }
  }

  time_pts <- c(time_pts, "end" = Sys.time())
  timer["total"] <- difftime(time_pts["end"], time_pts["start"], units = "secs")
  timer <- enframe(timer, name = "step", value = "secs")

  return(list("gamma" = gamma_out,
              "beta" = beta_out,
              "sigma_sq" = sigmasq_out,
              "tau_sq" = tausq_out,
              "lambda_basis" = lambda_basis_out,
              "lambda_demo" = lambda_demo_out,
              "burn_gamma" = gamma_burn,
              "burn_beta" = beta_burn,
              "burn_sigma_sq" = sigmasq_burn,
              "burn_tau_sq" = tausq_burn,
              "burn_lambda_basis" = lambda_basis_burn,
              "burn_lambda_demo" = lambda_demo_burn,
              "start-stop" = time_pts,
              "timing" = timer))
}
