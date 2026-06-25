




demo_s_obj <- s(age, bs = "tp", k = demo_age_knots)

demo_smooth_noshrink <- smooth.construct2(
  demo_s_obj, data = demo_test_df, knots = NULL
)

demo_s_diff_plot <- (demo_smooth_test$S[[1]] - demo_smooth_noshrink$S[[1]]) %>%
  reshape2::melt() %>%
  rename(row = Var1, col = Var2) %>%
  mutate(value = value * scale_demo_pen) %>%
  ggplot(aes(x = col, y = -row, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(breaks = scales::breaks_pretty(3)) +
  scale_x_continuous(breaks = scales::breaks_width(1)) +
  scale_y_continuous(breaks = scales::breaks_width(1)) +
  theme(legend.position = "bottom") +
  labs(title = "Demographics Penalty Shrinkage Difference")

demo_s_diff_plot


# Test penalty matrices rank ----------------------------------------------

# Steps to create total basis penalty matrix
num_pen_blocks_basis <- length(pen_basis)
num_penals_basis <- num_pen_blocks_basis - as.numeric(basis_has_unpenalized)

# non-zero columns in each penalty matrix
pen_block_cols_basis <- map(pen_basis, ~which(apply(., 2, sum) != 0))

# number of parameters for each penalty
pen_block_dims_basis <- map_dbl(pen_block_cols_basis, length)

lambda_basis <- rep_len(lambda_basis0, num_penals_basis)

lambda_basis_mult <- if(basis_has_unpenalized){
  c(1, lambda_basis)
} else { lambda_basis }

pen_basis_lambda <- map2(lambda_basis_mult, pen_basis, ~ .x * .y)

pen_basis_total <- Reduce("+", pen_basis_lambda)

# Steps to create total demographics penalty matrix
num_pen_blocks_demo <- length(pen_demo)
num_penals_demo <- num_pen_blocks_demo - as.numeric(demo_has_unpenalized)

lambda_demo <- rep_len(lambda_demo0, num_penals_demo)

lambda_demo_mult <- if(demo_has_unpenalized){
  c(1, lambda_demo)
} else { lambda_demo }

pen_demo_lambda <- map2(lambda_demo_mult, pen_demo, ~ .x * .y)

pen_demo_total <- Reduce("+", pen_demo_lambda)

# Calculation of demographics cross product
demo_sq <- crossprod(demo_mx)

# First matrix in the posterior variance
pen_basis_kron_demo_sq <- kronecker(pen_basis_total, demo_sq)

pen_demo_kron <- kronecker(Diagonal(K), pen_demo_total)

# Full posterior precision
penalties_kron <- pen_basis_kron_demo_sq + (sigmasq / tausq) * pen_demo_kron

# Determine dimension of posterior precision
post_prec_dim <- dim(penalties_kron)

# Rank of posterior precision
post_prec_rank <- qr2rankMatrix(qr(penalties_kron))

cat("Posterior Precision: Dimension - QR Rank =",
    post_prec_dim[2] - post_prec_rank, "\n")

cat("Is Invertible:",
    rcond(penalties_kron) >= .Machine$double.eps,
    "with rcond =", rcond(penalties_kron), "\n")

# Data piece of precision
post_data_dim <- dim(pen_basis_kron_demo_sq)

post_data_rank <- qr2rankMatrix(qr(pen_basis_kron_demo_sq))

cat("Posterior Data Kronecker: Dimension - QR Rank =",
    post_data_dim[2] - post_data_rank)

# Calculate dimensions
N <- nrow(outcome_mx)
M <- ncol(outcome_mx)
K <- ncol(basis_mx)
Q <- ncol(demo_mx)

# Check 2 diagonal blocks of full matrix

index <- 0

diag_block_df <- tibble(block_num = rep(NA, num_penals_basis),
                        block = rep(NA, num_penals_basis),
                        dim = rep(NA, num_penals_basis),
                        rank = rep(NA, num_penals_basis))

for(j in 1:num_penals_basis){

  first_index <- index + 1

  index <- index + pen_block_dims_basis[j] * Q

  block_j <- pen_basis_kron_demo_sq[first_index:index, first_index:index]

  row_j <- tibble_row(
    block_num <- j,
    block = list(block_j),
    dim = min(dim(block_j)),
    rank = qr2rankMatrix(qr(block_j))
  )

  diag_block_df[j,] <- row_j

}

diag_block_df %<>% mutate(diff = dim - rank)

kbl(diag_block_df %>% select(-block)) %>%
  kable_classic(full_width = F)

# Check segment corresponding to cyclic spline

cyclic_kron_demo <- diag_block_df$block[[1]] +
  kronecker(Diagonal(pen_block_dims_basis[1]), pen_demo_total)

cyclic_kron_dim <- dim(cyclic_kron_demo)

cyclic_kron_rank <- qr2rankMatrix(qr(cyclic_kron_demo))


cat("Posterior Data Cyclic Plus Demo: Dimension - QR Rank =",
    cyclic_kron_dim[2] - cyclic_kron_rank)
