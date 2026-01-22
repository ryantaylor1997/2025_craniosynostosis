
# Run Bayesian soap film on 1 individual -----------------------------

### Run Gibbs sampler for soap film

# Set data to train on, especially outcome
data_eg <- so_test_fit %>%
  filter(fusion_type == "Sagittal") %>%
  filter(fname == first(fname))

y_mx <- data_eg$diff

# Set penalty matrices
S_mx <- cranio_soap$S

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

set.seed(70)

# Run function
soap_bayes1 <- lm_penalized_gibbs(
  y = y_mx, X = so_X, S = S_mx,
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


# Print MCMC diagnostics --------------------------------------------------

# Check timing
mcmc_time <- soap_bayes1$timing %>%
  mutate(pct = as.numeric(secs) / last(secs))

# Make accessible all parameters together
soap_bayes1_tr <- imap(
  soap_bayes1 %>%
    list_modify("timing" = zap(),
                "start-stop" = zap()),
  function(x, xn){
    x <- as_tibble(x)
    if(ncol(x) == 1){
      names(x) <- xn
    } else {
      names(x) <- paste0(xn, ".", 1:ncol(x))
    }
    x %<>%
      mutate(iter = 1:n(), .before = everything()) %>%
      pivot_longer(-iter,
                   names_pattern = "(burn_)?(.*)",
                   names_to = c("draw", "name"))

    x
  }
) %>%
  unname() %>%
  list_rbind() %>%
  mutate(draw_cat = if_else(draw == "burn_",
                            "Burn-In", "Posterior"))

# Add default parameters to this table
trace_default <- soap_bayes1_tr %>%
  filter(draw_cat == "Burn-In") %>%
  distinct(draw_cat, name) %>%
  filter(!str_detect(name, "beta")) %>%
  mutate(
    iter = 0,
    value = case_when(
      str_detect(name, "sigma_sq") ~ def_pars["sigma_sq"],
      str_detect(name, "lambda") ~ def_pars["lambda"],
      T ~ NA
    ))

soap_bayes1_tr %<>% bind_rows(trace_default)

# Trace plot of sigma
sigma_trace <- ggplot(soap_bayes1_tr %>%
                        filter(str_detect(name, "sigma_sq"))) +
  geom_line(aes(x = iter, y = value)) +
  facet_wrap(~ name + draw_cat, ncol = 2, scales = "free",
             labeller = \(x) label_value(x, multi_line = F))

# Trace plot of lambdas
lambda_trace <- ggplot(soap_bayes1_tr %>%
                         filter(str_detect(name, "lambda"))) +
  geom_line(aes(x = iter, y = value), alpha = 0.5) +
  facet_wrap(~ name + draw_cat, ncol = 2, scales = "free",
             labeller = \(x) label_value(x, multi_line = F)) +
  theme(legend.position = "bottom")

# Trace plot of betas
beta_trace <- ggplot(soap_bayes1_tr %>%
                       filter(str_detect(name, "beta"))) +
  geom_line(aes(x = iter, y = value), alpha = 0.5) +
  theme(legend.position = "bottom") +
  facet_wrap(~ name + draw_cat, ncol = 2, scales = "free",
             labeller = \(x) label_value(x, multi_line = F))


# Visualize Bayesian soap film predictions --------------------------------

# Extract posterior mean coefficients
soap_bayes1_betas <- colMeans(soap_bayes1$beta)

# Calculate predictions from these coefficients
soap_bayes1_preds <- so_X %*% soap_bayes1_betas

# Add predictions to data
data_eg %<>% mutate(fit_bayes = as.numeric(soap_bayes1_preds))

data_toplot <- data_eg %>%
  select(row, col, diff, fit_gcv = fit, fit_bayes) %>%
  pivot_longer(-c(row, col))

# Plot truth and predictions from Bayesian example
soap_bayes1_plot <- ggplot(data_toplot) +
  geom_raster(aes(x = row, y = col, fill = value)) +
  scale_fill_viridis_c(option = "turbo", limits = so_diff_range) +
  facet_wrap(~name, ncol = 1) +
  labs(x = NULL, y = NULL,
       fill = "Growth",
       title = "Bayesian Prediction Comparison") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom")
