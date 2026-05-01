################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Summarize results from test Bayesian soap film on 1 individual
################################################################################

source(here::here("source", "000_definitions.R"))


# Load files --------------------------------------------------------------

# Load test bayesian soap film ("soap_bayes1")
load(file = here("data", "intermediate", "soap_test_single_obs_bayes.rda"))

# Load the single observation we fit the data on ("data_eg")
load(file = here("data", "intermediate", "soap_test_single_obs_data.rda"))

# Load soap film object ("cranio_soap")
load(file = here("data", "cleaned", "soap_object.rda"))

# Load consistent scale for plotting ("so_diff_range")
load(file = here("data", "intermediate", "soap_test_sample_range.rda"))

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
    x <- as_tibble(x, .name_repair = "unique")
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
                            "Burn-In", "Post-Burn-In"))

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
soap_bayes1_preds <- cranio_soap$X %*% soap_bayes1_betas

# Add predictions to data
data_eg %<>% mutate(fit_bayes = as.numeric(soap_bayes1_preds))

data_toplot <- data_eg %>%
  select(row, col, diff, fit_bayes) %>%
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
