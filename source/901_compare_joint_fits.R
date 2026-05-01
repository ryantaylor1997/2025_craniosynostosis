################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Compare Joint Model Fits
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load previous runs ("cranio_joint")
load(file = here("analysis", "intermediate",
                       "joint_model_fit_2026-04-21.rda"))

cranio_joint_recent <- cranio_joint

load(file = here("analysis", "intermediate",
                       "joint_model_fit 03.14.rda"))

cranio_joint_prev <- cranio_joint

rm(cranio_joint); gc()

# Load soap film penalty scales ("soap_pen_scale")
load(file = here("data", "intermediate", "soap_object.rda"))

# Load observation-level smooth scale terms ("obs_pen_scale")
load(file = here("data", "intermediate", "obs_level_penalty_scales.rda"))

# Summarize results -------------------------------------------------------

# Collect posterior mean values

model_summ_prev <- map(cranio_joint_prev[c("sigma_sq",
                                            "tau_sq",
                                            "lambda_basis",
                                            "lambda_demo")],
                       rowMeans) %>%
  enframe() %>%
  unnest_longer(value) %>%
  group_by(name) %>%
  mutate(num = row_number(), .after = name) %>%
  ungroup() %>%
  mutate(model = "previous")

model_summ_rec <- map(cranio_joint_recent[c("sigma_sq",
                                            "tau_sq",
                                            "lambda_basis",
                                            "lambda_demo")],
                      rowMeans) %>%
  enframe() %>%
  unnest_longer(value) %>%
  group_by(name) %>%
  mutate(num = row_number(), .after = name) %>%
  ungroup() %>%
  mutate(model = "recent")

# Bring in penalty scales
scales <- c(soap_pen_scale, obs_pen_scale)

model_summ_scale <- model_summ_prev %>%
  filter(str_detect(name, "lambda")) %>%
  mutate(model = "previous_scaled") %>%
  mutate(value = value / scales)

comp_tbl <- bind_rows(model_summ_prev, model_summ_rec, model_summ_scale) %>%
  arrange(name, num, model)
