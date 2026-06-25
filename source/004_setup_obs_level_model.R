################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Create subject-level characteristics modeling objects
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load cleaned / filtered subject-level data ("cranio_clean")
load(file = here("data", "cleaned", "obs_data_clean.rda"))

# Construct smooth --------------------------------------------------------

# Construct smooth without constraint or factor levels
age_smooth <- s(age, bs = age_bs, k = cranio_knots_age)

## Create design and penalty matrices using functions we wrote
obs_smooth_by_fusion <- construct_reference_smooth(
  sm = age_smooth, dat = cranio_clean,
  by_var = "fusion_type", param_formula = ~sex + fusion_type
)

# Combine suture fusion penalties so we estimate the same lambda
obs_S_combo <- list(obs_smooth_by_fusion$S[[1]],
                    Reduce("+", obs_smooth_by_fusion$S[-1]))

# Adjust penalty matrices so trace of inverse is similar to identity
obs_pen_eigen <- eigen(obs_S_combo[[2]],
                       symmetric = T, only.values = T)$values

obs_pen_scale <- sum(
  1 / obs_pen_eigen[which(abs(obs_pen_eigen) > 1e-10)]
) / sum(obs_pen_eigen != 0)

obs_S_combo[-1] <- map(obs_S_combo[-1], ~ .x * obs_pen_scale)

# Collect design and penalty back into one list
obs_smooth_list <- obs_smooth_by_fusion
obs_smooth_list[["S"]] <- obs_S_combo

# Save output -------------------------------------------------------------

# Save intermediate steps to share
save(obs_smooth_list,
     file = here::here("data", "cleaned", "obs_level_smooth.rda"))

# Save penalty scale terms
save(obs_pen_scale,
     file = here::here("data", "intermediate", "obs_level_penalty_scales.rda"))

# Save original reference smooth
save(obs_smooth_by_fusion,
     file = here("data", "intermediate", "obs_level_reference_smooth.rda"))

