################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Load packages, define constants, and source functions for work
################################################################################

rm(list = ls()); gc()

# Load packages -----------------------------------------------------------
pacman::p_load(
  concaveman, # Special function for concave hull
  mvtnorm, # Good multivariate normal
  gtsummary, # Summary tables
  broom, # Model tidying
  mgcv, # GAM Modeling
  Matrix, MASS, reshape2, # Useful matrix operations / generalized inverse
  ggridges, ggh4x, ggpubr, ggpmisc, # Useful extensions to ggplot
  gtools, readxl, # Extra programming tools
  tidyverse, magrittr, janitor, knitr, kableExtra, here # Basics
)

# Load functions ----------------------------------------------------------

# Load my Rcpp functions
library(DevCranio)

# Load my R functions
source(here("source", "load_matrix_data_fn.R"))
source(here("source", "obtain_predictions_fn.R"))
source(here("source", "lm_penalized_gibbs_fn.R"))
source(here("source", "construct_smooth_ref_code_fns.R"))
source(here("source", "hierarchical_penalized_gibbs_fn.R"))
source(here("source", "simulate_data_fn.R"))
source(here("source", "mv_normal_matrix_fn.R"))


# Set file paths ----------------------------------------------------------

# Set data folder
DATA_FOLDER <- here("../data/Lukemire, Joshua's files - data/")
COVARS_FOLDER <- here("../data/SphericalMaps_database/SphericalMaps_database/")

# Set subfolder with image data
MATRIX_FOLDER <- file.path(DATA_FOLDER, "unsigned_distances")

# Set constants -----------------------------------------------------------

# Define maximum age
cranio_max_age <- 365

# Define number of knots along each dimension of soap film
so_GRID_N <- 10

# Define number of knots along age effect
cranio_knots_age <- 10

# Set plotting constants --------------------------------------------------

# Define age category cutoffs
cranio_age_breaks <- c(seq(0, 300, 60), cranio_max_age)

# Label these with spaces so they can wrap in plots
cranio_age_break_labels <- map_chr(2:length(cranio_age_breaks),
                                   ~paste0(cranio_age_breaks[.-1],
                                           " - ",
                                           if_else(. == length(cranio_age_breaks),
                                                   "", "<"),
                                           cranio_age_breaks[.]))



# Set dictionaries for consistent merging ---------------------------------

# Dictionary of fusion indicator variables so we can keep 1 category column
fusion_dict <- tibble(
  fusion_type = c("Normative",
                  "Sagittal", "Metopic",
                  "Right Coronal", "Bicoronal", "Left Coronal"),
  fused_Sagittal = c(0, 1, 0, 0, 0, 0),
  fused_Metopic = c(0, 0, 1, 0, 0, 0),
  fused_RCoronal = c(0, 0, 0, 1, 1, 0),
  fused_LCoronal = c(0, 0, 0, 0, 1, 1)
)
