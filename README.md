---
editor_options: 
  markdown: 
    wrap: sentence
---

# 202403_porras

## Craniosynostosis Image Model Project Directory

This project contains analysis and simulations for modeling craniosynostosis imaging data.

The first iteration of the project analyzed volumes across individuals at the level of a region of the skull.

The second iteration analyzes individuals' growth at each of 52,800 registered points along the surface of the skull.

## Structure

-   `analysis/` contains Quarto markdown files that compile analyses into reports.
-   `data/` contains cleaned data, intermediate objects for the analysis, and a shortcut to the raw data in Josh's OneDrive.
-   `drafts/` contains presentations from meetings and collaborators as well as abstracts for conferences and grant proposals.
-   `results/` contains figures exported by the analysis files and model fit objects from the source folder.
-   `simulations/` will contain bare scripts and functions to run simulations for the project. This folder currently only contains the QMD file that we have used to run a simulation iteration and an unfinished R script written before that QMD file.
-   `source/` contains bare scripts and functions that perform the analysis for this project.

## Replicating this Analysis

To run the hierarchical soap film model in this analysis, use scripts in the `source` folder.

First, edit `source/000_definitions.R` so that variables ending in "\_FOLDER" point to the "unsigned distances" data and single covariates file created by Josh's preliminary data cleaning.

Second, run scripts prefixed with `00` in numerical order: 

- `001` cleans the dataset of covariates about each image.
- `002` loads and cleans the point-level data of growth distances that make up each image, each stored in a matrix of growth values by location in the spherical projection of the image.
This script also identifies images whose growth values are duplicated even though their labels are not identical.
These images are removed from the analysis, so the final subset of image-level covariates is created in this script as well.
- `003` uses coordinates from the cleaned pointwise data to construct a soap film smoother over this space.
- `004` uses the cleaned image-level covariate data to construct a smooth over age, separated by fusion type and "reference-coded," so that we can capture the growth patterns of a normative patient and also capture how the presence of each suture fusion is associated with deviations from this growth pattern varying over age.
- `005` uses the cleaned pointwise data to construct a large matrix of growth values at each point in each image.
- `006` collects cleaned data from previous scripts into one RDA file, so that this single file can be uploaded to the cluster and loaded in later scripts to run the main model for the analysis.

Third, run the script prefixed `201` to run the Gibbs sampler and save a list object of parameter estimates at each iteration.

Fourth, run the script prefixed `301` to generate trace plots and predictions from the estimated model.

Other scripts in the source folder are for additional analyses, but are not required to run the main model.

## Source Folder

The `source` folder is the most complicated of these folders, so we describe it in more detail below.

We prefix scripts in this folder to organize them: \* `000_definitions.R`: loads packages, functions, and constants we will use in the analysis.
\* `0xx`: clean raw data to prepare for the model.
\* `1xx`: summarize the cleaned data from `0xx`.
\* `2xx`: run the main hierarchical soap film model in our analysis.
\* `3xx`: summarize and visualize results from the main model in `2xx`.
\* `5xx`: run preliminary models (region-level regressions, pixel-level regressions, tensor product smooths) \* `6xx`: summarize and visualize results from preliminary models in `5xx`.
\* `7xx`: run simpler models we used to explore how the soap film model worked (models excluding boundary basis functions, frequentist models with GCV grid search, Bayesian soap film on one individual, etc.) \* `8xx`: summarize results from models in `8xx`.
Note: 801 summarizes results from 701, 802 summarizes results from 702, etc. \* `9xx`: "sandbox" scripts.
For example, 901 compares results from 2 iterations of results from 201 with different penalty matrix scales.
\* `fn_xxx`/`fns_xxx`: functions used in the analyses above.
All of these are sourced in `000`.

project structured using the [projectr](https://github.com/jeff-goldsmith/projectr) package.
