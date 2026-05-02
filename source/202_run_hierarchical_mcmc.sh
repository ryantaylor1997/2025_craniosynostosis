#!/bin/bash
#SBATCH --job-name=run_cranio_mcmc
#SBATCH --partition=wrobel
#SBATCH --output=run_cranio_mcmc.out
#SBATCH --error=run_cranio_mcmc.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript 201_run_hierarchical_mcmc.R $JOBID


