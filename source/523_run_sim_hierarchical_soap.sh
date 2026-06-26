#!/bin/bash
#SBATCH --job-name=sim_pos
#SBATCH --partition=wrobel
#SBATCH --output=sim_pos.out
#SBATCH --error=sim_pos.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript 522_run_sim_positive.R $JOBID


