#!/bin/bash
#SBATCH --job-name=sim_hier
#SBATCH --partition=wrobel
#SBATCH --output=sim_hier.out
#SBATCH --error=sim_hier.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript 512_run_sim_hierarchical_soap.R $JOBID


