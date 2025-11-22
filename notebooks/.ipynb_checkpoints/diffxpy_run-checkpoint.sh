#!/bin/bash
#SBATCH --job-name=adj_spec_rg_try
#SBATCH --account=power-miridan-users
#SBATCH --partition=power-general
#SBATCH --output=my_job_%j.out
#SBATCH --error=my_job_%j.err
#SBATCH --mem=500G
#SBATCH --cpus-per-task=20

module load mamba/mamba1.4.2-environmentally
mamba activate /scratch200/reutj/conda-envs/diffxpy_env

python -u diffxpy_test.py 
