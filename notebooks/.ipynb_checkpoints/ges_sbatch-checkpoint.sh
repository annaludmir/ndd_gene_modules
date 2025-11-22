#!/bin/bash
#SBATCH --job-name=search_enrichment
#SBATCH --account=power-miridan-users
#SBATCH --partition=power-general
#SBATCH --output=my_job_%j.out
#SBATCH --error=my_job_%j.err
#SBATCH --mem=300G
#SBATCH --cpus-per-task=10

module load mamba/mamba1.4.2-environmentally
mamba activate /scratch200/reutj/conda-envs/jupyter-scanpy

cd /scratch200/reutj/notebooks/

python -u ges_scores_forreal.py ['CellCycleStatus'] ['proliferating_non_cycling','differentiating_non_cycling'] 

wait

