#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=cc_conflicts
#SBATCH --mem=32G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-00:30:00
#SBATCH --partition=power-general-shared-pool
#SBATCH --qos=public

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Cells whose cycling score is <= 0.004 (Non-cycling by the rule) but whose
# CellCyclePhase is not Non-cycling, from each dataset, into one CSV.
mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/phase_score_conflicts.py \
    --dataset cortex=data/Cortex_EMX1_louvain3_passedQC_PostM_rev1.h5ad \
    --dataset human_dev=data/human_dev_without_week_5_cc_annotated_fraction.h5ad \
    --out results/additional_analyses/cell_cycle_comparison/phase_score_conflicts.csv

rc=$?
echo "Python exit code: $rc"
exit $rc
