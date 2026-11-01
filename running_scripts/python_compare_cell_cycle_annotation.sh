#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=cc_compare
#SBATCH --mem=64G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-01:00:00
#SBATCH --partition=power-general-shared-pool
#SBATCH --qos=public

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Edit these before submitting.
FILE_A="data/Cortex_EMX1_louvain3_passedQC_PostM_rev1.h5ad"
FILE_B="data/human_dev_without_week_5_cc_annotated_fraction.h5ad"
LABEL_A="cortex"
LABEL_B="human_dev"
OUT_DIR="results/additional_analyses/cell_cycle_comparison"

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/compare_cell_cycle_annotation.py \
    --a "$FILE_A" --b "$FILE_B" --label-a "$LABEL_A" --label-b "$LABEL_B" \
    --out "$OUT_DIR"

rc=$?
echo "Python exit code: $rc"
exit $rc
