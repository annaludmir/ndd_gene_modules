#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=tf_validation
#SBATCH --mem=8G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-00:30:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Edit these before submitting.
TF_NETWORK_DIR="results/tf_network/tf_network_cortex_v3_20260802"
ATAC_PER_CT_DIR="results/atac_analysis/atac_first_trimester_brain_20260816/3_motif_target_validation/per_cell_type"
SUBFOLDER_NAME="cortex_v3_20260802__vs__atac_20260816"

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/tf_validation_report.py \
    --tf-network-dir       "$TF_NETWORK_DIR" \
    --atac-per-cell-type-dir "$ATAC_PER_CT_DIR" \
    --subfolder-name       "$SUBFOLDER_NAME"

rc=$?
echo "Python exit code: $rc"
exit $rc
