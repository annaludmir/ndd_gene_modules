#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=gene_boxplots_umaps
#SBATCH --mem=250G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-06:00:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public
# # miris partition
# #SBATCH --partition=gpu-miridan-pool 
# #SBATCH --qos=owner
# #SBATCH --gres=gpu:0

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Edit these before submitting.
GENE_LIST="data/genes_for_analyses/ANX_only_no_autism.txt"
H5AD_PATH="data/Cortex_EMX1_louvain3_passedQC_PostM_rev1.h5ad"
COMPARISON_COLUMNS=(CellCyclePhase)
SUBFOLDER_NAME="ANX TFs only (not autism)"
CHEMISTRY="v3"

# Optional toggles — set to 1 to enable.
INCLUDE_ALL_CELLS=0   # boxplots include expression == 0 cells (default: expressing-only)
SKIP_UMAPS=0          # only produce boxplots

EXTRA_ARGS=()
[[ "$INCLUDE_ALL_CELLS" -eq 1 ]] && EXTRA_ARGS+=(--include-all-cells)
[[ "$SKIP_UMAPS"        -eq 1 ]] && EXTRA_ARGS+=(--skip-umaps)

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python modules/gene_boxplots_umaps.py \
    --h5ad-path "$H5AD_PATH" \
    --gene-list "$GENE_LIST" \
    --comparison-columns "${COMPARISON_COLUMNS[@]}" \
    --chemistry "$CHEMISTRY" \
    --subfolder-name "$SUBFOLDER_NAME" \
    "${EXTRA_ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
