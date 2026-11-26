#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=umap_no_cc
#SBATCH --account=miridan-users_v2
#SBATCH --partition=power-general-shared-pool
#SBATCH --qos=public
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --time=2-00:00:00
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err

# New UMAP of human_dev whose highly variable genes exclude the cell cycle
# genes (config_files/cell_cycle_genes_tirosh.json, all sets together).
# See modules/umap_without_cc_genes.py.
# Outputs: results/additional_analyses/umap_no_cc_genes/<label>/
#   umap_no_cc_genes.h5ad (obs + X_umap / X_pca; usable as --umap-from),
#   hvg_genes.csv, umap_<column>.png
#
# TEST=true: 100k random cells first, to check the run (minutes).

set -euo pipefail

module load mamba/mamba-1.5.8
cd /miridan-data/annaludmir/ndd_gene_modules

TEST=${TEST:-false}

ARGS=(
  --h5ad data/human_dev_without_week_5_cc_annotated_rescaled.h5ad
  --cc-genes config_files/cell_cycle_genes_tirosh.json
  --n-hvg 3000
  --hvg-cells 200000          # random cells used to pick the HVGs
  --batch-key Chemistry       # HVGs picked per chemistry
  --n-pcs 50
  --n-neighbors 15
  --min-dist 0.5
  --color CellCyclePhase CellClass Region RegionsGeneral Chemistry
  # --harmony                 # correct Chemistry in PCA space (needs harmonypy in the env)
)
if [[ "$TEST" == "true" ]]; then
  ARGS+=(--max-cells 100000 --label human_dev_test)
else
  ARGS+=(--label human_dev)
fi

export NUMBA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/umap_without_cc_genes.py "${ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
