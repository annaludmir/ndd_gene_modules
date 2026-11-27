#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=maturity_umaps
#SBATCH --mem=96G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-08:00:00
#SBATCH --partition=power-general-shared-pool
#SBATCH --qos=public

# Where are the neuron maturity markers expressed (cortex only or not)? Per
# marker (default and braun2023 sets): UMAP of regions | UMAP of expression |
# % detected and mean per region, on human_dev.
# Two runs: all cells, and neurons only (the cells the annotation was made on).
# See modules/maturity_marker_umaps.py.
# Outputs: results/additional_analyses/maturity_marker_umaps/<label>/

set -euo pipefail

module load mamba/mamba-1.5.8
cd /miridan-data/annaludmir/ndd_gene_modules

H5AD=data/human_dev_without_week_5_cc_annotated_neuron_type_braun2023.h5ad
COMMON=(
  --h5ad "$H5AD"
  --umap-from data/full_adata_with_umap.h5ad   # X_umap, matched by cell ID
  --marker-sets default braun2023
  --region-col Region
  # --chemistry v3
  # --subsample 300000                         # plot fewer points (numbers use all cells)
)

run() {
  mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
    python -u modules/maturity_marker_umaps.py "${COMMON[@]}" "$@"
}

rc=0
run --label human_dev_all_cells || rc=1
run --cell-class Neuron --label human_dev_neurons || rc=1

echo "Exit code: $rc"
exit $rc
