#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=cc_umap
#SBATCH --mem=64G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-02:00:00
#SBATCH --partition=power-general-shared-pool
#SBATCH --qos=public

# UMAP coloured by CellCyclePhase. See modules/cell_cycle_umap.py.
# Outputs: results/additional_analyses/cell_cycle_umap/<label>/

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Options: one per line; comment a line out (or in) to switch it.
ARGS=(
  --h5ad data/human_dev_without_week_5_cc_annotated_rescaled.h5ad   # rescaled cutoffs
  --label human_dev_rescaled
  # --h5ad data/human_dev_without_week_5_cc_annotated_fraction.h5ad # spec cutoffs
  # --label human_dev_fraction
  --umap-from data/full_adata_with_umap.h5ad   # X_umap, matched by cell ID
  --per-phase                    # also one UMAP per phase
  # --point-size 0.3
  --subsample 0.2               # random 20% of cells (or a number, e.g. 300000); remove for all
)

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/cell_cycle_umap.py "${ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
