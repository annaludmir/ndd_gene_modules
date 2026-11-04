#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=cell_class_timeline
#SBATCH --mem=32G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-01:00:00
#SBATCH --partition=power-general-shared-pool
#SBATCH --qos=public

# % of all cells of a cell class per region and age (one line per region).
# Outputs: results/additional_analyses/cell_class_timeline/<label>/

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Options: one per line; comment a line out (or in) to switch it.
ARGS=(
  --h5ad data/human_dev.h5ad
  --label human_dev
  --cell-class "Neuronal IPC"     # several allowed: --cell-class "Neuronal IPC" "Radial glia" Neuroblast
  # --by-region                   # also per Region as annotated
  # --chemistry v3                # one chemistry only (default: all cells)
)

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/cell_class_timeline.py "${ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
