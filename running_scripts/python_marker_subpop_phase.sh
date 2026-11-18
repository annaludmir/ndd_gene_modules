#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=marker_subpop_phase
#SBATCH --mem=64G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-02:00:00
#SBATCH --partition=power-general-shared-pool
#SBATCH --qos=public

# Neuronal IPC subpopulations by marker (SOX2+ proliferating, NEUROD6+
# differentiating) and the % of each in every cell cycle phase. See
# modules/marker_subpop_phase.py.
# Outputs: results/additional_analyses/marker_subpop_phase/<label>/

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Options: one per line; comment a line out (or in) to switch it.
ARGS=(
  --h5ad data/Cortex_EMX1_louvain3_passedQC_PostM_rev1.h5ad   # has its own CellCyclePhase
  --label cortex
  # --h5ad data/human_dev_without_week_5_cc_annotated_rescaled.h5ad   # human_dev instead
  # --label human_dev_without_week_5
  --cell-class "Neuronal IPC"
  --marker proliferating=SOX2
  --marker differentiating=NEUROD6
  --min-counts 1                # raw counts of the marker for a cell to be '+'
  # --by-region-general         # also per Forebrain / Midbrain / Hindbrain (human_dev)
  # --chemistry v3              # one chemistry only (default: all cells)
)

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/marker_subpop_phase.py "${ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
