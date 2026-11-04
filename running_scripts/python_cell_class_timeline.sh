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

# MODE=class: % of all cells of a class per region and age (default).
# MODE=phase: % of each class's cycling cells in S / G2M per age, a line per class.
#   MODE=phase sbatch running_scripts/python_cell_class_timeline.sh
# MODE=phase_lines: one figure per class with a line per phase (G1, S, G2M,
#   Post-M), % of that class's cycling cells per age.
MODE="${MODE:-class}"

# Options: one per line; comment a line out (or in) to switch it.
if [[ "$MODE" == "phase_lines" ]]; then
  ARGS=(
    --mode phase --lines phase
    --h5ad data/human_dev_without_week_5_cc_annotated_fraction.h5ad   # needs CellCyclePhase
    --label human_dev_without_week_5
    --cell-class Glioblast        # one figure per class
    --phases G1 S G2M Post-M
    --min-cycling-cells 20        # leave out ages with fewer cycling cells
    --by-region-general           # also one figure per Forebrain / Midbrain / Hindbrain
    # --chemistry v3              # one chemistry only (default: all cells)
  )
elif [[ "$MODE" == "phase" ]]; then
  ARGS=(
    --mode phase
    --h5ad data/human_dev_without_week_5_cc_annotated_fraction.h5ad   # needs CellCyclePhase
    --label human_dev_without_week_5
    --cell-class "Neuronal IPC" "Radial glia" Glioblast
    --phases S G2M                # % of cycling cells (G1, S, G2M, Post-M) in each phase
    --min-cycling-cells 20        # leave out class x age points with fewer cycling cells
    --by-region-general           # also one set of figures per Forebrain / Midbrain / Hindbrain
    # --chemistry v3              # one chemistry only (default: all cells)
  )
else
  ARGS=(
    --h5ad data/human_dev.h5ad
    --label human_dev
    --cell-class "Neuronal IPC"   # several allowed: --cell-class "Neuronal IPC" "Radial glia" Neuroblast
    # --by-region                 # also per Region as annotated
    # --chemistry v3              # one chemistry only (default: all cells)
  )
fi

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/cell_class_timeline.py "${ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
