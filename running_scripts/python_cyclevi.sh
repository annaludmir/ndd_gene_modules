#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=cyclevi
#SBATCH --account=miridan-users_v2
#SBATCH --partition=gpu-miridan-pool
#SBATCH --qos=owner
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=1-00:00:00
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err

# Cell cycle angle with CycleVI on the human_dev progenitors (all regions
# together, chemistry as batch covariate). See modules/cyclevi_annotation.py.
# Outputs: results/additional_analyses/cyclevi/<label>/
#
# One-time setup of the environment (as a job):
#   sbatch running_scripts/setup_cyclevi_env.sh
#
# A quick test first:  TEST=true sbatch running_scripts/python_cyclevi.sh
#   (50,000 random cells, 30 epochs, label <label>_test)

set -euo pipefail

module load mamba/mamba-1.5.8
cd /miridan-data/annaludmir/ndd_gene_modules
ENV=/miridan-data/annaludmir/conda-envs/cyclevi
export PYTHONNOUSERSITE=1   # ignore packages in ~/.local (they shadowed env packages before)
nvidia-smi || echo "[warning] nvidia-smi failed: no GPU visible"

LABEL="human_dev_progenitors"

# Options: one per line; comment a line out (or in) to switch it.
ARGS=(
  --h5ad data/human_dev_without_week_5_cc_annotated_rescaled.h5ad   # raw counts + CellCyclePhase (reference)
  --cell-class "Radial glia" "Neuronal IPC" Glioblast Oligo         # progenitor classes (Oligo = OPCs)
  --batch-key Chemistry          # batch covariate (v2 / v3)
  --n-hvg 4000                   # + all cell cycle marker genes
  --batch-size 256
  # --epochs 60                  # default: ~100k gradient steps (20..400 epochs)
  # --max-cells 300000           # random subsample
  --umap-from data/full_adata_with_umap.h5ad   # UMAP coloured by CycleVI phase
)
if [[ "${TEST:-false}" == "true" ]]; then
  ARGS+=(--max-cells 50000 --epochs 30)
  LABEL="${LABEL}_test"
fi

mamba run -p "$ENV" python -u modules/cyclevi_annotation.py --label "$LABEL" "${ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
