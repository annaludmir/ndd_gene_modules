#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=tf_validation_plot
#SBATCH --mem=8G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-00:30:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

set -euo pipefail

# Plots already-existing tf_validation CSVs. New tf_validation_report runs
# draw these figures themselves (see its --plot-style flag).

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Edit these before submitting.
# A single tf_validation subfolder, or results/tf_validation to plot them all.
VALIDATION_DIR="results/tf_validation"
PLOT_STYLE="connected"   # connected | hub_spoke
TOP_N_TFS=12
# Global per-pair scan of the matching ATAC run; stars unvalidated (grey)
# edges whose motif hit is in an accessible peak (all peaks). Leave empty to skip.
MOTIF_HITS_CSV="results/atac_analysis/atac_first_trimester_brain_20260816/3_motif_target_validation/motif_target_pair_scores.csv"
# Literature CSVs, pooled across gene lists: a pair is tagged direct / plausible
# in every figure it appears in, whatever the CSV is named.
LITERATURE_DIR="results/tf_validation/literature_evidence"

ARGS=(--validation-dir "$VALIDATION_DIR" --plot-style "$PLOT_STYLE" --top-n-tfs "$TOP_N_TFS")
if [[ -n "$MOTIF_HITS_CSV" ]]; then
  ARGS+=(--motif-hits-csv "$MOTIF_HITS_CSV")
fi
if [[ -n "$LITERATURE_DIR" ]]; then
  ARGS+=(--literature-dir "$LITERATURE_DIR")
fi

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/tf_validation_plot.py "${ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
