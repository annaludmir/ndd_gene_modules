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
ATAC_RUN="results/atac_analysis/atac_first_trimester_brain_20260816/3_motif_target_validation"
SUBFOLDER_NAME="cortex_v3_20260802_vs_atac_20260816"

# Either axis may be left empty. Cell phase needs per_cell_cycle_scoring to
# have been enabled in the ATAC config, which produces per_cell_cycle/.
ATAC_PER_CT_DIR="${ATAC_RUN}/per_cell_type"
ATAC_PER_PHASE_DIR="${ATAC_RUN}/per_cell_cycle"
# Global per-pair scan: stars unvalidated edges whose motif hit is in an
# accessible peak when all peaks are considered.
ATAC_PAIR_SCORES="${ATAC_RUN}/motif_target_pair_scores.csv"
# Literature CSVs, pooled across gene lists; tags direct / plausible evidence.
LITERATURE_DIR="results/tf_validation/literature_evidence"

# Validation network figures (figures/ in the output dir).
PLOT_STYLE="connected"   # connected | hub_spoke | none
TOP_N_TFS=12

ARGS=(--tf-network-dir "$TF_NETWORK_DIR" --subfolder-name "$SUBFOLDER_NAME"
      --plot-style "$PLOT_STYLE" --top-n-tfs "$TOP_N_TFS")
if [[ -n "$ATAC_PER_CT_DIR" && -d "$ATAC_PER_CT_DIR" ]]; then
  ARGS+=(--atac-per-cell-type-dir "$ATAC_PER_CT_DIR")
else
  echo "[note] no per_cell_type dir at ${ATAC_PER_CT_DIR} — skipping that axis."
fi
if [[ -n "$ATAC_PER_PHASE_DIR" && -d "$ATAC_PER_PHASE_DIR" ]]; then
  ARGS+=(--atac-per-cell-cycle-dir "$ATAC_PER_PHASE_DIR")
else
  echo "[note] no per_cell_cycle dir at ${ATAC_PER_PHASE_DIR} — skipping that axis."
  echo "       Run atac_seq_analysis with per_cell_cycle_scoring: true first."
fi
if [[ -n "$ATAC_PAIR_SCORES" && -f "$ATAC_PAIR_SCORES" ]]; then
  ARGS+=(--atac-pair-scores-csv "$ATAC_PAIR_SCORES")
else
  echo "[note] no pair scores at ${ATAC_PAIR_SCORES} — figures won't mark accessible motif hits."
fi
if [[ -n "$LITERATURE_DIR" && -d "$LITERATURE_DIR" ]]; then
  ARGS+=(--literature-dir "$LITERATURE_DIR")
else
  echo "[note] no literature dir at ${LITERATURE_DIR} — figures won't tag literature evidence."
fi

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/tf_validation_report.py "${ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
