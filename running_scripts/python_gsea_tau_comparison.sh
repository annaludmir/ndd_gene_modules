#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=python
#SBATCH --mem=500G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=2-00:00:00
# public
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

set -euo pipefail

# Usage: sbatch python_gsea_tau_comparison.sh <gene_list.csv|folder> [tau_percentile] [tau_score_cutoff]
#
# Runs GSEA enrichment in 4 variants (v2/v3 × with/without tau filtering)
# on all three scopes (cortex, cell_phase, all_layers).
# Accepts either a single gene-list CSV or a folder of CSVs.
#
# tau_percentile  (optional, arg 2): keep genes at or above the Nth percentile.
#   Default: read from each config file (typically 90).
#
# tau_score_cutoff (optional, arg 3): keep genes with tau >= this absolute value
#   (e.g. 0.5). Mutually exclusive with tau_percentile — pass "" as arg 2 to skip it.
#   Results folders will use 'tauscore{value}' instead of 'tau{pct}'.
#
# Produces one batch summary CSV per gene list:
#   results/enrichment_results/{gene_list_stem}_batch_summary_{tau_label}_vs_v2_v3_{date}.csv

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <gene_list.csv|folder> [tau_percentile] [tau_score_cutoff]" >&2
  exit 1
fi

INPUT="$1"
TAU_PCT="${2:-}"
TAU_SCORE="${3:-}"
NDD_ROOT="/miridan-data/annaludmir/ndd_gene_modules"

if [[ -n "$TAU_PCT" && -n "$TAU_SCORE" ]]; then
  echo "Error: tau_percentile and tau_score_cutoff are mutually exclusive." >&2
  exit 1
fi

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd "$NDD_ROOT"

echo "Input:            $INPUT"
echo "Tau percentile:   ${TAU_PCT:-from config}"
echo "Tau score cutoff: ${TAU_SCORE:-(none)}"
echo "NDD root:         $NDD_ROOT"

TAU_ARG=()
if [[ -n "$TAU_PCT" ]]; then
  TAU_ARG=(--tau-percentile "$TAU_PCT")
fi

TAU_SCORE_ARG=()
if [[ -n "$TAU_SCORE" ]]; then
  TAU_SCORE_ARG=(--tau-score-cutoff "$TAU_SCORE")
fi

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/enrichment_cal_lists_loop_tau_comparison.py \
  "$INPUT" \
  --ndd-root "$NDD_ROOT" \
  "${TAU_ARG[@]}" \
  "${TAU_SCORE_ARG[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
