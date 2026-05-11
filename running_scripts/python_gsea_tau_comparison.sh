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

# Usage: sbatch python_gsea_tau_comparison.sh <gene_list.csv|folder> [tau_percentile]
#
# Runs GSEA enrichment in 4 variants (v2/v3 × with/without tau filtering)
# on all three scopes (cortex, cell_phase, all_layers).
# Accepts either a single gene-list CSV or a folder of CSVs.
#
# tau_percentile (optional): overrides the tau_percentile in all tau-filtered configs.
#   Default: read from each config file (typically 90).
#
# Produces one batch summary CSV per gene list:
#   results/enrichment_results/{gene_list_stem}_batch_summary_tau_vs_v2_v3_{date}.csv

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <gene_list.csv|folder> [tau_percentile]" >&2
  exit 1
fi

INPUT="$1"
TAU_PCT="${2:-}"
NDD_ROOT="/miridan-data/annaludmir/ndd_gene_modules"

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd "$NDD_ROOT"

echo "Input:          $INPUT"
echo "Tau percentile: ${TAU_PCT:-from config}"
echo "NDD root:       $NDD_ROOT"

TAU_ARG=()
if [[ -n "$TAU_PCT" ]]; then
  TAU_ARG=(--tau-percentile "$TAU_PCT")
fi

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/enrichment_cal_lists_loop_tau_comparison.py \
  "$INPUT" \
  --ndd-root "$NDD_ROOT" \
  "${TAU_ARG[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
