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
# #SBATCH --partition=power-general-public-pool
# #SBATCH --qos=public
# # miris partition
#SBATCH --partition=gpu-miridan-pool
#SBATCH --qos=owner
#SBATCH --gres=gpu:0
set -euo pipefail

# Usage: sbatch python_gsea_tau_comparison.sh <gene_list.csv|folder> [tau_score_cutoff ...]
#
# Runs GSEA enrichment in variants (v2/v3 × with/without tau filtering)
# on all three scopes (cortex, cell_phase, all_layers).
# Accepts either a single gene-list CSV or a folder of CSVs.
#
# tau_score_cutoff (optional, args 2+): one or more absolute tau score thresholds.
#   Genes with tau >= cutoff are kept. Tau-filtered variants run once per cutoff.
#   Results folders will use 'tauscore{value}' in their name.
#   If omitted, only non-tau variants run.
#
# Examples:
#   sbatch python_gsea_tau_comparison.sh my_genes.csv
#   sbatch python_gsea_tau_comparison.sh my_genes.csv 0.5
#   sbatch python_gsea_tau_comparison.sh my_genes.csv 0.3 0.5 0.7
#
# Produces one batch summary CSV per gene list:
#   results/enrichment_results/{gene_list_stem}_batch_summary_{tau_label}_vs_v2_v3_{date}.csv

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <gene_list.csv|folder> [tau_score_cutoff ...]" >&2
  exit 1
fi

INPUT="$1"
shift
TAU_SCORES=("$@")  # all remaining args are tau score cutoffs (may be empty)
NDD_ROOT="/miridan-data/annaludmir/ndd_gene_modules"

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd "$NDD_ROOT"

echo "Input:            $INPUT"
echo "Tau score cutoffs: ${TAU_SCORES[*]:-(none)}"
echo "NDD root:         $NDD_ROOT"

TAU_SCORE_ARG=()
if [[ ${#TAU_SCORES[@]} -gt 0 ]]; then
  TAU_SCORE_ARG=(--tau-score-cutoff "${TAU_SCORES[@]}")
fi

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/enrichment_cal_lists_loop_tau_comparison.py \
  "$INPUT" \
  --ndd-root "$NDD_ROOT" \
  "${TAU_SCORE_ARG[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
