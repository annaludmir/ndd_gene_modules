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

# Usage: sbatch python_gsea_tau_comparison.sh <gene_lists_folder>
#
# Runs GSEA enrichment in 4 variants (v2/v3 × with/without tau filtering)
# on all three scopes (cortex, cell_phase, all_layers) for every gene list
# CSV in <gene_lists_folder>.
#
# Produces one batch summary CSV per gene list:
#   results/enrichment_results/{gene_list}_batch_summary_tau_vs_v2_v3_{date}.csv

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <gene_lists_folder>" >&2
  exit 1
fi

GENE_LISTS_FOLDER="$1"
NDD_ROOT="/miridan-data/annaludmir/ndd_gene_modules"

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd "$NDD_ROOT"

echo "Gene lists folder: $GENE_LISTS_FOLDER"
echo "NDD root:          $NDD_ROOT"

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/enrichment_cal_lists_loop_tau_comparison.py \
  "$GENE_LISTS_FOLDER" \
  "$NDD_ROOT"

rc=$?
echo "Python exit code: $rc"
exit $rc
