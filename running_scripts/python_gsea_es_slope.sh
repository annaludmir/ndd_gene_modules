#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=gsea_es_slope
#SBATCH --mem=200G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-12:00:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

set -euo pipefail

# Usage:
#   sbatch python_gsea_es_slope.sh <gene_list.csv>
#   sbatch python_gsea_es_slope.sh <gene_list.csv> <config.yaml> ...
#
# First argument is always the gene list.
# Remaining arguments are optional config files (default: combined v3).
#
# Results saved in:
#   results/gsea_es_slope/{dataset_name}_{gene_list_name}_{YYYYMMDD}/

NDD_ROOT="/miridan-data/annaludmir/ndd_gene_modules"

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd "$NDD_ROOT"

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch python_gsea_es_slope.sh <gene_list.csv> [config.yaml ...]"
  exit 1
fi

GENE_LIST="$1"
shift

if [[ $# -gt 0 ]]; then
  CONFIGS=("$@")
else
  CONFIGS=(
    config_files/gsea_es_slope_combined_v3_config.yaml
  )
fi

for CONFIG in "${CONFIGS[@]}"; do
  echo ""
  echo "============================================================"
  echo "GSEA ES Slope: $CONFIG"
  echo "Gene list: $GENE_LIST"
  echo "============================================================"
  mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
    python -u modules/gsea_es_slope.py "$CONFIG" "$GENE_LIST"
done

rc=$?
echo ""
echo "Python exit code: $rc"
exit $rc
