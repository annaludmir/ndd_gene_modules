#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=hypergeometric
#SBATCH --mem=32G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-01:00:00
#SBATCH --partition=gpu-miridan-pool
#SBATCH --qos=owner
#SBATCH --gres=gpu:0

set -euo pipefail

# Usage:
#   sbatch python_hypergeometric.sh <gene_list.csv>
#
#   Runs both v3 and v2 configs with the provided gene list.
#
# Examples:
#   sbatch python_hypergeometric.sh results/gene_clusters/microcephaly_leading_genes.csv
#   sbatch python_hypergeometric.sh gene_lists/my_genes.csv config_files/hypergeometric_combined_v3_config.yaml
#
# If a second argument is given it is treated as a specific config to run instead of the defaults.

NDD_ROOT="/miridan-data/annaludmir/ndd_gene_modules"

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd "$NDD_ROOT"

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch python_hypergeometric.sh <gene_list.csv> [config.yaml ...]"
  exit 1
fi

GENE_LIST="$1"
shift

if [[ $# -gt 0 ]]; then
  CONFIGS=("$@")
else
  CONFIGS=(
    config_files/hypergeometric_combined_v3_config.yaml
    config_files/hypergeometric_combined_v2_config.yaml
  )
fi

for CONFIG in "${CONFIGS[@]}"; do
  echo ""
  echo "============================================================"
  echo "Hypergeometric enrichment: $CONFIG"
  echo "Gene list: $GENE_LIST"
  echo "============================================================"
  mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
    python -u modules/hypergeometric_enrichment.py "$CONFIG" "$GENE_LIST"
done

rc=$?
echo ""
echo "Python exit code: $rc"
exit $rc
