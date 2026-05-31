#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=gene_clusters
#SBATCH --mem=200G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-06:00:00
# # public
# #SBATCH --partition=power-general-public-pool
# #SBATCH --qos=public
# # miris partition
#SBATCH --partition=gpu-miridan-pool
#SBATCH --qos=owner
#SBATCH --gres=gpu:0

set -euo pipefail

# Usage:
#   sbatch python_gene_clusters.sh                            # run combined v3 + v2 (default)
#   sbatch python_gene_clusters.sh <config.yaml> ...          # run specific config(s)
#   sbatch python_gene_clusters.sh --gene-list <genes.csv>    # highlight gene groups in all UMAPs
#   sbatch python_gene_clusters.sh config.yaml --gene-list genes.csv
#
# --gene-list  CSV with columns 'gene' and 'group'.  Each group gets one overlay
#              panel added to every UMAP figure: group genes in colour, all others gray.
#
# Default configs:
#   config_files/gene_clusters_combined_v3_config.yaml
#   config_files/gene_clusters_combined_v2_config.yaml
#
# Results are saved in:
#   results/gene_clusters/{dataset_name}_{YYYYMMDD}/

NDD_ROOT="/miridan-data/annaludmir/ndd_gene_modules"

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd "$NDD_ROOT"

# Separate --gene-list flag from config arguments
GENE_LIST_ARG=()
CONFIGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --gene-list)
      GENE_LIST_ARG=(--gene-list "$2")
      shift 2
      ;;
    *)
      CONFIGS+=("$1")
      shift
      ;;
  esac
done

if [[ ${#CONFIGS[@]} -eq 0 ]]; then
  CONFIGS=(
    config_files/gene_clusters_combined_v3_config.yaml
    config_files/gene_clusters_combined_v2_config.yaml
  )
fi

for CONFIG in "${CONFIGS[@]}"; do
  echo ""
  echo "============================================================"
  echo "Gene clustering: $CONFIG"
  echo "============================================================"
  mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
    python -u modules/gene_clustering.py "$CONFIG" "${GENE_LIST_ARG[@]}"
done

rc=$?
echo ""
echo "Python exit code: $rc"
exit $rc
