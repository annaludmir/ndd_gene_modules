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
#   sbatch python_gene_clusters.sh                   # run combined v3 + v2 (default)
#   sbatch python_gene_clusters.sh <config.yaml> ... # run specific config(s)
#
# Default configs:
#   config_files/gene_clusters_combined_v3_config.yaml
#   config_files/gene_clusters_combined_v2_config.yaml
#
# Single-dataset configs are also available:
#   config_files/gene_clusters_cortex_config.yaml
#   config_files/gene_clusters_all_layers_config.yaml
#
# Results are saved in:
#   results/gene_clusters/{dataset_name}_{YYYYMMDD}/

NDD_ROOT="/miridan-data/annaludmir/ndd_gene_modules"

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd "$NDD_ROOT"

if [[ $# -ge 1 ]]; then
  CONFIGS=("$@")
else
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
    python -u modules/gene_clustering.py "$CONFIG"
done

rc=$?
echo ""
echo "Python exit code: $rc"
exit $rc
