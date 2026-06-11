#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=gene_axes_scatter
#SBATCH --mem=32G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-02:00:00
#SBATCH --partition=gpu-miridan-pool
#SBATCH --qos=owner
#SBATCH --gres=gpu:0

set -euo pipefail

# Usage:
#   sbatch python_gene_axes_scatter.sh [config.yaml] [--gene-list genes.csv] [--gene-clusters gene_clusters.csv]
#
# Default config: gene_axes_scatter_v3_config.yaml
#
# Results saved in:
#   results/gene_axes_scatter/{dataset_name}_{YYYYMMDD}/

NDD_ROOT="/miridan-data/annaludmir/ndd_gene_modules"

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd "$NDD_ROOT"

CONFIG="${1:-config_files/gene_axes_scatter_v3_config.yaml}"
shift || true   # remaining args passed through to the python script

echo ""
echo "============================================================"
echo "Gene Axes Scatter: $CONFIG"
echo "============================================================"

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/gene_axes_scatter.py "$CONFIG" "$@"

rc=$?
echo ""
echo "Python exit code: $rc"
exit $rc
