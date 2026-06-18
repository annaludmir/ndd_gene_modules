#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=tf_network
#SBATCH --mem=500G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=1-00:00:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public
#SBATCH --cpus-per-task=8
# # miris partition
##SBATCH --partition=gpu-miridan-pool
##SBATCH --qos=owner
##SBATCH --gres=gpu:0

set -euo pipefail

# Usage:
#   # Full pipeline (aggregation + GRNBoost2 + optional cisTarget + AUCell):
#   sbatch python_tf_network.sh
#
#   # Full pipeline + query a gene set:
#   sbatch python_tf_network.sh --gene-list gene_lists/my_genes.csv
#
#   # Query only (pipeline already ran, just filter + plot):
#   sbatch python_tf_network.sh --gene-list gene_lists/my_genes.csv --query-only
#
#   # Custom config:
#   sbatch python_tf_network.sh --config config_files/tf_network_cortex_v3_config.yaml \
#                                --gene-list gene_lists/my_genes.csv
#
# Results saved in:
#   results/tf_network/{dataset_name}_{YYYYMMDD}/

NDD_ROOT="/miridan-data/annaludmir/ndd_gene_modules"

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd "$NDD_ROOT"

CONFIG="config_files/tf_network_cortex_v3_config.yaml"
GENE_LIST=""
QUERY_ONLY=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)     CONFIG="$2";     shift 2 ;;
    --gene-list)  GENE_LIST="$2";  shift 2 ;;
    --query-only) QUERY_ONLY="--query-only"; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

GENE_LIST_ARG=""
if [[ -n "$GENE_LIST" ]]; then
  GENE_LIST_ARG="--gene-list $GENE_LIST"
fi

echo ""
echo "============================================================"
echo "TF Network: $CONFIG"
[[ -n "$GENE_LIST" ]] && echo "Gene list:  $GENE_LIST"
[[ -n "$QUERY_ONLY" ]] && echo "Mode:       query-only"
echo "============================================================"

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/tf_network.py "$CONFIG" $GENE_LIST_ARG $QUERY_ONLY

rc=$?
echo ""
echo "Python exit code: $rc"
exit $rc
