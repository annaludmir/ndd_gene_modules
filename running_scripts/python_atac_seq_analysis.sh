#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=atac_seq
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

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

CONFIG="${1:-config_files/atac_seq_config.yaml}"
# Any additional args (--tasks, --run-dir, --force) are forwarded to Python.
# Usage examples:
#   sbatch python_atac_seq_analysis.sh                                      # run all
#   sbatch python_atac_seq_analysis.sh config.yaml --tasks motif_target_validation
#   sbatch python_atac_seq_analysis.sh config.yaml --run-dir results/atac_analysis/atac_first_trimester_brain_20260814 --tasks promoter
#   sbatch python_atac_seq_analysis.sh config.yaml --force
shift || true
EXTRA_ARGS=("$@")

echo "============================================================"
echo "ATAC-seq analysis: $CONFIG"
[[ ${#EXTRA_ARGS[@]} -gt 0 ]] && echo "Extra args:         ${EXTRA_ARGS[*]}"
echo "============================================================"

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/atac_seq_analysis.py "$CONFIG" "${EXTRA_ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
