#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=deg_analysis
#SBATCH --mem=150G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-06:00:00
#SBATCH --partition=gpu-miridan-pool
#SBATCH --qos=owner
#SBATCH --gres=gpu:0

set -euo pipefail

# Usage:
#   sbatch python_deg_analysis.sh                        # run v3 + v2 (default)
#   sbatch python_deg_analysis.sh <config.yaml> ...      # run specific config(s)
#
# Results saved in:
#   results/deg_analysis/{dataset_name}_{YYYYMMDD}/

NDD_ROOT="/miridan-data/annaludmir/ndd_gene_modules"

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd "$NDD_ROOT"

if [[ $# -gt 0 ]]; then
  CONFIGS=("$@")
else
  CONFIGS=(
    config_files/deg_combined_v3_config.yaml
    config_files/deg_combined_v2_config.yaml
  )
fi

for CONFIG in "${CONFIGS[@]}"; do
  echo ""
  echo "============================================================"
  echo "DEG analysis: $CONFIG"
  echo "============================================================"
  mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
    python -u modules/deg_analysis.py "$CONFIG"
done

rc=$?
echo ""
echo "Python exit code: $rc"
exit $rc
