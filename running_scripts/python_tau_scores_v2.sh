#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=python
#SBATCH --mem=500G
# #SBATCH --mem=250G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=1-00:00:00
# # public
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public
# # miris partition
# #SBATCH --partition=gpu-miridan-pool
# #SBATCH --qos=owner
# #SBATCH --gres=gpu:0

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

CONFIGS=(
  config_files/ges_score_cortex_config_v2.yaml
  config_files/ges_score_cortex_cell_phase_config_v2.yaml
  config_files/ges_score_all_layers_config_copy_v2.yaml
)

for CONFIG in "${CONFIGS[@]}"; do
  echo ""
  echo "=============================="
  echo "Running tau pipeline: $CONFIG"
  echo "=============================="
  mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
    python -u modules/tau_pipeline.py "$CONFIG"
done

rc=$?
echo "Python exit code: $rc"
exit $rc
