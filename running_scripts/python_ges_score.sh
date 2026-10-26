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
#SBATCH --partition=power-general-shared-pool
#SBATCH --qos=public
# # miris partition
# #SBATCH --partition=gpu-miridan-pool
# #SBATCH --qos=owner
# #SBATCH --gres=gpu:0
# # deprecated
# # SBATCH --account=public-users_v2

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Usage:
#   sbatch running_scripts/python_ges_score.sh [config_yaml]
# e.g.
#   sbatch running_scripts/python_ges_score.sh config_files/ges_score_all_layers_neuron_type_config.yaml
CONFIG="${1:-config_files/ges_score_all_layers_cell_phase_config.yaml}"

if [[ ! -f "$CONFIG" ]]; then
  echo "Error: config file not found: $CONFIG"
  exit 1
fi
echo "GES config: $CONFIG"

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/specificity_score_calculations.py "$CONFIG"

rc=$?
echo "Python exit code: $rc"
exit $rc
