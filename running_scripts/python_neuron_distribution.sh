#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=neuron_distribution
#SBATCH --mem=64G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-02:00:00
#SBATCH --partition=power-general-shared-pool
#SBATCH --qos=public

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Edit these before submitting: label -> h5ad with NeuronType / NeuronMaturity.
# Outputs: results/additional_analyses/neuron_distribution/<label>/
declare -A DATASETS=(
  [human_dev]="data/human_dev_without_week_5_cc_annotated_neuron_type.h5ad"
  [cortex]="data/Cortex_EMX1_louvain3_passedQC_PostM_rev1_neuron_type.h5ad"
)

rc=0
for label in "${!DATASETS[@]}"; do
  mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
    python -u modules/neuron_distribution.py --h5ad "${DATASETS[$label]}" --label "$label" || {
      echo "[error] $label failed"; rc=1; }
done

echo "Exit code: $rc"
exit $rc
