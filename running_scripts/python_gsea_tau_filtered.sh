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

# Usage: sbatch python_gsea_tau_filtered.sh <gene_list.csv>
# Example: sbatch python_gsea_tau_filtered.sh data/genes/microcephaly_genes.csv
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <gene_list_csv>" >&2
  exit 1
fi

GENE_LIST="$1"

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

CONFIGS=(
  config_files/enrichment_cortex_config_tau_filtered.yaml
  config_files/enrichment_cortex_cell_phase_config_tau_filtered.yaml
  config_files/enrichment_all_layers_config_tau_filtered.yaml
)

for CONFIG in "${CONFIGS[@]}"; do
  echo ""
  echo "=============================="
  echo "Running tau-filtered GSEA: $CONFIG"
  echo "Gene list: $GENE_LIST"
  echo "=============================="
  mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
    python -u modules/search_enrichment_gsea_tau_filtered.py \
    "$CONFIG" \
    --gene-list "$GENE_LIST"
done

rc=$?
echo "Python exit code: $rc"
exit $rc
