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
# # deprecated
# # SBATCH --account=public-users_v2

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

PIPELINE_SCRIPT="modules/enrichment_pipeline_for_gene_list.py"
GENE_LIST_PATH="${1:-data/genes/AH_MPRA_NPC_genes.csv}"
RUN_NAME_BASE_INPUT="${2:-}"

GENE_LIST_FILENAME="$(basename "$GENE_LIST_PATH")"
GENE_LIST_STEM="${GENE_LIST_FILENAME%.*}"
RUN_NAME_BASE="${RUN_NAME_BASE_INPUT:-$GENE_LIST_STEM}"

TMP_CONFIG_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_CONFIG_DIR"' EXIT

create_temp_config() {
  local source_config="$1"
  local run_name="$2"
  local temp_config="$TMP_CONFIG_DIR/$(basename "$source_config")"

  sed \
    -e 's|^gene_list_path:.*$|gene_list_path: "'"$GENE_LIST_PATH"'"|' \
    -e 's|^run_name: ".*"$|run_name: "'"$run_name"'"|' \
    "$source_config" > "$temp_config"

  printf '%s\n' "$temp_config"
}

run_config() {
  local base_config="$1"
  local scope_label="$2"
  local temp_config
  local run_name="${RUN_NAME_BASE} ${scope_label}"

  temp_config="$(create_temp_config "$base_config" "$run_name")"
  echo "Running enrichment for config: $base_config"
  echo "Gene list path: $GENE_LIST_PATH"
  echo "Run name: $run_name"
  echo "Using temporary config: $temp_config"
  mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
    python -u "$PIPELINE_SCRIPT" "$temp_config"
}

run_config "config_files/enrichment_cortex_config.yaml" "Cortex"
run_config "config_files/enrichment_cortex_cell_phase_config.yaml" "Cell Phase"
run_config "config_files/enrichment_all_layers_config.yaml" "All Layers"
