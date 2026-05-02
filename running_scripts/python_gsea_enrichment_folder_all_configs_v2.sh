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
GENE_LIST_DIR="${1:-}"
RUN_NAME_PREFIX="${2:-}"

if [[ -z "$GENE_LIST_DIR" ]]; then
  echo "Usage: sbatch running_scripts/python_gsea_enrichment_folder_all_configs.sh <gene_list_dir> [run_name_prefix]"
  exit 1
fi

if [[ ! -d "$GENE_LIST_DIR" ]]; then
  echo "Error: gene list directory not found: $GENE_LIST_DIR"
  exit 1
fi

TMP_CONFIG_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_CONFIG_DIR"' EXIT

create_temp_config() {
  local source_config="$1"
  local gene_list_path="$2"
  local run_name="$3"
  local temp_config="$TMP_CONFIG_DIR/$(basename "$source_config" .yaml)_$(basename "$gene_list_path" .csv).yaml"

  sed \
    -e 's|^gene_list_path:.*$|gene_list_path: "'"$gene_list_path"'"|' \
    -e 's|^run_name: ".*"$|run_name: "'"$run_name"'"|' \
    "$source_config" > "$temp_config"

  printf '%s\n' "$temp_config"
}

run_config() {
  local base_config="$1"
  local scope_label="$2"
  local gene_list_path="$3"
  local gene_list_stem="$4"
  local run_name_base="$gene_list_stem"
  local temp_config
  local run_name

  if [[ -n "$RUN_NAME_PREFIX" ]]; then
    run_name_base="${RUN_NAME_PREFIX} ${gene_list_stem}"
  fi

  run_name="${run_name_base} ${scope_label}"
  temp_config="$(create_temp_config "$base_config" "$gene_list_path" "$run_name")"

  echo "Running enrichment for config: $base_config"
  echo "Gene list path: $gene_list_path"
  echo "Run name: $run_name"
  echo "Using temporary config: $temp_config"

  mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
    python -u "$PIPELINE_SCRIPT" "$temp_config"
}

shopt -s nullglob
csv_files=("$GENE_LIST_DIR"/*.csv)
shopt -u nullglob

if [[ ${#csv_files[@]} -eq 0 ]]; then
  echo "Error: no CSV files found in directory: $GENE_LIST_DIR"
  exit 1
fi

echo "Found ${#csv_files[@]} gene list CSV files in: $GENE_LIST_DIR"

for gene_list_path in "${csv_files[@]}"; do
  gene_list_filename="$(basename "$gene_list_path")"
  gene_list_stem="${gene_list_filename%.*}"

  echo
  echo "============================================================"
  echo "Processing gene list: $gene_list_filename"
  echo "============================================================"

  run_config "config_files/enrichment_cortex_config_v2.yaml" "Cortex" "$gene_list_path" "$gene_list_stem"
  run_config "config_files/enrichment_cortex_cell_phase_config_v2.yaml" "Cell Phase" "$gene_list_path" "$gene_list_stem"
  run_config "config_files/enrichment_all_layers_config_v2.yaml" "All Layers" "$gene_list_path" "$gene_list_stem"
done

echo
echo "Completed all enrichment runs for directory: $GENE_LIST_DIR"
