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

# Run the enrichment pipeline with ONE config on every gene list CSV found in a
# folder (recursively, including subfolders). For each CSV a temporary copy of
# the config is created with gene_list_path and run_name replaced.
#
# run_name = [<prefix>_]<gene_list>_<config_label>
#   <gene_list>    : CSV path relative to <gene_list_dir>, "/" -> "_"
#                    (<gene_list_dir>/autism/strong.csv -> autism_strong)
#   <config_label> : derived from the config filename unless given as 4th arg
#                    (enrichment_all_layers_cell_phase_config_v2.yaml
#                     -> all_layers_cell_phase_v2)
# The pipeline itself then appends "_threshold_<min_ges_score_threshold>_<date>",
# so the final run dir is e.g. autism_strong_all_layers_cell_phase_v2_threshold_0_20260915
#
# Usage:
#   sbatch running_scripts/python_gsea_enrichment_folder_one_config.sh \
#       <gene_list_dir> <config_yaml> [run_name_prefix] [config_label]

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

PIPELINE_SCRIPT="modules/enrichment_pipeline_for_gene_list.py"
GENE_LIST_DIR="${1:-}"
BASE_CONFIG="${2:-}"
RUN_NAME_PREFIX="${3:-}"
CONFIG_LABEL="${4:-}"

if [[ -z "$GENE_LIST_DIR" || -z "$BASE_CONFIG" ]]; then
  echo "Usage: sbatch running_scripts/python_gsea_enrichment_folder_one_config.sh <gene_list_dir> <config_yaml> [run_name_prefix] [config_label]"
  exit 1
fi

if [[ ! -d "$GENE_LIST_DIR" ]]; then
  echo "Error: gene list directory not found: $GENE_LIST_DIR"
  exit 1
fi

if [[ ! -f "$BASE_CONFIG" ]]; then
  echo "Error: config file not found: $BASE_CONFIG"
  exit 1
fi

TMP_CONFIG_DIR="$(mktemp -d)"
MANIFEST_FILE="$(mktemp)"
trap 'rm -rf "$TMP_CONFIG_DIR" "$MANIFEST_FILE"' EXIT

TODAY="$(date +%Y%m%d)"
OUTPUT_ROOT="$(sed -nE 's|^output_folder:[[:space:]]*"?([^"]*[^"/])/?"?[[:space:]]*$|\1|p' "$BASE_CONFIG")"
OUTPUT_ROOT="${OUTPUT_ROOT:-results/enrichment_results}"

# Default config label: config filename without "enrichment_" prefix, "_config"
# and ".yaml", e.g. enrichment_all_layers_cell_phase_config_v2 -> all_layers_cell_phase_v2
if [[ -z "$CONFIG_LABEL" ]]; then
  CONFIG_LABEL="$(basename "$BASE_CONFIG" .yaml)"
  CONFIG_LABEL="${CONFIG_LABEL#enrichment_}"
  CONFIG_LABEL="${CONFIG_LABEL//_config/}"
fi

# Strip trailing slash so relative paths below are computed cleanly.
GENE_LIST_DIR="${GENE_LIST_DIR%/}"

create_temp_config() {
  local gene_list_path="$1"
  local run_name="$2"
  local safe_name="${run_name// /_}"
  local temp_config="$TMP_CONFIG_DIR/$(basename "$BASE_CONFIG" .yaml)_${safe_name}.yaml"

  sed \
    -e 's|^gene_list_path:.*$|gene_list_path: "'"$gene_list_path"'"|' \
    -e 's|^run_name:.*$|run_name: "'"$run_name"'"|' \
    "$BASE_CONFIG" > "$temp_config"

  printf '%s\n' "$temp_config"
}

run_gene_list() {
  local gene_list_path="$1"
  local rel_path="${gene_list_path#"$GENE_LIST_DIR"/}"
  local rel_stem="${rel_path%.*}"
  # Subfolder path becomes part of the run name: "sub/dir/list" -> "sub_dir_list"
  local run_name="${rel_stem//\//_}_${CONFIG_LABEL}"
  local temp_config

  if [[ -n "$RUN_NAME_PREFIX" ]]; then
    run_name="${RUN_NAME_PREFIX}_${run_name}"
  fi

  temp_config="$(create_temp_config "$gene_list_path" "$run_name")"

  echo "Running enrichment for config: $BASE_CONFIG"
  echo "Gene list path: $gene_list_path"
  echo "Run name: $run_name"
  echo "Using temporary config: $temp_config"

  mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
    python -u "$PIPELINE_SCRIPT" "$temp_config"

  # Record every run dir the pipeline just created (name follows the pattern
  #   "{run_name}_threshold_{min_ges_score_threshold}_{YYYYMMDD}")
  # so build_batch_summary.py can find them at the end.
  shopt -s nullglob
  for d in "$OUTPUT_ROOT/${run_name}_threshold_"*_"${TODAY}"; do
    [[ -d "$d" ]] && printf '%s\n' "$d" >> "$MANIFEST_FILE"
  done
  shopt -u nullglob
}

# Recursive, sorted list of CSV files (subfolders included).
mapfile -t csv_files < <(find "$GENE_LIST_DIR" -type f -name '*.csv' | sort)

if [[ ${#csv_files[@]} -eq 0 ]]; then
  echo "Error: no CSV files found under directory: $GENE_LIST_DIR"
  exit 1
fi

echo "Found ${#csv_files[@]} gene list CSV files under: $GENE_LIST_DIR"
echo "Config: $BASE_CONFIG"

for gene_list_path in "${csv_files[@]}"; do
  echo
  echo "============================================================"
  echo "Processing gene list: ${gene_list_path#"$GENE_LIST_DIR"/}"
  echo "============================================================"
  run_gene_list "$gene_list_path"
done

echo
echo "Completed all enrichment runs for directory: $GENE_LIST_DIR"

# --- Batch summary CSV ----------------------------------------------------
if [[ -s "$MANIFEST_FILE" ]]; then
  if [[ -n "$RUN_NAME_PREFIX" ]]; then
    suffix="${RUN_NAME_PREFIX// /_}"
  else
    suffix="$(basename "$GENE_LIST_DIR")"
  fi
  suffix="${suffix}_${CONFIG_LABEL}"
  SUMMARY_CSV="${OUTPUT_ROOT}/batch_summary_${TODAY}_${suffix}.csv"

  mapfile -t run_dirs < "$MANIFEST_FILE"
  echo
  echo "Building batch summary CSV → $SUMMARY_CSV"
  echo "  (${#run_dirs[@]} run dirs collected)"

  mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
    python -u modules/build_batch_summary.py \
      --run-dirs "${run_dirs[@]}" \
      --output-csv "$SUMMARY_CSV"
else
  echo
  echo "[warn] No run dirs were recorded — skipping batch summary."
fi
