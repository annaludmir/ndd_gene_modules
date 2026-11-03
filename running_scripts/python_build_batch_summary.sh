#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=batch_summary
#SBATCH --mem=16G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-01:00:00
# # public
#SBATCH --partition=power-general-shared-pool
#SBATCH --qos=public

# Rebuild the GSEA batch summary CSV for finished runs, without rerunning GSEA
# (e.g. after ai-adata-miner script 29 rewrote the gene lists' index.csv, to add
# the gene list descriptions).
#
# Usage:
#   sbatch running_scripts/python_build_batch_summary.sh <run_glob> <output_csv> [index.csv ...]
#
#   <run_glob>   : run dirs under OUTPUT_ROOT, e.g. '*_threshold_*_20261101'
#                  (quote it so the shell does not expand it at submit time)
#   <output_csv> : where to write the summary
#   [index.csv]  : gene list index files (script 29), only needed when the lists
#                  GSEA read are not next to their index.csv
#
# OUTPUT_ROOT defaults to results/enrichment_results/spectra; override with
#   OUTPUT_ROOT=... sbatch running_scripts/python_build_batch_summary.sh ...

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

RUN_GLOB="${1:-}"
OUTPUT_CSV="${2:-}"
if [[ -z "$RUN_GLOB" || -z "$OUTPUT_CSV" ]]; then
  echo "Usage: sbatch running_scripts/python_build_batch_summary.sh <run_glob> <output_csv> [index.csv ...]"
  exit 1
fi
shift 2
OUTPUT_ROOT="${OUTPUT_ROOT:-results/enrichment_results/spectra}"

INDEX_ARGS=()
if [[ $# -gt 0 ]]; then
  INDEX_ARGS=(--gene-list-index "$@")
fi

python -u modules/build_batch_summary.py \
  --base-folder "$OUTPUT_ROOT" \
  --run-glob "$RUN_GLOB" \
  --output-csv "$OUTPUT_CSV" \
  ${INDEX_ARGS[@]+"${INDEX_ARGS[@]}"}
