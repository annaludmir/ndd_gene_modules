#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=gsea_heatmap
#SBATCH --mem=8G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-00:30:00
#SBATCH --partition=power-general-shared-pool
#SBATCH --qos=public

# GSEA NES heatmap (gene lists x conditions) from batch summaries.
# See modules/gsea_heatmap.py.

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Edit these before submitting.
# Batch summaries; add =Label when the same gene list x condition is in several.
BATCH_SUMMARIES=(
  "results/enrichment_results/batch_summary_20261115_final_genes_to_run_on_all_layers.csv=All V3"
  "results/enrichment_results/batch_summary_20261115_final_genes_to_run_on_all_layers_v2.csv=All V2"
  "results/enrichment_results/batch_summary_20261115_final_genes_to_run_on_cortex.csv=Cortex V3"
  "results/enrichment_results/batch_summary_20261115_final_genes_to_run_on_cortex_v2.csv=Cortex V2"
)
# Gene lists, top to bottom: file name (with or without .csv), optionally =Row label.
# Leave empty to show every gene list in the summaries.
GENE_LISTS=(
  "ASD_HC"
  "ASD_SFARI"
  "autism_strong"
  "ASD_SFARI_2"
  "ASD_SFARI_3"
  "ASD_SFARI_syndromic"
  "ASD-DD"
  "ASD-predominant"
)
# Only these columns (column_condition_title), in this order; empty = all.
COLUMNS=(
  CellClass
  NeuronLineage
  NeuronSubtype
  NeuronMaturity
  Region
  RegionsGeneral
)
OUTPUT="results/additional_analyses/gsea_heatmap/gsea_heatmap.png"

ARGS=(--batch-summaries "${BATCH_SUMMARIES[@]}" --output "$OUTPUT")
[[ ${#GENE_LISTS[@]} -gt 0 ]] && ARGS+=(--gene-lists "${GENE_LISTS[@]}")
[[ ${#COLUMNS[@]}    -gt 0 ]] && ARGS+=(--columns "${COLUMNS[@]}")
# ARGS+=(--nes-cutoff 1.5 --trend-q 0.25 --vmax 3 --title "My title")

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/gsea_heatmap.py "${ARGS[@]}"
