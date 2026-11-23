#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=gsea_diagnostic
#SBATCH --mem=8G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-01:00:00
#SBATCH --partition=power-general-shared-pool
#SBATCH --qos=public

# What drives the GSEA NES of a gene list in each condition: ranked GES
# scores, running ES (weighted p=1 as GSEA ran, and classic p=0), and the
# null ES of random gene sets (NES = ES / mean null ES).
# See modules/gsea_diagnostic.py.
# Outputs: results/additional_analyses/gsea_diagnostic/<label>/

set -euo pipefail

module load mamba/mamba-1.5.8
cd /miridan-data/annaludmir/ndd_gene_modules

BATCH_SUMMARIES=(
  "results/enrichment_results/batch_summary_20261115_final_genes_to_run_on_all_layers.csv=All V3"
  "results/enrichment_results/batch_summary_20261115_final_genes_to_run_on_all_layers_v2.csv=All V2"
  "results/enrichment_results/batch_summary_20261115_final_genes_to_run_on_cortex.csv=Cortex V3"
  "results/enrichment_results/batch_summary_20261115_final_genes_to_run_on_cortex_v2.csv=Cortex V2"
)

ARGS=(
  --batch-summaries "${BATCH_SUMMARIES[@]}"
  --cases config_files/leading_edge_cases_neuron_neuroblast.csv   # gene_list, summary
  # or instead of --cases, every combination of:
  # --gene-lists ASD_SFARI_syndromic --summaries "All V3"
  --column CellClass
  --conditions Neuron Neuroblast
  --n-perm 1000
  --label autism_neuron_vs_neuroblast
)

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/gsea_diagnostic.py "${ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
