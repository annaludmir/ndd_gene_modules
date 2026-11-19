#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=lead_overlap
#SBATCH --mem=8G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-00:30:00
#SBATCH --partition=power-general-shared-pool
#SBATCH --qos=public

# Leading-edge overlap of the autism gene lists between Neuron and Neuroblast
# (Venn per case; ES curves of shared / exclusive genes for the "only" cases).
# See modules/leading_edge_overlap.py.
# Outputs: results/additional_analyses/leading_edge_overlap/<label>/

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
  --cases config_files/leading_edge_cases_neuron_neuroblast.csv   # gene_list, summary, category
  --column CellClass
  --condition-a Neuron
  --condition-b Neuroblast
  --label autism_neuron_vs_neuroblast
  # --es-plots all               # ES curves for the "both" cases too
)

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/leading_edge_overlap.py "${ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
