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

SUMMARY_CSV="results/enrichment_results/microcephaly/microcephaly_genes_all_layers_threshold_1_20260110/data/enrichment_results/GSEA/GSEA_final_summary.csv"
# Output directory is derived automatically from the run name in SUMMARY_CSV.
# Override with --outdir <path> if needed.
# To restrict to specific column conditions, add e.g.: --column-conditions CellCyclePhase Region

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/go_terms_pipeline.py \
  --summary-csv "$SUMMARY_CSV"

rc=$?
echo "Python exit code: $rc"
exit $rc
