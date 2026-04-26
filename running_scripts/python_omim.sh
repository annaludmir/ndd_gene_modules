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

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

SUMMARY_CSV="results/enrichment_results/AH MPRA NPC True Positive LogFC All Layers_threshold_1_20260412/data/enrichment_results/GSEA/GSEA_final_summary.csv"
# Output directory is derived automatically from the run name in SUMMARY_CSV.
# Override with --outdir <path> if needed.
# To restrict to specific column conditions, add e.g.: --column-conditions CellCyclePhase Region
# To use different OMIM gene sets, add e.g.: --gene-sets OMIM_Disease OMIM_Expanded

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/omim_pipeline.py \
  --summary-csv "$SUMMARY_CSV"
  --column-conditions Region \

rc=$?
echo "Python exit code: $rc"
exit $rc
