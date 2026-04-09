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

#module load mamba/mamba1.4.2-environmentally
module load mamba/mamba-1.5.8
#mamba activate /scratch200/reutj/conda-envs/jupyter-scanpy
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

GSEA_SUMMARY_FILE="results/enrichment_results/Microcephaly All Data (Without Week 5)_threshold_1_20260323/data/enrichment_results/GSEA/GSEA_final_summary.csv"

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new python modules/pseudotime_leading_genes_heatmap.py \
  --gsea-summary-file "$GSEA_SUMMARY_FILE" \
  --condition Forebrain \
  --h5ad-path data/human_dev_without_week_5.h5ad \
  --chemistry v3 \
  --subfolder-name microcephaly_forebrain_pseudotime

rc=$?
echo "Python exit code: $rc"
exit $rc

wait
