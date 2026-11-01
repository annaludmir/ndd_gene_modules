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
#SBATCH --partition=power-general-shared-pool
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

GSEA_SUMMARY_FILE="/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/no_ges_score_threshold/no_ges_theshold_microcephaly_genes All Layers_threshold_0_20260621/data/enrichment_results/GSEA/GSEA_final_summary.csv"
# LEADING_GENE_SUMMARY_FILE="/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/microcephaly/Microcephaly Cell Phase_threshold_1_20260222/data/enrichment_results/GSEA/GSEA_final_summary.csv"

# Options: one per line; comment a line out (or in) to switch it. Unlike
# options continued with "\", a commented line here can't cut the command.
ARGS=(
  --h5ad-path data/human_dev.h5ad
  --gsea-summary-file "$GSEA_SUMMARY_FILE"
  --condition Forebrain_General
  --chemistry v3
  --subfolder-name microcephaly_forebrain_gene_expression_summary
  --wilcoxon-region Forebrain
  --export-all-gene-boxplots
  --not-normalized-boxplots   # also <gene>_not_normalized.png boxplots (raw X)
  # --export-significant-gene-boxplots
  # --wilcoxon-compare-all-region-pairs
  # Dot plots over proliferating cells only (cortex file: "CellCycleFraction > 0.004"):
  # --dotplot-cell-filter "cell_cycle_score > 0.004"
  # --leading-gene-summary-file "$LEADING_GENE_SUMMARY_FILE"
  # --leading-gene-conditions S G2M G1 PostM Non-cycling
  # Raw-region pair mode instead (skips the Wilcoxon / meta-region pipeline):
  # --region-pairs Telencephalon:Midbrain
)

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python modules/gene_expression_summary.py "${ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc

wait
