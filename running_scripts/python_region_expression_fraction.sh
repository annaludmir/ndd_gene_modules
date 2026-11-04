#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=region_fraction
#SBATCH --mem=250G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-04:00:00
#SBATCH --partition=power-general-shared-pool
#SBATCH --qos=public

# Is a gene set expressed in more cells of one region than another because of
# the genes, or because that region's cells detect more genes overall? See
# modules/region_expression_fraction.py.
# Outputs: results/additional_analyses/region_expression_fraction/<subfolder>/

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

GSEA_SUMMARY_FILE="/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/no_ges_score_threshold/no_ges_theshold_microcephaly_genes All Layers_threshold_0_20260621/data/enrichment_results/GSEA/GSEA_final_summary.csv"

# Regions to compare; override at submit time, e.g.
#   REGION_B=Hindbrain sbatch running_scripts/python_region_expression_fraction.sh
REGION_A="${REGION_A:-Forebrain}"   # meta-regions: Forebrain = Forebrain + Telencephalon + Diencephalon,
REGION_B="${REGION_B:-Midbrain}"    # Hindbrain = Hindbrain + Cerebellum + Pons + Medulla
PAIR="$(echo "${REGION_A}_vs_${REGION_B}" | tr '[:upper:]' '[:lower:]')"

# Options: one per line; comment a line out (or in) to switch it.
ARGS=(
  --h5ad-path data/human_dev.h5ad
  --cell-class "Neuronal IPC"
  --region-a "$REGION_A"
  --region-b "$REGION_B"
  --chemistry v3
  --gsea-summary-file "$GSEA_SUMMARY_FILE"
  --condition Forebrain_General   # its Lead_genes are the gene set
  # --gene-list data/genes/microcephaly_genes.csv   # or a gene list instead of leading genes
  --subfolder-name "neuronal_ipc_${PAIR}_microcephaly_leading"
  # --downsample-to 2000          # UMIs per cell for the depth-matched part (default: 10% quantile)
  # --no-collapse                 # raw Region values (e.g. --region-a Telencephalon)
)

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/region_expression_fraction.py "${ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
