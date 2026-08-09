#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=expr_fraction_age_region
#SBATCH --mem=250G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-06:00:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Edit these before submitting.
H5AD_PATH="data/human_dev_without_week_5.h5ad"
GENE_LIST="data/genes_for_analyses/microcephaly_forebrain_leading_genes"
SUBFOLDER_NAME="microcephaly_forebrain_leading_genes"
CHEMISTRY="v3"
PLOT_METRIC="Proliferating_Cycling"   # Proliferating_Cycling | Differentiating_Cycling |
                                      # Proliferating_NonCycling | Differentiating_NonCycling

# Expression criterion. Leave both empty for "any gene" (≥1).
# --min-gene-fraction 0.5 → cell must express at least half of the found genes.
# --min-gene-count 5    → cell must express at least 5 target genes (overrides fraction).
MIN_GENE_FRACTION="0.5"
MIN_GENE_COUNT=""

# Optional overrides (keep defaults unless the all-layers h5ad uses different names).
CELL_CLASS_COL="CellClass"
CELL_CYCLE_SCORE_COL="cell_cycle_score"
CELL_CYCLE_THRESHOLD="0.004"
AGE_COL="Age"
REGION_COL="Region"
SYM_COL="Gene"

EXTRA_ARGS=()
[[ -n "$MIN_GENE_FRACTION" ]] && EXTRA_ARGS+=(--min-gene-fraction "$MIN_GENE_FRACTION")
[[ -n "$MIN_GENE_COUNT"    ]] && EXTRA_ARGS+=(--min-gene-count    "$MIN_GENE_COUNT")

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/dataset_analysis_helper.py \
    --task expression_fraction_by_age_region_cellcycle \
    --h5ad-path "$H5AD_PATH" \
    --gene-list "$GENE_LIST" \
    --subfolder-name "$SUBFOLDER_NAME" \
    --chemistry "$CHEMISTRY" \
    --plot-metric "$PLOT_METRIC" \
    --cell-class-col "$CELL_CLASS_COL" \
    --cell-cycle-score-col "$CELL_CYCLE_SCORE_COL" \
    --cell-cycle-threshold "$CELL_CYCLE_THRESHOLD" \
    --age-col "$AGE_COL" \
    --region-col "$REGION_COL" \
    --sym-col "$SYM_COL" \
    --exclude-regions Brain Head \
    "${EXTRA_ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
