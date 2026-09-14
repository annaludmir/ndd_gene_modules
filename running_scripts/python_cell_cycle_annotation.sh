#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=cc_annotation
#SBATCH --mem=250G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-04:00:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Edit these before submitting.
MODALITY="rna"   # rna | atac
H5AD_INPUT="data/human_dev_without_week_5.h5ad"
H5AD_OUTPUT="data/human_dev_without_week_5_cc_annotated.h5ad"

# RNA-specific
SYM_COL="Gene"

# ATAC-specific (used only when MODALITY=atac)
GTF_PATH=""             # leave empty → snapatac2 built-in
GENOME="hg38"
TSS_WINDOW="2000"

EXTRA_ARGS=(--modality "$MODALITY")
if [[ "$MODALITY" == "rna" ]]; then
  EXTRA_ARGS+=(--sym-col "$SYM_COL")
else
  EXTRA_ARGS+=(--genome "$GENOME" --tss-window "$TSS_WINDOW")
  [[ -n "$GTF_PATH" ]] && EXTRA_ARGS+=(--gtf-path "$GTF_PATH")
fi

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/cell_cycle_annotation.py \
    --h5ad-input  "$H5AD_INPUT" \
    --h5ad-output "$H5AD_OUTPUT" \
    "${EXTRA_ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
