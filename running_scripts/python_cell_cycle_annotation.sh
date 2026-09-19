#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=cc_annotation
#SBATCH --mem=500G
#SBATCH --cpus-per-task=8
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
<<<<<<< HEAD
#SBATCH --time=0-24:00:00
=======
#SBATCH --time=0-48:00:00
>>>>>>> 52e8d95 (Update notebook and running scripts)
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Edit these before submitting.
MODALITY="atac"   # rna | atac
H5AD_INPUT="/miridan-storage/annaludmir/atac-seq/e88a34d0-d28a-4d10-a8c5-d59f86ba621a.h5ad"
H5AD_OUTPUT="/miridan-storage/annaludmir/atac-seq/e88a34d0-d28a-4d10-a8c5-d59f86ba621a_cc_annotated_fragments.h5ad"

# RNA-specific
SYM_COL="Gene"

# ATAC-specific (used only when MODALITY=atac)
GTF_PATH=""             # leave empty → snapatac2 built-in
GENOME="hg38"
TSS_WINDOW="10000"

# ATAC only: fragments file → score by fragments in promoter windows instead of
# called peaks (adata.X). Leave empty to use the peak-based scoring.
FRAGMENTS_PATH="/miridan-storage/annaludmir/atac-seq/49a73aea-76d5-40a1-b5e4-fa6c0831e164-fragment.tsv.bgz"
FRAG_TEMPDIR="/miridan-storage/annaludmir/atac-seq/tmp_snapatac2"   # scratch for snapatac2 import (large)
FRAG_N_JOBS="${SLURM_CPUS_PER_TASK:-8}"

# Scoring scheme:
#   fraction   — share of the cell's total signal in the gene set. Thresholds
#                are absolute fractions and scale with how many features a set
#                has, so RNA-calibrated values do not transfer to ATAC.
#   background — set mean minus the mean of signal-matched control features,
#                as a z-score. Thresholds are in SDs and comparable across
#                sets. The run also prints a per-set 'signal' ratio; a value
#                near 1 means that set carries no real cell-to-cell structure.
SCORING="background"
N_BACKGROUND_GENES="5000"   # [background + fragments] size of the control pool
N_CTRL_PER_REGION="50"      # [background] controls drawn per set feature
N_BINS="25"                 # [background] signal bins used for matching
SEED="0"

# Optional threshold overrides. Leave empty to use the defaults for $SCORING
# (fraction: 0.004/0.002/0.002/0.03; background: 1.0 SD for all four).
CYCLING_THRESHOLD=""    # e.g. 1.0 under background
G1_THRESHOLD=""
S_THRESHOLD=""
G2M_THRESHOLD=""

EXTRA_ARGS=(--modality "$MODALITY" --scoring "$SCORING")
if [[ "$SCORING" == "background" ]]; then
  EXTRA_ARGS+=(--n-bins "$N_BINS" --n-ctrl-per-region "$N_CTRL_PER_REGION" --seed "$SEED")
fi
if [[ "$MODALITY" == "rna" ]]; then
  EXTRA_ARGS+=(--sym-col "$SYM_COL")
else
  EXTRA_ARGS+=(--genome "$GENOME" --tss-window "$TSS_WINDOW")
  [[ -n "$GTF_PATH" ]] && EXTRA_ARGS+=(--gtf-path "$GTF_PATH")
  if [[ -n "$FRAGMENTS_PATH" ]]; then
    EXTRA_ARGS+=(--fragments "$FRAGMENTS_PATH" --n-jobs "$FRAG_N_JOBS")
    [[ -n "$FRAG_TEMPDIR" ]] && EXTRA_ARGS+=(--tempdir "$FRAG_TEMPDIR")
    if [[ "$SCORING" == "background" ]]; then
      EXTRA_ARGS+=(--n-background-genes "$N_BACKGROUND_GENES")
    fi
  fi
fi

[[ -n "$CYCLING_THRESHOLD" ]] && EXTRA_ARGS+=(--cycling-threshold "$CYCLING_THRESHOLD")
[[ -n "$G1_THRESHOLD"      ]] && EXTRA_ARGS+=(--g1-threshold      "$G1_THRESHOLD")
[[ -n "$S_THRESHOLD"       ]] && EXTRA_ARGS+=(--s-threshold       "$S_THRESHOLD")
[[ -n "$G2M_THRESHOLD"     ]] && EXTRA_ARGS+=(--g2m-threshold     "$G2M_THRESHOLD")

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/cell_cycle_annotation.py \
    --h5ad-input  "$H5AD_INPUT" \
    --h5ad-output "$H5AD_OUTPUT" \
    "${EXTRA_ARGS[@]}"

rc=$?
echo "Python exit code: $rc"
exit $rc
