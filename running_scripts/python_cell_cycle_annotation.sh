#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=cc_annotation
#SBATCH --mem=500G
#SBATCH --cpus-per-task=8
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-48:00:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

set -euo pipefail

# ── End-of-job report ─────────────────────────────────────────────────────────
# Slurm's own --mail-type mail carries only the job status, never the output.
# This writes the interesting part of the log to a summary file and, when the
# compute node can relay mail, sends it too. The summary file is written either
# way, so a silently-dropped mail does not lose the report.
MAILTO="annaludmir@mail.tau.ac.il"
# Compute nodes accept mail (sendmail exits 0) but do not relay it, so this is
# off: the summary file is the reliable copy. Flip to "true" only if a test
# message from a compute node actually arrives.
EMAIL_OUTPUT="false"
JOBID="${SLURM_JOB_ID:-local}"
OUTFILE="/miridan-data/annaludmir/jobs_output/${JOBID}.out"
SUMMARY="/miridan-data/annaludmir/jobs_output/${JOBID}_summary.txt"

send_report() {
  rc=$?                          # must stay first: anything else overwrites it
  {
    echo "job ${JOBID}  host=$(hostname)  exit=${rc}"
    echo "log: ${OUTFILE}"
    echo
    if [[ -f "$OUTFILE" ]] && grep -q 'Per-cell signal distribution' "$OUTFILE"; then
      # the diagnostics, not a blind tail — a fixed tail would cut the
      # cell-type table short depending on how many cell types print
      sed -n '/Per-cell signal distribution/,$p' "$OUTFILE" | head -100
    elif [[ -f "$OUTFILE" ]]; then
      echo "(diagnostics not reached — last 40 lines of the log)"
      echo
      tail -n 40 "$OUTFILE"
    else
      echo "(no log file at ${OUTFILE})"
    fi
  } > "$SUMMARY" 2>/dev/null || true

  if [[ "$EMAIL_OUTPUT" == "true" && -x /usr/sbin/sendmail ]]; then
    {
      echo "To: ${MAILTO}"
      echo "Subject: [cc_annotation] job ${JOBID} exit=${rc}"
      echo "From: ${MAILTO}"
      echo
      cat "$SUMMARY"
    } | /usr/sbin/sendmail -t || true
  fi
}
trap send_report EXIT

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Edit these before submitting.
MODALITY="atac"   # rna | atac
H5AD_INPUT="/miridan-storage/annaludmir/atac-seq/e88a34d0-d28a-4d10-a8c5-d59f86ba621a.h5ad"
H5AD_OUTPUT="/miridan-storage/annaludmir/atac-seq/e88a34d0-d28a-4d10-a8c5-d59f86ba621a_cc_replication.h5ad"

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

# ATAC signal:
#   promoter    — score accessibility over cell-cycle gene promoters (gene sets).
#   replication — ignore gene sets. S phase from megabase-scale coverage
#                 overdispersion (a replicating cell has 2 copies of what it has
#                 copied and 1 of the rest), G2M from DNA content. Needs
#                 FRAGMENTS_PATH. Set CELL_TYPE_COL: the per-cell-type baseline
#                 matters for correctness, not just for reporting.
ATAC_SIGNAL="replication"
CELL_TYPE_COL="cell_type"
BIN_SIZE=""             # empty → chosen from the data's median depth
S_DISPERSION_THRESHOLD=""   # default 0.02 overdispersion above the type baseline
S_Z_THRESHOLD=""            # default 3 sd of the Poisson null
G2M_DNA_THRESHOLD=""        # default 0.6 log2 DNA above the type baseline
FORCE_G2M="false"           # G2M is auto-disabled when capture variation hides it

# Scoring scheme (ATAC_SIGNAL=promoter, or MODALITY=rna):
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
  EXTRA_ARGS+=(--genome "$GENOME" --tss-window "$TSS_WINDOW" --atac-signal "$ATAC_SIGNAL")
  [[ -n "$GTF_PATH" ]] && EXTRA_ARGS+=(--gtf-path "$GTF_PATH")
  if [[ "$ATAC_SIGNAL" == "replication" ]]; then
    [[ -n "$CELL_TYPE_COL"          ]] && EXTRA_ARGS+=(--cell-type-col "$CELL_TYPE_COL")
    [[ -n "$BIN_SIZE"               ]] && EXTRA_ARGS+=(--bin-size "$BIN_SIZE")
    [[ -n "$S_DISPERSION_THRESHOLD" ]] && EXTRA_ARGS+=(--s-dispersion-threshold "$S_DISPERSION_THRESHOLD")
    [[ -n "$S_Z_THRESHOLD"          ]] && EXTRA_ARGS+=(--s-z-threshold "$S_Z_THRESHOLD")
    [[ -n "$G2M_DNA_THRESHOLD"      ]] && EXTRA_ARGS+=(--g2m-dna-threshold "$G2M_DNA_THRESHOLD")
    [[ "$FORCE_G2M" == "true"       ]] && EXTRA_ARGS+=(--force-g2m)
  fi
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
