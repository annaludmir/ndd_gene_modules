#!/bin/bash
# mail_job_summary.sh <jobid> [email]
#
# Wait for a Slurm job to leave the queue, then mail its summary file.
#
# RUN THIS ON THE LOGIN NODE. Compute nodes accept mail (sendmail exits 0) but
# silently fail to relay it, so the job itself cannot send its own report; the
# login node can. The job writes <jobid>_summary.txt and this mails it.
#
#   JOBID=$(sbatch --parsable running_scripts/python_cell_cycle_annotation.sh)
#   nohup running_scripts/mail_job_summary.sh "$JOBID" >/dev/null 2>&1 &
#
# Tunable via env: POLL_SECONDS (default 60), MAX_HOURS (default 60).

set -uo pipefail   # not -e: squeue returns non-zero once the job is gone

JOBID="${1:?usage: mail_job_summary.sh <jobid> [email]}"
MAILTO="${2:-annaludmir@mail.tau.ac.il}"
OUTDIR="/miridan-data/annaludmir/jobs_output"
SUMMARY="${OUTDIR}/${JOBID}_summary.txt"
OUTFILE="${OUTDIR}/${JOBID}.out"
POLL="${POLL_SECONDS:-60}"
MAX_HOURS="${MAX_HOURS:-60}"

deadline=$(( SECONDS + MAX_HOURS * 3600 ))
while squeue -h -j "$JOBID" -o %T 2>/dev/null | grep -q .; do
  if (( SECONDS > deadline )); then
    echo "gave up waiting for job ${JOBID} after ${MAX_HOURS}h" >&2
    break
  fi
  sleep "$POLL"
done

sleep 10   # let the job's EXIT trap finish writing the summary

state=$(sacct -n -X -j "$JOBID" -o State%20 2>/dev/null | head -1 | tr -d ' ')

{
  echo "To: ${MAILTO}"
  echo "Subject: [cc_annotation] job ${JOBID} ${state:-finished}"
  echo "From: ${MAILTO}"
  echo
  if [[ -f "$SUMMARY" ]]; then
    cat "$SUMMARY"
  elif [[ -f "$OUTFILE" ]]; then
    echo "(no summary file — the job died before its exit trap ran)"
    echo "(last 40 lines of ${OUTFILE})"
    echo
    tail -n 40 "$OUTFILE"
  else
    echo "(no output found for job ${JOBID} in ${OUTDIR})"
  fi
} | /usr/sbin/sendmail -t
