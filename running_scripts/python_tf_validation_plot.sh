#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=tf_validation_plot
#SBATCH --mem=8G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-00:30:00
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

set -euo pipefail

# Plots already-existing tf_validation CSVs. New tf_validation_report runs
# draw these figures themselves (see its --plot-style flag).

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

# Edit these before submitting.
# A single tf_validation subfolder, or results/tf_validation to plot them all.
VALIDATION_DIR="results/tf_validation"
PLOT_STYLE="connected"   # connected | hub_spoke
TOP_N_TFS=12

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/tf_validation_plot.py \
    --validation-dir "$VALIDATION_DIR" \
    --plot-style "$PLOT_STYLE" \
    --top-n-tfs "$TOP_N_TFS"

rc=$?
echo "Python exit code: $rc"
exit $rc
