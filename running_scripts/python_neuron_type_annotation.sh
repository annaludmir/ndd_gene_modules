#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=neuron_type
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
H5AD_INPUT="/miridan-data/annaludmir/ndd_gene_modules/data/Cortex_EMX1_louvain3_passedQC_PostM_rev1.h5ad"
H5AD_OUTPUT="/miridan-data/annaludmir/ndd_gene_modules/data/Cortex_EMX1_louvain3_passedQC_PostM_rev1_neuron_type.h5ad"
SYM_COL="Gene"
CELL_CLASS_COL="CellClass"   # only cells with this == NEURON_CLASS are classified
NEURON_CLASS="Neuron"
MIN_RATIO="2"                # both marker sets detected: winner needs >= this x the other's UMIs

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python -u modules/neuron_type_annotation.py \
    --h5ad-input  "$H5AD_INPUT" \
    --h5ad-output "$H5AD_OUTPUT" \
    --sym-col "$SYM_COL" \
    --cell-class-col "$CELL_CLASS_COL" \
    --neuron-class "$NEURON_CLASS" \
    --min-ratio "$MIN_RATIO"

rc=$?
echo "Python exit code: $rc"
exit $rc
