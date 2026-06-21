#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=fix_pyscenic
#SBATCH --mem=4G
#SBATCH --account=miridan-users_v2
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=0-00:10:00
#SBATCH --qos=public

set -euo pipefail

NDD_ROOT="/miridan-data/annaludmir/ndd_gene_modules"

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd "$NDD_ROOT"

echo "Patching pyscenic for NumPy >= 1.24 and Dask compatibility..."

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new \
  python running_scripts/fix_pyscenic_numpy_compat.py

echo "Done."
