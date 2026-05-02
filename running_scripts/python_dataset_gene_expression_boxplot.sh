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
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public
# # miris partition
# #SBATCH --partition=gpu-miridan-pool
# #SBATCH --qos=owner
# #SBATCH --gres=gpu:0
# # deprecated
# # SBATCH --account=public-users_v2

set -euo pipefail

module load mamba/mamba-1.5.8
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new python modules/dataset_analysis_helper.py \
  --task gene_region_boxplot \
  --h5ad-path data/human_dev.h5ad \
  --chemistry v3 \
  --genes FOXG1 NDE1

rc=$?
echo "Python exit code: $rc"
exit $rc

wait
