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
# #SBATCH --partition=power-general-public-pool
# #SBATCH --qos=public
# # miris partition
#SBATCH --partition=gpu-miridan-pool 
#SBATCH --qos=owner
#SBATCH --gres=gpu:0
# # deprecated
# # SBATCH --account=public-users_v2

set -euo pipefail

#module load mamba/mamba1.4.2-environmentally
module load mamba/mamba-1.5.8
#mamba activate /scratch200/reutj/conda-envs/jupyter-scanpy
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules/modules

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new python -u enrichment_pipeline_for_gene_list.py /miridan-data/annaludmir/ndd_gene_modules/config_files/enrichment_all_layers_config_autism_strong_no_ges_treshold.yaml

echo "Python exit code: $rc"
exit $rc
