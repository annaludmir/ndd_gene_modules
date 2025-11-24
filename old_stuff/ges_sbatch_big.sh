#!/bin/bash
#SBATCH --job-name=search_all_ges
#SBATCH --account=power-miridan-users
#SBATCH --partition=power-general
#SBATCH --output=my_job_%j.out
#SBATCH --error=my_job_%j.err
#SBATCH --mem=800G


module load mamba/mamba1.4.2-environmentally
mamba activate /scratch200/reutj/conda-envs/jupyter-scanpy_new

cd /scratch200/reutj/notebooks/

python -u ges_score_corrected_no_permutations.py ['region_general'] ["b'Forebrain_general'","b'hindbrain_general'"] data_all /scratch200/reutj/data/spec_tables_new/

wait

