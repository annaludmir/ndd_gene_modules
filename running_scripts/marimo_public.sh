#!/bin/bash
#SBATCH --job-name=marimo
#SBATCH --mem=50G
#SBATCH --account=miridan-users_v2  
# #SBATCH --account=public-users_v2
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public
# #SBATCH --partition=gpu-miridan-pool 
# #SBATCH --qos=owner
# #SBATCH --gres=gpu:0
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=1-00:00:00

#module load mamba/mamba1.4.2-environmentally
module load mamba/mamba-1.5.8
#mamba activate /scratch200/reutj/conda-envs/jupyter-scanpy
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir

marimo edit --headless --port 8080 --host 132.66.114.22
wait
