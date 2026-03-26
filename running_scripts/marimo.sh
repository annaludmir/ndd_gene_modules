#!/bin/bash
#SBATCH --job-name=marimo
# #SBATCH --mem=500G
#SBATCH --mem=1G
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

#module load mamba/mamba1.4.2-environmentally
module load mamba/mamba-1.5.8
#mamba activate /scratch200/reutj/conda-envs/jupyter-scanpy
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir

HOST_IP=$(hostname -I | awk '{print $1}')
echo "Using host: $HOST_IP"

marimo edit --headless --port 8080 --host "$HOST_IP"
#marimo edit --headless --port 8080 --host 132.66.114.22
wait
