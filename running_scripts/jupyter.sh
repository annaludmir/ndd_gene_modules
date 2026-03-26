#!/bin/bash
#SBATCH --job-name=jupyter
#SBATCH --mem=500G
#SBATCH --account=miridan-users_v2  
# # SBATCH --account=public-users_v2
# # SBATCH --partition=power-general-public-pool
# #SBATCH --qos=public
#SBATCH --partition=gpu-miridan-pool 
#SBATCH --qos=owner
#SBATCH --gres=gpu:0
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=1-00:00:00

#module load mamba/mamba1.4.2-environmentally
module load mamba/mamba-1.5.8
#mamba activate /scratch200/reutj/conda-envs/jupyter-scanpy
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
echo "1"
cd /miridan-data/annaludmir

PORT=$(shuf -i 8000-8080 -n 1)

#jupyter lab --no-browser --ip=0.0.0.0 --port=8080 2>&1 | grep -m2 $(hostname) | sed -e 's/'"$(hostname)"'/'"$(ip address | grep '132.66' | awk '{print $2}'| awk -F'/' '{print $1}')"'/g' | sed '2!d'
jupyter lab --no-browser --ip=0.0.0.0 --port=8080 2>&1 | grep -m2 $(hostname) | sed -e 's/'"$(hostname)"'/'"$(ip -o -4 addr show | grep '132.66' | awk '{print $4}' | cut -d'/' -f1)"'/g' | sed '2!d'

wait
