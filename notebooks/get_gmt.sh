#!/bin/bash
#SBATCH --job-name=get_gmt
#SBATCH --account=power-miridan-users
#SBATCH --partition=power-general
#SBATCH --output=my_job_%j.out
#SBATCH --error=my_job_%j.err
#SBATCH --mem=300G
#SBATCH --cpus-per-task=10

module load mamba/mamba1.4.2-environmentally
mamba activate /scratch200/reutj/conda-envs/jupyter-scanpy_new

cd /scratch200/reutj/notebooks/

python -u get_gmt.py /scratch200/reutj/data/hsg_gene_lists/autism_genes.csv 

wait

