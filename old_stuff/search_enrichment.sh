#!/bin/bash
#SBATCH --job-name=search_enrichment
#SBATCH --account=power-miridan-users
#SBATCH --partition=power-general
#SBATCH --output=my_job_%j.out
#SBATCH --error=my_job_%j.err
#SBATCH --mem=200G
#SBATCH --cpus-per-task=10

module load mamba/mamba1.4.2-environmentally
mamba activate /scratch200/reutj/conda-envs/jupyter-scanpy

cd /scratch200/reutj/notebooks/

python -u search_enrichment.py ['radialglia','proliferating','differentiating','IPC','NPCs','G1','S','G2M','PostM','Neuron','Glioblast','Neuroblast','differentiating_non_cycling','proliferating_non_cycling'] /scratch200/reutj/data/hsg_gene_lists/har_associated.gmt /scratch200/reutj/data/ges_results/

wait

