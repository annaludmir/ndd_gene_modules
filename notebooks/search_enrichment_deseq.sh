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

python -u search_enrichment_deseq.py ['RG_other','IPC_other','differentiting_proliferating','NPCs_other','G1_other','S_other','G2M_other','SG2M_other','Neuron_other','Glioblast_other','Neuroblast_other','Non-cycling_other'] /scratch200/reutj/data/hsg_gene_lists/hsg_change_seq_gene_set.gmt /scratch200/reutj/data/deseq2/

wait


