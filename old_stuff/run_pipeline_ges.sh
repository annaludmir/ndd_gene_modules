#!/bin/bash
#SBATCH --job-name=search_enrichment
#SBATCH --account=power-miridan-users
#SBATCH --partition=power-general
#SBATCH --output=my_job_%j.out
#SBATCH --error=my_job_%j.err
#SBATCH --mem=300G
#SBATCH --cpus-per-task=10

module load mamba/mamba1.4.2-environmentally
mamba activate /scratch200/reutj/conda-envs/jupyter-scanpy_new

cd /scratch200/reutj/notebooks/

python -u gene_list_pipeline.py /scratch200/reutj/data/hsg_gene_lists/har_genes_neurons_df.csv /scratch200/reutj/enrichment_results/ ges_enrichment cortex /scratch200/reutj/data/hsg_gene_lists/   ['radialglia','proliferating','differentiating','IPC','NPCs','G1','S','G2M','PostM','Neuron','Glioblast','Neuroblast','differentiating_non_cycling','proliferating_non_cycling'] 

wait

