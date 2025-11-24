#!/bin/bash
#SBATCH --job-name=search_enrichment
#SBATCH --account=power-miridan-users
#SBATCH --partition=power-general
#SBATCH --output=my_job_%j.out
#SBATCH --error=my_job_%j.err
#SBATCH --mem=200G
#SBATCH --cpus-per-task=10

module load mamba/mamba1.4.2-environmentally
mamba activate /scratch200/reutj/conda-envs/diff_exp

cd /scratch200/reutj/notebooks/

python -u deseq_calculations.py /scratch200/reutj/data/deseq2/pseudobulk_cellclass_data.h5ad /scratch200/reutj/data/deseq2 /scratch200/reutj/data/hsg_gene_lists/har2000_df.csv
 

wait


