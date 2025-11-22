import pandas as pd
import numpy as np
import sys
import os
print('get_gmt script: imports done')

gene_list_df_path=sys.argv[1]

def save_to_gmt(gene_list, gene_set_name, output_file):
    """
    Save a list of gene symbols to a .gmt file for GSEA.
    
    Parameters:
    - gene_list: List of gene symbols.
    - gene_set_name: Name of the gene set.
    - output_file: Path to the output .gmt file.
    """
    with open(output_file, 'w') as f:
        # .gmt format: GeneSetName, Description, Genes...
        gene_list = gene_list_df.gene.to_list()
        genes_joined = "\t".join(gene_list)

        line = f"{gene_set_name}\tDescription\t{genes_joined}\n"
        f.write(line)

if __name__=="__main__":
    gene_list_df=pd.read_csv(gene_list_df_path)
    gene_set_name=gene_list_df_path.split("/")[-1][:-4]
    output_file=gene_list_df_path[:-4]+".gmt"
    if not os.path.exists(output_file):
        print("creating gmt file , path: " , output_file)
        save_to_gmt(gene_list_df, gene_set_name, output_file)
        print('finished gmt file creation')
    else:
        print("gmt file already exists, continue..")
