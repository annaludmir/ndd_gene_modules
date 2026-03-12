import numpy as np
import seaborn as sns
import pandas as pd
import scipy as sp
import scanpy as sc
import scipy.sparse as sp
from matplotlib import pyplot as plt
from matplotlib.pyplot import rc_context


def plot_genes_umap(adata, genes, sym_col="Gene"):
    sym2var = pd.Series(adata.var_names.values, index=adata.var[sym_col].astype(str)).to_dict()
    mapped = [sym2var[g] for g in genes if g in sym2var]
    print(mapped)

    sc.pl.umap(
        adata,
        color=mapped,
        title=genes,
        ncols=len(mapped),
        s=5,
        frameon=False,
        vmax="p99",
        save="_my_genes_umap.png"
    )


def main(genes):
    #upload data
    print('uploading data')
    adata=sc.read_h5ad("/miridan-data/annaludmir/ndd_gene_modules/data/full_adata_with_umap.h5ad")
    # sc.pp.filter_cells(adata, min_genes=200)
    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)
    # sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    # adata_highly_variable = adata[:, adata.var.highly_variable]
    # sc.pp.pca(adata_highly_variable)
    # sc.pp.neighbors(adata_highly_variable)
    # sc.tl.umap(adata_highly_variable)

    # adata_highly_variable.write("full_adata_with_umap.h5ad")
    plot_genes_umap(adata, genes)

if __name__ == '__main__':
	filtered_genes = ["FOXG1", "LHX2", "EOMES"]
	main(filtered_genes)
