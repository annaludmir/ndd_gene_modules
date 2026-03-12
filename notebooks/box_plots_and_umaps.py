import marimo

__generated_with = "0.17.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return


@app.cell
def _():
    import numpy as np
    import seaborn as sns
    import pandas as pd
    import scipy as sp
    import scanpy as sc
    import scipy.sparse as sp
    from matplotlib import pyplot as plt
    from matplotlib.pyplot import rc_context
    return np, pd, plt, sc, sns, sp


@app.cell
def _(sc):
    #upload data
    print('uploading data')
    adata=sc.read_h5ad("/miridan-data/annaludmir/ndd_gene_modules/data/human_dev.h5ad")
    return (adata,)


@app.cell
def _(adata):
    print('filtering data')
    adata_prolif = adata[adata.obs["cell_cycle_score"] > 0.004].copy()
    return


@app.cell
def _(adata):
    print('filtering data')
    adata_rad = adata[adata.obs["CellClass"] == "Radial glia"].copy()
    return (adata_rad,)


@app.cell
def _(adata_rad):
    adata_filtered = adata_rad
    return (adata_filtered,)


@app.cell
def _(adata_filtered, sc):
    print('normalizing data')
    sc.pp.normalize_total(adata_filtered)
    sc.pp.log1p(adata_filtered)
    return


@app.cell
def _(adata, adata_filtered, np, pd, plt, sns, sp):
    genes = ['CEP152','PLK4','STIL','SASS6','CEP135','RTTN','CEP63','CPAP','NINEIN','ASPM','WDR62','PCTN','CDK5RAP2','TUBGCP6','TUBGCP4','NDE1','KIF11','MAP11','CENPE','KNL1','MAP11','CIT','KIF14','ATR','TRAIP','RBBP8','MCPH1','GMNN','ORC1','ORC4','ORC6','CDC45','CDT1','NSMCE2','NCAPH','NACPD2','NCAPD3','NBS1','LIG4','XCCR4','LMNB1','LMNB2','U4ATAC','CENATAC','CDK6','ANKLE2','ZNF335','PHC1','MFSD2A','WDFY3','COPB2','RRP7A']
    region_col = "Region"
    sym_col = "Gene"   # column in adata.var that contains gene symbols

    # symbol -> var_name (with version)
    sym2varname = (
        pd.Series(adata.var_names.values, index=adata.var[sym_col].astype(str))
        .dropna()
        .to_dict()
    )

    filtered_varnames = [sym2varname[g] for g in genes if g in sym2varname]
    missing = [g for g in genes if g not in sym2varname]

    print("Missing symbols:", missing)
    print("Using var_names:", filtered_varnames)

    # Extract per-cell expression (n_cells x n_genes)
    X = adata_filtered[:, filtered_varnames].X
    if sp.issparse(X):
        X = X.toarray()

    # Build a tidy dataframe: one row per cell
    df_base = adata_filtered.obs[[region_col]].copy()

    # Plot one gene at a time
    for gene_sym, varname in zip([g for g in genes if g in sym2varname], filtered_varnames):
        y = np.asarray(X[:, filtered_varnames.index(varname)]).ravel()  # per-cell expression for this gene

        df = df_base.copy()
        df[gene_sym] = y

        plt.figure(figsize=(10, 4))
        sns.boxplot(data=df, x=region_col, y=gene_sym, showfliers=True)
        plt.xticks(rotation=45, ha="right")
        plt.title(f"{gene_sym} expression in Radial glia cells")
        plt.tight_layout()
        plt.show()
    return


@app.cell
def _(adata, sc):
    sc.pp.filter_cells(adata, min_genes=200)
    return


@app.cell
def _(adata, sc):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    adata_highly_variable = adata[:, adata.var.highly_variable]
    return (adata_highly_variable,)


@app.cell
def _(adata_highly_variable, sc):
    sc.pp.pca(adata_highly_variable)
    sc.pp.neighbors(adata_highly_variable)
    sc.tl.umap(adata_highly_variable)
    return


@app.cell
def _(adata_filtered, pd, sc):
    def plot_genes_umap(adata, genes, sym_col="Gene"):
        sym2var = pd.Series(adata.var_names.values,
                            index=adata.var[sym_col].astype(str)).to_dict()

        mapped = [sym2var[g] for g in genes if g in sym2var]

        sc.pl.umap(
            adata_filtered,
            color=mapped,
            title=genes,
            ncols=len(mapped),
            s=5,
            frameon=False,
            vmax="p99"
        )
    return (plot_genes_umap,)


@app.cell
def _(adata, plot_genes_umap):
    filtered_genes = ["FOXG1", "LHX2", "EOMES"]
    plot_genes_umap(adata, filtered_genes)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
