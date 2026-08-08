import marimo

__generated_with = "0.17.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return


@app.cell
def _():
    import scanpy as sc
    import numpy as np
    import pandas as pd

    # Load the AnnData metadata / peak matrix
    adata_atac = sc.read_h5ad("/miridan-storage/annaludmir/atac-seq/e88a34d0-d28a-4d10-a8c5-d59f86ba621a.h5ad")

    # Inspect cell annotations (e.g., cell types, brain regions like Telencephalon/Diencephalon)
    print("Cell Metadata Columns:", adata_atac.obs.columns.tolist())
    print(adata_atac.obs['cell_type'].value_counts().head(10))
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
