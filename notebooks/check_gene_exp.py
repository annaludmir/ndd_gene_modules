import marimo

__generated_with = "0.17.8"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    from __future__ import annotations

    import scanpy as sc
    import pandas as pd
    import numpy as np
    from matplotlib.pyplot import rc_context
    return np, pd, rc_context, sc


@app.cell
def _(pd):
    #upload_genes
    id_genes=pd.read_csv("/miridan-data/annaludmir/data/genes/sysndd_id_abnormal_facial_shape_abnormal_heart_morphology.csv")
    return (id_genes,)


@app.cell
def _(id_genes):
    id_genes
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Playing Around With The Data
    """)
    return


@app.cell
def _(adata):
    adata.obs
    return


@app.cell
def _(adata):
    adata_f = adata[
        adata.obs["Region"].isin(['Midbrain'])
        & (adata.obs["cell_cycle_score"] > 0.004)
    ].copy()
    return (adata_f,)


@app.cell
def _(adata_f):
    adata_f.obs
    return


@app.cell
def _(adata, sc):
    adata_in_memory = adata.to_memory()
    sc.pp.filter_cells(adata_in_memory, min_genes=200)
    sc.pp.filter_genes(adata_in_memory, min_cells=3)
    return


@app.cell
def _(adata):
    genes_adata=adata.var_names.tolist()
    return (genes_adata,)


@app.cell
def _(adata):
    ens_adata=adata.var["Accession"].to_list()
    return (ens_adata,)


@app.cell
def _(adata, ens_adata, np, random_genes):
    # Initialize results dictionary
    results_ens = {
        "gene": [],
        "num_cells_expressed": [],
        "mean_expression": [],
        "CellClass_expressed_in": [],
        "proliferation_expressed_in": [],
        "CellCycleStatus_expressed_in": []
    }

    # Loop through each gene in the list
    random_genes_ens = list(random_genes["ens"])

    for gene in random_genes_ens:
        if gene in ens_adata:  # <- safer than 'in ens_adata'
            print(f"Processing {gene}")
            ind= ens_adata.index(gene)
            print(ind)
            # Get expression data for the gene from counts layer
            gene_counts = adata.X[:,ind].toarray().flatten()

            # Find which cells express the gene
            expressed_cells = gene_counts > 0
            num_cells_expressed = np.sum(expressed_cells)
            mean_expression = gene_counts[expressed_cells].mean() if num_cells_expressed > 0 else 0

            # Get metadata for expressed cells
            if num_cells_expressed > 0:
                cell_classes = adata.obs.loc[expressed_cells, "CellClass"].unique().tolist()
                proliferations = adata.obs.loc[expressed_cells, "proliferation"].unique().tolist()
                cell_cycle_statuses = adata.obs.loc[expressed_cells, "CellCycleStatus"].unique().tolist()
            else:
                cell_classes, proliferations, cell_cycle_statuses = [], [], []

            # Store results
            results_ens["gene"].append(gene)
            results_ens["num_cells_expressed"].append(num_cells_expressed)
            results_ens["mean_expression"].append(mean_expression)
            results_ens["CellClass_expressed_in"].append(cell_classes)
            results_ens["proliferation_expressed_in"].append(proliferations)
            results_ens["CellCycleStatus_expressed_in"].append(cell_cycle_statuses)

        else:
            print(f"{gene} not found in ens_adata")
    return (results_ens,)


@app.cell
def _(pd, results_ens):
    genes_exp = pd.DataFrame(results_ens)
    genes_exp
    return (genes_exp,)


@app.cell
def _(genes_exp):
    genes_exp.sort_values(by="mean_expression")
    return


@app.cell
def _(genes_exp, pd, random_genes):
    genes_exp_1 = genes_exp.rename(columns={'gene': 'ens'})
    genes_exp_1 = pd.merge(random_genes, genes_exp_1, on='ens').iloc[:, 1:]
    return (genes_exp_1,)


@app.cell
def _(genes_exp_1):
    genes_exp_1.to_csv('genes_exp_in_adata_cortex.csv')
    return


@app.cell
def _(adata):
    adata.obs
    return


@app.cell
def _(sc):
    sc.set_figure_params(dpi=100, color_map="viridis")
    sc.settings.verbosity = 0
    sc.logging.print_header()
    return


@app.cell
def _():
    PECENTAGE_OF_CELLS_EXRESSING_GENE = 0.30
    return (PECENTAGE_OF_CELLS_EXRESSING_GENE,)


@app.cell
def _(
    PECENTAGE_OF_CELLS_EXRESSING_GENE,
    adata,
    color_vars,
    genes_adata,
    np,
    rc_context,
    sc,
):
    import scipy.sparse as sp

    # Extract expression matrix only for these genes
    X = adata[:, genes_adata].X

    # Convert to dense if sparse
    if sp.issparse(X):
        X = X.toarray()

    # Compute fraction of cells expressing each gene (>0)
    expr_frac = (X > 0).mean(axis=0)

    # Make sure it's a flat 1D numpy array
    expr_frac = np.ravel(genes_adata)

    # Filter: keep genes expressed in ≥5% of cells
    filtered_genes = [g for g, frac in zip(genes_adata, expr_frac) if frac >= PECENTAGE_OF_CELLS_EXRESSING_GENE]
    # remove duplicates
    filtered_genes = list(set(filtered_genes))

    print(f"Kept {len(filtered_genes)} of {len(color_vars)} genes (≥{PECENTAGE_OF_CELLS_EXRESSING_GENE*100}% of cells)")

    # Plot in batches of 16
    for i in range(0, len(filtered_genes), 16):
        subset = filtered_genes[i:i+16]
        print(f"Plotting genes {i+1}–{i+len(subset)} of {len(filtered_genes)}")

        with rc_context({"figure.figsize": (3, 3)}):
            sc.pl.umap(
                adata,
                color=subset,
                s=1,
                frameon=False,
                ncols=4,
                vmax="p99"
            )
    return expr_frac, sp


@app.cell
def _(color_vars, expr_frac):
    PECENTAGE_OF_CELLS_EXRESSING_GENE_1 = 0.5
    filtered_genes_1 = [g for g, frac in zip(color_vars, expr_frac) if frac >= PECENTAGE_OF_CELLS_EXRESSING_GENE_1]
    filtered_genes_1 = list(set(filtered_genes_1))
    return (filtered_genes_1,)


@app.cell
def _(adata, filtered_genes_1, sc):
    adata_in_memory_1 = adata.to_memory()
    top_genes = adata_in_memory_1[:, filtered_genes_1].to_df().var().sort_values(ascending=False).head(30).index
    sc.tl.dendrogram(adata_in_memory_1, groupby='CellClass')
    sc.pl.heatmap(adata_in_memory_1, var_names=top_genes, groupby='CellClass', standard_scale='var', swap_axes=True, cmap='viridis', dendrogram=True)
    return (adata_in_memory_1,)


@app.cell
def _(adata_in_memory_1, filtered_genes_1, sc):
    adata_in_memory_1.obs['Age'] = adata_in_memory_1.obs['Age'].astype('category')
    sc.tl.dendrogram(adata_in_memory_1, groupby='Age')
    sc.pl.heatmap(adata_in_memory_1, var_names=filtered_genes_1, groupby='Age', standard_scale='var', swap_axes=True, cmap='viridis', dendrogram=True)
    return


@app.cell
def _(adata_in_memory_1, filtered_genes_1, sc):
    adata_in_memory_2 = adata_in_memory_1[adata_in_memory_1.obs['CellCyclePhase'].notna(), :].copy()
    adata_in_memory_2.obs['CellCyclePhase'] = adata_in_memory_2.obs['CellCyclePhase'].astype('category').cat.remove_unused_categories()
    print(adata_in_memory_2.obs['CellCyclePhase'].value_counts())
    if 'dendrogram_CellCyclePhase' in adata_in_memory_2.uns:
        del adata_in_memory_2.uns['dendrogram_CellCyclePhase']
    sc.tl.dendrogram(adata_in_memory_2, groupby='CellCyclePhase')
    sc.pl.heatmap(adata_in_memory_2, var_names=filtered_genes_1, groupby='CellCyclePhase', standard_scale='var', swap_axes=True, cmap='viridis', dendrogram=True)
    return (adata_in_memory_2,)


@app.cell
def _(adata_in_memory_2, filtered_genes_1, sc):
    sc.tl.dendrogram(adata_in_memory_2, groupby='CellClass')
    sc.pl.dotplot(adata_in_memory_2, var_names=filtered_genes_1, groupby='CellClass', standard_scale='var', swap_axes=True, cmap='viridis', dendrogram=True)
    return


@app.cell
def _(adata_in_memory_2, filtered_genes_1, sc):
    adata_in_memory_2.obs['Age'] = adata_in_memory_2.obs['Age'].astype('category')
    sc.tl.dendrogram(adata_in_memory_2, groupby='Age')
    sc.pl.dotplot(adata_in_memory_2, var_names=filtered_genes_1, groupby='Age', standard_scale='var', swap_axes=True, cmap='viridis', dendrogram=True)
    return


@app.cell
def _(adata_in_memory_2, filtered_genes_1, sc):
    sc.tl.dendrogram(adata_in_memory_2, groupby='CellCyclePhase')
    sc.pl.dotplot(adata_in_memory_2, var_names=filtered_genes_1, groupby='CellCyclePhase', standard_scale='var', swap_axes=True, cmap='viridis', dendrogram=True)
    return


@app.cell
def _(mo):
    mo.md(r"""
    #Running GES Score and Pipeline
    """)
    return


@app.cell
def _(mo):
    import subprocess

    mo.md(
        """
        # Running Bash Commands with `subprocess`

        We are running the GES score calculation
        """
    )

    # Execute a Bash command and capture its output
    try:
        result = subprocess.run(
            [
                "python",
                "-u",
        "./notebooks/ges_score_corrected_no_permutations.py",
                "['region_general']",
                "['radialglia','Cerebellum','Diencephalon','Fibroblast','Glioblast','Forebrain_general','hindbrain_general','Immune','Midbrain_general','Neural_crest','IPC','Neuroblast','Neuron','Oligo','Forebrain','Hindbrain','Placodes','Pons','Medulla','Telencephalon','Vascular']",
                "data_all",
                "/miridan-data/annaludmir/data/ges_results/"
            ],
            capture_output=True,
            text=True,
            check=True
        )

    except subprocess.CalledProcessError as e:
        print("STDOUT:\n", e.stdout)
        print("\nSTDERR:\n", e.stderr)
        raise
    return


app._unparsable_cell(
    r"""
    !python -u gene_list_pipeline.py /miridan-data/annaludmir/data/genes/sysndd_id_abnormal_facial_shape_abnormal_heart_morphology_for_pipeline.csv /miridan-data/annaludmir/data/enrichment_results/ ges_enrichment data_all /miridan-data/annaludmir/data/genes/ ['radialglia','Cerebellum','Diencephalon','Fibroblast','Glioblast','Forebrain_general','hindbrain_general','Immune','Midbrain_general','Neural_crest','IPC','Neuroblast','Neuron','Oligo','Forebrain','Hindbrain','Placodes','Pons','Medulla','Telencephalon','Vascular']
    """,
    name="_"
)


@app.cell
def _(mo):
    mo.md(r"""
    # Analyzing Results
    """)
    return


@app.cell
def _(pd):
    full_enrichment_results = pd.read_csv('/miridan-data/annaludmir/data/enrichment_results/sysndd_id_abnormal_facial_shape_abnormal_heart_morphology_for_pipeline/ges_results_above1/cortex/sysndd_id_abnormal_facial_shape_abnormal_heart_morphology_for_pipeline/full_results.csv')
    return (full_enrichment_results,)


@app.cell
def _(full_enrichment_results):
    full_enrichment_results
    return


@app.cell
def _(full_enrichment_results):
    npc_genes_str = full_enrichment_results.loc[
        full_enrichment_results["condition"] == "NPCs",
        "lead_genes"
    ].iloc[0]
    npc_genes = set(npc_genes_str.split(';'))
    return (npc_genes,)


@app.cell
def _(full_enrichment_results):
    g1_genes_str = full_enrichment_results.loc[
        full_enrichment_results["condition"] == "G1",
        "lead_genes"
    ].iloc[0]
    g1_genes = set(g1_genes_str.split(';'))
    return (g1_genes,)


@app.cell
def _(g1_genes, npc_genes):
    shared_genes = npc_genes & g1_genes
    print(f"{len(shared_genes)} Shared genes in Neuroblastss and G1: {shared_genes}")

    only_npcs = npc_genes - g1_genes
    print(f"'\n{len(only_npcs)} Only in Neuroblasts:", only_npcs)

    only_g1 = g1_genes - npc_genes
    print(f"\n{len(only_g1)} Only in G1:", only_g1)
    return


@app.cell
def _(pd):
    ges_spec_v3_G1 = pd.read_csv('/miridan-data/annaludmir/data/ges_results/ges_spec_v3_G1.csv')
    return (ges_spec_v3_G1,)


@app.cell
def _(ges_spec_v3_G1):
    ges_spec_v3_G1
    return


@app.cell
def _(g1_genes, ges_spec_v3_G1):
    ges_spec_v3_G1_enriched_genes = ges_spec_v3_G1[ges_spec_v3_G1["gene"].isin(g1_genes)][["gene", "ges_score"]]
    ges_spec_v3_G1_enriched_genes.sort_values("ges_score", ascending=False)
    return


@app.cell
def _(g1_genes, ges_spec_v3_G1):
    low_ges_gene_names = ges_spec_v3_G1.loc[
        (ges_spec_v3_G1["gene"].isin(g1_genes)) & (ges_spec_v3_G1["ges_score"] < 1),
        "gene"
    ].tolist()

    low_ges_gene_names
    return


@app.cell
def _(pd):
    ges_spec_v3_NPCs = pd.read_csv('/miridan-data/annaludmir/data/ges_results/ges_spec_v3_NPCs.csv')
    return (ges_spec_v3_NPCs,)


@app.cell
def _(ges_spec_v3_NPCs, npc_genes):
    ges_spec_v3_NPCs_enriched_genes = ges_spec_v3_NPCs[ges_spec_v3_NPCs["gene"].isin(npc_genes)][["gene", "ges_score"]]
    ges_spec_v3_NPCs_enriched_genes.sort_values("ges_score", ascending=False)
    return


@app.cell
def _(full_enrichment_results):
    neuroblast_genes_str = full_enrichment_results.loc[
        full_enrichment_results["condition"] == "Neuroblast",
        "lead_genes"
    ].iloc[0]
    neuroblast_genes = set(neuroblast_genes_str.split(';'))
    return (neuroblast_genes,)


@app.cell
def _(pd):
    ges_spec_v3_neuroblast = pd.read_csv('/miridan-data/annaludmir/data/ges_results/ges_spec_v3_Neuroblast.csv')
    return (ges_spec_v3_neuroblast,)


@app.cell
def _(ges_spec_v3_neuroblast, neuroblast_genes):
    ges_spec_v3_Neuroblast_enriched_genes = ges_spec_v3_neuroblast[ges_spec_v3_neuroblast["gene"].isin(neuroblast_genes)][["gene", "ges_score"]]
    ges_spec_v3_Neuroblast_enriched_genes.sort_values("ges_score", ascending=False)
    return


@app.cell
def _(adata, neuroblast_genes, np, rc_context, sc, sp):
    X_1 = adata[:, list(neuroblast_genes)].X
    if sp.issparse(X_1):
        X_1 = X_1.toarray()
    expr_frac_1 = (X_1 > 0).mean(axis=0)
    expr_frac_1 = np.ravel(expr_frac_1)
    filtered_genes_2 = list(neuroblast_genes)
    for i_1 in range(0, len(filtered_genes_2), 16):
        subset_1 = filtered_genes_2[i_1:i_1 + 16]
        print(f'Plotting genes {i_1 + 1}–{i_1 + len(subset_1)} of {len(filtered_genes_2)}')
        with rc_context({'figure.figsize': (3, 3)}):
            sc.pl.umap(adata, color=subset_1, s=1, frameon=False, ncols=4, vmax='p99')
    return


@app.cell
def _(adata, g1_genes, np, rc_context, sc, sp):
    X_2 = adata[:, list(g1_genes)].X
    if sp.issparse(X_2):
        X_2 = X_2.toarray()
    expr_frac_2 = (X_2 > 0).mean(axis=0)
    expr_frac_2 = np.ravel(expr_frac_2)
    filtered_genes_3 = list(g1_genes)
    for i_2 in range(0, len(filtered_genes_3), 16):
        subset_2 = filtered_genes_3[i_2:i_2 + 16]
        print(f'Plotting genes {i_2 + 1}–{i_2 + len(subset_2)} of {len(filtered_genes_3)}')
        with rc_context({'figure.figsize': (3, 3)}):
            sc.pl.umap(adata, color=subset_2, s=1, frameon=False, ncols=4, vmax='p99')
    return


@app.cell
def _(adata, sc):
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    return


@app.cell
def _(adata, pd, rc_context, sc):
    import re

    def _strip_ensembl_version(x: str) -> str:
        return re.sub(r"\.\d+$", "", str(x))

    sym_col = "Gene"  # change if your symbols are in another column

    # symbol -> exact var_names (with version)
    sym2varname = pd.Series(adata.var_names.values, index=adata.var[sym_col].astype(str)).to_dict()

    filtered_symbols = ["FOXG1", "LHX2", "EOMES"]

    filtered_varnames = [sym2varname[s] for s in filtered_symbols if s in sym2varname]
    missing = sorted(set(filtered_symbols) - set(sym2varname))
    print("Missing symbols:", missing)
    print("Using var_names:", filtered_varnames)

    for sym in filtered_symbols:
        ens = sym2varname.get(sym)
        if ens is None:
            print(f"⚠️ {sym} not found in adata.var[{sym_col}]")
            continue

        with rc_context({"figure.figsize": (3, 3)}):
            sc.pl.umap(
                adata,
                color=ens,          # ✅ must be something in var_names
                s=1,
                frameon=False,
                vmax="p99",
                title=sym           # ✅ show symbol as title
            )
    return


@app.cell
def _(expr, np, pd, sc, sp):
    #upload data
    print('uploading data')
    adata=sc.read_h5ad("/miridan-data/annaludmir/ndd_gene_modules/data/human_dev.h5ad")

    #upload data
    import seaborn as sns
    from matplotlib import pyplot as plt

    print('filtering data')
    adata_prolif = adata[adata.obs["cell_cycle_score"] > 0.004].copy()
    print('normalizing data')
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    genes = ["FOXG1", "LHX2", "EOMES"]
    region_col = "Region"

    sym_col = "Gene"  # change if your symbols are in another column

    # symbol -> exact var_names (with version)
    sym2varname = pd.Series(adata.var_names.values, index=adata.var[sym_col].astype(str)).to_dict()

    filtered_varnames = [sym2varname[s] for s in genes if s in sym2varname]
    missing = sorted(set(genes) - set(sym2varname))
    print("Missing symbols:", missing)
    print("Using var_names:", filtered_varnames)

    X = adata_prolif[:, filtered_varnames].X
    if sp.issparse(X):
        X = X.toarray()
    expr_frac = (X > 0).mean(axis=0)
    expr_frac = np.ravel(expr_frac)

    for sym in genes:
        print(sym)
        ens = sym2varname.get(sym)
        df = pd.DataFrame({
            region_col: adata_prolif.obs[region_col].values,
            sym: expr
        })

        plt.figure(figsize=(10, 4))
        sns.boxplot(data=df, x=region_col, y=sym)
        plt.xticks(rotation=45, ha="right")
        plt.title(f"{sym} expression in proliferating cells")
        plt.tight_layout()
        plt.show()
    return adata, expr_frac


@app.cell
def _(adata):
    adata.obs
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
