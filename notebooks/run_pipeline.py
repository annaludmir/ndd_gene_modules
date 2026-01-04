import marimo

__generated_with = "0.17.8"
app = marimo.App()

with app.setup:
    # Initialization code that runs before all other cells
    import marimo as mo
    import sys
    from importlib import reload
    from types import ModuleType
    sys.path.append('/miridan-data/annaludmir/ndd_gene_modules')



@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ##Calculating GES score for all genes in data
    """)
    return


@app.cell
def _():
    import modules.ges_score_calculations as gsc
    reload(gsc)
    return (gsc,)


@app.cell
def _(gsc):
    gsc.run_ges_pipeline(config_path='/miridan-data/annaludmir/ndd_gene_modules/config_files/ges_score_cortex_config.yaml')
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Calculcating enrichment using GSEA for selected genes
    """)
    return


@app.cell
def _():
    import modules.enrichment_pipeline_for_gene_list as epfgl 
    reload(epfgl)

    epfgl.run_gene_list_pipeline(config_path='/miridan-data/annaludmir/ndd_gene_modules/config_files/enrichment_cortex_cell_phase_config.yaml')
    return (epfgl,)


@app.cell
def _(epfgl):
    epfgl.run_gene_list_pipeline(config_path='/miridan-data/annaludmir/ndd_gene_modules/config_files/enrichment_cortex_config.yaml')
    return


@app.cell
def _(epfgl):
    epfgl.run_gene_list_pipeline(config_path='/miridan-data/annaludmir/ndd_gene_modules/config_files/enrichment_all_layers_config.yaml')
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Plotting tSNE plots
    """)
    return


@app.cell
def _():
    import modules.tsne_plots as tsne
    reload(tsne)
    tsne.plot_tsne_for_gene_list("/miridan-data/annaludmir/ndd_gene_modules/data/Cortex_EMX1_louvain3_passedQC_PostM_rev1.h5ad", "/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/ID_AFS_HM_cortex_threshold_1_20251211/data/enrichment_results/GSEA/GSEA_final_summary.csv", "/miridan-data/annaludmir/ndd_gene_modules/results/ges_score_results/ges_score_for_cortex_20251211/","/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/ID_AFS_HM_cortex_threshold_1_20251211/data/additional_figures/","TSNE_for_GSEA")
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ##Plotting dot plots
    """)
    return


@app.cell
def _():
    import modules.dot_plots as dot
    reload(dot)
    dot.plot_dots_for_gene_list("/miridan-data/annaludmir/ndd_gene_modules/data/Cortex_EMX1_louvain3_passedQC_PostM_rev1.h5ad", "/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/ID_AFS_HM_cortex_threshold_1_20251211/data/enrichment_results/GSEA/GSEA_final_summary.csv","/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/ID_AFS_HM_cortex_threshold_1_20251211/data/additional_figures/","Dots_for_GSEA")
    return


@app.cell
def _():
    mo.md(r"""
    ##Playing around with the data
    """)
    return


@app.cell
def _():
    import scanpy as sc
    adata = sc.read_h5ad("/miridan-data/annaludmir/ndd_gene_modules/data/Cortex_EMX1_louvain3_passedQC_PostM_rev1.h5ad")
    adata.obs["CellCycleFraction"][0]
    return (adata,)


@app.cell
def _(adata):
    adata.obs["CellCycleFraction"] <= 0.004
    return


@app.cell
def _():
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ##Running in Loop for many gene lists
    """)
    return


@app.cell
def _():
    import modules.enrichment_cal_lists_loop as ecll
    reload(ecll)
    ecll.main('/miridan-data/annaludmir/ndd_gene_modules/config_files/enrichment_cortex_cell_phase_config.yaml','/miridan-data/annaludmir/ndd_gene_modules/data/genes/')
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
