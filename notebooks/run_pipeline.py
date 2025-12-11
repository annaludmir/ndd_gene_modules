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
    gsc.run_ges_pipeline(config_path='/miridan-data/annaludmir/ndd_gene_modules/config_files/ges_score_cortex_cell_phase_config.yaml')
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

    epfgl.run_gene_list_pipeline(config_path='/miridan-data/annaludmir/ndd_gene_modules/config_files/enrichment_all_layers_config.yaml')
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Plotting tSNE plots
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ##Plotting dot plots
    """)
    return


if __name__ == "__main__":
    app.run()
