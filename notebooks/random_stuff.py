import marimo

__generated_with = "0.17.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    mo.md(r"""
    #Check leading genes location in prerank list
    """)
    return


@app.cell
def _():
    gene_list = 'PDCD6IP','WDR4','DPP6','ZNHIT3','POGZ','CREBBP','TP53RK','CCND2','COPB1','STAMBP','TRAPPC12','WDR11','AP4E1','ARCN1','NIPBL','COPB2','KMT2B','ZNF335','ATRX','FRA10AC1','BRD4','CSNK2A1','RPL10','MCPH1','CHKA','TUBGCP4','PCNT','CRIPT','MECP2','EXOC7','TTC5','HDAC8','PTPN23','COG3','ORC4','TRIO','MINPP1','UBE3A','YIF1B','WDR37','NIN','PCLO','VPS50','TAF13','TRAPPC9','AP4M1','ARF3','MSMO1','PQBP1','RUSC2','AP4B1','UNC80','SMARCA5','PLAA','AKT3','ATP6V0A1','TRAPPC10','PDHA1','DYNC1I2','ZPR1','IGF1R','ATP9A','CASK','TRAPPC6B','CAMSAP1','ZNF668','SVBP'
    return (gene_list,)


@app.cell
def _():
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt


    def locate_genes_in_preranked(
        query_genes,
        preranked,
        gene_col=None,
        make_plot=True,
        use_percentile=False
    ):
        """
        Parameters
        ----------
        query_genes : list
            Genes to locate.

        preranked : list OR pandas Series OR pandas DataFrame
            Ordered gene list.

        gene_col : str, optional
            If preranked is a DataFrame, column containing gene names.

        make_plot : bool
            Whether to draw the position plot.

        use_percentile : bool
            Convert rank → percentile.

        Returns
        -------
        DataFrame with gene and rank.
        """

        # ---- extract ordered genes ----
        if isinstance(preranked, pd.DataFrame):
            if gene_col is None:
                raise ValueError("Provide gene_col for DataFrame input")
            ordered = preranked[gene_col].tolist()
        else:
            ordered = list(preranked)

        N = len(ordered)

        # ---- build rank lookup ----
        rank_dict = {g: i for i, g in enumerate(ordered)}

        ranks = []
        for g in query_genes:
            r = rank_dict.get(g, np.nan)
            if use_percentile and not np.isnan(r):
                r = r / N
            ranks.append(r)

        res = pd.DataFrame({
            "gene": query_genes,
            "rank": ranks
        })

        # ---- plot ----
        if make_plot:
            found = res.dropna()

            plt.figure(figsize=(12, 2))

            plt.scatter(
                found["rank"],
                np.zeros(len(found)),
                s=80
            )

            for _, row in found.iterrows():
                plt.text(
                    row["rank"],
                    0,
                    row["gene"],
                    rotation=90,
                    va="bottom",
                    ha="center",
                    fontsize=9
                )

            plt.yticks([])
            plt.xlabel("Rank in preranked list" if not use_percentile else "Percentile")
            plt.title("Gene locations in preranked list")

            plt.xlim(0, N if not use_percentile else 1)

            plt.tight_layout()
            plt.show()

        return res

    return locate_genes_in_preranked, pd


@app.cell
def _(pd):
    preranked = pd.read_csv('/miridan-data/annaludmir/ndd_gene_modules/results/ges_score_results/ges_score_for_cortex_20251221/data/ges_spec_CellClass_Glioblast.csv')
    return (preranked,)


@app.cell
def _(gene_list, locate_genes_in_preranked, preranked):
    locate_genes_in_preranked(gene_list,preranked,gene_col='gene')
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
