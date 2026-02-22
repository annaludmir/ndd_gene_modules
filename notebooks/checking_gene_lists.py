import marimo

__generated_with = "0.17.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return


@app.cell
def _():
    import pandas as pd
    return (pd,)


@app.cell
def _(df1, df2):
    def compare_gene_lists(df1_path, df2vcvb, gene_col="gene"):
        """
        Compare two gene lists.

        Parameters
        ----------
        df1, df2 : pd.DataFrame
            DataFrames containing gene symbols
        gene_col : str
            Column name with gene symbols

        Returns
        -------
        dict
            Dictionary with overlap, unique genes, and union size
        """
        genes1 = set(df1[gene_col].dropna())
        genes2 = set(df2[gene_col].dropna())

        overlap = genes1 & genes2
        unique_1 = genes1 - genes2
        unique_2 = genes2 - genes1
        union = genes1 | genes2

        return {
            "overlap_n": len(overlap),
            "overlap_genes": sorted(overlap),

            "unique_df1_n": len(unique_1),
            "unique_df1_genes": sorted(unique_1),

            "unique_df2_n": len(unique_2),
            "unique_df2_genes": sorted(unique_2),

            "union_n": len(union)
        }
    return (compare_gene_lists,)


@app.cell
def _(compare_gene_lists, pd):
    microcephaly_severe = pd.read_csv('/miridan-data/annaludmir/ndd_gene_modules/data/genes/microcephaly_genes.csv')
    microcephaly = pd.read_csv('/miridan-data/annaludmir/ndd_gene_modules/data/genes/ID_SYNDD_Microcephaly.csv')

    unique_severe_microcephaly_genes = compare_gene_lists(microcephaly_severe,microcephaly)["unique_df1_genes"]
    return (unique_severe_microcephaly_genes,)


@app.cell
def _(unique_severe_microcephaly_genes):
    unique_severe_microcephaly_genes
    return


@app.cell
def _(pd):
    def collect_leading_genes(
        summary_csv,
        run_prefix="microcephaly_genes_",
        run_col="run_name",
        lead_genes_col="lead_genes",
        sep=";"
    ):
        """
        Collect all unique leading genes for runs matching a prefix.

        Returns
        -------
        list[str]
            Sorted list of unique gene symbols
        """

        df = pd.read_csv(summary_csv)

        # filter runs
        df_filt = df[df[run_col].str.startswith(run_prefix, na=False)]
        df_filt = df_filt[df_filt["is_significant"] == True]

        # split + flatten + deduplicate
        genes = (
            df_filt[lead_genes_col]
            .dropna()
            .str.split(sep)
            .explode()
            .str.strip()
            .unique()
        )

        return sorted(genes)

    return (collect_leading_genes,)


@app.cell
def _(collect_leading_genes):
    genes = collect_leading_genes(
        "/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/batch_summary_20260110.csv",
        run_prefix="microcephaly_genes_"
    )

    print(len(genes))
    print(genes)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
