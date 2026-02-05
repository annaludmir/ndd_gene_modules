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
    import modules.specificity_score_calculations as ssc
    reload(ssc)
    return (ssc,)


@app.cell
def _(ssc):
    ssc.run_ges_pipeline(config_path='/miridan-data/annaludmir/ndd_gene_modules/config_files/ges_score_cortex_config.yaml')
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
    tsne.plot_tsne_for_gene_list("/miridan-data/annaludmir/ndd_gene_modules/results/ges_score_results/ges_score_for_all_layers_20251225/metadata/edited_adata.h5ad", "/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/ID_SYNDD_Neoplasm_all_layers_threshold_1_20260105/data/enrichment_results/GSEA/GSEA_final_summary.csv", "/miridan-data/annaludmir/ndd_gene_modules/results/ges_score_results/ges_score_for_all_layers_20251225/","/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/ID_SYNDD_Neoplasm_all_layers_threshold_1_20260105/data/additional_figures/","TSNE_for_GSEA")
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
    dot.plot_dots_for_gene_list("/miridan-data/annaludmir/ndd_gene_modules/results/ges_score_results/ges_score_for_cortex_20251221/metadata/edited_adata.h5ad", "/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/autism/autism_genes_cortex_threshold_1_20260105/data/enrichment_results/GSEA/GSEA_final_summary.csv","/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/autism/autism_genes_cortex_threshold_1_20260105/data/additional_figures/","Dots_for_GSEA")
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
    ecll.main('/miridan-data/annaludmir/ndd_gene_modules/config_files/enrichment_cortex_config_tau.yaml','/miridan-data/annaludmir/ndd_gene_modules/data/genes/')
    return


@app.cell
def _():
    import pandas as pd
    return (pd,)


@app.cell
def _(pd):
    summary = pd.read_csv('ndd_gene_modules/results/enrichment_results/batch_summary_20260110.csv')
    return (summary,)


@app.cell
def _(summary):
    summary.iloc()
    return


@app.cell
def _(summary):
    queried_summary = summary.query("is_significant == False & NES > 1.5 & FDR_qval_BH < 0.2")
    queried_summary
    return


@app.cell
def _(summary):
    queried_summary_1 = summary.query("run_name == 'updated_hsg_list_all_layers' & NES > 1.5")
    queried_summary_1
    return


@app.function
def check_similiarty(list_a, list_b):
    union_list = list(set(list_a) | set(list_b))
    intersection_list = set(list_a) & set(list_b)
    if len(union_list) in [len(list_a), len(list_b)]:
        return 1.0
    if len(intersection_list)/min(len(list_a), len(list_b)) > 0.5:
        return len(intersection_list)/min(len(list_a), len(list_b))
    return len(intersection_list)/len(union_list)


@app.cell
def _(pd):
    microcephaly_cortex = pd.read_csv('/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/microcephaly/microcephaly_genes_all_layers_threshold_1_20260110/data/enrichment_results/GSEA/GSEA_final_summary.csv')
    return (microcephaly_cortex,)


@app.cell
def _(pd):
    neoplasm_cortex = pd.read_csv('/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/ID/ID_SYNDD_Neoplasm_all_layers_threshold_1_20260110/data/enrichment_results/GSEA/GSEA_final_summary.csv')
    return (neoplasm_cortex,)


@app.cell
def _(microcephaly_cortex):
    microcephaly_cortex
    return


@app.cell
def _(neoplasm_cortex):
    neoplasm_cortex
    return


@app.cell
def _(pd):
    import re

    def _parse_lead_genes(x):
        """Turn a Lead_genes cell into a list of genes."""
        if pd.isna(x):
            return []
        s = str(x).strip()
        if not s:
            return []
        # robust split: handles ';' or ',' or tabs
        parts = re.split(r"[;,]\s*|\t+", s)
        return [p.strip() for p in parts if p.strip()]

    def compare_leading_genes_similarity(
        micro_df: pd.DataFrame,
        neo_df: pd.DataFrame,
        lead_col: str = "Lead_genes",
    ):
        # Keep only needed columns
        need = ["column", "condition", lead_col, "NES", "FDR q-val (BH corrected)"]
        for df_name, df in [("micro_df", micro_df), ("neo_df", neo_df)]:
            missing = [c for c in need if c not in df.columns]
            if missing:
                raise ValueError(f"{df_name} is missing columns: {missing}")

        m = micro_df[need].copy()
        n = neo_df[need].copy()

        # Parse leading genes
        m["lead_list_micro"] = m[lead_col].map(_parse_lead_genes)
        n["lead_list_neo"] = n[lead_col].map(_parse_lead_genes)
        m["NES_micro"] = m["NES"]
        n["NES_neo"] = n["NES"]
        m["FDR_micro"] = m["FDR q-val (BH corrected)"]
        n["FDR_neo"] = n["FDR q-val (BH corrected)"]
    

        # Merge on column+condition
        merged = m.merge(
            n[["column", "condition", "lead_list_neo","NES_neo","FDR_neo"]],
            on=["column", "condition"],
            how="inner",
        )

        if merged.empty:
            print("No overlapping (column, condition) pairs found between the two runs.")
            return

        # Compute + print
        for _, row in merged.iterrows():
            col = row["column"]
            cond = row["condition"]
            a = row["lead_list_micro"]
            b = row["lead_list_neo"]
            a_nes = round(row["NES_micro"], 2)
            b_nes = round(row["NES_neo"], 2)
            a_fdr = round(row["FDR_micro"], 3)
            b_fdr = round(row["FDR_neo"], 3)
        
            score = check_similiarty(a, b)  # <-- your function
            print(f"{col} - {cond}: similarity = {score:.3f}  (micro n={len(a)};NES={a_nes};FDR={a_fdr}, neo n={len(b)};NES:{b_nes};FDR:{b_fdr})")

    return compare_leading_genes_similarity, re


@app.cell
def _(compare_leading_genes_similarity, microcephaly_cortex, neoplasm_cortex):
    compare_leading_genes_similarity(microcephaly_cortex,neoplasm_cortex)
    return


@app.cell
def _(pd, re):
    import numpy as np
    from pathlib import Path

    def _parse_gene_list(x):
        """Parse Lead_genes string into a list (supports ; , or tab separators)."""
        if pd.isna(x):
            return []
        s = str(x).strip()
        if not s:
            return []
        parts = re.split(r"[;,]\s*|\t+", s)
        return [p.strip() for p in parts if p.strip()]

    def condition_x_condition_matrices_within_df(
        df: pd.DataFrame,
        out_dir: str | None = None,
        print_matrices: bool = True,
        include_self: bool = True,
    ) -> dict[str, pd.DataFrame]:
        """
        For each `column` in df, build a condition×condition similarity matrix
        based on df["Lead_genes"] and your `check_similiarty(listA, listB)`.

        Returns: {column_name: matrix_df}
        """
        required = {"column", "condition", "Lead_genes"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"df missing required columns: {sorted(missing)}")

        d = df[["column", "condition", "Lead_genes"]].copy()
        d["lead"] = d["Lead_genes"].map(_parse_gene_list)

        out_path = Path(out_dir) if out_dir else None
        if out_path:
            out_path.mkdir(parents=True, exist_ok=True)

        matrices = {}

        for col in sorted(d["column"].unique()):
            sub = d[d["column"] == col].copy()

            # If a condition appears multiple times, union its lead genes
            cond_to_genes: dict[str, set[str]] = {}
            for _, r in sub.iterrows():
                c = str(r["condition"])
                cond_to_genes.setdefault(c, set()).update(r["lead"])

            conds = list(cond_to_genes.keys())  # keep current order in df
            mat = pd.DataFrame(index=conds, columns=conds, dtype=float)

            for a in conds:
                for b in conds:
                    if (not include_self) and (a == b):
                        mat.loc[a, b] = np.nan
                    else:
                        mat.loc[a, b] = check_similiarty(
                            list(cond_to_genes[a]),
                            list(cond_to_genes[b]),
                        )

            matrices[col] = mat

            if print_matrices:
                print(f"\n=== {col}: condition × condition similarity ===")
                print(mat.round(3))

            if out_path:
                safe_col = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(col))
                mat.to_csv(out_path / f"{safe_col}_condition_x_condition_similarity.csv")

        return matrices

    return (condition_x_condition_matrices_within_df,)


@app.cell
def _(condition_x_condition_matrices_within_df, microcephaly_cortex):
    condition_x_condition_matrices_within_df(microcephaly_cortex)
    return


@app.cell
def _():
    micro = ['BUB1','KIF14','CENPE','ASPM','CENPF','CEP55','CKAP2L','NUF2','NDE1','KNL1','BUB1B','CIT','KIF11','NCAPH','NCAPD2','STIL','WDR62','POC1A','PLK4','LMNB2','PPP1R35','RAD21','TRAIP','CEP135','LMNB1','FILIP1','CDK5RAP2','XRCC4','FANCM','TRIP13','CEP152','SLX4','FANCD2','SASS6','BRCA2','TRMT10A','HMGB1','FANCB','DNA2','NSD2','NCAPD3','NIPBL','HIST1H4C','LARP7','PCNT','TNPO2','ORC6']
    return (micro,)


@app.cell
def _():
    hsg = ['NEK2','FAM72C','FAM72D','ASPM','KIF18A','FAM72B','ARHGAP11A','SPAG5','ZNF492','NOTCH2','CDK5RAP2','AR','NBPF14','GKAP1','PDCD4','CDH12','ANKRD20A4','ZFP36L1','SRGAP2B','SRGAP2C','TLN1']
    return (hsg,)


@app.cell
def _(hsg, micro):
    check_similiarty(micro,hsg)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
