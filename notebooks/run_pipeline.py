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
    ssc.run_ges_pipeline(config_path='/miridan-data/annaludmir/ndd_gene_modules/config_files/ges_score_all_layers_config_proliferating.yaml')
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


@app.cell
def _(epfgl):
    epfgl.run_gene_list_pipeline(config_path='/miridan-data/annaludmir/ndd_gene_modules/config_files/enrichment_all_layers_proliferating_config.yaml')
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
    ecll.main('/miridan-data/annaludmir/ndd_gene_modules/config_files/enrichment_all_layers_config.yaml','/miridan-data/annaludmir/ndd_gene_modules/data/genes/schizophrenia_bipolar/')
    return


@app.cell
def _():
    import pandas as pd
    return (pd,)


@app.cell
def _(pd):
    summary = pd.read_csv('/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/batch_summary_20260212_all_layers.csv')
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
    queried_summary_1 = summary.query("column_condition_value == 'Fibroblast' & is_significant == True")
    queried_summary_1
    return


@app.cell
def _(pd):
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns


    def plot_gsea_summary_relationships(summary_csv):
        df = pd.read_csv(summary_csv)

        # ensure numeric
        cols = [
            "FDR_qval_BH",
            "NES",
            "num_genes_in_gene_list",
            "num_of_lead_genes",
        ]
        for c in cols:
            df[c] = pd.to_numeric(df[c], errors="coerce")

        df = df.dropna(subset=cols)

        # nicer plotting style
        sns.set(style="whitegrid", context="talk")

        fig, axes = plt.subplots(2, 3, figsize=(18, 10))

        # 1️⃣ FDR vs gene list size
        sns.scatterplot(
            data=df,
            x="num_genes_in_gene_list",
            y="FDR_qval_BH",
            ax=axes[0, 0],
            alpha=0.6
        )
        axes[0, 0].set_yscale("log")
        axes[0, 0].set_title("FDR vs Gene List Size")
        axes[0, 0].set_xlabel("Number of genes in gene list")
        axes[0, 0].set_ylabel("FDR (log scale)")

        # 2️⃣ NES vs gene list size
        sns.scatterplot(
            data=df,
            x="num_genes_in_gene_list",
            y="NES",
            ax=axes[0, 1],
            alpha=0.6
        )
        axes[0, 1].axhline(0, ls="--", c="black", lw=1)
        axes[0, 1].set_title("NES vs Gene List Size")
        axes[0, 1].set_xlabel("Number of genes in gene list")
        axes[0, 1].set_ylabel("NES")

        # 3️⃣ FDR vs leading genes
        sns.scatterplot(
            data=df,
            x="num_of_lead_genes",
            y="FDR_qval_BH",
            ax=axes[0, 2],
            alpha=0.6
        )
        axes[0, 2].set_yscale("log")
        axes[0, 2].set_title("FDR vs Leading Genes")
        axes[0, 2].set_xlabel("Number of leading genes")
        axes[0, 2].set_ylabel("FDR (log scale)")

        # 4️⃣ NES vs leading genes
        sns.scatterplot(
            data=df,
            x="num_of_lead_genes",
            y="NES",
            ax=axes[1, 0],
            alpha=0.6
        )
        axes[1, 0].axhline(0, ls="--", c="black", lw=1)
        axes[1, 0].set_title("NES vs Leading Genes")
        axes[1, 0].set_xlabel("Number of leading genes")
        axes[1, 0].set_ylabel("NES")

        # 5️⃣ NEW: leading genes vs gene list size
        sns.scatterplot(
            data=df,
            x="num_genes_in_gene_list",
            y="num_of_lead_genes",
            ax=axes[1, 1],
            alpha=0.6
        )
        axes[1, 1].set_title("Leading Genes vs Gene List Size")
        axes[1, 1].set_xlabel("Number of genes in gene list")
        axes[1, 1].set_ylabel("Number of leading genes")

        # 6️⃣ empty panel → turn off
        axes[1, 2].axis("off")

        plt.tight_layout()
        plt.show()
    return (plot_gsea_summary_relationships,)


@app.cell
def _(plot_gsea_summary_relationships):
    plot_gsea_summary_relationships('/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/batch_summary_20260215_cell_phase.csv')
    return


@app.cell
def _():
    mo.md(r"""
    ##Check similarity between leading genes
    """)
    return


@app.cell
def _():
    from matplotlib import pyplot as plt
    from matplotlib_venn import venn2, venn3
    from upsetplot import UpSet, from_contents


    def plot_gene_overlap(
        gene_sets: dict,
        cluster,
        similiarity = 0):
        """
        Plot overlap for multiple gene lists.

        Parameters
        ----------
        gene_sets : dict
            {"label1": iterable_of_genes, "label2": iterable_of_genes, ...}
        title : str
        """

        # convert everything to sets
        gene_sets = {k: set(v) for k, v in gene_sets.items()}
        n = len(gene_sets)

        if n == 2:
            labels = list(gene_sets.keys())
            a, b = labels

            plt.figure(figsize=(5, 5))
            venn2([gene_sets[a], gene_sets[b]], set_labels=labels)
            plt.title(f"{cluster} Leading Genes - Similarity {similiarity}")
            plt.show()

        elif n == 3:
            labels = list(gene_sets.keys())
            a, b, c = labels

            plt.figure(figsize=(6, 6))
            venn3([gene_sets[a], gene_sets[b], gene_sets[c]], set_labels=labels)
            plt.title(f"{cluster} Leading Genes")
            plt.show()

        else:
            # UpSet plot (best for 4+ sets)
            contents = {k: list(v) for k, v in gene_sets.items()}
            upset_data = from_contents(contents)

            plt.figure(figsize=(10, 6))
            UpSet(upset_data, show_counts=True).plot()
            plt.suptitle(f"{cluster} Leading Genes")
            plt.show()
    return (plot_gene_overlap,)


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
def _(pd, plot_gene_overlap):
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
        label_a: str,
        label_b: str,
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

            if score > 0.5:
                plot_gene_overlap({label_a: a, label_b: b}, cond, round(score, 3))
    return compare_leading_genes_similarity, re


@app.cell
def _(compare_leading_genes_similarity, microcephaly_cortex, neoplasm_cortex):
    compare_leading_genes_similarity(microcephaly_cortex,neoplasm_cortex, "Microchepahly", "IDD Neoplasm")
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
    microcephaly_G1 = ["GINS2","CDC6","BRIP1","ORC1","GINS3","HPDL","CDT1","ATP1A2","MCM7","SLC38A3","BRCA2","GMNN","FANCA","CDK6","DDX11","ORC6","FANCB","BLM","RBBP8","PUS7","FANCI","FILIP1","DNA2","FANCE","MFSD2A","PALB2","SMO","FANCD2","PRIM1","MRE11","FANCG"]
    return (microcephaly_G1,)


@app.cell
def _():
    microcephaly_S = ["HIST1H4C","WDR62","BLM","FANCD2","FANCI","POC1A","CEP152","CDC6","STIL","CDT1","BRIP1","BRCA2","NCAPH","PLK4","FANCB","ORC6","KIF11","ORC1","DNA2","CKAP2L","GMNN","BUB1B","TRIP13","KNL1","NCAPD3","GINS3","TRAIP","FANCG","CIT","RBBP8","FANCA","NUF2","CEP55","NCAPD2","LMNB2","DDX11","CEP135","KIF14","CENPE","MCM7","ASPM","FILIP1","GINS2","BUB1","SASS6","PRIM1","CENPF","TUBG1","HMGB1","FANCM","TUBGCP3","NSD2","EOMES","CDK5RAP2","LMNB1","RMI1","RAD51C","VRK1","SMC1A","CDK6","NUP107","GPT2","MRE11","XRCC4","RTTN","SMO","FANCC","TRA2B","PSMC3","NDE1","PCNT","NUP188","LHX2","TTI1","SMC5","SMC3","FANCE","SLF2","RAD21","PPP1R35","CPSF3","PALB2","RRP7A","NSRP1","LARP7"]
    return (microcephaly_S,)


@app.cell
def _():
    microcephaly_G2M = ["BUB1","KIF14","CENPE","ASPM","CENPF","CEP55","CKAP2L","NUF2","NDE1","KNL1","BUB1B","CIT","KIF11","NCAPH","NCAPD2","STIL","WDR62","POC1A","PLK4","LMNB2","PPP1R35","RAD21","TRAIP","CEP135","LMNB1","FILIP1","CDK5RAP2","XRCC4","FANCM","TRIP13","CEP152","SLX4","FANCD2","SASS6","BRCA2","TRMT10A","HMGB1","FANCB","DNA2","NSD2","NCAPD3","NIPBL","HIST1H4C","LARP7","PCNT","TNPO2","ORC6"]
    return (microcephaly_G2M,)


@app.cell
def _():
    microcephaly_forebrain = ["FOXG1","CDK6","HPDL","BRIP1","CIT","KIF11","CEP152","KNL1","BLM","ORC1","FANCI","CDT1","BUB1B","POC1A","NCAPH","CDC6","CKAP2L","DNA2","NCAPD2","WDR62","BUB1","STIL","ZEB2","PLK4","CENPF","BRCA2","ASPM","CCDC88A","FANCB","FANCA","CENPE","LMNB2","HHAT","FANCD2","CEP55","NDE1","GINS2","MRE11","DDX11","RBBP8","TRAIP","FANCL","KIF14","GINS3","CEP135","FANCE","NUF2","TRIP13","GMNN","CDK5RAP2","TUBGCP3","CCND2","ORC6","VRK1","MCM7","LMNB1","DIAPH1","PCNT","SMC1A","NUP107","NUP188","SMO","ATP1A2","FANCG","GPT2","SASS6","PUS7","HIST1H4C","TCF4","TOP3A","NCAPD3","NUP214","NSD2","RMI1","SMC5","RAD21","FBRSL1","OSGEP","EFTUD2","TUBGCP4","ATRIP","TTI1","ZNF526","SLF2","SMC3","ANKLE2"]
    return (microcephaly_forebrain,)


@app.cell
def _(microcephaly_G1, microcephaly_forebrain):
    check_similiarty(microcephaly_G1,microcephaly_forebrain)
    return


@app.cell
def _(
    microcephaly_G1,
    microcephaly_G2M,
    microcephaly_S,
    microcephaly_forebrain,
    plot_gene_overlap,
):
    plot_gene_overlap(
        {'G1': microcephaly_G1, 'S': microcephaly_S, 'G2M': microcephaly_G2M, 'Forebrain': microcephaly_forebrain}, "Microcehpaly"
    )
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
