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
    from scipy.stats import hypergeom

    def overlap_pval(setA, setB, universe_size):
        A = set(setA)
        B = set(setB)

        N = universe_size
        K = len(A)
        M = len(B)
        k = len(A & B)

        # survival function = P(X >= k)
        pval = hypergeom.sf(k - 1, N, K, M)
        return pval, k
    return (overlap_pval,)


@app.cell
def _(overlap_pval, similarity):
    from matplotlib import pyplot as plt
    from matplotlib_venn import venn2, venn3
    from upsetplot import UpSet, from_contents


    def plot_gene_overlap(
        gene_sets: dict,
        cluster,
        similiarity = 0,
        universe_size = None):
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

            pval = None
            if universe_size is not None:
                pval, k = overlap_pval(gene_sets[a], gene_sets[b], universe_size)

            plt.figure(figsize=(5, 5))
            venn2([gene_sets[a], gene_sets[b]], set_labels=labels)

            title = f"{cluster} Leading Genes\nSimilarity={similarity}"
            if pval is not None:
                title += f"\nHypergeom p = {pval:.2e}"

            plt.title(title)
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
def _():
    microcephaly_NPCs = ['EOMES','KNL1','CEP55','BUB1B','WDR62','CKAP2L','KIF11','BUB1','CIT','STIL','NCAPH','ASPM','CEP152','NUF2','BRIP1','KIF14','DNA2','BLM','CENPE','POC1A','FANCD2','FANCI','PLK4','CENPF','NCAPD2','ORC1','BRCA2','CDT1','FANCB','FANCA','CDC6','ORC6','TRAIP','RBBP8','NDE1','NCAPD3','TRIP13','LMNB2','VRK1','CEP135','DDX11','GMNN','HIST1H4C','CDK6','CCND2','FANCG','FANCM','CDK5RAP2','LMNB1','MCM7','PRIM1','RTTN','GINS2','PCNT','TUBGCP3','SMC1A','RMI1','FOXG1','TOP3A','FANCL','GINS3','TUBG1','NSD2','FANCE','SASS6','MRE11','DONSON','XRCC4','NUP188']
    return (microcephaly_NPCs,)


@app.cell
def _():
    microcephaly_NPCs_cortex = ['EOMES','WDR62','CEP55','KNL1','CIT','KIF14','CENPE','BLM','CKAP2L','NUF2','KIF11','BUB1','NCAPH','FANCD2','BUB1B','ASPM','POC1A','FANCI','STIL','CDC6','CEP152','RBBP8','BRIP1','CENPF','DNA2','PLK4','BRCA2','NCAPD2','ORC6','ORC1','FANCB','CDT1','FANCG','NCAPD3','TRAIP','LMNB2','GINS3','GMNN','FANCA','DDX11','MFSD2A','CDK6','NDE1','MCM7','LMNB1','TUBG1','HIST1H4C','RMI1','PCNT','CDK5RAP2','NSD2','FANCM','CEP135','GINS2','NUP188','RTTN','TUBGCP3']
    return


@app.cell
def _():
    microcephaly_radial_glia = ['FILIP1','HPDL','CENPF','CEP55','NUF2','POC1A','CDC6','BUB1','BRCA2','CKAP2L','CENPE','ASPM','GINS2','NCAPH','TRIP13','CDT1','BUB1B','ORC1','FANCD2','KIF14','KIF11','ORC6','TRAIP','STIL','KNL1','GINS3','BLM','PLK4','FANCB','GMNN','MCM7','WDR62','BRIP1','FANCI','WLS','NCAPD2','CEP152','NDE1','MYCN','SLC25A19','CIT','CEP135','MFSD2A','TSEN15','DDX11','DNA2','RBBP8','SASS6','PPIL1','CDK6','FANCA','FANCG','TRMT10A','SMO']
    return (microcephaly_radial_glia,)


@app.cell
def _():
    microcephaly_radial_glia_cortex = ['FILIP1','BUB1B','BUB1','HPDL','ASPM','CKAP2L','CENPE','KIF14','CEP55','NUF2','POC1A','CENPF','NCAPH','CDT1','NDE1','KIF11','CIT','ORC1','GINS3','GINS2','TRAIP','KNL1','NCAPD2','CDC6','SASS6','ORC6','STIL','BLM','BRCA2','FANCB','GMNN','CEP152','MCM7','CEP135','TRMT10A','SMO','LMNB2','FANCD2','XRCC4','DDX11','PLK4','CDK6','DNA2','PNKP','PPFIBP1','WDR62','PRIM1','MFSD2A','CTNNB1','GPT2','METTL5']
    return


@app.cell
def _():
    microcephaly_glioblasts = ['ZEB2','HMGB1','DOHH','PDCD6IP','DYRK1A','NSD2','PTPN23','CKAP2L','EIF2S3','KIF14','UFM1','COPB1','UBA5','NUF2','PPP1R35','AKT3','ASPM','TPR','PSMC3','MCM7','PRUNE1','VPS50','BPTF','ANKLE2','RAD21','CTNNB1','RAD51C','INTS11','YIF1B','LMNB2','WDR73','CRIPT','WDFY3','CREBBP','TRAIP','RING1','TRMT1','POGZ','NAA20','HDAC8','CTU2','TUBGCP4','RTTN','NCAPD2','AP4E1','KNL1','SMC1A','QARS','OSGEP','EXOC7','ORC1','SASS6','NCAPH','TUBGCP2','NUP188','TUBGCP6','PUF60','DROSHA','PPP1R15B','BUB1B','TAF13','TRAPPC12','ERCC4','METTL5','LARP7','POC1A','STAMBP','TNPO2','TSEN54','TRIP13','DONSON','XRCC4','WDR62','TOP3A','NCAPD3','SLC9A6','HHAT','FANCG','STIL','RBBP8','FANCE','CDC6','NUP107','CIT','AP4S1','TTI1','FILIP1','DNA2','LAGE3','EFTUD2','GMNN','BLM','CDK5RAP2','NUP214','FANCL','UFC1','ORC6','NSMCE2','TUBGCP3','UGP2','FANCC','KIF1BP','COASY','NSRP1','PNKP','FOXG1','NBN','WDR11','NDE1','MYCN','KIF11','SMC5','TCF4','MSMO1','MRE11','PLK4','FANCB','GINS2','DDX11','CHKA']
    return (microcephaly_glioblasts,)


@app.cell
def _():
    microcephaly_glioblasts_cortex = ['ATP1A2','SLC38A3','BRIP1','FILIP1','PUS7','BRCA2','CDK6','FANCA','FANCI','SMO','GINS2','PLK4','FANCB','CDT1','ORC1','GINS3','FANCC','GMNN','FANCD2','CDC6','CEP135','CTSF','ATP11A','DDX11','ORC6','LHX2','MED11','MRE11','GPT2','FANCG','SMC5','AP4S1','KIF1BP','CTNNB1','SLF2','TUBGCP3','RBBP8','CEP152','STIL','ATR','FANCE','BLM','TRAIP','KIF11','WDR62','TTI1','DNA2','NCAPH','PNKP','NCAPD3','HPDL','MCM7','CCND2','NUP107','STAMBP','NSRP1','LARP7','WDR11','MSMO1','RING1','POC1A','COASY','CHKA','NDE1','NBN','TSEN54','SMC1A','ZNF526','NSD2','FANCM','TUBGCP4','QARS','EFTUD2','TAF13','PPP1R15B','ARCN1','SASS6','ANKLE2','PPFIBP1','TRMT1','RMI1','TRA2B','NUP188','TUBG1','CIT','FANCL','UFM1','MYCN','WDR73','LAGE3','PPIL1','XRCC4','PPP1R35','NAA20','PHC1','UGP2','TRMT10A','NUF2','RAD21','CKAP2L','CRIPT','NCAPD2']
    return


@app.cell
def _():
    microcephaly_fibroblasts = ['INTS11','SVBP','CKAP2L','NSRP1','RTTN','FANCA','TRAPPC12','ERCC8','WDR73','PALB2','TP53RK','HIKESHI','TRIO','RAD51C','FANCD2','CEP135','PTPN23','KIF11','HDAC8','TRAPPC10','WDR4','WDR11','DNMT3A','CCND2','FANCC','TRIP13','ATRIP','TRAIP','MPLKIP','RAD50','PDHA1','STIL','KMT2B','TUBG1','TCF4','FRA10AC1','SMC1A','TUBGCP6','PQBP1','NUP188','CDT1','TUBGCP2','ACBD6','BLM','CTNNB1','NCAPH','BRD4','ADARB1','LARP7','RPL10','PPIL1','NUP214','ORC6','KNL1','SLF2','ANKLE2','MRE11','TAF13','EFTUD2','COASY','AP4S1','CENPF','GINS2','NDE1','COG3','CENPE','TSEN54','ORC1','PSMC3','CREBBP','UGP2','SLC25A19','TRA2B','TRMT10A','DPM1','GTF2E2','DOHH','POC1A','YIF1B','PPP1R15B','DONSON','RING1','FANCL','EIF5A','SMC5','BRIP1','CDK6','DDX11','UFM1','LMNB2','NSMCE2','TRMT1','CTU2','NAA20','RRP7A','HIST1H4C','GINS3','ZNF335','EIF2S3','TPR','CHKA','METTL5','COPB2','PNKP','COPB1','CDC6','AP4M1','UFC1','ARPC4']
    return


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
    regions = {
        'Forebrain': ['FOXG1','CDK6','HPDL','BRIP1','CIT','KIF11','CEP152','KNL1','BLM','ORC1','FANCI','CDT1','BUB1B','POC1A','NCAPH','CDC6','CKAP2L','DNA2','NCAPD2','WDR62','BUB1','STIL','ZEB2','PLK4','CENPF','BRCA2','ASPM','CCDC88A','FANCB','FANCA','CENPE','LMNB2','HHAT','FANCD2','CEP55','NDE1','GINS2','MRE11','DDX11','RBBP8','TRAIP','FANCL','KIF14','GINS3','CEP135','FANCE','NUF2','TRIP13','GMNN','CDK5RAP2','TUBGCP3','CCND2','ORC6','VRK1','MCM7','LMNB1','DIAPH1','PCNT','SMC1A','NUP107','NUP188','SMO','ATP1A2','FANCG','GPT2','SASS6','PUS7','HIST1H4C','TCF4','TOP3A','NCAPD3','NUP214','NSD2','RMI1','SMC5','RAD21','FBRSL1','OSGEP','EFTUD2','TUBGCP4','ATRIP','TTI1','ZNF526','SLF2','SMC3','ANKLE2'],
        'Telecenphalon':
    ['EOMES','FOXG1','LHX2','VRK1','CCND2','TCF4','PUS7','BUB1B','FANCL','AKT3','PCNT','ASPM','BUB1','FANCA','DNA2','CKAP2L','CEP135','TRAPPC9','MYCN','KIF11','BRIP1','KNL1','DIAPH1','CEP152','TP53RK','CDK5RAP2','FANCI','PLK4','AP4B1','NCAPH','NUF2','CEP63','ATP11A','NDE1','SLF2','TUBGCP3','TRIP13','TRIO','ZEB2','NUP188','CENPF','PCLO','ACBD6','NCAPD2','HIST1H4C','FANCD2','CDT1','PDHA1','MCPH1','TOP3A','LMNB1','POC1A','ORC6','POGZ','PRIM1','CEP55','DROSHA','NUP107','WDR62','NCAPD3','KIF14','GINS2','BRCA2','LMNB2','ATRIP','CAMSAP1','STIL','ADARB1','ATRX','PALB2','FANCE','ATR','SMC1A','EXOC7','FANCB','TUBGCP6','KMT2B','ZNF526','EFTUD2','MCM7','DDX11','VPS50','SMC5','IGF1R','ZPR1','GPT2','CCDC88A','CENPE','MRE11','SMC3','WDFY3','FANCM','NSD2','OSGEP','PTPN23','CREBBP','CIT','NIPBL','TUBGCP2','FRA10AC1','RBBP8','WDR11','UBE3A','GMNN','UNC80','CTCF','RTTN','HDAC8','NUP214','TTI1','RMI1','ATP9A','RUSC2','TRAIP','FBRSL1','ERCC6','MORC2','TUBGCP4','PQBP1','AP4E1','TNPO2','PRUNE1','NSRP1','HMGB1'],
        'Diencephalon':
    ['ZNHIT3','PLK4','MFSD2A','HIKESHI','LMNB2','NCAPD2','SASS6','RMI1','TUBG1','CEP135','CRIPT','MSMO1','LARP7','RRP7A','SLC9A6','DPP6','LAGE3','KIF11','PRIM1','SVBP','NCAPH','FANCB','TPRKB','TRMT10A','MPLKIP','UNC80','HMGB1','CENPE','POC1A','FANCD2','ORC6','DYNC1I2','DHCR7','MED11','BUB1B','UFC1','CIT','KIF14','CTNNB1','ASPM','BUB1','KNL1','CEP55','HIST1H4C','CDK6','METTL5','CENPF','ATP1A2','NUF2','CKAP2L'],
        'Midbrain':
    ['EIF2S3','ARCN1','MINPP1','NSRP1','ERCC4','NBN','RAD21','MSMO1','ATP1A2','TUBG1','DHCR7','SMO','UBA5','UFC1','TRMT1','QARS','RPL10','AP4S1','RAD51C','STAMBP','CRIPT','ATP6V0A1','PUF60','WDR73','CEP63','NIN','NAA20','COPB2','PPFIBP1','ZNHIT3','TPRKB','GINS2','TRA2B','COASY','SVBP'],
        'Hindbrain':
    ['SMC1A','PPP1R35','FBRSL1','KMT2B','RPL10','NSMCE2','SMC3','TPR','MCPH1','ATP9A','CHAMP1','PDHA1','TUBGCP2','XRCC4','PTPN23','PUF60','PSMC3','ERCC6','ZNF335','COPB2','TTC5','MECP2','DPM1','CHKA','CAMSAP1','SLC1A4','COPB1','QARS','CRIPT','RAD51C','LMNB1','DYNC1I2','TUBG1','TRAPPC12','TUBGCP4','SLX4','CTCF','PQBP1','ZPR1','TRA2B','CPSF3','MSMO1','ERCC8','ARF3','ATP1A2','AP4E1','UFM1','WDR4','RING1','ANKLE2','COG3','SMARCA5','METTL5','MPLKIP','PLAA','TP53RK','MED11','LARP7','EIF5A','KIF1BP','WDR73','UBA5','GRM7','SVBP','TPRKB','DPP6','DOHH','YIF1B','PALB2','HIKESHI','STAMBP','ERCC4','LAGE3','ARCN1','MORC2','NIN','ZNF668','AP4S1','PPP1R15B','MINPP1','TRAPPC10','CTU2','ZNHIT3','PPIL1','TSEN15','TRMT10A','NBN'],
        'Pons':
        ['TSEN15','NAA20','YIF1B','MPLKIP','LIG4','MED11','TPRKB','UFC1','RPL10','NIN'],
        'Cerebellum':
    ['TNPO2','TRAPPC9','DIAPH1','UFM1','COPB2','WDFY3','DHCR7','MYCN','ATP1A2','MSMO1','TUBGCP2','HDAC8','WDR4','KMT2B','RPL10','SLC9A6','MECP2','ARF3','CAMSAP1','ATP6V0A1','ZNF335','PTPN23','TRAPPC12','FANCC','ZNHIT3','TUBGCP4','ATRIP','GINS3','TSEN54','ACBD6','PDHA1','CEP57','PALB2','ORC4','CHKA','TTC5','MCPH1','PPFIBP1','INTS11','TUBGCP6','UBE3A','PHC1','PDCD6IP','DYNC1I2','UNC80','ZPR1','WDR37','RUSC2','RING1','SLX4','AP4B1','SVBP','DOHH','FRA10AC1','EXOC7','CTU2','CASK','FANCG','MORC2','VPS50','EOMES','IGF1R','DPP6','TRIO','AP4M1','TRMT1','YIF1B','SLC1A4','AP4S1','NIN','CTSF','TRAPPC6B','ATP9A','ZNF668','ERCC6'],
        'Medulla':
    ['UNC80','HMGB1','NAA20','GPT2','UFC1','GRM7','ZPR1','DHCR7','KIF1BP','EIF5A','ARF3','COASY','RPL10','HIKESHI','CRIPT','TSEN15','TAF13','DPP6','MINPP1','MSMO1','WLS','MED11','METTL5','TCF4','YIF1B']
    }
    return (regions,)


@app.cell
def _(
    microcephaly_NPCs,
    microcephaly_glioblasts,
    microcephaly_radial_glia,
    plot_gene_overlap,
):
    plot_gene_overlap({'NPCs':microcephaly_NPCs, 'Radial glia':microcephaly_radial_glia, 'Glioblast':microcephaly_glioblasts}, 'Microcephaly (All layers)')
    return


@app.cell
def _(
    microcephaly_NPCs,
    microcephaly_glioblasts,
    microcephaly_radial_glia,
    plot_gene_overlap,
):
    plot_gene_overlap({'NPCs':microcephaly_NPCs, 'Radial glia':microcephaly_radial_glia, 'Glioblast':microcephaly_glioblasts}, 'Microcephaly (All layers)')
    return


@app.cell
def _(
    microcephaly_G1,
    microcephaly_G2M,
    microcephaly_S,
    plot_gene_overlap,
    regions,
):
    plot_gene_overlap({'G1':microcephaly_G1, 'S':microcephaly_S, 'G2M':microcephaly_G2M, 'Hindbrain': regions['Hindbrain']}, 'Microcephaly (All layers)')
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
