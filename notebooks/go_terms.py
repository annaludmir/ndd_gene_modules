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

    def get_microcephaly_gene_lists(summary_csv):
        df = pd.read_csv(summary_csv)

        # condition label
        df["condition"] = (
            df["column"] + " - " + df["condition"]
        )

        gene_lists = {}

        for cond, sub in df.groupby("condition"):
            genes = (
                sub["Lead_genes"]
                .dropna()
                .str.split(";")
                .explode()
                .str.strip()
                .unique()
            )
            gene_lists[cond] = list(genes)

        return gene_lists
    return get_microcephaly_gene_lists, pd


@app.cell
def _(pd):
    import re
    from pathlib import Path
    import gseapy as gp


    def run_go_enrichment_per_condition(
        gene_lists,
        outdir="microcephaly_GO_enrichment",
        organism="Human",
        gene_sets=None,
        cutoff=0.05,
        min_genes=5,
        save_all_results=True,
    ):
        """
        Runs Enrichr ORA per condition. Never crashes when no terms pass cutoff.

        Returns
        -------
        dict[str, pd.DataFrame]
            condition -> results dataframe (may be empty)
        """
        if gene_sets is None:
            gene_sets = [
                "GO_Biological_Process_2021",
                    "GO_Molecular_Function_2021",
                    "GO_Cellular_Component_2021",
                    "KEGG_2016",
                    "KEGG_2021_Human",
                    "WikiPathway_2021_Human"
            ]

        Path(outdir).mkdir(parents=True, exist_ok=True)

        results_by_condition = {}

        for condition, genes in gene_lists.items():
            genes = [g for g in genes if isinstance(g, str) and g.strip()]
            genes = list(dict.fromkeys(genes))  # de-dup, keep order

            safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", condition).strip("_")
            cond_outdir = Path(outdir) / safe_name
            cond_outdir.mkdir(parents=True, exist_ok=True)

            if len(genes) < min_genes:
                print(f"Skipping {condition} (too few genes: {len(genes)})")
                results_by_condition[condition] = pd.DataFrame()
                continue

            print(f"Running enrichment: {condition} (n_genes={len(genes)})")

            # IMPORTANT: outdir=None prevents gseapy from auto-plotting/saving figures
            enr = gp.enrichr(
                gene_list=genes,
                organism=organism,
                gene_sets=gene_sets,
                outdir=None,
                cutoff=cutoff,
                no_plot=True,   # <-- avoids the crash
            )

            res = enr.results if enr is not None and hasattr(enr, "results") else None
            if res is None:
                res = pd.DataFrame()

            results_by_condition[condition] = res

            # Save results even if empty
            if save_all_results:
                res.to_csv(cond_outdir / "enrichr_results.csv", index=False)

            # Also save the genes used (great for reproducibility)
            pd.Series(genes, name="gene").to_csv(cond_outdir / "leading_genes_used.csv", index=False)

            # If you want: save only significant hits (when they exist)
            if not res.empty and "Adjusted P-value" in res.columns:
                sig = res[res["Adjusted P-value"] <= cutoff].copy()
                sig.to_csv(cond_outdir / f"significant_hits_fdr_{cutoff}.csv", index=False)

                if sig.empty:
                    print(f"  No significant terms at cutoff={cutoff}")
            else:
                print("  No results returned (or unexpected format).")

        return results_by_condition
    return (run_go_enrichment_per_condition,)


@app.cell
def _(get_microcephaly_gene_lists, run_go_enrichment_per_condition):
    summary_csv = "/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/microcephaly/microcephaly_genes_all_layers_threshold_1_20260110/data/enrichment_results/GSEA/GSEA_final_summary.csv"
    gene_lists = get_microcephaly_gene_lists(summary_csv)
    go_results = run_go_enrichment_per_condition(gene_lists)
    return (go_results,)


@app.cell
def _(pd):

    def collect_top_terms(go_results, n=5, sort_col="Adjusted P-value", per_gene_set=False):
        """
        go_results: dict[str, pd.DataFrame]  # condition -> enrichr results
        n: top N terms to keep
        sort_col: column used for ranking (default: Adjusted P-value)
        per_gene_set: if True, take top n terms *per Gene_set* within each condition
        """
        rows = []

        for condition, df in go_results.items():
            if df is None or not isinstance(df, pd.DataFrame) or df.empty:
                continue

            # pick a usable sort column
            if sort_col not in df.columns:
                # fall back options
                fallback = next((c for c in ["Adjusted P-value", "P-value", "Combined Score"] if c in df.columns), None)
                if fallback is None:
                    continue
                use_sort = fallback
            else:
                use_sort = sort_col

            dfx = df.copy()

            # ensure sortable numeric
            dfx[use_sort] = pd.to_numeric(dfx[use_sort], errors="coerce")

            # sort direction: p-values ascending; combined score descending
            ascending = (use_sort in ["Adjusted P-value", "P-value"])

            if per_gene_set and "Gene_set" in dfx.columns:
                top_df = (
                    dfx.sort_values(use_sort, ascending=ascending)
                       .groupby("Gene_set", dropna=False)
                       .head(n)
                )
            else:
                top_df = dfx.sort_values(use_sort, ascending=ascending).head(n)

            for _, r in top_df.iterrows():
                rows.append({
                    "condition": condition,
                    "gene_set": r.get("Gene_set", None),
                    "term": r.get("Term", None),
                    "p_value": r.get("P-value", None),
                    "adj_p": r.get("Adjusted P-value", None),
                    "combined_score": r.get("Combined Score", None),
                    "overlap": r.get("Overlap", None),
                    "genes": r.get("Genes", None),
                    "ranked_by": use_sort
                })

        return pd.DataFrame(rows)
    return (collect_top_terms,)


@app.cell
def _(go_results):
    go_results
    return


@app.cell
def _(collect_top_terms, go_results):
    collect_top_terms(go_results)
    return


@app.cell
def _():
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    def load_enrichr_tables(phase_to_csv: dict, phase_col="Phase") -> pd.DataFrame:
        """phase_to_csv: dict like {'G1': 'g1.csv', 'S': 's.csv', ...}"""
        dfs = []
        for phase, path in phase_to_csv.items():
            d = pd.read_csv(path)
            d[phase_col] = phase
            dfs.append(d)
        return pd.concat(dfs, ignore_index=True)

    def bubble_go_by_phase(
        df: pd.DataFrame,
        phase_order: list,
        *,
        term_col: str = "Term",
        fdr_col: str = "Adjusted P-value",
        phase_col: str = "Phase",
        score_col_for_color: str | None = None,  # e.g. "Combined Score" or None
        top_terms: int = 30,
        figsize=(16, 5),
    ):
        """
        Bubble plot: y=phase, x=GO term, size=-log10(FDR).
        Terms are ordered by the phase where they are strongest (then by strength).
        """

        d = df.copy()
        d[fdr_col] = pd.to_numeric(d[fdr_col], errors="coerce")
        d = d.dropna(subset=[term_col, fdr_col, phase_col])

        # score used for size
        d["mlog10_fdr"] = -np.log10(np.clip(d[fdr_col].values, 1e-300, 1.0))

        # Keep only phases in desired order
        d = d[d[phase_col].isin(phase_order)].copy()
        d[phase_col] = pd.Categorical(d[phase_col], categories=phase_order, ordered=True)

        # Pick top terms globally (by best/strongest FDR across phases)
        term_strength = d.groupby(term_col)["mlog10_fdr"].max().sort_values(ascending=False)
        chosen_terms = term_strength.head(top_terms).index.tolist()
        d = d[d[term_col].isin(chosen_terms)].copy()

        # Order terms by "peak phase" then by peak strength within that phase
        peak = (
            d.sort_values("mlog10_fdr", ascending=False)
             .groupby(term_col)
             .first()[[phase_col, "mlog10_fdr"]]
        )
        # Convert phase to order index
        phase_rank = {p:i for i,p in enumerate(phase_order)}
        peak["phase_rank"] = peak[phase_col].map(phase_rank)

        term_order = (
            peak.sort_values(["phase_rank", "mlog10_fdr"], ascending=[True, False])
                .index.tolist()
        )

        # Prepare grid mapping to numeric x/y
        x_map = {t:i for i,t in enumerate(term_order)}
        y_map = {p:i for i,p in enumerate(phase_order)}

        d["x"] = d[term_col].map(x_map)
        d["y"] = d[phase_col].map(y_map)

        # Dot size scaling
        # (tune these if needed)
        min_size = 20
        max_size = 400
        s = d["mlog10_fdr"].values
        if np.nanmax(s) > 0:
            sizes = min_size + (s / np.nanmax(s)) * (max_size - min_size)
        else:
            sizes = np.full_like(s, min_size, dtype=float)

        # Dot color
        if score_col_for_color is None:
            c = d["mlog10_fdr"].values
            c_label = "-log10(FDR)"
        else:
            d[score_col_for_color] = pd.to_numeric(d[score_col_for_color], errors="coerce")
            c = d[score_col_for_color].values
            c_label = score_col_for_color

        plt.figure(figsize=figsize)
        sc = plt.scatter(d["x"], d["y"], s=sizes, c=c)

        plt.yticks(range(len(phase_order)), phase_order)
        plt.xticks(range(len(term_order)), term_order, rotation=75, ha="right")

        plt.gca().invert_yaxis()  # often nicer: top phase at top
        plt.colorbar(sc, label=c_label)
        plt.title("GO terms across cell cycle phases (dot size = -log10(FDR))")
        plt.tight_layout()
        plt.show()

        return {"data": d, "term_order": term_order}

    return bubble_go_by_phase, load_enrichr_tables, np, pd, plt


@app.cell
def _(bubble_go_by_phase, load_enrichr_tables):
    phase_to_csv = {
        "G1": "ndd_gene_modules/GO_terms/microcephaly_GO_enrichment/CellCyclePhase_-_G1/significant_hits_fdr_0.05.csv",
        "S": "ndd_gene_modules/GO_terms/microcephaly_GO_enrichment/CellCyclePhase_-_S/significant_hits_fdr_0.05.csv",
        "G2M": "ndd_gene_modules/GO_terms/microcephaly_GO_enrichment/CellCyclePhase_-_G2M/significant_hits_fdr_0.05.csv",
        "PostM": "ndd_gene_modules/GO_terms/microcephaly_GO_enrichment/CellCyclePhase_-_PostM/significant_hits_fdr_0.05.csv",
        "Non-cycling": "ndd_gene_modules/GO_terms/microcephaly_GO_enrichment/CellCyclePhase_-_Non-cycling/significant_hits_fdr_0.05.csv"
    }
    df_all = load_enrichr_tables(phase_to_csv)

    out = bubble_go_by_phase(df_all, phase_order=["G1","S","G2M","PostM","Non-cycling"], top_terms=25, score_col_for_color="Combined Score")
    return df_all, out


@app.cell
def _(np, pd, plt):
    def _parse_gene_list(s):
        if pd.isna(s):
            return set()
        return set(g.strip() for g in str(s).split(";") if g.strip())

    def go_gene_shift_heatmap(
        df: pd.DataFrame,
        phase_order: list,
        term_order: list,
        *,
        term_col="Term",
        phase_col="Phase",
        genes_col="Genes",
        figsize=(14, 6),
    ):
        """
        Heatmap: rows=GO terms, cols=phase transitions (G1→S, S→G2M, ...)
        values = Jaccard similarity of overlap genes for each term across phases.
        """

        d = df.copy()
        d = d[d[phase_col].isin(phase_order) & d[term_col].isin(term_order)].copy()

        # Build dict: (term, phase) -> gene set
        gene_map = {}
        for (term, phase), sub in d.groupby([term_col, phase_col]):
            # if multiple rows per term/phase, union them
            gs = set()
            for v in sub[genes_col].tolist():
                gs |= _parse_gene_list(v)
            gene_map[(term, phase)] = gs

        transitions = [(phase_order[i], phase_order[i+1]) for i in range(len(phase_order)-1)]
        trans_labels = [f"{a}→{b}" for a,b in transitions]

        mat = np.full((len(term_order), len(transitions)), np.nan, dtype=float)

        for ti, term in enumerate(term_order):
            for j, (a, b) in enumerate(transitions):
                A = gene_map.get((term, a), set())
                B = gene_map.get((term, b), set())
                if len(A) == 0 and len(B) == 0:
                    mat[ti, j] = np.nan
                else:
                    mat[ti, j] = len(A & B) / max(1, len(A | B))

        plt.figure(figsize=figsize)
        im = plt.imshow(mat, aspect="auto", interpolation="nearest")
        plt.colorbar(im, label="Jaccard similarity of overlap genes")

        plt.yticks(range(len(term_order)), term_order)
        plt.xticks(range(len(trans_labels)), trans_labels)

        plt.title("GO-term gene-driver stability across phases")
        plt.tight_layout()
        plt.show()

        return mat

    return (go_gene_shift_heatmap,)


@app.cell
def _(df_all, go_gene_shift_heatmap, out):
    go_gene_shift_heatmap(df_all, ["G1","S","G2M","PostM","Non-cycling"], out["term_order"])
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
