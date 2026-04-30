"""
GSEA enrichment with Tau pre-filtering.

Drop-in replacement for run_gsea() in search_enrichment_gsea.py.
For each condition:
  1. Load GES scores (prerank metric).
  2. Load global Tau scores for the matching column (from tau_pipeline.py output).
  3. Retain only genes at or above `tau_percentile` (default: 90th → top 10%).
  4. Apply the GES score threshold within that filtered set.
  5. Run prerank GSEA.

Results are written to GSEA_tau_filtered/ instead of GSEA/ so they coexist
with unfiltered runs.
"""

from pathlib import Path

import gseapy as gp
import pandas as pd
import statsmodels.stats.multitest as smm

from search_enrichment_gsea import plot_enhanced_gsea


# ---------------------------------------------------------------------------
# Tau helpers
# ---------------------------------------------------------------------------

def load_tau_scores(tau_scores_dir: Path, column: str) -> pd.DataFrame | None:
    """Load tau_scores_{column}.csv produced by tau_pipeline.py."""
    path = Path(tau_scores_dir) / f"tau_scores_{column}.csv"
    if not path.exists():
        print(f"⚠️ Tau scores not found for column '{column}': {path}")
        return None
    return pd.read_csv(path)


def top_tau_gene_set(tau_df: pd.DataFrame, percentile: float = 90.0) -> set:
    """
    Return the set of genes at or above `percentile` of the Tau distribution.
    E.g. percentile=90 keeps the top 10% most cell-type-specific genes.
    """
    threshold = tau_df["tau"].quantile(percentile / 100.0)
    return set(tau_df.loc[tau_df["tau"] >= threshold, "gene"])


# ---------------------------------------------------------------------------
# Main enrichment function
# ---------------------------------------------------------------------------

def run_gsea_tau_filtered(
    ges_score_path,
    tau_scores_dir,
    gmt_file,
    column_conditions,
    ges_score_threshold,
    out_folder,
    figs_folder,
    tau_percentile: float = 90.0,
) -> Path | None:
    """
    Run GSEA prerank enrichment, restricting the gene universe to those in the
    top (100 - tau_percentile)% by global Tau specificity before ranking.

    Parameters
    ----------
    ges_score_path : path-like
        Folder containing per-condition GES CSVs
        (data/ges_spec_{column}_{condition}.csv).
    tau_scores_dir : path-like
        Folder containing tau_scores_{column}.csv files
        (data/ subfolder of a tau_pipeline.py run).
    gmt_file : str
        Path to the .gmt gene-set file.
    column_conditions : dict
        {column: [condition, ...]} mapping.
    ges_score_threshold : float
        Minimum GES score; genes below this are excluded after the Tau filter.
    out_folder : path-like
        Root folder for enrichment result CSVs.
    figs_folder : path-like
        Root folder for figures.
    tau_percentile : float
        Tau quantile cutoff (default 90 → keep genes with tau >= 90th percentile).

    Returns
    -------
    Path to the summary CSV, or None if no results were produced.
    """
    tau_scores_dir = Path(tau_scores_dir)
    results_folder = Path(out_folder) / "GSEA_tau_filtered"
    fig_dir = Path(figs_folder) / "GSEA_tau_filtered"
    results_folder.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    gmt_name = Path(gmt_file).stem
    top_pct = 100.0 - tau_percentile
    print(f"\nGMT set: {gmt_name}")
    print(f"Tau filter: top {top_pct:.0f}% (tau >= {tau_percentile:.0f}th percentile)\n")

    # Cache top-tau gene sets per column to avoid repeated I/O and computation
    tau_cache: dict[str, set | None] = {}

    summary_rows = []

    for column, condition_list in column_conditions.items():
        if column not in tau_cache:
            tau_df = load_tau_scores(tau_scores_dir, column)
            if tau_df is None:
                tau_cache[column] = None
            else:
                gene_set = top_tau_gene_set(tau_df, percentile=tau_percentile)
                tau_cache[column] = gene_set
                print(
                    f"Column '{column}': {len(gene_set)} genes retained "
                    f"(tau >= {tau_percentile:.0f}th percentile)"
                )

        top_tau_set = tau_cache[column]

        for condition in condition_list:
            print(f"\nRunning Tau-filtered GSEA: {column} → {condition}")

            ges_path = Path(ges_score_path) / "data" / f"ges_spec_{column}_{condition}.csv"
            if not ges_path.exists():
                print(f"  ⚠️ Missing GES file: {ges_path} — skipping.")
                continue

            ges_df = pd.read_csv(ges_path)

            # --- Tau filter ---
            if top_tau_set is not None:
                n_before = len(ges_df)
                ges_df = ges_df[ges_df["gene"].isin(top_tau_set)].copy()
                print(f"  Tau filter:  {n_before} → {len(ges_df)} genes")
            else:
                print("  ⚠️ No Tau scores available for this column — Tau filter skipped.")

            # --- GES threshold ---
            n_before = len(ges_df)
            ges_df = ges_df[ges_df["ges_score"] > ges_score_threshold]
            print(f"  GES filter:  {n_before} → {len(ges_df)} genes (threshold={ges_score_threshold})")

            if ges_df.empty:
                print(f"  ⚠️ No genes pass both filters — skipping.")
                continue

            ranking = ges_df.set_index("gene")["ges_score"]

            cond_dir = results_folder / column / str(condition)
            raw_dir = cond_dir / "gsea_raw"
            cond_dir.mkdir(parents=True, exist_ok=True)

            try:
                gsea_res = gp.prerank(
                    rnk=ranking,
                    gene_sets=gmt_file,
                    outdir=str(raw_dir),
                    min_size=2,
                    max_size=2500,
                    seed=6,
                )
            except LookupError:
                print("  ⚠️ Not enough genes overlap with gene sets — skipping.")
                gsea_res = None

            if gsea_res is None:
                continue

            out_csv = cond_dir / "gsea_results.csv"
            gsea_res.res2d.to_csv(out_csv)
            print(f"  ✔ Saved results → {out_csv}")

            top = gsea_res.res2d.iloc[0]
            summary_rows.append({
                "column": column,
                "condition": condition,
                "term": top["Term"],
                "NES": top["NES"],
                "NOM p-val": top["NOM p-val"],
                "FDR q-val": top["FDR q-val"],
                "FWER p-val": top["FWER p-val"],
                "Tag %": top["Tag %"],
                "Gene %": top["Gene %"],
                "Lead_genes": top["Lead_genes"],
                "tau_percentile_cutoff": tau_percentile,
            })

            term = gsea_res.res2d["Term"].iloc[0]
            if term in gsea_res.results:
                plot_out = plot_enhanced_gsea(gsea_res, term, condition, raw_dir)
                print(f"  ✔ Saved plot → {plot_out}")
            else:
                print(f"  ⚠️ Term '{term}' not in result dict — skipping plot.")

    if not summary_rows:
        print("\n⚠️ No GSEA results produced — summary not written.")
        return None

    final_summary = pd.DataFrame(summary_rows)

    _, pvals_corr = smm.multipletests(
        final_summary["FDR q-val"].astype(float), method="fdr_bh"
    )[:2]
    final_summary["FDR q-val (BH corrected)"] = pvals_corr

    summary_path = results_folder / "GSEA_tau_filtered_summary.csv"
    final_summary.to_csv(summary_path, index=False)

    print(f"\n📄 Summary written: {summary_path}")
    print("🎉 Tau-filtered GSEA finished.")
    return summary_path
