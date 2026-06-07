"""
hypergeometric_enrichment.py

Standalone hypergeometric gene set enrichment — no GES scores, no tau, no GSEA.

For each dataset → column (e.g. CellClass, Region) → condition (e.g. Neuron):
  Background (N) = all unique genes across ALL conditions in that column.
  K             = input gene list ∩ background.
  n             = genes expressed in this specific condition.
  k             = input genes ∩ condition genes.
  p-value       = P(X ≥ k) = hypergeom.sf(k-1, N, K, n)   [one-sided enrichment]

BH correction applied within each column.

Outputs
-------
  results/hypergeometric/{dataset_name}_{YYYYMMDD}/
    metadata/pipeline_output.log
    {dataset}/
      {column}_results.csv
      {column}_barplot.png
    summary_all_results.csv

Usage
-----
  python modules/hypergeometric_enrichment.py <config.yaml> <gene_list.csv>

  Gene list: CSV with a 'gene' column, or plain text (one gene per line).
  Config:    same structure as gene_clustering configs
             (dataset_name, ndd_gene_modules_folder_root, output_folder, datasets).
"""

import argparse
import contextlib
import datetime
import sys
from pathlib import Path
from typing import NamedTuple

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.stats.multitest as smm
import yaml
from scipy.stats import hypergeom


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def _log_to_file(log_path: Path):
    class _Tee:
        def __init__(self, *streams):
            self._streams = streams
        def write(self, s):
            for st in self._streams:
                st.write(s)
        def flush(self):
            for st in self._streams:
                st.flush()
        @property
        def encoding(self):
            return getattr(self._streams[0], "encoding", "utf-8")

    orig = sys.stdout
    with open(log_path, "w", encoding="utf-8") as fh:
        sys.stdout = _Tee(orig, fh)
        try:
            yield
        finally:
            sys.stdout = orig


# ---------------------------------------------------------------------------
# Data structure
# ---------------------------------------------------------------------------

class DatasetSpec(NamedTuple):
    label: str
    ges_folder: Path
    column_conditions: dict


# ---------------------------------------------------------------------------
# Config loader
# ---------------------------------------------------------------------------

def _load_config(config_path: str):
    config_path = Path(config_path).resolve()
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    root = Path(cfg["ndd_gene_modules_folder_root"]).resolve()

    def resolve(p):
        p = Path(str(p))
        return p if p.is_absolute() else (root / p).resolve()

    cfg["output_folder"] = resolve(cfg.get("output_folder", "results/hypergeometric/"))

    if "datasets" in cfg:
        specs = []
        for label, ds in cfg["datasets"].items():
            specs.append(DatasetSpec(
                label=label,
                ges_folder=resolve(ds["ges_results_folder"]),
                column_conditions=ds["column_conditions"],
            ))
        return cfg, specs

    # Single-dataset fallback
    specs = [DatasetSpec(
        label=cfg["dataset_name"],
        ges_folder=resolve(cfg["ges_results_folder"]),
        column_conditions=cfg["column_conditions"],
    )]
    return cfg, specs


# ---------------------------------------------------------------------------
# Gene list loader
# ---------------------------------------------------------------------------

def _load_gene_list(path) -> set:
    path = Path(path)
    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path)
        col = "gene" if "gene" in df.columns else df.columns[0]
        return set(df[col].dropna().astype(str).str.strip())
    genes = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                genes.add(line.split("\t")[0].strip())
    return genes


# ---------------------------------------------------------------------------
# Load genes per condition from GES CSV, filtered by a minimum GES score
# ---------------------------------------------------------------------------

def _condition_genes(ges_folder: Path, column: str, condition, threshold: float) -> set:
    fname = ges_folder / "data" / f"ges_spec_{column}_{condition}.csv"
    if not fname.exists():
        return set()
    df = pd.read_csv(fname)
    if "gene" not in df.columns:
        return set()
    if "ges_score" in df.columns:
        df = df[df["ges_score"] >= threshold]
    return set(df["gene"].dropna().astype(str))


# ---------------------------------------------------------------------------
# Hypergeometric test for one column
# ---------------------------------------------------------------------------

def _test_column(
    dataset_label: str,
    column: str,
    conditions: list,
    ges_folder: Path,
    disease_genes: set,
    threshold: float,
) -> pd.DataFrame | None:
    """
    Run hypergeometric test for every condition in a column.

    Background = union of all genes (GES >= threshold) across all conditions in this column.
    Each condition is tested against that shared background.
    BH correction is applied across all conditions in the column.
    """
    cond_genes: dict[str, set] = {}
    for cond in conditions:
        g = _condition_genes(ges_folder, column, cond, threshold)
        if g:
            cond_genes[str(cond)] = g
        else:
            print(f"    ⚠️  No genes found for {column}/{cond} at threshold {threshold} — skipping.")

    if len(cond_genes) < 2:
        print(f"    ⚠️  Fewer than 2 conditions with data in {column} — skipping column.")
        return None

    background     = set().union(*cond_genes.values())
    disease_in_bg  = disease_genes & background
    N              = len(background)
    K              = len(disease_in_bg)

    print(f"    background={N:,} | disease in bg={K}")

    if K == 0:
        print(f"    ⚠️  No disease genes in background for column '{column}' — skipping.")
        return None

    rows = []
    for condition, genes in cond_genes.items():
        n       = len(genes)
        overlap = disease_in_bg & genes
        k       = len(overlap)

        expected = (K * n) / N if N > 0 else 0.0
        if expected > 0:
            fold_enrichment = k / expected
        elif k > 0:
            fold_enrichment = np.inf
        else:
            fold_enrichment = 1.0

        pval = float(hypergeom.sf(k - 1, N, K, n)) if k > 0 else 1.0

        rows.append({
            "dataset":                 dataset_label,
            "column":                  column,
            "condition":               condition,
            "n_background":            N,
            "n_disease_in_background": K,
            "n_condition_genes":       n,
            "n_overlap":               k,
            "expected_overlap":        round(expected, 4),
            "fold_enrichment":         fold_enrichment,
            "pvalue":                  pval,
            "FDR_qval_BH":             np.nan,
            "is_significant":          False,
            "overlap_genes":           ";".join(sorted(overlap)),
        })

    df = pd.DataFrame(rows)

    pvals        = df["pvalue"].astype(float).values
    _, corrected = smm.multipletests(pvals, method="fdr_bh")[:2]
    df["FDR_qval_BH"]  = corrected
    df["is_significant"] = corrected < 0.05

    return df


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def _fmt_pvals(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in ("pvalue", "FDR_qval_BH"):
        if col in df.columns:
            df[col] = df[col].apply(
                lambda v: f"{float(v):.6e}" if pd.notna(v) else v
            )
    return df


def _safe_fe(v) -> str:
    try:
        fv = float(v)
        return f"FE={fv:.2f}" if np.isfinite(fv) else "FE=∞"
    except (ValueError, TypeError):
        return ""


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _plot_column(df: pd.DataFrame, column: str, out_path: Path) -> None:
    df = df.copy()

    def _neg_log10(v):
        fv = float(v)
        return -np.log10(fv) if fv > 0 else -np.log10(1e-300)

    df["neg_log10_fdr"] = df["FDR_qval_BH"].apply(_neg_log10)
    df = df.sort_values("neg_log10_fdr", ascending=True)

    fig, ax = plt.subplots(figsize=(8, max(3, len(df) * 0.5 + 1.2)))
    colors = ["#d73027" if s else "#4575b4" for s in df["is_significant"]]
    bars = ax.barh(
        df["condition"].astype(str),
        df["neg_log10_fdr"],
        color=colors, edgecolor="white", height=0.65,
    )

    ax.axvline(-np.log10(0.05), color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("-log₁₀(FDR)", fontsize=11)
    ax.set_title(
        f"Hypergeometric enrichment — {column}\n"
        f"(n_disease={int(df['n_disease_in_background'].iloc[0])}, "
        f"background={int(df['n_background'].iloc[0])})",
        fontsize=11, fontweight="bold",
    )

    for bar, (_, row) in zip(bars, df.iterrows()):
        label = _safe_fe(row["fold_enrichment"])
        n_ov  = int(row["n_overlap"])
        ax.text(
            bar.get_width() + 0.05,
            bar.get_y() + bar.get_height() / 2,
            f"{label}  (k={n_ov})",
            va="center", ha="left", fontsize=8,
        )

    legend_handles = [
        mpatches.Patch(color="#d73027", label="FDR < 0.05"),
        mpatches.Patch(color="#4575b4", label="FDR ≥ 0.05"),
        plt.Line2D([0], [0], color="gray", linestyle="--", linewidth=1, label="FDR = 0.05"),
    ]
    ax.legend(handles=legend_handles, loc="lower right", fontsize=8)
    ax.set_xlim(left=0)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()


def _plot_summary(summary: pd.DataFrame, out_path: Path, gene_list_name: str) -> None:
    """Dot plot: datasets × conditions coloured by -log10(FDR), size by fold enrichment."""
    df = summary.copy()

    def _neg_log10(v):
        fv = float(v)
        return -np.log10(fv) if fv > 0 else -np.log10(1e-300)

    df["neg_log10_fdr"] = df["FDR_qval_BH"].apply(_neg_log10)
    df["fe_plot"] = df["fold_enrichment"].apply(
        lambda v: min(float(v), 20) if np.isfinite(float(v)) else 20
    )

    # y-axis: dataset__column__condition
    df["y_label"] = df["dataset"] + " | " + df["column"] + " | " + df["condition"].astype(str)
    df = df.sort_values(["dataset", "column", "neg_log10_fdr"], ascending=[True, True, False])

    fig_h = max(6, len(df) * 0.35 + 1)
    fig, ax = plt.subplots(figsize=(8, fig_h))

    sc = ax.scatter(
        df["neg_log10_fdr"],
        range(len(df)),
        c=df["neg_log10_fdr"],
        s=df["fe_plot"] * 30 + 20,
        cmap="YlOrRd",
        edgecolors="gray",
        linewidths=0.4,
        zorder=2,
    )

    ax.axvline(-np.log10(0.05), color="gray", linestyle="--", linewidth=1, zorder=1)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["y_label"], fontsize=7)
    ax.set_xlabel("-log₁₀(FDR)", fontsize=11)
    ax.set_title(
        f"Hypergeometric enrichment summary — {gene_list_name}",
        fontsize=11, fontweight="bold",
    )
    plt.colorbar(sc, ax=ax, label="-log₁₀(FDR)", shrink=0.5)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Per-dataset pipeline
# ---------------------------------------------------------------------------

def _run_dataset(
    spec: DatasetSpec,
    disease_genes: set,
    out_dir: Path,
    threshold: float,
) -> list[pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    all_dfs = []

    for column, conditions in spec.column_conditions.items():
        print(f"\n  Column: {column}  ({len(conditions)} conditions)")
        df = _test_column(spec.label, column, conditions, spec.ges_folder, disease_genes, threshold)
        if df is None:
            continue

        safe_col = column.replace(" ", "_")
        csv_path = out_dir / f"{safe_col}_results.csv"
        _fmt_pvals(df).to_csv(csv_path, index=False)
        print(f"    ✔ {csv_path.name}")

        fig_path = out_dir / f"{safe_col}_barplot.png"
        _plot_column(df, column, fig_path)
        print(f"    ✔ {fig_path.name}")

        all_dfs.append(df)

    return all_dfs


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run_hypergeometric_enrichment(config_path: str, gene_list_path: str) -> None:
    cfg, specs = _load_config(config_path)

    out_root       = Path(cfg["output_folder"])
    date_str       = datetime.datetime.now().strftime("%Y%m%d")
    gene_list_name = Path(gene_list_path).stem
    version        = cfg.get("version", "")
    folder_name    = f"{gene_list_name}_{version}_{date_str}" if version else f"{gene_list_name}_{date_str}"
    run_dir        = out_root / folder_name
    run_dir.mkdir(parents=True, exist_ok=True)

    meta_dir = run_dir / "metadata"
    meta_dir.mkdir(parents=True, exist_ok=True)

    with _log_to_file(meta_dir / "pipeline_output.log"):
        print("=" * 62)
        print(f"  Hypergeometric Enrichment — {gene_list_name}  [{version}]")
        print("=" * 62)
        print(f"  Gene list:   {gene_list_path}")
        print(f"  Output:      {run_dir}")
        print(f"  Datasets:")
        for s in specs:
            print(f"    · {s.label:<16} {s.ges_folder}")
        print("=" * 62)

        threshold = float(cfg.get("ges_score_threshold", 1.0))
        print(f"  GES threshold:   >= {threshold}  (defines specifically expressed genes per condition)")

        disease_genes = _load_gene_list(gene_list_path)
        print(f"\n🧬 Disease genes loaded: {len(disease_genes)}")

        all_results: list[pd.DataFrame] = []

        for spec in specs:
            print(f"\n{'─'*50}")
            print(f"  Dataset: {spec.label}")
            print(f"{'─'*50}")
            dfs = _run_dataset(spec, disease_genes, run_dir / spec.label, threshold)
            all_results.extend(dfs)

        if not all_results:
            print("\n⚠️  No results produced.")
            return

        summary = pd.concat(all_results, ignore_index=True)
        summary_path = run_dir / "summary_all_results.csv"
        _fmt_pvals(summary).to_csv(summary_path, index=False)
        print(f"\n📄 Summary saved: {summary_path}")

        summary_plot = run_dir / "summary_dotplot.png"
        _plot_summary(summary, summary_plot, gene_list_name)
        print(f"📊 Summary dot plot: {summary_plot}")

        n_sig = summary["is_significant"].sum()
        print(f"\n🎉 Done! {n_sig}/{len(summary)} conditions are significant (FDR < 0.05).")
        print(f"   Results: {run_dir}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Hypergeometric enrichment: test whether a gene list is enriched "
            "in each condition, compared to all conditions in the same column."
        )
    )
    parser.add_argument("config",     help="YAML config file (gene_clustering format).")
    parser.add_argument("gene_list",  help="Gene list: CSV with 'gene' column or plain text.")
    args = parser.parse_args()
    run_hypergeometric_enrichment(args.config, args.gene_list)
