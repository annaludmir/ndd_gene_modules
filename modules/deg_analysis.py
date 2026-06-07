"""
deg_analysis.py

Differential expression analysis (one condition vs rest) for the same
datasets and categories used in the GES pipeline.

For each dataset → column (e.g. CellClass) → condition (e.g. Neuron):
  Wilcoxon rank-sum test (scanpy rank_genes_groups, reference='rest').

Independent of GES scores, tau filtering, GSEA, and hypergeometric tests.

Preprocessing applied identically to the GES pipeline:
  sc.pp.filter_cells(min_genes=200)
  sc.pp.filter_genes(min_cells=3)
  sc.pp.normalize_total()
  sc.pp.log1p()

Usage
-----
  python modules/deg_analysis.py <config.yaml>

Output
------
  results/deg_analysis/{dataset_name}_{YYYYMMDD}/
    metadata/pipeline_output.log
    {dataset}/
      {column}_results.csv      all genes × all conditions (log2fc, pval_adj, …)
      {column}_volcano.png      grid of volcano plots, one panel per condition
    summary_significant.csv     significant DEGs across all datasets and columns
"""

import argparse
import contextlib
import datetime
import sys
import warnings
from pathlib import Path
from typing import NamedTuple

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import scanpy as sc
import yaml


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
    h5ad_path: Path
    column_conditions: dict
    chemistry: str | None      # filter adata.obs["Chemistry"] == chemistry if set


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

    cfg["output_folder"] = resolve(cfg.get("output_folder", "results/deg_analysis/"))

    if "datasets" in cfg:
        specs = []
        for label, ds in cfg["datasets"].items():
            specs.append(DatasetSpec(
                label=label,
                h5ad_path=resolve(ds["h5ad_path"]),
                column_conditions=ds["column_conditions"],
                chemistry=ds.get("chemistry"),
            ))
        return cfg, specs

    # Single-dataset fallback
    specs = [DatasetSpec(
        label=cfg["dataset_name"],
        h5ad_path=resolve(cfg["h5ad_path"]),
        column_conditions=cfg["column_conditions"],
        chemistry=cfg.get("chemistry"),
    )]
    return cfg, specs


# ---------------------------------------------------------------------------
# Data loader — same preprocessing as the GES pipeline
# ---------------------------------------------------------------------------

def _load_and_preprocess(spec: DatasetSpec) -> sc.AnnData:
    print(f"  Loading {spec.h5ad_path.name} ...")
    adata = sc.read_h5ad(spec.h5ad_path)
    print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Chemistry filter
    if spec.chemistry and "Chemistry" in adata.obs.columns:
        before = adata.n_obs
        adata = adata[adata.obs["Chemistry"] == spec.chemistry].copy()
        print(f"  Chemistry filter ({spec.chemistry}): {before:,} → {adata.n_obs:,} cells")

    # Same filters as GES pipeline
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"  After QC filters: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Normalize and log-transform (same as GES pipeline)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    print(f"  Preprocessing done (normalize_total + log1p)")

    return adata


# ---------------------------------------------------------------------------
# DEG for one column
# ---------------------------------------------------------------------------

def _run_column_deg(
    adata: sc.AnnData,
    dataset_label: str,
    column: str,
    conditions: list,
    out_dir: Path,
    logfc_threshold: float,
    pval_threshold: float,
) -> pd.DataFrame | None:

    if column not in adata.obs.columns:
        print(f"    ⚠️  Column '{column}' not in obs — skipping.")
        return None

    obs_vals = adata.obs[column].astype(str)
    available = [str(c) for c in conditions if str(c) in obs_vals.unique()]
    missing   = [str(c) for c in conditions if str(c) not in obs_vals.unique()]
    if missing:
        print(f"    ⚠️  Conditions not found in obs, skipping: {missing}")
    if len(available) < 2:
        print(f"    ⚠️  Fewer than 2 conditions available — skipping column.")
        return None

    # Subset to cells belonging to the specified conditions only
    adata_sub = adata[obs_vals.isin(available)].copy()
    print(f"    {adata_sub.n_obs:,} cells across {len(available)} conditions")

    # Wilcoxon rank-sum, each condition vs all others
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(
            adata_sub,
            groupby=column,
            groups=available,
            reference="rest",
            method="wilcoxon",
            n_genes=adata_sub.n_vars,
            pts=True,
        )

    # Extract results for every condition
    dfs = []
    for condition in available:
        df = sc.get.rank_genes_groups_df(adata_sub, group=condition)
        df = df.rename(columns={
            "names":          "gene",
            "scores":         "wilcoxon_score",
            "logfoldchanges": "log2fc",
            "pvals":          "pvalue",
            "pvals_adj":      "pvalue_adj",
            "pct_nz_group":   "pct_expressing_condition",
            "pct_nz_reference": "pct_expressing_rest",
        })
        df.insert(0, "dataset",   dataset_label)
        df.insert(1, "column",    column)
        df.insert(2, "condition", condition)
        df["is_significant"] = (
            (df["pvalue_adj"] < pval_threshold) &
            (df["log2fc"].abs() >= logfc_threshold)
        )
        dfs.append(df)

    result = pd.concat(dfs, ignore_index=True)

    # Save per-column CSV
    safe_col  = column.replace(" ", "_")
    csv_path  = out_dir / f"{safe_col}_results.csv"
    result.to_csv(csv_path, index=False)
    print(f"    ✔ {csv_path.name}  ({len(result):,} rows)")

    # Volcano plots
    fig_path = out_dir / f"{safe_col}_volcano.png"
    _plot_volcano_grid(result, column, available, logfc_threshold, pval_threshold, fig_path)
    print(f"    ✔ {fig_path.name}")

    return result


# ---------------------------------------------------------------------------
# Volcano plot grid
# ---------------------------------------------------------------------------

def _plot_volcano_grid(
    df: pd.DataFrame,
    column: str,
    conditions: list,
    logfc_threshold: float,
    pval_threshold: float,
    out_path: Path,
) -> None:
    n = len(conditions)
    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols * 4, nrows * 3.8),
                              squeeze=False)
    axes_flat = axes.flatten()

    neg_log_thresh = -np.log10(pval_threshold)

    for idx, condition in enumerate(conditions):
        ax  = axes_flat[idx]
        sub = df[df["condition"] == condition].copy()

        sub["neg_log10_pval"] = sub["pvalue_adj"].apply(
            lambda v: -np.log10(v) if v > 0 else -np.log10(1e-300)
        )

        # Color: up = red, down = blue, ns = gray
        colors = np.where(
            (sub["pvalue_adj"] < pval_threshold) & (sub["log2fc"] >= logfc_threshold),
            "#d73027",   # up-regulated
            np.where(
                (sub["pvalue_adj"] < pval_threshold) & (sub["log2fc"] <= -logfc_threshold),
                "#4575b4",   # down-regulated
                "#aaaaaa",   # not significant
            ),
        )

        ax.scatter(sub["log2fc"], sub["neg_log10_pval"],
                   c=colors, s=3, alpha=0.6, linewidths=0, rasterized=True)

        ax.axhline(neg_log_thresh, color="gray", linestyle="--", linewidth=0.8)
        ax.axvline( logfc_threshold, color="gray", linestyle="--", linewidth=0.8)
        ax.axvline(-logfc_threshold, color="gray", linestyle="--", linewidth=0.8)

        # Annotate top up-regulated genes
        top = (sub[sub["is_significant"] & (sub["log2fc"] > 0)]
               .nsmallest(5, "pvalue_adj"))
        for _, row in top.iterrows():
            ax.annotate(row["gene"],
                        (row["log2fc"], row["neg_log10_pval"]),
                        fontsize=5, ha="left", va="bottom",
                        xytext=(2, 2), textcoords="offset points")

        n_up   = ((sub["pvalue_adj"] < pval_threshold) & (sub["log2fc"] >=  logfc_threshold)).sum()
        n_down = ((sub["pvalue_adj"] < pval_threshold) & (sub["log2fc"] <= -logfc_threshold)).sum()
        ax.set_title(f"{condition}\n↑{n_up}  ↓{n_down}", fontsize=8)
        ax.set_xlabel("log₂FC", fontsize=7)
        ax.set_ylabel("-log₁₀(adj. p-value)", fontsize=7)
        ax.tick_params(labelsize=6)

    # Hide unused panels
    for idx in range(n, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    legend_handles = [
        mpatches.Patch(color="#d73027", label=f"Up  (log2FC≥{logfc_threshold}, FDR<{pval_threshold})"),
        mpatches.Patch(color="#4575b4", label=f"Down (log2FC≤-{logfc_threshold}, FDR<{pval_threshold})"),
        mpatches.Patch(color="#aaaaaa", label="Not significant"),
    ]
    fig.legend(handles=legend_handles, loc="lower right",
               bbox_to_anchor=(1.0, 0.0), fontsize=7)

    fig.suptitle(f"DEG — {column}", fontsize=11, fontweight="bold")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Per-dataset pipeline
# ---------------------------------------------------------------------------

def _run_dataset(
    spec: DatasetSpec,
    out_dir: Path,
    logfc_threshold: float,
    pval_threshold: float,
) -> list[pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    adata    = _load_and_preprocess(spec)
    all_dfs  = []

    for column, conditions in spec.column_conditions.items():
        print(f"\n  Column: {column}  ({len(conditions)} conditions)")
        df = _run_column_deg(
            adata, spec.label, column, conditions,
            out_dir, logfc_threshold, pval_threshold,
        )
        if df is not None:
            all_dfs.append(df)

    return all_dfs


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run_deg_analysis(config_path: str) -> None:
    cfg, specs = _load_config(config_path)

    dataset_name    = cfg["dataset_name"]
    out_root        = Path(cfg["output_folder"])
    date_str        = datetime.datetime.now().strftime("%Y%m%d")
    run_dir         = out_root / f"{dataset_name}_{date_str}"
    run_dir.mkdir(parents=True, exist_ok=True)

    meta_dir = run_dir / "metadata"
    meta_dir.mkdir(parents=True, exist_ok=True)

    logfc_threshold = float(cfg.get("logfc_threshold", 1.0))
    pval_threshold  = float(cfg.get("pval_threshold", 0.05))

    with _log_to_file(meta_dir / "pipeline_output.log"):
        print("=" * 62)
        print(f"  DEG Analysis — {dataset_name}")
        print("=" * 62)
        print(f"  Output:          {run_dir}")
        print(f"  log2FC threshold: {logfc_threshold}")
        print(f"  p-adj threshold:  {pval_threshold}")
        print(f"  Datasets:")
        for s in specs:
            print(f"    · {s.label:<16} {s.h5ad_path.name}"
                  + (f"  [chemistry={s.chemistry}]" if s.chemistry else ""))
        print("=" * 62)

        all_results = []

        for spec in specs:
            print(f"\n{'─'*50}")
            print(f"  Dataset: {spec.label}")
            print(f"{'─'*50}")
            dfs = _run_dataset(spec, run_dir / spec.label, logfc_threshold, pval_threshold)
            all_results.extend(dfs)

        if not all_results:
            print("\n⚠️  No results produced.")
            return

        # Summary: only significant DEGs
        full     = pd.concat(all_results, ignore_index=True)
        sig      = full[full["is_significant"]].copy()
        sig_path = run_dir / "summary_significant.csv"
        sig.to_csv(sig_path, index=False)

        n_sig = len(sig)
        n_conditions = full[["dataset", "column", "condition"]].drop_duplicates().shape[0]
        print(f"\n📄 Significant DEGs saved: {sig_path}")
        print(f"🎉 Done! {n_sig:,} significant DEGs across {n_conditions} conditions.")
        print(f"   Results: {run_dir}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Wilcoxon DEG analysis (one condition vs rest) for each column."
    )
    parser.add_argument("config", help="YAML config file.")
    args = parser.parse_args()
    run_deg_analysis(args.config)
