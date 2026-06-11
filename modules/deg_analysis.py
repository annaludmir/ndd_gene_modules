"""
deg_analysis.py

Two-step analysis:
  1. Compute DEGs per condition (Wilcoxon rank-sum, condition vs rest).
  2. Test whether an input gene list (risk genes) is enriched in each
     condition's upregulated DEG set — hypergeometric test.

Approach follows:
  "Risk genes were subjected to enrichment tests with cluster-specific DEGs
   with a threshold of FDR < 0.05 and log2-fold change > 0.2."
   (doi:10.1038/s12276-024-01328-6)

Preprocessing applied identically to the GES pipeline:
  sc.pp.filter_cells(min_genes=200)
  sc.pp.filter_genes(min_cells=3)
  sc.pp.normalize_total()
  sc.pp.log1p()

Independent of GES scores, tau filtering, and GSEA.

Usage
-----
  python modules/deg_analysis.py <config.yaml> <gene_list.csv>

  Gene list: CSV with a 'gene' column, or plain text (one gene per line).

Output
------
  results/deg_analysis/{dataset_name}_{gene_list_name}_{YYYYMMDD}/
    metadata/pipeline_output.log
    {dataset}/
      {column}_deg_results.csv         all genes × all conditions
      {column}_volcano.png             volcano grid (one panel per condition)
      {column}_enrichment_results.csv  hypergeometric enrichment per condition
      {column}_enrichment_barplot.png  -log10(FDR) bar chart per condition
    summary_enrichment.csv             enrichment across all datasets/columns
"""

import argparse
import contextlib
import datetime
import re
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, NamedTuple

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.stats.multitest as smm
import yaml
from scipy.stats import hypergeom

# Optional: mygene for Ensembl → symbol conversion (same as GES pipeline)
try:
    from mygene import MyGeneInfo
    _HAS_MYGENE = True
    _mg = MyGeneInfo()
except Exception:
    _HAS_MYGENE = False
    _mg = None


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
    chemistry: str | None


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
        specs = [
            DatasetSpec(
                label=label,
                h5ad_path=resolve(ds["h5ad_path"]),
                column_conditions=ds["column_conditions"],
                chemistry=ds.get("chemistry"),
            )
            for label, ds in cfg["datasets"].items()
        ]
        return cfg, specs

    specs = [DatasetSpec(
        label=cfg["dataset_name"],
        h5ad_path=resolve(cfg["h5ad_path"]),
        column_conditions=cfg["column_conditions"],
        chemistry=cfg.get("chemistry"),
    )]
    return cfg, specs


# ---------------------------------------------------------------------------
# Gene list loader
# ---------------------------------------------------------------------------

def _load_gene_list(path) -> set:
    path = Path(path)
    if path.suffix.lower() == ".csv":
        df  = pd.read_csv(path)
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
# Gene name helpers — identical logic to GES pipeline
# (all_layers stores Ensembl IDs in var_names; symbols are in a var column)
# ---------------------------------------------------------------------------

def _looks_like_ensembl(x: str) -> bool:
    return bool(re.match(r"^ENS[A-Z]*G\d+(\.\d+)?$", str(x)))


def _make_unique_symbols(symbols: List[str]) -> List[str]:
    seen: Dict[str, int] = {}
    out: List[str] = []
    for s in symbols:
        if s not in seen:
            seen[s] = 0
            out.append(s)
        else:
            seen[s] += 1
            out.append(f"{s}_dup{seen[s]}")
    return out


def _convert_ensembl_to_symbol(gene_list: List[str], species: str = "human") -> List[str]:
    if not _HAS_MYGENE:
        return gene_list
    ensembls = [g for g in gene_list if _looks_like_ensembl(g)]
    if not ensembls:
        return gene_list
    res = _mg.querymany(
        ensembls, scopes="ensembl.gene", fields="symbol",
        species=species, as_dataframe=True,
    )
    if "symbol" not in res.columns:
        return gene_list
    if "query" in res.columns:
        mapping = dict(zip(res["query"], res["symbol"].fillna(res["query"])))
    else:
        mapping = res["symbol"].fillna(res.index.to_series()).to_dict()
    converted = [mapping.get(g, g) for g in gene_list]
    return _make_unique_symbols([str(x) for x in converted])


def _get_gene_names(adata: sc.AnnData, species: str = "human") -> List[str]:
    """
    Extract gene symbols from AnnData — same logic as GES pipeline.
    Checks var columns first (gene_symbol, Gene, …), then var_names.
    Converts Ensembl IDs to symbols when detected.
    """
    for col in ("gene_symbol", "gene_symbols", "Gene", "gene", "name", "names"):
        if col in adata.var.columns:
            genes = adata.var[col].astype(str).tolist()
            break
    else:
        genes = adata.var_names.astype(str).tolist()

    if any(_looks_like_ensembl(g) for g in genes):
        print("  Detected Ensembl IDs — converting to gene symbols (mygene)...")
        genes = _convert_ensembl_to_symbol(genes, species=species)

    return genes


# ---------------------------------------------------------------------------
# Preprocessing — identical to GES pipeline
# ---------------------------------------------------------------------------

def _load_and_preprocess(spec: DatasetSpec) -> sc.AnnData:
    print(f"  Loading {spec.h5ad_path.name} ...")
    adata = sc.read_h5ad(spec.h5ad_path)
    print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    if spec.chemistry and "Chemistry" in adata.obs.columns:
        before = adata.n_obs
        adata  = adata[adata.obs["Chemistry"] == spec.chemistry].copy()
        print(f"  Chemistry={spec.chemistry}: {before:,} → {adata.n_obs:,} cells")

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Resolve gene names before normalization so var_names are symbols
    # (all_layers uses Ensembl IDs in var_names; symbols live in a var column)
    gene_names = _get_gene_names(adata)
    adata.var_names = pd.Index(gene_names)
    print(f"  Gene names resolved (e.g. {gene_names[:3]})")

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    print(f"  After QC + normalization: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    return adata


# ---------------------------------------------------------------------------
# Step 1 — DEG computation (Wilcoxon, condition vs rest)
# ---------------------------------------------------------------------------

def _compute_degs(
    adata: sc.AnnData,
    dataset_label: str,
    column: str,
    conditions: list,
) -> pd.DataFrame | None:
    """Run Wilcoxon rank-sum test for each condition vs all others."""
    if column not in adata.obs.columns:
        print(f"    ⚠️  Column '{column}' not in obs — skipping.")
        return None

    obs_vals  = adata.obs[column].astype(str)
    available = [str(c) for c in conditions if str(c) in obs_vals.unique()]
    missing   = [str(c) for c in conditions if str(c) not in obs_vals.unique()]
    if missing:
        print(f"    ⚠️  Conditions not in obs: {missing}")
    if len(available) < 2:
        print(f"    ⚠️  Fewer than 2 conditions available in '{column}' — skipping.")
        return None

    adata_sub = adata[obs_vals.isin(available)].copy()
    print(f"    {adata_sub.n_obs:,} cells  |  {len(available)} conditions")

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

    dfs = []
    for condition in available:
        df = sc.get.rank_genes_groups_df(adata_sub, group=condition)
        df = df.rename(columns={
            "names":            "gene",
            "scores":           "wilcoxon_score",
            "logfoldchanges":   "log2fc",
            "pvals":            "pvalue",
            "pvals_adj":        "pvalue_adj",
            "pct_nz_group":     "pct_expressing_condition",
            "pct_nz_reference": "pct_expressing_rest",
        })
        df.insert(0, "dataset",   dataset_label)
        df.insert(1, "column",    column)
        df.insert(2, "condition", condition)
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)


# ---------------------------------------------------------------------------
# Step 2 — Hypergeometric enrichment of risk genes in condition DEG sets
# ---------------------------------------------------------------------------

def _fmt_pvals(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in ("pvalue", "FDR_qval_BH"):
        if col in df.columns:
            df[col] = df[col].apply(
                lambda v: f"{float(v):.6e}" if pd.notna(v) else v
            )
    return df


def _enrich_gene_list(
    deg_df: pd.DataFrame,
    column: str,
    risk_genes: set,
    n_background: int,
    deg_logfc_threshold: float,
    deg_pval_threshold: float,
) -> pd.DataFrame | None:
    """
    For each condition, define the DEG set (upregulated: pval_adj < threshold
    and log2fc > logfc_threshold), then test hypergeometric enrichment of
    risk_genes in that DEG set.

    Background N = total genes in the analysis (adata.n_vars after QC).
    """
    risk_in_bg = risk_genes & set(deg_df["gene"].unique())
    N = n_background
    K = len(risk_in_bg)

    print(f"    Background={N:,} | risk genes in background={K}")

    if K == 0:
        print(f"    ⚠️  No risk genes overlap with background — skipping enrichment.")
        return None

    rows = []
    for condition, grp in deg_df.groupby("condition"):
        # Condition-specific upregulated DEGs
        deg_set = set(
            grp.loc[
                (grp["pvalue_adj"] < deg_pval_threshold) &
                (grp["log2fc"]     > deg_logfc_threshold),
                "gene",
            ]
        )
        n = len(deg_set)
        if n == 0:
            print(f"    ⚠️  No DEGs for {condition} at thresholds — enrichment p=1")

        overlap = risk_in_bg & deg_set
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
            "column":                  column,
            "condition":               str(condition),
            "n_background":            N,
            "n_risk_in_background":    K,
            "n_condition_degs":        n,
            "n_overlap":               k,
            "expected_overlap":        round(expected, 4),
            "fold_enrichment":         fold_enrichment,
            "pvalue":                  pval,
            "FDR_qval_BH":             np.nan,
            "is_significant":          False,
            "overlap_genes":           ";".join(sorted(overlap)),
        })

    if not rows:
        return None

    enr = pd.DataFrame(rows)
    _, corrected = smm.multipletests(enr["pvalue"].astype(float), method="fdr_bh")[:2]
    enr["FDR_qval_BH"]  = corrected
    enr["is_significant"] = corrected < 0.05
    return enr


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def _plot_volcano_grid(
    df: pd.DataFrame,
    column: str,
    deg_logfc_threshold: float,
    deg_pval_threshold: float,
    out_path: Path,
) -> None:
    conditions = df["condition"].unique().tolist()
    ncols      = min(4, len(conditions))
    nrows      = (len(conditions) + ncols - 1) // ncols
    fig, axes  = plt.subplots(nrows, ncols,
                               figsize=(ncols * 4, nrows * 3.8), squeeze=False)
    axes_flat  = axes.flatten()

    for idx, condition in enumerate(conditions):
        ax  = axes_flat[idx]
        sub = df[df["condition"] == condition].copy()
        sub["neg_log10_pval"] = sub["pvalue_adj"].apply(
            lambda v: -np.log10(v) if v > 0 else -np.log10(1e-300)
        )
        colors = np.where(
            (sub["pvalue_adj"] < deg_pval_threshold) & (sub["log2fc"] > deg_logfc_threshold),
            "#d73027",
            np.where(
                (sub["pvalue_adj"] < deg_pval_threshold) & (sub["log2fc"] < -deg_logfc_threshold),
                "#4575b4",
                "#aaaaaa",
            ),
        )
        ax.scatter(sub["log2fc"], sub["neg_log10_pval"],
                   c=colors, s=3, alpha=0.6, linewidths=0, rasterized=True)
        ax.axhline(-np.log10(deg_pval_threshold), color="gray", linestyle="--", linewidth=0.8)
        ax.axvline( deg_logfc_threshold, color="gray", linestyle="--", linewidth=0.8)
        ax.axvline(-deg_logfc_threshold, color="gray", linestyle="--", linewidth=0.8)

        top = (sub[(sub["pvalue_adj"] < deg_pval_threshold) & (sub["log2fc"] > deg_logfc_threshold)]
               .nsmallest(5, "pvalue_adj"))
        for _, row in top.iterrows():
            ax.annotate(row["gene"], (row["log2fc"], row["neg_log10_pval"]),
                        fontsize=5, ha="left", va="bottom",
                        xytext=(2, 2), textcoords="offset points")

        n_up   = ((sub["pvalue_adj"] < deg_pval_threshold) & (sub["log2fc"] >  deg_logfc_threshold)).sum()
        n_down = ((sub["pvalue_adj"] < deg_pval_threshold) & (sub["log2fc"] < -deg_logfc_threshold)).sum()
        ax.set_title(f"{condition}\n↑{n_up}  ↓{n_down}", fontsize=8)
        ax.set_xlabel("log₂FC", fontsize=7)
        ax.set_ylabel("-log₁₀(adj. p-value)", fontsize=7)
        ax.tick_params(labelsize=6)

    for idx in range(len(conditions), len(axes_flat)):
        axes_flat[idx].set_visible(False)

    fig.legend(
        handles=[
            mpatches.Patch(color="#d73027", label=f"Up  (log2FC>{deg_logfc_threshold}, FDR<{deg_pval_threshold})"),
            mpatches.Patch(color="#4575b4", label=f"Down (log2FC<-{deg_logfc_threshold}, FDR<{deg_pval_threshold})"),
            mpatches.Patch(color="#aaaaaa", label="Not significant"),
        ],
        loc="lower right", bbox_to_anchor=(1.0, 0.0), fontsize=7,
    )
    fig.suptitle(f"DEG — {column}", fontsize=11, fontweight="bold")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()


def _plot_enrichment_barplot(
    enr: pd.DataFrame,
    column: str,
    gene_list_name: str,
    out_path: Path,
) -> None:
    df = enr.copy()
    df["neg_log10_fdr"] = df["FDR_qval_BH"].apply(
        lambda v: -np.log10(float(v)) if float(v) > 0 else -np.log10(1e-300)
    )
    df = df.sort_values("neg_log10_fdr", ascending=True)

    fig, ax = plt.subplots(figsize=(8, max(3, len(df) * 0.5 + 1.2)))
    colors  = ["#d73027" if s else "#4575b4" for s in df["is_significant"]]
    bars    = ax.barh(df["condition"].astype(str), df["neg_log10_fdr"],
                      color=colors, edgecolor="white", height=0.65)

    ax.axvline(-np.log10(0.05), color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("-log₁₀(FDR)", fontsize=11)
    ax.set_title(
        f"Risk gene enrichment in DEGs — {column}\n"
        f"({gene_list_name}, n_risk={int(df['n_risk_in_background'].iloc[0])})",
        fontsize=11, fontweight="bold",
    )

    for bar, (_, row) in zip(bars, df.iterrows()):
        fe  = row["fold_enrichment"]
        fe_label = f"FE={float(fe):.2f}" if np.isfinite(float(fe)) else "FE=∞"
        ax.text(bar.get_width() + 0.05,
                bar.get_y() + bar.get_height() / 2,
                f"{fe_label}  (k={int(row['n_overlap'])}/{int(row['n_condition_degs'])})",
                va="center", ha="left", fontsize=8)

    ax.legend(
        handles=[
            mpatches.Patch(color="#d73027", label="FDR < 0.05"),
            mpatches.Patch(color="#4575b4", label="FDR ≥ 0.05"),
            plt.Line2D([0], [0], color="gray", linestyle="--", linewidth=1, label="FDR = 0.05"),
        ],
        loc="lower right", fontsize=8,
    )
    ax.set_xlim(left=0)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Per-dataset pipeline
# ---------------------------------------------------------------------------

def _run_dataset(
    spec: DatasetSpec,
    risk_genes: set,
    gene_list_name: str,
    out_dir: Path,
    deg_logfc_threshold: float,
    deg_pval_threshold: float,
) -> list[pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    adata    = _load_and_preprocess(spec)
    n_bg     = adata.n_vars          # background = all genes after QC
    enr_dfs  = []

    for column, conditions in spec.column_conditions.items():
        print(f"\n  Column: {column}  ({len(conditions)} conditions)")

        # Step 1: DEGs
        deg_df = _compute_degs(adata, spec.label, column, conditions)
        if deg_df is None:
            continue

        safe_col = column.replace(" ", "_")
        deg_csv  = out_dir / f"{safe_col}_deg_results.csv"
        deg_df.to_csv(deg_csv, index=False)
        print(f"    ✔ {deg_csv.name}")

        vol_path = out_dir / f"{safe_col}_volcano.png"
        _plot_volcano_grid(deg_df, column, deg_logfc_threshold, deg_pval_threshold, vol_path)
        print(f"    ✔ {vol_path.name}")

        # Step 2: Enrichment
        enr_df = _enrich_gene_list(
            deg_df, column, risk_genes, n_bg,
            deg_logfc_threshold, deg_pval_threshold,
        )
        if enr_df is None:
            continue

        enr_df.insert(0, "dataset", spec.label)
        enr_csv  = out_dir / f"{safe_col}_enrichment_results.csv"
        _fmt_pvals(enr_df).to_csv(enr_csv, index=False)
        print(f"    ✔ {enr_csv.name}")

        bar_path = out_dir / f"{safe_col}_enrichment_barplot.png"
        _plot_enrichment_barplot(enr_df, column, gene_list_name, bar_path)
        print(f"    ✔ {bar_path.name}")

        enr_dfs.append(enr_df)

    return enr_dfs


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run_deg_analysis(config_path: str, gene_list_path: str) -> None:
    cfg, specs = _load_config(config_path)

    dataset_name    = cfg["dataset_name"]
    out_root        = Path(cfg["output_folder"])
    date_str        = datetime.datetime.now().strftime("%Y%m%d")
    gene_list_name  = Path(gene_list_path).stem
    run_dir         = out_root / f"{dataset_name}_{gene_list_name}_{date_str}"
    run_dir.mkdir(parents=True, exist_ok=True)

    meta_dir = run_dir / "metadata"
    meta_dir.mkdir(parents=True, exist_ok=True)

    deg_logfc_threshold = float(cfg.get("deg_logfc_threshold", 0.2))
    deg_pval_threshold  = float(cfg.get("deg_pval_threshold",  0.05))

    with _log_to_file(meta_dir / "pipeline_output.log"):
        print("=" * 62)
        print(f"  DEG + Enrichment Analysis — {dataset_name}")
        print("=" * 62)
        print(f"  Risk gene list:    {gene_list_path}")
        print(f"  Output:            {run_dir}")
        print(f"  DEG thresholds:    log2FC > {deg_logfc_threshold}, FDR < {deg_pval_threshold}")
        print(f"  Datasets:")
        for s in specs:
            chem = f"  [chemistry={s.chemistry}]" if s.chemistry else ""
            print(f"    · {s.label:<16} {s.h5ad_path.name}{chem}")
        print("=" * 62)

        risk_genes = _load_gene_list(gene_list_path)
        print(f"\n🧬 Risk genes loaded: {len(risk_genes)}")

        all_enr = []
        for spec in specs:
            print(f"\n{'─'*50}")
            print(f"  Dataset: {spec.label}")
            print(f"{'─'*50}")
            enr_dfs = _run_dataset(
                spec, risk_genes, gene_list_name,
                run_dir / spec.label,
                deg_logfc_threshold, deg_pval_threshold,
            )
            all_enr.extend(enr_dfs)

        if not all_enr:
            print("\n⚠️  No enrichment results produced.")
            return

        summary = pd.concat(all_enr, ignore_index=True)
        summary_path = run_dir / "summary_enrichment.csv"
        _fmt_pvals(summary).to_csv(summary_path, index=False)

        n_sig = summary["is_significant"].sum()
        print(f"\n📄 Summary saved: {summary_path}")
        print(f"🎉 Done! {n_sig}/{len(summary)} conditions significantly enriched (FDR < 0.05).")
        print(f"   Results: {run_dir}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "DEG analysis (Wilcoxon, condition vs rest) followed by "
            "hypergeometric enrichment of a risk gene list in each condition's DEG set."
        )
    )
    parser.add_argument("config",    help="YAML config file.")
    parser.add_argument("gene_list", help="Risk gene list: CSV with 'gene' column or plain text.")
    args = parser.parse_args()
    run_deg_analysis(args.config, args.gene_list)
