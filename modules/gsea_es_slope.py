"""
gsea_es_slope.py

Compute GSEA-based enrichment metrics for a risk gene list ranked by GES scores.

For each condition the genes are ranked by GES score (descending) and
gseapy.prerank is applied with the user's gene set. Five metrics are computed:

  ES          enrichment score (maximum of the running enrichment curve)
  NES         normalized ES (permutation-based, from gseapy)
  pval_nes    nominal p-value for NES
  fdr_nes     FDR q-value for NES (BH-corrected across conditions per column)
  slope       ES / (peak_rank / n_total) — steepness of the enrichment curve
              up to the peak. High slope = gene set genes cluster at the very
              top of the ranked list; low slope = enrichment is spread out.
  pval_slope  empirical p-value for slope (permutation test)
  fdr_slope   FDR q-value for slope (BH-corrected across conditions per column)
  peak_rank   rank position (1-based) where ES is achieved
  n_hit       number of gene-set genes found in the ranked list

Outputs (in results/gsea_nes/{dataset_name}_{gene_list_stem}_{YYYYMMDD}/):
  data/
    {dataset}_{column}_results.csv
    all_results.csv
  figures/
    {dataset}_{column}_{condition}_gsea.png   (enrichment curve + slope line)
    {dataset}_{column}_nes_barplot.png
    {dataset}_{column}_slope_barplot.png
  metadata/pipeline_output.log

Usage:
  python modules/gsea_nes_module.py config.yaml gene_list.csv
"""

import argparse
import contextlib
import datetime
import shutil
import sys
from pathlib import Path
from typing import NamedTuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from statsmodels.stats.multitest import multipletests

try:
    import gseapy as gp
except ImportError:
    print("Error: gseapy not installed.  Run: pip install gseapy")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Data structure
# ---------------------------------------------------------------------------

class DatasetSpec(NamedTuple):
    label: str
    ges_folder: Path
    column_conditions: dict


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def _log_to_file(log_path: Path):
    class _Tee:
        def __init__(self, *streams): self._streams = streams
        def write(self, s):
            for st in self._streams: st.write(s)
        def flush(self):
            for st in self._streams: st.flush()
        @property
        def encoding(self): return getattr(self._streams[0], "encoding", "utf-8")
    orig = sys.stdout
    with open(log_path, "w", encoding="utf-8") as fh:
        sys.stdout = _Tee(orig, fh)
        try:
            yield
        finally:
            sys.stdout = orig


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

def load_config(config_path: str) -> tuple[dict, list[DatasetSpec]]:
    config_path = Path(config_path).resolve()
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    root = Path(cfg["ndd_gene_modules_folder_root"]).resolve()

    def resolve(p):
        p = Path(str(p))
        return p if p.is_absolute() else (root / p).resolve()

    cfg["output_folder"] = resolve(cfg.get("output_folder", "results/gsea_nes"))
    cfg["_config_path"]  = config_path

    if "datasets" in cfg:
        specs = []
        for label, ds in cfg["datasets"].items():
            specs.append(DatasetSpec(
                label=label,
                ges_folder=resolve(ds["ges_results_folder"]),
                column_conditions=ds["column_conditions"],
            ))
    else:
        specs = [DatasetSpec(
            label=cfg["dataset_name"],
            ges_folder=resolve(cfg["ges_results_folder"]),
            column_conditions=cfg["column_conditions"],
        )]

    return cfg, specs


# ---------------------------------------------------------------------------
# Gene set loader
# ---------------------------------------------------------------------------

def load_gene_set(gene_list_path: str) -> set[str]:
    """
    Load a gene list CSV.  Accepts:
      - CSV with a 'gene' column (possibly also 'group')
      - Single-column plain-text / headerless CSV (gene names, one per row)
    """
    df = pd.read_csv(gene_list_path)
    if "gene" in df.columns:
        genes = df["gene"].dropna().astype(str).tolist()
    else:
        genes = df.iloc[:, 0].dropna().astype(str).tolist()
    gene_set = {g.strip() for g in genes if g.strip()}
    print(f"Gene set loaded: {len(gene_set)} genes from {Path(gene_list_path).name}")
    return gene_set


# ---------------------------------------------------------------------------
# Running enrichment score (weighted GSEA, p=1)
# ---------------------------------------------------------------------------

def _running_es(ranked_genes: np.ndarray,
                scores: np.ndarray,
                gene_set_mask: np.ndarray) -> tuple[np.ndarray, float, int]:
    """
    Compute the weighted running enrichment score.

    Parameters
    ----------
    ranked_genes  : gene names in rank order (highest score first)
    scores        : corresponding GES scores
    gene_set_mask : boolean mask, True where gene is in the gene set

    Returns
    -------
    res      : running enrichment score array (length n)
    es       : enrichment score (value at peak)
    peak_idx : index (0-based) of the peak
    """
    n = len(scores)
    k = int(gene_set_mask.sum())
    if k == 0:
        return np.zeros(n), 0.0, 0

    n_r = float(np.sum(np.abs(scores[gene_set_mask])))
    if n_r == 0.0:
        return np.zeros(n), 0.0, 0

    miss_step = 1.0 / (n - k) if n > k else 0.0
    step = np.where(gene_set_mask,
                    np.abs(scores) / n_r,
                    -miss_step)
    res      = np.cumsum(step)
    peak_idx = int(np.argmax(np.abs(res)))
    es       = float(res[peak_idx])
    return res, es, peak_idx


# ---------------------------------------------------------------------------
# Slope permutation test
# ---------------------------------------------------------------------------

def _slope_permutation_pvalue(
    scores: np.ndarray,
    gene_set_mask: np.ndarray,
    obs_slope: float,
    n_perms: int,
    seed: int = 42,
) -> tuple[float, np.ndarray]:
    """
    Empirical p-value for the observed slope via permutation.

    Null: gene-set labels are randomly assigned to the ranked list
    (same gene-set size, random positions).

    Returns (p_slope, perm_slopes_array).
    """
    n = len(scores)
    k = int(gene_set_mask.sum())
    rng = np.random.default_rng(seed)

    perm_slopes = np.zeros(n_perms)
    for i in range(n_perms):
        perm_mask = np.zeros(n, dtype=bool)
        perm_mask[rng.choice(n, size=k, replace=False)] = True
        _, es_p, peak_p = _running_es(None, scores, perm_mask)
        peak_frac = (peak_p + 1) / n
        perm_slopes[i] = es_p / peak_frac if peak_frac > 0 else 0.0

    # One-sided: fraction of permutation |slopes| >= |observed slope|
    p_val = float(np.mean(np.abs(perm_slopes) >= abs(obs_slope)))
    return p_val, perm_slopes


# ---------------------------------------------------------------------------
# Single-condition analysis
# ---------------------------------------------------------------------------

def _analyze_condition(
    ges_path: Path,
    gene_set: set[str],
    gene_set_name: str,
    ges_score_threshold: float | None,
    permutation_num: int,
    min_size: int,
) -> dict | None:
    """
    Run gseapy.prerank + slope computation for one GES file.
    Returns a dict of metrics, or None if analysis cannot proceed.
    """
    df = pd.read_csv(ges_path)
    if "gene" not in df.columns or "ges_score" not in df.columns:
        print(f"  Warning: unexpected columns in {ges_path.name} — skipping.")
        return None

    if ges_score_threshold is not None:
        df = df[df["ges_score"] > ges_score_threshold]

    df = df.sort_values("ges_score", ascending=False).drop_duplicates("gene")
    ranking = df.set_index("gene")["ges_score"]

    n_hit = sum(1 for g in ranking.index if g in gene_set)
    if n_hit < min_size:
        print(f"  Warning: only {n_hit} gene-set genes in ranked list "
              f"(min_size={min_size}) — skipping.")
        return None

    # ── gseapy prerank ──────────────────────────────────────────────────────
    try:
        pre_res = gp.prerank(
            rnk=ranking,
            gene_sets={gene_set_name: sorted(gene_set)},
            min_size=min_size,
            max_size=len(ranking),
            permutation_num=permutation_num,
            ascending=False,
            outdir=None,
            no_plot=True,
            seed=42,
            verbose=False,
        )
    except Exception as exc:
        print(f"  gseapy error: {exc} — skipping.")
        return None

    if pre_res is None or pre_res.res2d is None or pre_res.res2d.empty:
        print("  gseapy returned empty results — skipping.")
        return None

    row = pre_res.res2d.iloc[0]
    es      = float(row.get("ES",        row.get("es",       np.nan)))
    nes     = float(row.get("NES",       row.get("nes",      np.nan)))
    pval    = float(row.get("NOM p-val", row.get("pval",     np.nan)))
    fdr_gp  = float(row.get("FDR q-val", row.get("fdr",     np.nan)))

    # ── RES curve from gseapy ───────────────────────────────────────────────
    term_key = gene_set_name
    if term_key not in pre_res.results:
        term_key = list(pre_res.results.keys())[0]

    term_data = pre_res.results[term_key]
    res_curve = np.array(term_data["RES"])
    hit_idx   = np.array(term_data.get("hits", []), dtype=int)

    peak_idx  = int(np.argmax(np.abs(res_curve)))
    n_total   = len(res_curve)
    peak_frac = (peak_idx + 1) / n_total
    slope     = es / peak_frac if peak_frac > 0 else 0.0

    # ── Slope permutation p-value ───────────────────────────────────────────
    ranked_genes_arr = ranking.index.to_numpy()
    scores_arr       = ranking.values.astype(float)
    gene_set_mask    = np.array([g in gene_set for g in ranked_genes_arr], dtype=bool)

    pval_slope, perm_slopes = _slope_permutation_pvalue(
        scores_arr, gene_set_mask, slope, permutation_num
    )

    return {
        "ES":          es,
        "NES":         nes,
        "pval_nes":    pval,
        "fdr_nes_gp":  fdr_gp,   # gseapy's own FDR (kept for reference)
        "slope":       slope,
        "pval_slope":  pval_slope,
        "peak_rank":   peak_idx + 1,
        "peak_frac":   peak_frac,
        "n_total":     n_total,
        "n_hit":       n_hit,
        # internal — used for plotting
        "_res_curve":  res_curve,
        "_hit_idx":    hit_idx,
        "_scores":     scores_arr,
        "_perm_slopes": perm_slopes,
    }


# ---------------------------------------------------------------------------
# Enrichment curve plot (per condition)
# ---------------------------------------------------------------------------

def _plot_gsea(
    result: dict,
    condition_label: str,
    x_label: str,
    fig_path: Path,
) -> None:
    """
    Three-panel GSEA plot with slope line overlay on the enrichment curve.

    Panel 1 — running enrichment score + slope line up to peak
    Panel 2 — gene hit barcode
    Panel 3 — GES ranking metric
    """
    res_curve  = result["_res_curve"]
    hit_idx    = result["_hit_idx"]
    scores     = result["_scores"]
    es         = result["ES"]
    nes        = result["NES"]
    pval_nes   = result["pval_nes"]
    slope      = result["slope"]
    pval_slope = result["pval_slope"]
    peak_rank  = result["peak_rank"]
    n_total    = result["n_total"]

    x_all  = np.arange(n_total)
    peak_i = peak_rank - 1

    fig, axes = plt.subplots(
        3, 1,
        figsize=(10, 8),
        gridspec_kw={"height_ratios": [4, 0.8, 1.2]},
        sharex=True,
    )

    # ── Panel 1: running ES ─────────────────────────────────────────────────
    ax = axes[0]
    ax.plot(x_all, res_curve, color="forestgreen", linewidth=1.8, label="RES")
    ax.axhline(0, color="gray", linewidth=0.6, linestyle="--")
    ax.axvline(peak_i, color="orange", linewidth=1.0, linestyle=":", alpha=0.8,
               label=f"Peak  (rank {peak_rank:,})")

    # Slope line: from (0, 0) to (peak_i, es), then dotted extension
    ax.plot([0, peak_i], [0, es], color="red", linewidth=1.2,
            linestyle="--", alpha=0.75, label=f"Slope = {slope:.3f}")

    ax.set_ylabel("Enrichment Score", fontsize=10)
    ax.set_title(
        f"{condition_label}\n"
        f"ES={es:.3f}   NES={nes:.3f}   p(NES)={pval_nes:.3e}   "
        f"slope={slope:.3f}   p(slope)={pval_slope:.3e}",
        fontsize=9,
    )
    ax.legend(fontsize=7, loc="upper right", frameon=False)
    ax.set_frame_on(False)

    # ── Panel 2: gene hit barcode ───────────────────────────────────────────
    ax2 = axes[1]
    ax2.axhline(0.5, color="#dddddd", linewidth=0.5)
    if len(hit_idx):
        ax2.scatter(hit_idx, np.ones(len(hit_idx)) * 0.5,
                    color="black", s=6, marker="|", alpha=0.7, linewidths=0.6)
    ax2.set_yticks([])
    ax2.set_ylabel("Hits", fontsize=8)
    ax2.set_frame_on(False)

    # ── Panel 3: ranking metric ─────────────────────────────────────────────
    ax3 = axes[2]
    pos_mask = scores >= 0
    ax3.fill_between(x_all[pos_mask],  scores[pos_mask], 0, color="#d62728", alpha=0.6)
    ax3.fill_between(x_all[~pos_mask], scores[~pos_mask], 0, color="#1f77b4", alpha=0.6)
    ax3.axhline(0, color="gray", linewidth=0.4)
    ax3.set_ylabel("GES score", fontsize=8)
    ax3.set_xlabel(x_label, fontsize=9)
    ax3.set_frame_on(False)

    plt.tight_layout()
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Summary bar plots (per column)
# ---------------------------------------------------------------------------

def _barplot(
    results_df: pd.DataFrame,
    value_col: str,
    pval_col: str,
    ylabel: str,
    title: str,
    fig_path: Path,
    sig_threshold: float = 0.05,
) -> None:
    df = results_df.sort_values(value_col, ascending=False)
    colors = ["#d62728" if p < sig_threshold else "#aaaaaa"
              for p in df[pval_col]]

    fig, ax = plt.subplots(figsize=(max(6, len(df) * 0.6 + 2), 5))
    bars = ax.bar(df["condition"].astype(str), df[value_col], color=colors, edgecolor="white")
    ax.axhline(0, color="gray", linewidth=0.6)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=10)
    ax.set_xticklabels(df["condition"].astype(str), rotation=45, ha="right", fontsize=8)

    from matplotlib.patches import Patch
    ax.legend(
        handles=[Patch(color="#d62728", label=f"p < {sig_threshold}"),
                 Patch(color="#aaaaaa", label=f"p ≥ {sig_threshold}")],
        fontsize=7, frameon=False,
    )
    ax.set_frame_on(False)
    plt.tight_layout()
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Per-dataset pipeline
# ---------------------------------------------------------------------------

def _process_dataset(
    spec: DatasetSpec,
    gene_set: set[str],
    gene_set_name: str,
    params: dict,
    out_dir: Path,
) -> list[dict]:
    """Run analysis for all columns/conditions in one dataset. Returns list of result rows."""
    ges_data_dir = spec.ges_folder / "data"
    fig_dir      = out_dir / "figures"
    data_dir     = out_dir / "data"
    fig_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict] = []

    for column, conditions in spec.column_conditions.items():
        col_rows: list[dict] = []

        for condition in conditions:
            ges_path = ges_data_dir / f"ges_spec_{column}_{condition}.csv"
            if not ges_path.exists():
                print(f"  Missing: {ges_path.name} — skipping.")
                continue

            print(f"  {column} / {condition} ...")
            metrics = _analyze_condition(
                ges_path=ges_path,
                gene_set=gene_set,
                gene_set_name=gene_set_name,
                ges_score_threshold=params["ges_score_threshold"],
                permutation_num=params["permutation_num"],
                min_size=params["min_size"],
            )
            if metrics is None:
                continue

            row = {
                "dataset":    spec.label,
                "column":     column,
                "condition":  str(condition),
                "ES":         metrics["ES"],
                "NES":        metrics["NES"],
                "pval_nes":   metrics["pval_nes"],
                "slope":      metrics["slope"],
                "pval_slope": metrics["pval_slope"],
                "peak_rank":  metrics["peak_rank"],
                "peak_frac":  metrics["peak_frac"],
                "n_total":    metrics["n_total"],
                "n_hit":      metrics["n_hit"],
            }
            col_rows.append(row)

            # Enrichment curve plot
            safe_cond = str(condition).replace(" ", "_").replace("/", "-")
            fig_path  = fig_dir / f"{spec.label}_{column}_{safe_cond}_gsea.png"
            _plot_gsea(
                result=metrics,
                condition_label=f"{spec.label} | {column} | {condition}",
                x_label=f"Rank in GES list (n={metrics['n_total']:,})",
                fig_path=fig_path,
            )
            print(f"    Saved: {fig_path.name}")

        if not col_rows:
            continue

        # BH correction within column
        col_df = pd.DataFrame(col_rows)
        for pval_col, fdr_col in [("pval_nes", "fdr_nes"), ("pval_slope", "fdr_slope")]:
            pvals  = col_df[pval_col].fillna(1.0).values
            _, fdr = multipletests(pvals, method="fdr_bh")[:2]
            col_df[fdr_col] = fdr

        # Save per-column CSV
        col_csv = data_dir / f"{spec.label}_{column}_results.csv"
        col_df.to_csv(col_csv, index=False)
        print(f"  Saved: {col_csv.name}")

        # Summary bar plots
        safe_col = column.replace(" ", "_")
        _barplot(col_df, "NES",   "fdr_nes",   "NES",
                 f"NES — {spec.label} / {column}",
                 fig_dir / f"{spec.label}_{safe_col}_nes_barplot.png")
        _barplot(col_df, "slope", "fdr_slope", "Slope",
                 f"Slope — {spec.label} / {column}",
                 fig_dir / f"{spec.label}_{safe_col}_slope_barplot.png")

        all_rows.extend(col_df.to_dict("records"))

    return all_rows


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def run_gsea_nes(config_path: str, gene_list_path: str) -> None:
    cfg, specs = load_config(config_path)

    dataset_name  = cfg["dataset_name"]
    gene_list_stem = Path(gene_list_path).stem
    date_str      = datetime.datetime.now().strftime("%Y%m%d")
    out_root      = cfg["output_folder"]
    run_dir       = out_root / f"{dataset_name}_{gene_list_stem}_{date_str}"
    run_dir.mkdir(parents=True, exist_ok=True)

    shutil.copy2(cfg["_config_path"], run_dir / Path(config_path).name)

    meta_dir = run_dir / "metadata"
    meta_dir.mkdir(exist_ok=True)

    with _log_to_file(meta_dir / "pipeline_output.log"):
        params = {
            "ges_score_threshold": cfg.get("ges_score_threshold"),
            "permutation_num":     int(cfg.get("permutation_num", 1000)),
            "min_size":            int(cfg.get("min_gene_set_size", 5)),
        }

        print(f"\n{'='*62}")
        print(f"  GSEA NES module — {dataset_name}")
        print(f"  Gene list:    {gene_list_path}")
        print(f"  Permutations: {params['permutation_num']}")
        print(f"  GES threshold:{params['ges_score_threshold']}")
        print(f"  Output:       {run_dir}")
        print(f"{'='*62}\n")

        gene_set      = load_gene_set(gene_list_path)
        gene_set_name = gene_list_stem

        all_rows: list[dict] = []

        for spec in specs:
            print(f"\n--- Dataset: {spec.label} ---")
            spec_dir = run_dir / spec.label if len(specs) > 1 else run_dir
            rows = _process_dataset(spec, gene_set, gene_set_name, params, spec_dir)
            all_rows.extend(rows)

        if all_rows:
            combined = pd.DataFrame(all_rows)
            combined_path = run_dir / "data" / "all_results.csv"
            combined_path.parent.mkdir(exist_ok=True)
            combined.to_csv(combined_path, index=False)
            print(f"\nCombined results: {combined_path}")

        print(f"\nDone → {run_dir}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="GSEA NES + slope analysis for a gene list ranked by GES scores."
    )
    parser.add_argument("config",     help="YAML config file.")
    parser.add_argument("gene_list",  help="Gene list CSV.")
    args = parser.parse_args()
    run_gsea_nes(args.config, args.gene_list)
