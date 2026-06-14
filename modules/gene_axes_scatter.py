"""
gene_axes_scatter.py

2-D scatter plot where each gene is positioned along two user-defined 'axis scores'
derived from GES specificity scores.

Each axis score = mean GES(numerator conditions) - mean GES(denominator conditions).

Typical use (mirrors the referenced figure):
  X: Differentiation score = mean GES(Neuron, Neuroblast) - mean GES(Radial glia, NPC)
  Y: Cell cycle score      = mean GES(S, G2M, PostM) - GES(Non-cycling)

Genes can be coloured by:
  - Gene group (CSV with 'gene','group' columns, passed via --gene-list)
  - Leiden cluster (gene_clusters.csv from the gene_clustering pipeline, via --gene-clusters)
  - Default: all gray

Top-N genes per quadrant are annotated on the plot.

Outputs (in results/gene_axes_scatter/{dataset_name}_{YYYYMMDD}/):
  data/gene_scores.csv              — gene × score table
  figures/gene_axes_scatter.png     — main scatter plot

Usage:
  python modules/gene_axes_scatter.py config_files/gene_axes_scatter_v3_config.yaml
  python modules/gene_axes_scatter.py config_files/gene_axes_scatter_v3_config.yaml \\
        --gene-list gene_lists/my_genes.csv
  python modules/gene_axes_scatter.py config_files/gene_axes_scatter_v3_config.yaml \\
        --gene-clusters results/gene_clusters/.../cortex/data/gene_clusters.csv
"""

import argparse
import datetime
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
import yaml


_GROUP_PALETTE = (
    list(sns.color_palette("Set1",  9))
    + list(sns.color_palette("Set2",  8))
    + list(sns.color_palette("tab20", 20))
)


# ---------------------------------------------------------------------------
# Config loader
# ---------------------------------------------------------------------------

def load_config(config_path: str) -> dict:
    config_path = Path(config_path).resolve()
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    if "dataset_name" not in cfg:
        raise ValueError("Config missing required key: 'dataset_name'")

    root = Path(cfg["ndd_gene_modules_folder_root"]).resolve()

    def resolve(p):
        p = Path(str(p))
        return p if p.is_absolute() else (root / p).resolve()

    cfg["output_folder"] = resolve(cfg.get("output_folder", "results/gene_axes_scatter"))
    cfg["_config_path"] = config_path
    cfg["_root"] = root

    for axis_key in ("x_axis", "y_axis"):
        if axis_key not in cfg:
            raise ValueError(f"Config missing required key: '{axis_key}'")
        cfg[axis_key]["ges_results_folder"] = resolve(
            cfg[axis_key]["ges_results_folder"]
        )

    if cfg.get("gene_clusters_csv"):
        cfg["gene_clusters_csv"] = resolve(cfg["gene_clusters_csv"])

    if cfg.get("gene_list_path"):
        cfg["gene_list_path"] = resolve(cfg["gene_list_path"])

    return cfg


# ---------------------------------------------------------------------------
# GES loading
# ---------------------------------------------------------------------------

def _load_ges_conditions(ges_folder: Path, column: str, conditions: list) -> pd.DataFrame:
    """
    Load GES files for given conditions from one GES results folder.
    Returns a DataFrame: genes as index, one column per condition (NaN if file missing).
    """
    ges_dir = ges_folder / "data"
    series_list = []
    for cond in conditions:
        fname = ges_dir / f"ges_spec_{column}_{cond}.csv"
        if not fname.exists():
            print(f"  Warning: {fname.name} not found — skipping condition '{cond}'.")
            continue
        df = pd.read_csv(fname)
        if "gene" not in df.columns or "ges_score" not in df.columns:
            print(f"  Warning: unexpected columns in {fname.name} — skipping.")
            continue
        s = df.set_index("gene")["ges_score"].rename(str(cond))
        series_list.append(s)
        print(f"  Loaded {fname.name} ({len(s):,} genes)")

    if not series_list:
        raise ValueError(
            f"No GES files found in {ges_dir} for column='{column}', conditions={conditions}"
        )
    return pd.concat(series_list, axis=1)


def compute_axis_score(axis_cfg: dict) -> pd.Series:
    """
    score = mean GES(numerator_conditions) - mean GES(denominator_conditions).
    """
    ges_folder = axis_cfg["ges_results_folder"]
    column     = axis_cfg["column"]
    num_conds  = [str(c) for c in axis_cfg["numerator_conditions"]]
    den_conds  = [str(c) for c in axis_cfg["denominator_conditions"]]

    df_num = _load_ges_conditions(ges_folder, column, num_conds)
    df_den = _load_ges_conditions(ges_folder, column, den_conds)

    # Fill missing values with 0 (gene not expressed in that condition)
    num_mean = df_num.reindex(df_num.index).fillna(0.0).mean(axis=1)
    den_mean = df_den.reindex(df_den.index).fillna(0.0).mean(axis=1)

    all_genes = num_mean.index.union(den_mean.index)
    num_mean  = num_mean.reindex(all_genes).fillna(0.0)
    den_mean  = den_mean.reindex(all_genes).fillna(0.0)

    return (num_mean - den_mean).rename("score")


# ---------------------------------------------------------------------------
# Optional colouring helpers
# ---------------------------------------------------------------------------

def load_gene_groups(path) -> dict[str, set[str]]:
    df = pd.read_csv(path)
    missing = [c for c in ("gene", "group") if c not in df.columns]
    if missing:
        raise ValueError(f"Gene-list CSV '{path}' missing columns: {missing}")
    groups: dict[str, set[str]] = {}
    for name, sub in df.groupby("group"):
        groups[str(name)] = set(sub["gene"].astype(str))
    print(f"Gene groups loaded:")
    for name, genes in groups.items():
        print(f"  '{name}': {len(genes)} genes")
    return groups


def load_gene_clusters(path) -> pd.Series:
    df = pd.read_csv(path, index_col=0)
    if "cluster" not in df.columns:
        raise ValueError(f"gene_clusters.csv '{path}' missing 'cluster' column")
    return df["cluster"].astype(str)


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _add_quadrant_lines(ax):
    ax.axhline(0, color="#bbbbbb", linewidth=0.7, linestyle="--", zorder=0)
    ax.axvline(0, color="#bbbbbb", linewidth=0.7, linestyle="--", zorder=0)


def _annotate_top_genes(ax, x, y, gene_names, top_n: int) -> None:
    """Annotate the top_n genes per quadrant ranked by distance from origin."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    gene_names = np.asarray(gene_names)
    dist = np.sqrt(x ** 2 + y ** 2)

    annotated: set[str] = set()
    for x_pos, y_pos in [(True, True), (True, False), (False, True), (False, False)]:
        mask = (x > 0) == x_pos
        mask &= (y > 0) == y_pos
        if mask.sum() == 0:
            continue
        idx_quad = np.where(mask)[0]
        top_idx  = idx_quad[np.argsort(dist[idx_quad])[::-1][:top_n]]
        for i in top_idx:
            name = str(gene_names[i])
            if name not in annotated:
                ax.annotate(
                    name,
                    (x[i], y[i]),
                    fontsize=5,
                    ha="center",
                    va="bottom",
                    xytext=(0, 3),
                    textcoords="offset points",
                    alpha=0.8,
                )
                annotated.add(name)


def _scatter_base(ax, x, y, c, s=5, alpha=0.6, **kwargs):
    return ax.scatter(x, y, c=c, s=s, alpha=alpha, linewidths=0, **kwargs)


# ---------------------------------------------------------------------------
# Scatter variants
# ---------------------------------------------------------------------------

def _make_scatter_plain(x, y, gene_names, x_label, y_label, title, fig_dir, top_n):
    fig, ax = plt.subplots(figsize=(8, 7))
    _scatter_base(ax, x, y, "#888888", s=5, alpha=0.5)
    _add_quadrant_lines(ax)
    _annotate_top_genes(ax, x, y, gene_names, top_n)
    ax.set_xlabel(x_label, fontsize=10)
    ax.set_ylabel(y_label, fontsize=10)
    ax.set_title(title, fontsize=11)
    plt.tight_layout()
    fig.savefig(fig_dir / "gene_axes_scatter.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved: gene_axes_scatter.png")


def _make_scatter_by_clusters(scores_df, x_label, y_label, title, fig_dir,
                               cluster_series, top_n):
    common  = scores_df.index.intersection(cluster_series.index)
    s_inner = scores_df.loc[common]
    clusters = cluster_series.loc[common]

    unique_clusters = sorted(clusters.unique(), key=lambda c: int(c) if c.isdigit() else c)
    n_clusters = len(unique_clusters)
    palette = (
        list(sns.color_palette("tab20",  min(20, n_clusters)))
        + list(sns.color_palette("tab20b", max(0, n_clusters - 20)))
    )
    color_map = {c: palette[i % len(palette)] for i, c in enumerate(unique_clusters)}

    fig, ax = plt.subplots(figsize=(9, 7))

    # Gray for genes not in cluster file
    others = scores_df.index.difference(common)
    if len(others):
        _scatter_base(ax, scores_df.loc[others, "x_score"],
                      scores_df.loc[others, "y_score"],
                      "#cccccc", s=3, alpha=0.3)

    for cluster in unique_clusters:
        mask = clusters == cluster
        _scatter_base(ax,
                      s_inner.loc[mask, "x_score"], s_inner.loc[mask, "y_score"],
                      [color_map[cluster]], s=6, alpha=0.75,
                      label=f"Cluster {cluster} ({mask.sum():,})")

    _add_quadrant_lines(ax)
    _annotate_top_genes(ax,
                        s_inner["x_score"].values, s_inner["y_score"].values,
                        s_inner.index.tolist(), top_n)
    ax.legend(markerscale=2, fontsize=6, frameon=False, loc="upper left",
              title="Cluster", title_fontsize=7)
    ax.set_xlabel(x_label, fontsize=10)
    ax.set_ylabel(y_label, fontsize=10)
    ax.set_title(title, fontsize=11)
    plt.tight_layout()
    fig.savefig(fig_dir / "gene_axes_scatter.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved: gene_axes_scatter.png")


def _make_scatter_by_groups(scores_df, x_label, y_label, title, fig_dir,
                             gene_groups, grp_colors, top_n):
    """Combined scatter: all genes in background, each group overlaid in its colour."""
    x = scores_df["x_score"].values
    y = scores_df["y_score"].values
    gene_names = np.array(scores_df.index.tolist())

    fig, ax = plt.subplots(figsize=(9, 7))

    # Gray background (ungrouped genes)
    all_set_genes: set[str] = set()
    for gs in gene_groups.values():
        all_set_genes.update(gs)
    gray_mask = np.array([g not in all_set_genes for g in gene_names], dtype=bool)
    if gray_mask.any():
        _scatter_base(ax, x[gray_mask], y[gray_mask], "#cccccc", s=3, alpha=0.25)

    # Coloured groups (on top)
    for i, (grp_name, gene_set) in enumerate(gene_groups.items()):
        mask = np.array([g in gene_set for g in gene_names], dtype=bool)
        if not mask.any():
            continue
        n_found = mask.sum()
        n_total = len(gene_set)
        _scatter_base(ax, x[mask], y[mask], [grp_colors[i]], s=18, alpha=0.85,
                      label=f"{grp_name} ({n_found}/{n_total})")

    _add_quadrant_lines(ax)
    _annotate_top_genes(ax, x, y, gene_names, top_n)
    ax.legend(markerscale=1.5, fontsize=7, frameon=False, title="Gene group",
              title_fontsize=7)
    ax.set_xlabel(x_label, fontsize=10)
    ax.set_ylabel(y_label, fontsize=10)
    ax.set_title(title, fontsize=11)
    plt.tight_layout()
    fig.savefig(fig_dir / "gene_axes_scatter.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved: gene_axes_scatter.png")


def _make_per_group_panels(scores_df, x_label, y_label, title, fig_dir,
                            gene_groups, grp_colors, top_n):
    """
    One subplot per gene group.
    Each panel: all background genes in gray, that group's genes in colour with labels.
    """
    n_grp  = len(gene_groups)
    ncols  = min(3, n_grp)
    nrows  = (n_grp + ncols - 1) // ncols

    x_all  = scores_df["x_score"].values
    y_all  = scores_df["y_score"].values
    names_all = np.array(scores_df.index.tolist())

    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols * 5.0, nrows * 4.5),
                              squeeze=False)
    axes_flat = axes.flatten()

    for g_idx, (grp_name, gene_set) in enumerate(gene_groups.items()):
        ax = axes_flat[g_idx]
        in_grp  = np.array([g in gene_set for g in names_all], dtype=bool)
        n_found = int(in_grp.sum())
        n_total = len(gene_set)

        # Gray background
        ax.scatter(x_all[~in_grp], y_all[~in_grp],
                   c="#cccccc", s=2, alpha=0.2, linewidths=0)
        # Group genes
        if n_found:
            ax.scatter(x_all[in_grp], y_all[in_grp],
                       c=[grp_colors[g_idx]], s=20, alpha=0.85,
                       linewidths=0, zorder=2)
            _annotate_top_genes(ax,
                                x_all[in_grp], y_all[in_grp],
                                names_all[in_grp], top_n)

        _add_quadrant_lines(ax)
        ax.set_title(f"{grp_name}\n({n_found}/{n_total} in scatter)", fontsize=8)
        ax.set_xlabel(x_label, fontsize=7)
        ax.set_ylabel(y_label, fontsize=7)

    for idx in range(n_grp, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    fig.suptitle(title, fontsize=11, y=1.01)
    plt.tight_layout()
    fig.savefig(fig_dir / "gene_axes_scatter_per_group.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved: gene_axes_scatter_per_group.png")


def _make_gene_set_only_panels(scores_df, x_label, y_label, title, fig_dir,
                                gene_groups, grp_colors, top_n):
    """
    One subplot per gene group, showing only gene-set genes (no background).
    Each panel: other gene-set genes in light gray, that group's genes in colour with labels.
    """
    all_set_genes: set[str] = set()
    for gs in gene_groups.values():
        all_set_genes.update(gs)

    set_gene_names = np.array([g for g in scores_df.index if g in all_set_genes])
    if len(set_gene_names) == 0:
        print("  Warning: none of the gene-set genes found in scatter data — "
              "skipping gene_axes_scatter_gene_set_only.png")
        return

    scores_set = scores_df.loc[set_gene_names]
    x_set  = scores_set["x_score"].values
    y_set  = scores_set["y_score"].values

    n_grp  = len(gene_groups)
    ncols  = min(3, n_grp)
    nrows  = (n_grp + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols * 5.0, nrows * 4.5),
                              squeeze=False)
    axes_flat = axes.flatten()

    for g_idx, (grp_name, gene_set) in enumerate(gene_groups.items()):
        ax = axes_flat[g_idx]
        in_grp  = np.array([g in gene_set for g in set_gene_names], dtype=bool)
        n_found = int(in_grp.sum())
        n_total = len(gene_set)

        # Other gene-set genes → light gray
        if (~in_grp).any():
            ax.scatter(x_set[~in_grp], y_set[~in_grp],
                       c="#cccccc", s=10, alpha=0.4, linewidths=0)
        # This group → coloured
        if n_found:
            ax.scatter(x_set[in_grp], y_set[in_grp],
                       c=[grp_colors[g_idx]], s=25, alpha=0.9,
                       linewidths=0, zorder=2)
            _annotate_top_genes(ax,
                                x_set[in_grp], y_set[in_grp],
                                set_gene_names[in_grp], top_n)

        _add_quadrant_lines(ax)
        ax.set_title(f"{grp_name}\n({n_found}/{n_total} in scatter)", fontsize=8)
        ax.set_xlabel(x_label, fontsize=7)
        ax.set_ylabel(y_label, fontsize=7)

    for idx in range(n_grp, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    fig.suptitle(f"{title} — gene set only", fontsize=11, y=1.01)
    plt.tight_layout()
    fig.savefig(fig_dir / "gene_axes_scatter_gene_set_only.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved: gene_axes_scatter_gene_set_only.png")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_gene_axes_scatter(
    config_path: str,
    gene_list_path=None,
    gene_clusters_path=None,
) -> None:
    cfg = load_config(config_path)

    dataset_name = cfg["dataset_name"]
    out_root     = cfg["output_folder"]
    date_str     = datetime.datetime.now().strftime("%Y%m%d")
    run_dir      = out_root / f"{dataset_name}_{date_str}"
    run_dir.mkdir(parents=True, exist_ok=True)

    src = Path(cfg["_config_path"])
    dst = run_dir / src.name
    if src.resolve() != dst.resolve():
        shutil.copy2(src, dst)

    print(f"\n{'='*62}")
    print(f"  Gene Axes Scatter — {dataset_name}")
    print(f"  Output: {run_dir}")
    print(f"{'='*62}\n")

    # ── Compute axis scores ──────────────────────────────────────────────────
    print("Computing X-axis scores...")
    x_scores = compute_axis_score(cfg["x_axis"])

    print("\nComputing Y-axis scores...")
    y_scores = compute_axis_score(cfg["y_axis"])

    common_genes = x_scores.index.intersection(y_scores.index)
    print(f"\nGenes present in both axes: {len(common_genes):,}")

    scores_df = pd.DataFrame({
        "x_score": x_scores.reindex(common_genes).fillna(0.0),
        "y_score": y_scores.reindex(common_genes).fillna(0.0),
    })
    scores_df.index.name = "gene"

    # ── Save scores CSV ──────────────────────────────────────────────────────
    data_dir = run_dir / "data"
    data_dir.mkdir(exist_ok=True)
    scores_df.to_csv(data_dir / "gene_scores.csv")
    print(f"  Saved: gene_scores.csv ({len(scores_df):,} genes)")

    # ── Load optional colouring ──────────────────────────────────────────────
    gene_groups = None
    gene_list_path = gene_list_path or cfg.get("gene_list_path")
    if gene_list_path:
        print("\nLoading gene groups...")
        gene_groups = load_gene_groups(gene_list_path)

    gene_cluster_series = None
    gene_clusters_path  = gene_clusters_path or cfg.get("gene_clusters_csv")
    if gene_clusters_path and not gene_groups:
        print("\nLoading gene clusters...")
        gene_cluster_series = load_gene_clusters(gene_clusters_path)
        print(f"  {gene_cluster_series.nunique()} clusters")

    # ── Make scatter plot ────────────────────────────────────────────────────
    fig_dir  = run_dir / "figures"
    fig_dir.mkdir(exist_ok=True)

    x_label = cfg["x_axis"].get("label", "X score")
    y_label = cfg["y_axis"].get("label", "Y score")
    title   = cfg.get("title", dataset_name)
    top_n   = int(cfg.get("top_genes_per_quadrant", 10))

    print("\nGenerating scatter plots...")
    if gene_groups:
        grp_colors = [_GROUP_PALETTE[i % len(_GROUP_PALETTE)]
                      for i in range(len(gene_groups))]
        # 1. Combined: all genes + all groups overlaid
        _make_scatter_by_groups(scores_df, x_label, y_label, title, fig_dir,
                                gene_groups, grp_colors, top_n)
        # 2. Per-group panels: background genes in gray + one group highlighted per panel
        _make_per_group_panels(scores_df, x_label, y_label, title, fig_dir,
                               gene_groups, grp_colors, top_n)
        # 3. Gene-set-only panels: only gene-set genes shown, one group highlighted per panel
        _make_gene_set_only_panels(scores_df, x_label, y_label, title, fig_dir,
                                   gene_groups, grp_colors, top_n)
    elif gene_cluster_series is not None:
        _make_scatter_by_clusters(scores_df, x_label, y_label, title, fig_dir,
                                  gene_cluster_series, top_n)
    else:
        _make_scatter_plain(
            scores_df["x_score"].values, scores_df["y_score"].values,
            scores_df.index.tolist(),
            x_label, y_label, title, fig_dir, top_n,
        )

    print(f"\nDone → {run_dir}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="2-D gene axes scatter plot derived from GES scores."
    )
    parser.add_argument("config", help="YAML config file.")
    parser.add_argument(
        "--gene-list",
        default=None,
        metavar="CSV",
        help="CSV with 'gene','group' columns — colours genes by group.",
    )
    parser.add_argument(
        "--gene-clusters",
        default=None,
        metavar="CSV",
        help="gene_clusters.csv from gene_clustering pipeline — colours by Leiden cluster.",
    )
    args = parser.parse_args()
    run_gene_axes_scatter(args.config, args.gene_list, args.gene_clusters)
