"""
gene_clustering.py

Discover gene modules by clustering genes based on their GES (Gene Expression
Specificity) score profiles.

Each gene is represented as a vector of GES scores across all conditions in a
dataset. PCA → neighbor graph → UMAP → Leiden clustering groups genes with
similar expression specificity patterns.

For multi-dataset configs (combined runs), the pipeline is executed
independently for each dataset and results are saved in separate subfolders:

  results/gene_clusters/combined_v3_{YYYYMMDD}/
    config_used.yaml
    cortex/
      data/gene_clusters.csv, cluster_profiles.csv
      figures/umap_clusters.png, umap_all_conditions.png,
              cluster_heatmap.png, umap_{column}.png …
      metadata/pipeline_output.log
    all_layers/ …
    cell_phase/  …

For single-dataset configs the output is flat (no dataset subfolder):

  results/gene_clusters/cortex_{YYYYMMDD}/
    data/ figures/ metadata/

Two config modes
----------------
Single-dataset:
    ges_results_folder: …
    column_conditions:  {column: [conditions]}

Multi-dataset (combined):
    datasets:
      cortex:
        ges_results_folder: …
        column_conditions: {…}
      all_layers: …
      cell_phase: …

Usage:
  python modules/gene_clustering.py config_files/gene_clusters_combined_v3_config.yaml
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
import scanpy as sc
import seaborn as sns
import yaml


# Palette used consistently for gene groups across all figures.
_GROUP_PALETTE = (
    list(sns.color_palette("Set1",  9))
    + list(sns.color_palette("Set2",  8))
    + list(sns.color_palette("tab20", 20))
)


# ---------------------------------------------------------------------------
# Internal data structure
# ---------------------------------------------------------------------------

class DatasetSpec(NamedTuple):
    label: str
    ges_folder: Path
    column_conditions: dict


# ---------------------------------------------------------------------------
# Logging helper
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
# Gene-group loader
# ---------------------------------------------------------------------------

def load_gene_groups(path: str | Path) -> dict[str, set[str]]:
    """
    Load a CSV with columns 'gene' and 'group'.
    Returns {group_name: set_of_gene_names}.
    """
    df = pd.read_csv(path)
    missing = [c for c in ("gene", "group") if c not in df.columns]
    if missing:
        raise ValueError(
            f"Gene-list CSV '{path}' is missing column(s): {missing}. "
            "Expected columns: 'gene', 'group'."
        )
    groups: dict[str, set[str]] = {}
    for group_name, sub in df.groupby("group"):
        groups[str(group_name)] = set(sub["gene"].astype(str))
    print(f"Gene groups loaded from {Path(path).name}:")
    for name, genes in groups.items():
        print(f"  '{name}': {len(genes)} genes")
    return groups


# ---------------------------------------------------------------------------
# Config loader
# ---------------------------------------------------------------------------

def load_config(config_path: str) -> tuple[dict, list[DatasetSpec], bool]:
    """
    Returns (cfg, dataset_specs, multi_dataset).
    multi_dataset=True when the config uses the 'datasets' block.
    """
    config_path = Path(config_path).resolve()
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    if "dataset_name" not in cfg:
        raise ValueError("Config missing required key: 'dataset_name'")

    root = Path(cfg["ndd_gene_modules_folder_root"]).resolve()

    def resolve(p):
        p = Path(str(p))
        return p if p.is_absolute() else (root / p).resolve()

    cfg["output_folder"] = resolve(cfg.get("output_folder", "results/gene_clusters"))
    cfg["_config_path"] = config_path

    if "datasets" in cfg:
        specs = []
        for label, ds in cfg["datasets"].items():
            if "ges_results_folder" not in ds or "column_conditions" not in ds:
                raise ValueError(
                    f"Dataset '{label}' must have 'ges_results_folder' and 'column_conditions'."
                )
            specs.append(DatasetSpec(
                label=label,
                ges_folder=resolve(ds["ges_results_folder"]),
                column_conditions=ds["column_conditions"],
            ))

        # Parse optional combined_runs: {run_label: [dataset_label, ...]}
        combined_runs: dict[str, list] = {}
        if "combined_runs" in cfg:
            spec_by_label = {s.label: s for s in specs}
            for run_label, ds_names in cfg["combined_runs"].items():
                missing = [n for n in ds_names if n not in spec_by_label]
                if missing:
                    raise ValueError(
                        f"combined_runs '{run_label}' references unknown datasets: {missing}"
                    )
                combined_runs[run_label] = [spec_by_label[n] for n in ds_names]

        return cfg, specs, True, combined_runs

    # Single-dataset (backward compatible)
    if "ges_results_folder" not in cfg or "column_conditions" not in cfg:
        raise ValueError(
            "Config must have either 'datasets' (multi-dataset) or "
            "'ges_results_folder' + 'column_conditions' (single-dataset)."
        )
    specs = [DatasetSpec(
        label=cfg["dataset_name"],
        ges_folder=resolve(cfg["ges_results_folder"]),
        column_conditions=cfg["column_conditions"],
    )]
    return cfg, specs, False, {}


# ---------------------------------------------------------------------------
# Data loading  (always called per-dataset; multi_dataset flag unused here)
# ---------------------------------------------------------------------------

def load_combined_ges_matrix(
    specs: list,
    ges_score_threshold: float | None,
) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    """
    Load GES CSVs from multiple DatasetSpecs and merge into one gene × condition matrix.

    Column labels: '{dataset}__{column}__{condition}' to avoid clashes across datasets.

    Returns:
      matrix     — gene × condition DataFrame (NaN filled with 0)
      n_present  — Series: number of GES files each gene appears in
      col_meta   — DataFrame with columns [dataset, column, condition]
                   indexed by col_label
    """
    series_list: list[pd.Series] = []
    col_meta_rows: list[dict] = []

    for spec in specs:
        ges_dir = spec.ges_folder / "data"
        for column, conditions in spec.column_conditions.items():
            for condition in conditions:
                fname = ges_dir / f"ges_spec_{column}_{condition}.csv"
                if not fname.exists():
                    print(f"  ⚠️  Missing: {fname.name} — skipping.")
                    continue
                df = pd.read_csv(fname)
                if "gene" not in df.columns or "ges_score" not in df.columns:
                    print(f"  ⚠️  Unexpected columns in {fname.name} — skipping.")
                    continue
                col_label = f"{spec.label}__{column}__{condition}"
                s = df.set_index("gene")["ges_score"].rename(col_label)
                series_list.append(s)
                col_meta_rows.append({
                    "col_label": col_label,
                    "dataset":   spec.label,
                    "column":    column,
                    "condition": str(condition),
                })
                print(f"  ✔  [{spec.label}] {fname.name}  ({len(s):,} genes)")

    if not series_list:
        raise ValueError("No GES files found for any of the combined specs.")

    col_meta   = pd.DataFrame(col_meta_rows).set_index("col_label")
    matrix_raw = pd.concat(series_list, axis=1)
    n_present  = matrix_raw.notna().sum(axis=1)
    matrix     = matrix_raw.fillna(0.0)

    print(
        f"\nCombined matrix: {matrix.shape[0]:,} genes × {matrix.shape[1]} conditions"
    )

    if ges_score_threshold is not None:
        before = len(matrix)
        keep      = (matrix >= ges_score_threshold).any(axis=1)
        matrix    = matrix[keep]
        n_present = n_present[keep]
        print(f"GES threshold (>= {ges_score_threshold}): {before:,} → {len(matrix):,} genes")

    return matrix, n_present, col_meta


def load_ges_matrix(
    spec: DatasetSpec,
    ges_score_threshold: float | None,
) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    """
    Load GES CSVs for one dataset and build a gene × condition matrix.

    Column labels: '{column}__{condition}'.

    Returns:
      matrix     — gene × condition DataFrame (NaN filled with 0)
      n_present  — Series: number of GES files each gene appears in
      col_meta   — DataFrame with columns [dataset, column, condition]
                   indexed by col_label
    """
    ges_dir = spec.ges_folder / "data"
    series_list: list[pd.Series] = []
    col_meta_rows: list[dict] = []

    for column, conditions in spec.column_conditions.items():
        for condition in conditions:
            fname = ges_dir / f"ges_spec_{column}_{condition}.csv"
            if not fname.exists():
                print(f"  ⚠️  Missing: {fname.name} — skipping.")
                continue
            df = pd.read_csv(fname)
            if "gene" not in df.columns or "ges_score" not in df.columns:
                print(f"  ⚠️  Unexpected columns in {fname.name} — skipping.")
                continue
            col_label = f"{column}__{condition}"
            s = df.set_index("gene")["ges_score"].rename(col_label)
            series_list.append(s)
            col_meta_rows.append({
                "col_label": col_label,
                "dataset":   spec.label,
                "column":    column,
                "condition": str(condition),
            })
            print(f"  ✔  {fname.name}  ({len(s):,} genes)")

    if not series_list:
        raise ValueError(
            f"No GES files found in {ges_dir}. "
            "Check ges_results_folder and column_conditions in the config."
        )

    col_meta   = pd.DataFrame(col_meta_rows).set_index("col_label")
    matrix_raw = pd.concat(series_list, axis=1)
    n_present  = matrix_raw.notna().sum(axis=1)
    matrix     = matrix_raw.fillna(0.0)

    print(
        f"\nMatrix: {matrix.shape[0]:,} genes × {matrix.shape[1]} conditions"
    )

    if ges_score_threshold is not None:
        before = len(matrix)
        keep      = (matrix >= ges_score_threshold).any(axis=1)
        matrix    = matrix[keep]
        n_present = n_present[keep]
        print(f"GES threshold (>= {ges_score_threshold}): {before:,} → {len(matrix):,} genes")

    return matrix, n_present, col_meta


# ---------------------------------------------------------------------------
# Gene filtering + scaling
# ---------------------------------------------------------------------------

def filter_genes(
    matrix: pd.DataFrame,
    n_present: pd.Series,
    min_conditions: int,
) -> pd.DataFrame:
    before = len(matrix)
    matrix = matrix[n_present >= min_conditions]
    print(
        f"Min-conditions filter (in >= {min_conditions} file(s)): "
        f"{before:,} → {len(matrix):,} genes"
    )
    return matrix


def scale_matrix(matrix: pd.DataFrame, mode: str) -> pd.DataFrame:
    """
    'gene'      — z-score each gene across conditions (default).
    'condition' — z-score each condition across genes.
    'none'      — raw GES scores.
    """
    if mode == "none":
        return matrix
    if mode == "condition":
        mean = matrix.mean(axis=0)
        std  = matrix.std(axis=0).replace(0, 1)
        return (matrix - mean) / std
    mean = matrix.mean(axis=1)
    std  = matrix.std(axis=1).replace(0, 1)
    return matrix.subtract(mean, axis=0).div(std, axis=0)


# ---------------------------------------------------------------------------
# UMAP + Leiden clustering
# ---------------------------------------------------------------------------

def cluster_genes(
    matrix_scaled: pd.DataFrame,
    pca_n_components: int,
    n_neighbors: int,
    umap_min_dist: float,
    leiden_resolution: float,
    random_state: int = 42,
) -> sc.AnnData:
    import anndata as ad

    adata = ad.AnnData(
        X=matrix_scaled.values.astype(np.float32),
        obs=pd.DataFrame(index=matrix_scaled.index),
        var=pd.DataFrame(index=matrix_scaled.columns),
    )

    n_comps = min(pca_n_components, adata.n_vars - 1, adata.n_obs - 1)
    print(f"\nPCA: {n_comps} components  ({adata.n_obs:,} genes × {adata.n_vars} conditions)")
    sc.pp.pca(adata, n_comps=n_comps, random_state=random_state)

    print(f"Building neighbor graph: k={n_neighbors}")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_pca",
                    random_state=random_state)

    print("UMAP...")
    sc.tl.umap(adata, min_dist=umap_min_dist, random_state=random_state)

    print(f"Leiden clustering (resolution={leiden_resolution})")
    sc.tl.leiden(adata, resolution=leiden_resolution, random_state=random_state,
                 flavor="igraph", n_iterations=2, directed=False)

    n_clusters = adata.obs["leiden"].nunique()
    print(f"  → {n_clusters} clusters found")
    for cluster, count in adata.obs["leiden"].value_counts().sort_index().items():
        print(f"     cluster {cluster:>3s}: {count:,} genes")

    return adata


# ---------------------------------------------------------------------------
# Save tabular results
# ---------------------------------------------------------------------------

def save_results(
    adata: sc.AnnData,
    matrix_raw: pd.DataFrame,
    out_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    data_dir = out_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    umap_df = pd.DataFrame(
        adata.obsm["X_umap"],
        index=adata.obs.index,
        columns=["UMAP_1", "UMAP_2"],
    )
    gene_clusters = (
        pd.DataFrame({"cluster": adata.obs["leiden"]}, index=adata.obs.index)
        .join(umap_df)
        .join(matrix_raw.reindex(adata.obs.index))
    )
    gene_clusters.index.name = "gene"
    gene_clusters.to_csv(data_dir / "gene_clusters.csv")
    print(f"  📄 gene_clusters.csv  ({len(gene_clusters):,} genes)")

    ges_cols = matrix_raw.columns.tolist()
    profiles = gene_clusters[["cluster"] + ges_cols].groupby("cluster", observed=True)[ges_cols].mean()
    profiles.to_csv(data_dir / "cluster_profiles.csv")
    print(f"  📄 cluster_profiles.csv  ({len(profiles)} clusters × {len(profiles.columns)} conditions)")

    return gene_clusters, profiles


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _draw_gene_group_panel(
    ax,
    umap_xy: np.ndarray,
    obs_index,          # iterable of gene names matching umap_xy rows
    gene_set: set[str],
    group_name: str,
    color,
) -> None:
    """
    Draw a single gene-group overlay panel:
      - All genes not in gene_set → light gray, tiny, translucent.
      - Genes in gene_set found in the UMAP → colored, larger, opaque.
    """
    in_group = np.array([g in gene_set for g in obs_index], dtype=bool)
    n_found  = int(in_group.sum())
    n_total  = len(gene_set)

    # Gray background
    ax.scatter(
        umap_xy[~in_group, 0], umap_xy[~in_group, 1],
        c="#cccccc", s=1, alpha=0.25, rasterized=True, linewidths=0,
    )
    # Coloured highlight
    if n_found:
        ax.scatter(
            umap_xy[in_group, 0], umap_xy[in_group, 1],
            c=[color], s=18, alpha=0.9, linewidths=0, zorder=2,
        )
    ax.set_title(f"{group_name}\n({n_found}/{n_total} in UMAP)", fontsize=7, pad=2)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_frame_on(False)

def make_figures(
    adata: sc.AnnData,
    profiles: pd.DataFrame,
    matrix_raw: pd.DataFrame,
    col_meta: pd.DataFrame,
    out_dir: Path,
    title_prefix: str,
    gene_groups: dict[str, set[str]] | None = None,
) -> None:
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    n_clusters   = adata.obs["leiden"].nunique()
    umap_xy      = adata.obsm["X_umap"]
    cluster_vals = adata.obs["leiden"].values
    ges_raw      = matrix_raw.reindex(adata.obs.index)
    obs_index    = list(adata.obs.index)

    groups       = gene_groups or {}
    n_grp        = len(groups)
    grp_names    = list(groups.keys())
    grp_colors   = [_GROUP_PALETTE[i % len(_GROUP_PALETTE)] for i in range(n_grp)]

    # ── 1. UMAP coloured by Leiden cluster (+ gene-group panels) ────────────
    cluster_palette = (
        sns.color_palette("tab20",  min(20, n_clusters))
        + sns.color_palette("tab20b", max(0, n_clusters - 20))
    )
    total_panels = 1 + n_grp
    ncols_fig    = min(4, total_panels)
    nrows_fig    = (total_panels + ncols_fig - 1) // ncols_fig
    fig, axes    = plt.subplots(
        nrows_fig, ncols_fig,
        figsize=(ncols_fig * 3.8, nrows_fig * 3.5),
        squeeze=False,
    )
    axes_flat = axes.flatten()

    ax0 = axes_flat[0]
    for i, cluster in enumerate(sorted(adata.obs["leiden"].unique(), key=int)):
        mask = cluster_vals == cluster
        ax0.scatter(
            umap_xy[mask, 0], umap_xy[mask, 1],
            s=4, alpha=0.7, color=cluster_palette[i % len(cluster_palette)],
            label=f"{cluster} ({mask.sum():,})",
        )
    ax0.set_title(
        f"Leiden clusters — {title_prefix}\n({adata.n_obs:,} genes, {n_clusters} clusters)",
        fontsize=7,
    )
    ax0.legend(
        markerscale=2, title="Cluster", title_fontsize=6,
        loc="upper left", fontsize=5, frameon=False,
    )
    ax0.set_frame_on(False)
    ax0.set_xticks([]); ax0.set_yticks([])

    for g_idx, (grp_name, gene_set) in enumerate(groups.items()):
        _draw_gene_group_panel(
            axes_flat[1 + g_idx], umap_xy, obs_index,
            gene_set, grp_name, grp_colors[g_idx],
        )

    for idx in range(total_panels, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    plt.tight_layout()
    fig.savefig(fig_dir / "umap_clusters.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  🖼  umap_clusters.png")

    # ── 2. All-conditions UMAP grid ──────────────────────────────────────────
    _make_all_conditions_umap(
        ges_raw, umap_xy, obs_index, col_meta, fig_dir, title_prefix, groups, grp_colors,
    )

    # ── 3. Cluster profile heatmap ───────────────────────────────────────────
    _make_cluster_heatmap(profiles, col_meta, fig_dir, title_prefix)

    # ── 4. Per-column UMAPs (+ gene-group panels) ────────────────────────────
    for column in col_meta["column"].unique():
        cond_cols  = col_meta[col_meta["column"] == column].index.tolist()
        conditions = col_meta.loc[cond_cols, "condition"].tolist()

        n_cond    = len(cond_cols)
        n_total   = n_cond + n_grp
        ncols_plot = min(4, n_total)
        nrows_plot = (n_total + ncols_plot - 1) // ncols_plot
        fig, axes = plt.subplots(
            nrows_plot, ncols_plot,
            figsize=(ncols_plot * 3.5, nrows_plot * 3.2),
            squeeze=False,
        )
        axes_flat = axes.flatten()

        for idx, (col_key, cond_label) in enumerate(zip(cond_cols, conditions)):
            ax = axes_flat[idx]
            vals = ges_raw[col_key].values
            nonzero = vals[vals != 0]
            vmax = float(np.percentile(np.abs(nonzero), 95)) if len(nonzero) else 1.0
            sc_plot = ax.scatter(
                umap_xy[:, 0], umap_xy[:, 1],
                c=vals, s=2, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
            )
            ax.set_title(cond_label, fontsize=8)
            ax.set_frame_on(False)
            ax.set_xticks([]); ax.set_yticks([])
            plt.colorbar(sc_plot, ax=ax, shrink=0.6, label="GES")

        for g_idx, (grp_name, gene_set) in enumerate(groups.items()):
            _draw_gene_group_panel(
                axes_flat[n_cond + g_idx], umap_xy, obs_index,
                gene_set, grp_name, grp_colors[g_idx],
            )

        for idx in range(n_total, len(axes_flat)):
            axes_flat[idx].set_visible(False)

        fig.suptitle(f"GES — {column} ({title_prefix})", fontsize=10)
        plt.tight_layout()
        safe_col = column.replace(" ", "_")
        fname = f"umap_{safe_col}.png"
        fig.savefig(fig_dir / fname, dpi=130, bbox_inches="tight")
        plt.close(fig)
        print(f"  🖼  {fname}")


def _make_all_conditions_umap(
    ges_raw: pd.DataFrame,
    umap_xy: np.ndarray,
    obs_index: list,
    col_meta: pd.DataFrame,
    fig_dir: Path,
    title_prefix: str,
    gene_groups: dict[str, set[str]],
    grp_colors: list,
) -> None:
    """Single figure with one subplot per condition + one per gene group."""
    all_cols = col_meta.index.tolist()
    n_cond   = len(all_cols)
    n_grp    = len(gene_groups)
    n_total  = n_cond + n_grp

    ncols_plot = min(6, n_total)
    nrows_plot = (n_total + ncols_plot - 1) // ncols_plot

    columns_in_meta = col_meta["column"].unique().tolist()
    col_palette     = dict(zip(columns_in_meta, sns.color_palette("tab10", len(columns_in_meta))))

    fig, axes = plt.subplots(
        nrows_plot, ncols_plot,
        figsize=(ncols_plot * 2.8, nrows_plot * 2.6),
        squeeze=False,
    )
    axes_flat = axes.flatten()

    # Condition panels
    for idx, col_key in enumerate(all_cols):
        ax = axes_flat[idx]
        vals = ges_raw[col_key].values
        nonzero = vals[vals != 0]
        vmax = float(np.percentile(np.abs(nonzero), 95)) if len(nonzero) else 1.0
        sc_plot = ax.scatter(
            umap_xy[:, 0], umap_xy[:, 1],
            c=vals, s=1, cmap="RdBu_r", vmin=-vmax, vmax=vmax, rasterized=True,
        )
        ax.set_title(col_meta.loc[col_key, "condition"], fontsize=6.5, pad=2)
        ax.set_xticks([]); ax.set_yticks([])

        border_color = col_palette[col_meta.loc[col_key, "column"]]
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_edgecolor(border_color)
            spine.set_linewidth(1.5)

        plt.colorbar(sc_plot, ax=ax, shrink=0.55, pad=0.02)

    # Gene-group panels
    for g_idx, (grp_name, gene_set) in enumerate(gene_groups.items()):
        _draw_gene_group_panel(
            axes_flat[n_cond + g_idx], umap_xy, obs_index,
            gene_set, grp_name, grp_colors[g_idx],
        )

    for idx in range(n_total, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    from matplotlib.patches import Patch
    legend_handles = [Patch(color=c, label=col) for col, c in col_palette.items()]
    if gene_groups:
        legend_handles += [Patch(facecolor="white", label="")]  # spacer
        legend_handles += [
            Patch(color=grp_colors[i], label=name)
            for i, name in enumerate(gene_groups)
        ]
    fig.legend(
        handles=legend_handles, title="Column / Gene group",
        loc="lower right", bbox_to_anchor=(1.0, 0.0),
        fontsize=7, frameon=False,
    )
    fig.suptitle(f"GES scores — all conditions ({title_prefix})", fontsize=11, y=1.01)
    plt.tight_layout()
    fig.savefig(fig_dir / "umap_all_conditions.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  🖼  umap_all_conditions.png")


def _make_cluster_heatmap(
    profiles: pd.DataFrame,
    col_meta: pd.DataFrame,
    fig_dir: Path,
    title_prefix: str,
) -> None:
    # Display labels: just the condition name (strip column prefix)
    short_labels = {c: col_meta.loc[c, "condition"] for c in profiles.columns}
    plot_profiles = profiles.rename(columns=short_labels)

    columns_in_meta = col_meta["column"].unique().tolist()
    col_palette     = dict(zip(columns_in_meta, sns.color_palette("Paired", len(columns_in_meta))))
    col_colors      = pd.Series(
        [col_palette[col_meta.loc[c, "column"]] for c in profiles.columns],
        index=plot_profiles.columns,
        name="Column",
    )

    fig_w = max(10, len(profiles.columns) * 0.45 + 3)
    fig_h = max(5,  len(profiles)         * 0.5  + 3)

    g = sns.clustermap(
        plot_profiles,
        col_colors=col_colors,
        row_cluster=True,
        col_cluster=False,
        cmap="RdBu_r",
        center=0,
        figsize=(fig_w, fig_h),
        xticklabels=True,
        yticklabels=[f"Cluster {c}" for c in plot_profiles.index],
        cbar_pos=(0.02, 0.82, 0.03, 0.12),
        linewidths=0.2,
        dendrogram_ratio=0.1,
    )
    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=7
    )
    g.ax_heatmap.set_title(f"Cluster expression profiles — {title_prefix}", pad=14)

    from matplotlib.patches import Patch
    handles = [Patch(color=c, label=col) for col, c in col_palette.items()]
    g.ax_col_colors.legend(
        handles=handles, loc="upper right",
        bbox_to_anchor=(1.18, 1.5),
        fontsize=7, title="Column", frameon=False,
    )

    g.savefig(fig_dir / "cluster_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  🖼  cluster_heatmap.png")


# ---------------------------------------------------------------------------
# Single-dataset pipeline (called once per dataset)
# ---------------------------------------------------------------------------

def _run_one_dataset(
    spec: DatasetSpec,
    params: dict,
    out_dir: Path,
    gene_groups: dict[str, set[str]] | None = None,
) -> None:
    metadata_dir = out_dir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    with _log_to_file(metadata_dir / "pipeline_output.log"):
        print(f"\n{'='*62}")
        print(f"  Dataset: {spec.label}")
        print(f"{'='*62}")
        print(f"• GES folder:          {spec.ges_folder}")
        print(f"• Output:              {out_dir}")
        print(f"• Scale mode:          {params['scale_mode']}")
        print(f"• Min conditions:      {params['min_conditions']}")
        ges_thr = params['ges_threshold']
        print(f"• GES threshold:       {ges_thr if ges_thr is not None else '(none)'}")
        print(f"• PCA components:      {params['pca_n_components']}")
        print(f"• UMAP neighbors:      {params['n_neighbors']}")
        print(f"• UMAP min_dist:       {params['umap_min_dist']}")
        print(f"• Leiden resolution:   {params['leiden_resolution']}")
        if gene_groups:
            print(f"• Gene groups:         {list(gene_groups.keys())}")
        print(f"{'='*62}\n")

        print("Loading GES scores...")
        matrix_raw, n_present, col_meta = load_ges_matrix(spec, ges_thr)

        matrix_raw = filter_genes(matrix_raw, n_present, params["min_conditions"])

        print(f"\nScaling (mode='{params['scale_mode']}')...")
        matrix_scaled = scale_matrix(matrix_raw, params["scale_mode"])

        adata = cluster_genes(
            matrix_scaled,
            params["pca_n_components"],
            params["n_neighbors"],
            params["umap_min_dist"],
            params["leiden_resolution"],
        )

        print("\nSaving results...")
        _, profiles = save_results(adata, matrix_raw, out_dir)

        print("\nGenerating figures...")
        make_figures(
            adata, profiles, matrix_raw, col_meta, out_dir, spec.label,
            gene_groups=gene_groups,
        )

        print(f"\n✔ {spec.label} complete → {out_dir}")


# ---------------------------------------------------------------------------
# Combined-datasets pipeline (merges multiple specs into one analysis)
# ---------------------------------------------------------------------------

def _run_combined_datasets(
    specs: list,
    label: str,
    params: dict,
    out_dir: Path,
    gene_groups: dict[str, set[str]] | None = None,
) -> None:
    """Run the full clustering pipeline on a merged matrix from multiple DatasetSpecs."""
    metadata_dir = out_dir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    log_path = metadata_dir / "pipeline_output.log"
    with _log_to_file(log_path):
        ges_thr = params["ges_threshold"]
        print(f"\n{'='*62}")
        print(f"  Combined dataset run: {label}")
        print(f"  Datasets: {[s.label for s in specs]}")
        print(f"• Output:              {out_dir}")
        print(f"• Scale mode:          {params['scale_mode']}")
        print(f"• Min conditions:      {params['min_conditions']}")
        print(f"• GES threshold:       {ges_thr if ges_thr is not None else '(none)'}")
        print(f"• PCA components:      {params['pca_n_components']}")
        print(f"• UMAP neighbors:      {params['n_neighbors']}")
        print(f"• UMAP min_dist:       {params['umap_min_dist']}")
        print(f"• Leiden resolution:   {params['leiden_resolution']}")
        if gene_groups:
            print(f"• Gene groups:         {list(gene_groups.keys())}")
        print(f"{'='*62}\n")

        print("Loading combined GES scores...")
        matrix_raw, n_present, col_meta = load_combined_ges_matrix(specs, ges_thr)

        matrix_raw = filter_genes(matrix_raw, n_present, params["min_conditions"])

        print(f"\nScaling (mode='{params['scale_mode']}')...")
        matrix_scaled = scale_matrix(matrix_raw, params["scale_mode"])

        adata = cluster_genes(
            matrix_scaled,
            params["pca_n_components"],
            params["n_neighbors"],
            params["umap_min_dist"],
            params["leiden_resolution"],
        )

        print("\nSaving results...")
        _, profiles = save_results(adata, matrix_raw, out_dir)

        print("\nGenerating figures...")
        make_figures(
            adata, profiles, matrix_raw, col_meta, out_dir, label,
            gene_groups=gene_groups,
        )

        print(f"\n✔ {label} complete → {out_dir}")


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def run_gene_clustering(
    config_path: str,
    gene_list_path: str | None = None,
) -> None:
    cfg, specs, multi_dataset, combined_runs = load_config(config_path)

    dataset_name = cfg["dataset_name"]
    out_root     = cfg["output_folder"]
    date_str     = datetime.datetime.now().strftime("%Y%m%d")
    run_dir      = out_root / f"{dataset_name}_{date_str}"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Copy config once to the top-level run folder
    src = Path(cfg["_config_path"])
    dst = run_dir / src.name
    if src.resolve() != dst.resolve():
        shutil.copy2(src, dst)

    params = {
        "min_conditions":    int(cfg.get("min_conditions_expressed", 2)),
        "ges_threshold":     float(cfg["ges_score_threshold"]) if cfg.get("ges_score_threshold") else None,
        "scale_mode":        cfg.get("scale", "gene"),
        "pca_n_components":  int(cfg.get("pca_n_components", 50)),
        "n_neighbors":       int(cfg.get("umap_n_neighbors", 15)),
        "umap_min_dist":     float(cfg.get("umap_min_dist", 0.3)),
        "leiden_resolution": float(cfg.get("leiden_resolution", 0.5)),
    }

    # Load optional gene groups (CLI arg takes priority over config key)
    gene_list_path = gene_list_path or cfg.get("gene_list_path")
    gene_groups: dict[str, set[str]] | None = None
    if gene_list_path:
        gene_groups = load_gene_groups(gene_list_path)

    mode_str = "multi-dataset" if multi_dataset else "single-dataset"
    print(f"\n{'='*62}")
    print(f"  Gene Clustering Pipeline — {dataset_name}  ({mode_str})")
    print(f"{'='*62}")
    for s in specs:
        print(f"  · {s.label:<16} {s.ges_folder}")
    if combined_runs:
        print(f"  Combined runs:")
        for run_label, cspecs in combined_runs.items():
            print(f"    · {run_label}: {[s.label for s in cspecs]}")
    print(f"• Output root: {run_dir}")
    if gene_groups:
        print(f"• Gene groups: {list(gene_groups.keys())}")
    print(f"{'='*62}\n")

    for spec in specs:
        out_dir = run_dir / spec.label if multi_dataset else run_dir
        _run_one_dataset(spec, params, out_dir, gene_groups=gene_groups)

    for run_label, cspecs in combined_runs.items():
        out_dir = run_dir / run_label
        _run_combined_datasets(cspecs, run_label, params, out_dir, gene_groups=gene_groups)

    print(f"\n🎉 All datasets complete → {run_dir}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Cluster genes by GES expression profile (UMAP + Leiden)."
    )
    parser.add_argument("config", help="Path to the YAML config file.")
    parser.add_argument(
        "--gene-list",
        default=None,
        metavar="CSV",
        help=(
            "Optional CSV with columns 'gene' and 'group'. "
            "For each group, a highlighted overlay panel is added to every UMAP figure "
            "showing those genes in colour against a gray background."
        ),
    )
    args = parser.parse_args()
    run_gene_clustering(args.config, gene_list_path=args.gene_list)
