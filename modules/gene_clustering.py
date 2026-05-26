"""
gene_clustering.py

Discover gene modules by clustering genes based on their GES (Gene Expression
Specificity) score profiles.

Each gene is represented as a vector of GES scores across all conditions from
one or more datasets. PCA → neighbor graph → UMAP → Leiden clustering groups
genes with similar expression specificity patterns.

Two config modes
----------------
Single-dataset (backward compatible):
    ges_results_folder: ...
    column_conditions:  {column: [conditions]}
    Column labels in output: '{column}__{condition}'

Multi-dataset (combined analysis):
    datasets:
      cortex:
        ges_results_folder: ...
        column_conditions: {...}
      all_layers:
        ges_results_folder: ...
        column_conditions: {...}
      cell_phase:
        ...
    Column labels in output: '{dataset}__{column}__{condition}'

Outputs in results/gene_clusters/{dataset_name}_{YYYYMMDD}/:
  data/
    gene_clusters.csv        — gene, cluster, UMAP_1, UMAP_2, all GES features
    cluster_profiles.csv     — mean GES per (cluster × condition)
  figures/
    umap_clusters.png        — UMAP coloured by Leiden cluster
    cluster_heatmap.png      — clustermap: clusters × conditions (mean GES)
    umap_{dataset}_{column}.png  — one per (dataset, column) coloured by GES
  metadata/
    config_used.yaml
    pipeline_output.log

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
        # Multi-dataset mode
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
        return cfg, specs, True

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
    return cfg, specs, False


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _col_label(dataset_label: str, column: str, condition, multi_dataset: bool) -> str:
    cond_str = str(condition)
    if multi_dataset:
        return f"{dataset_label}__{column}__{cond_str}"
    return f"{column}__{cond_str}"


def load_ges_matrix(
    specs: list[DatasetSpec],
    multi_dataset: bool,
    ges_score_threshold: float | None,
) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    """
    Load GES CSVs from all dataset specs and build a gene × condition matrix.

    Returns:
      matrix     — gene × condition DataFrame (NaN filled with 0)
      n_present  — Series: number of GES files each gene appears in
      col_meta   — DataFrame describing each matrix column:
                   columns: col_label, dataset, column, condition
    """
    series_list: list[pd.Series] = []
    col_meta_rows: list[dict] = []

    for spec in specs:
        ges_dir = spec.ges_folder / "data"
        print(f"\n  Dataset '{spec.label}'  ← {ges_dir}")

        for column, conditions in spec.column_conditions.items():
            for condition in conditions:
                fname = ges_dir / f"ges_spec_{column}_{condition}.csv"
                if not fname.exists():
                    print(f"    ⚠️  Missing: {fname.name} — skipping.")
                    continue
                df = pd.read_csv(fname)
                if "gene" not in df.columns or "ges_score" not in df.columns:
                    print(f"    ⚠️  Unexpected columns in {fname.name} — skipping.")
                    continue
                label = _col_label(spec.label, column, condition, multi_dataset)
                s = df.set_index("gene")["ges_score"].rename(label)
                series_list.append(s)
                col_meta_rows.append({
                    "col_label": label,
                    "dataset":   spec.label,
                    "column":    column,
                    "condition": str(condition),
                })
                print(f"    ✔  {fname.name}  ({len(s):,} genes)")

    if not series_list:
        raise ValueError("No GES files could be loaded. Check config paths and column_conditions.")

    col_meta = pd.DataFrame(col_meta_rows).set_index("col_label")
    matrix_raw = pd.concat(series_list, axis=1)
    n_present  = matrix_raw.notna().sum(axis=1)
    matrix     = matrix_raw.fillna(0.0)

    print(
        f"\nCombined matrix: {matrix.shape[0]:,} genes × {matrix.shape[1]} conditions"
        f" ({len(specs)} dataset(s))"
    )

    if ges_score_threshold is not None:
        before = len(matrix)
        keep   = (matrix >= ges_score_threshold).any(axis=1)
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
    'gene'      — z-score each gene across conditions (highlights which conditions
                  it is specific to, relative to its own average; default).
    'condition' — z-score each condition across genes.
    'none'      — raw GES scores.
    """
    if mode == "none":
        return matrix
    if mode == "condition":
        mean = matrix.mean(axis=0)
        std  = matrix.std(axis=0).replace(0, 1)
        return (matrix - mean) / std
    # gene-wise (default)
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
    sc.tl.leiden(adata, resolution=leiden_resolution, random_state=random_state)

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
    run_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    data_dir = run_dir / "data"
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
    profiles = gene_clusters[["cluster"] + ges_cols].groupby("cluster")[ges_cols].mean()
    profiles.to_csv(data_dir / "cluster_profiles.csv")
    print(f"  📄 cluster_profiles.csv  ({len(profiles)} clusters × {len(profiles.columns)} conditions)")

    return gene_clusters, profiles


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def make_figures(
    adata: sc.AnnData,
    profiles: pd.DataFrame,
    matrix_raw: pd.DataFrame,
    col_meta: pd.DataFrame,
    run_dir: Path,
    dataset_name: str,
    specs: list[DatasetSpec],
    multi_dataset: bool,
) -> None:
    fig_dir = run_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    n_clusters   = adata.obs["leiden"].nunique()
    umap_xy      = adata.obsm["X_umap"]
    cluster_vals = adata.obs["leiden"].values
    ges_raw      = matrix_raw.reindex(adata.obs.index)

    # ── 1. UMAP coloured by Leiden cluster ──────────────────────────────────
    palette = (
        sns.color_palette("tab20",  min(20, n_clusters))
        + sns.color_palette("tab20b", max(0, n_clusters - 20))
    )
    fig, ax = plt.subplots(figsize=(8, 7))
    for i, cluster in enumerate(sorted(adata.obs["leiden"].unique(), key=int)):
        mask = cluster_vals == cluster
        ax.scatter(
            umap_xy[mask, 0], umap_xy[mask, 1],
            s=4, alpha=0.7, color=palette[i % len(palette)],
            label=f"{cluster} ({mask.sum():,})",
        )
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_title(
        f"Gene modules — {dataset_name}\n"
        f"({adata.n_obs:,} genes, {n_clusters} Leiden clusters)"
    )
    ax.legend(
        markerscale=3, title="Cluster (n genes)",
        bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7,
    )
    ax.set_frame_on(False)
    plt.tight_layout()
    fig.savefig(fig_dir / "umap_clusters.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  🖼  umap_clusters.png")

    # ── 2. Cluster profile heatmap ───────────────────────────────────────────
    _make_cluster_heatmap(profiles, col_meta, fig_dir, dataset_name, multi_dataset)

    # ── 3. All-conditions UMAP grid ──────────────────────────────────────────
    _make_all_conditions_umap(ges_raw, umap_xy, col_meta, fig_dir, dataset_name)

    # ── 4. Per-(dataset, column) UMAPs coloured by GES ──────────────────────
    for spec in specs:
        for column, conditions in spec.column_conditions.items():
            cond_cols = [
                _col_label(spec.label, column, c, multi_dataset)
                for c in conditions
                if _col_label(spec.label, column, c, multi_dataset) in matrix_raw.columns
            ]
            if not cond_cols:
                continue

            n = len(cond_cols)
            ncols_plot = min(4, n)
            nrows_plot = (n + ncols_plot - 1) // ncols_plot
            fig, axes = plt.subplots(
                nrows_plot, ncols_plot,
                figsize=(ncols_plot * 3.5, nrows_plot * 3.2),
                squeeze=False,
            )
            axes_flat = axes.flatten()

            for idx, col_key in enumerate(cond_cols):
                ax = axes_flat[idx]
                vals = ges_raw[col_key].values
                nonzero = vals[vals != 0]
                vmax = float(np.percentile(np.abs(nonzero), 95)) if len(nonzero) else 1.0
                sc_plot = ax.scatter(
                    umap_xy[:, 0], umap_xy[:, 1],
                    c=vals, s=2, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                )
                cond_label = str(conditions[idx])
                ax.set_title(cond_label, fontsize=8)
                ax.set_frame_on(False)
                ax.set_xticks([]); ax.set_yticks([])
                plt.colorbar(sc_plot, ax=ax, shrink=0.6, label="GES")

            for idx in range(n, len(axes_flat)):
                axes_flat[idx].set_visible(False)

            title = f"GES — {spec.label} / {column} ({dataset_name})"
            fig.suptitle(title, fontsize=10)
            plt.tight_layout()
            safe_col = column.replace(" ", "_")
            fname = f"umap_{spec.label}_{safe_col}.png"
            fig.savefig(fig_dir / fname, dpi=130, bbox_inches="tight")
            plt.close(fig)
            print(f"  🖼  {fname}")


def _make_all_conditions_umap(
    ges_raw: pd.DataFrame,
    umap_xy: np.ndarray,
    col_meta: pd.DataFrame,
    fig_dir: Path,
    dataset_name: str,
) -> None:
    """
    Single figure with one subplot per condition (all datasets combined).
    Subplots are grouped and titled by dataset + column.
    """
    all_cols = col_meta.index.tolist()
    n = len(all_cols)
    ncols_plot = min(6, n)
    nrows_plot = (n + ncols_plot - 1) // ncols_plot

    fig, axes = plt.subplots(
        nrows_plot, ncols_plot,
        figsize=(ncols_plot * 2.8, nrows_plot * 2.6),
        squeeze=False,
    )
    axes_flat = axes.flatten()

    # Group colour: one colour per (dataset, column) combination for subplot borders
    groups = (col_meta["dataset"] + " / " + col_meta["column"]).unique().tolist()
    group_palette = dict(zip(groups, sns.color_palette("tab10", len(groups))))

    for idx, col_key in enumerate(all_cols):
        ax = axes_flat[idx]
        vals = ges_raw[col_key].values
        nonzero = vals[vals != 0]
        vmax = float(np.percentile(np.abs(nonzero), 95)) if len(nonzero) else 1.0
        sc_plot = ax.scatter(
            umap_xy[:, 0], umap_xy[:, 1],
            c=vals, s=1, cmap="RdBu_r", vmin=-vmax, vmax=vmax, rasterized=True,
        )
        condition_label = col_meta.loc[col_key, "condition"]
        ax.set_title(condition_label, fontsize=6.5, pad=2)
        ax.set_frame_on(False)
        ax.set_xticks([]); ax.set_yticks([])

        # Coloured border by (dataset, column) group
        group_key = col_meta.loc[col_key, "dataset"] + " / " + col_meta.loc[col_key, "column"]
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_edgecolor(group_palette[group_key])
            spine.set_linewidth(1.5)

        plt.colorbar(sc_plot, ax=ax, shrink=0.55, pad=0.02)

    for idx in range(n, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    # Legend for border colours
    from matplotlib.patches import Patch
    handles = [Patch(color=c, label=g) for g, c in group_palette.items()]
    fig.legend(
        handles=handles, title="Dataset / Column",
        loc="lower right", bbox_to_anchor=(1.0, 0.0),
        fontsize=7, frameon=False,
    )

    fig.suptitle(
        f"GES scores — all conditions ({dataset_name})",
        fontsize=11, y=1.01,
    )
    plt.tight_layout()
    out = fig_dir / "umap_all_conditions.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  🖼  umap_all_conditions.png")


def _make_cluster_heatmap(
    profiles: pd.DataFrame,
    col_meta: pd.DataFrame,
    fig_dir: Path,
    dataset_name: str,
    multi_dataset: bool,
) -> None:
    # Short display labels (strip dataset prefix if present)
    def short_label(col: str) -> str:
        parts = col.split("__")
        return parts[-1]  # just the condition name

    plot_profiles = profiles.rename(columns=short_label)

    # Build column annotation colours
    datasets_in_meta  = col_meta["dataset"].unique().tolist()
    columns_in_meta   = col_meta["column"].unique().tolist()
    dataset_palette   = dict(zip(datasets_in_meta, sns.color_palette("Set2",  len(datasets_in_meta))))
    col_grp_palette   = dict(zip(columns_in_meta,  sns.color_palette("Paired", len(columns_in_meta))))

    if multi_dataset:
        col_annotation = pd.DataFrame(
            {
                "Dataset": [dataset_palette[col_meta.loc[c, "dataset"]] for c in profiles.columns],
                "Column":  [col_grp_palette[col_meta.loc[c, "column"]]  for c in profiles.columns],
            },
            index=plot_profiles.columns,
        ).T
    else:
        col_annotation = pd.Series(
            [col_grp_palette[col_meta.loc[c, "column"]] for c in profiles.columns],
            index=plot_profiles.columns,
            name="Column",
        )

    fig_w = max(10, len(profiles.columns) * 0.45 + 3)
    fig_h = max(5,  len(profiles)         * 0.5  + 3)

    g = sns.clustermap(
        plot_profiles,
        col_colors=col_annotation,
        row_cluster=True,
        col_cluster=False,      # keep conditions in original dataset/column order
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
    g.ax_heatmap.set_title(f"Cluster expression profiles — {dataset_name}", pad=14)

    # Legend patches
    from matplotlib.patches import Patch
    legend_handles = []
    if multi_dataset:
        legend_handles += [Patch(color=c, label=d) for d, c in dataset_palette.items()]
        legend_handles += [Patch(facecolor="white", label="")]   # spacer
    legend_handles += [Patch(color=c, label=col) for col, c in col_grp_palette.items()]
    g.ax_col_colors.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(1.18, 2.0 if multi_dataset else 1.5),
        fontsize=7,
        title="Dataset / Column" if multi_dataset else "Column",
        frameon=False,
    )

    g.savefig(fig_dir / "cluster_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  🖼  cluster_heatmap.png")


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def run_gene_clustering(config_path: str) -> None:
    cfg, specs, multi_dataset = load_config(config_path)

    dataset_name      = cfg["dataset_name"]
    out_root          = cfg["output_folder"]
    min_conditions    = int(cfg.get("min_conditions_expressed", 2))
    ges_threshold     = cfg.get("ges_score_threshold", None)
    if ges_threshold is not None:
        ges_threshold = float(ges_threshold)
    scale_mode        = cfg.get("scale", "gene")
    pca_n_components  = int(cfg.get("pca_n_components", 50))
    n_neighbors       = int(cfg.get("umap_n_neighbors", 15))
    umap_min_dist     = float(cfg.get("umap_min_dist", 0.3))
    leiden_resolution = float(cfg.get("leiden_resolution", 0.5))

    date_str     = datetime.datetime.now().strftime("%Y%m%d")
    run_dir      = out_root / f"{dataset_name}_{date_str}"
    metadata_dir = run_dir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    with _log_to_file(metadata_dir / "pipeline_output.log"):
        print(f"📋 Log: {metadata_dir / 'pipeline_output.log'}")

        src = Path(cfg["_config_path"])
        dst = metadata_dir / src.name
        if src.resolve() != dst.resolve():
            shutil.copy2(src, dst)

        mode_str = "multi-dataset" if multi_dataset else "single-dataset"
        print(f"\n{'='*62}")
        print("  Gene Clustering Pipeline")
        print(f"{'='*62}")
        print(f"• Dataset name:        {dataset_name}  ({mode_str})")
        for s in specs:
            print(f"  · {s.label:<16} {s.ges_folder}")
        print(f"• Output:              {run_dir}")
        print(f"• Scale mode:          {scale_mode}")
        print(f"• Min conditions:      {min_conditions}")
        print(f"• GES threshold:       {ges_threshold if ges_threshold is not None else '(none)'}")
        print(f"• PCA components:      {pca_n_components}")
        print(f"• UMAP neighbors:      {n_neighbors}")
        print(f"• UMAP min_dist:       {umap_min_dist}")
        print(f"• Leiden resolution:   {leiden_resolution}")
        print(f"{'='*62}\n")

        print("Loading GES scores...")
        matrix_raw, n_present, col_meta = load_ges_matrix(specs, multi_dataset, ges_threshold)

        matrix_raw = filter_genes(matrix_raw, n_present, min_conditions)

        print(f"\nScaling (mode='{scale_mode}')...")
        matrix_scaled = scale_matrix(matrix_raw, scale_mode)

        adata = cluster_genes(
            matrix_scaled, pca_n_components, n_neighbors,
            umap_min_dist, leiden_resolution,
        )

        print("\nSaving results...")
        gene_clusters, profiles = save_results(adata, matrix_raw, run_dir)

        print("\nGenerating figures...")
        make_figures(
            adata, profiles, matrix_raw, col_meta,
            run_dir, dataset_name, specs, multi_dataset,
        )

        print(f"\n🎉 Gene clustering complete → {run_dir}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Cluster genes by GES expression profile (UMAP + Leiden)."
    )
    parser.add_argument("config", help="Path to the YAML config file.")
    args = parser.parse_args()
    run_gene_clustering(args.config)
