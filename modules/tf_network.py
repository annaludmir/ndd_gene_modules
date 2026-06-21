"""
tf_network.py

Transcription Factor regulatory network inference and gene-set query.

Pipeline (each step saves its output and is skipped on re-runs if the file exists):
  1. Aggregation   — SEACells meta-cell aggregation (or Leiden pseudo-bulk fallback)
  2. GRN           — GRNBoost2: infer TF → target importance scores
  3. cisTarget     — pySCENIC: prune adjacencies to motif-supported regulons
  4. AUCell        — compute per-cell regulon activity AUC scores
  5. Query / Plot  — filter network for a gene set, produce hub-and-spoke figure

Output structure:
  results/tf_network/{dataset_name}_{YYYYMMDD}/
    seacells/seacells_aggregated.h5ad
    grn/adjacencies.tsv
    ctx/regulons.csv
    aucell/auc_matrix.csv
    data/tf_targets_{gene_list_stem}.csv
    figures/tf_network_{gene_list_stem}.png
    metadata/pipeline_output.log

Usage:
  # Run full pipeline:
  python modules/tf_network.py config.yaml

  # Run full pipeline + query a gene set:
  python modules/tf_network.py config.yaml --gene-list gene_lists/my_genes.csv

  # Query only (pipeline results already exist):
  python modules/tf_network.py config.yaml --gene-list gene_lists/my_genes.csv --query-only
"""

import argparse
import contextlib
import datetime
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch
import numpy as np
import pandas as pd
import scanpy as sc
import yaml


# ---------------------------------------------------------------------------
# Colormap matching the reference figure (dark maroon → orange → yellow)
# ---------------------------------------------------------------------------

_TF_CMAP = LinearSegmentedColormap.from_list(
    "tf_score",
    ["#3d0000", "#8B0000", "#CC3300", "#FF6600", "#FFD700"],
)


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

def load_config(config_path: str) -> dict:
    config_path = Path(config_path).resolve()
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    root = Path(cfg["ndd_gene_modules_folder_root"]).resolve()

    def resolve(p):
        if p is None:
            return None
        p = Path(str(p))
        return p if p.is_absolute() else (root / p).resolve()

    cfg["output_folder"] = resolve(cfg.get("output_folder", "results/tf_network"))
    cfg["h5ad_path"]     = resolve(cfg["h5ad_path"])
    cfg["_config_path"]  = config_path
    cfg["_root"]         = root

    scenic = cfg.get("scenic", {})
    scenic["tf_list"]            = resolve(scenic.get("tf_list"))
    scenic["motif_annotations"]  = resolve(scenic.get("motif_annotations"))
    dbs = scenic.get("rankings_db", [])
    scenic["rankings_db"] = [resolve(d) for d in (dbs if isinstance(dbs, list) else [dbs])]
    cfg["scenic"] = scenic

    return cfg


# ---------------------------------------------------------------------------
# Gene set / TF list loaders
# ---------------------------------------------------------------------------

def load_gene_set(path: str | Path) -> set[str]:
    df = pd.read_csv(path)
    col = "gene" if "gene" in df.columns else df.columns[0]
    genes = {str(g).strip() for g in df[col].dropna() if str(g).strip()}
    print(f"Gene set: {len(genes)} genes from {Path(path).name}")
    return genes


def _load_tf_list(tf_list_path: Path) -> list[str]:
    with open(tf_list_path) as f:
        tfs = [line.strip() for line in f if line.strip()]
    print(f"TF list: {len(tfs)} transcription factors")
    return tfs


# ---------------------------------------------------------------------------
# Step 1 — Aggregation (SEACells or Leiden pseudo-bulk)
# ---------------------------------------------------------------------------

def _aggregate_seacells(adata: sc.AnnData, cfg: dict) -> sc.AnnData:
    """SEACells meta-cell aggregation."""
    try:
        import SEACells  # noqa: F401
    except ImportError:
        raise ImportError(
            "SEACells not installed. Install with: pip install SEACells\n"
            "Or set aggregation.method: 'leiden' in the config to use pseudo-bulk."
        )

    sc_cfg = cfg.get("seacells", {})
    n_cells_per_metacell = int(sc_cfg.get("cells_per_metacell", 75))
    n_metacells = max(10, adata.n_obs // n_cells_per_metacell)
    n_eigs = int(sc_cfg.get("n_eigs", 10))

    print(f"  SEACells: {adata.n_obs:,} cells → {n_metacells} meta-cells")

    # Ensure PCA is available (or use scVI embedding if present)
    use_rep = "X_scVI" if "X_scVI" in adata.obsm else "X_pca"
    if use_rep not in adata.obsm:
        print("  Computing PCA for SEACells kernel...")
        sc.pp.pca(adata, n_comps=50)
        use_rep = "X_pca"

    model = SEACells.core.SEACells(
        adata,
        build_kernel_on=use_rep,
        n_SEACells=n_metacells,
        n_waypoint_eigs=n_eigs,
        convergence_epsilon=float(sc_cfg.get("convergence_threshold", 0.01125)),
    )
    model.construct_kernel_matrix()
    model.initialize_archetypes()
    model.fit(
        min_iter=int(sc_cfg.get("min_iter", 10)),
        max_iter=int(sc_cfg.get("max_iter", 100)),
    )

    # Aggregate: sum log-normalised expression within each meta-cell
    meta_ad = SEACells.core.summarize_by_SEACell(adata, SEACells_label="SEACell", summarize_layer="X")
    print(f"  SEACells done: {meta_ad.n_obs} meta-cells × {meta_ad.n_vars} genes")
    return meta_ad


def _aggregate_leiden(adata: sc.AnnData, cfg: dict) -> sc.AnnData:
    """Pseudo-bulk: sum log-normalised expression per Leiden cluster."""
    leiden_col = cfg.get("aggregation", {}).get("leiden_column", "leiden")
    if leiden_col not in adata.obs.columns:
        print(f"  Column '{leiden_col}' not found — running Leiden clustering...")
        sc.pp.pca(adata, n_comps=50)
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata, resolution=0.5, flavor="igraph",
                     n_iterations=2, directed=False)
        leiden_col = "leiden"

    clusters = adata.obs[leiden_col].astype(str)
    unique_clusters = clusters.unique()
    print(f"  Pseudo-bulk: {len(unique_clusters)} Leiden clusters as meta-cells")

    rows = {}
    X = adata.X
    import scipy.sparse as sp
    if sp.issparse(X):
        X = X.toarray()
    for cl in unique_clusters:
        mask = (clusters == cl).values
        rows[cl] = X[mask].sum(axis=0)

    agg_df  = pd.DataFrame(rows, index=adata.var_names).T
    meta_ad = sc.AnnData(X=agg_df.values,
                         obs=pd.DataFrame(index=agg_df.index),
                         var=pd.DataFrame(index=adata.var_names))
    print(f"  Pseudo-bulk done: {meta_ad.n_obs} meta-cells × {meta_ad.n_vars} genes")
    return meta_ad


def run_aggregation(cfg: dict, out_dir: Path) -> sc.AnnData:
    out_file = out_dir / "seacells" / "seacells_aggregated.h5ad"
    if out_file.exists():
        print(f"[SKIP] Aggregation — loading existing: {out_file.name}")
        return sc.read_h5ad(out_file)

    print("\n[Step 1] Aggregation")
    out_file.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(cfg["h5ad_path"])
    print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    chemistry = cfg.get("chemistry")
    if chemistry and "Chemistry" in adata.obs.columns:
        adata = adata[adata.obs["Chemistry"] == chemistry].copy()
        print(f"  Filtered to chemistry={chemistry}: {adata.n_obs:,} cells")

    # Normalise if not already done (check if max value looks raw)
    import scipy.sparse as sp
    X_check = adata.X[:100].toarray() if sp.issparse(adata.X) else adata.X[:100]
    if X_check.max() > 50:
        print("  Normalising (log1p)...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # Optional: protein-coding genes only
    coding_path = cfg.get("protein_coding_genes")
    if coding_path:
        coding = set(pd.read_csv(coding_path, header=None)[0].astype(str))
        keep = [g for g in adata.var_names if g in coding]
        adata = adata[:, keep].copy()
        print(f"  Protein-coding filter: {len(keep):,} genes retained")

    method = cfg.get("aggregation", {}).get("method", "leiden")
    if method == "seacells":
        meta_ad = _aggregate_seacells(adata, cfg)
    else:
        meta_ad = _aggregate_leiden(adata, cfg)

    meta_ad.write_h5ad(out_file)
    print(f"  Saved: {out_file}")
    return meta_ad


# ---------------------------------------------------------------------------
# Step 2 — GRN inference (GENIE3 / ExtraTrees, no Dask dependency)
# ---------------------------------------------------------------------------
# arboreto/GRNBoost2 relies on a Dask distributed cluster that is broken in
# this environment (distributed API version mismatch).  We implement the same
# logic directly: for each target gene fit an ExtraTreesRegressor using TF
# expression values as features; feature importances become TF→target scores.
# This is equivalent to GENIE3 and produces the same TSV output format.
# ---------------------------------------------------------------------------

def _fit_one_target(target_name, y, X_tfs, tf_cols, params):
    """Fit one ExtraTrees model for a single target gene; return importance rows."""
    if float(np.std(y)) < 1e-6:
        return []
    from sklearn.ensemble import ExtraTreesRegressor
    model = ExtraTreesRegressor(**params)
    try:
        model.fit(X_tfs, y)
    except Exception:
        return []
    return [
        {"TF": tf, "target": target_name, "importance": float(imp)}
        for tf, imp in zip(tf_cols, model.feature_importances_)
        if imp > 0 and tf != target_name
    ]


def run_grn(meta_ad: sc.AnnData, cfg: dict, out_dir: Path) -> pd.DataFrame:
    out_file = out_dir / "grn" / "adjacencies.tsv"
    if out_file.exists():
        print(f"[SKIP] GRN — loading existing: {out_file.name}")
        return pd.read_csv(out_file, sep="\t")

    print("\n[Step 2] GRN inference (ExtraTrees / GENIE3-style, joblib parallel)")
    out_file.parent.mkdir(parents=True, exist_ok=True)

    scenic_cfg   = cfg.get("scenic", {})
    tf_list_path = scenic_cfg.get("tf_list")
    gene_names   = list(meta_ad.var_names)

    if tf_list_path and Path(str(tf_list_path)).exists():
        tf_names = _load_tf_list(Path(str(tf_list_path)))
        tf_names = [t for t in tf_names if t in set(gene_names)]
        print(f"  TFs found in expression matrix: {len(tf_names):,}")
    else:
        print("  Warning: no TF list — treating all genes as potential TFs")
        tf_names = gene_names

    import scipy.sparse as sp
    X = meta_ad.X
    if sp.issparse(X):
        X = X.toarray()
    X = X.astype(np.float32)

    tf_set  = set(tf_names)
    tf_cols = [g for g in gene_names if g in tf_set]
    tf_idx  = [gene_names.index(g) for g in tf_cols]
    X_tfs   = X[:, tf_idx]

    n_workers     = int(scenic_cfg.get("n_workers", 4))
    n_estimators  = int(scenic_cfg.get("n_estimators", 500))
    print(f"  {meta_ad.n_obs} meta-cells × {len(gene_names):,} target genes, "
          f"{len(tf_cols):,} TFs, n_estimators={n_estimators}, n_workers={n_workers}")

    params = dict(n_estimators=n_estimators, max_features="sqrt",
                  random_state=42, n_jobs=1)

    from joblib import Parallel, delayed
    results = Parallel(n_jobs=n_workers, verbose=5)(
        delayed(_fit_one_target)(
            gene_names[i],
            X[:, i],
            X_tfs,
            tf_cols,
            params,
        )
        for i in range(len(gene_names))
    )

    all_pairs = [row for pairs in results for row in pairs]
    adj = (pd.DataFrame(all_pairs)
           .sort_values("importance", ascending=False)
           .reset_index(drop=True))

    adj.to_csv(out_file, sep="\t", index=False)
    print(f"  Saved: {out_file.name}  ({len(adj):,} TF-target pairs)")
    return adj


# ---------------------------------------------------------------------------
# Step 3 — cisTarget regulon prediction (optional)
# ---------------------------------------------------------------------------

def run_ctx(adj: pd.DataFrame, meta_ad: sc.AnnData, cfg: dict, out_dir: Path) -> pd.DataFrame | None:
    out_file = out_dir / "ctx" / "regulons.csv"
    if out_file.exists():
        print(f"[SKIP] cisTarget — loading existing: {out_file.name}")
        return pd.read_csv(out_file)

    scenic_cfg = cfg.get("scenic", {})
    db_paths   = [p for p in scenic_cfg.get("rankings_db", []) if p and Path(p).exists()]
    motif_file = scenic_cfg.get("motif_annotations")

    if not db_paths:
        print("\n[Skip] cisTarget — no rankings_db files found in config.")
        print("  To enable: provide rankings_db paths pointing to .feather files from:")
        print("  https://resources.aertslab.org/cistarget/")
        return None

    if not motif_file or not Path(motif_file).exists():
        print("\n[Skip] cisTarget — motif_annotations file not found.")
        return None

    print("\n[Step 3] cisTarget regulon prediction")
    out_file.parent.mkdir(parents=True, exist_ok=True)

    # pyscenic uses removed NumPy aliases (np.object, np.bool, etc.) — restore them
    # before importing so the module loads correctly on NumPy >= 1.24
    import numpy as _np
    for _alias, _builtin in [("object", object), ("bool", bool), ("int", int),
                              ("float", float), ("complex", complex), ("str", str)]:
        if not hasattr(_np, _alias):
            setattr(_np, _alias, _builtin)

    try:
        from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
        from pyscenic.utils import modules_from_adjacencies
        from pyscenic.prune import prune2df, df2regulons
    except ImportError:
        print("\n[Skip] cisTarget — pyscenic / ctxcore not installed (pip install pyscenic).")
        return None

    import scipy.sparse as sp
    X = meta_ad.X
    if sp.issparse(X):
        X = X.toarray()
    ex_mtx = pd.DataFrame(X, columns=meta_ad.var_names)

    dbs = [RankingDatabase(fname=str(p), name=Path(p).stem) for p in db_paths]
    print(f"  Databases: {[db.name for db in dbs]}")

    modules = list(modules_from_adjacencies(adj, ex_mtx))
    print(f"  Modules from adjacencies: {len(modules)}")

    # Newer dask requires a list (with len()), but pyscenic passes a generator.
    # Patch from_delayed on both dask.dataframe and pyscenic.prune's local binding.
    import dask.dataframe as _dask_df
    import pyscenic.prune as _pyscenic_prune
    _orig_from_delayed = _dask_df.from_delayed
    def _compat_from_delayed(dfs, *args, **kwargs):
        if not hasattr(dfs, "__len__"):
            dfs = list(dfs)
        return _orig_from_delayed(dfs, *args, **kwargs)
    _dask_df.from_delayed = _compat_from_delayed
    if hasattr(_pyscenic_prune, "from_delayed"):
        _pyscenic_prune.from_delayed = _compat_from_delayed

    df = prune2df(
        dbs, modules, str(motif_file),
        num_workers=int(scenic_cfg.get("n_workers", 4)),
    )

    _dask_df.from_delayed = _orig_from_delayed  # restore after use
    regulons = df2regulons(df)
    print(f"  Regulons (motif-supported): {len(regulons)}")

    # Save as CSV: TF, targets (semicolon-separated), n_targets
    rows = [{"TF": r.name, "targets": ";".join(sorted(r.genes)), "n_targets": len(r.genes)}
            for r in regulons]
    reg_df = pd.DataFrame(rows)
    reg_df.to_csv(out_file, index=False)
    print(f"  Saved: {out_file.name}")
    return reg_df


# ---------------------------------------------------------------------------
# Step 4 — AUCell
# ---------------------------------------------------------------------------

def run_aucell(meta_ad: sc.AnnData, regulons_csv: pd.DataFrame | None,
               adj: pd.DataFrame, cfg: dict, out_dir: Path) -> pd.DataFrame | None:
    out_file = out_dir / "aucell" / "auc_matrix.csv"
    if out_file.exists():
        print(f"[SKIP] AUCell — loading existing: {out_file.name}")
        return pd.read_csv(out_file, index_col=0)

    if regulons_csv is None:
        print("\n[Skip] AUCell — requires cisTarget regulons (step 3).")
        return None

    print("\n[Step 4] AUCell — regulon activity scores")
    out_file.parent.mkdir(parents=True, exist_ok=True)

    try:
        from pyscenic.aucell import aucell
        from pyscenic.utils import modules_from_adjacencies
        from ctxcore.genesig import GeneSignature
    except ImportError:
        raise ImportError("pyscenic not installed.  Run: pip install pyscenic")

    import scipy.sparse as sp
    X = meta_ad.X
    if sp.issparse(X):
        X = X.toarray()
    ex_mtx = pd.DataFrame(X, index=meta_ad.obs_names, columns=meta_ad.var_names)

    # Rebuild regulon gene signatures from saved CSV
    regulons = [
        GeneSignature(name=row["TF"], gene2weight={g: 1.0 for g in row["targets"].split(";")})
        for _, row in regulons_csv.iterrows()
    ]

    auc_cfg = cfg.get("aucell", {})
    auc_mtx = aucell(
        ex_mtx,
        regulons,
        auc_threshold=float(auc_cfg.get("auc_threshold", 0.05)),
        num_workers=int(auc_cfg.get("num_workers", 4)),
    )
    auc_mtx.to_csv(out_file)
    print(f"  Saved: {out_file.name}  ({auc_mtx.shape[0]} cells × {auc_mtx.shape[1]} TFs)")
    return auc_mtx


# ---------------------------------------------------------------------------
# Step 5 — Gene set query
# ---------------------------------------------------------------------------

def query_gene_set(
    gene_set: set[str],
    adj: pd.DataFrame,
    regulons_csv: pd.DataFrame | None,
    cfg: dict,
) -> pd.DataFrame:
    """
    Filter TF-target network to keep only edges where target is in gene_set.

    Returns a DataFrame with columns:
      tf, target, importance, normalized_score
    where normalized_score ∈ [0,1] within each TF (importance / max_importance_for_TF).
    """
    q_cfg = cfg.get("query", {})

    if regulons_csv is not None:
        # Use motif-supported regulons (higher confidence)
        rows = []
        for _, row in regulons_csv.iterrows():
            targets = set(row["targets"].split(";")) & gene_set
            if not targets:
                continue
            # Look up importance scores from adjacencies
            tf_adj = adj[adj["TF"] == row["TF"]]
            for tgt in targets:
                imp_row = tf_adj[tf_adj["target"] == tgt]
                imp = float(imp_row["importance"].values[0]) if not imp_row.empty else 0.0
                rows.append({"tf": row["TF"], "target": tgt, "importance": imp})
        result = pd.DataFrame(rows)
    else:
        # Fall back to raw adjacencies
        col_map = {c.lower(): c for c in adj.columns}
        tf_col  = col_map.get("tf",     col_map.get("transcription_factor", adj.columns[0]))
        tgt_col = col_map.get("target", adj.columns[1])
        imp_col = col_map.get("importance", adj.columns[2])
        adj_filt = adj[adj[tgt_col].isin(gene_set)].copy()
        adj_filt = adj_filt.rename(columns={tf_col: "tf", tgt_col: "target", imp_col: "importance"})
        result   = adj_filt[["tf", "target", "importance"]].copy()

    if result.empty:
        print("  Warning: no TF-target pairs found for the provided gene set.")
        return result

    # Normalize importance to [0,1] within each TF
    result["normalized_score"] = result.groupby("tf")["importance"].transform(
        lambda x: (x - x.min()) / (x.max() - x.min() + 1e-9)
    )

    # Filter by minimum targets in gene set
    min_targets = int(q_cfg.get("min_targets", 3))
    counts = result.groupby("tf")["target"].count()
    valid_tfs = counts[counts >= min_targets].index
    result = result[result["tf"].isin(valid_tfs)]

    result = result.sort_values(["tf", "normalized_score"], ascending=[True, False])
    print(f"  Gene set query: {len(valid_tfs)} TFs with ≥{min_targets} targets in gene set")
    print(f"  Total TF-target pairs: {len(result):,}")
    return result


# ---------------------------------------------------------------------------
# Step 6 — Hub-and-spoke visualization
# ---------------------------------------------------------------------------

def _draw_hub_cluster(
    ax,
    tf: str,
    targets: list[str],
    scores: list[float],
    center: np.ndarray,
    radius: float,
    hub_radius: float,
    node_radius: float,
    cmap,
):
    """Draw one TF hub + target satellites onto ax, centred at `center`."""
    n = len(targets)
    if n == 0:
        return

    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    target_pos = center + radius * np.stack([np.cos(angles), np.sin(angles)], axis=1)

    for i, (pos, score) in enumerate(zip(target_pos, scores)):
        # Edge: TF → target
        direction = pos - center
        dist      = np.linalg.norm(direction)
        unit      = direction / dist

        # Start just outside TF circle, end just outside target circle
        start = center + unit * hub_radius
        end   = pos    - unit * node_radius

        ax.annotate(
            "",
            xy=end, xytext=start,
            arrowprops=dict(
                arrowstyle="-|>",
                color="#888888",
                lw=0.7,
                mutation_scale=8,
            ),
            zorder=2,
        )

        # Target node
        color = cmap(score)
        circle = plt.Circle(pos, node_radius, color=color, ec="none", zorder=3)
        ax.add_patch(circle)

        # Gene label — push outward
        label_pos = center + (radius + node_radius + 0.05) * (pos - center) / np.linalg.norm(pos - center)
        ha = "left" if pos[0] >= center[0] else "right"
        va = "center"
        ax.text(
            label_pos[0], label_pos[1],
            targets[i],
            ha=ha, va=va,
            fontsize=7,
            fontweight="bold",
            color="#222222",
            zorder=5,
        )

    # TF hub node (on top)
    hub = plt.Circle(center, hub_radius, color="white", ec="#333333", lw=1.8, zorder=4)
    ax.add_patch(hub)
    ax.text(
        center[0], center[1], tf,
        ha="center", va="center",
        fontsize=8, fontweight="bold",
        color="#111111", zorder=6,
    )


def plot_tf_network(
    query_df: pd.DataFrame,
    fig_path: Path,
    cfg: dict,
    title: str = "TF Regulatory Networks",
) -> None:
    """
    Hub-and-spoke figure matching the reference paper's style.

    Layout: top N TFs arranged in a grid, each with its gene-set target satellites.
    Node colour encodes normalized_score (dark maroon=0 → yellow=1).
    """
    if query_df.empty:
        print("  No data to plot.")
        return

    q_cfg   = cfg.get("query", {})
    top_n   = int(q_cfg.get("top_n_tfs", 8))
    min_tgt = int(q_cfg.get("min_targets", 3))

    # Select top TFs by number of targets (then mean importance)
    tf_stats = (
        query_df.groupby("tf")
        .agg(n_targets=("target", "count"), mean_score=("normalized_score", "mean"))
        .sort_values(["n_targets", "mean_score"], ascending=False)
        .head(top_n)
    )
    top_tfs = tf_stats.index.tolist()

    n_tfs = len(top_tfs)
    if n_tfs == 0:
        print("  No TFs to plot.")
        return

    ncols = min(3, n_tfs)
    nrows = (n_tfs + ncols - 1) // ncols

    # Compute per-cluster sizes to set axis limits
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(ncols * 5.5, nrows * 5.5),
        squeeze=False,
    )
    axes_flat = axes.flatten()

    cmap = _TF_CMAP

    for idx, tf in enumerate(top_tfs):
        ax  = axes_flat[idx]
        sub = query_df[query_df["tf"] == tf].sort_values("normalized_score", ascending=False)
        targets = sub["target"].tolist()
        scores  = sub["normalized_score"].tolist()
        n       = len(targets)

        # Scale radius with number of targets
        radius      = max(0.9, 0.15 * n)
        hub_radius  = min(0.22, radius * 0.22)
        node_radius = min(0.09, radius * 0.09)
        center      = np.array([0.0, 0.0])

        _draw_hub_cluster(
            ax, tf, targets, scores,
            center=center,
            radius=radius,
            hub_radius=hub_radius,
            node_radius=node_radius,
            cmap=cmap,
        )

        margin = radius + node_radius + 0.4
        ax.set_xlim(-margin, margin)
        ax.set_ylim(-margin, margin)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"{tf}  ({n} targets)", fontsize=8, pad=4)

    for idx in range(n_tfs, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    # Shared colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes_flat[:n_tfs].tolist(),
                        shrink=0.35, pad=0.04, aspect=20)
    cbar.set_label("Normalized\nPrediction score", fontsize=9)
    cbar.set_ticks([0, 0.5, 1])

    fig.suptitle(title, fontsize=12, y=1.01)
    plt.tight_layout()
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fig_path.name}")


# ---------------------------------------------------------------------------
# Connected-network visualization (shared targets appear once)
# ---------------------------------------------------------------------------

def plot_tf_network_connected(
    query_df: pd.DataFrame,
    fig_path: Path,
    cfg: dict,
    title: str = "TF Regulatory Networks",
) -> None:
    """
    Single connected graph where targets shared by multiple TFs appear once.
    Requires: pip install networkx (usually already present via scanpy).

    plot_style: "connected"  in config query section to enable.
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("networkx required: pip install networkx")

    if query_df.empty:
        print("  No data to plot.")
        return

    q_cfg = cfg.get("query", {})
    top_n = int(q_cfg.get("top_n_tfs", 8))

    tf_stats = (
        query_df.groupby("tf")
        .agg(n_targets=("target", "count"), mean_score=("normalized_score", "mean"))
        .sort_values(["n_targets", "mean_score"], ascending=False)
        .head(top_n)
    )
    top_tfs = tf_stats.index.tolist()
    plot_df = query_df[query_df["tf"].isin(top_tfs)].copy()

    # Build directed graph; for targets shared by multiple TFs use max score
    G = nx.DiGraph()
    for tf in top_tfs:
        G.add_node(tf, node_type="tf")

    target_max_score: dict[str, float] = {}
    for _, row in plot_df.iterrows():
        tgt   = row["target"]
        score = float(row["normalized_score"])
        if tgt not in target_max_score or score > target_max_score[tgt]:
            target_max_score[tgt] = score
        if tgt not in G:
            G.add_node(tgt, node_type="gene")
        G.add_edge(row["tf"], tgt, weight=score)

    # Spring layout — increase k so nodes are well-separated
    pos = nx.spring_layout(G, seed=42, k=3.0, iterations=150)

    fig, ax = plt.subplots(figsize=(14, 11))
    ax.set_aspect("equal")
    ax.axis("off")

    tf_nodes   = [n for n in G.nodes if G.nodes[n]["node_type"] == "tf"]
    gene_nodes = [n for n in G.nodes if G.nodes[n]["node_type"] == "gene"]

    # Node radii (spring layout coords are roughly in [-1, 1])
    tf_r   = 0.07
    gene_r = 0.035

    # ── Edges ────────────────────────────────────────────────────────────────
    for u, v in G.edges():
        pu = np.array(pos[u])
        pv = np.array(pos[v])
        d  = pv - pu
        dist = np.linalg.norm(d)
        if dist < 1e-6:
            continue
        unit  = d / dist
        r_u   = tf_r   if G.nodes[u]["node_type"] == "tf"   else gene_r
        r_v   = tf_r   if G.nodes[v]["node_type"] == "tf"   else gene_r
        start = pu + unit * r_u
        end   = pv - unit * r_v
        ax.annotate(
            "",
            xy=end, xytext=start,
            arrowprops=dict(arrowstyle="-|>", color="#888888", lw=0.8, mutation_scale=9),
            zorder=2,
        )

    # ── Gene (target) nodes ──────────────────────────────────────────────────
    for gene in gene_nodes:
        p     = np.array(pos[gene])
        score = target_max_score.get(gene, 0.0)
        color = _TF_CMAP(score)
        ax.add_patch(plt.Circle(p, gene_r, color=color, ec="none", zorder=3))

        # Push label away from the TF centroid
        tf_preds = [u for u in G.predecessors(gene) if G.nodes[u]["node_type"] == "tf"]
        if tf_preds:
            centroid  = np.mean([pos[t] for t in tf_preds], axis=0)
            push      = p - centroid
            push_norm = np.linalg.norm(push)
            push_dir  = push / push_norm if push_norm > 1e-6 else np.array([1.0, 0.0])
        else:
            push_dir = np.array([1.0, 0.0])

        lp = p + push_dir * (gene_r + 0.028)
        ha = "left" if push_dir[0] >= 0 else "right"
        ax.text(lp[0], lp[1], gene, ha=ha, va="center",
                fontsize=7, fontweight="bold", color="#222222", zorder=5)

    # ── TF hub nodes ─────────────────────────────────────────────────────────
    for tf in tf_nodes:
        p = np.array(pos[tf])
        ax.add_patch(plt.Circle(p, tf_r, color="white", ec="#333333", lw=2.0, zorder=4))
        ax.text(p[0], p[1], tf, ha="center", va="center",
                fontsize=9, fontweight="bold", color="#111111", zorder=6)

    # ── Colorbar ─────────────────────────────────────────────────────────────
    sm = plt.cm.ScalarMappable(cmap=_TF_CMAP, norm=plt.Normalize(0, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.28, pad=0.02, aspect=18, anchor=(1.0, 0.5))
    cbar.set_label("Normalized\nPrediction score", fontsize=9)
    cbar.set_ticks([0, 0.5, 1])

    all_pos = np.array(list(pos.values()))
    margin  = 0.30
    ax.set_xlim(all_pos[:, 0].min() - margin, all_pos[:, 0].max() + margin)
    ax.set_ylim(all_pos[:, 1].min() - margin, all_pos[:, 1].max() + margin)

    fig.suptitle(title, fontsize=12)
    plt.tight_layout()
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fig_path.name}")


# ---------------------------------------------------------------------------
# Pipeline orchestration
# ---------------------------------------------------------------------------

def run_tf_network(
    config_path: str,
    gene_list_path: str | None = None,
    query_only: bool = False,
) -> None:
    cfg = load_config(config_path)

    dataset_name = cfg["dataset_name"]
    date_str     = datetime.datetime.now().strftime("%Y%m%d")
    out_root     = cfg["output_folder"]
    run_dir      = out_root / f"{dataset_name}_{date_str}"
    run_dir.mkdir(parents=True, exist_ok=True)

    shutil.copy2(cfg["_config_path"], run_dir / Path(config_path).name)

    meta_dir = run_dir / "metadata"
    meta_dir.mkdir(exist_ok=True)

    with _log_to_file(meta_dir / "pipeline_output.log"):
        print(f"\n{'='*62}")
        print(f"  TF Network — {dataset_name}")
        print(f"  Output: {run_dir}")
        print(f"{'='*62}\n")

        # ── Pipeline steps ───────────────────────────────────────────────────
        if not query_only:
            meta_ad      = run_aggregation(cfg, run_dir)
            adj          = run_grn(meta_ad, cfg, run_dir)
            regulons_csv = run_ctx(adj, meta_ad, cfg, run_dir)
            run_aucell(meta_ad, regulons_csv, adj, cfg, run_dir)
        else:
            # Load existing results
            adj_path = run_dir / "grn" / "adjacencies.tsv"
            reg_path = run_dir / "ctx" / "regulons.csv"
            if not adj_path.exists():
                raise FileNotFoundError(
                    f"adjacencies.tsv not found at {adj_path}.\n"
                    "Run the full pipeline first (without --query-only)."
                )
            print("[query-only] Loading precomputed adjacencies...")
            adj          = pd.read_csv(adj_path, sep="\t")
            regulons_csv = pd.read_csv(reg_path) if reg_path.exists() else None

        # ── Gene set query + plot ────────────────────────────────────────────
        if gene_list_path:
            gene_set      = load_gene_set(gene_list_path)
            gene_list_stem = Path(gene_list_path).stem

            print(f"\n[Query] Filtering network for gene set: {gene_list_stem}")
            query_df = query_gene_set(gene_set, adj, regulons_csv, cfg)

            if not query_df.empty:
                data_dir = run_dir / "data"
                data_dir.mkdir(exist_ok=True)
                csv_path = data_dir / f"tf_targets_{gene_list_stem}.csv"
                query_df.to_csv(csv_path, index=False)
                print(f"  Saved: {csv_path.name}")

                fig_dir = run_dir / "figures"
                fig_dir.mkdir(exist_ok=True)
                plot_style = cfg.get("query", {}).get("plot_style", "hub_spoke")
                fig_path   = fig_dir / f"tf_network_{gene_list_stem}.png"
                plot_fn    = (plot_tf_network_connected
                              if plot_style == "connected"
                              else plot_tf_network)
                plot_fn(
                    query_df,
                    fig_path,
                    cfg,
                    title=f"TF Regulatory Networks — {gene_list_stem}",
                )

        print(f"\nDone → {run_dir}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="TF regulatory network: SEACells → GRNBoost2 → cisTarget → AUCell + hub-and-spoke plot."
    )
    parser.add_argument("config", help="YAML config file.")
    parser.add_argument(
        "--gene-list",
        default=None,
        metavar="CSV",
        help="Gene list CSV to query. Shows which TFs regulate these genes.",
    )
    parser.add_argument(
        "--query-only",
        action="store_true",
        help="Skip pipeline steps; load existing adjacencies.tsv and query only.",
    )
    args = parser.parse_args()
    run_tf_network(args.config, args.gene_list, args.query_only)
