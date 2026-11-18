"""
cell_cycle_umap.py

UMAP coloured by the cell cycle annotation (CellCyclePhase): Non-cycling cells
in light grey underneath, cycling cells on top coloured by phase
(G1 blue, S orange, G2M green, Post-M red), drawn in random order so no phase
hides another where they overlap.

UMAP coordinates come from the file's obsm["X_umap"] or, if it has none, from
--umap-from (matched by cell ID; cells missing there are left out).

--regions keeps only cells of those Region values (e.g. the forebrain /
cortex: Telencephalon Diencephalon Forebrain), at their positions on the
full-dataset UMAP; the plot zooms to them.

--subsample plots a random subset of cells (a fraction, or a number of
cells), so dense regions are not solid colour; each cell has the same chance
of being kept, so the phase proportions do not change.

Outputs (results/additional_analyses/cell_cycle_umap/<label>/):
  umap_cell_cycle.png
  umap_cell_cycle_<phase>.png   with --per-phase: one phase highlighted at a time

Usage:
  python modules/cell_cycle_umap.py \
    --h5ad data/human_dev_without_week_5_cc_annotated_rescaled.h5ad \
    --umap-from data/full_adata_with_umap.h5ad --label human_dev
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

DEFAULT_OUTPUT_DIR = Path("results/additional_analyses/cell_cycle_umap")
UMAP_KEYS = ("X_umap", "UMAP")
NON_CYCLING = "Non-cycling"
PHASE_COLORS = {             # as in the reference figure
    "G1": "#1f1fe0",
    "S": "#f5a623",
    "G2M": "#1e9e3c",
    "Post-M": "#e8231f",
    NON_CYCLING: "#d9d9d9",
}
LEGEND_NAMES = {"G2M": "G2/M"}


def load(h5ad: str, umap_from: str | None, phase_col: str,
         region_col: str = "Region", regions: list[str] | None = None):
    """(xy, phase) for the cells that have UMAP coordinates (and, with
    `regions`, whose `region_col` is one of them)."""
    import anndata as ad
    a = ad.read_h5ad(h5ad, backed="r")
    try:
        if phase_col not in a.obs.columns:
            raise KeyError(f"obs column '{phase_col}' not found in {h5ad}")
        phase = (a.obs[phase_col].astype(str).str.replace("PostM", "Post-M")
                 .replace({"nan": NON_CYCLING}).to_numpy())
        keep = np.ones(a.n_obs, dtype=bool)
        if regions:
            if region_col not in a.obs.columns:
                raise KeyError(f"obs column '{region_col}' not found in {h5ad}")
            reg = a.obs[region_col].astype(str)
            keep = reg.isin(regions).to_numpy()
            print(f"  {region_col} in {regions}: {int(keep.sum()):,}/{a.n_obs:,} cells")
            if not keep.any():
                raise ValueError(f"No cells in {regions}. {region_col} values: "
                                 f"{sorted(reg.unique())}")
        phase = phase[keep]
        key = next((k for k in UMAP_KEYS if k in a.obsm), None)
        if key is not None:
            print(f"  UMAP: obsm['{key}'] of {h5ad}")
            return np.asarray(a.obsm[key])[:, :2][keep], phase
        names = a.obs_names[keep]
    finally:
        a.file.close()
    if not umap_from:
        raise KeyError(f"No obsm['X_umap'] in {h5ad}; pass --umap-from")
    src = ad.read_h5ad(umap_from, backed="r")
    try:
        key = next((k for k in UMAP_KEYS if k in src.obsm), None)
        if key is None:
            raise KeyError(f"No obsm['X_umap'] in {umap_from} either")
        pos = src.obs_names.get_indexer(names)
        found = pos >= 0
        print(f"  UMAP from {umap_from} (obsm['{key}']): "
              f"{int(found.sum()):,}/{len(names):,} cells matched by ID"
              + (f"; {int((~found).sum()):,} left out" if (~found).any() else ""))
        if not found.any():
            raise ValueError("No cells matched by ID")
        xy = np.asarray(src.obsm[key])[:, :2][pos[found]]
        return xy, phase[found]
    finally:
        src.file.close()


def _base(xy, title):
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(title, fontsize=13, fontweight="bold", loc="left")
    return fig, ax


def _legend(ax, phases, counts):
    handles = [Line2D([], [], marker="o", linestyle="", markersize=9,
                      color=PHASE_COLORS[p],
                      label=f"{LEGEND_NAMES.get(p, p)} ({counts.get(p, 0):,})")
               for p in phases]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0),
              frameon=False, fontsize=11)


def plot_all(xy, phase, title, path, point_size, seed):
    counts = pd.Series(phase).value_counts().to_dict()
    fig, ax = _base(xy, title)
    bg = phase == NON_CYCLING
    ax.scatter(xy[bg, 0], xy[bg, 1], s=point_size, c=PHASE_COLORS[NON_CYCLING],
               linewidths=0, rasterized=True)
    cyc = np.flatnonzero(np.isin(phase, [p for p in PHASE_COLORS if p != NON_CYCLING]))
    cyc = np.random.default_rng(seed).permutation(cyc)       # no phase drawn on top of all
    ax.scatter(xy[cyc, 0], xy[cyc, 1], s=point_size * 1.5,
               c=[PHASE_COLORS[p] for p in phase[cyc]], linewidths=0, rasterized=True)
    _legend(ax, list(PHASE_COLORS), counts)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_one(xy, phase, target, title, path, point_size):
    """All cells grey, `target` phase coloured on top."""
    fig, ax = _base(xy, title)
    m = phase == target
    ax.scatter(xy[~m, 0], xy[~m, 1], s=point_size, c=PHASE_COLORS[NON_CYCLING],
               linewidths=0, rasterized=True)
    ax.scatter(xy[m, 0], xy[m, 1], s=point_size * 1.5, c=PHASE_COLORS[target],
               linewidths=0, rasterized=True)
    _legend(ax, [target], {target: int(m.sum())})
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


def main():
    p = argparse.ArgumentParser(description="UMAP coloured by cell cycle phase.")
    p.add_argument("--h5ad", required=True, help="Cell-cycle-annotated h5ad.")
    p.add_argument("--label", required=True, help="Output subfolder name.")
    p.add_argument("--umap-from", default=None,
                   help="h5ad to take obsm['X_umap'] from when --h5ad has none.")
    p.add_argument("--phase-col", default="CellCyclePhase")
    p.add_argument("--region-col", default="Region")
    p.add_argument("--regions", nargs="+", default=None,
                   help="Only cells whose --region-col is one of these, e.g. "
                        "Telencephalon Diencephalon Forebrain. Their positions on the "
                        "full-dataset UMAP are kept (not a new UMAP).")
    p.add_argument("--title", default="Cell cycle")
    p.add_argument("--point-size", type=float, default=0.3)
    p.add_argument("--per-phase", action="store_true",
                   help="Also one UMAP per phase, that phase highlighted.")
    p.add_argument("--subsample", type=float, default=None,
                   help="Plot a random subset of cells: a fraction (e.g. 0.2) or a "
                        "number of cells (e.g. 300000). Every cell has the same chance, "
                        "so phase proportions are kept. Default: all cells.")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    args = p.parse_args()

    print(f"Loading (backed): {args.h5ad}")
    xy, phase = load(args.h5ad, args.umap_from, args.phase_col,
                     args.region_col, args.regions)
    print("  " + ", ".join(f"{k} {v:,}" for k, v in pd.Series(phase).value_counts().items()))
    if args.subsample:
        n = len(phase)
        k = int(round(args.subsample * n)) if args.subsample <= 1 else int(args.subsample)
        if k < n:
            keep = np.sort(np.random.default_rng(args.seed).choice(n, k, replace=False))
            xy, phase = xy[keep], phase[keep]
            print(f"  Subsample: {k:,}/{n:,} random cells plotted (legend counts are of these)")
    out_dir = Path(args.output_dir) / args.label
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_all(xy, phase, args.title, out_dir / "umap_cell_cycle.png", args.point_size, args.seed)
    if args.per_phase:
        for ph in [p for p in PHASE_COLORS if p != NON_CYCLING]:
            plot_one(xy, phase, ph, f"{args.title}: {LEGEND_NAMES.get(ph, ph)}",
                     out_dir / f"umap_cell_cycle_{ph.replace('-', '')}.png", args.point_size)
    print(f"Outputs: {out_dir}")


if __name__ == "__main__":
    main()
