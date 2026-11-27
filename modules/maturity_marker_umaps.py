"""
maturity_marker_umaps.py

Where are the neuron maturity markers expressed? For each marker of the
NeuronMaturity annotation (neuron_type_annotation.py, MATURITY_MARKER_SETS:
default and/or braun2023), one figure on the human_dev UMAP:

  left    UMAP coloured by region (--region-col), the same colours in every
          figure
  middle  UMAP coloured by the gene's expression (log1p CP10k), cells with
          higher expression drawn on top, cells with 0 counts light grey
  right   per region: % of cells with >= 1 count, and mean log1p CP10k
          (so a gene seen only in the cortex stands out)

--cell-class Neuron restricts everything to neurons (the cells the maturity
annotation was made on), still at their positions on the full UMAP.

UMAP: obsm["X_umap"] of --h5ad or, if it has none, of --umap-from (matched by
cell ID; cells missing there are left out of the UMAP panels, not of the
per-region numbers).

Outputs (results/additional_analyses/maturity_marker_umaps/<label>/):
  <set>/<immature|mature>_<gene>.png
  expression_by_region.csv   gene, set, group, region: cells, % detected,
                             mean log1p CP10k

Usage:
  python modules/maturity_marker_umaps.py \
    --h5ad data/human_dev_without_week_5_cc_annotated_neuron_type_braun2023.h5ad \
    --umap-from data/full_adata_with_umap.h5ad --label human_dev
"""

import argparse
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from neuron_type_annotation import MATURITY_MARKER_SETS  # noqa: E402

DEFAULT_OUTPUT_DIR = Path("results/additional_analyses/maturity_marker_umaps")
UMAP_KEYS = ("X_umap", "UMAP")
ZERO_COLOR = "#e3e3e3"
EXPR_CMAP = "viridis"


# ---------------------------------------------------------------- data
def symbols(a, sym_col: str) -> np.ndarray:
    if sym_col in a.var.columns:
        return a.var[sym_col].astype(str).to_numpy()
    return np.asarray(a.var_names).astype(str)


def read_expression(a, rows: np.ndarray, cols: np.ndarray, chunk: int) -> np.ndarray:
    """log1p CP10k of the `cols` genes for `rows` (sorted), CP10k with the
    total counts over all genes; the file is read in blocks of rows."""
    out = np.zeros((len(rows), len(cols)), dtype=np.float32)
    t = time.time()
    for lo in range(0, a.n_obs, chunk):
        hi = min(lo + chunk, a.n_obs)
        k0, k1 = np.searchsorted(rows, [lo, hi])
        if k0 == k1:
            continue
        x = a.X[lo:hi]
        x = (sp.csr_matrix(x) if not sp.issparse(x) else x.tocsr())[rows[k0:k1] - lo]
        tot = np.asarray(x.sum(1)).ravel()
        m = x[:, cols].toarray().astype(np.float64)
        out[k0:k1] = np.log1p(m * (1e4 / np.maximum(tot, 1))[:, None])
        print(f"    {k1:,}/{len(rows):,} cells ({time.time() - t:.0f}s)", flush=True)
    return out


def umap_coords(a, names: pd.Index, umap_from: str | None):
    """(xy, found mask) for `names`: from `a` itself, else from `umap_from`."""
    import anndata as ad
    key = next((k for k in UMAP_KEYS if k in a.obsm), None)
    if key is not None:
        pos = a.obs_names.get_indexer(names)
        print(f"  UMAP: obsm['{key}'] of the input")
        return np.asarray(a.obsm[key])[:, :2][pos], np.ones(len(names), dtype=bool)
    if not umap_from:
        raise KeyError("No obsm['X_umap'] in --h5ad; pass --umap-from")
    src = ad.read_h5ad(umap_from, backed="r")
    try:
        key = next((k for k in UMAP_KEYS if k in src.obsm), None)
        if key is None:
            raise KeyError(f"No obsm['X_umap'] in {umap_from}")
        pos = src.obs_names.get_indexer(names)
        found = pos >= 0
        print(f"  UMAP from {umap_from}: {int(found.sum()):,}/{len(names):,} cells matched")
        if not found.any():
            raise ValueError("No cells matched by ID")
        xy = np.full((len(names), 2), np.nan)
        xy[found] = np.asarray(src.obsm[key])[:, :2][pos[found]]
        return xy, found
    finally:
        src.file.close()


# ---------------------------------------------------------------- plots
def region_colors(regions: list[str]) -> dict:
    pal = plt.get_cmap("tab20" if len(regions) > 10 else "tab10")
    return {r: pal(k % pal.N) for k, r in enumerate(regions)}


def scatter_regions(ax, xy, reg, colors, order, point_size):
    ax.scatter(xy[order, 0], xy[order, 1], s=point_size,
               c=[colors[r] for r in reg[order]], linewidths=0, rasterized=True)
    counts = pd.Series(reg).value_counts()
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=6, color=colors[r],
                              label=f"{r} ({counts.get(r, 0):,})") for r in colors],
              loc="upper left", bbox_to_anchor=(0.0, -0.02), ncol=2, frameon=False,
              fontsize=7.5)


def scatter_expression(fig, ax, xy, v, point_size):
    zero = v <= 0
    ax.scatter(xy[zero, 0], xy[zero, 1], s=point_size, c=ZERO_COLOR, linewidths=0,
               rasterized=True)
    pos = np.flatnonzero(~zero)
    pos = pos[np.argsort(v[pos])]                       # highest on top
    vmax = float(np.percentile(v[pos], 99)) if pos.size else 1.0
    sca = ax.scatter(xy[pos, 0], xy[pos, 1], s=point_size * 1.5, c=v[pos], cmap=EXPR_CMAP,
                     vmin=0, vmax=max(vmax, 1e-3), linewidths=0, rasterized=True)
    cb = fig.colorbar(sca, ax=ax, fraction=0.04, pad=0.01)
    cb.set_label("log1p CP10k (99th pct. cap)", fontsize=8)


def region_bars(ax, stats: pd.DataFrame, colors):
    s = stats.sort_values("pct_detected")
    y = np.arange(len(s))
    ax.barh(y, s["pct_detected"], color=[colors[r] for r in s["region"]])
    ax.set_yticks(y, s["region"], fontsize=8)
    ax.set_xlabel("% cells with ≥ 1 count")
    for k, (p, m) in enumerate(zip(s["pct_detected"], s["mean_log1p_cp10k"])):
        ax.text(p, k, f" {p:.0f}%  (mean {m:.2f})", va="center", fontsize=7)
    ax.set_xlim(0, 135)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Per region", fontsize=10, fontweight="bold", loc="left")


def gene_figure(gene, group, set_name, xy, reg, v, stats, colors, order, point_size,
                subset_txt, path):
    fig, axes = plt.subplots(1, 3, figsize=(21, 7.2),
                             gridspec_kw={"width_ratios": [1, 1.1, 0.75]})
    for ax in axes[:2]:
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
        ax.spines[["top", "right"]].set_visible(False)
    scatter_regions(axes[0], xy, reg, colors, order, point_size)
    axes[0].set_title(f"Region{subset_txt}", fontsize=11, fontweight="bold", loc="left")
    scatter_expression(fig, axes[1], xy, v, point_size)
    axes[1].set_title(f"{gene} expression ({int((v > 0).sum()):,} cells with counts)",
                      fontsize=11, fontweight="bold", loc="left")
    region_bars(axes[2], stats, colors)
    fig.suptitle(f"{gene}: {group} marker ({set_name} maturity marker set)",
                 fontsize=14, fontweight="bold", y=1.0)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# ---------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description="UMAPs of the neuron maturity markers.")
    p.add_argument("--h5ad", required=True, help="human_dev h5ad, raw counts in X.")
    p.add_argument("--umap-from", default=None,
                   help="h5ad with obsm['X_umap'] when --h5ad has none.")
    p.add_argument("--marker-sets", nargs="+", default=list(MATURITY_MARKER_SETS),
                   choices=list(MATURITY_MARKER_SETS))
    p.add_argument("--genes", nargs="+", default=None,
                   help="Extra genes to plot (group 'extra').")
    p.add_argument("--sym-col", default="Gene")
    p.add_argument("--region-col", default="Region")
    p.add_argument("--cell-class", default=None,
                   help="Only cells with --class-col == this (e.g. Neuron).")
    p.add_argument("--class-col", default="CellClass")
    p.add_argument("--chemistry", default=None, help="Only this --chemistry-col value.")
    p.add_argument("--chemistry-col", default="Chemistry")
    p.add_argument("--subsample", type=float, default=None,
                   help="Plot a random subset of cells (fraction or number); the "
                        "per-region numbers use all cells.")
    p.add_argument("--point-size", type=float, default=0.3)
    p.add_argument("--chunk", type=int, default=50000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--label", required=True)
    p.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    args = p.parse_args()

    import anndata as ad
    out_dir = Path(args.output_dir) / args.label
    out_dir.mkdir(parents=True, exist_ok=True)

    markers = []                                  # (set, group, gene)
    for s in args.marker_sets:
        imm, mat = MATURITY_MARKER_SETS[s]
        markers += [(s, "immature", g) for g in imm] + [(s, "mature", g) for g in mat]
    markers += [("extra", "extra", g) for g in args.genes or []]

    print(f"Loading (backed): {args.h5ad}")
    a = ad.read_h5ad(args.h5ad, backed="r")
    obs = a.obs
    keep = np.ones(a.n_obs, dtype=bool)
    subset_txt = ""
    if args.cell_class:
        keep &= (obs[args.class_col].astype(str) == args.cell_class).to_numpy()
        subset_txt += f", {args.cell_class} only"
    if args.chemistry:
        keep &= (obs[args.chemistry_col].astype(str) == args.chemistry).to_numpy()
        subset_txt += f", {args.chemistry}"
    rows = np.flatnonzero(keep)
    print(f"  {len(rows):,}/{a.n_obs:,} cells{subset_txt}")
    if not len(rows):
        raise ValueError("No cells left after the filters")

    sym = symbols(a, args.sym_col)
    first = {}
    for i, g in enumerate(sym):
        first.setdefault(g, i)
    missing = sorted({g for _, _, g in markers if g not in first})
    if missing:
        print(f"  [warning] not in the data, skipped: {missing}")
    markers = [m for m in markers if m[2] in first]
    genes = list(dict.fromkeys(g for _, _, g in markers))
    print(f"  {len(genes)} genes: {genes}")

    print("Expression (log1p CP10k)")
    expr = read_expression(a, rows, np.array([first[g] for g in genes]), args.chunk)
    reg_all = obs[args.region_col].astype(str).to_numpy()[rows]
    names = obs.index[rows]
    xy, found = umap_coords(a, names, args.umap_from)
    a.file.close()

    # per-region numbers, all selected cells
    recs = []
    for j, g in enumerate(genes):
        df = pd.DataFrame({"region": reg_all, "v": expr[:, j]})
        st = df.groupby("region")["v"].agg(cells="size",
                                           pct_detected=lambda x: 100 * (x > 0).mean(),
                                           mean_log1p_cp10k="mean").reset_index()
        st.insert(0, "gene", g)
        recs.append(st)
    stats = pd.concat(recs, ignore_index=True)
    tab = pd.DataFrame(markers, columns=["set", "group", "gene"]).merge(stats, on="gene")
    tab.to_csv(out_dir / "expression_by_region.csv", index=False)

    # cells plotted
    plot = np.flatnonzero(found)
    if args.subsample:
        k = (int(round(args.subsample * len(plot))) if args.subsample <= 1
             else int(args.subsample))
        if k < len(plot):
            plot = np.sort(np.random.default_rng(args.seed).choice(plot, k, replace=False))
            print(f"  plotting {k:,} random cells")
    regions = list(pd.Series(reg_all).value_counts().index)
    colors = region_colors(regions)
    order = np.random.default_rng(args.seed).permutation(len(plot))
    xy_p, reg_p = xy[plot], reg_all[plot]
    for s, grp, g in markers:
        j = genes.index(g)
        gene_figure(g, grp, s, xy_p, reg_p, expr[plot, j], stats[stats["gene"] == g],
                    colors, order, args.point_size, subset_txt,
                    out_dir / s / f"{grp}_{g}.png")
    print(f"Outputs: {out_dir}")


if __name__ == "__main__":
    main()
