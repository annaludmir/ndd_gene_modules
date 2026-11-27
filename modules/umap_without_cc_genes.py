"""
umap_without_cc_genes.py

A new UMAP of an h5ad (e.g. human_dev) whose highly variable genes do not
include cell cycle genes, so the embedding is not organised by cell cycle
phase.

--chemistry keeps one chemistry only (e.g. v3): v2 and v3 differ enough to
split a joint UMAP by chemistry, so each gets its own UMAP.

Steps (X = raw counts):
  1. Total counts per cell, over all genes (read in chunks, the file stays
     on disk).
  2. Highly variable genes on a random subset of cells (--hvg-cells),
     log1p CP10k, scanpy flavor "seurat" with --batch-key (Chemistry): the
     top --n-hvg genes among the genes that are NOT in --cc-genes (all gene
     sets of the JSON together). So there are still --n-hvg genes, none of
     them a cell cycle gene. For reference, the HVGs with the cell cycle
     genes included are computed too: hvg_genes.csv says which cell cycle
     genes would have been HVGs.
  3. All cells, HVG columns only: log1p CP10k (with the totals of step 1),
     scaled (max 10), PCA (--n-pcs), optional Harmony on --batch-key
     (--harmony, needs harmonypy), neighbours, UMAP.

Outputs (results/additional_analyses/umap_no_cc_genes/<label>/):
  umap_no_cc_genes.h5ad   obs + obsm X_umap / X_pca, no expression; use it
                          as --umap-from, e.g. for cell_cycle_umap.py
  hvg_genes.csv           every gene with its HVG rank without / with the
                          cell cycle genes, and whether it is a cell cycle gene
  umap_<column>.png       per --color column (CellCyclePhase: same colours as
                          cell_cycle_umap.py, Non-cycling grey underneath)

Usage:
  python modules/umap_without_cc_genes.py \
    --h5ad data/human_dev_without_week_5_cc_annotated_rescaled.h5ad \
    --cc-genes config_files/cell_cycle_genes_tirosh.json --chemistry v3 \
    --label human_dev_v3
"""

import argparse
import json
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
from cell_cycle_umap import NON_CYCLING, plot_all  # noqa: E402

DEFAULT_OUTPUT_DIR = Path("results/additional_analyses/umap_no_cc_genes")


# ---------------------------------------------------------------- inputs
def cc_gene_set(path: str, sets: list[str] | None) -> set[str]:
    d = json.load(open(path))
    keys = sets or list(d)
    missing = [k for k in keys if k not in d]
    if missing:
        raise KeyError(f"{missing} not in {path} (has {list(d)})")
    genes = set()
    for k in keys:
        v = d[k]
        genes |= set(v["geneSymbols"] if isinstance(v, dict) else v)
    print(f"  cell cycle genes: {len(genes)} from {keys}")
    return genes


def symbols(a, sym_col: str) -> np.ndarray:
    if sym_col in a.var.columns:
        return a.var[sym_col].astype(str).to_numpy()
    return np.asarray(a.var_names).astype(str)


def blocks(a, rows: np.ndarray, chunk: int):
    """Yield (positions in `rows`, CSR counts of those rows), reading the file
    in contiguous blocks of `chunk` rows (rows sorted)."""
    t = time.time()
    for lo in range(0, a.n_obs, chunk):
        hi = min(lo + chunk, a.n_obs)
        k0, k1 = np.searchsorted(rows, [lo, hi])
        if k0 == k1:
            continue
        x = a.X[lo:hi]
        x = (sp.csr_matrix(x) if not sp.issparse(x) else x.tocsr())[rows[k0:k1] - lo]
        yield slice(k0, k1), x
        print(f"    {k1:,}/{len(rows):,} cells ({time.time() - t:.0f}s)", flush=True)


# ---------------------------------------------------------------- steps
def totals_and_subset(a, rows: np.ndarray, in_sub: np.ndarray, chunk: int):
    """Total counts of `rows` (all genes), and the full count matrix of the
    rows flagged in `in_sub`, for the HVG selection."""
    tot = np.zeros(len(rows), dtype=np.float64)
    parts = []
    for sl, x in blocks(a, rows, chunk):
        tot[sl] = np.asarray(x.sum(1)).ravel()
        keep = in_sub[sl]
        if keep.any():
            parts.append(x[keep])
    return tot, sp.vstack(parts).tocsr()


def select_hvg(sub_x, tot_sub, obs_sub, genes, cc, n_hvg, batch_key):
    import anndata as ad
    import scanpy as sc
    x = sub_x.multiply(1e4 / np.maximum(tot_sub, 1)[:, None]).tocsr()
    x.data = np.log1p(x.data)
    h = ad.AnnData(X=x.astype(np.float32), obs=obs_sub.copy(),
                   var=pd.DataFrame(index=pd.Index(genes).astype(str)))
    h.var_names_make_unique()
    bk = batch_key if batch_key in h.obs and h.obs[batch_key].nunique() > 1 else None
    is_cc = np.isin(genes, list(cc))

    def ranks(adata):
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat", batch_key=bk)
        v = adata.var
        score = (v["highly_variable_nbatches"] if "highly_variable_nbatches" in v
                 else pd.Series(0, index=v.index))
        order = (pd.DataFrame({"hv": v["highly_variable"], "nb": score,
                               "disp": v["dispersions_norm"].fillna(-np.inf)})
                 .sort_values(["hv", "nb", "disp"], ascending=False))
        r = pd.Series(np.nan, index=v.index)
        hv = order.index[order["hv"]]
        r[hv] = np.arange(1, len(hv) + 1)
        return r

    with_cc = ranks(h.copy())
    no_cc = ranks(h[:, ~is_cc].copy())
    table = pd.DataFrame({"gene": h.var_names, "cell_cycle_gene": is_cc,
                          "hvg_rank_without_cc": no_cc.reindex(h.var_names).to_numpy(),
                          "hvg_rank_with_cc": with_cc.reindex(h.var_names).to_numpy()})
    n_cc_would = int((table["cell_cycle_gene"] & table["hvg_rank_with_cc"].notna()).sum())
    print(f"  HVGs: {int(table['hvg_rank_without_cc'].notna().sum()):,} without cell cycle "
          f"genes ({int(is_cc.sum())} cell cycle genes in the data left out; "
          f"{n_cc_would} of them would have been HVGs)")
    cols = np.flatnonzero(table["hvg_rank_without_cc"].notna().to_numpy())
    return cols, table


def hvg_matrix(a, rows, tot, cols, chunk):
    parts = []
    for sl, x in blocks(a, rows, chunk):
        x = x[:, cols].multiply(1e4 / np.maximum(tot[sl], 1)[:, None]).tocsr()
        x = x.astype(np.float32)
        x.data = np.log1p(x.data)
        parts.append(x)
    return sp.vstack(parts).tocsr()


# ---------------------------------------------------------------- plots
def plot_category(xy, values, title, path, point_size, seed):
    vals = pd.Series(values).astype(str).replace({"nan": "NA"})
    cats = vals.value_counts()
    pal = plt.get_cmap("tab20" if len(cats) > 10 else "tab10")
    colors = {c: pal(k % pal.N) for k, c in enumerate(cats.index)}
    order = np.random.default_rng(seed).permutation(len(vals))
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(xy[order, 0], xy[order, 1], s=point_size,
               c=[colors[v] for v in vals.to_numpy()[order]], linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(title, fontsize=13, fontweight="bold", loc="left")
    if len(cats) <= 40:
        ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=7, color=colors[c],
                                  label=f"{c} ({n:,})") for c, n in cats.items()],
                  loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, fontsize=8)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# ---------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description="UMAP without cell cycle genes in the HVGs.")
    p.add_argument("--h5ad", required=True, help="h5ad with raw counts in X.")
    p.add_argument("--cc-genes", default="config_files/cell_cycle_genes_tirosh.json",
                   help="JSON {set: {geneSymbols: [...]}} or {set: [...]}.")
    p.add_argument("--cc-sets", nargs="+", default=None,
                   help="Only these gene sets of the JSON. Default: all.")
    p.add_argument("--sym-col", default="Gene",
                   help="var column with gene symbols (else var_names).")
    p.add_argument("--n-hvg", type=int, default=3000)
    p.add_argument("--hvg-cells", type=int, default=200000,
                   help="Random cells used to pick the HVGs (0 = all).")
    p.add_argument("--batch-key", default="Chemistry",
                   help="obs column for HVG selection per batch (and Harmony).")
    p.add_argument("--n-pcs", type=int, default=50)
    p.add_argument("--harmony", action="store_true",
                   help="Harmony on --batch-key before the neighbours (needs harmonypy).")
    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--min-dist", type=float, default=0.5)
    p.add_argument("--chemistry", default=None,
                   help="Keep only cells with this --chemistry-col value (e.g. v3).")
    p.add_argument("--chemistry-col", default="Chemistry")
    p.add_argument("--max-cells", type=int, default=None,
                   help="Random subset of cells (for a test run).")
    p.add_argument("--color", nargs="+",
                   default=["CellCyclePhase", "CellClass", "Region", "RegionsGeneral",
                            "Chemistry"],
                   help="obs columns to colour the UMAP by (missing ones are skipped).")
    p.add_argument("--point-size", type=float, default=0.3)
    p.add_argument("--chunk", type=int, default=50000, help="Cells read at a time.")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--label", required=True)
    p.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    args = p.parse_args()

    import anndata as ad
    import scanpy as sc

    out_dir = Path(args.output_dir) / args.label
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    cc = cc_gene_set(args.cc_genes, args.cc_sets)

    print(f"Loading (backed): {args.h5ad}")
    a = ad.read_h5ad(args.h5ad, backed="r")
    genes = symbols(a, args.sym_col)
    print(f"  {a.n_obs:,} cells x {a.n_vars:,} genes; "
          f"{int(np.isin(genes, list(cc)).sum())} cell cycle genes found "
          f"(by {'var[' + repr(args.sym_col) + ']' if args.sym_col in a.var else 'var_names'})")
    rows = np.arange(a.n_obs)
    if args.chemistry:
        chem = a.obs[args.chemistry_col].astype(str)
        rows = np.flatnonzero((chem == args.chemistry).to_numpy())
        print(f"  {args.chemistry_col} == {args.chemistry}: {len(rows):,} cells")
        if not len(rows):
            raise ValueError(f"No cells; {args.chemistry_col} values: {sorted(chem.unique())}")
    if args.max_cells and args.max_cells < a.n_obs:
        rows = np.sort(rng.choice(rows, args.max_cells, replace=False))
        print(f"  test subset: {len(rows):,} cells")
    obs = a.obs.iloc[rows].copy()
    n_sub = len(rows) if not args.hvg_cells else min(args.hvg_cells, len(rows))
    in_sub = np.zeros(len(rows), dtype=bool)
    in_sub[rng.choice(len(rows), n_sub, replace=False)] = True

    print(f"Step 1: total counts per cell + {n_sub:,} random cells for the HVGs")
    tot, sub_x = totals_and_subset(a, rows, in_sub, args.chunk)
    sample = sub_x.data[:100000]
    if sample.size and not np.allclose(sample, np.round(sample)):
        print("  [warning] X does not look like raw counts")

    print("Step 2: highly variable genes (cell cycle genes left out)")
    cols, table = select_hvg(sub_x, tot[in_sub], obs.iloc[np.flatnonzero(in_sub)],
                             genes, cc, args.n_hvg, args.batch_key)
    del sub_x
    table.to_csv(out_dir / "hvg_genes.csv", index=False)

    print(f"Step 3: all cells, {len(cols):,} HVGs: log1p CP10k, scale, PCA, UMAP")
    x = hvg_matrix(a, rows, tot, cols, args.chunk)
    a.file.close()
    adata = ad.AnnData(X=x, obs=obs, var=pd.DataFrame(index=pd.Index(genes[cols])))
    adata.var_names_make_unique()
    t = time.time()
    adata.X = adata.X.toarray()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=args.n_pcs, random_state=args.seed)
    print(f"  PCA done ({time.time() - t:.0f}s)")
    rep = "X_pca"
    if args.harmony:
        import harmonypy
        ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, args.batch_key,
                                   random_state=args.seed)
        z = np.asarray(ho.Z_corr)
        adata.obsm["X_pca_harmony"] = z.T if z.shape[0] != adata.n_obs else z
        rep = "X_pca_harmony"
        print(f"  Harmony on {args.batch_key} ({time.time() - t:.0f}s)")
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep=rep,
                    random_state=args.seed)
    print(f"  neighbours done ({time.time() - t:.0f}s)")
    sc.tl.umap(adata, min_dist=args.min_dist, random_state=args.seed)
    print(f"  UMAP done ({time.time() - t:.0f}s)")

    res = ad.AnnData(obs=adata.obs.copy(),
                     obsm={k: adata.obsm[k] for k in adata.obsm})
    res.uns["umap_without_cc_genes"] = {
        "source": str(args.h5ad), "cc_genes": str(args.cc_genes),
        "n_hvg": int(len(cols)), "batch_key": str(args.batch_key),
        "chemistry": str(args.chemistry),
        "harmony": bool(args.harmony), "n_pcs": int(args.n_pcs),
        "n_neighbors": int(args.n_neighbors), "min_dist": float(args.min_dist)}
    res.write_h5ad(out_dir / "umap_no_cc_genes.h5ad")
    print(f"  Saved: {out_dir / 'umap_no_cc_genes.h5ad'}")

    xy = res.obsm["X_umap"]
    sfx = " (HVGs without cell cycle genes)"
    for col in args.color:
        if col not in res.obs:
            print(f"  [skip plot] obs['{col}'] not found")
            continue
        path = out_dir / f"umap_{col}.png"
        if col == "CellCyclePhase":
            ph = (res.obs[col].astype(str).str.replace("PostM", "Post-M")
                  .replace({"nan": NON_CYCLING}).to_numpy())
            plot_all(xy, ph, f"Cell cycle{sfx}", path, args.point_size, args.seed)
        else:
            plot_category(xy, res.obs[col].to_numpy(), f"{col}{sfx}", path,
                          args.point_size, args.seed)
    print(f"Outputs: {out_dir}")


if __name__ == "__main__":
    main()
