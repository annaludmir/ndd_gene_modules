"""
region_expression_fraction.py

Do the cells of one class (e.g. Neuronal IPC) express more of a gene set in
one region than another because of the genes, or because the cells of that
region detect more genes overall? Compares two (meta-)regions, Forebrain vs
Midbrain by default, in three steps:

  1. depth: genes detected and UMIs per cell, per region. If region A's cells
     detect many more genes, every gene looks more often "expressed" there.
  2. every gene: % of cells expressing it (UMI > 0) in each region, raw and
     after downsampling every cell to the same UMI depth (binomial thinning to
     --downsample-to UMIs; cells with fewer UMIs are left out of this part).
     The depth-matched difference is the one to read: what remains is not a
     sequencing-depth effect. Two-proportion z-test, BH FDR.
  3. the gene set (e.g. the GSEA leading genes): is its mean A - B difference
     larger than that of random gene sets with the same detection levels?
     Genes are binned by overall % expressing; each random set draws, per set
     gene, a gene from the same bin (--n-random sets). Each set gene also gets
     its percentile among the genes of its bin.

Cells are not independent (donors / samples differ between regions), so the
z-test p-values are optimistic; rank genes by the difference and use them as a
guide, not as proof.

Outputs (results/additional_analyses/region_expression_fraction/<subfolder>/):
  cell_depth_by_region.csv / .png      genes and UMIs per cell, per region
  gene_fraction_by_region.csv          every gene: pct_<A>, pct_<B>, diff, depth-
                                       matched pct_*_ds, diff_ds, fdr_ds, in_gene_set
  top_<A>_over_<B>.csv / top_<B>_over_<A>.csv
                                       --n-top genes each way (depth-matched diff)
  gene_set_vs_matched_random.csv / .png
  scatter_pct_<A>_vs_<B>.png           every gene (grey), set genes highlighted
  top_genes_<A>_vs_<B>.png             the top genes each way, as dot pairs

Usage:
  python modules/region_expression_fraction.py \
    --h5ad-path data/human_dev.h5ad --cell-class "Neuronal IPC" \
    --gsea-summary-file <GSEA_final_summary.csv> --condition Forebrain_General \
    --subfolder-name neuronal_ipc_microcephaly
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gene_expression_summary import (
    collapse_region, load_summary_row, parse_gene_list, read_table, sanitize_name,
)

DEFAULT_OUTPUT_DIR = Path("results/additional_analyses/region_expression_fraction")
DEFAULT_N_RANDOM = 10000
DEFAULT_N_BINS = 20
DEFAULT_N_TOP = 50
DEFAULT_MIN_PCT = 1.0
DEFAULT_DOWNSAMPLE_QUANTILE = 0.1
COLORS = {"A": "#c0392b", "B": "#2c7fb8", "set": "#e67e22"}


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values (NaN kept)."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(pv)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def two_proportion_p(k1, n1, k2, n2) -> np.ndarray:
    """Two-sided two-proportion z-test p-values, vectorized over genes."""
    from scipy.stats import norm
    p1, p2 = k1 / n1, k2 / n2
    pool = (k1 + k2) / (n1 + n2)
    se = np.sqrt(pool * (1 - pool) * (1 / n1 + 1 / n2))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(se > 0, (p1 - p2) / se, 0.0)
    return 2 * norm.sf(np.abs(z))


def detected_counts(X, mask: np.ndarray) -> np.ndarray:
    """Per gene: number of cells in `mask` with UMI > 0 (X csr)."""
    sub = X[mask]
    return np.asarray((sub > 0).sum(axis=0)).ravel()


def downsample(X, target: int, rng) -> tuple:
    """Binomial thinning of every cell with >= target UMIs to ~target UMIs.
    Returns (thinned csr of those cells, boolean mask of kept cells)."""
    import scipy.sparse as sp
    total = np.asarray(X.sum(axis=1)).ravel()
    keep = total >= target
    sub = X[keep].tocsr(copy=True).astype(np.int64)
    p_row = target / total[keep]
    row_of = np.repeat(np.arange(sub.shape[0]), np.diff(sub.indptr))
    sub.data = rng.binomial(sub.data, p_row[row_of])
    sub.eliminate_zeros()
    return sp.csr_matrix(sub), keep


def load_gene_set(args) -> tuple[list[str], str]:
    """The gene set and a label for it, from --gene-list or a GSEA summary row."""
    if args.gene_list:
        df = read_table(args.gene_list)
        col = "gene" if "gene" in df.columns else df.columns[0]
        genes = parse_gene_list(" ".join(df[col].dropna().astype(str)))
        return genes, Path(args.gene_list).stem
    if args.gsea_summary_file and args.condition:
        row = load_summary_row(args.gsea_summary_file, args.condition)
        return row["leading_genes"], f"leading genes ({row['condition']})"
    return [], ""


def load_cells(args):
    """Cells of --cell-class in the two regions (and chemistry), counts in csr."""
    import anndata as ad
    import scipy.sparse as sp
    print(f"Loading (backed): {args.h5ad_path}")
    a = ad.read_h5ad(args.h5ad_path, backed="r")
    obs = a.obs
    for col in (args.region_col, args.cell_class_col):
        if col not in obs.columns:
            raise KeyError(f"obs column '{col}' not found")
    region = (obs[args.region_col].map(collapse_region) if not args.no_collapse
              else obs[args.region_col].astype(str))
    mask = (obs[args.cell_class_col].astype(str) == args.cell_class).to_numpy()
    mask &= region.isin([args.region_a, args.region_b]).to_numpy()
    if args.chemistry:
        mask &= (obs[args.chemistry_col].astype(str) == args.chemistry).to_numpy()
    rows = np.flatnonzero(mask)
    print(f"  {len(rows):,} cells: {args.cell_class_col} == '{args.cell_class}', "
          f"{args.region_a} / {args.region_b}"
          + (f", {args.chemistry_col} == '{args.chemistry}'" if args.chemistry else ""))
    if len(rows) == 0:
        raise ValueError("No cells selected.")
    X = a[rows].to_memory().X
    X = sp.csr_matrix(X)
    labels = region.to_numpy()[rows].astype(str)
    sym = (a.var[args.sym_col].astype(str).to_numpy() if args.sym_col in a.var.columns
           else a.var_names.astype(str).to_numpy())
    return X, labels, sym


def cell_depth(X, labels, A, B, out_dir) -> pd.DataFrame:
    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    umis = np.asarray(X.sum(axis=1)).ravel()
    df = pd.DataFrame({"region": labels, "n_genes_detected": n_genes, "total_umis": umis})
    summ = df.groupby("region").agg(
        n_cells=("total_umis", "size"),
        median_genes=("n_genes_detected", "median"), mean_genes=("n_genes_detected", "mean"),
        median_umis=("total_umis", "median"), mean_umis=("total_umis", "mean")).loc[[A, B]]
    summ.to_csv(out_dir / "cell_depth_by_region.csv")
    print("\n  Depth per cell:\n" + summ.round(1).to_string())

    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    for ax, col, title in zip(axes, ["n_genes_detected", "total_umis"],
                              ["genes detected per cell", "UMIs per cell"]):
        vals = [df.loc[df.region == r, col] for r in (A, B)]
        bp = ax.boxplot(vals, tick_labels=[f"{A}\n(n={len(vals[0]):,})", f"{B}\n(n={len(vals[1]):,})"],
                        showfliers=False, patch_artist=True, widths=0.6)
        for patch, c in zip(bp["boxes"], (COLORS["A"], COLORS["B"])):
            patch.set_facecolor(c); patch.set_alpha(0.6)
        ax.set_title(title, fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_dir / "cell_depth_by_region.png", dpi=200)
    plt.close(fig)
    return summ


def gene_fractions(X, labels, sym, A, B, target, rng) -> pd.DataFrame:
    a, b = labels == A, labels == B
    nA, nB = int(a.sum()), int(b.sum())
    kA, kB = detected_counts(X, a), detected_counts(X, b)
    Xd, keep = downsample(X, target, rng)
    la = labels[keep]
    nAd, nBd = int((la == A).sum()), int((la == B).sum())
    print(f"\n  Depth-matched: every cell thinned to {target:,} UMIs; kept "
          f"{nAd:,}/{nA:,} {A} and {nBd:,}/{nB:,} {B} cells")
    kAd, kBd = detected_counts(Xd, la == A), detected_counts(Xd, la == B)
    df = pd.DataFrame({
        "gene": sym,
        f"pct_{A}": 100 * kA / nA, f"pct_{B}": 100 * kB / nB,
        f"pct_{A}_ds": 100 * kAd / max(nAd, 1), f"pct_{B}_ds": 100 * kBd / max(nBd, 1),
    })
    df["diff"] = df[f"pct_{A}"] - df[f"pct_{B}"]
    df["diff_ds"] = df[f"pct_{A}_ds"] - df[f"pct_{B}_ds"]
    df["pct_all_ds"] = 100 * (kAd + kBd) / max(nAd + nBd, 1)
    df["fdr_ds"] = bh(two_proportion_p(kAd, max(nAd, 1), kBd, max(nBd, 1)))
    df.attrs.update(nA=nA, nB=nB, nAd=nAd, nBd=nBd, target=target)
    return df


def gene_set_test(df, set_genes, n_bins, n_random, rng, A, B, out_dir, set_label):
    """Set genes' mean diff_ds vs random detection-matched sets."""
    df = df.copy()
    df["bin"] = pd.qcut(df["pct_all_ds"].rank(method="first"), n_bins, labels=False)
    in_set = df["gene"].isin(set_genes).to_numpy()
    df_set = df[in_set]
    bins = df.groupby("bin").indices
    diff = df["diff_ds"].to_numpy()
    obs_mean = float(df_set["diff_ds"].mean())
    draws = np.zeros(n_random)
    for b_ix, n in df_set["bin"].value_counts().items():
        pool = bins[b_ix]
        draws += diff[rng.choice(pool, size=(n_random, n), replace=True)].sum(axis=1)
    draws /= len(df_set)
    p_hi = (1 + (draws >= obs_mean).sum()) / (n_random + 1)
    p_lo = (1 + (draws <= obs_mean).sum()) / (n_random + 1)
    pct_rank = df.groupby("bin")["diff_ds"].rank(pct=True)
    per_gene = df_set[["gene", f"pct_{A}_ds", f"pct_{B}_ds", "diff_ds", "fdr_ds"]].assign(
        percentile_in_detection_bin=100 * pct_rank[in_set].to_numpy())
    summary = pd.DataFrame([{
        "gene_set": set_label, "n_genes_found": int(in_set.sum()),
        "mean_diff_ds_gene_set": obs_mean,
        "mean_diff_ds_all_genes": float(np.mean(diff)),
        "mean_diff_ds_random_matched": float(draws.mean()),
        "random_matched_p2.5": float(np.percentile(draws, 2.5)),
        "random_matched_p97.5": float(np.percentile(draws, 97.5)),
        f"p_set_more_{A}_than_random": p_hi, f"p_set_more_{B}_than_random": p_lo,
        "n_random": n_random,
    }])
    summary.to_csv(out_dir / "gene_set_vs_matched_random.csv", index=False)
    per_gene.sort_values("diff_ds", ascending=False).to_csv(
        out_dir / "gene_set_genes_percentiles.csv", index=False)
    print("\n  Gene set vs random detection-matched sets:\n"
          + summary.T.to_string(header=False))

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(draws, bins=60, color="#bbbbbb", edgecolor="white")
    ax.axvline(obs_mean, color=COLORS["set"], lw=2.5, label=f"{set_label} ({obs_mean:+.2f})")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel(f"mean % expressing, {A} − {B} (depth-matched)")
    ax.set_ylabel("random matched gene sets")
    ax.set_title(f"Gene set vs {n_random:,} random sets with the same detection levels",
                 fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_dir / "gene_set_vs_matched_random.png", dpi=200)
    plt.close(fig)


def scatter(df, set_genes, A, B, out_dir, set_label, suffix, label_max=30):
    x, y = df[f"pct_{B}{suffix}"], df[f"pct_{A}{suffix}"]
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, s=3, c="#cccccc", lw=0, rasterized=True, label="all genes")
    s = df["gene"].isin(set_genes)
    if s.any():
        ax.scatter(x[s], y[s], s=18, c=COLORS["set"], edgecolor="black", lw=0.4,
                   label=set_label, zorder=3)
        top = df[s].reindex(df[s][f"diff{suffix}"].abs().sort_values(ascending=False).index)
        for _, r in top.head(label_max).iterrows():
            ax.annotate(r["gene"], (r[f"pct_{B}{suffix}"], r[f"pct_{A}{suffix}"]),
                        fontsize=6, xytext=(2, 2), textcoords="offset points")
    lim = max(x.max(), y.max()) * 1.03
    ax.plot([0, lim], [0, lim], color="black", lw=0.8, ls="--")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel(f"% {B} cells expressing"); ax.set_ylabel(f"% {A} cells expressing")
    ax.set_title(f"{A} vs {B}" + (" (depth-matched)" if suffix else " (raw)"), fontsize=11)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_dir / f"scatter_pct_{A}_vs_{B}{suffix}.png", dpi=200)
    plt.close(fig)


def top_genes(df, A, B, n_top, min_pct, out_dir, set_genes):
    ok = df[[f"pct_{A}_ds", f"pct_{B}_ds"]].max(axis=1) >= min_pct
    d = df[ok].assign(in_gene_set=df.loc[ok, "gene"].isin(set_genes))
    up = d.sort_values("diff_ds", ascending=False).head(n_top)
    down = d.sort_values("diff_ds").head(n_top)
    up.to_csv(out_dir / f"top_{A}_over_{B}.csv", index=False)
    down.to_csv(out_dir / f"top_{B}_over_{A}.csv", index=False)

    show = min(n_top, 30)
    fig, axes = plt.subplots(1, 2, figsize=(9, 0.22 * show + 1.5), sharex=True)
    for ax, part, title in zip(axes, (up.head(show), down.head(show)),
                               (f"higher in {A}", f"higher in {B}")):
        yy = np.arange(len(part))[::-1]
        ax.hlines(yy, part[f"pct_{B}_ds"], part[f"pct_{A}_ds"], color="#999999", lw=1)
        ax.scatter(part[f"pct_{A}_ds"], yy, color=COLORS["A"], s=18, label=A, zorder=3)
        ax.scatter(part[f"pct_{B}_ds"], yy, color=COLORS["B"], s=18, label=B, zorder=3)
        ax.set_yticks(yy)
        ax.set_yticklabels(part["gene"], fontsize=7)
        for t, flag in zip(ax.get_yticklabels(), part["in_gene_set"]):
            if flag:
                t.set_color(COLORS["set"]); t.set_fontweight("bold")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("% cells expressing (depth-matched)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(out_dir / f"top_genes_{A}_vs_{B}.png", dpi=200)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.strip().splitlines()[2])
    p.add_argument("--h5ad-path", required=True)
    p.add_argument("--cell-class", default="Neuronal IPC")
    p.add_argument("--cell-class-col", default="CellClass")
    p.add_argument("--region-col", default="Region")
    p.add_argument("--region-a", default="Forebrain")
    p.add_argument("--region-b", default="Midbrain")
    p.add_argument("--no-collapse", action="store_true",
                   help="Compare raw --region-col values instead of the Forebrain / Midbrain / "
                        "Hindbrain meta-regions (Forebrain = Forebrain + Telencephalon + "
                        "Diencephalon).")
    p.add_argument("--chemistry", default="v3", help="Keep one chemistry; '' for all.")
    p.add_argument("--chemistry-col", default="Chemistry")
    p.add_argument("--sym-col", default="Gene")
    p.add_argument("--gene-list", default=None,
                   help="CSV/TSV with a 'gene' column (or genes in the first column).")
    p.add_argument("--gsea-summary-file", default=None,
                   help="GSEA_final_summary.csv; with --condition, its Lead_genes are the set.")
    p.add_argument("--condition", default=None)
    p.add_argument("--downsample-to", type=int, default=None,
                   help=f"UMIs per cell for the depth-matched comparison. Default: the "
                        f"{DEFAULT_DOWNSAMPLE_QUANTILE:.0%} quantile of the selected cells.")
    p.add_argument("--n-random", type=int, default=DEFAULT_N_RANDOM)
    p.add_argument("--n-bins", type=int, default=DEFAULT_N_BINS)
    p.add_argument("--n-top", type=int, default=DEFAULT_N_TOP)
    p.add_argument("--min-pct", type=float, default=DEFAULT_MIN_PCT,
                   help="Top-gene lists: gene expressed in >= this %% of cells in one region.")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--subfolder-name", default=None)
    p.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    args = p.parse_args()
    A, B = args.region_a, args.region_b
    rng = np.random.default_rng(args.seed)

    set_genes, set_label = load_gene_set(args)
    sub = args.subfolder_name or sanitize_name(f"{args.cell_class}_{A}_vs_{B}")
    out_dir = Path(args.output_dir) / sub
    out_dir.mkdir(parents=True, exist_ok=True)

    X, labels, sym = load_cells(args)
    cell_depth(X, labels, A, B, out_dir)
    total = np.asarray(X.sum(axis=1)).ravel()
    target = args.downsample_to or int(np.quantile(total, DEFAULT_DOWNSAMPLE_QUANTILE))
    df = gene_fractions(X, labels, sym, A, B, target, rng)
    df = df.groupby("gene", as_index=False, sort=False).first()   # duplicate symbols
    df["in_gene_set"] = df["gene"].isin(set_genes)
    df.sort_values("diff_ds", ascending=False).to_csv(
        out_dir / "gene_fraction_by_region.csv", index=False)

    expressed = df[df[["pct_" + A + "_ds", "pct_" + B + "_ds"]].max(axis=1) > 0]
    print(f"\n  Genes expressed in either region (depth-matched): {len(expressed):,}; "
          f"median diff_ds {expressed['diff_ds'].median():+.3f} pct points, "
          f"{(expressed['diff_ds'] > 0).mean():.1%} higher in {A}")

    if set_genes:
        missing = sorted(set(set_genes) - set(expressed["gene"]))
        print(f"  Gene set '{set_label}': {len(set_genes) - len(missing)}/{len(set_genes)} "
              f"genes expressed here" + (f"; not: {missing[:15]}" if missing else ""))
        gene_set_test(expressed, set_genes, args.n_bins, args.n_random, rng, A, B,
                      out_dir, set_label)
    for suffix in ("", "_ds"):
        scatter(expressed, set_genes, A, B, out_dir, set_label or "gene set", suffix)
    top_genes(expressed, A, B, args.n_top, args.min_pct, out_dir, set_genes)
    print(f"\nOutputs: {out_dir}")


if __name__ == "__main__":
    main()
