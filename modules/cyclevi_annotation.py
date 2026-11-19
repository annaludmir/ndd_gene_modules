"""
cyclevi_annotation.py

Cell cycle position with CycleVI (Mozdzanowski, Tarbier & Jeuken, Bioinformatics
2026; github.com/jeuken/CycleVI), as a second method next to the marker-fraction
annotation of cell_cycle_annotation.py.

CycleVI is an scVI-type model whose latent space has a 2-D circular part: the
angle of a cell in it is its position in the cycle; the other 8 dimensions hold
the rest of the variation, with the cell cycle removed. It gives an angle, not
phase labels, and puts every cell on the circle (non-cycling cells gather in
the G1/G0 region), so:
  - it is trained once, on all progenitor cell classes of all regions together
    (--cell-class), with the chemistry as batch covariate (--batch-key):
    angles from separately trained models are not comparable (each model
    picks its own start and direction), and the paper's benchmarks show it
    absorbs batch / heterogeneity through the covariate. Neurons and other
    post-mitotic classes are left out: their angle carries no information.
  - phases are assigned afterwards from the angle (see "Phases" below).

Steps
  1. Cells of --cell-class (raw counts, all genes), optional random subsample.
  2. Initial guess, as in cyclevi.prepare.compute_phase_initialization:
     Seurat S / G2M scores (Tirosh 2016 lists; scanpy score_genes_cell_cycle),
     angle = atan2(G2M, S) quantile-spread over the circle, phase G1/S/G2M.
     Difference: scored on log-normalised counts (as in Seurat / scanpy), not
     on raw counts as the package's prepare step does.
  3. Genes: --n-hvg highly variable genes (per batch) plus all cell cycle
     marker genes; raw counts of those go to the model.
  4. Train. The paper trains 400 epochs of batch 128 on <= 31,634 cells
     (~100k gradient steps). With hundreds of thousands of cells that is far
     more than needed, so by default the number of epochs is chosen to give
     ~--target-steps steps (between 20 and 400 epochs). The loss schedules
     (angle supervision, adversary, radius penalty) follow the KL warm-up,
     which is set to the number of epochs so they run their full course.
  5. Outputs (results/additional_analyses/cyclevi/<label>/):
       model/                     the trained model (scvi-tools format)
       cells.csv                  per cell: obs columns, initial scores,
                                  z_cycle_x/y, cyclevi_angle (0..2pi),
                                  cyclevi_radius, cyclevi_phase
       latent_other.npy           the 8-D cell-cycle-free latent (row order
                                  of cells.csv)
       phase_bins.csv             angle bins and the reference phase mix
       crosstab_cyclevi_vs_<ref>.csv / .png
       angle_hist_by_<ref>.png    angle histogram stacked by reference phase
       marker_curves.png          smoothed marker expression along the angle
       angle_by_cellclass.png     angle density per cell class
       z_cycle_by_<ref>.png       the circular latent coloured by phase
       umap_cyclevi_phase.png     the dataset UMAP (--umap-from): CycleVI phase
                                  and reference phase side by side

Phases
  The angle has no fixed origin or direction. cyclevi_phase labels each of
  --n-bins angle bins with the reference phase most enriched in it (fraction
  in the bin / fraction overall). The reference is --reference-col
  (CellCyclePhase from cell_cycle_annotation.py when the input has it, which
  includes Non-cycling and Post-M), else the initial Seurat phase. Check it
  against marker_curves.png: CCNE1 / PCNA peak at G1/S, histones in S,
  CCNB1 / TOP2A at G2/M, CDKN1A in G1/G0.

--analyze-only redoes step 5's plots / phases from an existing cells.csv.

Usage:
  python modules/cyclevi_annotation.py \
    --h5ad data/human_dev_without_week_5_cc_annotated_rescaled.h5ad \
    --label human_dev_progenitors \
    --cell-class "Radial glia" "Neuronal IPC" Glioblast Oligo
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_OUTPUT_DIR = Path("results/additional_analyses/cyclevi")
DEFAULT_CLASSES = ["Radial glia", "Neuronal IPC", "Glioblast", "Oligo"]

# Seurat cc.genes = Tirosh et al. 2016 (s.genes 43, g2m.genes 54); current
# HGNC symbols added next to the old ones so either naming matches.
SEURAT_S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "CENPU", "HELLS", "RFC2", "RPA2",
    "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2",
    "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]
SEURAT_G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
    "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "PIMREG", "SMC4",
    "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B",
    "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "JPT1", "CDC20", "TTK", "CDC25C",
    "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23",
    "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2",
    "G2E3", "GAS2L3", "CBX5", "CENPA",
]
# Plotted along the angle (as in the paper's Fig. 2); kept in the model genes
MARKERS = ["CCNE1", "CCNE2", "PCNA", "MCM2", "H4C3", "HIST1H4C", "H1-2", "HIST1H1C",
           "CCNA2", "CCNB1", "TOP2A", "MKI67", "CDC20", "CDKN1A"]
PHASE_COLORS = {"G1": "#1f1fe0", "S": "#f5a623", "G2M": "#1e9e3c", "Post-M": "#e8231f",
                "Non-cycling": "#bdbdbd"}
TARGET_STEPS = 100_000


# ---------------------------------------------------------------- data
def symbols(var: pd.DataFrame, var_names, sym_col: str) -> np.ndarray:
    return (var[sym_col].astype(str).to_numpy() if sym_col in var.columns
            else np.asarray(var_names).astype(str))


def load_cells(args):
    import anndata as ad
    import scipy.sparse as sp
    print(f"Loading (backed): {args.h5ad}")
    a = ad.read_h5ad(args.h5ad, backed="r")
    obs = a.obs
    counts = obs[args.class_col].astype(str).value_counts()
    print(f"  {args.class_col} values: " + ", ".join(f"{k} {v:,}" for k, v in counts.items()))
    missing = [c for c in args.cell_class if c not in counts.index]
    if missing:
        raise ValueError(f"Not in obs['{args.class_col}']: {missing}")
    mask = obs[args.class_col].astype(str).isin(args.cell_class).to_numpy()
    if args.chemistry:
        mask &= (obs[args.batch_key].astype(str) == args.chemistry).to_numpy()
    rows = np.flatnonzero(mask)
    print(f"  {len(rows):,} cells in {args.cell_class}")
    if args.max_cells and len(rows) > args.max_cells:
        rows = np.sort(np.random.default_rng(args.seed).choice(rows, args.max_cells,
                                                                replace=False))
        print(f"  random subsample: {len(rows):,} cells")
    t = time.time()
    sub = a[rows].to_memory()
    a.file.close()
    sub.X = sp.csr_matrix(sub.X)
    print(f"  loaded {sub.n_obs:,} x {sub.n_vars:,} in {time.time() - t:.0f}s")
    sub.var_names = pd.Index(symbols(sub.var, sub.var_names, args.sym_col))
    sub.var_names_make_unique()
    vals = sub.X.data[:100000]
    if vals.size and not np.allclose(vals, np.round(vals)):
        print("  [warning] X is not integer counts; CycleVI expects raw counts")
    return sub


def initial_phase(adata, seed: int):
    """Seurat S / G2M scores on log-normalised data, angle atan2(G2M, S)
    quantile-spread to [0, 2pi), phase G1/S/G2M (cyclevi.prepare logic)."""
    import scanpy as sc
    from sklearn.preprocessing import QuantileTransformer
    lognorm = adata.copy()
    sc.pp.normalize_total(lognorm, target_sum=1e4)
    sc.pp.log1p(lognorm)
    genes = set(lognorm.var_names)
    s = [g for g in SEURAT_S_GENES if g in genes]
    g2m = [g for g in SEURAT_G2M_GENES if g in genes]
    print(f"  Seurat genes found: S {len(s)}, G2M {len(g2m)}")
    sc.tl.score_genes_cell_cycle(lognorm, s_genes=s, g2m_genes=g2m, random_state=seed)
    for c in ("S_score", "G2M_score", "phase"):
        adata.obs[c] = lognorm.obs[c].values
    adata.obs["phase"] = pd.Categorical(adata.obs["phase"].astype(str),
                                        categories=["G1", "S", "G2M"])
    raw = np.mod(np.arctan2(adata.obs["G2M_score"], adata.obs["S_score"]), 2 * np.pi)
    adata.obs["cycle_angle"] = raw
    qt = QuantileTransformer(output_distribution="uniform", random_state=0,
                             n_quantiles=min(len(raw), 1000), subsample=len(raw))
    adata.obs["cycle_angle_uniform"] = qt.fit_transform(
        np.asarray(raw).reshape(-1, 1)).ravel() * 2 * np.pi
    print("  initial phase: " + ", ".join(
        f"{k} {v:,}" for k, v in adata.obs["phase"].value_counts().items()))
    return lognorm


def select_genes(adata, lognorm, n_hvg: int, batch_key: str | None) -> list[str]:
    import scanpy as sc
    bk = batch_key if batch_key and lognorm.obs[batch_key].nunique() > 1 else None
    sc.pp.highly_variable_genes(lognorm, n_top_genes=n_hvg, flavor="seurat", batch_key=bk)
    hvg = set(lognorm.var_names[lognorm.var["highly_variable"].to_numpy()])
    keep = hvg | {g for g in SEURAT_S_GENES + SEURAT_G2M_GENES + MARKERS
                  if g in adata.var_names}
    genes = [g for g in adata.var_names if g in keep]
    print(f"  genes: {len(hvg):,} HVGs + cell cycle markers = {len(genes):,}")
    return genes


# ---------------------------------------------------------------- train
def train(adata, args, out_dir: Path):
    import torch
    from cyclevi import CycleVI
    bk = args.batch_key if args.batch_key in adata.obs.columns else None
    if args.batch_key and bk is None:
        print(f"  [note] no obs '{args.batch_key}': no batch covariate")
    CycleVI.setup_anndata(adata, batch_key=bk,
                          cycle_initiation_label_key="phase",
                          cycle_initiation_angle_key="cycle_angle_uniform")
    model = CycleVI(adata, n_latent=args.n_latent, n_hidden=args.n_hidden,
                    n_layers=args.n_layers)
    epochs = args.epochs or int(np.clip(round(args.target_steps * args.batch_size
                                              / adata.n_obs), 20, 400))
    steps = int(np.ceil(adata.n_obs / args.batch_size)) * epochs
    gpu = torch.cuda.is_available()
    print(f"  device: {'GPU ' + torch.cuda.get_device_name(0) if gpu else 'CPU (slow!)'}")
    print(f"  training {epochs} epochs x batch {args.batch_size} = {steps:,} steps "
          f"(KL warm-up / loss schedules over {epochs} epochs)")
    t = time.time()
    model.train(max_epochs=epochs, batch_size=args.batch_size,
                accelerator="gpu" if gpu else "cpu", devices=1,
                train_size=1.0, early_stopping=False,
                plan_kwargs={"lr": args.lr, "n_epochs_kl_warmup": epochs})
    print(f"  trained in {(time.time() - t) / 60:.1f} min")
    model.save(str(out_dir / "model"), overwrite=True, save_anndata=False)
    hist = {k: v.iloc[:, 0].tolist() for k, v in model.history.items()
            if hasattr(v, "iloc")}
    (out_dir / "training_history.json").write_text(json.dumps(hist))
    z = model.get_latent_representation(adata, batch_size=4096)
    return z, epochs


# ---------------------------------------------------------------- phases / plots
def assign_phases(df: pd.DataFrame, ref_col: str, n_bins: int, out_dir: Path) -> None:
    """cyclevi_phase = reference phase most enriched in the cell's angle bin."""
    edges = np.linspace(0, 2 * np.pi, n_bins + 1)
    df["angle_bin"] = np.clip(np.digitize(df["cyclevi_angle"], edges) - 1, 0, n_bins - 1)
    tab = pd.crosstab(df["angle_bin"], df[ref_col]).reindex(range(n_bins), fill_value=0)
    frac = tab.div(tab.sum(axis=1).replace(0, np.nan), axis=0)
    overall = df[ref_col].value_counts(normalize=True)
    enrich = frac.div(overall[frac.columns], axis=1)
    label = enrich.idxmax(axis=1)
    out = tab.add_prefix("n_").join(enrich.add_prefix("enrich_"))
    out.insert(0, "angle_from", edges[:-1]); out.insert(1, "angle_to", edges[1:])
    out["cyclevi_phase"] = label
    out.to_csv(out_dir / "phase_bins.csv")
    df["cyclevi_phase"] = df["angle_bin"].map(label)
    ct = pd.crosstab(df["cyclevi_phase"], df[ref_col])
    ct.to_csv(out_dir / f"crosstab_cyclevi_vs_{ref_col}.csv")
    print(f"\n  cyclevi_phase (rows) vs {ref_col} (columns):\n" + ct.to_string())


def _ordered(values) -> list[str]:
    vals = [v for v in PHASE_COLORS if v in set(values)]
    return vals + sorted(set(values) - set(vals))


def plots(df: pd.DataFrame, ref_col: str, marker_expr: pd.DataFrame | None,
          class_col: str, out_dir: Path, seed: int) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ang = df["cyclevi_angle"].to_numpy()
    ticks = ([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi],
             ["0", "π/2", "π", "3π/2", "2π"])
    phases = _ordered(df[ref_col].astype(str).unique())
    bins = np.linspace(0, 2 * np.pi, 73)

    # 1. angle histogram stacked by reference phase, + cyclevi_phase strip
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.hist([ang[df[ref_col].astype(str) == p] for p in phases], bins=bins, stacked=True,
            color=[PHASE_COLORS.get(p, "#777") for p in phases], label=phases)
    ax.set_xticks(*ticks); ax.set_xlim(0, 2 * np.pi)
    ax.set_xlabel("CycleVI angle"); ax.set_ylabel("Cells")
    ax.set_title(f"CycleVI angle, coloured by {ref_col}", fontweight="bold")
    ax.legend(title=ref_col, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(out_dir / f"angle_hist_by_{ref_col}.png", dpi=200)
    plt.close(fig)

    # 2. smoothed marker expression along the angle
    if marker_expr is not None and marker_expr.shape[1]:
        # mean per angle bin, then a circular moving average over bins
        nb = 120
        b = np.clip((ang / (2 * np.pi) * nb).astype(int), 0, nb - 1)
        centers = (np.arange(nb) + 0.5) * 2 * np.pi / nb
        k = 5
        fig, ax = plt.subplots(figsize=(10, 5))
        for g in marker_expr.columns:
            means = (pd.Series(marker_expr[g].to_numpy()).groupby(b).mean()
                       .reindex(range(nb)).interpolate(limit_direction="both").to_numpy())
            pad = np.r_[means[-k:], means, means[:k]]
            sm = np.convolve(pad, np.ones(2 * k + 1) / (2 * k + 1), mode="same")[k:-k]
            rng = sm.max() - sm.min()
            ax.plot(centers, (sm - sm.min()) / rng if rng > 0 else sm * 0, lw=2, label=g)
        ax.set_xticks(*ticks); ax.set_xlim(0, 2 * np.pi)
        ax.set_xlabel("CycleVI angle")
        ax.set_ylabel("Expression, log1p CP10k\n(smoothed, min-max scaled)")
        ax.set_title("Cell cycle markers along the CycleVI angle", fontweight="bold")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        # cyclevi_phase strip along the top
        edges = np.linspace(0, 2 * np.pi, 200)
        mids = (edges[:-1] + edges[1:]) / 2
        lab = (df.assign(b=np.digitize(ang, edges) - 1).groupby("b")["cyclevi_phase"]
                 .agg(lambda s: s.mode().iat[0]))
        for i, m in enumerate(mids):
            if i in lab.index:
                ax.axvspan(edges[i], edges[i + 1], ymin=0.97, ymax=1.0, lw=0,
                           color=PHASE_COLORS.get(lab[i], "#777"))
        fig.tight_layout(); fig.savefig(out_dir / "marker_curves.png", dpi=200)
        plt.close(fig)

    # 3. angle density per cell class
    fig, ax = plt.subplots(figsize=(10, 4.5))
    for c in sorted(df[class_col].astype(str).unique()):
        m = df[class_col].astype(str) == c
        h, _ = np.histogram(ang[m], bins=bins, density=True)
        ax.plot((bins[:-1] + bins[1:]) / 2, h, lw=2, label=f"{c} (n={int(m.sum()):,})")
    ax.set_xticks(*ticks); ax.set_xlim(0, 2 * np.pi)
    ax.set_xlabel("CycleVI angle"); ax.set_ylabel("Density")
    ax.set_title("CycleVI angle by cell class", fontweight="bold")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(out_dir / "angle_by_cellclass.png", dpi=200)
    plt.close(fig)

    # 4. the circular latent
    idx = np.random.default_rng(seed).permutation(len(df))[:100_000]
    d = df.iloc[idx]
    fig, ax = plt.subplots(figsize=(7, 6.5))
    for p in phases:
        m = d[ref_col].astype(str) == p
        ax.scatter(d.loc[m, "z_cycle_x"], d.loc[m, "z_cycle_y"], s=1, lw=0,
                   c=PHASE_COLORS.get(p, "#777"), label=p, rasterized=True)
    ax.set_aspect("equal"); ax.set_xlabel("z_cycle 1"); ax.set_ylabel("z_cycle 2")
    ax.set_title(f"CycleVI cell cycle latent ({len(d):,} cells)", fontweight="bold")
    ax.legend(title=ref_col, markerscale=8, bbox_to_anchor=(1.02, 1), loc="upper left",
              frameon=False)
    fig.tight_layout(); fig.savefig(out_dir / f"z_cycle_by_{ref_col}.png", dpi=200)
    plt.close(fig)

    # 5. crosstab heatmap (% of each cyclevi_phase row)
    ct = pd.crosstab(df["cyclevi_phase"], df[ref_col], normalize="index") * 100
    ct = ct.reindex(index=_ordered(ct.index), columns=_ordered(ct.columns))
    fig, ax = plt.subplots(figsize=(1.2 * ct.shape[1] + 2.5, 0.8 * ct.shape[0] + 1.5))
    im = ax.imshow(ct.to_numpy(), cmap="Blues", vmin=0, vmax=100)
    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            v = ct.iat[i, j]
            ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                    color="white" if v > 60 else "black", fontsize=10)
    ax.set_xticks(range(ct.shape[1]), ct.columns, rotation=45, ha="right")
    ax.set_yticks(range(ct.shape[0]), ct.index)
    ax.set_xlabel(ref_col); ax.set_ylabel("cyclevi_phase")
    ax.set_title("% of each CycleVI phase", fontweight="bold")
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout(); fig.savefig(out_dir / f"crosstab_cyclevi_vs_{ref_col}.png", dpi=200)
    plt.close(fig)
    print(f"  Saved plots to {out_dir}")


def umap_plot(df: pd.DataFrame, ref_col: str, umap_from: str, out_dir: Path,
              seed: int, max_background: int = 400_000) -> None:
    """The dataset UMAP (obsm['X_umap'] of `umap_from`, matched by cell ID):
    all cells light grey, the trained cells on top coloured by cyclevi_phase
    (left) and by the reference phase (right), in random order."""
    import anndata as ad
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    src = ad.read_h5ad(umap_from, backed="r")
    try:
        key = next((k for k in ("X_umap", "UMAP") if k in src.obsm), None)
        if key is None:
            print(f"  [skip] UMAP — no obsm['X_umap'] in {umap_from}")
            return
        xy_all = np.asarray(src.obsm[key])[:, :2]
        pos = src.obs_names.get_indexer(df.index.astype(str))
    finally:
        src.file.close()
    found = pos >= 0
    print(f"  UMAP: {int(found.sum()):,}/{len(df):,} trained cells matched by ID in {umap_from}")
    if not found.any():
        return
    rng = np.random.default_rng(seed)
    bg = rng.permutation(len(xy_all))[:max_background]
    d = df[found]
    xy = xy_all[pos[found]]
    order = rng.permutation(len(d))
    fig, axes = plt.subplots(1, 2, figsize=(16, 7.5))
    for ax, col, title in ((axes[0], "cyclevi_phase", "CycleVI phase"),
                           (axes[1], ref_col, f"{ref_col} (reference)")):
        ax.scatter(xy_all[bg, 0], xy_all[bg, 1], s=0.2, c="#e6e5e0", lw=0, rasterized=True)
        lab = d[col].astype(str).to_numpy()[order]
        ax.scatter(xy[order, 0], xy[order, 1], s=0.4, lw=0, rasterized=True,
                   c=[PHASE_COLORS.get(v, "#777") for v in lab])
        counts = pd.Series(lab).value_counts()
        handles = [Line2D([], [], marker="o", ls="", ms=8, color=PHASE_COLORS.get(v, "#777"),
                          label=f"{v} ({counts.get(v, 0):,})") for v in _ordered(counts.index)]
        handles.append(Line2D([], [], marker="o", ls="", ms=8, color="#e6e5e0",
                              label="not in the model"))
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False)
        ax.set_title(title, fontsize=13, fontweight="bold", loc="left")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_dir / "umap_cyclevi_phase.png", dpi=250, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_dir / 'umap_cyclevi_phase.png'}")


def marker_expression(adata) -> pd.DataFrame:
    """log1p CP10k of the marker genes present (from the raw counts)."""
    import scipy.sparse as sp
    genes = [g for g in MARKERS if g in adata.var_names]
    total = np.asarray(adata.X.sum(axis=1)).ravel()
    sub = adata[:, genes].X
    sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
    return pd.DataFrame(np.log1p(sub * (1e4 / np.where(total > 0, total, 1))[:, None]),
                        columns=genes, index=adata.obs_names)


# ---------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description="Cell cycle angle with CycleVI.")
    p.add_argument("--h5ad", required=True, help="h5ad with raw counts in X.")
    p.add_argument("--label", required=True, help="Output subfolder name.")
    p.add_argument("--cell-class", nargs="+", default=DEFAULT_CLASSES,
                   help=f"Classes to train on. Default: {DEFAULT_CLASSES}.")
    p.add_argument("--class-col", default="CellClass")
    p.add_argument("--sym-col", default="Gene")
    p.add_argument("--batch-key", default="Chemistry",
                   help="obs column used as batch covariate ('' for none).")
    p.add_argument("--chemistry", default=None, help="Keep only this --batch-key value.")
    p.add_argument("--max-cells", type=int, default=None,
                   help="Random subsample (e.g. 50000 for a quick test run).")
    p.add_argument("--n-hvg", type=int, default=4000)
    p.add_argument("--n-latent", type=int, default=10)
    p.add_argument("--n-hidden", type=int, default=128)
    p.add_argument("--n-layers", type=int, default=1)
    p.add_argument("--epochs", type=int, default=None,
                   help="Default: ~--target-steps gradient steps, 20..400 epochs.")
    p.add_argument("--target-steps", type=int, default=TARGET_STEPS)
    p.add_argument("--batch-size", type=int, default=256)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--reference-col", default="CellCyclePhase",
                   help="obs column the angle bins are named from; falls back to the "
                        "initial Seurat phase if absent.")
    p.add_argument("--n-bins", type=int, default=36)
    p.add_argument("--umap-from", default="data/full_adata_with_umap.h5ad",
                   help="h5ad whose obsm['X_umap'] is used for umap_cyclevi_phase.png "
                        "(matched by cell ID); '' to skip.")
    p.add_argument("--analyze-only", action="store_true",
                   help="Redo phases and plots from an existing cells.csv.")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    args = p.parse_args()
    args.batch_key = args.batch_key or None

    out_dir = Path(args.output_dir) / args.label
    out_dir.mkdir(parents=True, exist_ok=True)
    marker_expr = None
    if args.analyze_only:
        df = pd.read_csv(out_dir / "cells.csv", index_col=0)
        mpath = out_dir / "marker_expression.csv.gz"
        marker_expr = pd.read_csv(mpath, index_col=0) if mpath.exists() else None
    else:
        import scvi
        scvi.settings.seed = args.seed
        adata = load_cells(args)
        print("Initial phase guess (Seurat scoring)")
        lognorm = initial_phase(adata, args.seed)
        genes = select_genes(adata, lognorm, args.n_hvg, args.batch_key)
        del lognorm
        marker_expr = marker_expression(adata)
        marker_expr.to_csv(out_dir / "marker_expression.csv.gz")
        adata = adata[:, genes].copy()
        print("Training CycleVI")
        z, epochs = train(adata, args, out_dir)
        np.save(out_dir / "latent_other.npy", z[:, 2:].astype(np.float32))
        keep = [c for c in (args.class_col, "Region", "Age", args.batch_key,
                            args.reference_col, "S_score", "G2M_score", "phase",
                            "cycle_angle_uniform") if c and c in adata.obs.columns]
        df = adata.obs[list(dict.fromkeys(keep))].copy()
        df = df.rename(columns={"phase": "seurat_phase"})
        df["z_cycle_x"], df["z_cycle_y"] = z[:, 0], z[:, 1]
        df["cyclevi_angle"] = np.mod(np.arctan2(z[:, 1], z[:, 0]), 2 * np.pi)
        df["cyclevi_radius"] = np.hypot(z[:, 0], z[:, 1])
        (out_dir / "run_info.json").write_text(json.dumps({
            "h5ad": args.h5ad, "cell_class": args.cell_class, "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars), "batch_key": args.batch_key, "epochs": epochs,
            "batch_size": args.batch_size, "n_latent": args.n_latent}, indent=2))

    ref = args.reference_col if args.reference_col in df.columns else "seurat_phase"
    if ref != args.reference_col:
        print(f"  [note] no '{args.reference_col}': phases named from the Seurat phase")
    df[ref] = df[ref].astype(str).str.replace("PostM", "Post-M")
    assign_phases(df, ref, args.n_bins, out_dir)
    df.to_csv(out_dir / "cells.csv")
    plots(df, ref, marker_expr, args.class_col, out_dir, args.seed)
    if args.umap_from:
        if Path(args.umap_from).exists():
            umap_plot(df, ref, args.umap_from, out_dir, args.seed)
        else:
            print(f"  [skip] UMAP — {args.umap_from} not found")
    print(f"Outputs: {out_dir}")


if __name__ == "__main__":
    main()
