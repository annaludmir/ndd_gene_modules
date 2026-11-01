"""
neuron_distribution.py

How neuron labels from neuron_type_annotation.py (NeuronType, NeuronMaturity)
are distributed over age and region, and where the uncertain calls sit.

Per dataset (one h5ad), writes to
results/additional_analyses/neuron_distribution/<label>/:

  1. timeline_NeuronMaturity_<grouping>_<label>.png, + one .csv per grouping
  2. timeline_NeuronType_<grouping>_<label>.png,     + one .csv per grouping
       One figure per label (e.g. Excitatory): a line per region over age
       (weeks), y = % of that region's neurons at that age with the label.
       <grouping> is Region (as annotated) and RegionsGeneral (Forebrain /
       Midbrain / Hindbrain, as in the GES configs; cells in other regions are
       left out of it). Chemistries are pooled. Without a region column
       (e.g. the cortex) there is a single line over all cells.
  3. umap_mixed_unassigned_neurons.png
       All cells in light grey, other neurons darker grey, Mixed and
       Unassigned neurons coloured on top. Uses the stored obsm["X_umap"];
       skipped with a note when the file has none.
  4. umap_neuron_type.png
       Neurons coloured by NeuronType (Excitatory / Inhibitory / Mixed /
       Unassigned) over all cells in light grey; the smallest groups are drawn
       last so rare ones (e.g. cortex Inhibitory) stay visible.
  5. umap_neuron_maturity.png
       The same, coloured by NeuronMaturity (Immature / Intermediate /
       Mature / Unassigned).

Maturity diagnostics (is there an age trend under the class calls?):

  6. maturity_score_<grouping>.png / .csv
       The continuous neuron_maturity_score (log2 mature / immature markers)
       per age: median line per region with the interquartile range shaded.
  7. maturity_split/<subset>/
       Items 1 and 6 again within one neuron type (NeuronType_Excitatory,
       NeuronType_Inhibitory) or one chemistry (Chemistry_v2, Chemistry_v3),
       to tell an age trend apart from changes in neuron-type or chemistry mix.
  8. marker_expression_<grouping>.png / .csv
       Each maturity marker gene (immature row, mature row): mean
       log1p(counts per 10k) over neurons per age, a line per region; the CSV
       also has the % of neurons expressing it. Reads X in row chunks (raw
       counts expected, as in neuron_type_annotation.py).

Items 1–7 read only obs and obsm (backed mode); item 8 reads X once.
If the file has no UMAP, --umap-from takes obsm["X_umap"] from another h5ad,
matched by cell ID (cells missing there are left out of the UMAPs).

Usage:
  python modules/neuron_distribution.py \
    --h5ad data/human_dev_without_week_5_cc_annotated_neuron_type.h5ad --label human_dev
"""

import argparse
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

from cell_cycle_annotation import _symbol_to_var_map
from neuron_type_annotation import IMMATURE_MARKERS, MATURE_MARKERS, MATURITY_SCORE_COL


DEFAULT_OUTPUT_DIR = Path("results/additional_analyses/neuron_distribution")

# Same grouping as the RegionsGeneral rules in the GES configs
REGIONS_GENERAL = {
    "Forebrain_General": ["Forebrain", "Telencephalon", "Diencephalon"],
    "Midbrain_General":  ["Midbrain"],
    "Hindbrain_General": ["Hindbrain", "Cerebellum", "Pons", "Medulla"],
}

LABELS = {
    "NeuronMaturity": {
        "Immature":     "#9ec5f4",
        "Intermediate": "#3987e5",
        "Mature":       "#104281",
        "Unassigned":   "#c8c7c2",
    },
    "NeuronType": {
        "Excitatory": "#2a78d6",
        "Inhibitory": "#e34948",
        "Mixed":      "#eda100",
        "Unassigned": "#c8c7c2",
    },
}

UMAP_HIGHLIGHT = {"Mixed": "#eda100", "Unassigned": "#4a3aa7"}   # violet: not reused by the timelines


def regions_general(region: pd.Series) -> pd.Series:
    lookup = {r: g for g, rs in REGIONS_GENERAL.items() for r in rs}
    return region.astype(str).map(lookup)


def age_label(a) -> str:
    """6.599999904 -> "6.6", 7.0 -> "7"."""
    return f"{round(float(a), 1):g}"


def timeline(obs: pd.DataFrame, label_col: str, group_col: str, age_col: str,
             dataset: str, out_dir: Path, subset: str | None = None) -> None:
    """For each label, one line plot: % of each region's neurons with that
    label at each age (one line per region)."""
    labels = list(LABELS[label_col])
    df = obs.loc[obs[label_col].notna() & obs[group_col].notna(), [group_col, age_col, label_col]].copy()
    if df.empty:
        print(f"  [skip] {label_col} × {group_col}: no labelled neurons with a region")
        return
    df[age_col] = pd.to_numeric(df[age_col], errors="coerce").round(1)
    df = df.dropna(subset=[age_col])
    df[label_col] = df[label_col].astype(str)
    df[group_col] = df[group_col].astype(str)

    counts = (df.groupby([group_col, age_col, label_col], observed=True).size()
                .unstack(label_col, fill_value=0)
                .reindex(columns=labels, fill_value=0))
    totals = counts.sum(axis=1)
    pct = counts.div(totals, axis=0) * 100
    table = counts.join(totals.rename("n_neurons")).join(pct.add_prefix("pct_"))
    table.to_csv(out_dir / f"timeline_{label_col}_{group_col}.csv")

    groups = sorted(df[group_col].unique())
    legend_name = {g: g.removesuffix("_General") for g in groups}
    ages = sorted(df[age_col].unique())
    lo, hi = np.floor(min(ages) * 2) / 2, np.ceil(max(ages) * 2) / 2

    for lab in labels:
        fig, ax = plt.subplots(figsize=(10, 5))
        for g in groups:
            series = pct.loc[g, lab].sort_index()
            ax.plot(series.index, series.values, marker="o", markersize=6, linewidth=2,
                    markeredgecolor="white", label=legend_name[g])
        ax.set_xticks(np.arange(lo, hi + 0.01, 0.5))
        ax.tick_params(axis="x", labelrotation=45)
        ax.grid(True, linestyle="--", alpha=0.5)
        ax.set_xlabel(f"{age_col} (PCW)", fontsize=12)
        ax.set_ylabel("Percentage of neurons (%)", fontsize=12)
        ax.set_title(f"Neuron {label_col.removeprefix('Neuron')} Proportions: {lab}"
                     + (f" ({subset})" if subset else ""), fontsize=14, fontweight="bold")
        if group_col != "AllCells":
            ax.legend(title=group_col, loc="upper left", bbox_to_anchor=(1.04, 1),
                      fontsize=11, title_fontsize=11)
        fig.tight_layout()
        path = out_dir / f"timeline_{label_col}_{group_col}_{lab}.png"
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(fig)
    print(f"  Saved: timeline_{label_col}_{group_col}_<label>.png ({len(labels)} figures)")


def _age_ticks(ax, ages) -> None:
    lo, hi = np.floor(min(ages) * 2) / 2, np.ceil(max(ages) * 2) / 2
    ax.set_xticks(np.arange(lo, hi + 0.01, 0.5))
    ax.tick_params(axis="x", labelrotation=45)


def maturity_score_by_age(obs: pd.DataFrame, group_col: str, age_col: str,
                          out_dir: Path, subset: str | None = None) -> None:
    """Median neuron_maturity_score per age (IQR shaded), a line per region."""
    if MATURITY_SCORE_COL not in obs:
        print(f"  [skip] maturity score plot — no '{MATURITY_SCORE_COL}' column")
        return
    df = obs.loc[obs[group_col].notna(), [group_col, age_col]].copy()
    df["score"] = pd.to_numeric(obs[MATURITY_SCORE_COL], errors="coerce")
    df[age_col] = pd.to_numeric(df[age_col], errors="coerce").round(1)
    df = df.dropna()
    if df.empty:
        print(f"  [skip] maturity score × {group_col}: no scored neurons")
        return
    df[group_col] = df[group_col].astype(str)
    stats = (df.groupby([group_col, age_col])["score"]
               .agg(n_neurons="size", median="median", mean="mean",
                    q25=lambda v: v.quantile(0.25), q75=lambda v: v.quantile(0.75)))
    stats.to_csv(out_dir / f"maturity_score_{group_col}.csv")

    fig, ax = plt.subplots(figsize=(10, 5))
    for g in sorted(df[group_col].unique()):
        st = stats.loc[g].sort_index()
        line, = ax.plot(st.index, st["median"], marker="o", markersize=6, linewidth=2,
                        markeredgecolor="white", label=g.removesuffix("_General"))
        ax.fill_between(st.index, st["q25"], st["q75"], color=line.get_color(), alpha=0.12,
                        linewidth=0)
    ax.axhline(0, color="#8f8e88", linewidth=1)
    _age_ticks(ax, df[age_col].unique())
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.set_xlabel(f"{age_col} (PCW)", fontsize=12)
    ax.set_ylabel("Maturity score, log2(mature / immature)\nmedian, IQR shaded", fontsize=12)
    ax.set_title("Neuron maturity score" + (f" ({subset})" if subset else ""),
                 fontsize=14, fontweight="bold")
    if group_col != "AllCells":
        ax.legend(title=group_col, loc="upper left", bbox_to_anchor=(1.04, 1),
                  fontsize=11, title_fontsize=11)
    fig.tight_layout()
    fig.savefig(out_dir / f"maturity_score_{group_col}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: maturity_score_{group_col}.png" + (f" ({subset})" if subset else ""))


def maturity_split(obs: pd.DataFrame, groupings, age_col: str, dataset: str,
                   out_dir: Path, chemistry_col: str) -> None:
    """Maturity timelines and score plots within each neuron type and chemistry."""
    subsets = {}
    if "NeuronType" in obs:
        for t in ("Excitatory", "Inhibitory"):
            subsets[f"NeuronType_{t}"] = obs["NeuronType"].astype(object) == t
    if chemistry_col in obs:
        for c in sorted(obs.loc[obs["NeuronMaturity"].notna(), chemistry_col].dropna().astype(str).unique()):
            subsets[f"{chemistry_col}_{c}"] = obs[chemistry_col].astype(str) == c
    else:
        print(f"  [note] no '{chemistry_col}' column — maturity not split by chemistry")
    for name, mask in subsets.items():
        sub = obs[mask.to_numpy() & obs["NeuronMaturity"].notna().to_numpy()]
        print(f"  maturity_split/{name}: {len(sub):,} neurons")
        if sub.empty:
            continue
        sub_dir = out_dir / "maturity_split" / name
        sub_dir.mkdir(parents=True, exist_ok=True)
        for group_col in groupings:
            timeline(sub, "NeuronMaturity", group_col, age_col, dataset, sub_dir, subset=name)
            maturity_score_by_age(sub, group_col, age_col, sub_dir, subset=name)


def marker_expression_by_age(a, obs: pd.DataFrame, groupings, age_col: str,
                             sym_col: str, out_dir: Path, chunk: int = 20000) -> None:
    """Mean log1p(CP10k) and % expressing of each maturity marker per age and
    region, over neurons (cells with a NeuronMaturity call)."""
    neuron = obs["NeuronMaturity"].notna().to_numpy()
    sym2var = _symbol_to_var_map(a, sym_col)
    var_pos = {v: i for i, v in enumerate(a.var_names.astype(str))}
    genes = [g for g in IMMATURE_MARKERS + MATURE_MARKERS if sym2var.get(g) in var_pos]
    missing = [g for g in IMMATURE_MARKERS + MATURE_MARKERS if g not in genes]
    if missing:
        print(f"  [note] marker genes not in the file: {missing}")
    if not genes or not neuron.any():
        print("  [skip] marker expression — no marker genes or no neurons")
        return
    gene_ix = np.array([var_pos[sym2var[g]] for g in genes])

    print(f"  Reading X for {int(neuron.sum()):,} neurons × {len(genes)} marker genes ...")
    parts = []
    for start in range(0, a.n_obs, chunk):
        end = min(start + chunk, a.n_obs)
        keep = neuron[start:end]
        if not keep.any():
            continue
        X = a.X[start:end]
        total = np.asarray(X.sum(axis=1)).ravel()[keep]
        sub = X[:, gene_ix]
        sub = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
        sub = sub[keep]
        parts.append(np.log1p(sub * (1e4 / np.where(total > 0, total, 1.0))[:, None]))
    expr = pd.DataFrame(np.vstack(parts).astype(np.float32), columns=genes,
                        index=obs.index[neuron])

    base = obs.loc[neuron].copy()
    base[age_col] = pd.to_numeric(base[age_col], errors="coerce").round(1)
    for group_col in groupings:
        df = base[[group_col, age_col]].join(expr)
        df = df.dropna(subset=[group_col, age_col])
        if df.empty:
            continue
        df[group_col] = df[group_col].astype(str)
        g = df.groupby([group_col, age_col])
        mean = g[genes].mean()
        pct = g[genes].agg(lambda v: 100 * (v > 0).mean())
        table = pd.concat({"mean_log1p_cp10k": mean.stack(),
                           "pct_expressing": pct.stack()}, axis=1)
        table.index.names = [group_col, age_col, "gene"]
        table = table.join(g.size().rename("n_neurons"))
        table.to_csv(out_dir / f"marker_expression_{group_col}.csv")

        regions = sorted(df[group_col].unique())
        ncol = max(len(IMMATURE_MARKERS), len(MATURE_MARKERS))
        fig, axes = plt.subplots(2, ncol, figsize=(3.2 * ncol, 6.4), squeeze=False)
        for r, (kind, markers) in enumerate((("immature", IMMATURE_MARKERS),
                                             ("mature", MATURE_MARKERS))):
            for c in range(ncol):
                ax = axes[r, c]
                gene = markers[c] if c < len(markers) else None
                if gene not in genes:
                    ax.set_visible(False)
                    continue
                for reg in regions:
                    st = mean.loc[reg, gene].sort_index()
                    ax.plot(st.index, st.values, marker="o", markersize=3, linewidth=1.5,
                            label=reg.removesuffix("_General"))
                ax.set_title(f"{gene} ({kind})", fontsize=10)
                ax.grid(True, linestyle="--", alpha=0.4)
                ax.tick_params(labelsize=8)
                ax.spines[["top", "right"]].set_visible(False)
                if c == 0:
                    ax.set_ylabel("mean log1p(CP10k)", fontsize=9)
                if r == 1:
                    ax.set_xlabel(f"{age_col} (PCW)", fontsize=9)
        handles, names = next(ax for ax in axes.ravel() if ax.get_visible()).get_legend_handles_labels()
        if group_col != "AllCells":
            fig.legend(handles, names, title=group_col, loc="lower center",
                       ncol=min(len(names), 8), frameon=False, bbox_to_anchor=(0.5, -0.06))
        fig.suptitle("Maturity marker expression in neurons by age", fontsize=13, fontweight="bold")
        fig.tight_layout()
        fig.savefig(out_dir / f"marker_expression_{group_col}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: marker_expression_{group_col}.png")


UMAP_KEYS = ("X_umap", "UMAP")


def load_umap(a, umap_from: str | None):
    """(xy, labels) for the cells that have UMAP coordinates, or None; labels
    holds the NeuronType / NeuronMaturity columns the file has. Coordinates
    come from the file itself or, if it has none, from obsm of `umap_from`,
    matched by cell ID; cells absent there are left out."""
    cols = [c for c in LABELS if c in a.obs]
    if not cols:
        print(f"  [skip] UMAPs — no {' / '.join(LABELS)} column")
        return None
    labels = a.obs[cols].astype(object)
    key = next((k for k in UMAP_KEYS if k in a.obsm), None)
    if key is not None:
        return np.asarray(a.obsm[key])[:, :2], labels
    if not umap_from:
        print("  [skip] UMAPs — no obsm['X_umap'] in this file (pass --umap-from)")
        return None
    src = ad.read_h5ad(umap_from, backed="r")
    try:
        key = next((k for k in UMAP_KEYS if k in src.obsm), None)
        if key is None:
            print(f"  [skip] UMAPs — no obsm['X_umap'] in {umap_from} either")
            return None
        pos = src.obs_names.get_indexer(a.obs_names)
        found = pos >= 0
        print(f"  UMAP from {umap_from} (obsm['{key}']): "
              f"{int(found.sum()):,}/{a.n_obs:,} cells matched by ID"
              + (f"; {int((~found).sum()):,} left out" if (~found).any() else ""))
        if not found.any():
            return None
        xy = np.asarray(src.obsm[key])[:, :2][pos[found]]
        return xy, labels[found]
    finally:
        src.file.close()


def umap_uncertain(xy: np.ndarray, nt: np.ndarray, dataset: str, out_dir: Path) -> None:
    """All cells grey, other neurons darker grey, Mixed / Unassigned on top."""
    is_neuron = pd.notna(nt)

    fig, ax = plt.subplots(figsize=(7.5, 7))
    ax.scatter(xy[~is_neuron, 0], xy[~is_neuron, 1], s=0.3, c="#e6e5e0",
               linewidths=0, rasterized=True)
    other = is_neuron & ~np.isin(nt, list(UMAP_HIGHLIGHT))
    ax.scatter(xy[other, 0], xy[other, 1], s=0.3, c="#a3a29c", linewidths=0, rasterized=True)
    handles = [Patch(color="#e6e5e0", label=f"non-neurons ({int((~is_neuron).sum()):,})"),
               Patch(color="#a3a29c", label=f"Excitatory / Inhibitory neurons ({int(other.sum()):,})")]
    for lab, col in UMAP_HIGHLIGHT.items():
        m = nt == lab
        ax.scatter(xy[m, 0], xy[m, 1], s=3, c=col, linewidths=0, rasterized=True)
        handles.append(Patch(color=col, label=f"{lab} neurons ({int(m.sum()):,})"))
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.06),
              ncol=2, frameon=False, fontsize=9)
    ax.set_title(f"{dataset}: Mixed and Unassigned neurons", fontsize=12)
    path = out_dir / "umap_mixed_unassigned_neurons.png"
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path.name}")


UMAP_LABEL_PLOTS = {   # label column -> (file name, title)
    "NeuronType":     ("umap_neuron_type.png",     "excitatory / inhibitory neurons"),
    "NeuronMaturity": ("umap_neuron_maturity.png", "neuron maturity"),
}


def umap_by_label(xy: np.ndarray, nt: np.ndarray, label_col: str,
                  dataset: str, out_dir: Path) -> None:
    """Neurons coloured by `label_col` over all cells in light grey."""
    is_neuron = pd.notna(nt)
    # Timeline colours, but Unassigned darker so it shows against the background
    colors = dict(LABELS[label_col], Unassigned="#8f8e88")

    fig, ax = plt.subplots(figsize=(7.5, 7))
    ax.scatter(xy[~is_neuron, 0], xy[~is_neuron, 1], s=0.3, c="#e6e5e0",
               linewidths=0, rasterized=True)
    handles = [Patch(color="#e6e5e0", label=f"non-neurons ({int((~is_neuron).sum()):,})")]
    counts = {lab: int((nt == lab).sum()) for lab in colors}
    for lab in sorted(colors, key=lambda l: -counts[l]):        # largest first, rarest on top
        m = nt == lab
        size = 0.6 if counts[lab] > 20000 else 3
        ax.scatter(xy[m, 0], xy[m, 1], s=size, c=colors[lab], linewidths=0, rasterized=True)
    for lab in colors:                                          # legend in a fixed order
        handles.append(Patch(color=colors[lab], label=f"{lab} ({counts[lab]:,})"))
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.06),
              ncol=3, frameon=False, fontsize=9)
    file_name, title = UMAP_LABEL_PLOTS[label_col]
    ax.set_title(f"{dataset}: {title}", fontsize=12)
    path = out_dir / file_name
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path.name}")


def run(h5ad: str, label: str, output_dir: str = str(DEFAULT_OUTPUT_DIR),
        region_col: str = "Region", age_col: str = "Age",
        umap_from: str | None = None, chemistry_col: str = "Chemistry",
        sym_col: str = "Gene", marker_expression: bool = True) -> None:
    out_dir = Path(output_dir) / label
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n{label}: {h5ad}")
    a = ad.read_h5ad(h5ad, backed="r")
    try:
        obs = a.obs.copy()
        print(f"  {a.n_obs:,} cells; neurons with NeuronType: "
              f"{int(obs['NeuronType'].notna().sum()) if 'NeuronType' in obs else 0:,}")

        if age_col not in obs:
            print(f"  [skip] timelines — '{age_col}' not in obs")
        else:
            if region_col in obs:
                obs["RegionsGeneral"] = regions_general(obs[region_col])
                groupings = (region_col, "RegionsGeneral")
            else:
                # e.g. a single-region dataset like the cortex: one panel, all cells
                print(f"  [note] no '{region_col}' column — timelines over all cells")
                obs["AllCells"] = label
                groupings = ("AllCells",)
            for label_col in ("NeuronMaturity", "NeuronType"):
                if label_col not in obs:
                    print(f"  [skip] {label_col} timelines — column missing (rerun the annotation)")
                    continue
                for group_col in groupings:
                    timeline(obs, label_col, group_col, age_col, label, out_dir)
            if "NeuronMaturity" in obs:
                for group_col in groupings:
                    maturity_score_by_age(obs, group_col, age_col, out_dir)
                maturity_split(obs, groupings, age_col, label, out_dir, chemistry_col)
                if marker_expression:
                    marker_expression_by_age(a, obs, groupings, age_col, sym_col, out_dir)

        umap = load_umap(a, umap_from)
        if umap is not None:
            xy, labels = umap
            if "NeuronType" in labels:
                umap_uncertain(xy, labels["NeuronType"].to_numpy(), label, out_dir)
            for label_col in UMAP_LABEL_PLOTS:
                if label_col in labels:
                    umap_by_label(xy, labels[label_col].to_numpy(), label_col, label, out_dir)
                else:
                    print(f"  [skip] {label_col} UMAP — column missing (rerun the annotation)")
    finally:
        a.file.close()
    print(f"  Outputs: {out_dir}")


def main():
    p = argparse.ArgumentParser(description="Neuron type / maturity distribution by age, region and on the UMAP.")
    p.add_argument("--h5ad", required=True, help="h5ad with NeuronType / NeuronMaturity.")
    p.add_argument("--label", required=True, help="Dataset name; also the output subfolder.")
    p.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    p.add_argument("--region-col", default="Region")
    p.add_argument("--age-col", default="Age")
    p.add_argument("--umap-from", default=None,
                   help="h5ad to take obsm['X_umap'] from (matched by cell ID) when "
                        "--h5ad has none, e.g. data/full_adata_with_umap.h5ad.")
    p.add_argument("--chemistry-col", default="Chemistry",
                   help="obs column to split the maturity diagnostics by (skipped if absent).")
    p.add_argument("--sym-col", default="Gene", help="var column with gene symbols.")
    p.add_argument("--no-marker-expression", action="store_true",
                   help="Skip the marker-gene expression plots (the only step that reads X).")
    args = p.parse_args()
    run(args.h5ad, args.label, args.output_dir, args.region_col, args.age_col,
        umap_from=args.umap_from, chemistry_col=args.chemistry_col, sym_col=args.sym_col,
        marker_expression=not args.no_marker_expression)


if __name__ == "__main__":
    main()
