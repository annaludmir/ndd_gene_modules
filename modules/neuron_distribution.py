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

Only obs and obsm are read (backed mode), so even large files are cheap.
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
             dataset: str, out_dir: Path) -> None:
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
        ax.set_title(f"Neuron {label_col.removeprefix('Neuron')} Proportions: {lab}",
                     fontsize=14, fontweight="bold")
        if group_col != "AllCells":
            ax.legend(title=group_col, loc="upper left", bbox_to_anchor=(1.04, 1),
                      fontsize=11, title_fontsize=11)
        fig.tight_layout()
        path = out_dir / f"timeline_{label_col}_{group_col}_{lab}.png"
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(fig)
    print(f"  Saved: timeline_{label_col}_{group_col}_<label>.png ({len(labels)} figures)")


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
        umap_from: str | None = None) -> None:
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
    args = p.parse_args()
    run(args.h5ad, args.label, args.output_dir, args.region_col, args.age_col,
        umap_from=args.umap_from)


if __name__ == "__main__":
    main()
