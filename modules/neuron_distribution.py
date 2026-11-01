"""
neuron_distribution.py

How neuron labels from neuron_type_annotation.py (NeuronType, NeuronMaturity)
are distributed over age and region, and where the uncertain calls sit.

Per dataset (one h5ad), writes to
results/additional_analyses/neuron_distribution/<label>/:

  1. timeline_NeuronMaturity_<grouping>.png / .csv
  2. timeline_NeuronType_<grouping>.png / .csv
       One panel per region, a stacked bar per age: the share of each label
       among that region's neurons at that age, with the neuron count above
       the bar. <grouping> is Region (as annotated) and RegionsGeneral
       (Forebrain / Midbrain / Hindbrain _General, as in the GES configs;
       cells in other regions are left out of it). Chemistries are pooled.
       Skipped when the dataset has no region column.
  3. umap_mixed_unassigned_neurons.png
       All cells in light grey, other neurons darker grey, Mixed and
       Unassigned neurons coloured on top. Uses the stored obsm["X_umap"];
       skipped with a note when the file has none.
  4. umap_neuron_type.png
       Neurons coloured by NeuronType (Excitatory / Inhibitory / Mixed /
       Unassigned) over all cells in light grey; the smallest groups are drawn
       last so rare ones (e.g. cortex Inhibitory) stay visible.

Only obs and obsm are read (backed mode), so even large files are cheap.

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
    """Stacked share of each label per age, one panel per region."""
    colors = LABELS[label_col]
    df = obs.loc[obs[label_col].notna() & obs[group_col].notna(), [group_col, age_col, label_col]].copy()
    if df.empty:
        print(f"  [skip] {label_col} × {group_col}: no labelled neurons with a region")
        return
    df[age_col] = pd.to_numeric(df[age_col], errors="coerce").round(1)
    df = df.dropna(subset=[age_col])
    df[label_col] = df[label_col].astype(str)

    counts = (df.groupby([group_col, age_col, label_col], observed=True).size()
                .unstack(label_col, fill_value=0)
                .reindex(columns=list(colors), fill_value=0))
    counts.to_csv(out_dir / f"timeline_{label_col}_{group_col}.csv")

    groups = [g for g in (list(REGIONS_GENERAL) if group_col == "RegionsGeneral"
                          else sorted(df[group_col].astype(str).unique()))
              if g in counts.index.get_level_values(0)]
    ages = sorted(df[age_col].unique())
    ncols = min(3, len(groups))
    nrows = int(np.ceil(len(groups) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 3.4 * nrows),
                             squeeze=False, sharey=True)
    for ax, g in zip(axes.flat, groups):
        sub = counts.loc[g].reindex(ages, fill_value=0)
        totals = sub.sum(axis=1)
        frac = sub.div(totals.replace(0, np.nan), axis=0).fillna(0)
        x = np.arange(len(ages))
        bottom = np.zeros(len(ages))
        for lab, col in colors.items():
            ax.bar(x, frac[lab].to_numpy(), bottom=bottom, color=col, width=0.8,
                   edgecolor="white", linewidth=0.5)
            bottom += frac[lab].to_numpy()
        for xi, n in zip(x, totals):
            if n:
                ax.text(xi, 1.01, f"{int(n):,}", ha="center", va="bottom",
                        fontsize=6, rotation=90, color="#52514e")
        ax.set_title(str(g), fontsize=10, pad=26)
        ax.set_xticks(x)
        ax.set_xticklabels([age_label(a) for a in ages], fontsize=7, rotation=90)
        ax.set_ylim(0, 1)
        ax.set_xlabel(f"{age_col} (weeks post-conception)", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes.flat[len(groups):]:
        ax.set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("share of neurons", fontsize=9)

    fig.legend(handles=[Patch(color=c, label=l) for l, c in colors.items()],
               loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=len(colors),
               frameon=False, fontsize=9)
    fig.suptitle(f"{dataset}: {label_col} by age and {group_col} "
                 f"(chemistries v2 + v3 pooled; numbers = neurons per bar)", fontsize=11)
    fig.tight_layout()
    path = out_dir / f"timeline_{label_col}_{group_col}.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path.name}")


def umap_uncertain(a, dataset: str, out_dir: Path) -> None:
    """All cells grey, other neurons darker grey, Mixed / Unassigned on top."""
    key = next((k for k in ("X_umap", "UMAP") if k in a.obsm), None)
    if key is None:
        print("  [skip] UMAP — no obsm['X_umap'] in this file")
        return
    if "NeuronType" not in a.obs:
        print("  [skip] UMAP — no NeuronType column")
        return
    xy = np.asarray(a.obsm[key])[:, :2]
    nt = a.obs["NeuronType"].astype(object).to_numpy()
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


def umap_neuron_type(a, dataset: str, out_dir: Path) -> None:
    """Neurons coloured by NeuronType over all cells in light grey."""
    key = next((k for k in ("X_umap", "UMAP") if k in a.obsm), None)
    if key is None or "NeuronType" not in a.obs:
        print("  [skip] neuron-type UMAP — no obsm['X_umap'] or NeuronType column")
        return
    xy = np.asarray(a.obsm[key])[:, :2]
    nt = a.obs["NeuronType"].astype(object).to_numpy()
    is_neuron = pd.notna(nt)
    # Timeline colours, but Unassigned darker so it shows against the background
    colors = dict(LABELS["NeuronType"], Unassigned="#8f8e88")

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
    ax.set_title(f"{dataset}: excitatory / inhibitory neurons", fontsize=12)
    path = out_dir / "umap_neuron_type.png"
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path.name}")


def run(h5ad: str, label: str, output_dir: str = str(DEFAULT_OUTPUT_DIR),
        region_col: str = "Region", age_col: str = "Age") -> None:
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

        umap_uncertain(a, label, out_dir)
        umap_neuron_type(a, label, out_dir)
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
    args = p.parse_args()
    run(args.h5ad, args.label, args.output_dir, args.region_col, args.age_col)


if __name__ == "__main__":
    main()
