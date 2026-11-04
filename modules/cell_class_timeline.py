"""
cell_class_timeline.py

What fraction of all cells belongs to a cell class (e.g. Neuronal IPC), per
region and age. Same layout as the cell cycle phase and neuron timelines: one
figure per class, a line per region, x = age (PCW), y = % of all of that
region's cells at that age that are of the class.

Only obs is read (backed mode), so it is quick on the full h5ad.

Outputs (results/additional_analyses/cell_class_timeline/<label>/):
  timeline_<class>_<grouping>.png   one per class and grouping
  timeline_CellClass_<grouping>.csv counts, total cells and % per region x age

<grouping> is RegionsGeneral (Forebrain / Midbrain / Hindbrain, as in the GES
configs) and, with --by-region, Region as annotated.

Usage:
  python modules/cell_class_timeline.py --h5ad data/human_dev.h5ad \
    --cell-class "Neuronal IPC" --label human_dev
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gene_expression_summary import sanitize_name
from neuron_distribution import _age_ticks, regions_general

DEFAULT_OUTPUT_DIR = Path("results/additional_analyses/cell_class_timeline")
GROUP_COLORS = {"Forebrain": "#1f77b4", "Hindbrain": "#ff7f0e", "Midbrain": "#2ca02c"}


def class_timeline(obs: pd.DataFrame, classes: list[str], class_col: str, group_col: str,
                   age_col: str, out_dir: Path, subtitle: str = "") -> None:
    """% of all cells of each region x age that are of each class."""
    df = obs.loc[obs[group_col].notna(), [group_col, age_col, class_col]].copy()
    df[age_col] = pd.to_numeric(df[age_col], errors="coerce").round(1)
    df = df.dropna(subset=[age_col])
    df[group_col] = df[group_col].astype(str).str.removesuffix("_General")
    df[class_col] = df[class_col].astype(str)

    totals = df.groupby([group_col, age_col]).size().rename("n_cells")
    counts = (df[df[class_col].isin(classes)]
              .groupby([group_col, age_col, class_col]).size()
              .unstack(class_col, fill_value=0)
              .reindex(index=totals.index, columns=classes, fill_value=0))
    pct = counts.div(totals, axis=0) * 100
    counts.add_prefix("n_").join(totals).join(pct.add_prefix("pct_")).to_csv(
        out_dir / f"timeline_{class_col}_{group_col}.csv")

    groups = sorted(df[group_col].unique())
    ages = sorted(df[age_col].unique())
    for cls in classes:
        fig, ax = plt.subplots(figsize=(10, 5))
        for g in groups:
            series = pct.loc[g, cls].sort_index()
            ax.plot(series.index, series.values, marker="o", markersize=6, linewidth=2,
                    markeredgecolor="white", label=g, color=GROUP_COLORS.get(g))
        _age_ticks(ax, ages)
        ax.set_ylim(bottom=0)
        ax.grid(True, linestyle="--", alpha=0.5)
        ax.set_xlabel(f"{age_col} (PCW)", fontsize=12)
        ax.set_ylabel("Percentage of all cells (%)", fontsize=12)
        ax.set_title(f"Cell Class Proportions: {cls}" + (f" ({subtitle})" if subtitle else ""),
                     fontsize=14, fontweight="bold")
        ax.legend(title=group_col, loc="upper left", bbox_to_anchor=(1.04, 1),
                  fontsize=11, title_fontsize=11)
        fig.tight_layout()
        path = out_dir / f"timeline_{sanitize_name(cls)}_{group_col}.png"
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {path}")


def main():
    p = argparse.ArgumentParser(description="% of all cells of a cell class per region and age.")
    p.add_argument("--h5ad", required=True)
    p.add_argument("--label", required=True, help="Output subfolder name, e.g. human_dev.")
    p.add_argument("--cell-class", nargs="+", default=["Neuronal IPC"],
                   help="Classes to plot, one figure each. Default 'Neuronal IPC'.")
    p.add_argument("--class-col", default="CellClass")
    p.add_argument("--region-col", default="Region")
    p.add_argument("--age-col", default="Age")
    p.add_argument("--by-region", action="store_true",
                   help="Also plot per Region as annotated, not only RegionsGeneral.")
    p.add_argument("--chemistry", default=None,
                   help="Keep one chemistry (e.g. v3). Default: all cells.")
    p.add_argument("--chemistry-col", default="Chemistry")
    p.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    args = p.parse_args()

    import anndata as ad
    print(f"Loading obs: {args.h5ad}")
    obs = ad.read_h5ad(args.h5ad, backed="r").obs.copy()
    print(f"  {len(obs):,} cells")
    if args.chemistry:
        obs = obs[obs[args.chemistry_col].astype(str) == args.chemistry]
        print(f"  {args.chemistry_col} == '{args.chemistry}': {len(obs):,} cells")
    found = set(obs[args.class_col].astype(str))
    missing = [c for c in args.cell_class if c not in found]
    if missing:
        raise ValueError(f"Not in obs['{args.class_col}']: {missing}. Values: {sorted(found)}")

    obs["RegionsGeneral"] = regions_general(obs[args.region_col])
    out_dir = Path(args.output_dir) / args.label
    out_dir.mkdir(parents=True, exist_ok=True)
    subtitle = f"{args.chemistry_col} {args.chemistry}" if args.chemistry else ""
    groupings = ["RegionsGeneral"] + ([args.region_col] if args.by_region else [])
    for group_col in groupings:
        class_timeline(obs, args.cell_class, args.class_col, group_col, args.age_col,
                       out_dir, subtitle)
    print(f"Outputs: {out_dir}")


if __name__ == "__main__":
    main()
