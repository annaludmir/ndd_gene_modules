"""
marker_subpop_phase.py

Split one cell class (e.g. Neuronal IPC) into subpopulations by marker gene
expression, and show what % of each subpopulation is in each cell cycle phase.

Default subpopulations (--marker NAME=GENE, repeatable):
  proliferating   = SOX2+     (SOX2 counts >= --min-counts)
  differentiating = NEUROD6+
A cell is "+" for a marker when its raw count of that gene is >= --min-counts
(default 1). The marker subpopulations can overlap, so with two markers the
table also has the exclusive groups: SOX2+ only, NEUROD6+ only, both, neither
(and "all" = the whole class).

For each subpopulation, % of its cells in each CellCyclePhase:
  - of all cells (Non-cycling included)
  - of cycling cells only (G1, S, G2M, Post-M; same denominator as
    cell_class_timeline.py --mode phase)

Outputs (results/additional_analyses/marker_subpop_phase/<label>/):
  phase_by_subpop.csv          counts and % per subpopulation x phase
                               (all regions, and per RegionsGeneral with
                               --by-region-general)
  phase_by_subpop[_<region>].png   stacked bars, one bar per subpopulation,
                                   left: of all cells, right: of cycling cells
  marker_detection.csv         % of the class's cells detecting each marker

Only the marker genes' columns are read, in row chunks (backed mode).

Usage:
  python modules/marker_subpop_phase.py \
    --h5ad data/Cortex_EMX1_louvain3_passedQC_PostM_rev1.h5ad \
    --cell-class "Neuronal IPC" --label cortex \
    --marker proliferating=SOX2 --marker differentiating=NEUROD6
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from neuron_distribution import regions_general

DEFAULT_OUTPUT_DIR = Path("results/additional_analyses/marker_subpop_phase")
PHASES = ["Non-cycling", "G1", "S", "G2M", "Post-M"]
CYCLING = PHASES[1:]
PHASE_COLORS = {"Non-cycling": "#bdbdbd", "G1": "#1f77b4", "S": "#2ca02c",
                "G2M": "#d62728", "Post-M": "#9467bd"}
CHUNK = 20000


def parse_markers(items: list[str]) -> dict[str, str]:
    out = {}
    for it in items:
        if "=" not in it:
            raise ValueError(f"--marker must be NAME=GENE, got '{it}'")
        name, gene = it.split("=", 1)
        out[name.strip()] = gene.strip()
    return out


def marker_counts(a, rows: np.ndarray, genes: list[str], sym_col: str) -> pd.DataFrame:
    """Raw counts of `genes` for `rows` (sorted), read in row chunks."""
    import scipy.sparse as sp
    sym = (a.var[sym_col].astype(str).to_numpy() if sym_col in a.var.columns
           else a.var_names.astype(str).to_numpy())
    cols = []
    for g in genes:
        hit = np.flatnonzero(sym == g)
        if len(hit) == 0:
            raise KeyError(f"Gene '{g}' not found in var['{sym_col}'] / var_names")
        if len(hit) > 1:
            print(f"  [note] {g}: {len(hit)} var entries, summing them")
        cols.append(hit)
    flat = np.concatenate(cols)
    parts = []
    for i in range(0, len(rows), CHUNK):
        X = a.X[rows[i:i + CHUNK]]
        X = X[:, flat] if sp.issparse(X) else np.asarray(X)[:, flat]
        parts.append(X.toarray() if sp.issparse(X) else np.asarray(X))
    M = np.vstack(parts)
    out, j = {}, 0
    for g, hit in zip(genes, cols):
        out[g] = M[:, j:j + len(hit)].sum(axis=1)
        j += len(hit)
    return pd.DataFrame(out, index=a.obs_names[rows])


def subpop_masks(pos: pd.DataFrame, markers: dict[str, str]) -> dict[str, np.ndarray]:
    """Marker subpopulations (overlapping), plus exclusive groups for 2 markers.
    pos: one bool column per marker name; markers: name -> gene."""
    names = list(pos.columns)
    masks = {f"{markers[n]}+ ({n})": pos[n].to_numpy() for n in names}
    if len(names) == 2:
        a, b = (pos[n].to_numpy() for n in names)
        masks[f"{markers[names[0]]}+ only"] = a & ~b
        masks[f"{markers[names[1]]}+ only"] = b & ~a
        masks["both"] = a & b
        masks["neither"] = ~a & ~b
    masks["all"] = np.ones(len(pos), dtype=bool)
    return masks


def phase_table(phase: pd.Series, masks: dict[str, np.ndarray]) -> pd.DataFrame:
    rows = []
    for name, m in masks.items():
        ph = phase[m]
        counts = ph.value_counts().reindex(PHASES, fill_value=0)
        n, n_cyc = int(m.sum()), int(counts[CYCLING].sum())
        row = {"subpopulation": name, "n_cells": n, "n_cycling": n_cyc}
        for p in PHASES:
            row[f"n_{p}"] = int(counts[p])
        for p in PHASES:
            row[f"pct_all_{p}"] = 100 * counts[p] / n if n else np.nan
        for p in CYCLING:
            row[f"pct_cycling_{p}"] = 100 * counts[p] / n_cyc if n_cyc else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def plot(tab: pd.DataFrame, title: str, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 0.6 * len(tab) + 2), sharey=True)
    y = np.arange(len(tab))
    labels = [f"{s} (n={n:,})" for s, n in zip(tab["subpopulation"], tab["n_cells"])]
    for ax, kind, phases, xlab in (
            (axes[0], "all", PHASES, "% of cells"),
            (axes[1], "cycling", CYCLING, "% of cycling cells")):
        left = np.zeros(len(tab))
        for p in phases:
            v = tab[f"pct_{kind}_{p}"].fillna(0).to_numpy()
            ax.barh(y, v, left=left, color=PHASE_COLORS[p], label=p,
                    edgecolor="white", linewidth=1)
            for yi, (l, w) in enumerate(zip(left, v)):
                if w >= 6:
                    ax.text(l + w / 2, yi, f"{w:.0f}", ha="center", va="center",
                            fontsize=9, color="white" if p != "Non-cycling" else "black")
            left += v
        ax.set_xlim(0, 100)
        ax.set_xlabel(xlab, fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_yticks(y, labels, fontsize=10)
    axes[0].invert_yaxis()
    axes[1].set_title("Cycling cells only", fontsize=11)
    axes[0].set_title("All cells", fontsize=11)
    handles, lbls = axes[0].get_legend_handles_labels()
    fig.legend(handles, lbls, title="CellCyclePhase", loc="upper left",
               bbox_to_anchor=(1.0, 0.9), fontsize=10)
    fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


def main():
    p = argparse.ArgumentParser(description="Cell cycle phases of marker-defined subpopulations.")
    p.add_argument("--h5ad", required=True, help="Cell-cycle-annotated h5ad (CellCyclePhase).")
    p.add_argument("--label", required=True, help="Output subfolder name.")
    p.add_argument("--cell-class", default="Neuronal IPC")
    p.add_argument("--class-col", default="CellClass")
    p.add_argument("--marker", action="append", default=None,
                   help="NAME=GENE, repeatable. Default proliferating=SOX2 "
                        "differentiating=NEUROD6.")
    p.add_argument("--min-counts", type=float, default=1,
                   help="Raw count of a marker for a cell to be '+'. Default 1.")
    p.add_argument("--sym-col", default="Gene")
    p.add_argument("--phase-col", default="CellCyclePhase")
    p.add_argument("--region-col", default="Region")
    p.add_argument("--by-region-general", action="store_true",
                   help="Also one table/figure per Forebrain / Midbrain / Hindbrain.")
    p.add_argument("--chemistry", default=None, help="Keep one chemistry (e.g. v3).")
    p.add_argument("--chemistry-col", default="Chemistry")
    p.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    args = p.parse_args()
    markers = parse_markers(args.marker or ["proliferating=SOX2", "differentiating=NEUROD6"])

    import anndata as ad
    print(f"Loading (backed): {args.h5ad}")
    a = ad.read_h5ad(args.h5ad, backed="r")
    obs = a.obs
    if args.phase_col not in obs.columns:
        raise KeyError(f"obs column '{args.phase_col}' not found; use a cell-cycle-annotated h5ad.")
    mask = (obs[args.class_col].astype(str) == args.cell_class).to_numpy()
    if args.chemistry:
        mask &= (obs[args.chemistry_col].astype(str) == args.chemistry).to_numpy()
    rows = np.flatnonzero(mask)
    print(f"  {len(rows):,} cells: {args.class_col} == '{args.cell_class}'"
          + (f", {args.chemistry_col} == '{args.chemistry}'" if args.chemistry else ""))
    if len(rows) == 0:
        raise ValueError(f"No cells selected. {args.class_col} values: "
                         f"{sorted(obs[args.class_col].astype(str).unique())}")

    counts = marker_counts(a, rows, list(markers.values()), args.sym_col)
    sub = obs.iloc[rows]
    phase = sub[args.phase_col].astype(str).str.replace("PostM", "Post-M").reset_index(drop=True)
    pos = pd.DataFrame({n: counts[g].to_numpy() >= args.min_counts for n, g in markers.items()})
    if args.by_region_general and args.region_col not in sub.columns:
        raise KeyError(f"--by-region-general: obs column '{args.region_col}' not found")
    region = (regions_general(sub[args.region_col]).reset_index(drop=True)
              if args.by_region_general else None)
    a.file.close()

    out_dir = Path(args.output_dir) / args.label
    out_dir.mkdir(parents=True, exist_ok=True)

    det = pd.DataFrame([{"subpopulation": n, "gene": g, "n_pos": int(pos[n].sum()),
                         "pct_pos": 100 * pos[n].mean(),
                         "mean_counts": float(counts[g].mean())}
                        for n, g in markers.items()])
    det.to_csv(out_dir / "marker_detection.csv", index=False)
    print("\n  Marker detection:")
    for r in det.itertuples():
        print(f"    {r.subpopulation:16s} {r.gene:10s} {r.n_pos:>9,} cells "
              f"({r.pct_pos:.1f}%)  mean counts {r.mean_counts:.3f}")
        if r.pct_pos < 1:
            print(f"    [warning] {r.gene} is detected in < 1% of cells; check that its "
                  f"counts are not lost in this file (e.g. gene overlap / annotation).")

    gene_desc = ", ".join(f"{n}={g}" for n, g in markers.items())
    tables = []
    scopes = [("all regions", "", np.ones(len(pos), dtype=bool))]
    if args.by_region_general:
        for r in sorted(region.dropna().unique()):
            scopes.append((r.removesuffix("_General"), f"_{r.removesuffix('_General')}",
                           (region == r).to_numpy()))
    for scope, tag, m in scopes:
        tab = phase_table(phase[m].reset_index(drop=True),
                          {k: v[m] for k, v in subpop_masks(pos, markers).items()})
        tab.insert(0, "region", scope)
        tables.append(tab)
        print(f"\n  {scope}:")
        show = ["subpopulation", "n_cells"] + [f"pct_all_{p}" for p in PHASES]
        print(tab[show].round(1).to_string(index=False))
        plot(tab, f"{args.cell_class}: cell cycle phase by marker ({scope})\n"
                  f"{gene_desc}, counts >= {args.min_counts:g}",
             out_dir / f"phase_by_subpop{tag}.png")
    pd.concat(tables).to_csv(out_dir / "phase_by_subpop.csv", index=False)
    print(f"\nOutputs: {out_dir}")


if __name__ == "__main__":
    main()
