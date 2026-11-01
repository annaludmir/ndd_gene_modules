"""
compare_cell_cycle_annotation.py

Compare the cell-cycle annotation of two h5ad files on the cells they share,
to explain why their CellCyclePhase calls differ:
  1. how each annotation was made (uns["cell_cycle_scoring"])
  2. gene sets: total genes, overlap, and how many G1 / S / G2M genes each has
  3. total UMIs of shared cells (the phase-fraction denominator)
  4. correlation of the per-cell scores (cycling score and phase scores)
  5. a phase cross-table (% of each file-A phase going to each file-B phase)

Files are opened in backed mode; only shared cells' rows of X are read, in a
sample, for the UMI totals.

Usage:
  python modules/compare_cell_cycle_annotation.py --a cortex.h5ad --b human_dev.h5ad \
      --out results/additional_analyses/cell_cycle_comparison
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

import cell_cycle_annotation as C


PHASE_SETS = {"G1": C.G1_GENES, "S": C.S_GENES, "G2M": C.G2M_GENES}
# Per-cell scores to correlate when present in both files
SCORE_COLS = ["cc_score_g1", "cc_score_s", "cc_score_g2m"]
CYCLING_COLS = ["cell_cycle_score", "CellCycleFraction", "cycling_score"]


def symbols(a) -> set[str]:
    return set(a.var["Gene"].astype(str)) if "Gene" in a.var else set(a.var_names.astype(str))


def main():
    p = argparse.ArgumentParser(description="Compare cell-cycle annotation of two h5ad files.")
    p.add_argument("--a", required=True, help="First h5ad (rows of the cross-table).")
    p.add_argument("--b", required=True, help="Second h5ad (columns of the cross-table).")
    p.add_argument("--label-a", default="A")
    p.add_argument("--label-b", default="B")
    p.add_argument("--phase-col", default="CellCyclePhase")
    p.add_argument("--n-umi-cells", type=int, default=2000,
                   help="Shared cells sampled for the total-UMI comparison.")
    p.add_argument("--out", default=None, help="Folder for the CSV outputs (optional).")
    args = p.parse_args()
    la, lb = args.label_a, args.label_b

    a = ad.read_h5ad(args.a, backed="r")
    b = ad.read_h5ad(args.b, backed="r")
    try:
        common = a.obs_names.intersection(b.obs_names)
        print(f"{la}: {a.n_obs:,} cells × {a.n_vars:,} genes  ({args.a})")
        print(f"{lb}: {b.n_obs:,} cells × {b.n_vars:,} genes  ({args.b})")
        print(f"shared cells: {len(common):,}")

        print("\n== 1. How each annotation was made (uns['cell_cycle_scoring'])")
        print(f"  {la}: {a.uns.get('cell_cycle_scoring', '(none: not made by cell_cycle_annotation.py)')}")
        print(f"  {lb}: {b.uns.get('cell_cycle_scoring', '(none: not made by cell_cycle_annotation.py)')}")

        print("\n== 2. Gene sets")
        ga, gb = symbols(a), symbols(b)
        print(f"  {la}: {len(ga):,} | {lb}: {len(gb):,} | shared: {len(ga & gb):,} | "
              f"only {la}: {len(ga - gb):,} | only {lb}: {len(gb - ga):,}")
        for name, genes in PHASE_SETS.items():
            genes = set(genes)
            print(f"  {name:4s} genes found: {la} {len(genes & ga)}/{len(genes)}, "
                  f"{lb} {len(genes & gb)}/{len(genes)}"
                  f" | missing in {la}: {sorted(genes - ga)[:10]}"
                  f" | missing in {lb}: {sorted(genes - gb)[:10]}")

        print("\n== 3. Total UMIs of the same cells (phase-fraction denominator)")
        rng = np.random.default_rng(0)
        ids = (common if len(common) <= args.n_umi_cells
               else pd.Index(rng.choice(common, args.n_umi_cells, replace=False)))
        ids = ids.sort_values()
        ta = np.asarray(a[ids].X.sum(axis=1)).ravel()
        tb = np.asarray(b[ids].X.sum(axis=1)).ravel()
        ratio = ta / np.where(tb > 0, tb, np.nan)
        print(f"  {len(ids):,} shared cells | median total {la}: {np.median(ta):,.1f}, "
              f"{lb}: {np.median(tb):,.1f} | median ratio {la}/{lb}: {np.nanmedian(ratio):.3f}")
        print(f"  (ratio far from 1 → different genes or different count processing; "
              f"non-integer totals → X is normalized, not raw counts)")

        print("\n== 4. Correlation of per-cell scores on shared cells (Pearson)")
        for col in SCORE_COLS + CYCLING_COLS:
            if col in a.obs and col in b.obs:
                x = pd.to_numeric(a.obs.loc[common, col], errors="coerce")
                y = pd.to_numeric(b.obs.loc[common, col], errors="coerce")
                print(f"  {col:20s} r = {x.corr(y):.3f}")
        ca = [c for c in CYCLING_COLS if c in a.obs]
        cb = [c for c in CYCLING_COLS if c in b.obs]
        if ca and cb and ca[0] != cb[0]:
            x = pd.to_numeric(a.obs.loc[common, ca[0]], errors="coerce")
            y = pd.to_numeric(b.obs.loc[common, cb[0]], errors="coerce")
            print(f"  {la}.{ca[0]} vs {lb}.{cb[0]}  r = {x.corr(y):.3f}")
        missing = [c for c in SCORE_COLS if c not in a.obs or c not in b.obs]
        if missing:
            print(f"  (not in both files, skipped: {missing})")

        print(f"\n== 5. Phase cross-table, % of each {la} phase (rows) going to each {lb} phase")
        pa = a.obs.loc[common, args.phase_col].astype(str)
        pb = b.obs.loc[common, args.phase_col].astype(str)
        tab = pd.crosstab(pa, pb)
        pct = (pd.crosstab(pa, pb, normalize="index") * 100).round(1)
        print(pct.to_string())
        same = sum(tab.loc[r, c] for r in tab.index for c in tab.columns
                   if r.replace("-", "") == c.replace("-", ""))
        print(f"  overall agreement (PostM == Post-M): {100 * same / tab.values.sum():.1f}%")

        if args.out:
            out = Path(args.out)
            out.mkdir(parents=True, exist_ok=True)
            tab.to_csv(out / "phase_crosstab_counts.csv")
            pct.to_csv(out / "phase_crosstab_percent.csv")
            pd.DataFrame({"cell": ids, f"total_{la}": ta, f"total_{lb}": tb}).to_csv(
                out / "total_umis_shared_cells.csv", index=False)
            print(f"\nSaved CSVs to {out}")
    finally:
        a.file.close()
        b.file.close()


if __name__ == "__main__":
    main()
