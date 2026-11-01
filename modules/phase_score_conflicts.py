"""
phase_score_conflicts.py

List the cells whose cycling score says Non-cycling (score <= cutoff, compared
at the score's stored precision, as in cell_cycle_annotation.py) but whose
CellCyclePhase is something else. One combined CSV over all given files, with
a `dataset` column, the cell's phase and score, and context columns when
present (CellClass, Region, Age, Chemistry, phase scores).

Usage:
  python modules/phase_score_conflicts.py \
    --dataset cortex=data/Cortex_EMX1_louvain3_passedQC_PostM_rev1.h5ad \
    --dataset human_dev=data/human_dev_without_week_5_cc_annotated_fraction.h5ad \
    --out results/additional_analyses/cell_cycle_comparison/phase_score_conflicts.csv
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from cell_cycle_annotation import CYCLING_THRESHOLD, _above_cutoff

SCORE_COLS   = ["CellCycleFraction", "cell_cycle_score"]
CONTEXT_COLS = ["CellClass", "Region", "Age", "Chemistry",
                "cc_score_g1", "cc_score_s", "cc_score_g2m"]


def conflicts(path: str, label: str, phase_col: str, cutoff: float) -> pd.DataFrame:
    a = ad.read_h5ad(path, backed="r")
    try:
        obs = a.obs
        score_col = next((c for c in SCORE_COLS if c in obs), None)
        if score_col is None or phase_col not in obs:
            print(f"{label}: [skip] needs one of {SCORE_COLS} and '{phase_col}'")
            return pd.DataFrame()
        score = pd.to_numeric(obs[score_col], errors="coerce").to_numpy()
        phase = obs[phase_col].astype(object)
        non_cycling_by_score = ~np.isnan(score) & ~_above_cutoff(score, cutoff)
        mask = non_cycling_by_score & phase.notna().to_numpy() & (phase != "Non-cycling").to_numpy()

        out = pd.DataFrame({
            "dataset": label,
            "cell": obs.index[mask],
            "phase": phase[mask].to_numpy(),
            "score_column": score_col,
            "score": score[mask],
        })
        for c in CONTEXT_COLS:
            if c in obs:
                out[c] = obs[c].to_numpy()[mask]
        print(f"{label}: {int(mask.sum()):,} cells with {score_col} <= {cutoff:g} "
              f"(of {int(non_cycling_by_score.sum()):,}) but {phase_col} != Non-cycling")
        if mask.any():
            print("  by phase: " + ", ".join(f"{k} {v:,}" for k, v in out["phase"].value_counts().items()))
            print(f"  score range: {out['score'].min():.6g} – {out['score'].max():.6g}")
        return out
    finally:
        a.file.close()


def main():
    p = argparse.ArgumentParser(description="Cells Non-cycling by score but with another phase.")
    p.add_argument("--dataset", action="append", required=True, metavar="LABEL=PATH",
                   help="Dataset to check; repeat for several.")
    p.add_argument("--phase-col", default="CellCyclePhase")
    p.add_argument("--cutoff", type=float, default=CYCLING_THRESHOLD)
    p.add_argument("--out", required=True, help="Output CSV path.")
    args = p.parse_args()

    frames = []
    for spec in args.dataset:
        label, _, path = spec.partition("=")
        if not path:
            p.error(f"--dataset expects LABEL=PATH, got: {spec}")
        frames.append(conflicts(path, label, args.phase_col, args.cutoff))

    result = pd.concat([f for f in frames if not f.empty], ignore_index=True) \
        if any(not f.empty for f in frames) else pd.DataFrame()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out, index=False)
    print(f"\nSaved {len(result):,} rows → {out}")


if __name__ == "__main__":
    main()
