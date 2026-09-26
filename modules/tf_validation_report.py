"""
tf_validation_report.py

Cross-reference tf_network's per-gene-list query outputs against the ATAC-seq
MOODS-scan CSVs to produce, per gene list, a version of the tf_targets CSV
annotated with validation flags along one or both axes.

Cell type (--atac-per-cell-type-dir):
  - `validated_<cell_type>` — one boolean column per cell type, True if that
    (TF, target) pair has an ATAC-supported motif hit in that cell type
    (source: `validated_by_motif_and_accessibility` in the per-CT CSV).
  - `validated_in_any_cell_type`, `n_cell_types_validated`,
    `validated_cell_types` (semicolon-separated list).

Cell-cycle phase (--atac-per-cell-cycle-dir):
  - `validated_phase_<phase>` — one boolean column per phase. Prefixed so the
    columns stay distinguishable from the cell-type ones in the same CSV.
  - `validated_in_any_cell_phase`, `n_cell_phases_validated`,
    `validated_cell_phases` (semicolon-separated list of every phase in which
    the pair was validated).

Phase columns are emitted in cell-cycle order (Non-cycling -> G1 -> S -> G2M
-> Post-M) rather than alphabetically, so the CSV reads in cycle order.

Either axis may be omitted; supplying both annotates against both.

Also writes a `validation_summary.csv` reporting what percentage of pairs were
validated on each axis, per input gene list, plus per-group counts.

Inputs:
  --tf-network-dir <path>
     A tf_network run directory. Its `data/` subfolder holds
     `tf_targets_<gene_list>.csv` files with columns (tf, target, importance,
     normalized_score).

  --atac-per-cell-type-dir <path>
     The `per_cell_type/` folder from a matching atac_seq_analysis run, e.g.
     `results/atac_analysis/<run>/3_motif_target_validation/per_cell_type/`
     with `motif_target_pair_scores_<cell_type>.csv` files.

  --atac-per-cell-cycle-dir <path>
     The `per_cell_cycle/` folder from the same run, same file naming.

Usage:
  python modules/tf_validation_report.py \
    --tf-network-dir results/tf_network/tf_network_cortex_v3_20260802 \
    --atac-per-cell-type-dir  results/atac_analysis/<run>/3_motif_target_validation/per_cell_type \
    --atac-per-cell-cycle-dir results/atac_analysis/<run>/3_motif_target_validation/per_cell_cycle \
    --subfolder-name cortex_v3_vs_atac_20260816
"""

import argparse
import re
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd


DEFAULT_OUTPUT_DIR = Path("results/tf_validation")

TF_KEY_CANDIDATES     = ("tf", "TF", "transcription_factor")
TARGET_KEY_CANDIDATES = ("target", "Target", "gene", "Gene")
VALIDATED_COL         = "validated_by_motif_and_accessibility"

# One entry per validation axis. Cell type keeps the column names it has
# always had; phase columns carry a prefix so both axes can sit in one CSV.
AXES = {
    "cell_type": dict(
        col_prefix="validated_",
        any_col="validated_in_any_cell_type",
        count_col="n_cell_types_validated",
        list_col="validated_cell_types",
        summary_any="n_validated_any_cell_type",
        summary_pct="pct_validated_any_cell_type",
        label="cell type",
        label_plural="cell types",
        short="CT",
    ),
    "cell_phase": dict(
        col_prefix="validated_phase_",
        any_col="validated_in_any_cell_phase",
        count_col="n_cell_phases_validated",
        list_col="validated_cell_phases",
        summary_any="n_validated_any_cell_phase",
        summary_pct="pct_validated_any_cell_phase",
        label="cell phase",
        label_plural="cell phases",
        short="phase",
    ),
}

# Cycle order for phase columns. Covers both annotation vocabularies: the
# promoter path emits Non-cycling/G1/S/G2M/Post-M, the replication path
# G1-G0/S/G2M. Anything unrecognized sorts alphabetically after these.
PHASE_ORDER_HINT = ("Non-cycling", "G0", "G1/G0", "G1", "S", "G2M", "Post-M")


def _phase_sort_key(slug: str):
    """Order phase slugs by position in the cell cycle, unknown ones last."""
    hint = [_sanitize(p) for p in PHASE_ORDER_HINT]
    s = str(slug)
    return (hint.index(s), "") if s in hint else (len(hint), s)


def _order_groups(axis: str, slugs) -> list[str]:
    slugs = list(slugs)
    return (sorted(slugs, key=_phase_sort_key) if axis == "cell_phase"
            else sorted(slugs))


def _sanitize(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_") or "unnamed"


def _pick_col(df: pd.DataFrame, candidates) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(
        f"None of {candidates} found in columns: {list(df.columns)}"
    )


# ---------------------------------------------------------------------------
# Load per-cell-type ATAC validation
# ---------------------------------------------------------------------------

def load_group_validation(
    group_dir: Path, label: str = "cell type",
) -> Dict[str, Dict[Tuple[str, str], bool]]:
    """Return {group_slug: {(TF_upper, target_upper): bool}} for each
    `motif_target_pair_scores_<group>.csv` under `group_dir`."""
    if not group_dir.is_dir():
        raise FileNotFoundError(f"per-{label} dir not found: {group_dir}")

    result: Dict[str, Dict[Tuple[str, str], bool]] = {}
    csvs = sorted(group_dir.glob("motif_target_pair_scores_*.csv"))
    if not csvs:
        raise FileNotFoundError(
            f"No motif_target_pair_scores_*.csv files under {group_dir}"
        )

    for path in csvs:
        # `motif_target_pair_scores_radial_glial_cell.csv` -> `radial_glial_cell`
        ct_slug = path.stem[len("motif_target_pair_scores_"):]

        df = pd.read_csv(path)
        try:
            tf_col     = _pick_col(df, TF_KEY_CANDIDATES)
            target_col = _pick_col(df, TARGET_KEY_CANDIDATES)
        except KeyError as e:
            print(f"  [skip] {path.name}: {e}")
            continue
        if VALIDATED_COL not in df.columns:
            print(f"  [skip] {path.name}: missing '{VALIDATED_COL}' column")
            continue

        keys   = zip(df[tf_col].astype(str).str.upper(),
                     df[target_col].astype(str).str.upper())
        values = df[VALIDATED_COL].astype(bool)
        result[ct_slug] = dict(zip(keys, values))
        print(f"  {ct_slug}: {len(df):,} pairs  "
              f"({int(values.sum()):,} validated)")
    return result


# ---------------------------------------------------------------------------
# Annotate one tf_network CSV
# ---------------------------------------------------------------------------

def annotate_tf_network_csv(
    tfnet_csv: Path,
    axis_validation: Dict[str, Dict[str, Dict[Tuple[str, str], bool]]],
    out_csv: Path,
) -> dict:
    """Read a `tf_targets_*.csv` and add, for every axis in `axis_validation`,
    per-group boolean columns plus any/count/list summary columns. Save to
    `out_csv` and return stats for the aggregate report.

    `axis_validation` maps axis name -> {group_slug: {(TF, target): bool}}.
    """
    df = pd.read_csv(tfnet_csv)

    try:
        tf_col     = _pick_col(df, TF_KEY_CANDIDATES)
        target_col = _pick_col(df, TARGET_KEY_CANDIDATES)
    except KeyError as e:
        raise KeyError(f"{tfnet_csv.name}: {e}") from e

    # Build lookup keys once.
    keys = list(zip(df[tf_col].astype(str).str.upper(),
                    df[target_col].astype(str).str.upper()))

    stats = {"n_pairs": int(len(df)), "axes": {}}

    for axis, validation in axis_validation.items():
        spec  = AXES[axis]
        slugs = _order_groups(axis, validation.keys())
        cols  = []
        for g in slugs:
            d = validation[g]
            colname = f"{spec['col_prefix']}{g}"
            cols.append(colname)
            df[colname] = [d.get(k, False) for k in keys]

        df[spec["any_col"]]   = df[cols].any(axis=1) if cols else False
        df[spec["count_col"]] = df[cols].sum(axis=1).astype(int) if cols else 0

        # Semicolon-separated list of the groups where the pair passed.
        # Empty string when the pair wasn't validated anywhere on this axis.
        if cols:
            matrix = df[cols].to_numpy(dtype=bool)
            df[spec["list_col"]] = [
                ";".join(slugs[i] for i in np.where(row)[0]) for row in matrix
            ]
        else:
            df[spec["list_col"]] = ""

        n_any = int(df[spec["any_col"]].sum())
        stats["axes"][axis] = {
            "slugs": slugs,
            "n_validated_any": n_any,
            "pct_validated_any": (100.0 * n_any / len(df)) if len(df) else 0.0,
            "per_group_validated": {
                g: int(df[f"{spec['col_prefix']}{g}"].sum()) for g in slugs
            },
        }

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    return stats


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def _find_tf_network_csvs(tf_network_dir: Path) -> list[Path]:
    """Find tf_targets_*.csv files. Checks the dir itself and its `data/` sub."""
    candidates = list((tf_network_dir / "data").glob("tf_targets_*.csv")) \
                 if (tf_network_dir / "data").is_dir() else []
    if not candidates:
        candidates = list(tf_network_dir.glob("tf_targets_*.csv"))
    return sorted(candidates)


def run(
    tf_network_dir: str,
    atac_per_cell_type_dir: str | None = None,
    atac_per_cell_cycle_dir: str | None = None,
    output_dir: str = str(DEFAULT_OUTPUT_DIR),
    subfolder_name: str | None = None,
) -> Path:
    tf_network_dir = Path(tf_network_dir).resolve()
    output_dir     = Path(output_dir)

    axis_dirs = {}
    if atac_per_cell_type_dir:
        axis_dirs["cell_type"] = Path(atac_per_cell_type_dir).resolve()
    if atac_per_cell_cycle_dir:
        axis_dirs["cell_phase"] = Path(atac_per_cell_cycle_dir).resolve()
    if not axis_dirs:
        raise ValueError(
            "Give at least one of --atac-per-cell-type-dir / --atac-per-cell-cycle-dir."
        )

    if subfolder_name:
        output_dir = output_dir / _sanitize(subfolder_name)
    else:
        # Default: encode both sources so different comparisons don't collide.
        # a group dir looks like .../<atac_run>/3_motif_target_validation/<group>
        any_dir = next(iter(axis_dirs.values()))
        atac_run_name = any_dir.parent.parent.name
        output_dir = output_dir / f"{tf_network_dir.name}__vs__{atac_run_name}"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\ntf_network dir:   {tf_network_dir}")
    for axis, d in axis_dirs.items():
        print(f"atac per-{AXES[axis]['label']} dir: {d}")
    print(f"output dir:       {output_dir}")

    axis_validation = {}
    for axis, d in axis_dirs.items():
        spec = AXES[axis]
        print(f"\nLoading per-{spec['label']} ATAC validation:")
        axis_validation[axis] = load_group_validation(d, spec["label"])
        print(f"-> {len(axis_validation[axis])} {spec['label_plural']} loaded")

    tfnet_csvs = _find_tf_network_csvs(tf_network_dir)
    if not tfnet_csvs:
        raise FileNotFoundError(
            f"No tf_targets_*.csv files found under {tf_network_dir} or its data/ subdir."
        )
    print(f"\nFound {len(tfnet_csvs)} tf_network CSV(s):")
    for p_ in tfnet_csvs:
        print(f"  {p_.name}")

    all_summaries = []
    for csv_path in tfnet_csvs:
        gene_list = csv_path.stem.removeprefix("tf_targets_")
        out_csv = output_dir / csv_path.name
        try:
            stats = annotate_tf_network_csv(csv_path, axis_validation, out_csv)
        except KeyError as e:
            print(f"  [skip] {csv_path.name}: {e}")
            continue
        stats["gene_list"]  = gene_list
        stats["source_csv"] = str(csv_path)
        stats["output_csv"] = str(out_csv)
        all_summaries.append(stats)

        parts = [
            f"{stats['axes'][a]['n_validated_any']:,}/{stats['n_pairs']:,} "
            f"({stats['axes'][a]['pct_validated_any']:.1f}%) in >=1 {AXES[a]['label']}"
            for a in axis_validation
        ]
        print(f"  {csv_path.name}: " + "; ".join(parts))

    _write_summary(output_dir, all_summaries, axis_validation)
    print(f"\nDone -> {output_dir}")
    return output_dir


def _write_summary(output_dir: Path, summaries: list[dict], axis_validation: dict) -> None:
    if not summaries:
        print("[warn] no per-list summaries — skipping validation_summary.csv")
        return

    axes_present = list(axis_validation.keys())
    slugs_by_axis = {
        a: _order_groups(a, axis_validation[a].keys()) for a in axes_present
    }

    rows = []
    for s in summaries:
        n = s["n_pairs"]
        row = {
            "gene_list": s["gene_list"],
            "source_csv": s["source_csv"],
            "output_csv": s["output_csv"],
            "n_pairs": n,
        }
        for axis in axes_present:
            spec = AXES[axis]
            ax = s["axes"][axis]
            row[spec["summary_any"]] = ax["n_validated_any"]
            row[spec["summary_pct"]] = ax["pct_validated_any"]
            for g in slugs_by_axis[axis]:
                n_valid = ax["per_group_validated"].get(g, 0)
                row[f"n_validated_{spec['col_prefix'][len('validated_'):]}{g}"] = n_valid
                row[f"pct_validated_{spec['col_prefix'][len('validated_'):]}{g}"] = (
                    (100.0 * n_valid / n) if n else 0.0)
        rows.append(row)

    sort_col = AXES[axes_present[0]]["summary_pct"]
    df = pd.DataFrame(rows).sort_values(sort_col, ascending=False)
    summary_path = output_dir / "validation_summary.csv"
    df.to_csv(summary_path, index=False)
    print(f"\nSaved: {summary_path}")

    # Human-readable text report.
    txt = output_dir / "validation_report.txt"
    with open(txt, "w") as f:
        f.write("TF-network vs ATAC-seq validation report\n")
        f.write("=" * 60 + "\n\n")
        header = f"{'Gene list':40s}  {'Pairs':>7s}"
        for axis in axes_present:
            header += f"  {'Val >=1 ' + AXES[axis]['short']:>12s}  {'%':>6s}"
        f.write(header + "\n")
        f.write("-" * len(header) + "\n")
        for r in rows:
            line = f"{r['gene_list']:40s}  {r['n_pairs']:>7,}"
            for axis in axes_present:
                spec = AXES[axis]
                line += f"  {r[spec['summary_any']]:>12,}  {r[spec['summary_pct']]:>5.1f}%"
            f.write(line + "\n")

        for axis in axes_present:
            spec = AXES[axis]
            key  = spec["col_prefix"][len("validated_"):]
            f.write(f"\nPer-{spec['label']} validation counts "
                    f"(rows = gene lists, cols = {spec['label_plural']}):\n\n")
            for r in rows:
                f.write(f"[{r['gene_list']}]  total pairs = {r['n_pairs']:,}\n")
                for g in slugs_by_axis[axis]:
                    n_valid = r[f"n_validated_{key}{g}"]
                    pct     = r[f"pct_validated_{key}{g}"]
                    f.write(f"    {g:45s}  {n_valid:>6,}  ({pct:5.1f}%)\n")
                f.write("\n")
    print(f"Saved: {txt}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_arg_parser():
    p = argparse.ArgumentParser(
        description=(
            "Cross-check tf_network's per-gene-list outputs against "
            "atac_seq_analysis's per-cell-type MOODS-scan CSVs, and produce "
            "annotated CSVs + a validation summary."
        )
    )
    p.add_argument("--tf-network-dir", required=True,
                   help="tf_network run dir (holds data/tf_targets_*.csv).")
    p.add_argument("--atac-per-cell-type-dir", default=None,
                   help="per_cell_type/ folder under an atac_seq_analysis run's "
                        "3_motif_target_validation/.")
    p.add_argument("--atac-per-cell-cycle-dir", default=None,
                   help="per_cell_cycle/ folder under the same run. Adds "
                        "validated_phase_<phase> columns plus "
                        "validated_cell_phases listing every phase in which "
                        "each pair was validated.")
    p.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR),
                   help="Base output folder (default: results/tf_validation).")
    p.add_argument("--subfolder-name", default=None,
                   help="Subfolder under --output-dir. Defaults to "
                        "'<tf_net_dir_name>__vs__<atac_run_name>'.")
    return p


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    run(
        tf_network_dir=args.tf_network_dir,
        atac_per_cell_type_dir=args.atac_per_cell_type_dir,
        atac_per_cell_cycle_dir=args.atac_per_cell_cycle_dir,
        output_dir=args.output_dir,
        subfolder_name=args.subfolder_name,
    )
