"""
tf_validation_report.py

Cross-reference tf_network's per-gene-list query outputs against the ATAC-seq
per-cell-type MOODS-scan CSVs to produce, per gene list, a version of the
tf_targets CSV annotated with:

  - `validated_<cell_type>` — one boolean column per cell type, True if that
    (TF, target) pair has an ATAC-supported motif hit in that cell type
    (source: `validated_by_motif_and_accessibility` in the per-CT CSV).
  - `validated_in_any_cell_type` — True if any of the above are True.

Also writes a `validation_summary.csv` reporting what percentage of pairs were
validated in any cell type, per input gene list, plus per-cell-type counts.

Inputs:
  --tf-network-dir <path>
     A tf_network run directory. Its `data/` subfolder holds
     `tf_targets_<gene_list>.csv` files with columns (tf, target, importance,
     normalized_score).

  --atac-per-cell-type-dir <path>
     The `per_cell_type/` folder from a matching atac_seq_analysis run, e.g.
     `results/atac_analysis/<run>/3_motif_target_validation/per_cell_type/`
     with `motif_target_pair_scores_<cell_type>.csv` files.

Usage:
  python modules/tf_validation_report.py \
    --tf-network-dir results/tf_network/tf_network_cortex_v3_20260802 \
    --atac-per-cell-type-dir results/atac_analysis/atac_first_trimester_brain_20260816/3_motif_target_validation/per_cell_type \
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

def load_ct_validation(
    per_ct_dir: Path,
) -> Dict[str, Dict[Tuple[str, str], bool]]:
    """Return {cell_type_slug: {(TF_upper, target_upper): bool}} for each
    `motif_target_pair_scores_<ct>.csv` under `per_ct_dir`."""
    if not per_ct_dir.is_dir():
        raise FileNotFoundError(f"per-cell-type dir not found: {per_ct_dir}")

    result: Dict[str, Dict[Tuple[str, str], bool]] = {}
    csvs = sorted(per_ct_dir.glob("motif_target_pair_scores_*.csv"))
    if not csvs:
        raise FileNotFoundError(
            f"No motif_target_pair_scores_*.csv files under {per_ct_dir}"
        )

    for path in csvs:
        # `motif_target_pair_scores_radial_glial_cell.csv` → `radial_glial_cell`
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
    ct_validation: Dict[str, Dict[Tuple[str, str], bool]],
    out_csv: Path,
) -> dict:
    """Read a `tf_targets_*.csv`, add per-CT `validated_<ct>` columns and an
    `validated_in_any_cell_type` column, save to `out_csv`. Returns summary
    stats used by the aggregate report."""
    df = pd.read_csv(tfnet_csv)

    try:
        tf_col     = _pick_col(df, TF_KEY_CANDIDATES)
        target_col = _pick_col(df, TARGET_KEY_CANDIDATES)
    except KeyError as e:
        raise KeyError(f"{tfnet_csv.name}: {e}") from e

    # Build lookup keys once.
    keys = list(zip(df[tf_col].astype(str).str.upper(),
                    df[target_col].astype(str).str.upper()))

    ct_slugs = sorted(ct_validation.keys())
    ct_cols  = []
    for ct in ct_slugs:
        d      = ct_validation[ct]
        colname = f"validated_{ct}"
        ct_cols.append(colname)
        df[colname] = [d.get(k, False) for k in keys]

    df["validated_in_any_cell_type"] = df[ct_cols].any(axis=1) if ct_cols else False

    # Also expose how many CTs each pair was validated in — useful for sorting.
    df["n_cell_types_validated"] = df[ct_cols].sum(axis=1).astype(int) if ct_cols else 0

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)

    n = int(len(df))
    n_any = int(df["validated_in_any_cell_type"].sum())
    per_ct = {ct: int(df[f"validated_{ct}"].sum()) for ct in ct_slugs}
    return {
        "n_pairs": n,
        "n_validated_any": n_any,
        "pct_validated_any": (100.0 * n_any / n) if n else 0.0,
        "per_ct_validated": per_ct,
    }


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
    atac_per_cell_type_dir: str,
    output_dir: str = str(DEFAULT_OUTPUT_DIR),
    subfolder_name: str | None = None,
) -> Path:
    tf_network_dir         = Path(tf_network_dir).resolve()
    atac_per_cell_type_dir = Path(atac_per_cell_type_dir).resolve()
    output_dir             = Path(output_dir)

    if subfolder_name:
        output_dir = output_dir / _sanitize(subfolder_name)
    else:
        # Default: encode both sources so different comparisons don't collide.
        # atac_per_ct parent path: .../<atac_run>/3_motif_target_validation/per_cell_type
        atac_run_name = atac_per_cell_type_dir.parent.parent.name
        output_dir = output_dir / f"{tf_network_dir.name}__vs__{atac_run_name}"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\ntf_network dir:   {tf_network_dir}")
    print(f"atac per-CT dir:  {atac_per_cell_type_dir}")
    print(f"output dir:       {output_dir}")

    print("\nLoading per-cell-type ATAC validation:")
    ct_validation = load_ct_validation(atac_per_cell_type_dir)
    print(f"→ {len(ct_validation)} cell types loaded")

    tfnet_csvs = _find_tf_network_csvs(tf_network_dir)
    if not tfnet_csvs:
        raise FileNotFoundError(
            f"No tf_targets_*.csv files found under {tf_network_dir} or its data/ subdir."
        )
    print(f"\nFound {len(tfnet_csvs)} tf_network CSV(s):")
    for p in tfnet_csvs:
        print(f"  {p.name}")

    all_summaries = []
    for csv_path in tfnet_csvs:
        gene_list = csv_path.stem.removeprefix("tf_targets_")
        out_csv = output_dir / csv_path.name
        try:
            stats = annotate_tf_network_csv(csv_path, ct_validation, out_csv)
        except KeyError as e:
            print(f"  [skip] {csv_path.name}: {e}")
            continue
        stats["gene_list"]  = gene_list
        stats["source_csv"] = str(csv_path)
        stats["output_csv"] = str(out_csv)
        all_summaries.append(stats)
        print(f"  {csv_path.name}: "
              f"{stats['n_validated_any']:,}/{stats['n_pairs']:,} "
              f"({stats['pct_validated_any']:.1f}%) validated in ≥1 cell type")

    _write_summary(output_dir, all_summaries, sorted(ct_validation.keys()))
    print(f"\nDone → {output_dir}")
    return output_dir


def _write_summary(output_dir: Path, summaries: list[dict], ct_slugs: list[str]) -> None:
    if not summaries:
        print("[warn] no per-list summaries — skipping validation_summary.csv")
        return

    rows = []
    for s in summaries:
        n = s["n_pairs"]
        row = {
            "gene_list": s["gene_list"],
            "source_csv": s["source_csv"],
            "output_csv": s["output_csv"],
            "n_pairs": n,
            "n_validated_any_cell_type": s["n_validated_any"],
            "pct_validated_any_cell_type": s["pct_validated_any"],
        }
        for ct in ct_slugs:
            n_valid = s["per_ct_validated"].get(ct, 0)
            row[f"n_validated_{ct}"]  = n_valid
            row[f"pct_validated_{ct}"] = (100.0 * n_valid / n) if n else 0.0
        rows.append(row)

    df = pd.DataFrame(rows).sort_values("pct_validated_any_cell_type", ascending=False)
    summary_path = output_dir / "validation_summary.csv"
    df.to_csv(summary_path, index=False)
    print(f"\nSaved: {summary_path}")

    # Human-readable text report.
    txt = output_dir / "validation_report.txt"
    with open(txt, "w") as f:
        f.write("TF-network vs ATAC-seq validation report\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"{'Gene list':40s}  {'Pairs':>7s}  {'Val≥1CT':>7s}  {'%≥1CT':>6s}\n")
        f.write("-" * 66 + "\n")
        for r in rows:
            f.write(
                f"{r['gene_list']:40s}  "
                f"{r['n_pairs']:>7,}  "
                f"{r['n_validated_any_cell_type']:>7,}  "
                f"{r['pct_validated_any_cell_type']:>5.1f}%\n"
            )
        f.write("\nPer-cell-type validation counts (rows = gene lists, cols = cell types):\n\n")
        for r in rows:
            f.write(f"[{r['gene_list']}]  total pairs = {r['n_pairs']:,}\n")
            for ct in ct_slugs:
                n_valid = r[f"n_validated_{ct}"]
                pct     = r[f"pct_validated_{ct}"]
                f.write(f"    {ct:45s}  {n_valid:>6,}  ({pct:5.1f}%)\n")
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
    p.add_argument("--atac-per-cell-type-dir", required=True,
                   help="per_cell_type/ folder under an atac_seq_analysis run's "
                        "3_motif_target_validation/.")
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
        output_dir=args.output_dir,
        subfolder_name=args.subfolder_name,
    )
