"""
build_batch_summary.py

Collect GSEA_final_summary.csv files from a set of enrichment run folders and
produce a single batch_summary CSV with the same schema as the one produced by
`enrichment_cal_lists_loop.py`.

Columns:
  run_name, enrichment_config_path, ges_config_path, gene_list_file_name,
  num_genes_in_gene_list, column_condition_title, column_condition_value,
  condition_compared_to, NES, FDR_qval_BH, is_significant, lead_genes,
  num_of_lead_genes

Two invocation styles:

  # 1) Explicit list of run dirs (recommended — used by the bash launcher):
  python modules/build_batch_summary.py \
    --run-dirs run_dir_1 run_dir_2 ... \
    --output-csv results/enrichment_results/batch_summary_YYYYMMDD_myrun.csv

  # 2) Glob a base folder by pattern:
  python modules/build_batch_summary.py \
    --base-folder results/enrichment_results \
    --run-glob '*BDgene*_threshold_*_20260723' \
    --output-csv results/enrichment_results/batch_summary_20260723.csv
"""

import argparse
import sys
from pathlib import Path

# Make sibling modules importable regardless of CWD.
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

from enrichment_cal_lists_loop import (
    _conditions_compared_to,
    _find_nes_column,
    _find_qval_column,
    _locate_ges_config_in_folder,
    _read_gene_list_size,
    _safe_load_yaml,
)


def _resolve_relative(path_str: str, root: str) -> Path:
    p = Path(str(path_str))
    if p.is_absolute():
        return p
    return (Path(root) / p).resolve()


def build_rows_from_run_dir(run_dir: Path) -> list[dict]:
    """Read one run dir's config + GSEA summary and return per-row records."""
    cfg_path = run_dir / "metadata" / "config_used.yaml"
    if not cfg_path.exists():
        yamls = sorted((run_dir / "metadata").glob("*.yaml"))
        if not yamls:
            print(f"  [skip] no config yaml under {run_dir}/metadata")
            return []
        cfg_path = yamls[0]

    cfg      = _safe_load_yaml(cfg_path)
    root     = cfg.get("ndd_gene_modules_folder_root", ".")
    run_name = cfg.get("run_name", run_dir.name)

    gene_list_path = _resolve_relative(cfg.get("gene_list_path", ""), root)
    gene_list_n    = (
        _read_gene_list_size(gene_list_path)
        if gene_list_path.exists() else np.nan
    )

    ges_folder   = _resolve_relative(cfg.get("ges_results_folder", ""), root)
    ges_cfg_path = _locate_ges_config_in_folder(ges_folder)

    summary_path = run_dir / "data" / "enrichment_results" / "GSEA" / "GSEA_final_summary.csv"
    if not summary_path.exists():
        print(f"  [skip] no GSEA_final_summary.csv in {run_dir}")
        return []

    df      = pd.read_csv(summary_path)
    nes_col = _find_nes_column(df)
    q_col   = _find_qval_column(df)
    if "column" not in df.columns or "condition" not in df.columns or nes_col is None:
        print(f"  [skip] required columns missing in {summary_path}")
        return []

    rows = []
    for _, r in df.iterrows():
        col_name    = str(r["column"])
        cond_val    = str(r["condition"])
        nes_val     = r[nes_col]
        q_val       = r[q_col] if q_col is not None and q_col in df.columns else np.nan
        lead_genes  = r.get("Lead_genes", "")
        num_lead    = (
            len(str(lead_genes).split(";"))
            if pd.notna(lead_genes) and str(lead_genes).strip() else 0
        )

        is_sig = bool(
            pd.notna(nes_val) and (nes_val > 1.5 or nes_val < -1.5)
            and pd.notna(q_val) and q_val <= 0.05
        )

        rows.append({
            "run_name": run_name,
            "enrichment_config_path": str(cfg_path),
            "ges_config_path": ges_cfg_path,
            "gene_list_file_name": gene_list_path.name,
            "num_genes_in_gene_list": gene_list_n,
            "column_condition_title": col_name,
            "column_condition_value": cond_val,
            "condition_compared_to": ";".join(_conditions_compared_to(cfg, col_name, cond_val)),
            "NES": nes_val,
            "FDR_qval_BH": q_val,
            "is_significant": is_sig,
            "lead_genes": lead_genes,
            "num_of_lead_genes": num_lead,
        })
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__.strip().splitlines()[1])
    p.add_argument("--run-dirs", nargs="+", default=None,
                   help="Explicit list of run directory paths.")
    p.add_argument("--base-folder", default=None,
                   help="Base enrichment output folder for --run-glob mode.")
    p.add_argument("--run-glob", default=None,
                   help="Glob pattern under --base-folder, e.g. '*_threshold_*_20260723'.")
    p.add_argument("--output-csv", required=True,
                   help="Where to write the batch summary CSV.")
    args = p.parse_args()

    run_dirs: list[Path] = []
    if args.run_dirs:
        run_dirs.extend(Path(x).resolve() for x in args.run_dirs)
    if args.base_folder and args.run_glob:
        base = Path(args.base_folder).resolve()
        run_dirs.extend(sorted(base.glob(args.run_glob)))

    if not run_dirs:
        print("No run dirs provided (use --run-dirs or --base-folder + --run-glob).")
        sys.exit(1)

    seen: set[Path] = set()
    unique_dirs = []
    for d in run_dirs:
        if d not in seen:
            seen.add(d)
            unique_dirs.append(d)

    all_rows = []
    for run_dir in unique_dirs:
        if not run_dir.is_dir():
            print(f"  [skip] not a directory: {run_dir}")
            continue
        print(f"Processing: {run_dir}")
        all_rows.extend(build_rows_from_run_dir(run_dir))

    if not all_rows:
        print("\nNo summary rows collected; nothing written.")
        sys.exit(1)

    out_path = Path(args.output_csv).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(all_rows).to_csv(out_path, index=False)
    print(f"\nWrote batch summary CSV: {out_path}  ({len(all_rows)} rows)")


if __name__ == "__main__":
    main()
