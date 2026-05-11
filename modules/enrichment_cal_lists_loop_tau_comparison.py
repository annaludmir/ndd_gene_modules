"""
enrichment_cal_lists_loop_tau_comparison.py

Run GSEA enrichment in 4 variants (v2/v3 chemistry × with/without tau filtering)
across all three scopes (cortex, cell_phase, all_layers).

This allows direct comparison of:
  1. v3 chemistry, no tau filtering
  2. v2 chemistry, no tau filtering
  3. v3 chemistry, tau filtered
  4. v2 chemistry, tau filtered

Accepts either a single gene-list CSV or a folder of CSVs. One batch-summary CSV
is produced per gene list, named:
  {gene_list_stem}_batch_summary_tau_vs_v2_v3_{YYYYMMDD}.csv

Usage:
  python modules/enrichment_cal_lists_loop_tau_comparison.py <gene_list.csv|folder>
         [--ndd-root PATH] [--tau-percentile N]

  --ndd-root        defaults to /miridan-data/annaludmir/ndd_gene_modules
  --tau-percentile  overrides the tau_percentile field in all tau-filtered configs
                    (default: read from each config file, typically 90)
"""

import argparse
import datetime
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


# ---------------------------------------------------------------------------
# Config registry — (scope, chemistry, tau_filtered) -> path relative to ndd_root
# ---------------------------------------------------------------------------

_ENRICHMENT_CONFIGS = {
    ("cortex",     "v3", False): "config_files/enrichment_cortex_config.yaml",
    ("cortex",     "v2", False): "config_files/enrichment_cortex_config_v2.yaml",
    ("cortex",     "v3", True):  "config_files/enrichment_cortex_config_tau_filtered.yaml",
    ("cortex",     "v2", True):  "config_files/enrichment_cortex_config_tau_filtered_v2.yaml",
    ("cell_phase", "v3", False): "config_files/enrichment_cortex_cell_phase_config.yaml",
    ("cell_phase", "v2", False): "config_files/enrichment_cortex_cell_phase_config_v2.yaml",
    ("cell_phase", "v3", True):  "config_files/enrichment_cortex_cell_phase_config_tau_filtered.yaml",
    ("cell_phase", "v2", True):  "config_files/enrichment_cortex_cell_phase_config_tau_filtered_v2.yaml",
    ("all_layers", "v3", False): "config_files/enrichment_all_layers_config.yaml",
    ("all_layers", "v2", False): "config_files/enrichment_all_layers_config_v2.yaml",
    ("all_layers", "v3", True):  "config_files/enrichment_all_layers_config_tau_filtered.yaml",
    ("all_layers", "v2", True):  "config_files/enrichment_all_layers_config_tau_filtered_v2.yaml",
}

# Summary CSV paths relative to run_dir, for each pipeline type
_SUMMARY_REL = {
    False: "data/enrichment_results/GSEA/GSEA_final_summary.csv",
    True:  "data/enrichment_results/GSEA_tau_filtered/GSEA_tau_filtered_summary.csv",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_gene_list_size(gene_list_csv: Path) -> int:
    df = pd.read_csv(gene_list_csv)
    col = "gene" if "gene" in df.columns else df.columns[0]
    genes = (
        df[col].astype(str).str.strip()
        .replace({"nan": np.nan, "None": np.nan})
        .dropna().tolist()
    )
    return len(set(genes))


def _find_nes_column(df: pd.DataFrame) -> str | None:
    for c in ["NES", "nes", "nes_score"]:
        if c in df.columns:
            return c
    return None


def _find_qval_column(df: pd.DataFrame) -> str | None:
    for c in ["FDR q-val (BH corrected)", "FDR q-val", "FDR_qval_BH"]:
        if c in df.columns:
            return c
    return None


def _conditions_compared_to(cfg: dict, column: str, condition_value: str) -> list[str]:
    vals = [str(v) for v in cfg.get("column_conditions_for_gsea", {}).get(column, [])]
    return [v for v in vals if v != str(condition_value)]


def _locate_ges_config(ges_results_folder: Path) -> str | None:
    for d in [ges_results_folder / "metadata", ges_results_folder]:
        if d.is_dir():
            ymls = sorted(list(d.glob("*.yaml")) + list(d.glob("*.yml")))
            if ymls:
                return str(ymls[0])
    return None


def _load_yaml(path: Path) -> dict:
    with path.open("r") as f:
        return yaml.safe_load(f)


def _dump_yaml(cfg: dict, path: Path) -> None:
    with path.open("w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)


# ---------------------------------------------------------------------------
# Summary-row collector
# ---------------------------------------------------------------------------

def _collect_summary_rows(
    run_dir: Path,
    tau_filtered: bool,
    run_name: str,
    scope: str,
    chemistry: str,
    gene_list_path: Path,
    gene_list_n: int,
    base_cfg: dict,
    ges_cfg_path: str | None,
) -> list[dict]:
    summary_path = run_dir / _SUMMARY_REL[tau_filtered]

    if not summary_path.exists():
        print(f"  ⚠️ Summary not found: {summary_path}")
        return []

    df_sum = pd.read_csv(summary_path)
    nes_col = _find_nes_column(df_sum)
    q_col   = _find_qval_column(df_sum)

    if "column" not in df_sum.columns or "condition" not in df_sum.columns or nes_col is None:
        print(f"  ⚠️ Missing required columns in {summary_path.name}")
        return []

    rows = []
    for _, r in df_sum.iterrows():
        col_name   = str(r["column"])
        cond_val   = str(r["condition"])
        nes_val    = r[nes_col]
        lead_genes = r.get("Lead_genes", "")
        q_val      = r[q_col] if q_col and q_col in df_sum.columns else np.nan

        try:
            is_sig = bool(
                (float(nes_val) > 1.5 or float(nes_val) < -1.5)
                and float(q_val) <= 0.05
            )
        except (TypeError, ValueError):
            is_sig = False

        n_lead = (
            len(str(lead_genes).split(";"))
            if pd.notna(lead_genes) and str(lead_genes).strip()
            else 0
        )

        rows.append({
            "run_name":               run_name,
            "scope":                  scope,
            "chemistry":              chemistry,
            "tau_filtered":           tau_filtered,
            "enrichment_config_path": str(run_dir / "metadata" / "config_used.yaml"),
            "ges_config_path":        ges_cfg_path,
            "gene_list_file_name":    gene_list_path.name,
            "num_genes_in_gene_list": gene_list_n,
            "column_condition_title": col_name,
            "column_condition_value": cond_val,
            "condition_compared_to":  ";".join(
                _conditions_compared_to(base_cfg, col_name, cond_val)
            ),
            "NES":           nes_val,
            "FDR_qval_BH":   q_val,
            "is_significant": is_sig,
            "lead_genes":     lead_genes,
            "num_of_lead_genes": n_lead,
        })

    return rows


# ---------------------------------------------------------------------------
# Per-gene-list runner
# ---------------------------------------------------------------------------

def _run_one(gene_list_path: Path, ndd_root: Path, today: str, tau_percentile: int | None = None) -> None:
    gene_list_stem = gene_list_path.stem
    gene_list_n    = _read_gene_list_size(gene_list_path)

    print(f"\n{'='*64}")
    print(f"Gene list: {gene_list_path.name}  ({gene_list_n} genes)")
    print(f"Running 12 enrichment variants.")
    print(f"{'='*64}")

    all_rows: list[dict] = []

    for (scope, chemistry, tau_filtered), cfg_rel in _ENRICHMENT_CONFIGS.items():
        cfg_path = ndd_root / cfg_rel
        if not cfg_path.exists():
            print(f"  ⚠️ Config not found: {cfg_path} — skipping.")
            continue

        base_cfg = _load_yaml(cfg_path)
        min_thr  = base_cfg.get("gsea", {}).get("min_ges_score_threshold", 1)
        tau_pct  = tau_percentile if tau_percentile is not None else int(base_cfg.get("tau_percentile", 90))
        out_root = (ndd_root / base_cfg["output_folder"]).resolve()

        # run_name encodes gene list + scope + chemistry;
        # the tau pipeline appends _tau{pct}_ in the folder name automatically.
        run_name = f"{gene_list_stem}_{scope}_{chemistry}"

        if tau_filtered:
            run_dir_name = f"{run_name}_tau{tau_pct}_threshold_{min_thr}_{today}"
        else:
            run_dir_name = f"{run_name}_threshold_{min_thr}_{today}"

        run_dir      = out_root / run_dir_name
        tmp_cfg_path = run_dir / "metadata" / "config_used.yaml"
        tmp_cfg_path.parent.mkdir(parents=True, exist_ok=True)

        cfg = dict(base_cfg)
        cfg["run_name"]       = run_name
        cfg["gene_list_path"] = str(gene_list_path)
        if tau_filtered and tau_percentile is not None:
            cfg["tau_percentile"] = tau_percentile
        _dump_yaml(cfg, tmp_cfg_path)

        print(f"\n  ▶ scope={scope}  chemistry={chemistry}  tau={tau_filtered}")
        print(f"    run_dir: {run_dir}")

        pipeline_ok = True
        try:
            if tau_filtered:
                from search_enrichment_gsea_tau_filtered import run_tau_filtered_pipeline
                run_tau_filtered_pipeline(config_path=str(tmp_cfg_path))
            else:
                from enrichment_pipeline_for_gene_list import run_gene_list_pipeline
                run_gene_list_pipeline(config_path=str(tmp_cfg_path))
        except Exception:
            print(f"  ✗ Pipeline raised an exception — skipping summary collection.")
            traceback.print_exc()
            pipeline_ok = False

        expected_summary = run_dir / _SUMMARY_REL[tau_filtered]
        print(f"    Expected summary: {expected_summary}")
        print(f"    Summary exists:   {expected_summary.exists()}")

        if not pipeline_ok:
            continue

        ges_results_dir = (ndd_root / base_cfg.get("ges_results_folder", "")).resolve()
        ges_cfg_path    = _locate_ges_config(ges_results_dir)

        rows = _collect_summary_rows(
            run_dir        = run_dir,
            tau_filtered   = tau_filtered,
            run_name       = run_name,
            scope          = scope,
            chemistry      = chemistry,
            gene_list_path = gene_list_path,
            gene_list_n    = gene_list_n,
            base_cfg       = base_cfg,
            ges_cfg_path   = ges_cfg_path,
        )
        all_rows.extend(rows)
        print(f"    ✔ Collected {len(rows)} rows.")

    if all_rows:
        out_root_default = (ndd_root / "results/enrichment_results").resolve()
        batch_path = out_root_default / f"{gene_list_stem}_batch_summary_tau_vs_v2_v3_{today}.csv"
        pd.DataFrame(all_rows).to_csv(batch_path, index=False)
        print(f"\n📌 Batch summary: {batch_path}")
    else:
        print(f"\n⚠️ No results collected for {gene_list_path.name}.")


# ---------------------------------------------------------------------------
# Main — accepts a single CSV or a folder of CSVs
# ---------------------------------------------------------------------------

def main(input_path: str, ndd_root: str, tau_percentile: int | None = None) -> None:
    p        = Path(input_path).resolve()
    ndd_root = Path(ndd_root).resolve()
    today    = datetime.datetime.now().strftime("%Y%m%d")

    if not p.exists():
        raise FileNotFoundError(f"Input not found: {p}")

    if tau_percentile is not None:
        print(f"Tau percentile override: {tau_percentile}")

    if p.is_dir():
        csv_files = sorted(p.glob("*.csv"))
        if not csv_files:
            raise ValueError(f"No .csv files found in: {p}")
        print(f"Found {len(csv_files)} gene list(s) in {p}.")
        for gene_list_path in csv_files:
            _run_one(gene_list_path, ndd_root, today, tau_percentile)
    else:
        _run_one(p, ndd_root, today, tau_percentile)

    print("\nDONE.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run 12-variant GSEA enrichment (v2/v3 × tau/no-tau × all scopes)."
    )
    parser.add_argument("input", help="Gene-list CSV or folder of CSVs.")
    parser.add_argument(
        "--ndd-root",
        default="/miridan-data/annaludmir/ndd_gene_modules",
        help="Project root directory (default: /miridan-data/annaludmir/ndd_gene_modules).",
    )
    parser.add_argument(
        "--tau-percentile",
        type=int,
        default=None,
        metavar="N",
        help="Override tau_percentile for all tau-filtered runs (e.g. 90). "
             "If omitted, each config file's own value is used.",
    )
    args = parser.parse_args()
    main(args.input, args.ndd_root, args.tau_percentile)
