import os
import re
import sys
import shutil
import datetime
from pathlib import Path

import yaml
import pandas as pd

def _extract_context_suffix(run_name: str) -> str:
    """
    Extract and normalize biological context from run_name.

    Order rule:
        [layer scope]_[analysis context]

    Where layer scope is:
        - cortex
        - all_layers

    Examples:
        ID_microcephaly_behvioral_cortex_cell_phase
          -> cortex_cell_phase

        ID_microcephaly_behvioral_all_layers
          -> all_layers

        ID_microcephaly_behvioral_cortex_cell_phase_all_layers
          -> all_layers_cell_phase
    """

    # Detect layer scope
    layer = None
    if "all_layers" in run_name:
        layer = "all_layers"
    elif "cortex" in run_name:
        layer = "cortex"
    else:
        raise ValueError(f"No layer scope (cortex/all_layers) in run_name: {run_name}")

    # Remove gene list prefix
    tail = run_name

    # Remove everything before first cortex / all_layers
    tail = re.sub(r"^.*?(cortex|all_layers)", r"\1", tail)

    # Remove layer words from tail to get biological context
    context = tail.replace("cortex", "").replace("all_layers", "")
    context = context.strip("_")

    # Rebuild in canonical order
    if context:
        return f"{layer}_{context}"
    else:
        return layer


def _safe_load_yaml(path: Path) -> dict:
    with path.open("r") as f:
        return yaml.safe_load(f)


def _safe_dump_yaml(cfg: dict, path: Path) -> None:
    with path.open("w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)


def _find_nes_column(df: pd.DataFrame) -> str | None:
    # common possibilities we’ve seen
    candidates = ["NES", "nes", "nes_score", "NES score", "Normalized Enrichment Score"]
    for c in candidates:
        if c in df.columns:
            return c

    # fallback: case-insensitive match
    low = {c.lower(): c for c in df.columns}
    for key in ["nes", "nes_score"]:
        if key in low:
            return low[key]
    return None


def _is_significant_run(run_dir: Path, nes_threshold: float) -> bool:
    """
    Decide whether to keep results:
    - If summary exists and ANY NES > threshold -> keep
    - If we can't find summary / NES column -> keep (fail-open so we don't delete useful outputs)
    """
    # Adjust these if your pipeline writes to a different location
    summary_path = run_dir / "data" / "enrichment_results" / "GSEA" / "GSEA_final_summary.csv"
  
    if summary_path is None:
        print(f"⚠️ Could not find GSEA summary CSV under {run_dir} -> keeping results.")
        return True

    df = pd.read_csv(summary_path)
    nes_col = _find_nes_column(df)
    if nes_col is None:
        print(f"⚠️ Could not find NES column in {summary_path.name} -> keeping results.")
        return True

    # coerce to numeric; ignore non-numeric values
    nes = pd.to_numeric(df[nes_col], errors="coerce")
    max_nes = nes.max(skipna=True)

    print(f"   ↳ Max NES found: {max_nes}")
    if pd.isna(max_nes):
        print("⚠️ NES values are all NaN -> keeping results.")
        return True

    return bool(max_nes > nes_threshold)


def main(base_cfg_path, gene_lists_folder):
    """
    Usage:
      python -u enrichment_cal_lists_loop.py <base_config.yaml> <gene_lists_folder>

    gene_lists_folder should contain many *.csv files (each with a 'gene' column, etc.)
    """
    base_cfg_path = Path(base_cfg_path).resolve()
    gene_lists_folder = Path(gene_lists_folder).resolve()
  
    if not base_cfg_path.exists():
        raise FileNotFoundError(f"Base config not found: {base_cfg_path}")
    if not gene_lists_folder.exists():
        raise FileNotFoundError(f"Gene lists folder not found: {gene_lists_folder}")

    base_cfg = _safe_load_yaml(base_cfg_path)

    # Flags (default behavior if missing)
    keep_only_significant = bool(base_cfg.get("keep_only_significant", False))
    nes_keep_threshold = float(base_cfg.get("nes_keep_threshold", 1))
    run_name_from_cfg = base_cfg.get("run_name")
    context_suffix = _extract_context_suffix(run_name_from_cfg)

    # Import your pipeline (adjust module path/name to your project)
    # Example assumes: modules/enrichment_pipeline_for_gene_list.py defines run_gene_list_pipeline()
    from modules import enrichment_pipeline_for_gene_list as epfgl

    # The pipeline’s run folder naming convention (based on what you requested earlier):
    #   {run_name}_{min_threshold}_{YYYYMMDD}
    today = datetime.datetime.now().strftime("%Y%m%d")

    out_root = Path(base_cfg["ndd_gene_modules_folder_root"]) / base_cfg["output_folder"]
    out_root = out_root.resolve()

    min_thr = base_cfg.get("gsea", {}).get("min_ges_score_threshold", None)
    if min_thr is None:
        raise ValueError("Base config missing gsea.min_ges_score_threshold")

    # iterate gene list csvs
    csv_files = sorted(gene_lists_folder.glob("*.csv"))
    if not csv_files:
        raise ValueError(f"No .csv files found in {gene_lists_folder}")

    print(f"Found {len(csv_files)} gene lists in: {gene_lists_folder}")
    print(f"keep_only_significant={keep_only_significant} | nes_keep_threshold={nes_keep_threshold}")

    for csv_path in csv_files:
        gene_list_name = csv_path.stem  # filename only, no extension

        run_name = f"{gene_list_name}_{context_suffix}"
        # you asked to set gene_list_path = "data/genes/{gene_list_file_name}.csv"
        # BUT in batch mode we’ll point to the *actual* csv_path so it always works.
        # If you *must* force "data/genes/...", copy the file there first.
        gene_list_path = str(csv_path)

        cfg = dict(base_cfg)
        cfg["run_name"] = run_name
        cfg["gene_list_path"] = gene_list_path

        # Write a temp config next to the outputs (so you can inspect later if needed)
        run_dir = out_root / f"{run_name}_threshold_{min_thr}_{today}"
        tmp_cfg_path = run_dir / "metadata" / "config_used.yaml"
        tmp_cfg_path.parent.mkdir(parents=True, exist_ok=True)
        _safe_dump_yaml(cfg, tmp_cfg_path)

        print("\n============================================================")
        print(f"▶ Running gene list: {csv_path.name}")
        print(f"   run_name: {run_name}")
        print(f"   run_dir:  {run_dir}")
        print("============================================================")

        # Run the pipeline using the saved temp config
        try:
          epfgl.run_gene_list_pipeline(config_path=str(tmp_cfg_path))

          # Decide keep/delete
          if keep_only_significant:
              keep = _is_significant_run(run_dir, nes_keep_threshold)
              if not keep:
                  print(f"🧹 No NES > {nes_keep_threshold}. Deleting: {run_dir}")
                  shutil.rmtree(run_dir, ignore_errors=True)
              else:
                  print(f"✅ Significant (NES > {nes_keep_threshold}). Keeping: {run_dir}")
          else:
              print("ℹ️ keep_only_significant=false -> keeping run folder")
        except LookupError:
          print("Probably not enough relevant genes were found for encrichment.")

    print("\nDONE.")


if __name__ == "__main__":
    main()
