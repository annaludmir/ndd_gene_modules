"""
Computes global Tau specificity scores for all genes, per condition column,
from an AnnData (.h5ad) file. Uses the same YAML config as the GES pipeline.

Tau_g = sum_k(1 - x_norm_{k,g}) / (n_groups - 1)
where x_norm_{k,g} = mean_expr_{k,g} / max_k(mean_expr_{k,g})

Output per column: tau_scores_{column}.csv
  gene | tau | max_group | mean_expr_{group1} | mean_expr_{group2} | ...
"""

import argparse
import datetime
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from specificity_score_calculations import (
    apply_derived_columns,
    compute_file_hash,
    get_gene_names,
    get_obs_columns,
    load_and_preprocess_adata,
    load_config,
    normalize_label,
)


def _count_nonzero_per_gene(X) -> np.ndarray:
    if sp.issparse(X):
        return np.asarray((X > 0).sum(axis=0)).ravel()
    return np.asarray((X > 0).sum(axis=0)).ravel()


def compute_global_tau(
    adata,
    condition_col: str,
    genes: list,
    min_cells_global: int = 5,
    agg: str = "mean",
) -> pd.DataFrame:
    """
    Compute global Tau specificity for every gene across all groups in condition_col.

    Tau is in [0, 1]: 1 = expressed in exactly one group; 0 = uniform across all groups.
    The max_group column indicates which group has the highest mean expression per gene.

    Returns
    -------
    pd.DataFrame
        Columns: gene, tau, max_group, mean_expr_{group} for each group.
    """
    t0 = time.time()

    global_counts = _count_nonzero_per_gene(adata.X)
    keep_mask = global_counts >= min_cells_global
    n_keep = int(keep_mask.sum())

    adata_f = adata[:, keep_mask].copy()
    genes_f = [g for g, k in zip(genes, keep_mask) if k]

    labels = adata_f.obs[condition_col].astype(str).fillna("NA")
    adata_f.obs[condition_col] = labels.values
    group_names = list(pd.unique(adata_f.obs[condition_col]))
    n_groups = len(group_names)

    if n_groups < 2:
        raise ValueError(f"Need >=2 groups in '{condition_col}' to compute Tau (got {n_groups})")

    expr_by_group = np.zeros((n_groups, adata_f.n_vars), dtype=float)
    for i, grp in enumerate(group_names):
        idx = (adata_f.obs[condition_col] == grp).to_numpy()
        if idx.sum() == 0:
            continue
        sub_X = adata_f[idx].X
        if agg == "mean":
            expr_by_group[i] = np.asarray(sub_X.mean(axis=0)).ravel()
        elif agg == "median":
            dense = sub_X.toarray() if sp.issparse(sub_X) else np.asarray(sub_X)
            expr_by_group[i] = np.median(dense, axis=0)
        else:
            raise ValueError(f"Unknown agg='{agg}'. Use 'mean' or 'median'.")

    eps = 1e-12
    x_max = expr_by_group.max(axis=0)
    x_max_safe = np.where(x_max > 0, x_max, eps)
    x_norm = expr_by_group / x_max_safe
    tau = np.sum(1.0 - x_norm, axis=0) / (n_groups - 1)

    max_group_idx = np.argmax(expr_by_group, axis=0)
    max_group = [group_names[j] for j in max_group_idx]

    result = pd.DataFrame({"gene": genes_f, "tau": tau, "max_group": max_group})
    for i, grp in enumerate(group_names):
        result[f"mean_expr_{grp}"] = expr_by_group[i]

    print(
        f"compute_global_tau | col='{condition_col}' | "
        f"genes_kept={n_keep}/{adata.n_vars} | "
        f"n_groups={n_groups} | elapsed={time.time() - t0:.2f}s"
    )
    return result


def run_tau_pipeline(config_path: str) -> Path:
    """
    Run the Tau pipeline using a GES-format YAML config.
    For each column in column_conditions, compute global Tau and save to CSV.

    Returns the run output directory.
    """
    config = load_config(config_path)

    name_of_run = config["name_of_run"]
    data_path = config["data_path"]
    output_root = Path(config["output_folder"])
    column_conditions = config["column_conditions"]

    normalize_data = bool(config.get("normalize_data", True))
    chemistry = config.get("chemistry", None)
    min_cells_global = int(config.get("min_cells_global", 5))
    species = config.get("species", "human")
    agg = config.get("tau_agg", "mean")

    date_tag = datetime.datetime.now().strftime("%Y%m%d")
    run_dir = output_root / f"{name_of_run}_{date_tag}"
    metadata_dir = run_dir / "metadata"
    data_dir = run_dir / "data"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    shutil.copy2(config["config_path"], metadata_dir / Path(config["config_path"]).name)

    hash_value = compute_file_hash(data_path)
    with open(metadata_dir / "adata_hash.txt", "w") as hf:
        hf.write(f"file: {data_path}\nalgorithm: sha256\nhash: {hash_value}\n")

    print(f"\n{'='*50}")
    print("  Tau Specificity Pipeline")
    print(f"{'='*50}")
    print(f"• Run name:        {name_of_run}")
    print(f"• Data path:       {data_path}")
    print(f"• Output dir:      {run_dir}")
    print(f"• Normalize data:  {normalize_data}")
    print(f"• Aggregation:     {agg}")
    print(f"• Min cells global:{min_cells_global}")
    print(f"{'='*50}\n")

    print("Available .obs columns (fast h5py scan):")
    print(get_obs_columns(data_path))

    adata = load_and_preprocess_adata(data_path, chemistry=chemistry, normalize_data=normalize_data)
    apply_derived_columns(adata, column_conditions)
    genes = get_gene_names(adata, species=species)

    if len(genes) != adata.n_vars:
        raise ValueError(
            f"Gene list length ({len(genes)}) does not match adata.n_vars ({adata.n_vars})."
        )

    for condition_col in column_conditions:
        print(f"\n=== Column: {condition_col} ===")

        if condition_col not in adata.obs.columns:
            print(f"  ⚠️ '{condition_col}' not in adata.obs — skipping")
            continue

        adata.obs[condition_col] = adata.obs[condition_col].apply(normalize_label)

        tau_df = compute_global_tau(
            adata,
            condition_col,
            genes,
            min_cells_global=min_cells_global,
            agg=agg,
        )

        out_csv = data_dir / f"tau_scores_{condition_col}.csv"
        tau_df.to_csv(out_csv, index=False)
        print(f"  ✔ Saved → {out_csv}")

    print(f"\n✅ DONE. All outputs under: {run_dir}\n")
    return run_dir


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compute global Tau specificity scores for all genes per condition column. "
            "Accepts the same YAML config format as the GES pipeline."
        )
    )
    parser.add_argument(
        "config",
        help="Path to GES-format YAML config file.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    run_tau_pipeline(args.config)
    sys.exit(0)
