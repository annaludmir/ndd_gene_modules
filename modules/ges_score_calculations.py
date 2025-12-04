import time
import yaml
import os
import h5py
from tqdm import tqdm
import datetime
import hashlib
import shutil
from pathlib import Path
from typing import Iterable, List, Optional, Any

import numpy as np
import pandas as pd
import scanpy as sc
from statsmodels.stats.multitest import multipletests


def load_config(config_path: str) -> dict:
    """
    Load a YAML configuration file for the GES pipeline.

    The YAML must define:
        - name_of_run
        - ndd_gene_modules_folder_root
        - data_path          (relative or absolute)
        - output_folder      (relative or absolute)
        - column_conditions  (mapping: column -> spec)

    data_type is optional; if missing, defaults to 'cortex'.
    """
    config_path = Path(config_path).resolve()

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    required = [
        "name_of_run",
        "ndd_gene_modules_folder_root",
        "data_path",
        "output_folder",
        "column_conditions",
    ]
    for key in required:
        if key not in config:
            raise ValueError(f"Missing required parameter in config: {key}")

    root = Path(config["ndd_gene_modules_folder_root"]).expanduser().resolve()

    def resolve(p: str | Path) -> str:
        p = Path(p)
        return str(p if p.is_absolute() else (root / p).resolve())

    # Resolve paths
    config["data_path"] = resolve(config["data_path"])
    config["output_folder"] = resolve(config["output_folder"])
    config["config_path"] = str(config_path)  # keep for metadata copy

    # Optional: allow data_type to be omitted
    if "data_type" not in config:
        config["data_type"] = "cortex"

    # We only make sure the base output folder exists;
    # the per-run subfolder is created inside run_ges_pipeline.
    os.makedirs(config["output_folder"], exist_ok=True)

    # Sanity check for data_path
    if not os.path.exists(config["data_path"]):
        raise FileNotFoundError(f"Data file not found: {config['data_path']}")

    return config


def get_obs_columns(h5ad_path: str) -> List[str]:
    """Return a list of column names stored in adata.obs using h5py only."""
    with h5py.File(h5ad_path, "r") as f:
        if "obs" not in f:
            return []

        obs_group = f["obs"]

        # Columns appear as datasets within /obs
        col_names = [
            key for key, value in obs_group.items()
            if (isinstance(value, h5py.Dataset) or isinstance(value, h5py.Group)) and key != "_index"
        ]

        return col_names

def get_obs_unique_values(h5ad_path: str, column: str):
    """
    Robust extraction of obs[column] values from any AnnData H5AD file.
    Supports both old and new storage formats.
    """

    with h5py.File(h5ad_path, "r") as f:
        if "obs" not in f:
            raise KeyError("This file contains no 'obs' group.")

        obs_group = f["obs"]

        if column not in obs_group:
            raise KeyError(f"Column '{column}' not found in obs.")

        node = obs_group[column]

        # ----------------------------------------------------------
        # Case 1: Old-style AnnData where obs[column] is a dataset
        # ----------------------------------------------------------
        if isinstance(node, h5py.Dataset):
            data = node[:]

        # ----------------------------------------------------------
        # Case 2: AnnData >= 0.8 where obs[column] is a GROUP
        # ----------------------------------------------------------
        elif isinstance(node, h5py.Group):

            # 2A — String categories encoding (most common)
            if "categories" in node:
                # e.g., string columns encoded as categorical
                data = node["categories"][:]

            # 2B — Numerical/boolean encoded as data + mask
            elif "data" in node:
                data = node["data"][:]

            else:
                raise TypeError(
                    f"obs/{column} is a group but does not contain 'categories' or 'data'. "
                    f"Keys present: {list(node.keys())}"
                )
        else:
            raise TypeError(f"Unsupported obs format for column '{column}'.")

    # ---- Decode bytes to strings ----
    if isinstance(data[0], bytes):
        data = [d.decode("utf-8") for d in data]

    # Return unique sorted values
    return sorted(pd.unique(data).tolist())


def compute_file_hash(path: str, algorithm: str = "sha256") -> str:
    """
    Compute a hash of the given file (default: sha256).
    """
    h = hashlib.new(algorithm)
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


# -------------------------------------------------------------------
# Utility helpers
# -------------------------------------------------------------------

def normalize_label(label) -> str:
    """
    Normalize cell-type labels for safe file naming.

    Handles:
    - bytes -> utf-8 string
    - strings like \"b'Neuron'\" -> \"Neuron\"
    """
    if isinstance(label, bytes):
        label = label.decode("utf-8")

    label_str = str(label)
    if label_str.startswith("b'") and label_str.endswith("'"):
        label_str = label_str[2:-1]

    return label_str


def load_and_preprocess_adata(data_path: str, chemistry: str | None = "v3"):
    """
    Generic adata loader + preprocessing.
    Works for any dataset — no cortex/all-layers branching.

    Parameters
    ----------
    data_path : str
        Path to the .h5ad file.
    chemistry : str or None
        If provided and a 'Chemistry' column exists in adata.obs,
        keeps only cells where Chemistry == chemistry.
        If None, no filtering is done.
    """
    import time
    import scanpy as sc

    print(f"Uploading data from: {data_path}")
    start_time = time.time()

    # Load any h5ad
    adata = sc.read_h5ad(data_path)
    print(f"Data loaded ({adata.n_obs} cells, {adata.n_vars} genes) "
          f"in {time.time() - start_time:.2f} seconds.")

    # -------------------------
    # Optional Chemistry filter
    # -------------------------
    if chemistry is not None and "Chemistry" in adata.obs.columns:
        print(f"Filtering only {chemistry} chemistry.")
        before = adata.n_obs
        adata = adata[adata.obs["Chemistry"] == chemistry, :]
        after = adata.n_obs
        print(f"Filtered Chemistry == '{chemistry}': {before} → {after} cells")

    # -------------------------
    # Preprocessing
    # -------------------------
    prepro_time = time.time()

    print('Starting data preprocessing...')

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    print(f"Preprocessing done in {time.time() - prepro_time:.2f} seconds.")
    print(f"Final adata: {adata.n_obs} cells, {adata.n_vars} genes")

    return adata



def compute_derived_column(adata, column_name, derived_spec):
    """
    Create a new column in adata.obs from logical rules defined in the YAML.
    Example:
       derived_spec = {
         "type": "categorical",
         "rules": {
             "RG": "CellType == 'Radial glia'",
             "IPCs": "CellType == 'Neuronal IPC'"
         }
       }
    """
    rules = derived_spec["rules"]
    new_col = []

    # Build the new column row-by-row
    for idx in range(adata.n_obs):
        row = adata.obs.iloc[idx]

        assigned = None
        for label, expr in rules.items():
            # Safely evaluate expression using row namespace
            if eval(expr, {}, dict(row)):
                assigned = label
                break

        new_col.append(assigned)

    adata.obs[column_name] = new_col

# -------------------------------------------------------------------
# Core GES functions (logic kept the same as original)
# -------------------------------------------------------------------

def calculate_ges(
    adata,
    condition_col: str,
    target_cell_type,
    chemistry: str,
    data_type: str,
):
    """
    Calculate the Generalized Expression Specificity (GES) score.

    Parameters
    ----------
    adata : AnnData
        Single-cell data (already preprocessed).
    condition_col : str
        Column in `adata.obs` that contains cell-type labels.
    target_cell_type : str or bytes
        Cell type to calculate specificity for.
    chemistry : str
        A label used only for file naming/logging (kept for compatibility).
    data_type : {'cortex', 'data_all', ...}
        Used to decide how to get gene names.

    Returns
    -------
    results_df : pd.DataFrame
        DataFrame with gene names, GES scores, mean expressions,
        and percentages of expressing cells.
    """
    ges_time = time.time()
    total_cells = adata.shape[0]

    target_idx = adata.obs[condition_col] == target_cell_type
    target_fraction = target_idx.sum() / total_cells
    filtered_data = adata[target_idx]

    # Mean expression in the target cell type
    mean_target_expr = filtered_data.X.mean(axis=0).A1
    weighted_target_expr = (1 - target_fraction) * mean_target_expr

    # Weighted mean expression in other cell types
    cell_types = adata.obs[condition_col].unique()
    other_cell_types = [ct for ct in cell_types if ct != target_cell_type]
    weighted_sum = np.zeros(adata.shape[1])
    per_exp_other = {}

    for ct in other_cell_types:
        ct_idx = adata.obs[condition_col] == ct
        ct_idx = np.array(ct_idx)
        ct_fraction = ct_idx.sum() / total_cells
        ct_mean_expr = adata[ct_idx].X.mean(axis=0).A1
        weighted_sum += ct_fraction * ct_mean_expr
        if len(other_cell_types) > 1:
            per_exp_other[f"per_exp_{ct}"] = (adata.X[ct_idx] > 0).mean(axis=0).A1

    # Calculate GES scores
    epsilon = 1e-10
    ges_scores = weighted_target_expr / (weighted_sum + epsilon)

    # Recompute masks for filtered data - other
    filtered_other_data = adata[~target_idx]

    # Calculate mean expression
    mean_expression_other = filtered_other_data.X.mean(axis=0).A1

    # Calculate percentages of expressing cells
    percent_expressed_target = (filtered_data.X > 0).mean(axis=0).A1
    percent_expressed_other = (filtered_other_data.X > 0).mean(axis=0).A1

    # TODO: maybe make it more generic
    if data_type == "cortex":
        genes = adata.var_names
    else:
        # in the big data gene symbols are in var["Gene"]
        genes = adata.var["Gene"]

    # Create DataFrame with results
    results_df = pd.DataFrame(
        {
            "gene": genes,
            "ges_score": ges_scores,
            "mean_expression_target": mean_target_expr,
            "mean_expression_other": mean_expression_other,
            "per_expressed_target": percent_expressed_target,
            "per_expressed_other": percent_expressed_other,
        }
    )

    if len(per_exp_other) > 0:
        per_expressed_df = pd.DataFrame.from_dict(per_exp_other)
        per_expressed_df["gene"] = adata.var_names
        results_df = results_df.merge(per_expressed_df, on="gene")

    print(f"finished_ges_caculations in {time.time() - ges_time:.2f} seconds")

    return results_df


def calculate_ges_with_permutations(
    adata,
    condition_col: str,
    target_cell_type,
    chemistry: str,
    data_type: str,
    expression_threshold: float = 0.05,
    n_permutations: int = 500,
):
    """
    Calculate GES scores and perform permutation-based p-value estimation.

    Parameters
    ----------
    adata : AnnData
    condition_col : str
    target_cell_type : str or bytes
    chemistry : str
    data_type : str
    expression_threshold : float, optional
        Minimum fraction of cells in target population expressing each gene.
    n_permutations : int, optional
        Number of label permutations.

    Returns
    -------
    actual_results : pd.DataFrame
        GES results with p-values and FDR-adjusted p-values.
    """

    # Subset the data for the target cell type
    target_mask = adata.obs[condition_col] == target_cell_type
    target_data = adata[target_mask]

    # Calculate the number of cells expressing each gene in the target cell type
    gene_expression_num = (target_data.X > 0).sum(axis=0).A1

    # Filter genes with at least `expression_threshold` proportion of expression
    expressed_genes_mask = gene_expression_num >= (target_data.shape[0] * expression_threshold)
    filtered_adata = adata[:, expressed_genes_mask]

    actual_results = calculate_ges(
        filtered_adata, condition_col, target_cell_type, chemistry, data_type
    )
    actual_results.to_csv(f"ges_only_{chemistry}_{normalize_label(target_cell_type)}.csv")
    actual_scores = actual_results["ges_score"].values

    # Initialize array to store permuted scores
    permuted_scores = np.zeros((n_permutations, len(actual_scores)))
    p_values = np.zeros(len(actual_scores))

    # Perform permutations
    for i in range(n_permutations):
        print("permutation ", i + 1)
        # Shuffle the condition column
        shuffled_conditions = np.random.permutation(filtered_adata.obs[condition_col])
        filtered_adata.obs[condition_col] = shuffled_conditions

        # Recalculate GES scores for the shuffled data
        permuted_results = calculate_ges(
            filtered_adata, condition_col, target_cell_type, chemistry, data_type
        )
        permuted_scores[i] = permuted_results["ges_score"].values
        p_values += np.array((actual_scores < permuted_scores[i]).astype(int).tolist())

    # Add p-values to results
    actual_results["p_values"] = p_values / n_permutations

    # subset the data for genes with more than expression_threshold in the cell population
    actual_results = actual_results.loc[
        actual_results["per_expressed_target"] >= expression_threshold, :
    ]

    # FDR correction
    actual_results["adj_p_val"] = multipletests(
        actual_results["p_values"], method="fdr_bh"
    )[1]

    return actual_results


# -------------------------------------------------------------------
# High-level driver
# -------------------------------------------------------------------
def run_ges_pipeline(config_path: str):
    """
    Entry point: run the GES pipeline using a YAML configuration file.

    Layout of outputs:

    <output_folder>/
      <name_of_run>_<YYYYMMDD_HHMMSS>/
        metadata/
          config.yaml           # copy of the config used
          adata_hash.txt        # hash of the .h5ad file in data_path
        data/
          ges_spec_<chemistry>_<column>_<target>.csv
          ges_spec_<chemistry>_<column>_<target>_perm.csv (optional)

    The YAML must define:
        - name_of_run
        - ndd_gene_modules_folder_root
        - data_path
        - output_folder
        - column_conditions

    Optional:
        - data_type
        - chemistry
        - expression_threshold
        - permutations
        - n_permutations
    """

    # ------------------------------
    # Load YAML config
    # ------------------------------
    config = load_config(config_path)

    name_of_run = config["name_of_run"]
    data_type = config["data_type"]
    data_path = config["data_path"]
    output_root = Path(config["output_folder"])
    column_conditions = config["column_conditions"]  # dict: {column: spec}

    chemistry = config.get("chemistry")
    expression_threshold = config.get("expression_threshold", 0.05)
    permutations = config.get("permutations", False)
    n_permutations = config.get("n_permutations", 500)

    # ------------------------------
    # Create run-specific folder structure
    # ------------------------------
    timestamp = datetime.datetime.now().strftime("%Y%m%d")
    run_dir = output_root / f"{name_of_run}_{timestamp}"
    metadata_dir = run_dir / "metadata"
    data_dir = run_dir / "data"

    metadata_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    # Copy configuration file into metadata
    cfg_dest = metadata_dir / Path(config["config_path"]).name
    shutil.copy2(config["config_path"], cfg_dest)

    # Compute and store hash of the h5ad data file
    hash_value = compute_file_hash(data_path, algorithm="sha256")
    hash_file_path = metadata_dir / "adata_hash.txt"
    with open(hash_file_path, "w") as hf:
        hf.write(f"file: {data_path}\n")
        hf.write("algorithm: sha256\n")
        hf.write(f"hash: {hash_value}\n")

    print("\n==============================")
    print("  Loading GES YAML config")
    print("==============================")
    print(f"• Name of run:          {name_of_run}")
    print(f"• Data path:            {data_path}")
    print(f"• Output root:          {output_root}")
    print(f"• Run directory:        {run_dir}")
    print(f"• Column_conditions:    {column_conditions}")
    if chemistry is not None:
        print(f"• Chemistry:            {chemistry}")
    print(f"• Expression threshold: {expression_threshold}")
    print(f"• Permutations:         {permutations}")
    if permutations:
        print(f"• Num permutations:     {n_permutations}")
    print("==============================\n")

    print("\nAvailable columns in .obs:")
    obs_columns = get_obs_columns(data_path)
    print(obs_columns)
    for column in column_conditions:
      if not (isinstance(column_conditions[column], dict) and "derived" in column_conditions[column]):
        column_unique_values= get_obs_unique_values(data_path, column)
        print(f'In column {column}, there are the following unique values:')
        print(column_unique_values)


    # ------------------------------
    # Load & preprocess data
    # ------------------------------
    adata = load_and_preprocess_adata(
    data_path = data_path,
    chemistry = chemistry,
)

    # ------------------------------
    # Process each (possibly derived) column
    # ------------------------------
    for condition_col, spec in tqdm(column_conditions.items()):
        print(f"\n=== Column: {condition_col} ===")

        # Support both list-of-conditions and dict-spec
        if isinstance(spec, dict):
            # Derived column? (boolean_expr / rules etc. handled earlier in your code)
            if "derived" in spec and "boolean_expr" in spec["derived"]:
                expr = spec["derived"]["boolean_expr"]
                print(f"  Detected BOOLEAN derived column '{condition_col}' with expr: {expr}")
                mask = adata.obs.eval(expr)
                mask = mask.fillna(False).astype(bool)
                adata.obs[condition_col] = mask
                condition_list = spec.get("conditions", [True])
            elif "conditions" in spec:
                condition_list = spec["conditions"]
            else:
                raise ValueError(
                    f"Column '{condition_col}' has a dict spec but no 'conditions' or 'derived'."
                )
        else:
            # Backwards compatibility: spec is just a list of conditions
            condition_list = spec

        column_values = adata.obs[condition_col].unique().tolist()
        print(f"  Available values:   {column_values}")
        print(f"  Requested targets:  {condition_list}")

        for target in condition_list:
            if target not in column_values:
                print(f"  ⚠️ '{target}' NOT FOUND in column '{condition_col}' → skipping")
                continue

            print(f"\n  → Processing target: {target}")

            target_mask = adata.obs[condition_col] == target
            target_data = adata[target_mask]

            if target_data.n_obs == 0:
                print(f"  ⚠️ target '{target}' has zero cells → skipping")
                continue

            # Filter genes by expression in this target group
            gene_expression_num = (target_data.X > 0).sum(axis=0).A1
            expressed_mask = gene_expression_num >= (
                target_data.shape[0] * expression_threshold
            )

            if expressed_mask.sum() == 0:
                print(f"  ⚠️ No genes pass expression threshold → skipping")
                continue

            filtered_adata = adata[:, expressed_mask]

            # Compute GES
            ges_results = calculate_ges(
                filtered_adata,
                condition_col,
                target,
                chemistry if chemistry is not None else "",
                data_type,
            )

            ges_results = ges_results.sort_values("ges_score", ascending=False)
            target_name = normalize_label(target)

            # NOTE: we **do not** include data_type in the file name anymore
            ges_filename = f"ges_spec_{condition_col}_{target_name}.csv"
            ges_path = data_dir / ges_filename

            ges_results.to_csv(ges_path, index=False)
            print(f"  ✔ Saved GES results → {ges_path}")

            # --------- permutations (optional) ----------
            if permutations:
                print(f"  🔁 Running {n_permutations} permutations for {target_name}")

                perm_results = calculate_ges_with_permutations(
                    filtered_adata,
                    condition_col,
                    target,
                    chemistry if chemistry is not None else "",
                    data_type,
                    expression_threshold=expression_threshold,
                    n_permutations=n_permutations,
                )

                perm_filename = (
f"ges_spec_{condition_col}_{target_name}_perm.csv"
                )
                perm_path = data_dir / perm_filename

                perm_results.to_csv(perm_path, index=False)
                print(f"  ✔ Saved permutation results → {perm_path}")

    print("\n🎉 DONE — GES pipeline completed.\n")
    print(f"All outputs for this run are under:\n  {run_dir}\n")



# -------------------------------------------------------------------
# Optional: allow script execution with defaults
# -------------------------------------------------------------------

def _main_():
    # Runs with default settings from constants above.
    # No CLI args are used; this is just a convenience.
    run_ges_pipeline()

if __name__ == '__main__':
  _main_()