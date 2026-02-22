import time
import yaml
import os
import re
import h5py
import datetime
import hashlib
import shutil
from pathlib import Path
from typing import List, Optional, Any, Dict

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests

# Optional dependency: mygene for Ensembl->symbol conversion
try:
    from mygene import MyGeneInfo
    _HAS_MYGENE = True
    mg = MyGeneInfo()
except Exception:
    _HAS_MYGENE = False
    mg = None


# -------------------------------------------------------------------
# Config loading
# -------------------------------------------------------------------

def load_config(config_path: str) -> dict:
    """
    Load a YAML configuration file for the GES pipeline.

    Required keys:
        - name_of_run
        - ndd_gene_modules_folder_root
        - data_path          (relative or absolute)
        - output_folder      (relative or absolute)
        - column_conditions  (mapping: column -> spec)
    """
    config_path = Path(config_path).expanduser().resolve()

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

    def resolve_path(p: str | Path) -> str:
        p = Path(p).expanduser()
        return str(p if p.is_absolute() else (root / p).resolve())

    # Resolve paths to absolute
    config["data_path"] = resolve_path(config["data_path"])
    config["output_folder"] = resolve_path(config["output_folder"])
    config["config_path"] = str(config_path)  # keep for metadata copy

    os.makedirs(config["output_folder"], exist_ok=True)

    if not os.path.exists(config["data_path"]):
        raise FileNotFoundError(f"Data file not found: {config['data_path']}")

    if not isinstance(config["column_conditions"], dict):
        raise TypeError("column_conditions must be a dict mapping column_name -> spec")

    return config


# -------------------------------------------------------------------
# Fast H5AD inspection helpers (h5py-only)
# -------------------------------------------------------------------

def get_obs_columns(h5ad_path: str) -> List[str]:
    """Return a list of column names stored in adata.obs using h5py only."""
    with h5py.File(h5ad_path, "r") as f:
        if "obs" not in f:
            return []
        obs_group = f["obs"]
        return [
            key for key, value in obs_group.items()
            if (isinstance(value, h5py.Dataset) or isinstance(value, h5py.Group)) and key != "_index"
        ]


def get_obs_unique_values(h5ad_path: str, column: str) -> List[Any]:
    """
    Robust extraction of obs[column] values from AnnData H5AD using h5py only.
    Supports both old and new storage formats.
    """
    with h5py.File(h5ad_path, "r") as f:
        if "obs" not in f:
            raise KeyError("This file contains no 'obs' group.")
        obs_group = f["obs"]
        if column not in obs_group:
            raise KeyError(f"Column '{column}' not found in obs.")

        node = obs_group[column]

        # Case 1: dataset
        if isinstance(node, h5py.Dataset):
            data = node[:]

        # Case 2: group (categorical / sparse encoding)
        elif isinstance(node, h5py.Group):
            if "categories" in node:
                data = node["categories"][:]
            elif "data" in node:
                data = node["data"][:]
            else:
                raise TypeError(
                    f"obs/{column} is a group but does not contain 'categories' or 'data'. "
                    f"Keys present: {list(node.keys())}"
                )
        else:
            raise TypeError(f"Unsupported obs format for column '{column}'.")

    # Decode bytes
    if len(data) and isinstance(data[0], bytes):
        data = [d.decode("utf-8") for d in data]

    # Normalize "b'X'" style strings too
    data = [normalize_label(d) for d in data]

    return sorted(pd.unique(pd.Series(data)).tolist())


def compute_file_hash(path: str, algorithm: str = "sha256") -> str:
    """Compute a hash of the given file (default: sha256)."""
    h = hashlib.new(algorithm)
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


# -------------------------------------------------------------------
# Gene helpers: Ensembl -> symbol (optional)
# -------------------------------------------------------------------

def make_unique_gene_symbols(symbols: List[str]) -> List[str]:
    """Ensure gene symbols are unique by appending _dupN where necessary."""
    seen: Dict[str, int] = {}
    out: List[str] = []
    for s in symbols:
        if s not in seen:
            seen[s] = 0
            out.append(s)
        else:
            seen[s] += 1
            out.append(f"{s}_dup{seen[s]}")
    return out


def looks_like_ensembl(x: str) -> bool:
    # Typical human Ensembl gene IDs: ENSG..., ENSG... with version sometimes
    return bool(re.match(r"^ENS[A-Z]*G\d+(\.\d+)?$", str(x)))


def convert_ensembl_to_symbol(gene_list: List[str], species: str = "human") -> List[str]:
    """
    Convert Ensembl IDs -> gene symbols using mygene.
    If mygene is unavailable, returns the original list.
    """
    if not _HAS_MYGENE:
        return gene_list

    ensembls = [g for g in gene_list if looks_like_ensembl(g)]
    if not ensembls:
        return gene_list

    res = mg.querymany(
        ensembls,
        scopes="ensembl.gene",
        fields="symbol",
        species=species,
        as_dataframe=True,
    )

    # mygene dataframe: index is query by default when as_dataframe=True,
    # but depending on version it might include 'query' column. Handle both.
    if "symbol" not in res.columns:
        return gene_list

    # If 'query' column exists use it; else use index
    if "query" in res.columns:
        query_keys = res["query"]
        symbols = res["symbol"]
        mapping = dict(zip(query_keys, symbols.fillna(query_keys)))
    else:
        mapping = res["symbol"].fillna(res.index.to_series()).to_dict()

    converted = [mapping.get(g, g) for g in gene_list]
    return make_unique_gene_symbols([str(x) for x in converted])


def get_gene_names(adata: sc.AnnData, species: str = "human") -> List[str]:
    """
    Safely extract gene names from an AnnData object.
    If names look like Ensembl IDs, convert to symbols (optional).
    """
    possible_columns = ["gene_symbol", "gene_symbols", "Gene", "gene", "name", "names"]

    genes: List[str]
    for col in possible_columns:
        if col in adata.var.columns:
            genes = adata.var[col].astype(str).tolist()
            break
    else:
        genes = adata.var_names.astype(str).tolist()

    # convert if needed
    if any(looks_like_ensembl(g) for g in genes):
        print("Detected Ensembl IDs in genes — converting to gene symbols...")
        genes = convert_ensembl_to_symbol(genes, species=species)

    return genes


# -------------------------------------------------------------------
# Label normalization
# -------------------------------------------------------------------

def normalize_label(label: Any) -> str:
    """
    Normalize labels for safe comparisons and file naming.

    Handles:
    - bytes -> utf-8 string
    - strings like "b'Neuron'" -> "Neuron"
    - strips surrounding whitespace
    """
    if isinstance(label, bytes):
        label = label.decode("utf-8")

    s = str(label).strip()

    # unwrap literal "b'X'" or 'b"X"'
    if (s.startswith("b'") and s.endswith("'")) or (s.startswith('b"') and s.endswith('"')):
        s = s[2:-1]

    return s


# -------------------------------------------------------------------
# Load + preprocess
# -------------------------------------------------------------------

def load_and_preprocess_adata(
    data_path: str,
    chemistry: Optional[str] = None,
    normalize_data: bool = True,
) -> sc.AnnData:
    """
    Generic adata loader + preprocessing.

    If chemistry is not None and adata.obs contains "Chemistry",
    filter to Chemistry == chemistry.
    """
    print(f"Uploading data from: {data_path}")
    start_time = time.time()
    adata = sc.read_h5ad(data_path)
    print(f"Data loaded ({adata.n_obs} cells, {adata.n_vars} genes) in {time.time() - start_time:.2f}s")

    print("Filtering chemistry...")
    if chemistry is not None and "Chemistry" in adata.obs.columns:
        before = adata.n_obs
        adata = adata[adata.obs["Chemistry"] == chemistry, :].copy()
        print(f"Filtered Chemistry == '{chemistry}': {before} → {adata.n_obs} cells")
    else:
      print("Chemistry is none or Chemistry column is not in data.")

    # preprocessing
    print("Preprocessing...")
    prepro_time = time.time()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    if normalize_data:
        print("Normalizing data...")
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

    print(f"Preprocessing done in {time.time() - prepro_time:.2f}s")
    return adata


# -------------------------------------------------------------------
# Derived column evaluation (vectorized)
# -------------------------------------------------------------------

def apply_derived_columns(adata: sc.AnnData, column_conditions: dict) -> None:
    """
    Mutates adata.obs in-place by creating any derived columns found in YAML spec.
    Supports:
      - derived.boolean_expr -> boolean column
      - derived.rules -> categorical column (ordered rules)
    """
    for col_name, spec in column_conditions.items():
        if not isinstance(spec, dict):
            continue
        if "derived" not in spec or not isinstance(spec["derived"], dict):
            continue

        d = spec["derived"]

        # Boolean derived column
        if "boolean_expr" in d:
            expr = d["boolean_expr"]
            print(f"Creating derived boolean column '{col_name}' from expr: {expr}")
            mask = adata.obs.eval(expr)
            # ensure boolean dtype
            mask = mask.fillna(False).astype(bool)
            adata.obs[col_name] = mask
            continue

        # Categorical derived column via rules
        if "rules" in d:
            rules = d["rules"]
            if not isinstance(rules, dict) or len(rules) == 0:
                raise ValueError(f"Derived rules for '{col_name}' must be a non-empty mapping")

            print(f"Creating derived categorical column '{col_name}' from {len(rules)} rules (ordered)")
            # Evaluate each rule vectorized
            masks = []
            labels = []
            for label, expr in rules.items():
                labels.append(label)
                m = adata.obs.eval(expr)
                m = m.fillna(False).astype(bool)
                masks.append(m.to_numpy())

            # np.select assigns first True rule; default None
            out = np.select(masks, labels, default=None)
            adata.obs[col_name] = out
            continue

# -------------------------------------------------------------------
# Core tau functions
# -------------------------------------------------------------------
def _count_nonzero_per_gene(X) -> np.ndarray:
    """#cells expressing gene (X>0) for each gene. Dense or sparse."""
    if sp.issparse(X):
        return np.asarray((X > 0).sum(axis=0)).ravel()
    return np.asarray((X > 0).sum(axis=0)).ravel()

def calculate_tau(
    adata: sc.AnnData,
    condition_col: str,
    target_cell_type: Any,
    genes: List[str],
    min_cells_global: int = 5,
    target_expr_frac: float = 0.05,
    agg: str = "mean",  # "mean" (default), or "median" (slower), or "pct90" (slower)
) -> pd.DataFrame:
    """
    Drop-in replacement for calculate_ges().

    Returns a DataFrame with the same columns as your GES function, but:
      - results_df["ges_score"] is replaced by a tau-based target specificity score:
          tau_target_score = tau * (x_target / x_max)

    Filtering (same spirit as your GES filtering):
      1) keep genes expressed in >= min_cells_global cells globally
      2) keep genes expressed in >= target_expr_frac of TARGET cells
    """
    t0 = time.time()

    # ---- basic checks ----
    if condition_col not in adata.obs.columns:
        raise ValueError(f"Column '{condition_col}' not found in adata.obs")

    total_cells = adata.n_obs
    target_mask = (adata.obs[condition_col] == target_cell_type).to_numpy()
    n_target = int(target_mask.sum())
    if n_target == 0:
        raise ValueError(f"Target '{target_cell_type}' has 0 cells in '{condition_col}'")

    # ---- filtering ----
    global_counts = _count_nonzero_per_gene(adata.X)
    global_keep = global_counts >= int(min_cells_global)

    target_X = adata.X[target_mask, :]
    target_counts = _count_nonzero_per_gene(target_X)
    target_keep = target_counts >= (float(target_expr_frac) * float(n_target))

    keep_mask = global_keep & target_keep
    if int(keep_mask.sum()) == 0:
        return pd.DataFrame(
            columns=[
                "gene",
                "ges_score",
                "tau",
                "tau_target_score",
                "mean_expression_target",
                "mean_expression_other",
                "per_expressed_target",
                "per_expressed_other",
                "max_group",
                "x_target",
                "x_max",
            ]
        )

    adata_f = adata[:, keep_mask].copy()
    genes_f = [g for g, keep in zip(genes, keep_mask) if keep]

    # ---- group definitions ----
    # normalize labels to strings so grouping is stable
    groups = pd.Series(adata_f.obs[condition_col].astype(str)).fillna("NA")
    adata_f.obs[condition_col] = groups.values

    group_names = list(pd.unique(adata_f.obs[condition_col]))
    n_groups = len(group_names)
    if n_groups < 2:
        raise ValueError(f"Need >=2 groups in '{condition_col}' to compute tau (got {n_groups})")

    # ---- compute x_{group,gene} matrix (aggregated expression per group) ----
    # We'll compute mean/median/pct90 per gene per group.
    # Mean is fast and aligns with pseudo-bulk thinking.
    X = adata_f.X

    def _agg_vec(Xsub) -> np.ndarray:
        if agg == "mean":
            return np.asarray(Xsub.mean(axis=0)).ravel()
        # For median/pct, convert to dense per-group; slower but workable if groups small.
        Xdense = Xsub.toarray() if sp.issparse(Xsub) else np.asarray(Xsub)
        if agg == "median":
            return np.median(Xdense, axis=0)
        if agg == "pct90":
            return np.percentile(Xdense, 90, axis=0)
        raise ValueError(f"Unknown agg='{agg}'. Use 'mean', 'median', or 'pct90'.")

    expr_by_group = np.zeros((n_groups, adata_f.n_vars), dtype=float)
    for i, g in enumerate(group_names):
        idx = (adata_f.obs[condition_col] == g).to_numpy()
        if int(idx.sum()) == 0:
            # keep zeros; avoids divide-by-zero later
            continue
        expr_by_group[i, :] = _agg_vec(adata_f[idx, :].X)

    # ---- tau per gene ----
    x_max = expr_by_group.max(axis=0)
    eps = 1e-12
    x_max_safe = np.where(x_max > 0, x_max, eps)
    x_norm = expr_by_group / x_max_safe  # each in [0,1] if x_max>0
    tau = (np.sum(1.0 - x_norm, axis=0)) / float(n_groups - 1)

    # ---- make it target-specific (so it behaves like your per-target GES) ----
    target_label = str(target_cell_type)
    if target_label not in group_names:
        # This can happen if labels got normalized differently outside
        raise ValueError(f"Target '{target_label}' not found among groups for '{condition_col}': {group_names}")

    target_i = group_names.index(target_label)
    x_target = expr_by_group[target_i, :]
    tau_target_score = tau * (x_target / x_max_safe)  # in [0,1], equals tau if target is the max group

    # ---- the same helper fields you used in GES ----
    filtered_target = adata_f[target_mask, :]
    filtered_other = adata_f[~target_mask, :]

    mean_target_expr = np.asarray(filtered_target.X.mean(axis=0)).ravel()
    mean_other_expr = np.asarray(filtered_other.X.mean(axis=0)).ravel()

    per_target = np.asarray((filtered_target.X > 0).mean(axis=0)).ravel()
    per_other = np.asarray((filtered_other.X > 0).mean(axis=0)).ravel()

    max_group_idx = np.argmax(expr_by_group, axis=0)
    max_group = [group_names[j] for j in max_group_idx]

    results_df = pd.DataFrame(
        {
            "gene": genes_f,
            # keep the same column name your pipeline expects:
            "ges_score": tau_target_score,
            # extra useful columns:
            "tau_target_score": tau_target_score,
            "tau": tau,
            "max_group": max_group,
            "x_target": x_target,
            "x_max": x_max,
            "mean_expression_target": mean_target_expr,
            "mean_expression_other": mean_other_expr,
            "per_expressed_target": per_target,
            "per_expressed_other": per_other,
        }
    )

    print(
        f"finished_tau_calculations in {time.time() - t0:.2f}s | "
        f"kept {adata_f.n_vars}/{adata.n_vars} genes "
        f"(global>= {min_cells_global} cells AND target>= {target_expr_frac:.3f} of target cells) | "
        f"agg={agg} | groups={n_groups}"
    )
    return results_df


# -------------------------------------------------------------------
# Core GES functions
# -------------------------------------------------------------------
def _count_nonzero_per_gene(X) -> np.ndarray:
    """
    Return (#cells expressing gene) for each gene.
    Works for dense or sparse matrices.
    Expression is defined as X > 0.
    """
    if sp.issparse(X):
        # (X > 0) stays sparse boolean; sum over rows gives counts per gene
        return np.asarray((X > 0).sum(axis=0)).ravel()
    else:
        return np.asarray((X > 0).sum(axis=0)).ravel()

def normalize_ges_zscore(df):
    """
    Takes a GES result DataFrame with a column 'ges_score'
    and adds a standardized z-score column.
    """
    mu = df["ges_score"].mean()
    sigma = df["ges_score"].std() or 1e-9  # avoid division by zero
    df["ges_zscore"] = (df["ges_score"] - mu) / sigma
    return df

def calculate_ges(
    adata: sc.AnnData,
    condition_col: str,
    target_cell_type: Any,
    genes: List[str],
    min_cells_global: int = 5,
    target_expr_frac: float = 0.05,
) -> pd.DataFrame:
    """
    Compute GES scores for a specific target cell type, with built-in filtering:

    1) Keep genes expressed in >= min_cells_global cells globally.
    2) From those, keep genes expressed in >= target_expr_frac fraction of TARGET cells.

    Notes:
    - "expressed" means X > 0.
    - Returns a DataFrame only for genes that pass both filters.
    """
    ges_time = time.time()

    total_cells = adata.n_obs
    target_idx = (adata.obs[condition_col] == target_cell_type).to_numpy()
    n_target = int(target_idx.sum())

    if n_target == 0:
        raise ValueError(f"Target '{target_cell_type}' has 0 cells in column '{condition_col}'.")

    target_fraction = float(n_target) / float(total_cells)

    # -------------------------------
    # Filtering
    # -------------------------------
    global_counts = _count_nonzero_per_gene(adata.X)  # length = adata.n_vars
    global_mask = global_counts >= min_cells_global

    target_X = adata.X[target_idx, :]
    target_counts = _count_nonzero_per_gene(target_X)
    target_mask = target_counts >= (target_expr_frac * n_target)

    keep_mask = global_mask & target_mask

    n_keep = int(keep_mask.sum())
    if n_keep == 0:
        # Return empty, but with correct columns
        return pd.DataFrame(
            columns=[
                "gene",
                "ges_score",
                "mean_expression_target",
                "mean_expression_other",
                "per_expressed_target",
                "per_expressed_other",
            ]
        )

    # Subset adata and genes to kept genes
    adata_f = adata[:, keep_mask]
    genes_f = [g for g, keep in zip(genes, keep_mask) if keep]

    # -------------------------------
    # Compute GES on filtered genes
    # -------------------------------
    filtered_data = adata_f[target_idx, :]

    # Mean expression in the target cell type
    mean_target_expr = np.asarray(filtered_data.X.mean(axis=0)).ravel()
    weighted_target_expr = (1.0 - target_fraction) * mean_target_expr

    # Weighted mean expression in other cell types
    cell_types = adata_f.obs[condition_col].unique()
    other_cell_types = [ct for ct in cell_types if ct != target_cell_type]

    weighted_sum = np.zeros(adata_f.n_vars, dtype=float)
    per_exp_other = {}
    
    print("Group sizes:")
    print(adata_f.obs[condition_col].value_counts(dropna=False).head(30))

    for ct in other_cell_types:
        ct_idx = (adata_f.obs[condition_col] == ct).to_numpy()
        ct_fraction = float(ct_idx.sum()) / float(total_cells)

        subset = adata_f[ct_idx]

        if subset.n_obs == 0:
            print(f"⚠️ No cells for {cell_type} – skipping")
            continue
        
        ct_mean_expr = np.asarray(adata_f[ct_idx].X.mean(axis=0)).ravel()
        weighted_sum += ct_fraction * ct_mean_expr

        if len(other_cell_types) > 1:
            per_exp_other[f"per_exp_{normalize_label(ct)}"] = np.asarray(
                (adata_f.X[ct_idx] > 0).mean(axis=0)
            ).ravel()

    epsilon = 1e-10
    ges_scores = weighted_target_expr / (weighted_sum + epsilon)

    filtered_other_data = adata_f[~target_idx, :]
    mean_expression_other = np.asarray(filtered_other_data.X.mean(axis=0)).ravel()

    percent_expressed_target = np.asarray((filtered_data.X > 0).mean(axis=0)).ravel()
    percent_expressed_other = np.asarray((filtered_other_data.X > 0).mean(axis=0)).ravel()

    results_df = pd.DataFrame(
        {
            "gene": genes_f,
            "ges_score": ges_scores,
            "mean_expression_target": mean_target_expr,
            "mean_expression_other": mean_expression_other,
            "per_expressed_target": percent_expressed_target,
            "per_expressed_other": percent_expressed_other,
        }
    )

    if per_exp_other:
        per_expressed_df = pd.DataFrame(per_exp_other)
        per_expressed_df["gene"] = genes_f
        results_df = results_df.merge(per_expressed_df, on="gene", how="left")

    print(
        f"finished_ges_calculations in {time.time() - ges_time:.2f}s | "
        f"kept {n_keep}/{adata.n_vars} genes "
        f"(global>= {min_cells_global} cells AND target>= {target_expr_frac:.3f} of target cells)"
    )
    return results_df


# NOTE: left here for compatibility with your previous workflow.
# It is not called unless permutations=True.
def calculate_ges_with_permutations(
    adata: sc.AnnData,
    condition_col: str,
    target_cell_type: Any,
    genes: List[str],
    expression_threshold: float = 0,
    n_permutations: int = 500,
    min_cells_global = 5
) -> pd.DataFrame:
    """Permutation-based p-values for GES scores."""
    global_counts = np.asarray((adata.X > 0).sum(axis=0)).ravel()
    mask = global_counts >= min_cells_global   # e.g. 20
    adata = adata[:, mask].copy()

    target_mask = adata.obs[condition_col] == target_cell_type
    target_data = adata[target_mask]

    gene_expression_num = np.asarray((target_data.X > 0).sum(axis=0)).ravel()
    expressed_genes_mask = gene_expression_num >= (target_data.n_obs * expression_threshold)
    filtered_adata = adata[:, expressed_genes_mask].copy()

    # Recompute genes list for filtered adata
    filtered_genes = [g for g, keep in zip(genes, expressed_genes_mask) if keep]

    actual_results = calculate_ges(filtered_adata, condition_col, target_cell_type, filtered_genes)
    actual_scores = actual_results["ges_score"].values

    permuted_scores = np.zeros((n_permutations, len(actual_scores)))
    p_values = np.zeros(len(actual_scores))

    for i in range(n_permutations):
        print("permutation", i + 1)
        shuffled = np.random.permutation(filtered_adata.obs[condition_col].to_numpy())
        filtered_adata.obs[condition_col] = shuffled
        permuted = calculate_ges(filtered_adata, condition_col, target_cell_type, filtered_genes)
        permuted_scores[i] = permuted["ges_score"].values
        p_values += (actual_scores < permuted_scores[i]).astype(int)

    actual_results["p_values"] = p_values / float(n_permutations)
    actual_results = actual_results.loc[actual_results["per_expressed_target"] >= expression_threshold, :].copy()
    actual_results["adj_p_val"] = multipletests(actual_results["p_values"], method="fdr_bh")[1]
    return actual_results


# -------------------------------------------------------------------
# Pipeline runner
# -------------------------------------------------------------------

def run_ges_pipeline(config_path: str):
    """
    Run the GES pipeline using a YAML configuration file.
    """
    config = load_config(config_path)

    name_of_run = config["name_of_run"]
    data_path = config["data_path"]
    output_root = Path(config["output_folder"])
    column_conditions = config["column_conditions"]

    normalize_data = bool(config.get("normalize_data", True))
    chemistry = config.get("chemistry", None)
    expression_threshold = float(config.get("expression_threshold", 0.05))
    permutations = bool(config.get("permutations", False))
    n_permutations = int(config.get("n_permutations", 500))
    species = config.get("species", "human")

    # Run folder structure
    date_tag = datetime.datetime.now().strftime("%Y%m%d")
    run_dir = output_root / f"{name_of_run}_{date_tag}"
    metadata_dir = run_dir / "metadata"
    data_dir = run_dir / "data"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    # Copy config file into metadata
    cfg_dest = metadata_dir / Path(config["config_path"]).name
    shutil.copy2(config["config_path"], cfg_dest)

    # Hash the data file
    hash_value = compute_file_hash(data_path, algorithm="sha256")
    with open(metadata_dir / "adata_hash.txt", "w") as hf:
        hf.write(f"file: {data_path}\n")
        hf.write("algorithm: sha256\n")
        hf.write(f"hash: {hash_value}\n")

    # Pretty print
    print("\n==============================")
    print("  Loading GES YAML config")
    print("==============================")
    print(f"• Name of run:          {name_of_run}")
    print(f"• Data path:            {data_path}")
    print(f"• Output root:          {output_root}")
    print(f"• Run directory:        {run_dir}")
    print(f"• Normalize data:       {normalize_data}")
    print(f"• Chemistry:            {chemistry}")
    print(f"• Expression threshold: {expression_threshold}")
    print(f"• Permutations:         {permutations}")
    if permutations:
        print(f"• Num permutations:     {n_permutations}")
    print("==============================\n")

    # Quick obs inspection
    print("Available columns in .obs (fast h5py check):")
    obs_cols = get_obs_columns(data_path)
    print(obs_cols)

    # Load & preprocess full data
    adata = load_and_preprocess_adata(
        data_path=data_path,
        chemistry=chemistry,
        normalize_data=normalize_data,
    )

    # Apply derived columns FIRST (so they exist for GES calculation and for saving edited_adata)
    apply_derived_columns(adata, column_conditions)

    # Save edited adata (with derived columns) into metadata
    print("\nSaving edited adata file with derived columns → metadata/edited_adata.h5ad")
    edited_adata_path = metadata_dir / "edited_adata.h5ad"
    adata.write_h5ad(str(edited_adata_path))

    # Precompute gene names once, aligned to adata.var
    genes = get_gene_names(adata, species=species)
    if len(genes) != adata.n_vars:
        raise ValueError(
            f"Gene list length ({len(genes)}) does not match adata.n_vars ({adata.n_vars}). "
            "Check get_gene_names() and your AnnData var annotation."
        )

    # Process each column spec
    for condition_col, spec in column_conditions.items():
        print(f"\n=== Column: {condition_col} ===")

        # Determine condition_list
        if isinstance(spec, dict):
            condition_list = spec.get("conditions", None)
            if condition_list is None:
                # If user provided only derived.* and forgot conditions
                raise ValueError(f"Column '{condition_col}' dict spec must include 'conditions'")
        else:
            condition_list = spec

        if condition_col not in adata.obs.columns:
            raise ValueError(
                f"Column '{condition_col}' not found in adata.obs after derived processing.\n"
                f"Available columns: {list(adata.obs.columns)}"
            )

        # Normalize labels in the column for robust comparisons
        # (handles bytes and "b'X'" strings)
        adata.obs[condition_col] = adata.obs[condition_col].apply(normalize_label)

        col_series = adata.obs[condition_col]

        # remove NaNs before sorting
        column_values = sorted(pd.unique(col_series.dropna()).tolist())
        print(f"  Available values (sample): {column_values[:25]}{' ...' if len(column_values) > 25 else ''}")
        print(f"  Requested targets:         {condition_list}")

        # Make sure targets are normalized the same way
        norm_targets = [normalize_label(t) for t in condition_list]

        for target_raw, target in zip(condition_list, norm_targets):
            if target not in set(column_values):
                print(f"  ⚠️ '{target_raw}' NOT FOUND in '{condition_col}' → skipping")
                continue

            print(f"\n  → Processing target: {target_raw}")

            target_mask = (adata.obs[condition_col] == target).to_numpy()
            target_data = adata[target_mask]

            if target_data.n_obs == 0:
                print(f"  ⚠️ target '{target_raw}' has zero cells → skipping")
                continue

            # Filter genes by expression threshold in target group
            gene_expression_num = np.asarray((target_data.X > 0).sum(axis=0)).ravel()
            expressed_mask = gene_expression_num >= (target_data.n_obs * expression_threshold)

            if expressed_mask.sum() == 0:
                print("  ⚠️ No genes pass expression threshold → skipping")
                continue

            filtered_adata = adata[:, expressed_mask].copy()
            filtered_genes = [g for g, keep in zip(genes, expressed_mask) if keep]

            # Compute GES
            ges_results = calculate_ges(
                filtered_adata,
                condition_col,
                target,
                filtered_genes
            ).sort_values("ges_score", ascending=False)

            # Compute TAU
            # ges_results = calculate_tau(
            #   filtered_adata,
            #   condition_col,
            #   target,
            #   filtered_genes
            # ).sort_values("tau_target_score", ascending=False)

            target_name = normalize_label(target_raw)
            out_csv = data_dir / f"ges_spec_{condition_col}_{target_name}.csv"
            ges_results.to_csv(out_csv, index=False)
            print(f"  ✔ Saved GES results → {out_csv}")

            if permutations:
                print(f"  🔁 Running {n_permutations} permutations for {target_raw}")
                perm_results = calculate_ges_with_permutations(
                    filtered_adata,
                    condition_col,
                    target,
                    filtered_genes,
                    expression_threshold=expression_threshold,
                    n_permutations=n_permutations,
                )
                out_perm = data_dir / f"ges_spec_{condition_col}_{target_name}_perm.csv"
                perm_results.to_csv(out_perm, index=False)
                print(f"  ✔ Saved permutation results → {out_perm}")

    print("\n🎉 DONE — GES pipeline completed.")
    print(f"All outputs are under:\n  {run_dir}\n")


if __name__ == "__main__":
    # Example:
    #   python ges_score_calculations.py /path/to/config.yaml
    import sys
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python ges_score_calculations.py <config.yaml>")
    run_ges_pipeline(sys.argv[1])
