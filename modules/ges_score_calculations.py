import marimo

__generated_with = "0.17.8"
app = marimo.App()


@app.cell
def _():
    """
    Refactored GES (Generalized Expression Specificity) pipeline.

    Key points:
    - No work is done at import time (no data loaded, no loops run).
    - Main entry point for running the whole pipeline: `run_ges_pipeline(...)`.
    - Core logic of `calculate_ges` and `calculate_ges_with_permutations` is preserved.
    - Data paths, output folder, column names, and conditions are configurable via arguments.
    """

    import time
    import yaml
    import os
    from typing import Iterable, List, Optional

    import numpy as np
    import pandas as pd
    import scanpy as sc
    from statsmodels.stats.multitest import multipletests


    def load_config(config_path: str) -> dict:
        """
        Load a YAML configuration file for the GES pipeline.

        The YAML must define:
            - data_type
            - data_path
            - output_folder
            - column_conditions  (mapping: column -> list of conditions)
        """
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)

        # Basic validation
        required = ["data_type", "data_path", "output_folder", "column_conditions"]
        for key in required:
            if key not in config:
                raise ValueError(f"Missing required parameter in config: {key}")

        if not os.path.exists(config["data_path"]):
            raise FileNotFoundError(f"Data file not found: {config['data_path']}")

        os.makedirs(config["output_folder"], exist_ok=True)

        return config


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


    def load_and_preprocess_adata(data_type, path_all, path_cortex):
        """
        Load and preprocess AnnData according to the chosen data_type.

        Parameters
        ----------
        data_type : {'data_all', 'cortex'}
        path_all : path to full dataset .h5ad
        path_cortex : path to cortex dataset .h5ad

        Returns
        -------
        adata : AnnData
            Preprocessed AnnData object.
        """
        # TODO: make it more generic, there is no real difference between cortex and all data, at least not supposed to be
        print("Uploading data...")
        start_time = time.time()

        if data_type == "cortex":
            adata = sc.read_h5ad(path_cortex)
            print("cortex data uploaded")
            # extract only v3 chemistry
            adata = adata[adata.obs["Chemistry"] == "v3", :]
            print("data copied to v3 array")
        else:
            # all data
            adata = sc.read_h5ad(path_all)
            # Originally there was an option to remove 'Head' and 'Brain'
            # # adata = adata[~adata.obs.Region.isin([\"b'Head'\",\"b'Brain'\"])]
            adata = adata[adata.obs["Chemistry"] == "v3", :]
            print("all-data uploaded")

        print(f"data loaded in {time.time() - start_time:.2f} seconds.")

        # data preprocessing
        prepro_time = time.time()
        # filter cells with less than 200 genes expressed and genes that express in less than 3 cells
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

        print(f"preprocessing done in {time.time() - prepro_time:.2f} seconds")
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

        YAML structure supports two column types:
    
        1. Normal columns:
            column_name:
              conditions: ["RG", "IPC", ...]

        2. Derived boolean columns:
            column_name:
              derived:
                boolean_expr: "CellClass == 'RG' or CellClass == 'IPC'"
              conditions:
                - True
        """

        # -------------------------------------------------------
        # Helper for boolean derived columns
        # -------------------------------------------------------
        def compute_boolean_column(adata, column_name: str, expr: str):
            """
            Create adata.obs[column_name] as a boolean mask evaluated from expr.
            Example expr:
                "(CellClass == 'RG') or (CellClass == 'IPC')"
            """
            print(f"    Evaluating boolean expression for '{column_name}': {expr}")

            mask = adata.obs.eval(expr)

            # Ensure boolean dtype, no NaN
            mask = mask.fillna(False).astype(bool)

            adata.obs[column_name] = mask
            print(f"    → Derived boolean column '{column_name}' created ({mask.sum()} True cells)")

        # -------------------------------------------------------
        # Load config
        # -------------------------------------------------------
        config = load_config(config_path)

        data_type = config["data_type"]
        data_path = config["data_path"]
        output_folder = config["output_folder"]
        column_conditions = config["column_conditions"]

        chemistry = config.get("chemistry", "v3")
        expression_threshold = config.get("expression_threshold", 0.05)
        permutations = config.get("permutations", False)
        n_permutations = config.get("n_permutations", 500)

        print("\n==============================")
        print("  Loading GES YAML config")
        print("==============================")
        print(f"• Data type:            {data_type}")
        print(f"• Data path:            {data_path}")
        print(f"• Output folder:        {output_folder}")
        print(f"• Chemistry:            {chemistry}")
        print(f"• Expression threshold: {expression_threshold}")
        print(f"• Permutations:         {permutations}")
        print("==============================\n")

        os.makedirs(output_folder, exist_ok=True)

        # -------------------------------------------------------
        # Load ANNData
        # -------------------------------------------------------
        adata = load_and_preprocess_adata(
            data_type=data_type,
            path_all=data_path,
            path_cortex=data_path,
        )

        # -------------------------------------------------------
        # Process each column (normal or derived)
        # -------------------------------------------------------
        for condition_col, spec in column_conditions.items():
            print(f"\n=== Column: {condition_col} ===")

            # Handle dictionary-based spec
            if isinstance(spec, dict):

                # -------- BOOLEAN DERIVED COLUMN -------
                if "derived" in spec:

                    if "boolean_expr" not in spec["derived"]:
                        raise ValueError(
                            f"Derived column '{condition_col}' requires 'boolean_expr' in YAML."
                        )

                    boolean_expr = spec["derived"]["boolean_expr"]

                    print(f"  Detected BOOLEAN derived column: {condition_col}")
                    compute_boolean_column(adata, condition_col, boolean_expr)

                    condition_list = [True]   # always run on True only

                # -------- NORMAL COLUMN -------
                else:
                    condition_list = spec.get("conditions", [])
                    if not condition_list:
                        raise ValueError(
                            f"Column '{condition_col}' has no 'conditions' list in YAML."
                        )

            else:
                # Legacy support: spec is just a list
                condition_list = spec

            # Extract present values
            column_values = adata.obs[condition_col].unique().tolist()
            print(f"  Available values:   {column_values}")
            print(f"  Requested targets:  {condition_list}")

            # -------------------------------------------------------
            # Run GES for each requested condition
            # -------------------------------------------------------
            for target in condition_list:

                # Boolean columns: target = True/False
                if target not in column_values:
                    print(f"  ⚠️ '{target}' NOT FOUND in column '{condition_col}' → skipping")
                    continue

                print(f"\n  → Processing target: {target}")

                target_mask = adata.obs[condition_col] == target
                target_data = adata[target_mask]

                # Empty population? Skip
                if target_data.n_obs == 0:
                    print(f"  ⚠️ target '{target}' has zero cells → skipping")
                    continue

                # Filter genes by minimal expression
                gene_expression_num = (target_data.X > 0).sum(axis=0).A1
                expressed_mask = gene_expression_num >= (
                    target_data.shape[0] * expression_threshold
                )

                # Avoid zero genes issue
                if expressed_mask.sum() == 0:
                    print(f"  ⚠️ No genes pass expression threshold → skipping")
                    continue

                filtered_adata = adata[:, expressed_mask]

                # Run GES
                ges_results = calculate_ges(
                    filtered_adata,
                    condition_col,
                    target,
                    chemistry,
                    data_type,
                )

                ges_results = ges_results.sort_values("ges_score", ascending=False)

                target_name = normalize_label(target)

                out_path = os.path.join(
                    output_folder,
                    f"ges_spec_{data_type}_{chemistry}_{condition_col}_{target_name}.csv",
                )

                ges_results.to_csv(out_path)
                print(f"  ✔ Saved GES results → {out_path}")

                # Optional permutations
                if permutations:
                    print(f"  🔁 Running {n_permutations} permutations for {target_name}")

                    perm_results = calculate_ges_with_permutations(
                        filtered_adata,
                        condition_col,
                        target,
                        chemistry,
                        data_type,
                        expression_threshold=expression_threshold,
                        n_permutations=n_permutations,
                    )

                    perm_path = os.path.join(
                        output_folder,
                        f"ges_spec_{data_type}_{chemistry}_{condition_col}_{target_name}_perm.csv",
                    )

                    perm_results.to_csv(perm_path)
                    print(f"  ✔ Saved permutation results → {perm_path}")

        print("\n🎉 DONE — GES pipeline completed.\n")


    # -------------------------------------------------------------------
    # Optional: allow script execution with defaults
    # -------------------------------------------------------------------

    def _main_():
        # Runs with default settings from constants above.
        # No CLI args are used; this is just a convenience.
        run_ges_pipeline()
    return


if __name__ == "__main__":
    app.run()
