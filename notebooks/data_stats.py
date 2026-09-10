import marimo

__generated_with = "0.17.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return


@app.cell
def _():
    import mygene
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    return mygene, np, pd, plt, sc, sns


@app.cell
def _(sc):
    all_layers=sc.read_h5ad("/miridan-data/annaludmir/ndd_gene_modules/data/full_adata_with_umap_without_week_5.h5ad")
    return (all_layers,)


@app.cell
def _(all_layers):
    all_layers.obs
    return


@app.cell
def _(all_layers, np, pd):
    # Extract metadata DataFrame
    obs_general = (
        all_layers.obs.copy()
        if hasattr(all_layers, "obs")
        else all_layers.copy()
    )
    df_general = obs_general.copy()

    # Dynamically identify regional metadata column
    region_col_key = next(
        (
            c
            for c in ["region", "Region", "top_region", "tissue", "brain_region"]
            if c in df_general.columns
        ),
        "region",
    )

    # ---------------------------------------------------------
    # 1. Map Detailed Regions to Top-Level General Regions
    # ---------------------------------------------------------
    forebrain_set = {"Diencephalon", "Forebrain", "Telencephalon"}
    hindbrain_set = {"Hindbrain", "Cerebellum", "Pons", "Medulla"}
    midbrain_set = {"Midbrain"}


    def assign_general_region(region_str):
        if region_str in forebrain_set:
            return "Forebrain_General"
        elif region_str in hindbrain_set:
            return "Hindbrain_General"
        elif region_str in midbrain_set:
            return "Midbrain_General"
        return np.nan


    df_general["RegionsGeneral"] = (
        df_general[region_col_key].astype(str).map(assign_general_region)
    )

    # Drop rows that don't belong to any general region (e.g. 'Brain', 'Head', NaNs)
    df_general = df_general.dropna(subset=["RegionsGeneral"]).copy()

    # ---------------------------------------------------------
    # 2. Define Cell State Categories
    # ---------------------------------------------------------
    is_proliferating = df_general["CellClass"].isin(
        ["Radial glia", "Neuronal IPC", "Glioblast"]
    )
    is_differentiating = df_general["CellClass"].isin(["Neuroblast", "Neuron"])
    is_cycling_cell = df_general["cell_cycle_score"] > 0.004

    conditions_list = [
        is_proliferating & is_cycling_cell,
        is_differentiating & is_cycling_cell,
        is_proliferating & (~is_cycling_cell),
        is_differentiating & (~is_cycling_cell),
    ]

    category_labels = [
        "Proliferating_Cycling",
        "Differentiating_Cycling",
        "Proliferating_NonCycling",
        "Differentiating_NonCycling",
    ]

    df_general["State_Category"] = np.select(
        conditions_list, category_labels, default="Other"
    )

    # ---------------------------------------------------------
    # 3. Aggregate Counts & Calculate Percentages
    # ---------------------------------------------------------
    general_counts = (
        df_general.groupby(["Age", "RegionsGeneral", "State_Category"])
        .size()
        .unstack(fill_value=0)
    )

    target_columns = [
        "Proliferating_Cycling",
        "Differentiating_Cycling",
        "Proliferating_NonCycling",
        "Differentiating_NonCycling",
    ]
    general_counts = general_counts.reindex(columns=target_columns, fill_value=0)

    # Calculate percentages normalized per (Age, RegionsGeneral) row
    general_totals = general_counts.sum(axis=1)
    general_percentages = general_counts.div(general_totals, axis=0) * 100

    summary_general_df = general_percentages.reset_index().copy()

    # Clean numeric values
    summary_general_df["Proliferating_Cycling"] = pd.to_numeric(
        summary_general_df["Proliferating_Cycling"], errors="coerce"
    )
    summary_general_df["Age"] = pd.to_numeric(
        summary_general_df["Age"], errors="coerce"
    )
    return (summary_general_df,)


@app.cell
def _(pd, plt, sns, summary_general_df):
    # ---------------------------------------------------------
    # 4. Apply Custom Ordered Categorical (Top Regions Only)
    # ---------------------------------------------------------
    general_desired_order = [
        "Forebrain_General",
        "Midbrain_General",
        "Hindbrain_General",
    ]

    summary_general_df["RegionsGeneral"] = pd.Categorical(
        summary_general_df["RegionsGeneral"].astype(str).str.strip(),
        categories=general_desired_order,
        ordered=True,
    )

    # Clean missing entries and sort by Age
    plot_general_df = (
        summary_general_df.dropna(
            subset=["Age", "RegionsGeneral", "Proliferating_Cycling"]
        )
        .sort_values(["Age", "RegionsGeneral"])
        .copy()
    )

    # ---------------------------------------------------------
    # 5. Plotting (Fixed for Marimo)
    # ---------------------------------------------------------
    fig_gen, ax_gen = plt.subplots(figsize=(9, 5.5), dpi=300)

    sns.lineplot(
        data=plot_general_df,
        x="Age",
        y="Proliferating_Cycling",
        hue="RegionsGeneral",
        style="RegionsGeneral",
        marker="o",
        markersize=8,
        linewidth=2.5,
        dashes=False,
        ax=ax_gen,
    )

    ax_gen.set_title(
        "Proliferating Cycling Cell Percentage Across Ages & Top Regions",
        fontsize=12,
        fontweight="bold",
        pad=12,
    )
    ax_gen.set_xlabel("Age", fontsize=11, fontweight="bold")
    ax_gen.set_ylabel(
        "Proliferating Cycling Cells (%)", fontsize=11, fontweight="bold"
    )
    ax_gen.legend(
        title="General Region",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=True,
    )

    plt.tight_layout()

    # In Marimo, evaluating the figure object as the final line renders it automatically
    fig_gen
    return


@app.cell
def _(sns):
    CANONICAL_REGIONS = [
        "Forebrain",
        "Telencephalon",
        "Diencephalon",
        "Midbrain",
        "Hindbrain",
        "Cerebellum",
        "Pons",
        "Medulla",
    ]

    sns.set_theme(style="whitegrid")
    return


@app.cell
def _(all_layers):
    # Replace 'all_layers' with your AnnData variable if different
    adata_dataset = all_layers

    # Extract raw IDs from var_names
    raw_var_ids = adata_dataset.var_names.tolist()

    # CLEANING: Strip whitespace and remove version suffixes (e.g., '.2')
    cleaned_ensembl_ids = [str(gid).strip().split(".")[0] for gid in raw_var_ids]

    # Map cleaned IDs to original var_names to maintain an exact index reference
    id_clean_to_raw_map = dict(zip(cleaned_ensembl_ids, raw_var_ids))
    return adata_dataset, cleaned_ensembl_ids, id_clean_to_raw_map


@app.cell
def _(cleaned_ensembl_ids, id_clean_to_raw_map, mygene):
    mg_client = mygene.MyGeneInfo()

    # Query mygene using the cleaned ENSEMBL IDs
    query_response = mg_client.querymany(
        cleaned_ensembl_ids,
        scopes="ensembl.gene",
        fields="symbol",
        species="human",
        as_dataframe=True,
        silent=True,
    )

    # Build a lookup table from Symbol -> Original var_name
    symbol_to_varname_dict = {}

    if query_response is not None and not query_response.empty:
        # Ensure both 'query' and 'symbol' columns exist in response
        if "symbol" in query_response.columns:
            # Drop rows where symbol is NaN/missing
            valid_mappings = query_response.dropna(subset=["symbol"]).copy()

            # Iterate through rows cleanly using itertuples to handle duplicate query hits
            for row in valid_mappings.itertuples():
                # Get the query ENSEMBL ID (handles case where query is column vs index)
                cleaned_id = getattr(row, "query", getattr(row, "Index", None))
                sym = getattr(row, "symbol", None)

                if cleaned_id and sym and isinstance(sym, str):
                    raw_id = id_clean_to_raw_map.get(str(cleaned_id))
                    if raw_id:
                        # Store mapping (e.g. 'FOXG1': 'ENSG00000176887.12')
                        symbol_to_varname_dict[sym.strip().upper()] = raw_id

    print(
        f"✅ Successfully mapped {len(symbol_to_varname_dict)} unique gene symbols to dataset var_names."
    )
    return (symbol_to_varname_dict,)


app._unparsable_cell(
    r"""
    symbol_to_varname_dict[]
    """,
    name="_"
)


@app.cell
def _(adata_dataset, symbol_to_varname_dict):
    # ------------------------------------------------------------------------------
    # Cell 4 (FIXED): Process Target Gene List & Extract Expression Matrix
    # ------------------------------------------------------------------------------
    input_gene_symbols = ["FOXG1", "TCF4", "ZEB2", "NDE1"]

    # 1. Look up target symbols using uppercase keys
    target_var_ids = []
    found_symbols = []
    missing_symbols = []

    for g in input_gene_symbols:
        sym_clean = g.strip().upper()
        if sym_clean in symbol_to_varname_dict:
            target_var_ids.append(symbol_to_varname_dict[sym_clean])
            found_symbols.append(g)
        else:
            missing_symbols.append(g)

    print(f"✅ Found target genes in mapping: {found_symbols}")
    if missing_symbols:
        print(f"⚠️ Missing from mapping: {missing_symbols}")

    if not target_var_ids:
        raise ValueError("None of the target genes were mapped to dataset var_names.")

    # 2. Filter metadata for Chemistry == 'v3'
    v3_mask = adata_dataset.obs["Chemistry"].astype(str).str.lower() == "v3"
    v3_obs_df = adata_dataset.obs[v3_mask].copy()

    # Dynamically find region column name
    region_column_name = next(
        (c for c in ["region", "Region", "top_region", "tissue", "brain_region"] if c in v3_obs_df.columns),
        "region"
    )

    # 3. Subset expression matrix for Chemistry v3 & target genes
    v3_adata_subset = adata_dataset[v3_obs_df.index, target_var_ids]
    expr_matrix = v3_adata_subset.X
    if hasattr(expr_matrix, "toarray"):
        expr_matrix = expr_matrix.toarray()

    # Calculate union detection (expresses ANY target gene > 0)
    v3_obs_df["expresses_gene_list_union"] = (expr_matrix > 0).any(axis=1)
    return


app._unparsable_cell(
    r"""
    summary_percentage_df = (
            v3_obs_df.groupby([\"Age\", region_column_name])[
                \"expresses_gene_list_union\"
            ]
            .mean()
            .reset_index()
        )

        summary_percentage_df[\"expresses_gene_list_union\"] *= 100
        summary_percentage_df.rename(
            columns={\"expresses_gene_list_union\": \"Pct_Expressing\"}, inplace=True
        )

        # Clean numeric Age and filter out excluded regions ('Brain', 'Head')
        summary_percentage_df[\"Age\"] = pd.to_numeric(
            summary_percentage_df[\"Age\"], errors=\"coerce\"
        )
        excluded_regions = [\"brain\", \"head\"]
        summary_percentage_df = summary_percentage_df[
            ~summary_percentage_df[region_column_name]
            .astype(str)
            .str.strip()
            .str.lower()
            .isin(excluded_regions)
        ].copy()

        # Apply ordered Categorical dtype to enforce exact legend order
        summary_percentage_df[region_column_name] = pd.Categorical(
            summary_percentage_df[region_column_name].astype(str).str.strip(),
            categories=CANONICAL_REGIONS,
            ordered=True,
        )

        # Clean NaN values and sort by Age
        summary_percentage_df = summary_percentage_df.dropna(
            subset=[\"Age\", region_column_name, \"Pct_Expressing\"]
        ).sort_values(\"Age\")
    """,
    name="_"
)


app._unparsable_cell(
    r"""
    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=300)

        sns.lineplot(
            data=summary_percentage_df,
            x=\"Age\",
            y=\"Pct_Expressing\",
            hue=region_column_name,
            marker=\"o\",
            markersize=8,
            linewidth=2.5,
            ax=ax,
        )

        found_symbols_str = \", \".join(
            [
                g
                for g in input_gene_symbols
                if g.upper() in symbol_to_varname_dict
            ]
        )

        ax.set_title(
            f\"Fraction of Cells Expressing Gene List Union Across Ages & Regions\n({found_symbols_str})\",
            fontsize=12,
            fontweight=\"bold\",
            pad=12,
        )
        ax.set_xlabel(\"Age\", fontsize=11, fontweight=\"bold\")
        ax.set_ylabel(
            \"Cells Expressing Gene Union (%)\", fontsize=11, fontweight=\"bold\"
        )
        ax.legend(
            title=\"Region\",
            bbox_to_anchor=(1.02, 1),
            loc=\"upper left\",
            frameon=True,
        )

        plt.tight_layout()
        plt.show()
    """,
    name="_"
)


if __name__ == "__main__":
    app.run()
