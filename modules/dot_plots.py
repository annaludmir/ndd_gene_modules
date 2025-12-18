import pandas as pd
import scanpy as sc
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

from pathlib import Path
from matplotlib import rc_context

GROUP_ORDERS = {
    "CellClass": ["Glioblast", "Radial glia", "Neuronal IPC", "Neuron", "Neuroblast"],
    "CellCyclePhase": ["G1", "S", "G2M", "PostM", "Non-cycling"],
    # add more groupby columns here if you want
}

def apply_group_order(adata, groupby: str, desired_order: list[str]) -> None:
    """
    Force a specific y-axis order for scanpy dotplot by making adata.obs[groupby]
    a pandas Categorical with an explicit category order.

    - Keeps only categories that actually exist in the data
    - Appends any missing/unknown categories at the end (so nothing disappears)
    """
    s = adata.obs[groupby].astype(str)

    present = list(pd.unique(s))
    desired_present = [x for x in desired_order if x in set(present)]
    leftovers = [x for x in present if x not in set(desired_present)]

    final_order = desired_present + leftovers

    adata.obs[groupby] = pd.Categorical(s, categories=final_order, ordered=True)


# ----------------------------------------------------------
# 1. Extract leading genes column from enrichment results
# ----------------------------------------------------------
def extract_leading_genes_from_enrichment_results(enrichment_results_path):
    print("\n==============================")
    print(f"📄 Loading enrichment results from:\n   {enrichment_results_path}")
    print("==============================")

    df = pd.read_csv(enrichment_results_path)

    if "Lead_genes" not in df.columns:
        raise ValueError("❌ Column 'Lead_genes' not found in enrichment results")

    print(f"✔ Loaded {len(df)} enrichment rows")
    print("🔍 Parsing 'Lead_genes' column... (splitting by ';')")

    df["Lead_genes_list"] = df["Lead_genes"].apply(
        lambda x: [g.strip() for g in str(x).split(";") if g.strip()]
    )

    print("✔ Extracted leading gene lists\n")
    return df


def plot_dots_for_gene_list(
    h5ad_path,
    enrichment_results_path,
    output_path,
    parent_folder,
):
    """
    Creates:
      → ONE dot plot per column
    showing all leading genes across all conditions of that column.

    Saves to:
      {output_path}/{parent_folder}/{column}_leading_genes.png
    """

    print("\n==============================================================")
    print("🔧 STARTING DOT PLOT FOR LEADING GENES (per column)")
    print("==============================================================")
    print(f"📁 Output directory: {output_path}")
    print(f"📁 Parent folder:    {parent_folder}")
    print("--------------------------------------------------------------\n")

    # ----------------------------------------------------------
    # 1. Load enrichment results + parse leading genes
    # ----------------------------------------------------------
    enr_df = extract_leading_genes_from_enrichment_results(enrichment_results_path)

    if "column" not in enr_df.columns:
        raise ValueError("❌ 'column' column is missing in enrichment_results CSV")

    output_path = Path(output_path)
    save_dir = output_path / parent_folder
    save_dir.mkdir(parents=True, exist_ok=True)

    # ----------------------------------------------------------
    # 2. Load AnnData (full, not backed, for dotplot)
    # ----------------------------------------------------------
    print(f"📥 Loading dataset:\n   {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print("✔ Dataset loaded")
    print(f"   🔹 Number of cells: {adata.n_obs}")
    print(f"   🔹 Number of genes: {adata.n_vars}\n")

    sc.set_figure_params(dpi=100, color_map="viridis")
    sc.settings.verbosity = 0

    # ----------------------------------------------------------
    # 3. For each column: gather all leading genes and dot-plot
    # ----------------------------------------------------------
    unique_columns = enr_df["column"].unique()
    print(f"📊 Found {len(unique_columns)} columns in enrichment results: {unique_columns}\n")

        # ----------------------------------------------------------
    # 3. For each column: gather all leading genes and dot-plots
    # ----------------------------------------------------------
    unique_columns = enr_df["column"].unique()
    print(f"📊 Found {len(unique_columns)} columns in enrichment results: {unique_columns}\n")

    for column in unique_columns:
        print("--------------------------------------------------------------")
        print(f"▶ Processing column: {column}")

        sub = enr_df[enr_df["column"] == column]

        # Collect union of all leading genes for this column
        all_genes = set()
        for genes in sub["Lead_genes_list"]:
            all_genes.update(genes)

        all_genes = sorted(all_genes)
        print(f"🔍 Total unique leading genes for column '{column}': {len(all_genes)}")

        # Filter to genes present in adata
        genes_in_adata = [g for g in all_genes if g in adata.var_names]
        print(f"🔎 Genes present in dataset: {len(genes_in_adata)}/{len(all_genes)}")

        if len(genes_in_adata) == 0:
            print(f"⚠️ No valid genes found in adata for column '{column}' — skipping.\n")
            continue

        # Check that the column exists in adata.obs for grouping
        if column not in adata.obs.columns:
            print(f"⚠️ Column '{column}' not found in adata.obs — skipping.\n")
            continue

        # ------------------------------------------------------
        # 3a. Compute fraction of expressing cells per gene
        #     (used to sort genes by "dot size")
        # ------------------------------------------------------
        adata_sub = adata[:, genes_in_adata]

        if sp.issparse(adata_sub.X):
            frac_expr = (adata_sub.X > 0).mean(axis=0)
            frac_expr = np.asarray(frac_expr).ravel()
        else:
            frac_expr = (adata_sub.X > 0).mean(axis=0)
            frac_expr = np.asarray(frac_expr).ravel()

        # Sort genes by fraction expressed, descending (biggest dots first)
        order_idx = np.argsort(-frac_expr)
        genes_sorted = [genes_in_adata[i] for i in order_idx]
        print("📐 Sorted genes by fraction of expressing cells (descending)")

        # ------------------------------------------------------
        # 3b. Split into chunks of N genes per dot plot
        # ------------------------------------------------------
        chunk_size = 30  # you can change to 20/30 if you like
        num_chunks = int(np.ceil(len(genes_sorted) / chunk_size))
        print(f"📦 Splitting into {num_chunks} chunks of up to {chunk_size} genes\n")

        for i in range(num_chunks):
            chunk_genes = genes_sorted[i * chunk_size : (i + 1) * chunk_size]
            if not chunk_genes:
                continue

            out_png = save_dir / f"{column}_leading_genes_part{i+1}.png"
            print(f"🖼 Creating dot plot {i+1}/{num_chunks} with {len(chunk_genes)} genes")
            print(f"📂 Saving to: {out_png}")

            # Width scales with number of genes so labels are readable
            fig_width = max(6, len(chunk_genes) * 0.35)
            fig_height = 5

            with rc_context({"figure.figsize": (fig_width, fig_height)}):
              # Enforce y-axis order if we have a desired order for this column
              if column in GROUP_ORDERS:
                apply_group_order(adata, column, GROUP_ORDERS[column])
                dp = sc.pl.dotplot(
                    adata,
                    var_names=chunk_genes,
                    groupby=column,
                    standard_scale="var",
                    show=False,
                    return_fig=True,   # <-- IMPORTANT!
                )
              
                dp.make_figure()
                fig = dp.fig
                
                # Adjust margins manually instead of tight_layout()
                fig.subplots_adjust(
                    left=0.15,
                    right=0.95,
                    bottom=0.25,
                    top=0.95
                )
                
                fig.savefig(out_png, dpi=200, bbox_inches="tight")
                plt.close(fig)

                print(f"✔ Saved {column} chunk {i+1}/{num_chunks} → {out_png}\n")

    print("==============================================================")
    print("🎉 FINISHED: Dot plots for leading genes generated (per column)!")
    print("==============================================================\n")

