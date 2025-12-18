import pandas as pd
import scanpy as sc
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

from pathlib import Path
from matplotlib import rc_context

# ----------------------------------------------------------
# 1. Extract leading genes column from enrichment results
# ----------------------------------------------------------
def extract_leading_genes_from_enrichemnt_results(enrichment_results_path):
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


# ----------------------------------------------------------
# 2. Plot multiple images per condition (16 genes per image)
# ----------------------------------------------------------
def plot_tsne_for_gene_list(h5ad_path, enrichment_results_path, ges_scores_folder, output_path, parent_folder):

    print("\n==============================================================")
    print("🔧 STARTING TSNE/UMAP LEADING-GENE PLOTTING")
    print("==============================================================")
    print(f"📁 Output directory: {output_path}")
    print(f"📁 Parent folder:    {parent_folder}")
    print("--------------------------------------------------------------\n")

    enr_df = extract_leading_genes_from_enrichemnt_results(enrichment_results_path)
    
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load dataset
    print(f"📥 Loading dataset:\n   {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path, backed="r")
    print("✔ Dataset loaded")
    print(f"   🔹 Number of cells: {adata.n_obs}")
    print(f"   🔹 Number of genes: {adata.n_vars}\n")

    sc.set_figure_params(dpi=100, color_map="viridis")
    sc.settings.verbosity = 0

    # Iterate over enrichment rows
    for idx, row in enr_df.iterrows():

        print("--------------------------------------------------------------")
        print(f"▶ Processing row {idx+1}/{len(enr_df)}")

        condition = row.get("condition", f"row{idx}")
        column = row.get("column")
        print(f"📌 Condition: {condition}")

        gene_list = row["Lead_genes_list"]
        print(f"🔍 Total leading genes found: {len(gene_list)}")

        if len(gene_list) == 0:
            print(f"⚠️ No leading genes for {parent_folder}/{condition} — skipping.\n")
            continue

        # Filter genes present in dataset
        filtered_genes = [g for g in gene_list if g in adata.var_names]
        print(f"🔎 Genes present in dataset: {len(filtered_genes)}/{len(gene_list)}")

        if len(filtered_genes) == 0:
            print(f"⚠️ No valid genes found in dataset for {condition} — skipping.\n")
            continue

        # Save location
        save_dir = output_path / parent_folder
        save_dir.mkdir(parents=True, exist_ok=True)

        ges_score_path = ges_scores_folder + f'/data/ges_spec_{column}_{condition}.csv'
        print(f"📥 Loading GES score file:\n   {ges_score_path}")
        ges_df = pd.read_csv(ges_score_path)
        
        if "gene" not in ges_df.columns or "ges_score" not in ges_df.columns:
            raise ValueError("GES results file must contain 'gene' and 'ges_score' columns")
        
        # Create a lookup dictionary
        ges_dict = dict(zip(ges_df["gene"], ges_df["ges_zscore"]))
        
        print(f"✔ Loaded GES table ({len(ges_dict)} genes)\n")

        # -----------------------------------------
        # Sort genes by GES score (descending)
        # -----------------------------------------
        # Keep only genes that exist in the GES file (otherwise score is missing)
        scored_genes = [g for g in filtered_genes if g in ges_dict]
        missing_genes = [g for g in filtered_genes if g not in ges_dict]

        # Sort by score descending
        scored_genes = sorted(scored_genes, key=lambda g: ges_dict[g], reverse=True)

        # Final ordered list:
        #   1) scored genes (high→low)
        #   2) then any genes missing from GES (optional; keep them at the end)
        filtered_genes = scored_genes + missing_genes

        print("✅ Sorted genes by GES score (descending)")
        if missing_genes:
            print(f"⚠️ {len(missing_genes)} genes missing from GES file; appended at the end")

        # ---- Split into batches of 16 genes ----
        batch_size = 12
        total_batches = int(np.ceil(len(filtered_genes) / batch_size))

        print(f"📦 Splitting into {total_batches} batches (each up to {batch_size} genes)\n")

        for batch_index in range(total_batches):

            batch_genes = filtered_genes[batch_index * batch_size : (batch_index + 1) * batch_size]

            print(f"🖼 Creating batch {batch_index+1}/{total_batches} → {len(batch_genes)} genes")

            save_file = save_dir / f"{column}_{condition}_part{batch_index+1}.png"
            print(f"📂 Saving to: {save_file}")

            # Compute grid
            n = len(batch_genes)
            ncols = 4
            nrows = int(np.ceil(n / ncols))

            fig, axes = plt.subplots(
                nrows=nrows,
                ncols=ncols,
                figsize=(4 * ncols, 4 * nrows),
                squeeze=False
            )

            for ax, gene in zip(axes.flat, batch_genes):
                sc.pl.umap(
                    adata,
                    color=gene,
                    ax=ax,
                    show=False,
                    s=8,
                    frameon=False,
                    vmax="p99"
                )
                ges_val = ges_dict.get(gene)
                ges_str = f"{ges_val:.2f}" if isinstance(ges_val, (int, float)) else "NA"
                ax.set_title(f"{gene}\nGES={ges_str}")

            # Hide unused axes
            for ax in axes.flat[n:]:
                ax.axis("off")

            plt.tight_layout()
            plt.savefig(save_file, dpi=200)
            plt.close()

            print(f"✔ Saved batch {batch_index+1}/{total_batches}\n")

    print("==============================================================")
    print("🎉 FINISHED: All TSNE/UMAP leading-gene plots generated!")
    print("==============================================================\n")
