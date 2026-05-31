import scanpy as sc
import numpy as np
import pandas as pd
import sys
import os
import shutil
import gseapy as gp
import statsmodels.stats.multitest as smm
import matplotlib.pyplot as plt
from pathlib import Path


_PVAL_COLS = ["NOM p-val", "FDR q-val", "FWER p-val", "FDR q-val (BH corrected)"]

def _fmt_pvals(df: pd.DataFrame) -> pd.DataFrame:
    """Format p-value/FDR columns as exact scientific notation strings (e.g. 1.23e-08)."""
    df = df.copy()
    for col in _PVAL_COLS:
        if col in df.columns:
            df[col] = df[col].apply(
                lambda v: f"{float(v):.6e}" if pd.notna(v) and v != "" else v
            )
    return df



def plot_enhanced_gsea(gsea_res, term, cell_type, output_dir):
    """
    Create an enhanced GSEA plot and save it to a file.
    
    Parameters:
    -----------
    gsea_res : GSEApy result object
        The result object from GSEApy prerank analysis
    term : str
        The gene set term to plot
    output_dir : str
        Directory to save the output figure
        
    Returns:
    --------
    fig_path : str
        Path to the saved figure
    """
    # Get the enrichment results for the term
    enrichment_results = gsea_res.results[term]
    
    # Create an enhanced custom GSEA plot
    plt.figure(figsize=(10, 8))
    
    # Plot the enrichment score profile with better styling
    plt.subplot(211)
    plt.plot(enrichment_results.get('RES'), color='forestgreen', linewidth=2)
    plt.axhline(y=0, color='gray', linestyle='--', alpha=0.7)
    plt.ylabel('Enrichment Score', fontsize=12)
    plt.title(f"{term} enrichment in {cell_type} ges scores\nNES={enrichment_results.get('nes'):.3f} FDR={enrichment_results.get('fdr'):.3e}", 
              fontsize=14, fontweight='bold')
    plt.grid(alpha=0.3)
    
    # Add vertical range to match GSEA standard plots
    max_es = max(abs(min(enrichment_results.get('RES'))), abs(max(enrichment_results.get('RES'))))
    plt.ylim(-max_es*1.1, max_es*1.1)
    
    # Plot the hits with a more visible representation
    plt.subplot(212)
    hit_indices = enrichment_results.get('hits')
    # Create a rug plot for the hits
    y = [1] * len(hit_indices)
    plt.plot([0, len(enrichment_results.get('RES'))], [1, 1], color='black', linewidth=0.5, alpha=0.5)
    plt.scatter(hit_indices, y, color='red', s=15, marker='|', alpha=0.8)
    
    # Add a heatmap-style gradient for the rank metric
    ax2 = plt.twinx()
    xs = range(len(gsea_res.ranking))
    ys = [0] * len(xs)
    ax2.scatter(xs, ys, c=gsea_res.ranking, cmap='coolwarm', s=10, marker='_', alpha=0.7)
    ax2.set_ylim(-0.5, 0.5)
    ax2.set_yticks([])
    
    plt.xlabel('Rank in Ordered Dataset', fontsize=12)
    plt.yticks([])
    
    # Add some labels - with adjusted position for "Hits" to avoid the frame
    plt.annotate('Hits', xy=(0.01, 0.95), xycoords='axes fraction', fontsize=10, 
                 fontweight='bold', color='red')
    plt.annotate('Ranked list', xy=(0.01, 0.02), xycoords='axes fraction', fontsize=10, 
                 fontweight='bold', color='blue')
    
    # Improve overall appearance
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.3)
    
    # Save figure with higher quality
    fig = plt.gcf()
    fig_path = os.path.join(output_dir, f"custom_gseaplot_{term}.png")
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return fig_path

def run_gsea(ges_score_path,
             gmt_file,
             column_conditions,
             ges_score_threshold,
             out_folder,
             figs_folder):

    # create resuls folder
    results_folder = out_folder / "GSEA"
    fig_dir = figs_folder / "GSEA"
    
    results_folder.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    # storage for final combined results
    summary_rows = []

  
    # --------------------------
    # Extract GMT name
    # --------------------------
    gmt_name = Path(gmt_file).stem
    print(f"\nUsing GMT set: {gmt_name}\n")

    # --------------------------
    # Run GSEA per column + condition
    # --------------------------
    for column, condition_list in column_conditions.items():

        for condition in condition_list:
            print(f"\nRunning GSEA for {column} → {condition}")

            # locate the matching GES CSV file
            ges_path = Path(ges_score_path) / "data" / f"ges_spec_{column}_{condition}.csv"
            if not ges_path.exists():
                print(f"⚠️ Missing GES file for {column}/{condition}, skipping.")
                continue

            ges_df = pd.read_csv(ges_path)
            ges_df = ges_df[ges_df["ges_score"] > ges_score_threshold]

            if ges_df.empty:
                print(f"⚠️ No genes pass threshold {ges_score_threshold} for {column}/{condition}")
                continue

            # Prepare ranking
            ranking = ges_df["ges_score"]
            ranking.index = ges_df["gene"]

            # Folder for this exact GSEA run
            cond_dir = results_folder / column / str(condition)
            raw_dir = cond_dir / "gsea_raw"
            cond_dir.mkdir(parents=True, exist_ok=True)

            try:
              # Run GSEA
              gsea_res = gp.prerank(
                  rnk=ranking,
                  gene_sets=gmt_file,
                  outdir=str(raw_dir),      # keep GSEApy default raw output
                  min_size=2,
                  max_size=2500,
                  seed=6,
                  verbose=True
              )
            except LookupError:
              print("Probably not enough relevant genes were found for encrichment.")
              gsea_res = None

            if gsea_res:
              # Save main results
              out_csv = cond_dir / "gsea_results.csv"
              _fmt_pvals(gsea_res.res2d).to_csv(out_csv)
              print(f"  ✔ Saved GSEA table to {out_csv}")
  
              # Collect summary row (first term only, like old script)
              top = gsea_res.res2d.iloc[0]
  
              summary_rows.append({
                  "column": column,
                  "condition": condition,
                  "term": top["Term"],
                  "NES": top["NES"],
                  "NOM p-val": top["NOM p-val"],
                  "FDR q-val": top["FDR q-val"],
                  "FWER p-val": top["FWER p-val"],
                  "Tag %": top["Tag %"],
                  "Gene %": top["Gene %"],
                  "Lead_genes": top["Lead_genes"]
              })
            
              # Save plot (first term)
              term = gsea_res.res2d.Term.iloc[0]
              if term in gsea_res.results:
                  plot_out = plot_enhanced_gsea(gsea_res, term, condition, raw_dir)
                  print(f"  ✔ Saved plot to {plot_out}")
                  # fig = gsea_res.plot(term)
                  # plot_out = fig_dir / f"GSEA_{column}_{condition}.png"
                  # fig.savefig(plot_out, dpi=300, bbox_inches="tight")
                  # print(f"  ✔ Saved plot to {plot_out}")
              else:
                  print(f"⚠️ Term '{term}' not found in GSEA result dictionary; skipping plot.")


      # --------------------------------
    # SAVE FINAL SUMMARY
    # --------------------------------
    if summary_rows:
        final_summary = pd.DataFrame(summary_rows)

        # adjust p-values (as old script did)
        reject, pvals_corr = smm.multipletests(
            final_summary["FDR q-val"].astype(float),
            method="fdr_bh"
        )[:2]

        final_summary["FDR q-val (BH corrected)"] = pvals_corr

        final_summary_path = results_folder / "GSEA_final_summary.csv"
        _fmt_pvals(final_summary).to_csv(final_summary_path, index=False)

        print("\n📄 Final summary file created:")
        print(final_summary_path)
        print("\n🎉 GSEA pipeline finished successfully!")
        print(f"All files saved under:\n{results_folder}\n")
        return final_summary_path
    else:
        print("\n⚠️ No GSEA results found — skipping final summary table.")
    
    
