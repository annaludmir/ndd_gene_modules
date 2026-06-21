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



def plot_enhanced_gsea(gsea_res, term, cell_type, output_dir, low_ges_warning=False):
    """
    Create an enhanced GSEA plot and save it to a file.

    low_ges_warning : bool
        When True (ges_score_threshold is None and all leading genes have GES < 1),
        the enrichment score panel is covered with diagonal stripes and a warning label.
    """
    enrichment_results = gsea_res.results[term]

    fig = plt.figure(figsize=(10, 8))

    # ── Panel 1: Enrichment score ──────────────────────────────────────────
    ax1 = fig.add_subplot(211)
    ax1.plot(enrichment_results.get('RES'), color='forestgreen', linewidth=2)
    ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.7)
    ax1.set_ylabel('Enrichment Score', fontsize=12)
    ax1.set_title(
        f"{term} enrichment in {cell_type} ges scores\n"
        f"NES={enrichment_results.get('nes'):.3f}  FDR={enrichment_results.get('fdr'):.3e}",
        fontsize=14, fontweight='bold'
    )
    ax1.grid(alpha=0.3)

    max_es = max(abs(min(enrichment_results.get('RES'))), abs(max(enrichment_results.get('RES'))))
    ax1.set_ylim(-max_es * 1.1, max_es * 1.1)

    # Diagonal-stripe overlay when all leading genes have GES < 1
    if low_ges_warning:
        from matplotlib.patches import Rectangle
        ax1.add_patch(Rectangle(
            (0, 0), 1, 1,
            transform=ax1.transAxes,
            hatch='///', facecolor='#fffacc', edgecolor='#aaaaaa',
            alpha=0.45, zorder=0, linewidth=0,
        ))
        ax1.text(
            0.98, 0.97,
            "All leading genes GES < 1",
            transform=ax1.transAxes, ha='right', va='top',
            fontsize=9, color='#555555', style='italic',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                      edgecolor='#aaaaaa', alpha=0.85),
            zorder=5,
        )

    # ── Panel 2: Hit rug + ranked-list gradient ───────────────────────────
    ax2 = fig.add_subplot(212)
    hit_indices = enrichment_results.get('hits')
    y = [1] * len(hit_indices)
    ax2.plot([0, len(enrichment_results.get('RES'))], [1, 1], color='black', linewidth=0.5, alpha=0.5)
    ax2.scatter(hit_indices, y, color='red', s=15, marker='|', alpha=0.8)

    ax3 = ax2.twinx()
    xs = range(len(gsea_res.ranking))
    ys = [0] * len(xs)
    ax3.scatter(xs, ys, c=gsea_res.ranking, cmap='coolwarm', s=10, marker='_', alpha=0.7)
    ax3.set_ylim(-0.5, 0.5)
    ax3.set_yticks([])

    ax2.set_xlabel('Rank in Ordered Dataset', fontsize=12)
    ax2.set_yticks([])
    ax2.annotate('Hits', xy=(0.01, 0.95), xycoords='axes fraction', fontsize=10,
                 fontweight='bold', color='red')
    ax2.annotate('Ranked list', xy=(0.01, 0.02), xycoords='axes fraction', fontsize=10,
                 fontweight='bold', color='blue')

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.3)

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

            ges_df_full = pd.read_csv(ges_path)
            if ges_score_threshold is not None:
                ges_df = ges_df_full[ges_df_full["ges_score"] > ges_score_threshold]
            else:
                ges_df = ges_df_full.copy()

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

              # When no GES threshold was applied, warn if all leading genes have GES < 1
              low_ges_warning = False
              if ges_score_threshold is None:
                  lead_str = str(top.get("Lead_genes", "") or "")
                  if lead_str and lead_str.lower() != "nan":
                      lead_genes = [g.strip() for g in lead_str.replace(";", ",").split(",") if g.strip()]
                      if lead_genes:
                          lead_scores = ges_df_full[ges_df_full["gene"].isin(lead_genes)]["ges_score"]
                          if len(lead_scores) > 0 and (lead_scores < 1).all():
                              low_ges_warning = True
                              print(f"  ⚠️ All leading genes have GES < 1 — marking plot with stripes")

              # Save plot (first term)
              term = gsea_res.res2d.Term.iloc[0]
              if term in gsea_res.results:
                  plot_out = plot_enhanced_gsea(gsea_res, term, condition, raw_dir,
                                                low_ges_warning=low_ges_warning)
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
    
    
