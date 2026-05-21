import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

def plot_bar_chart(ges_results_df, output_file, gene_list_name, comparison):
    # ---- Find q-value column robustly ----
    q_candidates = [
        "FDR q-val (BH corrected)",
        "FDR q-val (BH corrected) ",
        "FDR q-val",
        "FDR q-val (BH)",
        "FDR q-val (BH corrected)\n",  # just in case of weird headers
    ]
    q_col = next((c for c in q_candidates if c in ges_results_df.columns), None)

    # Ensure NES numeric (safe if strings exist)
    nes_vals = pd.to_numeric(ges_results_df["NES"], errors="coerce").to_numpy()

    # Bar colors (blue for >=0, red for <0)
    colors = ['#1f77b4' if v >= 0 else '#d62728' for v in nes_vals]
    x_pos = np.arange(len(ges_results_df)) * 0.9

    # Decide legend
    add_legend = float(np.nanmax(nes_vals)) >= 2 if np.isfinite(np.nanmax(nes_vals)) else False
    fig_width = 4 if not add_legend else 6
    fig, ax = plt.subplots(figsize=(fig_width, 4))

    bars = ax.bar(
        x_pos,
        nes_vals,
        color=colors,
        edgecolor="black",
        linewidth=0.7,
        width=0.7,
        alpha=0.8,
    )

    # X labels
    ax.set_xticks(x_pos)
    ax.set_xticklabels(ges_results_df["condition"].astype(str), fontsize=10, rotation=80)

    # Y limits
    ymin = float(np.nanmin(nes_vals)) - 0.2 if np.isfinite(np.nanmin(nes_vals)) else -1
    ymax = float(np.nanmax(nes_vals)) + 0.2 if np.isfinite(np.nanmax(nes_vals)) else 1
    ax.set_ylim(ymin, ymax)

    # Baseline
    ax.axhline(y=0, color="black", linestyle="-", linewidth=1)

    # Significance lines
    if add_legend:
        ax.axhline(y=1.5, color="#555555", linestyle="--", linewidth=1.5, label="Significant enrichment")
        ax.axhline(y=2, color="darkorange", linestyle="-", linewidth=2, label="Very Significant enrichment")
        ax.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize=8)
    else:
        ax.axhline(y=1.5, color="#555555", linewidth=1.5, linestyle=(0, (5, 3)))

    ax.set_title(f"GSEA enrichment results for {gene_list_name} - {comparison}", pad=14, fontsize=11)
    ax.set_ylabel("NES Score", fontsize=10)

    max_q_text_y = -np.inf

    # ---- Add q-value annotations per bar (black if significant, gray otherwise) ----
    if q_col is not None:
        y_range = max(1e-9, (ymax - ymin))
        offset = 0.03 * y_range

        q_series = pd.to_numeric(ges_results_df[q_col], errors="coerce").to_numpy()

        for i, rect in enumerate(bars):
          nes = float(rect.get_height())
          qv = q_series[i] if i < len(q_series) else np.nan
      
          if np.isfinite(qv):
              q_txt = "q≈0" if qv == 0.0 else f"q={qv:.2e}"
              sig = (qv < 0.05)
          else:
              q_txt = "q=NA"
              sig = False
      
          txt_color = "black" if sig else "#7a7a7a"
      
          x = rect.get_x() + rect.get_width() / 2
          if nes >= 0:
              y = nes + offset
              va = "bottom"
          else:
              y = 2 * offset
              va = "top"
      
          # Track highest label
          if y > max_q_text_y:
              max_q_text_y = y
      
          ax.text(
              x, y, q_txt,
              ha="center", va=va,
              fontsize=8,
              color=txt_color,
              fontweight="bold" if sig else "normal",
              clip_on=False
          )

    else:
        print("⚠️ Could not find q-value column. Skipping q-value annotations.")

    # Ensure q-labels never touch the top border
    if max_q_text_y > ymax:
        padding = 0.05 * (ymax - ymin)
        ax.set_ylim(ymin, max_q_text_y + padding)
  
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__=="__main__":
    gene_list_name=gene_list.split("/")[-1][:-4]
    print("starting plotting ges results for gene list: ", gene_list_name)
    #upload ges_results for the gene list
    ges_results_path=output_path+"ges_results_above1/"+data_type+"/"+gene_list_name+"/full_results.csv"
    ges_results_df=pd.read_csv(ges_results_path)
    print("uploaded ges results from: ", ges_results_path)
    #make figure output folder
    figs_outpath=output_path+"figs/"+gene_list_name+"/"+data_type+"/ges_results/"
    os.makedirs(figs_outpath,exist_ok=True)#make sure path exists..
    #divide into cell class and cell cycle and run each file saparatly with the barplot
    if data_type=="cortex":
        ges_results_cell_class=ges_results_df.iloc[[10,0,3,11,9,1,4,2],:]
        output_cell_class=figs_outpath+"cell_class_enrichment.png"
        print("plotting cell class gsea results, output: ", output_cell_class)
        plot_bar_chart(ges_results_cell_class, output_cell_class, gene_list_name,"cell class")
        ges_results_cell_cycle=ges_results_df.iloc[[5,6,7,8,13,12],:]
        output_cell_cycle=figs_outpath+"cell_cycle_enrichment.png"
        print("plotting cell cycle gsea results, output: " , output_cell_cycle)
        plot_bar_chart(ges_results_cell_cycle, output_cell_cycle, gene_list_name,"cell cycle")
    else:
        ges_results_regions=ges_results_df.iloc[[5,8,6,14,19,2,15,1,18,17],:]
        output_regions=figs_outpath+"regions_enrichment.png"
        print("plotting regions gsea results , output: ", output_regions)
        plot_bar_chart(ges_results_regions,output_regions,gene_list_name,"region")
        ges_results_cell_class=ges_results_df.iloc[[7,20,3,16,9,4,13,0,10,11,12],:]
        output_cell_class=figs_outpath+"cell_class_enrichment.png"
        print("plotting cell class gsea results, output: ", output_cell_class)
        plot_bar_chart(ges_results_cell_class, output_cell_class, gene_list_name,"cell class")

    print('DONE')

    
    
