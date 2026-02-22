import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_bar_chart(
    ges_results_df: pd.DataFrame,
    output_file: str,
    gene_list_name: str,
    comparison: str,
    *,
    fdr_sig: float = 0.05,
    sig_nes: float = 1.5,
    very_sig_nes: float = 2.0,
):
    """
    Bar chart of NES per condition.
    - Bars colored to match a red↔white↔purple heatmap aesthetic.
    - Annotates FDR only when significant (FDR < fdr_sig), and ALWAYS places it above the bar.
      For negative NES bars, the label is placed just above y=0 (still "above the bar" visually).
    - Removes top/right spines (keeps only left/bottom).
    - Replaces the previous orange "very significant" line with a purple-toned line.

    Expects columns:
      - 'condition'
      - 'NES'
      - an FDR column (tries several common names)
    """

    # ---- robust FDR column detection ----
    fdr_candidates = [
        "FDR_qval_BH",
        "FDR q-val (BH corrected)",
        "FDR q-val (BH corrected) ",
        "FDR q-val",
        "FDR q-val (BH)",
        "FDR q-val (BH corrected)\n",
        "Adjusted P-value",
        "Adjusted P-value ",
        "Adjusted P-value\n",
    ]
    fdr_col = next((c for c in fdr_candidates if c in ges_results_df.columns), None)

    # ---- basic checks ----
    if "NES" not in ges_results_df.columns:
        raise KeyError("ges_results_df must include a 'NES' column.")
    if "condition" not in ges_results_df.columns:
        raise KeyError("ges_results_df must include a 'condition' column for x-axis labels.")

    # ---- numeric NES ----
    nes_vals = pd.to_numeric(ges_results_df["NES"], errors="coerce").to_numpy()

    # ---- colors to match heatmap (red↔white↔purple) ----
    # strong neg (red), weak neg (light red), weak pos (light purple), strong pos (purple)
    def _nes_color(v: float) -> str:
        if not np.isfinite(v):
            return "#E0E0E0"  # NA
        if v <= -sig_nes:
            return "#ab4d4d"  # strong red
        if v < 0:
            return "#ab4d4d"  # light red
        if v >= sig_nes:
            return "#825ea6"  # strong purple
        return "#825ea6"      # light purple

    colors = [_nes_color(v) for v in nes_vals]

    # ---- x positions ----
    x_pos = np.arange(len(ges_results_df)) * 0.9

    # ---- decide legend + figure width ----
    max_nes = np.nanmax(nes_vals) if np.any(np.isfinite(nes_vals)) else np.nan
    add_legend = bool(np.isfinite(max_nes) and max_nes >= very_sig_nes)
    fig_width = 4 if not add_legend else 6
    fig, ax = plt.subplots(figsize=(fig_width, 4))

    # ---- bars ----
    bars = ax.bar(
        x_pos,
        nes_vals,
        color=colors,
        edgecolor="black",
        linewidth=0.7,
        width=0.7,
        alpha=0.9,
    )

    # ---- x labels ----
    ax.set_xticks(x_pos)
    ax.set_xticklabels(ges_results_df["condition"].astype(str), fontsize=10, rotation=80)

    # ---- y limits ----
    min_nes = np.nanmin(nes_vals) if np.any(np.isfinite(nes_vals)) else -1.0
    max_nes = np.nanmax(nes_vals) if np.any(np.isfinite(nes_vals)) else 1.0
    ymin = float(min_nes) - 0.2
    ymax = float(max_nes) + 0.2
    ax.set_ylim(ymin, ymax)

    # ---- baseline ----
    ax.axhline(0, color="black", linewidth=1)

    # ---- guide lines (significance by NES magnitude) ----
    # Significant enrichment: grey dashed
    ax.axhline(sig_nes, color="#555555", linestyle="--", linewidth=1.5,
               label="Significant enrichment" if add_legend else None)
    # Very significant enrichment: purple-toned (replaces orange)
    ax.axhline(very_sig_nes, color="#542788", linestyle="-", linewidth=2.0,
               label="Very significant enrichment" if add_legend else None)

    if add_legend:
        ax.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize=8)

    # ---- title/labels ----
    ax.set_title(f"GSEA enrichment results for {gene_list_name} - {comparison}", pad=14, fontsize=11)
    ax.set_ylabel("NES Score", fontsize=10)

    # ---- remove top/right spines ----
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(top=False, right=False)

    # ---- annotate significant FDR only (FDR < fdr_sig), always above bar ----
    max_text_y = -np.inf
    if fdr_col is not None:
        fdr_vals = pd.to_numeric(ges_results_df[fdr_col], errors="coerce").to_numpy()

        y_range = max(1e-9, (ymax - ymin))
        offset = 0.04 * y_range

        for i, rect in enumerate(bars):
            fdr = fdr_vals[i] if i < len(fdr_vals) else np.nan
            if not (np.isfinite(fdr) and fdr < fdr_sig):
                continue

            bar_h = float(rect.get_height())
            x = rect.get_x() + rect.get_width() / 2

            # ALWAYS above the bar visually:
            # - for positive bars: above bar top
            # - for negative bars: just above y=0 (still above the bar)
            y = (bar_h + offset) if (np.isfinite(bar_h) and bar_h >= 0) else (0 + offset)

            max_text_y = max(max_text_y, y)

            ax.text(
                x, y, f"q={fdr:.2g}",
                ha="center", va="bottom",
                fontsize=8,
                color="black",
                fontweight="bold",
                clip_on=False,
            )
    else:
        print("⚠️ Could not find an FDR/q-value column. Skipping FDR annotations.")

    # ---- expand ylim if needed so text doesn't touch border ----
    if np.isfinite(max_text_y) and max_text_y > ymax:
        padding = 0.06 * (ymax - ymin)
        ax.set_ylim(ymin, max_text_y + padding)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    print("DONE")
