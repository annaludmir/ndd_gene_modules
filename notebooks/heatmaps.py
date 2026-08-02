import marimo

__generated_with = "0.17.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return


@app.cell
def _():
    import re
    from pathlib import Path
    import textwrap

    import pandas as pd
    import matplotlib.pyplot as plt


    def truncate_run_name(gene_list_file_name: str) -> str:
        """
        Examples:
          ".../ID_SYNDD_Intellectual_disability.csv" -> "Intellectual_disability"
          ".../ID_SYNDD_Intellectual_disability_cortex.csv" -> "Intellectual_disability"
          ".../ID_SYNDD_Microcephaly_all.csv" -> "Microcephaly"
        """
        stem = Path(str(gene_list_file_name)).stem  # filename without extension

        # remove common prefix
        stem = re.sub(r"^ID_SYNDD_", "", stem)

        # remove common trailing dataset tags (edit/extend if you have more)
        stem = re.sub(r"_(cortex|all|all_layers|data_cortex)$", "", stem)

        return stem


    def _wrap_labels(labels, width=18):
        """Wrap long tick labels into multiple lines for readability."""
        return ['\n'.join(textwrap.wrap(str(l), width=width)) for l in labels]


    def make_fdr_heatmap(
        summary_csv: str,
        data_type: str,
        *,
        run_col: str = "gene_list_file_name",
        title_col: str = "column_condition_title",
        value_col: str = "column_condition_value",
        fdr_col: str = "FDR_qval_BH",
        agg: str = "min",          # if duplicates exist, take min FDR (most significant)
        cmap: str = "YlGnBu_r",      # calmer colors
        wrap_x: bool = True,       # wrap x labels to avoid overlap
        wrap_width: int = 18,
        rotate_x: int = 45,         # if wrap_x=True, rotation=0 is best
        x_fontsize: int = 9,
        y_fontsize: int = 10,
        figsize=None,              # if None, auto-scale based on matrix size
        vmin: float = 0.0,
        vmax: float = 1.0,
    ):
        df = pd.read_csv(summary_csv)

        # axis labels
        df["run_trunc"] = (
            df[run_col]
            .apply(truncate_run_name)
            .str.replace("_", " ", regex=False)
            .str.replace(r"\s+", " ", regex=True)
            .str.strip()
        )

        df["colcond"] = (
            df[title_col].astype(str)
            + " - "
            + df[value_col].astype(str)
        ).str.replace("_", " ", regex=False).str.replace(r"\s+", " ", regex=True).str.strip()

        # ensure numeric
        df[fdr_col] = pd.to_numeric(df[fdr_col], errors="coerce")

        # pivot to matrix
        mat = df.pivot_table(
            index="colcond",
            columns="run_trunc",
            values=fdr_col,
            aggfunc=agg
        ).sort_index()

        # auto figure sizing: wide x, shorter y
        if figsize is None:
            fig_w = max(12, len(mat.columns) * 0.5)   # widen with more runs
            fig_h = max(4,  len(mat.index) * 0.45)    # keep relatively short
            figsize = (fig_w, fig_h)

        # plot
        plt.figure(figsize=figsize)
        im = plt.imshow(
            mat.values,
            aspect="auto",
            interpolation="nearest",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax
        )
        plt.colorbar(im, label="FDR", shrink=0.85)

        plt.yticks(range(len(mat.index)), mat.index, fontsize=y_fontsize)

        xlabels = mat.columns.tolist()

        plt.xticks(
            range(len(mat.columns)),
            xlabels,
            rotation=rotate_x,
            ha="right" if rotate_x else "center",
            fontsize=x_fontsize
        )

        plt.title(f"Heatmap of FDR for SYNDD data - {data_type}")
        plt.tight_layout()
        # extra bottom margin if wrapped labels are still tight
        plt.subplots_adjust(bottom=0.25 if wrap_x else 0.15)
        plt.show()

        return mat  # handy if you want to save/use later
    return make_fdr_heatmap, pd, plt, textwrap, truncate_run_name


@app.cell
def _(make_fdr_heatmap):
    mat = make_fdr_heatmap("/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/batch_summary_20260726_ID_SYNDD.csv", 'All Layers')
    return


@app.cell
def _(pd):
    summary = pd.read_csv('/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/batch_summary_20260212.csv')
    return (summary,)


@app.cell
def _(summary):
    summary
    return


@app.cell
def _(summary):
    queried_summary = summary.query("column_condition_value == 'Neural crest' & NES > 1.5")
    queried_summary
    return


@app.cell
def _():
    from matplotlib.colors import LinearSegmentedColormap

    def make_red_white_purple():
        return LinearSegmentedColormap.from_list(
            "red_white_purple",
            [
                "#6A3D9A",   # deep red (negative NES)
                "#FFFFFF",   # white center
                "#FFFFFF",
                "#B2182B"    # deep purple (positive NES)
            ],
            N=256
        )
    return (make_red_white_purple,)


@app.cell
def _():
    from matplotlib.colors import Normalize
    import numpy as np

    class NonlinearMidpointNorm(Normalize):
        """
        Non-linear norm that expands the region around zero
        so small NES values appear near-white.
        """

        def __init__(self, vmin, vmax, mid=0, softness=1.2, clip=False):
            super().__init__(vmin, vmax, clip)
            self.mid = mid
            self.softness = softness  # bigger = wider white band

        def __call__(self, value, clip=None):
            v = np.asarray(value, dtype=float)

            # scale to [-1, 1]
            x = (v - self.mid) / max(abs(self.vmax - self.mid), abs(self.mid - self.vmin))

            # non-linear transform
            x = np.sign(x) * (np.abs(x) ** self.softness)

            # back to [0, 1]
            return 0.5 * (x + 1)
    return Normalize, np


@app.function
def simplify_id_label(s: str) -> str:
    s = str(s).strip()
    s_low = s.lower()

    # Only map things that actually mention intellectual disability
    if "intellectual disability" not in s_low:
        return s  # keep e.g. "Progressive", "Microcephaly", etc.

    if "borderline" in s_low:
        return "Borderline"
    if "mild" in s_low:
        return "Mild"
    if "moderate" in s_low:
        return "Moderate"
    if "profound" in s_low:
        return "Profound"
    if "severe" in s_low:
        return "Severe"

    # plain "intellectual disability" (no modifier)
    return "ID"


@app.cell
def _(np):
    from matplotlib.colors import ListedColormap, BoundaryNorm

    def make_manual_nes_cmap_and_norm_red_white_purple(clip_nes=3.0):
        """
        Discrete NES bins with manual colors: red (neg) -> white -> purple (pos).

        Returns:
          cmap, norm, bounds
        """
        # Bin edges (edit freely). Must be sorted.
        bounds = np.array([
            -clip_nes, -2.5, -2.0, -1.5, -1.25, -1.0,
             1.0, 1.25, 1.5,  2.0,  2.5,  clip_nes
        ], dtype=float)

        # One color per interval (len(bounds) - 1).
        colors = [
            "#67001F",  # [-clip, -2.5]  darkest red
            "#B2182B",  # [-2.5, -2.0]   strong red
            "#D6604D",  # [-2.0, -1.5]   medium red
            "#F4A582",  # [-1.5, -1.0]   light red
            "#FFFFFF",
            "#FFFFFF",  # [-1.0,  1.0]   white band
            "#FFFFFF",
            "#DADAEB",  # [ 1.0,  1.5]   very light purple
            "#B2ABD2",  # [ 1.5,  2.0]   light/medium purple
            "#8073AC",  # [ 2.0,  2.5]   strong purple
            "#542788",  # [ 2.5, clip]   darkest purple
        ]

        cmap = ListedColormap(colors, name="manual_RWP")
        norm = BoundaryNorm(bounds, ncolors=cmap.N, clip=True)
        return cmap, norm, bounds
    return


@app.cell
def _(
    Normalize,
    make_red_white_purple,
    np,
    pd,
    plt,
    textwrap,
    truncate_run_name,
):
    desired_order = [
        "Intellectual disability",
        "Borderline",
        "Mild",
        "Moderate",
        "Severe",
        "Profound",
        "Progressive"
    ]

    cellclass_y_order = [
        "Radial glia",
        "Neuronal IPC",
        "Neuroblast",
        "Neuron",
        "Fibroblast",
        "Glioblast",
        "Neural crest",
        "Oligo"
    ]

    region_y_order = [
        "Forebrain",
        "Telencephalon",
        "Diencephalon",
        "Midbrain",
        "Hindbrain",
        "Pons",
        "Cerebellum",
        "Medulla"
    ]

    def reorder_columns(mat, desired_order):
        cols_present = [c for c in desired_order if c in mat.columns]
        cols_rest = [c for c in mat.columns if c not in cols_present]
        return mat[cols_present + cols_rest]

    def reorder_rows(mat, desired_order):
        rows_present = [r for r in desired_order if r in mat.index]
        rows_rest = [r for r in mat.index if r not in rows_present]
        return mat.reindex(index=rows_present + rows_rest)

    def _clean_label(s: str) -> str:
        return str(s).replace("_", " ").strip()

    def _wrap_labels(labels, width=18):
        return ['\n'.join(textwrap.wrap(str(l), width=width)) for l in labels]


    def make_nes_heatmaps_with_fdr_text(
        summary_csv: str,
        data_type: str,
        *,
        run_col: str = "gene_list_file_name",
        title_col: str = "column_condition_title",
        value_col: str = "column_condition_value",
        nes_col: str = "NES",
        fdr_col: str = "FDR_qval_BH",
        agg_nes: str = "mean",         # or "mean"/"median" depending on duplicates
        agg_fdr: str = "min",          # best (most significant) in case of duplicates
        fdr_text_cutoff: float = 0.1,  # annotate only if FDR < cutoff
        cmap: str = "RdBu_r",          # red -> yellow -> green
        wrap_x: bool = True,
        wrap_width: int = 18,
        rotate_x: int = 45,
        x_fontsize: int = 9,
        y_fontsize: int = 10,
        figsize=None,
        clip_nes: float = 2.5,         # clip color scale to [-clip_nes, clip_nes]
        separate_titles: bool = True,  # 1 heatmap per column_condition_title
        id_subset_runs=None,           # optionally pass set/list of run_trunc names to plot separately
    ):
        """
        NES shown as a continuous diverging heatmap (red negative -> green positive).
        FDR is written inside a cell only when FDR < fdr_text_cutoff.

        If separate_titles=True: plots one heatmap per column_condition_title.
        If id_subset_runs is provided: also plots a second heatmap per title with only those runs.
        Returns matrices for downstream use.
        """

        df = pd.read_csv(summary_csv)

        # -------------------------------------------------------------------------
        # FILTER: Keep only rows where data_type matches the value in "run_name"
        # -------------------------------------------------------------------------
        # Option A: Exact string match
        # df = df[df["run_name"].astype(str) == str(data_type)]

        # (Optional) Option B: If data_type is a substring inside "run_name", use this instead:
        df = df[df["run_name"].astype(str).str.contains(str(data_type), na=False)]

        if df.empty:
            print(f"Warning: No rows found matching run_name == '{data_type}'. Returning empty result.")
            return {"nes": {}, "fdr": {}, "nes_subset": {}, "fdr_subset": {}}

        # labels
        df["run_trunc"] = (
            df[run_col]
            .apply(truncate_run_name)
            .map(_clean_label)
            # remove ONLY lowercase "intellectual disability"
            .str.replace("intellectual disability", "", regex=False)
            # cleanup spacing (multiple spaces created by removal)
            .str.replace(r"\s+", " ", regex=True)
            .str.strip()
        )
        df["colcond"] = (df[value_col].astype(str)).map(_clean_label)

        # numeric
        df[nes_col] = pd.to_numeric(df[nes_col], errors="coerce")
        df[fdr_col] = pd.to_numeric(df[fdr_col], errors="coerce")

        # which titles to plot
        groups = [(None, df)] if not separate_titles else list(df.groupby(title_col, dropna=False))

        out = {"nes": {}, "fdr": {}, "nes_subset": {}, "fdr_subset": {}}

        for title, dft in groups:
            title_clean = _clean_label(title) if title is not None else "ALL"

            # pivot matrices
            mat_nes = (
                dft.pivot_table(index="colcond", columns="run_trunc", values=nes_col, aggfunc=agg_nes)
                   .sort_index()
            )
            mat_fdr = (
                dft.pivot_table(index="colcond", columns="run_trunc", values=fdr_col, aggfunc=agg_fdr)
                   .reindex(index=mat_nes.index, columns=mat_nes.columns)
            )

            # ✅ enforce x-axis order HERE (for the full plot)
            mat_nes = reorder_columns(mat_nes, desired_order)
            mat_fdr = reorder_columns(mat_fdr, desired_order)

            # ✅ enforce y-axis order depending on which title we're plotting
            if title_clean.lower() == "cellclass":
                mat_nes = reorder_rows(mat_nes, cellclass_y_order)
                mat_fdr = reorder_rows(mat_fdr, cellclass_y_order)

            if title_clean.lower() == "region":
                mat_nes = reorder_rows(mat_nes, region_y_order)
                mat_fdr = reorder_rows(mat_fdr, region_y_order)

            out["nes"][title_clean] = mat_nes
            out["fdr"][title_clean] = mat_fdr

            def _plot_one(mat_nes_plot, mat_fdr_plot, plot_title):
                # auto figure sizing
                if figsize is None:
                    fig_w = max(12, len(mat_nes_plot.columns) * 0.5)
                    fig_h = max(4,  len(mat_nes_plot.index) * 0.45)
                    _figsize = (fig_w, fig_h)
                else:
                    _figsize = figsize

                plt.figure(figsize=_figsize)

                # clamp NES for stable color scaling
                vals = np.clip(mat_nes_plot.values.astype(float), -clip_nes, clip_nes)

                # symmetric linear norm around 0
                norm = Normalize(vmin=-clip_nes, vmax=clip_nes)

                im = plt.imshow(
                    vals,
                    aspect="auto",
                    interpolation="nearest",
                    cmap=make_red_white_purple(),   # or any diverging cmap you like
                    norm=norm
                )

                cbar = plt.colorbar(im, shrink=0.85)
                cbar.set_label("NES")

                # symmetric, evenly spaced ticks
                ticks = [-2, -1, 0, 1, 2]
                ticks = [t for t in ticks if -clip_nes <= t <= clip_nes]  # keep only ticks in range
                cbar.set_ticks(ticks)
                cbar.set_ticklabels([str(t) for t in ticks])

                # ticks
                plt.yticks(range(len(mat_nes_plot.index)), mat_nes_plot.index, fontsize=y_fontsize)

                xlabels = mat_nes_plot.columns.tolist()

                plt.xticks(
                    range(len(mat_nes_plot.columns)),
                    xlabels,
                    rotation=rotate_x,
                    ha="right" if rotate_x else "center",
                    fontsize=x_fontsize
                )

                # annotate FDR only if < cutoff
                for i in range(mat_fdr_plot.shape[0]):
                    for j in range(mat_fdr_plot.shape[1]):
                        fdr_val = mat_fdr_plot.iat[i, j]
                        if pd.notna(fdr_val) and fdr_val < fdr_text_cutoff:
                            # compact formatting similar to papers
                            txt = f"{fdr_val:.2g}"
                            plt.text(j, i, txt, ha="center", va="center", fontsize=7, color="black")

                plt.title(plot_title, pad=15)
                plt.tight_layout()
                plt.subplots_adjust(bottom=0.30 if wrap_x else 0.15)
                plt.show()

            # plot full
            _plot_one(
                mat_nes, mat_fdr,
                plot_title=f"NES heatmap - {data_type} - {title_clean}"
            )

            # optional subset plot
            if id_subset_runs is not None:
                wanted = set(id_subset_runs)
                cols_present = [c for c in mat_nes.columns if c in wanted]
                if cols_present:
                    mat_nes_sub = mat_nes[cols_present].copy()
                    mat_fdr_sub = mat_fdr[cols_present].copy()
                    out["nes_subset"][title_clean] = mat_nes_sub
                    out["fdr_subset"][title_clean] = mat_fdr_sub

                    cols_present = [c for c in desired_order if c in mat_nes.columns]
                    cols_rest = [c for c in mat_nes.columns if c not in cols_present]

                    new_order = cols_present + cols_rest  # keeps anything unexpected at the end

                    mat_nes = mat_nes[new_order]
                    mat_fdr = mat_fdr[new_order]

                    _plot_one(
                        mat_nes_sub, mat_fdr_sub,
                        plot_title=f"NES heatmap - {data_type} - {title_clean}"
                    )
                else:
                    out["nes_subset"][title_clean] = pd.DataFrame()
                    out["fdr_subset"][title_clean] = pd.DataFrame()

        return out
    return (make_nes_heatmaps_with_fdr_text,)


@app.cell
def _(make_nes_heatmaps_with_fdr_text):
    id_subset_runs = {
        "Intellectual disability",
        "Borderline",
        "Mild",
        "Moderate",
        "Profound",
        "Progressive",
        "Severe",
    }


    res = make_nes_heatmaps_with_fdr_text(
        summary_csv="/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/batch_summary_20260726_ID_SYNDD.csv",
        data_type="All Layers",
        separate_titles=True,
        id_subset_runs=id_subset_runs,
        fdr_text_cutoff=0.1,
        cmap="RdBu_r",      # red (negative) -> green (positive)
        clip_nes=2.5,       # adjust if you want stronger/weaker contrast
        wrap_x=True,
        wrap_width=18
    )
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
