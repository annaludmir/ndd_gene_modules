import marimo

__generated_with = "0.17.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return


@app.cell
def _():
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import seaborn as sns
    from scipy.stats import pearsonr
    return np, pd, pearsonr, plt, sc, sns


@app.cell
def _(sc):
    a=sc.read_h5ad('/miridan-storage/annaludmir/atac-seq/e88a34d0-d28a-4d10-a8c5-d59f86ba621a.h5ad', backed='r')
    print(a.obs_names[:3].tolist())
    return


@app.cell
def _(sc):
    # 1. Load Data
    h5ad_file_path = "/miridan-data/annaludmir/ndd_gene_modules/data/human_dev_without_week_5.h5ad"

    adata_raw_input = sc.read_h5ad(h5ad_file_path)
    return (adata_raw_input,)


@app.cell
def _(adata_raw_input, np, sc):
    # 2. Filter Chemistry == 'v3'
    chem_col = "Chemistry"
    if chem_col in adata_raw_input.obs.columns:
        adata_v3_only = adata_raw_input[
            adata_raw_input.obs[chem_col].astype(str).str.lower() == "v3"
        ].copy()
    else:
        adata_v3_only = adata_raw_input.copy()

    # 3. Filter CellClass in Proliferating List
    proliferating_cell_types = ["Radial glia", "Neuronal IPC", "Glioblast"]
    cell_class_key = (
        "CellClass" if "CellClass" in adata_v3_only.obs.columns else "cell_class"
    )

    adata_proliferating = adata_v3_only[
        adata_v3_only.obs[cell_class_key].isin(proliferating_cell_types)
    ].copy()

    # 4. Derive RegionsGeneral from Region Column
    region_obs_key = (
        "Region" if "Region" in adata_proliferating.obs.columns else "region"
    )

    forebrain_members = {"Diencephalon", "Forebrain", "Telencephalon"}
    midbrain_members = {"Midbrain"}
    hindbrain_members = {"Hindbrain", "Cerebellum", "Pons", "Medulla"}


    def assign_general_region(region_name):
        if region_name in forebrain_members:
            return "Forebrain_General"
        elif region_name in midbrain_members:
            return "Midbrain_General"
        elif region_name in hindbrain_members:
            return "Hindbrain_General"
        return np.nan


    adata_proliferating.obs["RegionsGeneral"] = (
        adata_proliferating.obs[region_obs_key]
        .astype(str)
        .map(assign_general_region)
    )

    # 5. Keep target regions
    target_derived_regions = ["Forebrain_General", "Midbrain_General"]
    adata_filtered_final = adata_proliferating[
        adata_proliferating.obs["RegionsGeneral"].isin(target_derived_regions)
    ].copy()

    # 6. Normalize and log-transform
    sc.pp.normalize_total(adata_filtered_final)
    sc.pp.log1p(adata_filtered_final)
    return (adata_filtered_final,)


@app.cell
def _(adata_filtered_final, np, pd):
    # ---------------------------------------------------------
    # 2. Calculate Median Expression & Fraction Expressing
    # ---------------------------------------------------------
    # Define your target gene list
    target_gene_list = [
        "FOXG1",
        "CDK6",
        "HPDL",
        "BRIP1",
        "CEP152",
        "CIT",
        "BLM",
        "RAD51",
        "KIF11",
        "KNL1",
        "FANCI",
        "NCAPH",
        "CDT1",
        "ORC1",
        "POC1A",
        "CDC6",
        "BUB1B",
        "DNA2",
        "CKAP2L",
        "NCAPD2",
        "ZEB2",
        "BUB1",
        "STIL",
        "CENPF",
        "PLK4",
        "ASPM",
        "FANCB",
        "BRCA2",
        "WDR62",
        "CENPE",
        "LMNB2",
        "CEP55",
        "FANCD2",
        "HHAT",
        "FANCA",
        "CCDC88A",
        "GINS2",
        "MRE11",
        "NDE1",
        "CCND2",
        "DDX11",
        "TRAIP",
        "RBBP8",
        "FANCL",
        "TRIP13",
        "GINS3",
        "CEP135",
        "NUF2",
        "FANCE",
        "GMNN",
        "TUBGCP3",
        "CDK5RAP2",
        "ORC6",
        "KIF14",
        "VRK1",
        "MCM7",
        "LMNB1",
        "SMC1A",
        "DIAPH1",
        "PCNT",
        "NUP188",
        "NUP107",
        "FANCG",
        "GPT2",
        "NCAPD3",
        "TCF4",
        "NUP214",
        "SASS6",
        "PUS7",
        "SMO",
        "RAD50",
        "TTI1",
        "SMC5",
        "ZNF526",
        "RAD21",
        "ATP1A2",
        "TOP3A",
        "NSD2",
        "FBRSL1",
        "SLF2",
        "ERCC5",
        "OSGEP",
        "SMC3",
        "EFTUD2",
        "ANKLE2",
        "TPR",
        "POGZ",
        "MCPH1",
        "RMI1",
        "TUBGCP4",
        "CPSF3",
        "IARS1",
        "CAMSAP1",
        "PSMC3",
        "CTCF",
        "UGP2",
        "PPP1R35",
        "UFC1",
        "DROSHA",
        "ATRIP",
        "RTTN",
        "ATRX",
        "ERCC6",
        "SLC9A6",
        "FRA10AC1",
        "CEP57",
        "CREBBP",
        "BPTF",
        "YIPF5",
        "PTPN23",
        "CHAMP1",
        "NIPBL",
        "PDCD6IP",
        "PALB2",
        "TUBG1",
        "ACBD6",
        "DYRK1A",
        "EIF2S3",
        "DONSON",
        "LARP7",
        "PNKP",
        "NSRP1",
        "EXOC7",
        "SMARCA5",
        "XRCC4",
        "BRD4",
        "AARS1",
        "CEP63",
        "PHC1",
        "ARCN1",
        "FANCC",
        "TSEN54",
        "FANCM",
        "AKT3",
        "TRIO",
        "PPIL1",
        "PRIM1",
        "CASK",
        "HDAC8",
        "NSMCE2",
        "MED17",
        "CTU2",
        "ATP11A",
        "UBE3A",
        "IGF1R",
        "COPB1",
        "ERCC8",
        "EIF5A",
    ]

    # Map gene symbols from adata.var['Gene'] if present
    sym_col = "Gene" if "Gene" in adata_filtered_final.var.columns else None
    if sym_col:
        extracted_gene_names = (
            adata_filtered_final.var[sym_col].astype(str).values
        )
    else:
        extracted_gene_names = adata_filtered_final.var_names.values

    # Convert sparse matrix X to dense numpy array
    expression_matrix_x = adata_filtered_final.X
    if hasattr(expression_matrix_x, "toarray"):
        expression_matrix_x = expression_matrix_x.toarray()

    # Regional masks
    mask_forebrain_cells = (
        adata_filtered_final.obs["RegionsGeneral"] == "Forebrain_General"
    ).values
    mask_midbrain_cells = (
        adata_filtered_final.obs["RegionsGeneral"] == "Midbrain_General"
    ).values

    matrix_forebrain = expression_matrix_x[mask_forebrain_cells, :]
    matrix_midbrain = expression_matrix_x[mask_midbrain_cells, :]

    # Calculate Median & Fraction Expressing
    forebrain_medians = np.median(matrix_forebrain, axis=0)
    midbrain_medians = np.median(matrix_midbrain, axis=0)

    # Calculate Fraction Expressing (across all cells)
    forebrain_fractions = (matrix_forebrain > 0).mean(axis=0) * 100
    midbrain_fractions = (matrix_midbrain > 0).mean(axis=0) * 100

    # Calculate Median Expression ONLY in Expressing Cells (> 0)
    n_genes = expression_matrix_x.shape[1]
    forebrain_medians = np.zeros(n_genes)
    midbrain_medians = np.zeros(n_genes)

    for i in range(n_genes):
        fb_expressed = matrix_forebrain[matrix_forebrain[:, i] > 0, i]
        forebrain_medians[i] = np.median(fb_expressed) if len(fb_expressed) > 0 else 0.0

        mb_expressed = matrix_midbrain[matrix_midbrain[:, i] > 0, i]
        midbrain_medians[i] = np.median(mb_expressed) if len(mb_expressed) > 0 else 0.0

    # Assemble Summary DataFrame
    regional_metrics_df = pd.DataFrame(
        {
            "gene": extracted_gene_names,
            "Forebrain_median_expr": forebrain_medians,
            "Midbrain_median_expr": midbrain_medians,
            "Forebrain_fraction_expressing": forebrain_fractions,
            "Midbrain_fraction_expressing": midbrain_fractions,
        }
    )

    # ==============================================================================
    # INSERT FILTER HERE: Keep ONLY your list of genes (case-insensitive match)
    # ==============================================================================
    target_genes_upper = [g.upper() for g in target_gene_list]

    regional_metrics_df = regional_metrics_df[
        regional_metrics_df["gene"].astype(str).str.upper().isin(target_genes_upper)
    ].copy()
    # ==============================================================================

    # Keep expressed genes among the target list
    regional_metrics_df = regional_metrics_df[
        (regional_metrics_df["Forebrain_fraction_expressing"] > 0)
        | (regional_metrics_df["Midbrain_fraction_expressing"] > 0)
    ].copy()

    target_gene_symbol = "NDE1"
    return regional_metrics_df, target_gene_symbol


@app.cell
def _(np, pearsonr, plt, regional_metrics_df, sns, target_gene_symbol):
    plot1_clean_df = regional_metrics_df.dropna(
        subset=["Forebrain_median_expr", "Midbrain_median_expr", "gene"]
    ).copy()

    fig1, ax1_instance = plt.subplots(figsize=(8, 6.5), dpi=300)

    x1_data = plot1_clean_df["Midbrain_median_expr"]
    y1_data = plot1_clean_df["Forebrain_median_expr"]
    r1_stat, p1_stat = pearsonr(x1_data, y1_data)

    poly1_fit = np.polyfit(x1_data, y1_data, 1)
    y1_predicted = np.polyval(poly1_fit, x1_data)
    residuals_p1 = np.abs(y1_data - y1_predicted)

    top_outliers_p1 = plot1_clean_df.loc[residuals_p1.nlargest(5).index]

    sns.regplot(
        x=x1_data,
        y=y1_data,
        data=plot1_clean_df,
        ax=ax1_instance,
        color="#2563eb",
        scatter_kws={"alpha": 0.5, "s": 50},
        line_kws={"color": "#dc2626", "linewidth": 2, "label": "Regression fit"},
    )

    # Identity Line
    axis1_limits = [
        min(ax1_instance.get_xlim()[0], ax1_instance.get_ylim()[0]),
        max(ax1_instance.get_xlim()[1], ax1_instance.get_ylim()[1]),
    ]
    ax1_instance.plot(
        axis1_limits,
        axis1_limits,
        color="#64748b",
        linestyle="--",
        linewidth=1.8,
        label="f(x) = x",
        zorder=2,
    )
    ax1_instance.set_xlim(axis1_limits)
    ax1_instance.set_ylim(axis1_limits)

    # Outliers
    ax1_instance.scatter(
        top_outliers_p1["Midbrain_median_expr"],
        top_outliers_p1["Forebrain_median_expr"],
        color="#f59e0b",
        edgecolor="#b45309",
        s=70,
        zorder=5,
        label="Outlier Genes",
    )

    for idx_1, row_p1 in top_outliers_p1.reset_index().iterrows():
        gene_p1 = row_p1["gene"]
        gx1, gy1 = row_p1["Midbrain_median_expr"], row_p1["Forebrain_median_expr"]
        off_x1 = 15 if idx_1 % 2 == 0 else -45
        off_y1 = 12 if idx_1 % 3 == 0 else -18
        ax1_instance.annotate(
            gene_p1,
            xy=(gx1, gy1),
            xytext=(gx1 + off_x1 * 0.003, gy1 + off_y1 * 0.003),
            fontsize=9,
            fontweight="bold",
            color="#0f172a",
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor="#fef3c7",
                edgecolor="#f59e0b",
                alpha=0.9,
            ),
            arrowprops=dict(
                arrowstyle="->",
                connectionstyle="arc3,rad=0.2",
                color="#b45309",
                lw=1.2,
            ),
        )

    # Target Gene (NDE1)
    nde1_data_p1 = plot1_clean_df[
        plot1_clean_df["gene"].astype(str).str.upper()
        == target_gene_symbol.upper()
    ]
    if not nde1_data_p1.empty:
        gx_nde1 = nde1_data_p1["Midbrain_median_expr"].values[0]
        gy_nde1 = nde1_data_p1["Forebrain_median_expr"].values[0]

        ax1_instance.scatter(
            gx_nde1,
            gy_nde1,
            color="#10b981",
            edgecolor="#047857",
            s=110,
            zorder=6,
            label=f"Target: {target_gene_symbol}",
        )
        ax1_instance.annotate(
            target_gene_symbol,
            xy=(gx_nde1, gy_nde1),
            xytext=(gx_nde1 - 0.05, gy_nde1 + 0.08),
            fontsize=10,
            fontweight="bold",
            color="#064e3b",
            bbox=dict(
                boxstyle="round,pad=0.35",
                facecolor="#d1fae5",
                edgecolor="#10b981",
                alpha=0.95,
            ),
            arrowprops=dict(
                arrowstyle="->",
                connectionstyle="arc3,rad=-0.2",
                color="#047857",
                lw=1.5,
            ),
        )

    ax1_instance.set_title(
        "Proliferating Cells: Forebrain vs Midbrain (Median Expression)",
        fontsize=12,
        fontweight="bold",
        pad=12,
    )
    ax1_instance.set_xlabel(
        "Midbrain Median Expression", fontsize=11, fontweight="bold"
    )
    ax1_instance.set_ylabel(
        "Forebrain Median Expression", fontsize=11, fontweight="bold"
    )

    p_text_format1 = f"p < 0.001" if p1_stat < 0.001 else f"p = {p1_stat:.4f}"
    ax1_instance.text(
        0.05,
        0.85,
        f"Pearson r = {r1_stat:.3f}\n{p_text_format1}",
        transform=ax1_instance.transAxes,
        fontsize=11,
        fontweight="bold",
        bbox=dict(
            boxstyle="round,pad=0.5",
            facecolor="white",
            alpha=0.9,
            edgecolor="#cbd5e1",
        ),
    )

    ax1_instance.legend(loc="lower right", frameon=True)
    plt.tight_layout()

    # Returning the fig object directly renders it cleanly in Marimo
    fig1
    return


@app.cell
def _(ax2, np, pearsonr, plt, regional_metrics_df, sns, target_gene_symbol):
    plot2_clean_df = regional_metrics_df.dropna(
        subset=[
            "Midbrain_fraction_expressing",
            "Forebrain_fraction_expressing",
            "gene",
        ]
    ).copy()

    fig2, ax2_instance = plt.subplots(figsize=(8, 6.5), dpi=300)

    x2_data = plot2_clean_df["Midbrain_fraction_expressing"]
    y2_data = plot2_clean_df["Forebrain_fraction_expressing"]
    r2_stat, p2_stat = pearsonr(x2_data, y2_data)

    poly2_fit = np.polyfit(x2_data, y2_data, 1)
    y2_predicted = np.polyval(poly2_fit, x2_data)
    residuals_p2 = np.abs(y2_data - y2_predicted)

    top_outliers_p2 = plot2_clean_df.loc[residuals_p2.nlargest(5).index]

    sns.regplot(
        x=x2_data,
        y=y2_data,
        data=plot2_clean_df,
        ax=ax2_instance,
        color="#0d9488",
        scatter_kws={"alpha": 0.5, "s": 50},
        line_kws={"color": "#dc2626", "linewidth": 2, "label": "Regression fit"},
    )

    axis2_limits = [
        min(ax2_instance.get_xlim()[0], ax2_instance.get_ylim()[0]),
        max(ax2_instance.get_xlim()[1], ax2.get_ylim()[1] if 'ax2' in locals() else ax2_instance.get_ylim()[1]),
    ]
    ax2_instance.plot(
        axis2_limits,
        axis2_limits,
        color="#64748b",
        linestyle="--",
        linewidth=1.8,
        label="f(x) = x",
        zorder=2,
    )
    ax2_instance.set_xlim(axis2_limits)
    ax2_instance.set_ylim(axis2_limits)

    # Outliers
    ax2_instance.scatter(
        top_outliers_p2["Midbrain_fraction_expressing"],
        top_outliers_p2["Forebrain_fraction_expressing"],
        color="#f59e0b",
        edgecolor="#b45309",
        s=70,
        zorder=5,
        label="Outlier Genes",
    )

    for idx_2, row_p2 in top_outliers_p2.reset_index().iterrows():
        gene_p2 = row_p2["gene"]
        gx2, gy2 = (
            row_p2["Midbrain_fraction_expressing"],
            row_p2["Forebrain_fraction_expressing"],
        )
        off_x2 = 20 if idx_2 % 2 == 0 else -45
        off_y2 = 15 if idx_2 % 3 == 0 else -20
        ax2_instance.annotate(
            gene_p2,
            xy=(gx2, gy2),
            xytext=(gx2 + off_x2 * 0.003, gy2 + off_y2 * 0.003),
            fontsize=9,
            fontweight="bold",
            color="#0f172a",
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor="#fef3c7",
                edgecolor="#f59e0b",
                alpha=0.9,
            ),
            arrowprops=dict(
                arrowstyle="->",
                connectionstyle="arc3,rad=0.2",
                color="#b45309",
                lw=1.2,
            ),
        )

    # Target Gene (NDE1)
    nde1_data_p2 = plot2_clean_df[
        plot2_clean_df["gene"].astype(str).str.upper()
        == target_gene_symbol.upper()
    ]
    if not nde1_data_p2.empty:
        gx_nde1_2 = nde1_data_p2["Midbrain_fraction_expressing"].values[0]
        gy_nde1_2 = nde1_data_p2["Forebrain_fraction_expressing"].values[0]

        ax2_instance.scatter(
            gx_nde1_2,
            gy_nde1_2,
            color="#10b981",
            edgecolor="#047857",
            s=110,
            zorder=6,
            label=f"Target: {target_gene_symbol}",
        )
        ax2_instance.annotate(
            target_gene_symbol,
            xy=(gx_nde1_2, gy_nde1_2),
            xytext=(gx_nde1_2 - 5.0, gy_nde1_2 + 6.0),
            fontsize=10,
            fontweight="bold",
            color="#064e3b",
            bbox=dict(
                boxstyle="round,pad=0.35",
                facecolor="#d1fae5",
                edgecolor="#10b981",
                alpha=0.95,
            ),
            arrowprops=dict(
                arrowstyle="->",
                connectionstyle="arc3,rad=-0.2",
                color="#047857",
                lw=1.5,
            ),
        )

    ax2_instance.set_title(
        "Proliferating Cells: Forebrain vs Midbrain (Fraction Expressing)",
        fontsize=12,
        fontweight="bold",
        pad=12,
    )
    ax2_instance.set_xlabel(
        "Midbrain Fraction Expressing (%)", fontsize=11, fontweight="bold"
    )
    ax2_instance.set_ylabel(
        "Forebrain Fraction Expressing (%)", fontsize=11, fontweight="bold"
    )

    p_text_format2 = f"p < 0.001" if p2_stat < 0.001 else f"p = {p2_stat:.4f}"
    ax2_instance.text(
        0.05,
        0.85,
        f"Pearson r = {r2_stat:.3f}\n{p_text_format2}",
        transform=ax2_instance.transAxes,
        fontsize=11,
        fontweight="bold",
        bbox=dict(
            boxstyle="round,pad=0.5",
            facecolor="white",
            alpha=0.9,
            edgecolor="#cbd5e1",
        ),
    )

    ax2_instance.legend(loc="lower right", frameon=True)
    plt.tight_layout()

    fig2
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
