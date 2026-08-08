import marimo

__generated_with = "0.17.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return


@app.cell
def _():
    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.stats import pearsonr
    import numpy as np

    # ---------------------------------------------------------
    # 1. Scanpy Preprocessing (Total Normalization & log1p)
    # ---------------------------------------------------------
    # Load single-cell AnnData object (or CSV if pre-extracted)
    # adata = sc.read_h5ad('/miridan-data/annaludmir/autism_strong_adata.h5ad')

    # Set aesthetic theme
    sns.set_theme(style="whitegrid")

    # Configuration toggle for v3 chemistry
    filter_v3_only = True

    # Normalize total counts per cell & log-transform expression
    # sc.pp.normalize_total(adata)
    # sc.pp.log1p(adata)

    # Load / Extract aggregated DataFrame
    df = pd.read_csv('/miridan-data/annaludmir/autism_strong_telencephalon_vs_midbrain.csv')

    # ---------------------------------------------------------
    # 2. Data Cleaning & Chemistry Filtering
    # ---------------------------------------------------------
    def clean_numeric_column(series):
        if series.dtype == object:
            series = series.astype(str).str.rstrip('%')
        return pd.to_numeric(series, errors='coerce')


    # Map expression columns dynamically
    col_gene = 'gene'
    col_med_tel = 'Telencephalon_median_expr' if 'Telencephalon_median_expr' in df.columns else 'median_expr_a'
    col_med_die = 'Midbrain_median_expr' if 'Midbrain_median_expr' in df.columns else 'median_expr_b'
    col_frac_tel = 'Telencephalon_fraction_expressing' if 'Telencephalon_fraction_expressing' in df.columns else 'Telencephalon %'
    col_frac_die = 'Midbrain_fraction_expressing' if 'Midbrain_fraction_expressing' in df.columns else 'Midbrain %'

    # Clean numeric data
    df[col_med_tel] = clean_numeric_column(df[col_med_tel])
    df[col_med_die] = clean_numeric_column(df[col_med_die])
    df[col_frac_tel] = clean_numeric_column(df[col_frac_tel])
    df[col_frac_die] = clean_numeric_column(df[col_frac_die])


    # ---------------------------------------------------------
    # 3. Plot 1: Median Expression with Outliers & Identity Line f(x)=x
    # ---------------------------------------------------------
    p1_df = df.dropna(subset=[col_med_die, col_med_tel, col_gene]).copy()

    fig, ax1 = plt.subplots(figsize=(8, 6.5))

    x1 = p1_df[col_med_die]
    y1 = p1_df[col_med_tel]
    r1, p1 = pearsonr(x1, y1)

    poly1 = np.polyfit(x1, y1, 1)
    y1_pred = np.polyval(poly1, x1)
    res1 = np.abs(y1 - y1_pred)

    outliers1_df = p1_df.loc[res1.nlargest(5).index]

    # Main scatter plot + regression line
    sns.regplot(
        x=x1, y=y1, data=p1_df, ax=ax1,
        color='#2563eb',
        scatter_kws={'alpha': 0.6, 's': 55},
        line_kws={'color': '#dc2626', 'linewidth': 2, 'label': 'Regression fit'}
    )

    # Identity Line f(x) = x
    lims1 = [min(ax1.get_xlim()[0], ax1.get_ylim()[0]), max(ax1.get_xlim()[1], ax1.get_ylim()[1])]
    ax1.plot(lims1, lims1, color='#64748b', linestyle='--', linewidth=1.8, label='f(x) = x', zorder=2)
    ax1.set_xlim(lims1)
    ax1.set_ylim(lims1)

    # Outlier points
    ax1.scatter(
        outliers1_df[col_med_die],
        outliers1_df[col_med_tel],
        color='#f59e0b', edgecolor='#b45309', s=70, zorder=5, label='Outlier Genes'
    )

    # Annotations
    for i, (idx, row) in enumerate(outliers1_df.iterrows()):
        gene = row[col_gene]
        gx, gy = row[col_med_die], row[col_med_tel]
        offset_x = 15 if i % 2 == 0 else -45
        offset_y = 12 if i % 3 == 0 else -18
        ax1.annotate(
            gene,
            xy=(gx, gy),
            xytext=(gx + offset_x * 0.003, gy + offset_y * 0.003),
            fontsize=9, fontweight='bold', color='#0f172a',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#fef3c7', edgecolor='#f59e0b', alpha=0.9),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2', color='#b45309', lw=1.2)
        )

    ax1.set_title(f'Gene Expression: Telencephalon vs Midbrain (Median Expression)', fontsize=12, fontweight='bold', pad=12)
    ax1.set_xlabel('Midbrain Median Expression', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Telencephalon Median Expression', fontsize=11, fontweight='bold')

    p_text1 = f"p < 0.001" if p1 < 0.001 else f"p = {p1:.4f}"
    ax1.text(0.05, 0.88, f'Pearson r = {r1:.3f}\n{p_text1}', transform=ax1.transAxes,
             fontsize=11, fontweight='bold', bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='#cbd5e1'))

    ax1.legend(loc='lower right', frameon=True)
    plt.tight_layout()
    plt.savefig('plot1_median_expr_tagged_v3.png', dpi=300)
    plt.show()


    # ---------------------------------------------------------
    # 4. Plot 2: Fraction Expressing with Outliers & Identity Line f(x)=x
    # ---------------------------------------------------------
    p2_df = df.dropna(subset=[col_frac_die, col_frac_tel, col_gene]).copy()

    fig, ax2 = plt.subplots(figsize=(8, 6.5))

    x2 = p2_df[col_frac_die]
    y2 = p2_df[col_frac_tel]
    r2, p2 = pearsonr(x2, y2)

    poly2 = np.polyfit(x2, y2, 1)
    y2_pred = np.polyval(poly2, x2)
    res2 = np.abs(y2 - y2_pred)

    outliers2_df = p2_df.loc[res2.nlargest(5).index]

    # Main scatter plot + regression line
    sns.regplot(
        x=x2, y=y2, data=p2_df, ax=ax2,
        color='#0d9488',
        scatter_kws={'alpha': 0.6, 's': 55},
        line_kws={'color': '#dc2626', 'linewidth': 2, 'label': 'Regression fit'}
    )

    # Identity Line f(x) = x
    lims2 = [min(ax2.get_xlim()[0], ax2.get_ylim()[0]), max(ax2.get_xlim()[1], ax2.get_ylim()[1])]
    ax2.plot(lims2, lims2, color='#64748b', linestyle='--', linewidth=1.8, label='f(x) = x', zorder=2)
    ax2.set_xlim(lims2)
    ax2.set_ylim(lims2)

    # Outlier points
    ax2.scatter(
        outliers2_df[col_frac_die],
        outliers2_df[col_frac_tel],
        color='#f59e0b', edgecolor='#b45309', s=70, zorder=5, label='Outlier Genes'
    )

    # Annotations
    for i, (idx, row) in enumerate(outliers2_df.iterrows()):
        gene = row[col_gene]
        gx, gy = row[col_frac_die], row[col_frac_tel]
        offset_x = 20 if i % 2 == 0 else -45
        offset_y = 15 if i % 3 == 0 else -20
        ax2.annotate(
            gene,
            xy=(gx, gy),
            xytext=(gx + offset_x * 0.003, gy + offset_y * 0.003),
            fontsize=9, fontweight='bold', color='#0f172a',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#fef3c7', edgecolor='#f59e0b', alpha=0.9),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2', color='#b45309', lw=1.2)
        )

    ax2.set_title(f'Gene Expression: Telencephalon vs Midbrain (Fraction Expressing)', fontsize=12, fontweight='bold', pad=12)
    ax2.set_xlabel('Midbrain Fraction Expressing', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Telencephalon Fraction Expressing', fontsize=11, fontweight='bold')

    p_text2 = f"p < 0.001" if p2 < 0.001 else f"p = {p2:.4f}"
    ax2.text(0.05, 0.88, f'Pearson r = {r2:.3f}\n{p_text2}', transform=ax2.transAxes,
             fontsize=11, fontweight='bold', bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='#cbd5e1'))

    ax2.legend(loc='lower right', frameon=True)
    plt.tight_layout()
    plt.savefig('plot2_fraction_expressing_tagged_v3.png', dpi=300)
    plt.show()
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
