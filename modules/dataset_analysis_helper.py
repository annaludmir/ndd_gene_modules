import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


DEFAULT_OUTPUT_DIR = Path("results/dataset_analysis")
DEFAULT_CHEMISTRY_COL = "Chemistry"
DEFAULT_REGION_COL = "Region"


def create_sample_count_by_age_plot(
    h5ad_path,
    output_dir=DEFAULT_OUTPUT_DIR,
    age_col="Age",
    cell_id_col="CellID",
    chemistry=None,
    chemistry_col=DEFAULT_CHEMISTRY_COL,
    region=None,
    region_col=DEFAULT_REGION_COL,
):
    """Create a bar plot showing the number of unique CellIDs at each age."""
    h5ad_path = Path(h5ad_path)
    if not h5ad_path.exists():
        raise FileNotFoundError(f"AnnData file not found: {h5ad_path}")

    print(f"Loading AnnData from: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)

    if age_col not in adata.obs.columns:
        raise KeyError(f"Column '{age_col}' not found in adata.obs")
    if chemistry is not None:
        if chemistry_col not in adata.obs.columns:
            raise KeyError(f"Column '{chemistry_col}' not found in adata.obs")
        adata = adata[adata.obs[chemistry_col].astype(str) == str(chemistry)].copy()
    if region is not None:
        if region_col not in adata.obs.columns:
            raise KeyError(f"Column '{region_col}' not found in adata.obs")
        adata = adata[adata.obs[region_col].astype(str) == str(region)].copy()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if cell_id_col in adata.obs.columns:
        cell_id_series = adata.obs[cell_id_col].astype(str)
        resolved_cell_id_source = f"obs['{cell_id_col}']"
    elif adata.obs.index.name == cell_id_col:
        cell_id_series = pd.Series(adata.obs.index.astype(str), index=adata.obs.index)
        resolved_cell_id_source = "obs.index"
    else:
        raise KeyError(
            f"Cell identifier '{cell_id_col}' was not found in adata.obs columns "
            f"or as the named obs index."
        )

    summary_df = (
        pd.DataFrame(
            {
                age_col: adata.obs[age_col],
                "cell_id_value": cell_id_series,
            },
            index=adata.obs.index,
        )
        .dropna()
        .assign(age_numeric=lambda df: pd.to_numeric(df[age_col], errors="coerce"))
        .dropna(subset=["age_numeric"])
        .groupby("age_numeric")["cell_id_value"]
        .nunique()
        .reset_index(name="n_cells")
        .sort_values("age_numeric")
    )

    if summary_df.empty:
        raise ValueError("No rows remained after filtering with valid age and CellID values.")

    csv_path = output_dir / "cell_counts_by_age.csv"
    plot_path = output_dir / "cell_counts_by_age_barplot.png"
    summary_df.to_csv(csv_path, index=False)

    fig_width = max(8, min(16, len(summary_df) * 0.55 + 3))
    fig, ax = plt.subplots(figsize=(fig_width, 5.5), constrained_layout=False)
    ax.bar(summary_df["age_numeric"].astype(str), summary_df["n_cells"], color="#3a6ea5")
    ax.set_title("Number of Cells at Each Age", fontsize=14, pad=12)
    ax.set_xlabel(age_col)
    ax.set_ylabel("Number of cells")
    ax.tick_params(axis="x", rotation=45)

    for idx, value in enumerate(summary_df["n_cells"]):
        ax.text(idx, value, str(value), ha="center", va="bottom", fontsize=9)

    fig.subplots_adjust(left=0.1, right=0.96, top=0.88, bottom=0.24)
    fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return {
        "output_dir": str(output_dir),
        "summary_csv_path": str(csv_path),
        "barplot_path": str(plot_path),
        "h5ad_path": str(h5ad_path),
        "age_column": age_col,
        "cell_id_column": cell_id_col,
        "cell_id_source": resolved_cell_id_source,
        "chemistry": chemistry,
        "chemistry_column": chemistry_col,
        "region": region,
        "region_column": region_col,
        "n_ages": int(summary_df.shape[0]),
    }


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Create helper dataset-analysis plots and tables."
    )
    parser.add_argument(
        "--h5ad-path",
        required=True,
        help="Path to the input AnnData .h5ad file.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Folder where the dataset-analysis outputs will be saved.",
    )
    parser.add_argument(
        "--age-col",
        default="Age",
        help="Column in adata.obs containing age values.",
    )
    parser.add_argument(
        "--cell-id-col",
        default="CellID",
        help="Column in adata.obs containing cell identifiers.",
    )
    parser.add_argument(
        "--chemistry",
        default=None,
        help="Optional chemistry value to keep from adata.obs before plotting.",
    )
    parser.add_argument(
        "--chemistry-col",
        default=DEFAULT_CHEMISTRY_COL,
        help="Column in adata.obs containing chemistry labels.",
    )
    parser.add_argument(
        "--region",
        default=None,
        help="Optional region value to keep from adata.obs before plotting.",
    )
    parser.add_argument(
        "--region-col",
        default=DEFAULT_REGION_COL,
        help="Column in adata.obs containing region labels.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    result = create_sample_count_by_age_plot(
        h5ad_path=args.h5ad_path,
        output_dir=args.output_dir,
        age_col=args.age_col,
        cell_id_col=args.cell_id_col,
        chemistry=args.chemistry,
        chemistry_col=args.chemistry_col,
        region=args.region,
        region_col=args.region_col,
    )

    print("Saved dataset analysis results:")
    for key, value in result.items():
        print(f"  {key}: {value}")
