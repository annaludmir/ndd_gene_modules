import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc


DEFAULT_OUTPUT_DIR = Path("results/dataset_analysis")
DEFAULT_CHEMISTRY_COL = "Chemistry"
DEFAULT_REGION_COL = "Region"


def create_sample_count_by_age_plot(
    h5ad_path,
    output_dir=DEFAULT_OUTPUT_DIR,
    age_col="Age",
    chemistry="v3",
    chemistry_col=DEFAULT_CHEMISTRY_COL,
    region=None,
    region_col=DEFAULT_REGION_COL,
):
    """Export per-age cell counts after optional filtering."""
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

    summary_df = (
        adata.obs[age_col]
        .dropna()
        .astype(str)
        .value_counts()
        .rename_axis(age_col)
        .reset_index(name="n_cells")
    )

    if summary_df.empty:
        raise ValueError("No rows remained after filtering with valid age values.")

    csv_path = output_dir / "cell_counts_by_age.csv"
    summary_df.to_csv(csv_path, index=False)

    return {
        "output_dir": str(output_dir),
        "summary_csv_path": str(csv_path),
        "h5ad_path": str(h5ad_path),
        "age_column": age_col,
        "chemistry": chemistry,
        "chemistry_column": chemistry_col,
        "region": region,
        "region_column": region_col,
        "n_ages": int(summary_df.shape[0]),
        "age_cell_counts": summary_df.to_dict(orient="records"),
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
        "--chemistry",
        default="v3",
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
        chemistry=None if args.chemistry == "None" else args.chemistry,
        chemistry_col=args.chemistry_col,
        region=args.region,
        region_col=args.region_col,
    )

    print("Saved dataset analysis results:")
    for key, value in result.items():
        print(f"  {key}: {value}")
