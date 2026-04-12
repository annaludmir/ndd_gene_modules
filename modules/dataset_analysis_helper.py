import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


DEFAULT_OUTPUT_DIR = Path("results/dataset_analysis")


def create_sample_count_by_age_plot(
    h5ad_path,
    output_dir=DEFAULT_OUTPUT_DIR,
    age_col="Age",
    sample_col="sample",
):
    """Create a bar plot showing the number of unique samples at each age."""
    h5ad_path = Path(h5ad_path)
    if not h5ad_path.exists():
        raise FileNotFoundError(f"AnnData file not found: {h5ad_path}")

    print(f"Loading AnnData from: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)

    if age_col not in adata.obs.columns:
        raise KeyError(f"Column '{age_col}' not found in adata.obs")
    if sample_col not in adata.obs.columns:
        raise KeyError(f"Column '{sample_col}' not found in adata.obs")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_df = (
        adata.obs[[age_col, sample_col]]
        .dropna()
        .assign(
            age_numeric=lambda df: pd.to_numeric(df[age_col], errors="coerce"),
            sample_value=lambda df: df[sample_col].astype(str),
        )
        .dropna(subset=["age_numeric"])
        .groupby("age_numeric")["sample_value"]
        .nunique()
        .reset_index(name="n_samples")
        .sort_values("age_numeric")
    )

    if summary_df.empty:
        raise ValueError("No rows with valid age and sample values were found.")

    csv_path = output_dir / "sample_counts_by_age.csv"
    plot_path = output_dir / "sample_counts_by_age_barplot.png"
    summary_df.to_csv(csv_path, index=False)

    fig_width = max(8, min(16, len(summary_df) * 0.55 + 3))
    fig, ax = plt.subplots(figsize=(fig_width, 5.5), constrained_layout=False)
    ax.bar(summary_df["age_numeric"].astype(str), summary_df["n_samples"], color="#3a6ea5")
    ax.set_title("Number of Samples at Each Age", fontsize=14, pad=12)
    ax.set_xlabel(age_col)
    ax.set_ylabel("Number of samples")
    ax.tick_params(axis="x", rotation=45)

    for idx, value in enumerate(summary_df["n_samples"]):
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
        "sample_column": sample_col,
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
        "--sample-col",
        default="sample",
        help="Column in adata.obs containing sample identifiers.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    result = create_sample_count_by_age_plot(
        h5ad_path=args.h5ad_path,
        output_dir=args.output_dir,
        age_col=args.age_col,
        sample_col=args.sample_col,
    )

    print("Saved dataset analysis results:")
    for key, value in result.items():
        print(f"  {key}: {value}")
