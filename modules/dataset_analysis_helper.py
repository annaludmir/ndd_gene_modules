import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


DEFAULT_OUTPUT_DIR = Path("results/dataset_analysis")
DEFAULT_CHEMISTRY_COL = "Chemistry"
DEFAULT_REGION_COL = "Region"
DEFAULT_REGION_GROUP_COL = "Region_group"
REGION_MAP = {
    "Forebrain": ["Forebrain", "Diencephalon", "Telencephalon"],
    "Midbrain": ["Midbrain"],
    "Hindbrain": ["Hindbrain", "Cerebellum", "Pons", "Medulla"],
}


def map_region_to_group(region_value):
    """Collapse detailed regions into Forebrain, Midbrain, or Hindbrain."""
    if pd.isna(region_value):
        return None
    for region_group, region_values in REGION_MAP.items():
        if region_value in region_values:
            return region_group
    return None


def load_filtered_adata(
    h5ad_path,
    chemistry="v3",
    chemistry_col=DEFAULT_CHEMISTRY_COL,
):
    """Load AnnData and optionally filter by chemistry."""
    h5ad_path = Path(h5ad_path)
    if not h5ad_path.exists():
        raise FileNotFoundError(f"AnnData file not found: {h5ad_path}")

    print(f"Loading AnnData from: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)

    if chemistry is not None:
        if chemistry_col not in adata.obs.columns:
            raise KeyError(f"Column '{chemistry_col}' not found in adata.obs")
        adata = adata[adata.obs[chemistry_col].astype(str) == str(chemistry)].copy()

    return adata, h5ad_path


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
    adata, h5ad_path = load_filtered_adata(
        h5ad_path,
        chemistry=chemistry,
        chemistry_col=chemistry_col,
    )

    if age_col not in adata.obs.columns:
        raise KeyError(f"Column '{age_col}' not found in adata.obs")
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


def create_gene_expression_boxplot_by_region(
    h5ad_path,
    genes,
    output_dir=DEFAULT_OUTPUT_DIR,
    chemistry="v3",
    chemistry_col=DEFAULT_CHEMISTRY_COL,
    region_col=DEFAULT_REGION_COL,
    region_group_col=DEFAULT_REGION_GROUP_COL,
):
    """Create one boxplot showing requested genes across Forebrain, Midbrain, and Hindbrain."""
    adata, h5ad_path = load_filtered_adata(
        h5ad_path,
        chemistry=chemistry,
        chemistry_col=chemistry_col,
    )

    if region_col not in adata.obs.columns:
        raise KeyError(f"Column '{region_col}' not found in adata.obs")
    if "Gene" not in adata.var.columns:
        raise KeyError("Column 'Gene' not found in adata.var")

    cleaned_genes = []
    seen = set()
    for gene in genes:
        gene = str(gene).strip()
        if gene and gene not in seen:
            seen.add(gene)
            cleaned_genes.append(gene)
    if not cleaned_genes:
        raise ValueError("At least one gene must be provided.")

    adata = adata.copy()
    adata.obs[region_group_col] = adata.obs[region_col].map(map_region_to_group)
    adata = adata[adata.obs[region_group_col].notna(), :].copy()
    if adata.n_obs == 0:
        raise ValueError("No cells remained after collapsing regions to Forebrain/Midbrain/Hindbrain.")

    gene_to_var = (
        pd.Series(adata.var_names.values, index=adata.var["Gene"].astype(str))
        .dropna()
        .to_dict()
    )
    found_genes = [gene for gene in cleaned_genes if gene in gene_to_var]
    missing_genes = [gene for gene in cleaned_genes if gene not in gene_to_var]
    if not found_genes:
        raise ValueError("None of the requested genes were found in adata.var['Gene'].")

    expr_data = {}
    for gene in found_genes:
        gene_matrix = adata[:, gene_to_var[gene]].X
        # AnnData can return dense arrays, scipy sparse matrices, or sparse view wrappers.
        if hasattr(gene_matrix, "toarray"):
            gene_values = gene_matrix.toarray()
        elif hasattr(gene_matrix, "A"):
            gene_values = gene_matrix.A
        else:
            gene_values = np.asarray(gene_matrix)
        expr_data[gene] = np.asarray(gene_values).reshape(-1)

    expr_df = pd.DataFrame(expr_data)
    expr_df[region_group_col] = adata.obs[region_group_col].to_numpy()

    expr_long = expr_df.melt(
        id_vars=region_group_col,
        var_name="Gene",
        value_name="Expression",
    )
    expr_long[region_group_col] = pd.Categorical(
        expr_long[region_group_col],
        categories=["Forebrain", "Midbrain", "Hindbrain"],
        ordered=True,
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    genes_stem = "_".join(found_genes)
    csv_path = output_dir / f"{genes_stem}_expression_by_region.csv"
    png_path = output_dir / f"{genes_stem}_expression_by_region.png"
    expr_long.to_csv(csv_path, index=False)

    plt.figure(figsize=(8, 5))
    sns.boxplot(
        data=expr_long,
        x=region_group_col,
        y="Expression",
        hue="Gene",
    )
    plt.title(f"Expression of {', '.join(found_genes)} across brain regions")
    plt.xlabel("Region")
    plt.ylabel("Expression")
    plt.tight_layout()
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close()

    return {
        "output_dir": str(output_dir),
        "plot_png_path": str(png_path),
        "expression_long_csv_path": str(csv_path),
        "h5ad_path": str(h5ad_path),
        "chemistry": chemistry,
        "chemistry_column": chemistry_col,
        "region_column": region_col,
        "region_group_column": region_group_col,
        "genes_requested": cleaned_genes,
        "genes_found": found_genes,
        "genes_missing": missing_genes,
        "n_cells": int(adata.n_obs),
    }


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Create helper dataset-analysis plots and tables."
    )
    parser.add_argument(
        "--task",
        default="age_counts",
        choices=["age_counts", "gene_region_boxplot"],
        help="Which dataset-analysis helper task to run.",
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
    parser.add_argument(
        "--genes",
        nargs="+",
        default=None,
        help="Gene symbols to plot for --task gene_region_boxplot.",
    )
    parser.add_argument(
        "--region-group-col",
        default=DEFAULT_REGION_GROUP_COL,
        help="Temporary/ad hoc column name to store collapsed region groups.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    chemistry = None if args.chemistry == "None" else args.chemistry
    if args.task == "age_counts":
        result = create_sample_count_by_age_plot(
            h5ad_path=args.h5ad_path,
            output_dir=args.output_dir,
            age_col=args.age_col,
            chemistry=chemistry,
            chemistry_col=args.chemistry_col,
            region=args.region,
            region_col=args.region_col,
        )
    else:
        result = create_gene_expression_boxplot_by_region(
            h5ad_path=args.h5ad_path,
            genes=args.genes or [],
            output_dir=args.output_dir,
            chemistry=chemistry,
            chemistry_col=args.chemistry_col,
            region_col=args.region_col,
            region_group_col=args.region_group_col,
        )

    print("Saved dataset analysis results:")
    for key, value in result.items():
        print(f"  {key}: {value}")
