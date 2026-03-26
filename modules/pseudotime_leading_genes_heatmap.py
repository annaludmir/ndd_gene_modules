import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns


DEFAULT_OUTPUT_DIR = Path("results/time_analysis/pseudotime")


def parse_gene_list(genes_value):
    """Parse Lead_genes text into a unique, ordered list of gene symbols."""
    if pd.isna(genes_value):
        return []

    text = str(genes_value)
    for sep in [";", ",", "\n", "\t", "|", "/"]:
        text = text.replace(sep, " ")

    seen = set()
    genes = []
    for gene in text.split():
        gene = gene.strip()
        if gene and gene not in seen:
            seen.add(gene)
            genes.append(gene)
    return genes


def load_summary_row(summary_file, condition):
    """Load the GSEA summary row for the requested condition."""
    summary_file = Path(summary_file)
    if not summary_file.exists():
        raise FileNotFoundError(f"GSEA summary file not found: {summary_file}")

    if summary_file.suffix.lower() in {".tsv", ".txt"}:
        summary_df = pd.read_csv(summary_file, sep="\t")
    else:
        summary_df = pd.read_csv(summary_file)

    required_cols = {"column", "condition", "Lead_genes"}
    missing_cols = required_cols - set(summary_df.columns)
    if missing_cols:
        raise ValueError(
            f"GSEA summary file is missing required columns: {sorted(missing_cols)}"
        )

    matched = summary_df.loc[summary_df["condition"].astype(str) == str(condition)].copy()
    if matched.empty:
        raise ValueError(f"Condition '{condition}' was not found in the summary file.")

    if len(matched) > 1:
        if "NES" in matched.columns:
            matched["NES_abs"] = pd.to_numeric(matched["NES"], errors="coerce").abs()
            matched = matched.sort_values("NES_abs", ascending=False)
        print(
            f"Warning: found {len(matched)} rows for condition '{condition}'. "
            "Using the first match."
        )

    row = matched.iloc[0]
    genes = parse_gene_list(row["Lead_genes"])
    if not genes:
        raise ValueError(f"No leading genes found for condition '{condition}'.")

    return {
        "column": str(row["column"]),
        "condition": str(row["condition"]),
        "leading_genes": genes,
    }


def build_condition_time_matrices(
    adata,
    filter_column,
    filter_value,
    leading_genes,
    age_col="Age",
    sym_col="Gene",
):
    """Compute mean expression per age for the requested condition and genes."""
    if filter_column not in adata.obs.columns:
        raise KeyError(f"Column '{filter_column}' not found in adata.obs")
    if age_col not in adata.obs.columns:
        raise KeyError(f"Column '{age_col}' not found in adata.obs")
    if sym_col not in adata.var.columns:
        raise KeyError(f"Column '{sym_col}' not found in adata.var")

    obs_condition = adata.obs[filter_column].astype(str)
    cond_mask = obs_condition == str(filter_value)
    adata_f = adata[cond_mask].copy()
    if adata_f.n_obs == 0:
        raise ValueError(
            f"No cells found for condition '{filter_value}' in adata.obs['{filter_column}']."
        )

    adata_f.obs["_age_numeric"] = pd.to_numeric(adata_f.obs[age_col], errors="coerce")
    adata_f = adata_f[adata_f.obs["_age_numeric"].notna()].copy()
    if adata_f.n_obs == 0:
        raise ValueError("No cells with numeric age values were found after filtering.")

    sym2var = (
        pd.Series(adata_f.var_names.values, index=adata_f.var[sym_col].astype(str))
        .dropna()
        .to_dict()
    )

    found_genes = [gene for gene in leading_genes if gene in sym2var]
    missing_genes = [gene for gene in leading_genes if gene not in sym2var]
    if not found_genes:
        raise ValueError("None of the leading genes were found in adata.var[sym_col].")

    var_names = [sym2var[gene] for gene in found_genes]
    X = adata_f[:, var_names].X
    if sp.issparse(X):
        X = X.toarray()

    expr_df = pd.DataFrame(X, columns=found_genes)
    expr_df["_age_numeric"] = adata_f.obs["_age_numeric"].to_numpy()

    mean_expr = expr_df.groupby("_age_numeric")[found_genes].mean().sort_index()
    mean_expr.index.name = "Age"

    zscore_expr = mean_expr.copy()
    for gene in found_genes:
        values = mean_expr[gene].astype(float)
        std = values.std(ddof=0)
        if pd.isna(std) or std == 0:
            zscore_expr[gene] = 0.0
        else:
            zscore_expr[gene] = (values - values.mean()) / std

    return mean_expr, zscore_expr, found_genes, missing_genes, adata_f.n_obs


def plot_pseudotime_heatmap(zscore_expr, condition, output_path):
    """Plot genes x age heatmap resembling a pseudotime progression."""
    plot_df = zscore_expr.T

    fig, ax = plt.subplots(
        figsize=(max(5.5, len(zscore_expr.index) * 0.55), max(6, len(plot_df.index) * 0.34)),
        constrained_layout=True,
    )

    sns.heatmap(
        plot_df,
        ax=ax,
        cmap="Spectral_r",
        center=0,
        linewidths=0,
        cbar_kws={"label": "Gene expression z-score"},
    )

    ax.set_title(str(condition), fontsize=16, pad=10)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_yticklabels(plot_df.index.tolist(), rotation=0, fontsize=9, fontweight="bold")

    age_labels = [f"{age:g}" for age in plot_df.columns.to_list()]
    if len(age_labels) <= 12:
        ax.set_xticklabels(age_labels, rotation=0, fontsize=9)
    else:
        tick_idx = np.linspace(0, len(age_labels) - 1, num=min(8, len(age_labels)), dtype=int)
        ax.set_xticks(tick_idx + 0.5)
        ax.set_xticklabels([age_labels[i] for i in tick_idx], rotation=0, fontsize=9)

    ax.text(0.0, -0.12, "Early", transform=ax.transAxes, ha="left", va="center", fontsize=14)
    ax.text(1.0, -0.12, "Late", transform=ax.transAxes, ha="right", va="center", fontsize=14)
    ax.annotate(
        "",
        xy=(0.98, -0.12),
        xytext=(0.02, -0.12),
        xycoords="axes fraction",
        textcoords="axes fraction",
        arrowprops={"arrowstyle": "->", "lw": 1.8, "color": "black"},
    )

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_outputs(
    output_dir,
    mean_expr,
    zscore_expr,
    found_genes,
    missing_genes,
    condition,
    filter_column,
):
    """Save supporting output tables for the pseudotime plot."""
    output_dir.mkdir(parents=True, exist_ok=True)

    mean_expr.T.reset_index(names="gene").to_csv(
        output_dir / "leading_genes_mean_expression_by_age.csv",
        index=False,
    )
    zscore_expr.T.reset_index(names="gene").to_csv(
        output_dir / "leading_genes_zscore_by_age.csv",
        index=False,
    )
    pd.DataFrame({"gene": found_genes}).to_csv(output_dir / "genes_found_in_data.csv", index=False)
    pd.DataFrame({"gene": missing_genes}).to_csv(output_dir / "genes_missing_from_data.csv", index=False)
    pd.DataFrame(
        [{"column": filter_column, "condition": condition, "n_genes": len(found_genes)}]
    ).to_csv(output_dir / "selection_metadata.csv", index=False)


def create_pseudotime_leading_gene_plot(
    gsea_summary_file,
    condition,
    h5ad_path,
    output_dir=DEFAULT_OUTPUT_DIR,
    subfolder_name=None,
    age_col="Age",
    sym_col="Gene",
):
    """
    Create a pseudotime-style heatmap for leading genes from a GSEA condition.
    """
    summary_info = load_summary_row(gsea_summary_file, condition)
    h5ad_path = Path(h5ad_path)
    if not h5ad_path.exists():
        raise FileNotFoundError(f"AnnData file not found: {h5ad_path}")

    output_dir = Path(output_dir)
    if subfolder_name:
        output_dir = output_dir / str(subfolder_name).strip()

    adata = sc.read_h5ad(h5ad_path)
    mean_expr, zscore_expr, found_genes, missing_genes, n_cells = build_condition_time_matrices(
        adata=adata,
        filter_column=summary_info["column"],
        filter_value=summary_info["condition"],
        leading_genes=summary_info["leading_genes"],
        age_col=age_col,
        sym_col=sym_col,
    )

    heatmap_path = output_dir / "leading_genes_pseudotime_heatmap.png"
    plot_pseudotime_heatmap(zscore_expr, summary_info["condition"], heatmap_path)
    save_outputs(
        output_dir=output_dir,
        mean_expr=mean_expr,
        zscore_expr=zscore_expr,
        found_genes=found_genes,
        missing_genes=missing_genes,
        condition=summary_info["condition"],
        filter_column=summary_info["column"],
    )

    return {
        "heatmap_path": str(heatmap_path),
        "output_dir": str(output_dir),
        "column": summary_info["column"],
        "condition": summary_info["condition"],
        "n_cells_used": int(n_cells),
        "n_found_genes": len(found_genes),
        "n_missing_genes": len(missing_genes),
    }


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Create a pseudotime-style heatmap for GSEA leading genes."
    )
    parser.add_argument(
        "--gsea-summary-file",
        required=True,
        help="GSEA summary CSV/TSV containing 'column', 'condition', and 'Lead_genes' columns.",
    )
    parser.add_argument(
        "--condition",
        required=True,
        help="Condition to select from the GSEA summary file.",
    )
    parser.add_argument(
        "--h5ad-path",
        required=True,
        help="Path to the AnnData file used for expression summaries.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Base folder where the heatmap and summary tables will be saved.",
    )
    parser.add_argument(
        "--subfolder-name",
        default=None,
        help="Optional subfolder name created under the output directory for this run.",
    )
    parser.add_argument("--age-col", default="Age", help="Age column in adata.obs.")
    parser.add_argument("--sym-col", default="Gene", help="Gene symbol column in adata.var.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    result = create_pseudotime_leading_gene_plot(
        gsea_summary_file=args.gsea_summary_file,
        condition=args.condition,
        h5ad_path=args.h5ad_path,
        output_dir=args.output_dir,
        subfolder_name=args.subfolder_name,
        age_col=args.age_col,
        sym_col=args.sym_col,
    )
    print("Saved pseudotime leading-gene heatmap results:")
    for key, value in result.items():
        print(f"  {key}: {value}")
