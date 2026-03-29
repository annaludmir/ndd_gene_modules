import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist


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


def parse_go_genes(genes_text):
    """Parse the GO enrichment 'Genes' column into a clean gene list."""
    if pd.isna(genes_text):
        return []

    text = str(genes_text)
    for sep in [";", ",", "\n", "\t", "|", "/"]:
        text = text.replace(sep, " ")
    return [g.strip() for g in text.split() if g.strip()]


def clean_go_term_name(term):
    """Remove trailing GO accession text from a term label."""
    return re.sub(r"\s*\(GO:[^)]+\)\s*$", "", str(term).strip())


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


def load_go_term_file(go_term_file):
    """Load a GO-term enrichment table if provided."""
    if go_term_file is None:
        return None

    go_term_file = Path(go_term_file)
    if not go_term_file.exists():
        raise FileNotFoundError(f"GO term file not found: {go_term_file}")

    if go_term_file.suffix.lower() in {".tsv", ".txt"}:
        go_df = pd.read_csv(go_term_file, sep="\t")
    else:
        go_df = pd.read_csv(go_term_file)

    required_cols = {"Genes", "Term"}
    missing_cols = required_cols - set(go_df.columns)
    if missing_cols:
        raise ValueError(f"GO term file is missing required columns: {sorted(missing_cols)}")

    return go_df


def build_go_grouped_gene_order(leading_genes, go_df):
    """
    Assign genes to GO terms in GO-file order.
    Each gene is assigned once, to the first matching term.
    """
    if go_df is None:
        return [("All genes", list(leading_genes))], list(leading_genes)

    remaining = set(leading_genes)
    grouped = []

    for _, row in go_df.iterrows():
        term = clean_go_term_name(row["Term"])
        term_genes = parse_go_genes(row["Genes"])
        matched = [gene for gene in leading_genes if gene in remaining and gene in term_genes]
        if matched:
            grouped.append((term, matched))
            remaining -= set(matched)

    if remaining:
        grouped.append(("Unassigned", [gene for gene in leading_genes if gene in remaining]))

    ordered_genes = [gene for _, genes in grouped for gene in genes]
    return grouped, ordered_genes


def cluster_genes_by_pattern(zscore_expr):
    """Order genes by similarity of their pseudotime expression profiles."""
    genes = zscore_expr.columns.tolist()
    if len(genes) <= 2:
        return genes

    gene_matrix = zscore_expr[genes].T.to_numpy()
    if np.allclose(gene_matrix, gene_matrix[0]):
        return genes

    distances = pdist(gene_matrix, metric="euclidean")
    if np.allclose(distances, 0):
        return genes

    linkage_matrix = linkage(distances, method="average", optimal_ordering=True)
    return [genes[i] for i in leaves_list(linkage_matrix)]


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


def _add_pseudotime_arrow(ax):
    """Add a centered two-headed Early/Late arrow just below the heatmap."""
    y_pos = -0.065
    ax.text(0.12, y_pos, "Early", transform=ax.transAxes, ha="center", va="center", fontsize=12)
    ax.text(0.88, y_pos, "Late", transform=ax.transAxes, ha="center", va="center", fontsize=12)
    ax.annotate(
        "",
        xy=(0.78, y_pos),
        xytext=(0.22, y_pos),
        xycoords="axes fraction",
        textcoords="axes fraction",
        arrowprops={"arrowstyle": "<->", "lw": 1.6, "color": "black"},
    )


def plot_pseudotime_heatmap(
    zscore_expr,
    condition,
    output_path,
    gene_order,
    grouped_terms=None,
):
    """Plot a pseudotime heatmap, optionally with GO-term labels on the far left."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    plot_df = zscore_expr[gene_order].T
    show_group_labels = bool(grouped_terms) and any(term != "All genes" for term, _ in grouped_terms)

    fig_width = max(5.8, len(zscore_expr.index) * 0.55)
    fig_height = max(6, len(plot_df.index) * 0.34)
    left_margin = 0.36 if show_group_labels else 0.18

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=False)

    sns.heatmap(
        plot_df,
        ax=ax,
        cmap="Spectral_r",
        center=0,
        linewidths=0,
        cbar_kws={"label": "Gene expression z-score", "shrink": 0.9},
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

    if show_group_labels:
        boundary = 0
        for term, genes in grouped_terms:
            if not genes:
                continue

            start = boundary
            end = start + len(genes)
            center = (start + end) / 2

            if start > 0:
                ax.hlines(start, *ax.get_xlim(), colors="white", linewidth=2.2)

            ax.text(
                -1.15,
                center,
                term,
                ha="right",
                va="center",
                fontsize=10,
                fontweight="bold",
                clip_on=False,
            )
            boundary = end

    _add_pseudotime_arrow(ax)
    fig.subplots_adjust(left=left_margin, right=0.92, top=0.94, bottom=0.11)
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
    go_grouped_terms,
    clustered_genes,
):
    """Save supporting output tables for the pseudotime plots."""
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

    go_rows = []
    for term, genes in go_grouped_terms:
        for order, gene in enumerate(genes, start=1):
            go_rows.append({"Term": term, "gene": gene, "term_order": order})
    pd.DataFrame(go_rows).to_csv(output_dir / "go_term_gene_order.csv", index=False)

    pd.DataFrame(
        {"gene": clustered_genes, "expression_pattern_order": range(1, len(clustered_genes) + 1)}
    ).to_csv(output_dir / "expression_pattern_gene_order.csv", index=False)


def create_pseudotime_leading_gene_plot(
    gsea_summary_file,
    condition,
    h5ad_path,
    output_dir=DEFAULT_OUTPUT_DIR,
    subfolder_name=None,
    age_col="Age",
    sym_col="Gene",
    go_term_file=None,
):
    """
    Create pseudotime-style heatmaps for leading genes from a GSEA condition.
    """
    summary_info = load_summary_row(gsea_summary_file, condition)
    h5ad_path = Path(h5ad_path)
    if not h5ad_path.exists():
        raise FileNotFoundError(f"AnnData file not found: {h5ad_path}")

    output_dir = Path(output_dir)
    if subfolder_name:
        output_dir = output_dir / str(subfolder_name).strip()

    go_df = load_go_term_file(go_term_file)
    go_grouped_terms, go_ordered_genes = build_go_grouped_gene_order(
        summary_info["leading_genes"], go_df
    )

    adata = sc.read_h5ad(h5ad_path)
    mean_expr, zscore_expr, found_genes, missing_genes, n_cells = build_condition_time_matrices(
        adata=adata,
        filter_column=summary_info["column"],
        filter_value=summary_info["condition"],
        leading_genes=go_ordered_genes,
        age_col=age_col,
        sym_col=sym_col,
    )

    go_grouped_terms = [
        (term, [gene for gene in genes if gene in found_genes])
        for term, genes in go_grouped_terms
    ]
    go_grouped_terms = [(term, genes) for term, genes in go_grouped_terms if genes]
    go_gene_order = [gene for _, genes in go_grouped_terms for gene in genes]
    clustered_genes = cluster_genes_by_pattern(zscore_expr[go_gene_order])

    go_heatmap_path = output_dir / "leading_genes_pseudotime_heatmap_go_terms.png"
    pattern_heatmap_path = output_dir / "leading_genes_pseudotime_heatmap_expression_patterns.png"

    plot_pseudotime_heatmap(
        zscore_expr=zscore_expr,
        condition=summary_info["condition"],
        output_path=go_heatmap_path,
        gene_order=go_gene_order,
        grouped_terms=go_grouped_terms,
    )
    plot_pseudotime_heatmap(
        zscore_expr=zscore_expr,
        condition=summary_info["condition"],
        output_path=pattern_heatmap_path,
        gene_order=clustered_genes,
        grouped_terms=None,
    )

    save_outputs(
        output_dir=output_dir,
        mean_expr=mean_expr,
        zscore_expr=zscore_expr,
        found_genes=found_genes,
        missing_genes=missing_genes,
        condition=summary_info["condition"],
        filter_column=summary_info["column"],
        go_grouped_terms=go_grouped_terms,
        clustered_genes=clustered_genes,
    )

    return {
        "go_heatmap_path": str(go_heatmap_path),
        "expression_pattern_heatmap_path": str(pattern_heatmap_path),
        "output_dir": str(output_dir),
        "column": summary_info["column"],
        "condition": summary_info["condition"],
        "n_cells_used": int(n_cells),
        "n_found_genes": len(found_genes),
        "n_missing_genes": len(missing_genes),
    }


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Create pseudotime-style heatmaps for GSEA leading genes."
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
        "--go-term-file",
        default=None,
        help="Optional CSV/TSV enrichment file containing 'Genes' and 'Term' columns.",
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
        go_term_file=args.go_term_file,
    )
    print("Saved pseudotime leading-gene heatmap results:")
    for key, value in result.items():
        print(f"  {key}: {value}")
