import argparse
import re
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist


DEFAULT_OUTPUT_DIR = Path("results/time_analysis/early_late")
DEFAULT_DATA_PATH = Path("data/human_dev.h5ad")
STAGE_ORDER = ["Early", "Mid", "Late"]
DEFAULT_CHEMISTRY = "v3"
DEFAULT_CHEMISTRY_COL = "Chemistry"


def parse_gene_list(genes_value):
    """Accept a Python list or a delimited string and return clean gene symbols."""
    if isinstance(genes_value, (list, tuple, pd.Series, np.ndarray)):
        genes = [str(g).strip() for g in genes_value]
    else:
        text = str(genes_value).strip()
        for sep in [";", ",", "\n", "\t", "|"]:
            text = text.replace(sep, " ")
        genes = [g.strip() for g in text.split() if g.strip()]

    seen = set()
    cleaned = []
    for gene in genes:
        if gene not in seen:
            seen.add(gene)
            cleaned.append(gene)
    return cleaned


def filter_adata_by_chemistry(adata, chemistry=DEFAULT_CHEMISTRY, chemistry_col=DEFAULT_CHEMISTRY_COL):
    """Optionally filter AnnData to one chemistry value."""
    if chemistry is None:
        return adata
    if chemistry_col not in adata.obs.columns:
        raise KeyError(f"Chemistry column '{chemistry_col}' not found in adata.obs")

    before = adata.n_obs
    adata = adata[adata.obs[chemistry_col].astype(str) == str(chemistry)].copy()
    print(f"Filtered {chemistry_col} == '{chemistry}': {before} -> {adata.n_obs} cells")
    if adata.n_obs == 0:
        raise ValueError(
            f"No cells remained after filtering {chemistry_col} == '{chemistry}'."
        )
    return adata


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


def format_go_term_label(term, span_width):
    """Shorten and wrap a GO term label when there is enough room to show it."""
    term = str(term).strip()
    if not term:
        return None

    max_chars = max(10, int(span_width * 4.5))
    shortened = textwrap.shorten(term, width=max_chars, placeholder="...")

    if span_width < 3.0:
        return None
    if span_width < 5.5:
        return shortened if len(shortened) <= max(10, int(span_width * 3.0)) else None

    wrap_width = max(10, int(span_width * 2.6))
    wrapped = textwrap.wrap(shortened, width=wrap_width, max_lines=2, break_long_words=False)
    if not wrapped:
        return None
    if len(wrapped) == 2 and len(wrapped[1]) >= wrap_width:
        wrapped[1] = textwrap.shorten(wrapped[1], width=wrap_width, placeholder="...")
    return "\n".join(wrapped[:2])


def get_leading_genes_from_summary(summary_file, condition):
    """Extract leading genes from the summary row matching the requested condition."""
    summary_file = Path(summary_file)
    if not summary_file.exists():
        raise FileNotFoundError(f"GSEA summary file not found: {summary_file}")

    if summary_file.suffix.lower() in {".tsv", ".txt"}:
        summary_df = pd.read_csv(summary_file, sep="\t")
    else:
        summary_df = pd.read_csv(summary_file)

    required_cols = {"condition", "Lead_genes"}
    missing_cols = required_cols - set(summary_df.columns)
    if missing_cols:
        raise ValueError(
            f"GSEA summary file is missing required columns: {sorted(missing_cols)}"
        )

    matched = summary_df.loc[summary_df["condition"].astype(str) == str(condition)]
    if matched.empty:
        raise ValueError(
            f"Condition '{condition}' was not found in the summary file column 'condition'."
        )
    if len(matched) > 1:
        print(
            f"Warning: found {len(matched)} rows for condition '{condition}'. "
            "Using the first match."
        )

    lead_genes_value = matched.iloc[0]["Lead_genes"]
    leading_genes = parse_gene_list(str(lead_genes_value).replace(";", " "))
    if not leading_genes:
        raise ValueError(
            f"No leading genes were parsed from Lead_genes for condition '{condition}'."
        )
    return leading_genes


def classify_stage(age_value):
    """
    Map age into the requested developmental bins.

    Assumption to avoid overlap at week 9:
    - Early: 5.5 <= age <= 7.0
    - Mid:   7.0 < age < 9.0
    - Late:  9.0 <= age <= 14.0
    """
    if pd.isna(age_value):
        return np.nan

    age = float(age_value)

    if 5.5 <= age <= 7.0:
        return "Early"
    if 7.0 < age < 9.0:
        return "Mid"
    if 9.0 <= age <= 14.0:
        return "Late"
    return np.nan


def load_go_term_file(go_term_file):
    """Load the GO-term enrichment table."""
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
    """Order genes by similarity of their Early/Mid/Late expression profiles."""
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


def compute_stage_expression(
    adata,
    genes,
    age_col="Age",
    sym_col="Gene",
):
    """Compute mean expression per gene for Early/Mid/Late stages."""
    if age_col not in adata.obs.columns:
        raise KeyError(f"Column '{age_col}' not found in adata.obs")
    if sym_col not in adata.var.columns:
        raise KeyError(f"Column '{sym_col}' not found in adata.var")

    adata_f = adata.copy()
    adata_f.obs["time_stage"] = pd.to_numeric(adata_f.obs[age_col], errors="coerce").map(classify_stage)
    adata_f = adata_f[adata_f.obs["time_stage"].notna()].copy()

    sym2var = (
        pd.Series(adata_f.var_names.values, index=adata_f.var[sym_col].astype(str))
        .dropna()
        .to_dict()
    )

    found_genes = [gene for gene in genes if gene in sym2var]
    missing_genes = [gene for gene in genes if gene not in sym2var]

    if not found_genes:
        raise ValueError("None of the requested genes were found in adata.var[sym_col].")

    var_names = [sym2var[gene] for gene in found_genes]
    X = adata_f[:, var_names].X
    if sp.issparse(X):
        X = X.toarray()

    expr_df = pd.DataFrame(X, columns=found_genes)
    expr_df["time_stage"] = adata_f.obs["time_stage"].to_numpy()

    mean_expr = pd.DataFrame(index=STAGE_ORDER)
    for gene in found_genes:
        gene_expr = expr_df.loc[expr_df[gene] > 0, ["time_stage", gene]].copy()
        gene_means = gene_expr.groupby("time_stage")[gene].mean()
        mean_expr[gene] = gene_means.reindex(STAGE_ORDER)

    zscore_expr = mean_expr.copy()
    for gene in found_genes:
        values = mean_expr[gene].astype(float)
        std = values.std(ddof=0)
        if pd.isna(std) or std == 0:
            zscore_expr[gene] = 0.0
        else:
            zscore_expr[gene] = (values - values.mean()) / std

    return mean_expr, zscore_expr, found_genes, missing_genes


def plot_stage_heatmap(zscore_expr, gene_order, output_path, grouped_terms=None):
    """Plot a 3-row Early/Mid/Late heatmap with optional GO-term brackets below."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    plot_df = zscore_expr[gene_order]
    show_group_labels = bool(grouped_terms)

    if show_group_labels:
        fig, (ax, label_ax) = plt.subplots(
            2,
            1,
            figsize=(max(12, len(gene_order) * 0.34), 4.9),
            gridspec_kw={"height_ratios": [9, 2.7]},
            constrained_layout=False,
            sharex=True,
        )
    else:
        fig, ax = plt.subplots(
            figsize=(max(12, len(gene_order) * 0.34), 3.8),
            constrained_layout=False,
        )
        label_ax = None

    sns.heatmap(
        plot_df,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-1.5,
        vmax=1.5,
        linewidths=0.35,
        linecolor="#9f9f9f",
        cbar_kws={"label": "Mean z-score", "shrink": 0.7},
    )

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_yticklabels(STAGE_ORDER, rotation=0, fontsize=11)
    ax.set_xticklabels(gene_order, rotation=90, fontsize=9)
    ax.tick_params(axis="x", length=0)

    if show_group_labels:
        boundary = 0
        for term, genes in grouped_terms:
            if not genes:
                continue

            start = boundary
            end = start + len(genes)
            span_width = end - start
            center = (start + end) / 2
            term_label = format_go_term_label(term, span_width)

            if start > 0:
                ax.axvline(start, color="#6f6f6f", linewidth=1.0)

            label_ax.plot([start + 0.08, end - 0.08], [0.8, 0.8], color="#505050", linewidth=1.2)
            if term_label:
                label_ax.text(
                    center,
                    0.15,
                    term_label,
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    linespacing=1.0,
                )
            boundary = end

        label_ax.set_xlim(0, len(gene_order))
        label_ax.set_ylim(0, 1)
        label_ax.axis("off")
        fig.subplots_adjust(left=0.08, right=0.95, top=0.96, bottom=0.23, hspace=0.42)
    else:
        fig.subplots_adjust(left=0.08, right=0.95, top=0.96, bottom=0.32)

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_outputs(
    output_dir,
    grouped_terms,
    clustered_genes,
    mean_expr,
    zscore_expr,
    found_genes,
    missing_genes,
):
    """Write supporting tables for the heatmap output."""
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_order_rows = []
    for term, genes in grouped_terms:
        for order, gene in enumerate(genes, start=1):
            gene_order_rows.append({"Term": term, "gene": gene, "term_order": order})

    pd.DataFrame(gene_order_rows).to_csv(output_dir / "go_term_gene_order.csv", index=False)
    pd.DataFrame(
        {"gene": clustered_genes, "expression_pattern_order": range(1, len(clustered_genes) + 1)}
    ).to_csv(output_dir / "expression_pattern_gene_order.csv", index=False)
    mean_expr.T.reset_index(names="gene").to_csv(
        output_dir / "leading_genes_mean_expression_by_stage.csv",
        index=False,
    )
    zscore_expr.T.reset_index(names="gene").to_csv(
        output_dir / "leading_genes_stage_zscore.csv",
        index=False,
    )
    pd.DataFrame({"gene": found_genes}).to_csv(output_dir / "genes_found_in_data.csv", index=False)
    pd.DataFrame({"gene": missing_genes}).to_csv(output_dir / "genes_missing_from_data.csv", index=False)


def create_early_late_go_heatmap(
    gsea_summary_file,
    condition,
    go_term_file,
    h5ad_path=DEFAULT_DATA_PATH,
    output_dir=DEFAULT_OUTPUT_DIR,
    subfolder_name=None,
    age_col="Age",
    sym_col="Gene",
    chemistry=DEFAULT_CHEMISTRY,
    chemistry_col=DEFAULT_CHEMISTRY_COL,
):
    """
    Build Early/Mid/Late heatmaps for a list of leading genes.
    """
    leading_genes = get_leading_genes_from_summary(gsea_summary_file, condition)
    h5ad_path = Path(h5ad_path)
    output_dir = Path(output_dir)
    if subfolder_name:
        output_dir = output_dir / str(subfolder_name).strip()

    if not h5ad_path.exists():
        raise FileNotFoundError(f"AnnData file not found: {h5ad_path}")

    go_df = load_go_term_file(go_term_file)
    grouped_terms, ordered_genes = build_go_grouped_gene_order(leading_genes, go_df)

    adata = sc.read_h5ad(h5ad_path)
    adata = filter_adata_by_chemistry(
        adata,
        chemistry=chemistry,
        chemistry_col=chemistry_col,
    )
    mean_expr, zscore_expr, found_genes, missing_genes = compute_stage_expression(
        adata=adata,
        genes=ordered_genes,
        age_col=age_col,
        sym_col=sym_col,
    )

    grouped_terms = [
        (term, [gene for gene in genes if gene in found_genes])
        for term, genes in grouped_terms
    ]
    grouped_terms = [(term, genes) for term, genes in grouped_terms if genes]
    go_gene_order = [gene for _, genes in grouped_terms for gene in genes]
    clustered_genes = cluster_genes_by_pattern(zscore_expr[go_gene_order])

    go_heatmap_path = output_dir / "leading_genes_early_mid_late_heatmap_go_terms.png"
    pattern_heatmap_path = output_dir / "leading_genes_early_mid_late_heatmap_expression_patterns.png"

    plot_stage_heatmap(
        zscore_expr=zscore_expr,
        gene_order=go_gene_order,
        output_path=go_heatmap_path,
        grouped_terms=grouped_terms,
    )
    plot_stage_heatmap(
        zscore_expr=zscore_expr,
        gene_order=clustered_genes,
        output_path=pattern_heatmap_path,
        grouped_terms=None,
    )

    save_outputs(
        output_dir=output_dir,
        grouped_terms=grouped_terms,
        clustered_genes=clustered_genes,
        mean_expr=mean_expr,
        zscore_expr=zscore_expr,
        found_genes=found_genes,
        missing_genes=missing_genes,
    )

    return {
        "go_heatmap_path": str(go_heatmap_path),
        "expression_pattern_heatmap_path": str(pattern_heatmap_path),
        "output_dir": str(output_dir),
        "condition": str(condition),
        "chemistry": chemistry,
        "chemistry_column": chemistry_col,
        "n_found_genes": len(found_genes),
        "n_missing_genes": len(missing_genes),
    }


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Create Early/Mid/Late heatmaps using leading genes from a GSEA summary file."
    )
    parser.add_argument(
        "--gsea-summary-file",
        required=True,
        help="GSEA summary CSV/TSV containing 'condition' and 'Lead_genes' columns.",
    )
    parser.add_argument(
        "--condition",
        required=True,
        help="Condition to select from the GSEA summary file.",
    )
    parser.add_argument(
        "--go-term-file",
        required=True,
        help="CSV/TSV enrichment file containing 'Genes' and 'Term' columns.",
    )
    parser.add_argument(
        "--h5ad-path",
        default=str(DEFAULT_DATA_PATH),
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
    parser.add_argument(
        "--chemistry",
        default=DEFAULT_CHEMISTRY,
        help="Chemistry value to keep from adata.obs before analysis. Use 'None' to disable.",
    )
    parser.add_argument(
        "--chemistry-col",
        default=DEFAULT_CHEMISTRY_COL,
        help="Column in adata.obs containing chemistry labels.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    result = create_early_late_go_heatmap(
        gsea_summary_file=args.gsea_summary_file,
        condition=args.condition,
        go_term_file=args.go_term_file,
        h5ad_path=args.h5ad_path,
        output_dir=args.output_dir,
        subfolder_name=args.subfolder_name,
        age_col=args.age_col,
        sym_col=args.sym_col,
        chemistry=None if args.chemistry == "None" else args.chemistry,
        chemistry_col=args.chemistry_col,
    )
    print("Saved Early/Mid/Late GO heatmap results:")
    for key, value in result.items():
        print(f"  {key}: {value}")
