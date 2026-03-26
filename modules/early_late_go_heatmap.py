import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns


DEFAULT_OUTPUT_DIR = Path("results/time_analysis/early_late")
DEFAULT_DATA_PATH = Path("data/human_dev.h5ad")
STAGE_ORDER = ["Early", "Mid", "Late"]


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


def parse_go_genes(genes_text):
    """Parse the GO enrichment 'Genes' column into a clean gene list."""
    if pd.isna(genes_text):
        return []

    text = str(genes_text)
    for sep in [";", ",", "\n", "\t", "|", "/"]:
        text = text.replace(sep, " ")
    return [g.strip() for g in text.split() if g.strip()]


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


def build_gene_order(leading_genes, go_df):
    """
    Assign genes to GO terms in GO-file order.
    Each gene is assigned once, to the first matching term.
    """
    remaining = set(leading_genes)
    grouped = []

    for _, row in go_df.iterrows():
        term = str(row["Term"]).strip()
        term_genes = parse_go_genes(row["Genes"])
        matched = [gene for gene in leading_genes if gene in remaining and gene in term_genes]
        if matched:
            grouped.append((term, matched))
            remaining -= set(matched)

    if remaining:
        grouped.append(("Unassigned", [gene for gene in leading_genes if gene in remaining]))

    ordered_genes = [gene for _, genes in grouped for gene in genes]
    return grouped, ordered_genes


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

    mean_expr = expr_df.groupby("time_stage")[found_genes].mean().reindex(STAGE_ORDER)

    # Gene-wise z-score across the 3 stages to match the reference style.
    zscore_expr = mean_expr.copy()
    for gene in found_genes:
        values = mean_expr[gene].astype(float)
        std = values.std(ddof=0)
        if pd.isna(std) or std == 0:
            zscore_expr[gene] = 0.0
        else:
            zscore_expr[gene] = (values - values.mean()) / std

    return mean_expr, zscore_expr, found_genes, missing_genes


def plot_stage_heatmap(zscore_expr, grouped_terms, output_path):
    """Plot a 3-row early/mid/late heatmap with GO-term group labels."""
    ordered_genes = [gene for _, genes in grouped_terms for gene in genes]
    plot_df = zscore_expr[ordered_genes]

    fig, (ax, label_ax) = plt.subplots(
        2,
        1,
        figsize=(max(12, len(ordered_genes) * 0.32), 4.8),
        gridspec_kw={"height_ratios": [10, 2]},
        constrained_layout=True,
    )

    sns.heatmap(
        plot_df,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-1.5,
        vmax=1.5,
        linewidths=0.5,
        linecolor="#d0d0d0",
        cbar_kws={"label": "Mean z-score"},
    )

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_yticklabels(STAGE_ORDER, rotation=0, fontsize=11)
    ax.set_xticklabels(ordered_genes, rotation=90, fontsize=9)

    boundary = 0
    for term, genes in grouped_terms:
        start = boundary
        end = start + len(genes)
        center = (start + end) / 2

        if start > 0:
            ax.axvline(start, color="black", linewidth=1.2)

        label_ax.plot([start, end], [0.8, 0.8], color="black", linewidth=1.0)
        label_ax.text(center, 0.15, term, ha="center", va="bottom", fontsize=10)
        boundary = end

    label_ax.set_xlim(0, len(ordered_genes))
    label_ax.set_ylim(0, 1)
    label_ax.axis("off")

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_outputs(
    output_dir,
    grouped_terms,
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
):
    """
    Build an Early/Mid/Late heatmap for a list of leading genes grouped by GO term.

    Parameters
    ----------
    gsea_summary_file
        GSEA summary CSV/TSV containing 'condition' and 'Lead_genes' columns.
    condition
        Value to match in the summary file 'condition' column.
    go_term_file
        CSV/TSV file containing at least 'Genes' and 'Term' columns.
    h5ad_path
        AnnData file used to compute expression.
    output_dir
        Base output folder. Defaults to results/time_analysis/early_late.
    subfolder_name
        Optional subfolder created under output_dir for this run.
    """
    leading_genes = get_leading_genes_from_summary(gsea_summary_file, condition)
    gsea_summary_file = Path(gsea_summary_file)
    go_term_file = Path(go_term_file)
    h5ad_path = Path(h5ad_path)
    output_dir = Path(output_dir)
    if subfolder_name:
        output_dir = output_dir / str(subfolder_name).strip()

    if not go_term_file.exists():
        raise FileNotFoundError(f"GO term file not found: {go_term_file}")
    if not h5ad_path.exists():
        raise FileNotFoundError(f"AnnData file not found: {h5ad_path}")

    if go_term_file.suffix.lower() in {".tsv", ".txt"}:
        go_df = pd.read_csv(go_term_file, sep="\t")
    else:
        go_df = pd.read_csv(go_term_file)

    required_cols = {"Genes", "Term"}
    missing_cols = required_cols - set(go_df.columns)
    if missing_cols:
        raise ValueError(f"GO term file is missing required columns: {sorted(missing_cols)}")

    grouped_terms, ordered_genes = build_gene_order(leading_genes, go_df)
    adata = sc.read_h5ad(h5ad_path)
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

    heatmap_path = output_dir / "leading_genes_early_mid_late_heatmap.png"
    plot_stage_heatmap(zscore_expr, grouped_terms, heatmap_path)
    save_outputs(output_dir, grouped_terms, mean_expr, zscore_expr, found_genes, missing_genes)

    return {
        "heatmap_path": str(heatmap_path),
        "output_dir": str(output_dir),
        "condition": str(condition),
        "n_found_genes": len(found_genes),
        "n_missing_genes": len(missing_genes),
    }


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Create an Early/Mid/Late heatmap using leading genes from a GSEA summary file."
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
    )
    print("Saved Early/Mid/Late GO heatmap results:")
    for key, value in result.items():
        print(f"  {key}: {value}")
