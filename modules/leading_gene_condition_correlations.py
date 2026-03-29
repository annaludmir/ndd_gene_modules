import argparse
import itertools
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


DEFAULT_OUTPUT_DIR = Path("results/correlations")
DEFAULT_TOP_N = 25


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


def sanitize_name(text):
    """Convert free text to a filesystem-friendly stem."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(text).strip())
    cleaned = cleaned.strip("._")
    return cleaned or "correlations"


def discover_summary_files(input_folder):
    """Find summary-like tables that contain leading genes."""
    input_folder = Path(input_folder)
    if not input_folder.exists():
        raise FileNotFoundError(f"Input folder not found: {input_folder}")
    if not input_folder.is_dir():
        raise NotADirectoryError(f"Input path is not a directory: {input_folder}")

    candidate_files = sorted(
        path
        for path in input_folder.rglob("*")
        if path.is_file() and path.suffix.lower() in {".csv", ".tsv", ".txt"}
    )

    summary_files = []
    for path in candidate_files:
        try:
            if path.suffix.lower() in {".tsv", ".txt"}:
                preview = pd.read_csv(path, sep="\t", nrows=5)
            else:
                preview = pd.read_csv(path, nrows=5)
        except Exception:
            continue

        required_cols = {"column", "condition", "Lead_genes"}
        if required_cols.issubset(preview.columns):
            summary_files.append(path)

    if not summary_files:
        raise ValueError(
            f"No summary files with columns {sorted(required_cols)} were found under {input_folder}."
        )
    return summary_files


def load_condition_entries(summary_files):
    """Load all valid condition rows from the discovered summary files."""
    entries = []

    for summary_path in summary_files:
        if summary_path.suffix.lower() in {".tsv", ".txt"}:
            df = pd.read_csv(summary_path, sep="\t")
        else:
            df = pd.read_csv(summary_path)

        for _, row in df.iterrows():
            genes = parse_gene_list(row["Lead_genes"])
            if not genes:
                continue

            column_name = str(row["column"]).strip()
            condition_value = str(row["condition"]).strip()
            label = f"{column_name}: {condition_value}"

            entries.append(
                {
                    "column": column_name,
                    "condition": condition_value,
                    "condition_label": label,
                    "summary_file": str(summary_path.resolve()),
                    "leading_genes": genes,
                    "leading_gene_set": set(genes),
                }
            )

    if not entries:
        raise ValueError("Summary files were found, but none contained parseable leading genes.")
    return entries


def compute_jaccard_score(genes_a, genes_b):
    """Measure overlap between two leading-gene sets."""
    union = genes_a | genes_b
    if not union:
        return 0.0
    return len(genes_a & genes_b) / len(union)


def build_pairwise_correlation_table(entries):
    """Compute all unique condition-pair correlations."""
    rows = []

    for left_entry, right_entry in itertools.combinations_with_replacement(entries, 2):
        overlap_genes = sorted(left_entry["leading_gene_set"] & right_entry["leading_gene_set"])
        score = compute_jaccard_score(
            left_entry["leading_gene_set"], right_entry["leading_gene_set"]
        )

        rows.append(
            {
                "left_column": left_entry["column"],
                "left_condition": left_entry["condition"],
                "left_condition_label": left_entry["condition_label"],
                "left_summary_file": left_entry["summary_file"],
                "left_n_leading_genes": len(left_entry["leading_genes"]),
                "right_column": right_entry["column"],
                "right_condition": right_entry["condition"],
                "right_condition_label": right_entry["condition_label"],
                "right_summary_file": right_entry["summary_file"],
                "right_n_leading_genes": len(right_entry["leading_genes"]),
                "correlation_score": score,
                "n_overlap_genes": len(overlap_genes),
                "overlap_genes": ";".join(overlap_genes),
            }
        )

    pairwise_df = pd.DataFrame(rows).sort_values(
        by=["correlation_score", "n_overlap_genes", "left_condition_label", "right_condition_label"],
        ascending=[False, False, True, True],
    )
    return pairwise_df.reset_index(drop=True)


def build_full_correlation_matrix(entries):
    """Build a symmetric condition-by-condition matrix."""
    labels = []
    label_counts = {}
    for entry in entries:
        base = entry["condition_label"]
        count = label_counts.get(base, 0) + 1
        label_counts[base] = count
        label = base if count == 1 else f"{base} [{count}]"
        entry["matrix_label"] = label
        labels.append(label)

    matrix = pd.DataFrame(np.eye(len(entries)), index=labels, columns=labels, dtype=float)

    for i, left_entry in enumerate(entries):
        for j in range(i + 1, len(entries)):
            right_entry = entries[j]
            score = compute_jaccard_score(
                left_entry["leading_gene_set"], right_entry["leading_gene_set"]
            )
            matrix.iloc[i, j] = score
            matrix.iloc[j, i] = score

    return matrix


def select_top_conditions(matrix, top_n):
    """Pick conditions participating in the strongest non-self correlations."""
    if matrix.shape[0] <= top_n:
        return matrix

    scores = []
    for idx in matrix.index:
        row = matrix.loc[idx].drop(index=idx, errors="ignore")
        scores.append((idx, float(row.max()) if not row.empty else 0.0))

    top_labels = [
        label
        for label, _ in sorted(scores, key=lambda item: (-item[1], item[0]))[:top_n]
    ]
    return matrix.loc[top_labels, top_labels]


def plot_top_correlation_heatmap(matrix, output_path):
    """Save a heatmap of the strongest condition correlations."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig_width = max(8, min(18, matrix.shape[1] * 0.48 + 4))
    fig_height = max(7, min(18, matrix.shape[0] * 0.38 + 3))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=False)

    sns.heatmap(
        matrix,
        ax=ax,
        cmap="mako",
        vmin=0,
        vmax=1,
        square=False,
        linewidths=0.35,
        linecolor="#d0d0d0",
        cbar_kws={"label": "Leading-gene correlation (Jaccard)", "shrink": 0.8},
    )

    ax.set_title("Highest Correlations of Leading Genes", fontsize=15, pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)

    fig.subplots_adjust(left=0.33, bottom=0.28, right=0.93, top=0.92)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def create_leading_gene_condition_correlations(
    input_folder,
    output_dir=DEFAULT_OUTPUT_DIR,
    top_n=DEFAULT_TOP_N,
):
    """Build pairwise leading-gene correlations across summary files in a folder tree."""
    input_folder = Path(input_folder)
    output_dir = Path(output_dir) / sanitize_name(input_folder.resolve().name)
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_files = discover_summary_files(input_folder)
    entries = load_condition_entries(summary_files)

    pairwise_df = build_pairwise_correlation_table(entries)
    full_matrix = build_full_correlation_matrix(entries)
    top_matrix = select_top_conditions(full_matrix, top_n=top_n)

    csv_path = output_dir / "leading_gene_condition_correlations.csv"
    matrix_csv_path = output_dir / "leading_gene_condition_correlation_matrix.csv"
    heatmap_path = output_dir / "leading_gene_condition_correlation_heatmap.png"

    pairwise_df.to_csv(csv_path, index=False)
    full_matrix.to_csv(matrix_csv_path, index=True)
    plot_top_correlation_heatmap(top_matrix, heatmap_path)

    return {
        "output_dir": str(output_dir),
        "pairwise_csv_path": str(csv_path),
        "matrix_csv_path": str(matrix_csv_path),
        "heatmap_path": str(heatmap_path),
        "n_summary_files": len(summary_files),
        "n_conditions": len(entries),
        "heatmap_n_conditions": int(top_matrix.shape[0]),
    }


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Compute pairwise leading-gene correlations across summary files in a folder tree."
    )
    parser.add_argument(
        "input_folder",
        help="Folder to scan recursively for summary CSV/TSV/TXT files with columns 'column', 'condition', and 'Lead_genes'.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Base folder where correlation outputs will be saved.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=DEFAULT_TOP_N,
        help="Number of top-correlated conditions to include in the heatmap.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    result = create_leading_gene_condition_correlations(
        input_folder=args.input_folder,
        output_dir=args.output_dir,
        top_n=args.top_n,
    )
    print("Saved leading-gene correlation results:")
    for key, value in result.items():
        print(f"  {key}: {value}")
