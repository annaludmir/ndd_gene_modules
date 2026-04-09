import argparse
import itertools
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


DEFAULT_OUTPUT_DIR = Path("results/correlations")
DEFAULT_TOP_N = 25
DEFAULT_MIN_CORRELATION = 0.5
DEFAULT_HEATMAP_PVALUE_THRESHOLD = 0.3


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


def parse_run_name_from_config(config_path):
    """Extract run_name from a metadata YAML file with a simple line parser."""
    try:
        for line in config_path.read_text(encoding="utf-8").splitlines():
            match = re.match(r"^\s*run_name\s*:\s*(.+?)\s*$", line)
            if not match:
                continue

            value = match.group(1).strip()
            if value and value[0] == value[-1] and value[0] in {'"', "'"}:
                value = value[1:-1]
            return value.strip()
    except Exception:
        return None
    return None


def find_run_name_for_summary(summary_path):
    """Locate the nearest metadata config for a summary file and read its run name."""
    summary_path = Path(summary_path).resolve()

    for parent in [summary_path.parent, *summary_path.parents]:
        metadata_dir = parent / "metadata"
        if not metadata_dir.is_dir():
            continue

        preferred_paths = [metadata_dir / "config_used.yaml"]
        preferred_paths.extend(sorted(metadata_dir.glob("*_config.yaml")))

        for config_path in preferred_paths:
            if not config_path.exists():
                continue
            run_name = parse_run_name_from_config(config_path)
            if run_name:
                return run_name, str(config_path.resolve())

    return summary_path.parent.parent.name, ""


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
    seen_keys = set()

    for summary_path in summary_files:
        run_name, metadata_config_path = find_run_name_for_summary(summary_path)

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
            label = f"{run_name} | {column_name}: {condition_value}"
            entry_key = (run_name, column_name, condition_value)
            if entry_key in seen_keys:
                continue
            seen_keys.add(entry_key)

            entries.append(
                {
                    "run_name": run_name,
                    "column": column_name,
                    "condition": condition_value,
                    "condition_label": label,
                    "summary_file": str(summary_path.resolve()),
                    "metadata_config_path": metadata_config_path,
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


def log_comb(n, k):
    """Compute log(n choose k) safely."""
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def compute_overlap_pvalue(universe_size, left_size, right_size, overlap_size):
    """Hypergeometric one-sided p-value for observing at least this overlap."""
    if universe_size <= 0:
        return np.nan

    min_overlap = max(0, left_size + right_size - universe_size)
    max_overlap = min(left_size, right_size)
    if overlap_size < min_overlap or overlap_size > max_overlap:
        return np.nan

    log_denom = log_comb(universe_size, right_size)
    tail_probability = 0.0
    for k in range(overlap_size, max_overlap + 1):
        log_num = log_comb(left_size, k) + log_comb(universe_size - left_size, right_size - k)
        tail_probability += math.exp(log_num - log_denom)

    return min(1.0, tail_probability)


def build_pairwise_correlation_table(entries):
    """Compute all unique condition-pair correlations."""
    rows = []
    universe_genes = set()
    for entry in entries:
        universe_genes.update(entry["leading_gene_set"])
    universe_size = len(universe_genes)

    for left_entry, right_entry in itertools.combinations_with_replacement(entries, 2):
        overlap_genes = sorted(left_entry["leading_gene_set"] & right_entry["leading_gene_set"])
        overlap_size = len(overlap_genes)
        score = compute_jaccard_score(
            left_entry["leading_gene_set"], right_entry["leading_gene_set"]
        )
        overlap_pvalue = compute_overlap_pvalue(
            universe_size=universe_size,
            left_size=len(left_entry["leading_gene_set"]),
            right_size=len(right_entry["leading_gene_set"]),
            overlap_size=overlap_size,
        )

        rows.append(
            {
                "left_run_name": left_entry["run_name"],
                "left_column": left_entry["column"],
                "left_condition": left_entry["condition"],
                "left_condition_label": left_entry["condition_label"],
                "left_summary_file": left_entry["summary_file"],
                "left_metadata_config_path": left_entry["metadata_config_path"],
                "left_n_leading_genes": len(left_entry["leading_genes"]),
                "right_run_name": right_entry["run_name"],
                "right_column": right_entry["column"],
                "right_condition": right_entry["condition"],
                "right_condition_label": right_entry["condition_label"],
                "right_summary_file": right_entry["summary_file"],
                "right_metadata_config_path": right_entry["metadata_config_path"],
                "right_n_leading_genes": len(right_entry["leading_genes"]),
                "correlation_score": score,
                "overlap_pvalue": overlap_pvalue,
                "n_overlap_genes": overlap_size,
                "overlap_genes": ";".join(overlap_genes),
            }
        )

    pairwise_df = pd.DataFrame(rows).sort_values(
        by=["correlation_score", "n_overlap_genes", "left_condition_label", "right_condition_label"],
        ascending=[False, False, True, True],
    )
    return pairwise_df.reset_index(drop=True)


def build_filtered_correlation_export(pairwise_df, min_correlation=DEFAULT_MIN_CORRELATION):
    """Build a simple condition-pair export for strong, non-self correlations."""
    filtered_df = pairwise_df.loc[
        (pairwise_df["correlation_score"] > min_correlation)
        & (pairwise_df["left_condition_label"] != pairwise_df["right_condition_label"]),
        [
            "left_condition_label",
            "right_condition_label",
            "correlation_score",
            "overlap_pvalue",
        ],
    ].copy()

    filtered_df.columns = [
        "Condition 1",
        "Condition 2",
        "Correlation Value",
        "P-value",
    ]
    return filtered_df.reset_index(drop=True)


def build_full_correlation_matrices(entries):
    """Build symmetric condition-by-condition matrices for Jaccard and p-values."""
    labels = []
    label_counts = {}
    for entry in entries:
        base = entry["condition_label"]
        count = label_counts.get(base, 0) + 1
        label_counts[base] = count
        label = base if count == 1 else f"{base} [{count}]"
        entry["matrix_label"] = label
        labels.append(label)

    universe_genes = set()
    for entry in entries:
        universe_genes.update(entry["leading_gene_set"])
    universe_size = len(universe_genes)

    correlation_matrix = pd.DataFrame(np.eye(len(entries)), index=labels, columns=labels, dtype=float)
    pvalue_matrix = pd.DataFrame(np.nan, index=labels, columns=labels, dtype=float)

    for i, left_entry in enumerate(entries):
        left_size = len(left_entry["leading_gene_set"])
        pvalue_matrix.iloc[i, i] = compute_overlap_pvalue(
            universe_size=universe_size,
            left_size=left_size,
            right_size=left_size,
            overlap_size=left_size,
        )

        for j in range(i + 1, len(entries)):
            right_entry = entries[j]
            overlap_size = len(left_entry["leading_gene_set"] & right_entry["leading_gene_set"])
            score = compute_jaccard_score(
                left_entry["leading_gene_set"], right_entry["leading_gene_set"]
            )
            pvalue = compute_overlap_pvalue(
                universe_size=universe_size,
                left_size=left_size,
                right_size=len(right_entry["leading_gene_set"]),
                overlap_size=overlap_size,
            )
            correlation_matrix.iloc[i, j] = score
            correlation_matrix.iloc[j, i] = score
            pvalue_matrix.iloc[i, j] = pvalue
            pvalue_matrix.iloc[j, i] = pvalue

    return correlation_matrix, pvalue_matrix


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


def format_pvalue_text(pvalue):
    """Format a p-value compactly for heatmap annotations."""
    if pd.isna(pvalue):
        return ""
    if pvalue < 1e-3:
        return f"p={pvalue:.1e}"
    return f"p={pvalue:.3f}"


def plot_top_correlation_heatmap(
    matrix,
    pvalue_matrix,
    output_path,
    pvalue_annotation_threshold=DEFAULT_HEATMAP_PVALUE_THRESHOLD,
):
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

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            score = matrix.iloc[i, j]
            if score <= pvalue_annotation_threshold:
                continue
            pvalue = pvalue_matrix.iloc[i, j]
            pvalue_text = format_pvalue_text(pvalue)
            if not pvalue_text:
                continue
            ax.text(
                j + 0.5,
                i + 0.5,
                pvalue_text,
                ha="center",
                va="center",
                fontsize=6.5,
                color="white" if score >= 0.55 else "black",
            )

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
    filtered_export_df = build_filtered_correlation_export(pairwise_df)
    full_matrix, full_pvalue_matrix = build_full_correlation_matrices(entries)
    top_matrix = select_top_conditions(full_matrix, top_n=top_n)
    top_pvalue_matrix = full_pvalue_matrix.loc[top_matrix.index, top_matrix.columns]

    csv_path = output_dir / "leading_gene_condition_correlations.csv"
    heatmap_path = output_dir / "leading_gene_condition_correlation_heatmap.png"

    filtered_export_df.to_csv(csv_path, index=False)
    plot_top_correlation_heatmap(top_matrix, top_pvalue_matrix, heatmap_path)

    return {
        "output_dir": str(output_dir),
        "filtered_csv_path": str(csv_path),
        "heatmap_path": str(heatmap_path),
        "n_summary_files": len(summary_files),
        "n_conditions": len(entries),
        "n_exported_condition_pairs": len(filtered_export_df),
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
