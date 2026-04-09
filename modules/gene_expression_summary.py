import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_OUTPUT_DIR = Path("results/additional_analyses/gene_expression_summary")
META_REGIONS = ["Forebrain", "Midbrain", "Hindbrain"]


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
    return cleaned or "gene_expression_summary"


def read_table(path):
    """Read CSV/TSV/TXT tables using the file suffix to choose a separator."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    if path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def load_summary_row(summary_file, condition):
    """Load the GSEA summary row for the requested condition."""
    summary_df = read_table(summary_file)

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
    leading_genes = parse_gene_list(row["Lead_genes"])
    if not leading_genes:
        raise ValueError(f"No leading genes found for condition '{condition}'.")

    return {
        "column": str(row["column"]),
        "condition": str(row["condition"]),
        "leading_genes": leading_genes,
    }


def load_condition_leading_gene_sets(summary_file, conditions=None):
    """Load leading-gene sets for requested summary conditions."""
    summary_df = read_table(summary_file)

    required_cols = {"condition", "Lead_genes"}
    missing_cols = required_cols - set(summary_df.columns)
    if missing_cols:
        raise ValueError(
            f"Leading-gene summary file is missing required columns: {sorted(missing_cols)}"
        )

    if conditions is None:
        ordered_conditions = summary_df["condition"].astype(str).drop_duplicates().tolist()
    else:
        ordered_conditions = [str(condition) for condition in conditions]

    condition_gene_sets = {}
    for condition in ordered_conditions:
        matched = summary_df.loc[summary_df["condition"].astype(str) == condition].copy()
        if matched.empty:
            raise ValueError(f"Condition '{condition}' was not found in {summary_file}.")

        if len(matched) > 1 and "NES" in matched.columns:
            matched["NES_abs"] = pd.to_numeric(matched["NES"], errors="coerce").abs()
            matched = matched.sort_values("NES_abs", ascending=False)

        genes = parse_gene_list(matched.iloc[0]["Lead_genes"])
        condition_gene_sets[condition] = set(genes)

    return condition_gene_sets


def build_leading_gene_flag_columns(genes, condition_gene_sets):
    """Create per-gene leading-gene membership flags for additional conditions."""
    flag_df = pd.DataFrame({"gene": list(genes)})
    for condition, gene_set in condition_gene_sets.items():
        column_name = f"{sanitize_name(condition)}_leading_gene"
        flag_df[column_name] = flag_df["gene"].isin(gene_set)
    return flag_df


def collapse_region(region):
    """Map detailed regions to meta-regions."""
    if pd.isna(region):
        return np.nan

    forebrain = {"Forebrain", "Telencephalon", "Diencephalon"}
    hindbrain = {"Hindbrain", "Cerebellum", "Pons", "Medulla"}
    midbrain = {"Midbrain"}

    if region in forebrain:
        return "Forebrain"
    if region in hindbrain:
        return "Hindbrain"
    if region in midbrain:
        return "Midbrain"
    return np.nan


def compute_top_vs_second_wilcoxon(final_summary):
    """Run a one-sided paired Wilcoxon signed-rank test across genes."""
    from scipy.stats import wilcoxon

    valid = final_summary.loc[
        final_summary["top_region_score"].notna() & final_summary["second_region_score"].notna()
    ].copy()

    result = {
        "wilcoxon_top_vs_second_statistic": np.nan,
        "wilcoxon_top_vs_second_pvalue": np.nan,
        "wilcoxon_top_vs_second_n_pairs": int(len(valid)),
        "wilcoxon_top_vs_second_alternative": "greater",
        "wilcoxon_top_vs_second_note": "",
    }

    if valid.empty:
        result["wilcoxon_top_vs_second_note"] = "No genes had both top and second region scores."
        return result

    differences = valid["top_region_score"] - valid["second_region_score"]
    nonzero_mask = differences != 0
    if not nonzero_mask.any():
        result["wilcoxon_top_vs_second_note"] = (
            "All top-vs-second region score differences were zero."
        )
        return result

    valid_nonzero = valid.loc[nonzero_mask]
    result["wilcoxon_top_vs_second_n_pairs"] = int(len(valid_nonzero))

    test = wilcoxon(
        valid_nonzero["top_region_score"],
        valid_nonzero["second_region_score"],
        alternative="greater",
        zero_method="wilcox",
    )
    result["wilcoxon_top_vs_second_statistic"] = float(test.statistic)
    result["wilcoxon_top_vs_second_pvalue"] = float(test.pvalue)
    return result


def list_leading_gene_flag_columns(df):
    """Return boolean flag columns created from optional leading-gene condition sets."""
    return [column for column in df.columns if column.endswith("_leading_gene")]


def list_region_score_columns(df):
    """Return region-score columns and their corresponding region names."""
    return [
        (region_name, f"{region_name}_region_score")
        for region_name in META_REGIONS
        if f"{region_name}_region_score" in df.columns
    ]


def compute_gene_set_top_region(gene_set_summary):
    """Pick the top region for a gene set using mean region scores across genes."""
    region_columns = list_region_score_columns(gene_set_summary)
    if not region_columns or gene_set_summary.empty:
        return np.nan

    mean_scores = {
        region_name: pd.to_numeric(gene_set_summary[column], errors="coerce").mean()
        for region_name, column in region_columns
    }
    mean_scores = pd.Series(mean_scores).dropna()
    if mean_scores.empty:
        return np.nan
    return str(mean_scores.idxmax())


def run_selected_region_wilcoxon(gene_set_summary, selected_region, other_region):
    """Run a one-sided paired Wilcoxon test for one selected-vs-other region comparison."""
    from scipy.stats import wilcoxon

    selected_col = f"{selected_region}_region_score"
    other_col = f"{other_region}_region_score"
    valid = gene_set_summary.loc[
        gene_set_summary[selected_col].notna() & gene_set_summary[other_col].notna()
    ].copy()

    result = {
        "selected_region": selected_region,
        "comparison_region": other_region,
        "wilcoxon_statistic": np.nan,
        "wilcoxon_pvalue": np.nan,
        "wilcoxon_n_pairs": int(len(valid)),
        "wilcoxon_alternative": "greater",
        "wilcoxon_note": "",
    }

    if valid.empty:
        result["wilcoxon_note"] = (
            f"No genes had both {selected_region} and {other_region} region scores."
        )
        return result

    differences = valid[selected_col] - valid[other_col]
    nonzero_mask = differences != 0
    if not nonzero_mask.any():
        result["wilcoxon_note"] = (
            f"All {selected_region}-vs-{other_region} region score differences were zero."
        )
        return result

    valid_nonzero = valid.loc[nonzero_mask]
    result["wilcoxon_n_pairs"] = int(len(valid_nonzero))

    test = wilcoxon(
        valid_nonzero[selected_col],
        valid_nonzero[other_col],
        alternative="greater",
        zero_method="wilcox",
    )
    result["wilcoxon_statistic"] = float(test.statistic)
    result["wilcoxon_pvalue"] = float(test.pvalue)
    return result


def build_leading_gene_group_wilcoxon_summary(final_summary, selected_region):
    """Create one summary table for all genes and each optional leading-gene subgroup."""
    available_regions = [region_name for region_name, _ in list_region_score_columns(final_summary)]
    if selected_region not in available_regions:
        raise ValueError(
            f"Selected Wilcoxon region '{selected_region}' is not available. "
            f"Choose from: {available_regions}"
        )

    group_frames = [("all", final_summary)]
    for flag_column in list_leading_gene_flag_columns(final_summary):
        group_name = flag_column[: -len("_leading_gene")]
        group_df = final_summary.loc[final_summary[flag_column] == True].copy()
        group_frames.append((group_name, group_df))

    rows = []
    other_regions = [region for region in available_regions if region != selected_region]

    for group_name, group_df in group_frames:
        top_region = compute_gene_set_top_region(group_df)
        for other_region in other_regions:
            row = {
                "leading_gene_group": group_name,
                "top_region_for_gene_set": top_region,
            }
            row.update(run_selected_region_wilcoxon(group_df, selected_region, other_region))
            rows.append(row)

    return pd.DataFrame(rows)


def summarize_gene_expression_top_region(
    adata,
    genes,
    region_col="Region",
    cell_filter_col="CellClass",
    cell_filter_val="Radial glia",
    sym_col="Gene",
    expr_threshold=0,
    leading_gene_flags=None,
):
    """
    Returns:
      1) long_summary: one row per gene x meta-region
      2) final_summary: one row per gene with top region and region-specific stats
      3) wilcoxon_result: paired one-sided Wilcoxon top-vs-second score test across genes
    """

    if region_col not in adata.obs.columns:
        raise KeyError(f"Column '{region_col}' not found in adata.obs")
    if cell_filter_col not in adata.obs.columns:
        raise KeyError(f"Column '{cell_filter_col}' not found in adata.obs")
    if sym_col not in adata.var.columns:
        raise KeyError(f"Column '{sym_col}' not found in adata.var")

    print("Filtering cells...")
    adata_f = adata[adata.obs[cell_filter_col].astype(str) == str(cell_filter_val)].copy()
    if adata_f.n_obs == 0:
        raise ValueError(
            f"No cells found for condition '{cell_filter_val}' in adata.obs['{cell_filter_col}']."
        )

    print("Normalizing...")
    import scanpy as sc
    import scipy.sparse as sp

    sc.pp.normalize_total(adata_f)
    sc.pp.log1p(adata_f)

    print("Collapsing regions...")
    adata_f.obs["meta_region"] = adata_f.obs[region_col].map(collapse_region)
    adata_f = adata_f[adata_f.obs["meta_region"].notna()].copy()
    if adata_f.n_obs == 0:
        raise ValueError("No cells remained after collapsing to Forebrain/Midbrain/Hindbrain.")

    sym2var = (
        pd.Series(adata_f.var_names.values, index=adata_f.var[sym_col].astype(str))
        .dropna()
        .to_dict()
    )

    found_genes = [g for g in genes if g in sym2var]
    missing_genes = [g for g in genes if g not in sym2var]

    print(f"Found {len(found_genes)} genes")
    print(f"Missing {len(missing_genes)} genes")
    if missing_genes:
        print("Missing symbols:", missing_genes)
    if not found_genes:
        raise ValueError("None of the requested genes were found in adata.var[sym_col].")

    varnames = [sym2var[g] for g in found_genes]
    X = adata_f[:, varnames].X
    if sp.issparse(X):
        X = X.toarray()

    meta_regions = adata_f.obs["meta_region"].values
    long_results = []

    print("Calculating summaries...")
    for i, gene in enumerate(found_genes):
        vals = X[:, i]

        tmp = pd.DataFrame({"meta_region": meta_regions, "expr": vals})
        grouped_all = tmp.groupby("meta_region")["expr"]

        n_total = grouped_all.size().reset_index(name="n_total_cells")
        frac_expr = grouped_all.apply(
            lambda x: (x > expr_threshold).mean()
        ).reset_index(name="fraction_expressing")

        tmp_expr = tmp[tmp["expr"] > expr_threshold].copy()
        grouped_expr = tmp_expr.groupby("meta_region")["expr"]

        summary_expr = grouped_expr.agg(
            mean_expr="mean",
            median_expr="median",
            std_expr="std",
            min_expr="min",
            max_expr="max",
            n_expressing_cells="count",
        ).reset_index()

        summary = n_total.merge(frac_expr, on="meta_region", how="left")
        summary = summary.merge(summary_expr, on="meta_region", how="left")
        summary["region_score"] = summary["mean_expr"] * summary["fraction_expressing"]
        summary["gene"] = gene
        long_results.append(summary)

    long_summary = pd.concat(long_results, ignore_index=True)

    final_rows = []

    for gene, sub in long_summary.groupby("gene"):
        sub = sub.set_index("meta_region").reindex(META_REGIONS)
        scores = sub["region_score"].fillna(0)

        top_region = scores.idxmax()
        top_score = scores.max()

        sorted_scores = scores.sort_values(ascending=False)
        second_region = sorted_scores.index[1] if len(sorted_scores) > 1 else np.nan
        second_score = sorted_scores.iloc[1] if len(sorted_scores) > 1 else np.nan

        enrichment_score = (
            top_score / second_score
            if pd.notna(second_score) and second_score > 0
            else np.nan
        )

        row = {
            "gene": gene,
            "top_region": top_region,
            "top_region_score": top_score,
            "second_region": second_region,
            "second_region_score": second_score,
            "enrichment_score_top_vs_second": enrichment_score,
        }

        for region in META_REGIONS:
            row[f"{region}_mean_expr"] = sub.loc[region, "mean_expr"]
            row[f"{region}_median_expr"] = sub.loc[region, "median_expr"]
            row[f"{region}_std_expr"] = sub.loc[region, "std_expr"]
            row[f"{region}_fraction_expressing"] = sub.loc[region, "fraction_expressing"]
            row[f"{region}_n_total_cells"] = sub.loc[region, "n_total_cells"]
            row[f"{region}_n_expressing_cells"] = sub.loc[region, "n_expressing_cells"]
            row[f"{region}_region_score"] = sub.loc[region, "region_score"]

        final_rows.append(row)

    final_summary = pd.DataFrame(final_rows).sort_values(
        ["top_region", "enrichment_score_top_vs_second"],
        ascending=[True, False],
    )

    if leading_gene_flags is not None:
        long_summary = long_summary.merge(leading_gene_flags, on="gene", how="left")
        final_summary = final_summary.merge(leading_gene_flags, on="gene", how="left")

    wilcoxon_result = compute_top_vs_second_wilcoxon(final_summary)
    for column, value in wilcoxon_result.items():
        final_summary[column] = value

    return long_summary, final_summary, wilcoxon_result


def create_gene_expression_summary(
    h5ad_path,
    gsea_summary_file,
    condition,
    output_dir=DEFAULT_OUTPUT_DIR,
    subfolder_name=None,
    leading_gene_summary_file=None,
    leading_gene_conditions=None,
    region_col="Region",
    cell_filter_col="CellClass",
    cell_filter_val="Radial glia",
    sym_col="Gene",
    expr_threshold=0,
    wilcoxon_region="Forebrain",
):
    """Create long and final gene-expression summaries for a GSEA leading-gene set."""
    summary_info = load_summary_row(gsea_summary_file, condition)

    output_dir = Path(output_dir)
    if subfolder_name:
        output_dir = output_dir / sanitize_name(subfolder_name)
    else:
        output_dir = output_dir / sanitize_name(condition)
    output_dir.mkdir(parents=True, exist_ok=True)

    leading_gene_flags = None
    if leading_gene_summary_file is not None:
        condition_gene_sets = load_condition_leading_gene_sets(
            leading_gene_summary_file,
            conditions=leading_gene_conditions,
        )
        leading_gene_flags = build_leading_gene_flag_columns(
            summary_info["leading_genes"],
            condition_gene_sets,
        )

    print("Loading AnnData...")
    import scanpy as sc

    adata = sc.read_h5ad(h5ad_path)

    long_summary, final_summary, wilcoxon_result = summarize_gene_expression_top_region(
        adata,
        summary_info["leading_genes"],
        region_col=region_col,
        cell_filter_col=cell_filter_col,
        cell_filter_val=cell_filter_val,
        sym_col=sym_col,
        expr_threshold=expr_threshold,
        leading_gene_flags=leading_gene_flags,
    )

    long_csv_path = output_dir / "gene_expression_region_long_summary.csv"
    final_csv_path = output_dir / "gene_expression_top_region_summary.csv"
    group_wilcoxon_csv_path = output_dir / "gene_expression_group_wilcoxon_summary.csv"
    long_summary.to_csv(long_csv_path, index=False)
    final_summary.to_csv(final_csv_path, index=False)
    group_wilcoxon_summary = build_leading_gene_group_wilcoxon_summary(
        final_summary,
        selected_region=wilcoxon_region,
    )
    group_wilcoxon_summary.to_csv(group_wilcoxon_csv_path, index=False)

    return {
        "output_dir": str(output_dir),
        "long_summary_csv_path": str(long_csv_path),
        "final_summary_csv_path": str(final_csv_path),
        "group_wilcoxon_summary_csv_path": str(group_wilcoxon_csv_path),
        "gsea_summary_file": str(Path(gsea_summary_file)),
        "condition": summary_info["condition"],
        "condition_column": summary_info["column"],
        "n_input_leading_genes": len(summary_info["leading_genes"]),
        "n_found_genes": int(final_summary.shape[0]),
        "wilcoxon_region": wilcoxon_region,
        "wilcoxon_top_vs_second_statistic": wilcoxon_result[
            "wilcoxon_top_vs_second_statistic"
        ],
        "wilcoxon_top_vs_second_pvalue": wilcoxon_result[
            "wilcoxon_top_vs_second_pvalue"
        ],
        "wilcoxon_top_vs_second_n_pairs": wilcoxon_result[
            "wilcoxon_top_vs_second_n_pairs"
        ],
    }


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize top regional expression for leading genes taken from a GSEA summary condition."
        )
    )
    parser.add_argument(
        "--h5ad-path",
        required=True,
        help="Path to the input AnnData .h5ad file.",
    )
    parser.add_argument(
        "--gsea-summary-file",
        required=True,
        help="GSEA summary CSV/TSV/TXT containing 'column', 'condition', and 'Lead_genes'.",
    )
    parser.add_argument(
        "--condition",
        required=True,
        help="Condition value to match in the GSEA summary and use for leading genes.",
    )
    parser.add_argument(
        "--leading-gene-summary-file",
        default=None,
        help=(
            "Optional separate summary CSV/TSV/TXT used to add per-condition leading-gene "
            "membership flags."
        ),
    )
    parser.add_argument(
        "--leading-gene-conditions",
        nargs="+",
        default=None,
        help=(
            "Optional list of condition names to annotate from --leading-gene-summary-file. "
            "If omitted, all conditions in that file are used."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Base output folder for the generated CSV files.",
    )
    parser.add_argument(
        "--subfolder-name",
        default=None,
        help="Optional subfolder name under the output directory.",
    )
    parser.add_argument(
        "--region-col",
        default="Region",
        help="Column in adata.obs describing anatomical regions.",
    )
    parser.add_argument(
        "--cell-filter-col",
        default="CellClass",
        help="Column in adata.obs used to filter cells before summarizing.",
    )
    parser.add_argument(
        "--cell-filter-val",
        default="Radial glia",
        help="Value in --cell-filter-col to keep before summarizing.",
    )
    parser.add_argument(
        "--sym-col",
        default="Gene",
        help="Column in adata.var containing gene symbols.",
    )
    parser.add_argument(
        "--expr-threshold",
        type=float,
        default=0,
        help="Expression threshold used for fraction-expressing and expressing-cell summaries.",
    )
    parser.add_argument(
        "--wilcoxon-region",
        default="Forebrain",
        help=(
            "Region to compare against each other available region in the group Wilcoxon summary CSV."
        ),
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    result = create_gene_expression_summary(
        h5ad_path=args.h5ad_path,
        gsea_summary_file=args.gsea_summary_file,
        condition=args.condition,
        output_dir=args.output_dir,
        subfolder_name=args.subfolder_name,
        leading_gene_summary_file=args.leading_gene_summary_file,
        leading_gene_conditions=args.leading_gene_conditions,
        region_col=args.region_col,
        cell_filter_col=args.cell_filter_col,
        cell_filter_val=args.cell_filter_val,
        sym_col=args.sym_col,
        expr_threshold=args.expr_threshold,
        wilcoxon_region=args.wilcoxon_region,
    )

    print("Saved gene expression summary results:")
    for key, value in result.items():
        print(f"  {key}: {value}")
