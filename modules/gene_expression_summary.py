import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_OUTPUT_DIR = Path("results/additional_analyses/gene_expression_summary")
META_REGIONS = ["Forebrain", "Midbrain", "Hindbrain"]
REGION_COLLAPSE_MAP = {
    "Forebrain": {"Forebrain", "Telencephalon", "Diencephalon"},
    "Midbrain": {"Midbrain"},
    "Hindbrain": {"Hindbrain", "Cerebellum", "Pons", "Medulla"},
}
DEFAULT_CHEMISTRY = "v3"
DEFAULT_CHEMISTRY_COL = "Chemistry"


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


def add_fdr_column(df, pvalue_col, fdr_col):
    """Add BH/FDR-corrected p-values to a DataFrame."""
    import statsmodels.stats.multitest as smm

    df = df.copy()
    df[fdr_col] = np.nan

    numeric_pvalues = pd.to_numeric(df[pvalue_col], errors="coerce")
    valid_mask = numeric_pvalues.notna()
    if valid_mask.any():
        df.loc[valid_mask, fdr_col] = smm.multipletests(
            numeric_pvalues.loc[valid_mask], method="fdr_bh"
        )[1]
    return df


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

    for meta_region in META_REGIONS:
        if region in REGION_COLLAPSE_MAP[meta_region]:
            return meta_region
    return np.nan


def prepare_expression_adata(
    adata,
    region_col="Region",
    cell_filter_col="CellClass",
    cell_filter_val="Radial glia",
    normalize=True,
):
    """Filter cells, optionally normalize/log-transform, and collapse to shared meta-regions."""
    if region_col not in adata.obs.columns:
        raise KeyError(f"Column '{region_col}' not found in adata.obs")
    if cell_filter_col not in adata.obs.columns:
        raise KeyError(f"Column '{cell_filter_col}' not found in adata.obs")

    if cell_filter_val is not None:
        print(f"Filtering cells: {cell_filter_col} == '{cell_filter_val}'")
        adata_f = adata[adata.obs[cell_filter_col].astype(str) == str(cell_filter_val)].copy()
        if adata_f.n_obs == 0:
            raise ValueError(
                f"No cells found for condition '{cell_filter_val}' in adata.obs['{cell_filter_col}']."
            )
    else:
        print("No cell-type filter applied — using all cells.")
        adata_f = adata.copy()

    if normalize:
        print("Normalizing...")
        import scanpy as sc

        sc.pp.normalize_total(adata_f)
        sc.pp.log1p(adata_f)

    print("Collapsing regions...")
    adata_f.obs["meta_region"] = adata_f.obs[region_col].map(collapse_region)
    adata_f = adata_f[adata_f.obs["meta_region"].notna()].copy()
    if adata_f.n_obs == 0:
        raise ValueError("No cells remained after collapsing to Forebrain/Midbrain/Hindbrain.")

    return adata_f


def build_expression_dataframe(adata, genes, sym_col="Gene"):
    """Return a dense expression table for the requested genes plus meta-region labels."""
    import scipy.sparse as sp

    if sym_col not in adata.var.columns:
        raise KeyError(f"Column '{sym_col}' not found in adata.var")
    if "meta_region" not in adata.obs.columns:
        raise KeyError("Column 'meta_region' not found in adata.obs")

    sym2var = (
        pd.Series(adata.var_names.values, index=adata.var[sym_col].astype(str))
        .dropna()
        .to_dict()
    )

    found_genes = [gene for gene in genes if gene in sym2var]
    if not found_genes:
        raise ValueError("None of the requested genes were found in adata.var[sym_col].")

    varnames = [sym2var[gene] for gene in found_genes]
    X = adata[:, varnames].X
    if sp.issparse(X):
        X = X.toarray()

    expr_df = pd.DataFrame(X, columns=found_genes)
    expr_df["meta_region"] = adata.obs["meta_region"].to_numpy()
    invalid_regions = set(expr_df["meta_region"].dropna().astype(str)) - set(META_REGIONS)
    if invalid_regions:
        raise ValueError(
            "Found unexpected meta-region values after region collapsing: "
            f"{sorted(invalid_regions)}"
        )
    return expr_df, found_genes


def get_unordered_meta_region_pairs():
    """Return all unordered pairs of collapsed meta-regions."""
    return [
        (META_REGIONS[i], META_REGIONS[j])
        for i in range(len(META_REGIONS))
        for j in range(i + 1, len(META_REGIONS))
    ]


def compute_top_vs_second_wilcoxon(final_summary):
    """Run a one-sided paired Wilcoxon signed-rank test across genes."""
    from scipy.stats import wilcoxon
    import statsmodels.stats.multitest as smm

    valid = final_summary.loc[
        final_summary["top_region_score"].notna() & final_summary["second_region_score"].notna()
    ].copy()

    result = {
        "wilcoxon_top_vs_second_statistic": np.nan,
        "wilcoxon_top_vs_second_pvalue": np.nan,
        "wilcoxon_top_vs_second_fdr_bh_pvalue": np.nan,
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
    result["wilcoxon_top_vs_second_fdr_bh_pvalue"] = float(
        smm.multipletests([float(test.pvalue)], method="fdr_bh")[1][0]
    )
    return result


def get_gene_flag_column(condition_name):
    """Return the flag-column name used for a leading-gene condition."""
    return f"{sanitize_name(condition_name)}_leading_gene"


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


def build_gene_group_frames(final_summary):
    """Build the requested group summaries for the top-region group CSV."""
    group_frames = [("all", final_summary)]

    g1_col = get_gene_flag_column("G1")
    s_col = get_gene_flag_column("S")
    g2m_col = get_gene_flag_column("G2M")

    if g1_col in final_summary.columns and s_col in final_summary.columns:
        group_frames.append(
            ("G1 & S", final_summary.loc[final_summary[g1_col] & final_summary[s_col]].copy())
        )
    if s_col in final_summary.columns and g2m_col in final_summary.columns:
        group_frames.append(
            ("S & G2M", final_summary.loc[final_summary[s_col] & final_summary[g2m_col]].copy())
        )

    return group_frames


def get_combined_gene_group_label(row):
    """Label a gene as belonging to G1 & S, S & G2M, or none."""
    g1_flag = bool(row.get(get_gene_flag_column("G1"), False))
    s_flag = bool(row.get(get_gene_flag_column("S"), False))
    g2m_flag = bool(row.get(get_gene_flag_column("G2M"), False))

    if g1_flag and s_flag:
        return "G1 & S"
    if s_flag and g2m_flag:
        return "S & G2M"
    return "none"


def build_leading_gene_group_wilcoxon_summary(final_summary, selected_region):
    """Create one summary table for all genes and each optional leading-gene subgroup."""
    available_regions = [region_name for region_name, _ in list_region_score_columns(final_summary)]
    if selected_region not in available_regions:
        raise ValueError(
            f"Selected Wilcoxon region '{selected_region}' is not available. "
            f"Choose from: {available_regions}"
        )

    group_frames = build_gene_group_frames(final_summary)

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

    return add_fdr_column(pd.DataFrame(rows), "wilcoxon_pvalue", "wilcoxon_fdr_bh_pvalue")


def build_per_gene_region_wilcoxon_summary(
    adata,
    genes,
    selected_region,
    sym_col="Gene",
    leading_gene_flags=None,
    compare_all_region_pairs=False,
):
    """Create per-gene region-vs-region Wilcoxon rank-sum summaries using expressing cells only."""
    from scipy.stats import mannwhitneyu
    if not compare_all_region_pairs and selected_region not in META_REGIONS:
        raise ValueError(f"Selected region '{selected_region}' must be one of {META_REGIONS}")

    expr_df, found_genes = build_expression_dataframe(adata, genes, sym_col=sym_col)

    rows = []
    if compare_all_region_pairs:
        region_pairs = [
            (region_of_interest, region_of_comparison)
            for region_of_interest in META_REGIONS
            for region_of_comparison in META_REGIONS
            if region_of_interest != region_of_comparison
        ]
    else:
        region_pairs = [
            (selected_region, region_of_comparison)
            for region_of_comparison in META_REGIONS
            if region_of_comparison != selected_region
        ]
    gene_group_labels = {}
    if leading_gene_flags is not None and not leading_gene_flags.empty:
        gene_group_labels = (
            leading_gene_flags.assign(
                combined_leading_gene_group=leading_gene_flags.apply(
                    get_combined_gene_group_label,
                    axis=1,
                )
            )
            .set_index("gene")["combined_leading_gene_group"]
            .to_dict()
        )

    for gene in found_genes:
        for region_of_interest, comparison_region in region_pairs:
            selected_values = expr_df.loc[
                (expr_df["meta_region"] == region_of_interest) & (expr_df[gene] > 0),
                gene,
            ].astype(float)
            comparison_values = expr_df.loc[
                (expr_df["meta_region"] == comparison_region) & (expr_df[gene] > 0),
                gene,
            ].astype(float)

            row = {
                "gene": gene,
                "combined_leading_gene_group": gene_group_labels.get(gene, "none"),
                "region_of_interest": region_of_interest,
                "region_of_comparison": comparison_region,
                "selected_region": region_of_interest,
                "comparison_region": comparison_region,
                "n_selected_region_expressing_cells": int(len(selected_values)),
                "n_comparison_region_expressing_cells": int(len(comparison_values)),
                "selected_region_mean_expr": float(selected_values.mean()) if len(selected_values) else np.nan,
                "comparison_region_mean_expr": float(comparison_values.mean()) if len(comparison_values) else np.nan,
                "wilcoxon_statistic": np.nan,
                "wilcoxon_pvalue": np.nan,
                "wilcoxon_alternative": "greater",
                "wilcoxon_note": "",
            }

            if len(selected_values) == 0 or len(comparison_values) == 0:
                row["wilcoxon_note"] = (
                    "No expressing cells were available for one or both regions."
                )
            else:
                test = mannwhitneyu(
                    selected_values,
                    comparison_values,
                    alternative="greater",
                )
                row["wilcoxon_statistic"] = float(test.statistic)
                row["wilcoxon_pvalue"] = float(test.pvalue)

            rows.append(row)

    return add_fdr_column(pd.DataFrame(rows), "wilcoxon_pvalue", "wilcoxon_fdr_bh_pvalue")


def get_combined_leading_gene_group_order(per_gene_wilcoxon_summary):
    """Return a stable display order for combined leading-gene groups."""
    preferred_order = ["G1 & S", "S & G2M", "none"]
    observed_groups = (
        per_gene_wilcoxon_summary["combined_leading_gene_group"]
        .fillna("none")
        .astype(str)
        .drop_duplicates()
        .tolist()
    )
    ordered_groups = [group for group in preferred_order if group in observed_groups]
    ordered_groups.extend(sorted(group for group in observed_groups if group not in preferred_order))
    return ordered_groups


def create_region_pair_significant_gene_count_plot(per_gene_wilcoxon_summary, output_path):
    """Plot significant per-gene Wilcoxon counts for each unordered region pair."""
    required_cols = {
        "combined_leading_gene_group",
        "region_of_interest",
        "region_of_comparison",
        "wilcoxon_fdr_bh_pvalue",
    }
    missing_cols = required_cols - set(per_gene_wilcoxon_summary.columns)
    if missing_cols:
        raise ValueError(
            "Per-gene Wilcoxon summary is missing required columns for plotting: "
            f"{sorted(missing_cols)}"
        )

    group_order = get_combined_leading_gene_group_order(per_gene_wilcoxon_summary)
    region_pairs = get_unordered_meta_region_pairs()
    colors = {
        "left": "#ff5a52",
        "right": "#16b3b1",
    }

    fig, axes = plt.subplots(
        1,
        len(region_pairs),
        figsize=(5 * len(region_pairs), max(3.5, 1.1 * len(group_order) + 1.5)),
        sharey=True,
    )
    if len(region_pairs) == 1:
        axes = [axes]

    max_abs_count = 0
    pair_counts = {}
    significant_mask = pd.to_numeric(
        per_gene_wilcoxon_summary["wilcoxon_fdr_bh_pvalue"], errors="coerce"
    ) < 0.05
    significant_summary = per_gene_wilcoxon_summary.loc[significant_mask].copy()

    for left_region, right_region in region_pairs:
        left_counts = (
            significant_summary.loc[
                (significant_summary["region_of_interest"] == left_region)
                & (significant_summary["region_of_comparison"] == right_region)
            ]
            .groupby("combined_leading_gene_group")["gene"]
            .nunique()
            .reindex(group_order, fill_value=0)
        )
        right_counts = (
            significant_summary.loc[
                (significant_summary["region_of_interest"] == right_region)
                & (significant_summary["region_of_comparison"] == left_region)
            ]
            .groupby("combined_leading_gene_group")["gene"]
            .nunique()
            .reindex(group_order, fill_value=0)
        )
        pair_counts[(left_region, right_region)] = (left_counts, right_counts)
        max_abs_count = max(
            max_abs_count,
            int(left_counts.max()) if len(left_counts) else 0,
            int(right_counts.max()) if len(right_counts) else 0,
        )

    x_limit = max(1, max_abs_count)
    y_positions = np.arange(len(group_order))

    for ax, (left_region, right_region) in zip(axes, region_pairs):
        left_counts, right_counts = pair_counts[(left_region, right_region)]
        ax.barh(
            y_positions,
            -left_counts.to_numpy(),
            color=colors["left"],
            edgecolor="white",
            height=0.7,
        )
        ax.barh(
            y_positions,
            right_counts.to_numpy(),
            color=colors["right"],
            edgecolor="white",
            height=0.7,
        )
        ax.axvline(0, color="#666666", linewidth=1)
        ax.set_xlim(-x_limit - 0.5, x_limit + 0.5)
        ax.set_title(f"{left_region} vs {right_region}", fontsize=11, pad=18)
        ax.set_yticks(y_positions)
        ax.set_yticklabels(group_order, fontsize=10)
        ax.grid(axis="x", alpha=0.3)
        ax.set_axisbelow(True)
        ax.invert_yaxis()
        ax.text(
            0.02,
            1.04,
            left_region,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=10,
            color=colors["left"],
        )
        ax.text(
            0.98,
            1.04,
            right_region,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=10,
            color=colors["right"],
        )

    axes[0].set_ylabel("combined_leading_gene_group", fontsize=10)
    for ax in axes:
        ax.set_xlabel("significant gene count", fontsize=10)

    fig.suptitle(
        "Genes with region-enriched expression by Wilcoxon FDR < 0.05",
        fontsize=13,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def get_requested_region_pairs(selected_region, compare_all_region_pairs):
    """Return ordered region pairs requested for per-gene comparisons."""
    if compare_all_region_pairs:
        return [
            (region_of_interest, region_of_comparison)
            for region_of_interest in META_REGIONS
            for region_of_comparison in META_REGIONS
            if region_of_interest != region_of_comparison
        ]
    return [
        (selected_region, region_of_comparison)
        for region_of_comparison in META_REGIONS
        if region_of_comparison != selected_region
    ]


def create_significant_gene_boxplots(
    adata,
    genes,
    per_gene_wilcoxon_summary,
    output_dir,
    selected_region,
    compare_all_region_pairs=False,
    sym_col="Gene",
):
    """Save per-gene boxplots across all meta-regions for strong requested enrichments."""
    expr_df, found_genes = build_expression_dataframe(adata, genes, sym_col=sym_col)
    requested_region_pairs = get_requested_region_pairs(
        selected_region=selected_region,
        compare_all_region_pairs=compare_all_region_pairs,
    )
    significant_summary = per_gene_wilcoxon_summary.loc[
        pd.to_numeric(per_gene_wilcoxon_summary["wilcoxon_fdr_bh_pvalue"], errors="coerce") < 0.05
    ].copy()
    significant_summary = significant_summary.loc[
        significant_summary.apply(
            lambda row: (row["region_of_interest"], row["region_of_comparison"])
            in requested_region_pairs,
            axis=1,
        )
    ]
    requested_regions_of_interest = sorted(
        {region_of_interest for region_of_interest, _ in requested_region_pairs},
        key=META_REGIONS.index,
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    colors = {
        "Forebrain": "#ff5a52",
        "Midbrain": "#16b3b1",
        "Hindbrain": "#4c78a8",
    }
    n_saved_plots = 0

    total_cells_by_region = expr_df["meta_region"].value_counts()

    for region_of_interest in requested_regions_of_interest:
        region_rows = significant_summary.loc[
            significant_summary["region_of_interest"] == region_of_interest
        ].copy()
        if region_rows.empty:
            continue

        successful_other_regions = (
            region_rows.groupby("gene")["region_of_comparison"].agg(lambda values: sorted(set(values)))
        )
        required_other_regions = sorted(
            [region for region in META_REGIONS if region != region_of_interest]
        )
        qualifying_genes = [
            gene for gene, comparison_regions in successful_other_regions.items()
            if comparison_regions == required_other_regions
        ]

        region_dir = output_dir / sanitize_name(region_of_interest)
        region_dir.mkdir(parents=True, exist_ok=True)

        for gene in sorted(qualifying_genes):
            if gene not in found_genes:
                continue

            row_match = region_rows.loc[region_rows["gene"] == gene].iloc[0]
            region_values = []
            region_labels = []
            region_sample_sizes = []
            for region_name in META_REGIONS:
                values = expr_df.loc[
                    (expr_df["meta_region"] == region_name) & (expr_df[gene] > 0),
                    gene,
                ].astype(float)
                if len(values) == 0:
                    region_values = []
                    break
                region_values.append(values.to_numpy())
                region_labels.append(region_name)
                n_total = int(total_cells_by_region.get(region_name, 0))
                region_sample_sizes.append(f"{region_name} n={len(values)}/{n_total}")

            if not region_values:
                continue

            plot_path = region_dir / f"{sanitize_name(gene)}.png"

            fig, ax = plt.subplots(figsize=(5.2, 4.5))
            box = ax.boxplot(
                region_values,
                patch_artist=True,
                labels=region_labels,
                widths=0.6,
            )
            for patch, region_name in zip(box["boxes"], region_labels):
                patch.set_facecolor(colors.get(region_name, "#999999"))
                patch.set_alpha(0.75)
            for median in box["medians"]:
                median.set_color("#222222")
                median.set_linewidth(1.5)

            ax.set_title(gene, fontsize=12, pad=12)
            ax.set_ylabel("log-normalized expression (expressing cells only)", fontsize=10)
            ax.set_xlabel("meta_region", fontsize=10)
            ax.grid(axis="y", alpha=0.3)
            ax.set_axisbelow(True)
            ax.text(
                0.5,
                1.02,
                (
                    f"{region_of_interest} enriched vs both other regions "
                    f"(FDR < 0.05)"
                ),
                transform=ax.transAxes,
                ha="center",
                va="bottom",
                fontsize=9,
            )
            ax.text(
                0.5,
                0.98,
                (
                    "log-normalized, expressing cells only | "
                    + " | ".join(region_sample_sizes)
                    + f" | group={row_match.combined_leading_gene_group}"
                ),
                transform=ax.transAxes,
                ha="center",
                va="top",
                fontsize=8,
            )

            fig.tight_layout()
            fig.savefig(plot_path, dpi=300, bbox_inches="tight")
            plt.close(fig)
            n_saved_plots += 1

    return {
        "significant_gene_boxplots_dir": str(output_dir),
        "n_significant_gene_boxplots": int(n_saved_plots),
    }


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

    if sym_col not in adata.var.columns:
        raise KeyError(f"Column '{sym_col}' not found in adata.var")
    import scipy.sparse as sp

    adata_f = prepare_expression_adata(
        adata,
        region_col=region_col,
        cell_filter_col=cell_filter_col,
        cell_filter_val=cell_filter_val,
    )

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
    cell_filter_val=None,
    sym_col="Gene",
    expr_threshold=0,
    wilcoxon_region="Forebrain",
    wilcoxon_compare_all_region_pairs=False,
    export_significant_gene_boxplots=False,
    chemistry=DEFAULT_CHEMISTRY,
    chemistry_col=DEFAULT_CHEMISTRY_COL,
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
    adata = filter_adata_by_chemistry(
        adata,
        chemistry=chemistry,
        chemistry_col=chemistry_col,
    )

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
    prepared_adata = prepare_expression_adata(
        adata,
        region_col=region_col,
        cell_filter_col=cell_filter_col,
        cell_filter_val=cell_filter_val,
    )
    all_pairs_per_gene_wilcoxon_summary = build_per_gene_region_wilcoxon_summary(
        prepared_adata,
        summary_info["leading_genes"],
        selected_region=wilcoxon_region,
        sym_col=sym_col,
        leading_gene_flags=leading_gene_flags,
        compare_all_region_pairs=True,
    )
    per_gene_wilcoxon_summary = build_per_gene_region_wilcoxon_summary(
        prepared_adata,
        summary_info["leading_genes"],
        selected_region=wilcoxon_region,
        sym_col=sym_col,
        leading_gene_flags=leading_gene_flags,
        compare_all_region_pairs=wilcoxon_compare_all_region_pairs,
    )

    long_csv_path = output_dir / "gene_expression_region_long_summary.csv"
    final_csv_path = output_dir / "gene_expression_top_region_summary.csv"
    group_wilcoxon_csv_path = output_dir / "gene_expression_group_wilcoxon_summary.csv"
    per_gene_wilcoxon_csv_path = output_dir / "gene_expression_per_gene_wilcoxon_summary.csv"
    region_pair_counts_plot_path = output_dir / "gene_expression_region_pair_significant_gene_counts.png"
    long_summary.to_csv(long_csv_path, index=False)
    final_summary.to_csv(final_csv_path, index=False)
    group_wilcoxon_summary = build_leading_gene_group_wilcoxon_summary(
        final_summary,
        selected_region=wilcoxon_region,
    )
    group_wilcoxon_summary.to_csv(group_wilcoxon_csv_path, index=False)
    per_gene_wilcoxon_summary.to_csv(per_gene_wilcoxon_csv_path, index=False)
    create_region_pair_significant_gene_count_plot(
        all_pairs_per_gene_wilcoxon_summary,
        region_pair_counts_plot_path,
    )
    boxplot_export_result = {
        "significant_gene_boxplots_dir": None,
        "n_significant_gene_boxplots": 0,
    }
    if export_significant_gene_boxplots:
        significant_gene_boxplots_dir = output_dir / "significant_gene_boxplots"
        boxplot_export_result = create_significant_gene_boxplots(
            prepared_adata,
            summary_info["leading_genes"],
            per_gene_wilcoxon_summary=all_pairs_per_gene_wilcoxon_summary,
            output_dir=significant_gene_boxplots_dir,
            selected_region=wilcoxon_region,
            compare_all_region_pairs=wilcoxon_compare_all_region_pairs,
            sym_col=sym_col,
        )

    return {
        "output_dir": str(output_dir),
        "long_summary_csv_path": str(long_csv_path),
        "final_summary_csv_path": str(final_csv_path),
        "group_wilcoxon_summary_csv_path": str(group_wilcoxon_csv_path),
        "per_gene_wilcoxon_summary_csv_path": str(per_gene_wilcoxon_csv_path),
        "region_pair_significant_gene_count_plot_path": str(region_pair_counts_plot_path),
        "gsea_summary_file": str(Path(gsea_summary_file)),
        "condition": summary_info["condition"],
        "condition_column": summary_info["column"],
        "chemistry": chemistry,
        "chemistry_column": chemistry_col,
        "n_input_leading_genes": len(summary_info["leading_genes"]),
        "n_found_genes": int(final_summary.shape[0]),
        "wilcoxon_region": wilcoxon_region,
        "wilcoxon_compare_all_region_pairs": bool(wilcoxon_compare_all_region_pairs),
        "export_significant_gene_boxplots": bool(export_significant_gene_boxplots),
        "wilcoxon_top_vs_second_statistic": wilcoxon_result[
            "wilcoxon_top_vs_second_statistic"
        ],
        "wilcoxon_top_vs_second_pvalue": wilcoxon_result[
            "wilcoxon_top_vs_second_pvalue"
        ],
        "wilcoxon_top_vs_second_fdr_bh_pvalue": wilcoxon_result[
            "wilcoxon_top_vs_second_fdr_bh_pvalue"
        ],
        "wilcoxon_top_vs_second_n_pairs": wilcoxon_result[
            "wilcoxon_top_vs_second_n_pairs"
        ],
        "n_per_gene_wilcoxon_rows": int(per_gene_wilcoxon_summary.shape[0]),
        **boxplot_export_result,
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
        default=None,
        help=(
            "Value in --cell-filter-col to keep before summarizing. "
            "If omitted (default), all cell types are included."
        ),
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
    parser.add_argument(
        "--wilcoxon-compare-all-region-pairs",
        action="store_true",
        help=(
            "For the per-gene Wilcoxon summary CSV, run all ordered meta-region comparisons "
            "(for example, 6 rows per gene for 3 regions) instead of only the selected region "
            "versus the others."
        ),
    )
    parser.add_argument(
        "--export-significant-gene-boxplots",
        action="store_true",
        help=(
            "Create a folder of per-gene boxplots for significant region comparisons, "
            "using only cells with expression > 0 for each gene."
        ),
    )
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
        wilcoxon_compare_all_region_pairs=args.wilcoxon_compare_all_region_pairs,
        export_significant_gene_boxplots=args.export_significant_gene_boxplots,
        chemistry=None if args.chemistry == "None" else args.chemistry,
        chemistry_col=args.chemistry_col,
    )

    print("Saved gene expression summary results:")
    for key, value in result.items():
        print(f"  {key}: {value}")
