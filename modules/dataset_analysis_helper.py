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


def sanitize_filename_component(value):
    """Convert a value to a filesystem-friendly filename component."""
    return (
        str(value)
        .strip()
        .replace("/", "_")
        .replace("\\", "_")
        .replace(" ", "_")
    )


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

    regions_all = adata.obs[region_group_col].to_numpy()
    gene_records = []
    for gene in found_genes:
        gene_matrix = adata[:, gene_to_var[gene]].X
        # AnnData can return dense arrays, scipy sparse matrices, or sparse view wrappers.
        if hasattr(gene_matrix, "toarray"):
            gene_values = gene_matrix.toarray()
        elif hasattr(gene_matrix, "A"):
            gene_values = gene_matrix.A
        else:
            gene_values = np.asarray(gene_matrix)
        gene_values = np.asarray(gene_values).reshape(-1)

        expressed_mask = gene_values > 0
        gene_records.append(
            pd.DataFrame(
                {
                    region_group_col: regions_all[expressed_mask],
                    "Gene": gene,
                    "Expression": gene_values[expressed_mask],
                }
            )
        )

    expr_long = pd.concat(gene_records, ignore_index=True)
    expr_long[region_group_col] = pd.Categorical(
        expr_long[region_group_col],
        categories=["Forebrain", "Midbrain", "Hindbrain"],
        ordered=True,
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    genes_stem = "_".join(found_genes)
    chemistry_label = sanitize_filename_component(chemistry if chemistry is not None else "all_chemistries")
    h5ad_stem = sanitize_filename_component(h5ad_path.stem)
    file_suffix = f"{chemistry_label}_{h5ad_stem}"

    summary_df = (
        expr_long.groupby([region_group_col, "Gene"], observed=False)["Expression"]
        .agg(
            n_cells="size",
            mean_expression="mean",
            std_expression="std",
            median_expression="median",
            min_expression="min",
            max_expression="max",
        )
        .reset_index()
        .sort_values(["Gene", region_group_col])
    )

    csv_path = output_dir / f"{genes_stem}_expression_by_region_summary_{file_suffix}.csv"
    png_path = output_dir / f"{genes_stem}_expression_by_region_{file_suffix}.png"
    summary_df.to_csv(csv_path, index=False)

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
        "summary_csv_path": str(csv_path),
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


# ---------------------------------------------------------------------------
# Expressing-fraction by (Age, Region) × cell-cycle × prolif/diff class
# ---------------------------------------------------------------------------

CELL_BUCKETS = [
    "Proliferating_Cycling",
    "Differentiating_Cycling",
    "Proliferating_NonCycling",
    "Differentiating_NonCycling",
]

DEFAULT_REGION_ORDER_DETAILED = (
    "Forebrain", "Telencephalon", "Diencephalon",
    "Midbrain",
    "Hindbrain", "Cerebellum", "Pons", "Medulla",
)
DEFAULT_EXCLUDE_REGIONS = ("Brain", "Head")

DEFAULT_PROLIFERATING_CLASSES   = ("Radial glia", "Neuronal IPC", "Glioblast")
DEFAULT_DIFFERENTIATING_CLASSES = ("Neuroblast", "Neuron")


def _cell_bucket_series(
    obs, cell_class_col, cell_cycle_score_col, cell_cycle_threshold,
    proliferating_classes, differentiating_classes,
):
    """Return a Series (index = obs.index) with one of CELL_BUCKETS or NA."""
    is_prolif = obs[cell_class_col].isin(list(proliferating_classes))
    is_diff   = obs[cell_class_col].isin(list(differentiating_classes))
    is_cyc    = pd.to_numeric(obs[cell_cycle_score_col], errors="coerce") > float(cell_cycle_threshold)

    bucket = pd.Series(pd.NA, index=obs.index, dtype="object")
    bucket[is_prolif &  is_cyc] = "Proliferating_Cycling"
    bucket[is_diff   &  is_cyc] = "Differentiating_Cycling"
    bucket[is_prolif & ~is_cyc] = "Proliferating_NonCycling"
    bucket[is_diff   & ~is_cyc] = "Differentiating_NonCycling"
    return bucket


def _cells_expressing_enough_genes(adata, genes, sym_col, min_gene_count):
    """Per-cell boolean: does the cell express ≥ `min_gene_count` of the target
    genes (X > 0)? Also returns the found/missing lists and the actual count used."""
    import scipy.sparse as sp

    if sym_col in adata.var.columns:
        sym2var = (
            pd.Series(adata.var_names.values, index=adata.var[sym_col].astype(str))
            .dropna().to_dict()
        )
    else:
        sym2var = {g: g for g in adata.var_names.astype(str)}

    found   = [g for g in genes if g in sym2var]
    missing = [g for g in genes if g not in sym2var]
    if missing:
        preview = ", ".join(missing[:5]) + (" ..." if len(missing) > 5 else "")
        print(f"  [warn] {len(missing):,} gene(s) not found: {preview}")
    if not found:
        raise ValueError("None of the requested genes were found in adata.")

    threshold = max(1, min(int(min_gene_count), len(found)))
    print(f"  Expression criterion: ≥ {threshold}/{len(found)} genes with X > 0")

    varnames = [sym2var[g] for g in found]
    X_sub    = adata[:, varnames].X
    if sp.issparse(X_sub):
        n_nonzero = np.asarray(X_sub.getnnz(axis=1)).ravel()
    else:
        n_nonzero = (np.asarray(X_sub) > 0).sum(axis=1)
    n_nonzero = np.asarray(n_nonzero).ravel()
    expresses_enough = n_nonzero >= threshold
    return expresses_enough, found, missing, threshold


def create_expression_fraction_by_age_region_cellcycle(
    h5ad_path,
    genes,
    output_dir=DEFAULT_OUTPUT_DIR,
    subfolder_name=None,
    chemistry="v3",
    chemistry_col=DEFAULT_CHEMISTRY_COL,
    age_col="Age",
    region_col=DEFAULT_REGION_COL,
    cell_class_col="CellClass",
    cell_cycle_score_col="cell_cycle_score",
    cell_cycle_threshold=0.004,
    sym_col="Gene",
    proliferating_classes=DEFAULT_PROLIFERATING_CLASSES,
    differentiating_classes=DEFAULT_DIFFERENTIATING_CLASSES,
    exclude_regions=DEFAULT_EXCLUDE_REGIONS,
    region_order=DEFAULT_REGION_ORDER_DETAILED,
    plot_metric="Proliferating_Cycling",
    min_gene_fraction=None,
    min_gene_count=None,
):
    """
    Fraction of cells (in each of 4 prolif/diff × cycling buckets) that express
    a sufficient subset of the target gene list, aggregated by (Age, Region).

    Expression criterion (choose one; both None = "any gene" i.e. ≥ 1):
      - `min_gene_fraction`: fraction of the *found* gene list (0-1); rounded up.
        e.g. 0.5 with 20 found genes → cell must express ≥ 10 of them.
      - `min_gene_count`:    absolute integer count of genes required.

    Writes:
      - a long CSV: Age, Region, bucket, n_cells, n_expressing, fraction
      - a wide CSV: Age, Region, one column per bucket (fraction)
      - a plot: x=Age, y=fraction for `plot_metric`, one line per region
        (in `region_order`, `exclude_regions` skipped)
    """
    import math
    adata, h5ad_path = load_filtered_adata(
        h5ad_path, chemistry=chemistry, chemistry_col=chemistry_col,
    )

    for col in (age_col, region_col, cell_class_col, cell_cycle_score_col):
        if col not in adata.obs.columns:
            raise KeyError(
                f"Column '{col}' not found in adata.obs. "
                f"Available: {list(adata.obs.columns)}"
            )

    # 1. bucket each cell
    obs = adata.obs.copy()
    obs["_bucket"] = _cell_bucket_series(
        obs, cell_class_col, cell_cycle_score_col, cell_cycle_threshold,
        proliferating_classes, differentiating_classes,
    )
    print(f"  Cell-bucket counts: {obs['_bucket'].value_counts(dropna=False).to_dict()}")

    # 2. per-cell "expresses enough genes" flag.
    # Resolve the required count from either min_gene_fraction or min_gene_count.
    # If both None → 1 (any gene). Fraction uses ceil so 0.5 * 3 = 2.
    n_input_genes = len(list(genes))
    if min_gene_count is not None:
        required = int(min_gene_count)
    elif min_gene_fraction is not None:
        required = max(1, math.ceil(float(min_gene_fraction) * n_input_genes))
    else:
        required = 1

    expresses_enough, found_genes, missing_genes, threshold_used = (
        _cells_expressing_enough_genes(
            adata, genes, sym_col=sym_col, min_gene_count=required,
        )
    )
    obs["_expresses_any"] = expresses_enough

    # 3. drop excluded regions + cells with no bucket
    excluded_set = set(map(str, exclude_regions or ()))
    keep_mask = (
        obs["_bucket"].notna()
        & ~obs[region_col].astype(str).isin(excluded_set)
    )
    obs_kept = obs.loc[keep_mask].copy()
    print(f"  Kept cells after excluding regions {sorted(excluded_set)}: "
          f"{len(obs_kept):,} / {len(obs):,}")

    # 4. aggregate: per (Age, Region, bucket) → n_cells + n_expressing → fraction
    grouped = (
        obs_kept.groupby([age_col, region_col, "_bucket"], observed=True)
        .agg(n_cells=("_expresses_any", "size"),
             n_expressing=("_expresses_any", "sum"))
        .reset_index()
    )
    grouped["fraction"] = grouped["n_expressing"] / grouped["n_cells"].replace(0, np.nan)
    grouped = grouped.rename(columns={"_bucket": "cell_bucket"})

    # 5. wide format
    wide = grouped.pivot_table(
        index=[age_col, region_col], columns="cell_bucket", values="fraction",
    ).reset_index()
    # Ensure all four bucket columns exist (some buckets may be missing entirely).
    for b in CELL_BUCKETS:
        if b not in wide.columns:
            wide[b] = np.nan
    wide = wide[[age_col, region_col, *CELL_BUCKETS]]

    # 6. output paths
    output_dir = Path(output_dir)
    if subfolder_name:
        output_dir = output_dir / sanitize_filename_component(subfolder_name)
    output_dir.mkdir(parents=True, exist_ok=True)

    chemistry_label = sanitize_filename_component(chemistry if chemistry is not None else "all")
    h5ad_stem       = sanitize_filename_component(h5ad_path.stem)
    threshold_tag   = f"min{threshold_used}of{len(found_genes)}"
    long_csv_path   = output_dir / f"expression_fraction_long_{threshold_tag}_{chemistry_label}_{h5ad_stem}.csv"
    wide_csv_path   = output_dir / f"expression_fraction_wide_{threshold_tag}_{chemistry_label}_{h5ad_stem}.csv"
    plot_png_path   = output_dir / f"expression_fraction_{sanitize_filename_component(plot_metric)}_{threshold_tag}_{chemistry_label}_{h5ad_stem}.png"

    grouped.to_csv(long_csv_path, index=False)
    wide.to_csv(wide_csv_path, index=False)
    print(f"  Saved long CSV: {long_csv_path.name}")
    print(f"  Saved wide CSV: {wide_csv_path.name}")

    # 7. plot
    if plot_metric not in CELL_BUCKETS:
        raise ValueError(f"plot_metric must be one of {CELL_BUCKETS}")
    _plot_fraction_by_age_and_region(
        grouped[grouped["cell_bucket"] == plot_metric],
        age_col=age_col, region_col=region_col,
        region_order=region_order, exclude_regions=exclude_regions,
        title=(f"Fraction of {plot_metric.replace('_', ' ')} cells "
               f"expressing ≥ {threshold_used} of {len(found_genes)} genes"),
        out_path=plot_png_path,
    )
    print(f"  Saved plot:    {plot_png_path.name}")

    return {
        "output_dir":         str(output_dir),
        "long_csv_path":      str(long_csv_path),
        "wide_csv_path":      str(wide_csv_path),
        "plot_png_path":      str(plot_png_path),
        "n_cells_after_filter": int(len(obs_kept)),
        "n_input_genes":      n_input_genes,
        "n_found_genes":      len(found_genes),
        "n_missing_genes":    len(missing_genes),
        "min_gene_count_used": threshold_used,
        "plot_metric":        plot_metric,
    }


def _plot_fraction_by_age_and_region(
    df, age_col, region_col, region_order, exclude_regions, title, out_path,
):
    """Line plot: x=age, y=fraction, one line per region (in `region_order`)."""
    if df.empty:
        print("  [warn] no data to plot for the selected metric.")
        return

    excluded_set = set(map(str, exclude_regions or ()))
    ordered_regions = [r for r in region_order if r not in excluded_set]
    # Add any observed regions not in the fixed order (kept but plotted after).
    observed = df[region_col].dropna().astype(str).unique().tolist()
    extras   = [r for r in observed if r not in ordered_regions and r not in excluded_set]
    plot_regions = ordered_regions + sorted(extras)

    # Convert Age to numeric where possible for a proper x-axis.
    df = df.copy()
    df["_age_num"] = pd.to_numeric(df[age_col], errors="coerce")
    sort_key = "_age_num" if df["_age_num"].notna().any() else age_col

    palette = plt.get_cmap("tab10").colors + plt.get_cmap("Set2").colors
    color_map = {r: palette[i % len(palette)] for i, r in enumerate(plot_regions)}

    fig, ax = plt.subplots(figsize=(9, 5.5))
    for region in plot_regions:
        sub = df[df[region_col].astype(str) == region].sort_values(sort_key)
        if sub.empty:
            continue
        ax.plot(sub[sort_key], sub["fraction"],
                marker="o", linewidth=1.5, color=color_map[region], label=region)

    ax.set_xlabel(age_col)
    ax.set_ylabel("Fraction of cells expressing ≥1 gene")
    ax.set_title(title, fontsize=11)
    ax.set_ylim(-0.02, min(1.02, ax.get_ylim()[1] * 1.02) if ax.get_ylim()[1] > 0 else 1.02)
    ax.grid(alpha=0.3)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False,
              title="Region")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Create helper dataset-analysis plots and tables."
    )
    parser.add_argument(
        "--task",
        default="age_counts",
        choices=["age_counts", "gene_region_boxplot", "expression_fraction_by_age_region_cellcycle"],
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

    # -- expression_fraction_by_age_region_cellcycle -----------------------
    parser.add_argument(
        "--gene-list",
        default=None,
        help=("CSV/text file listing genes for --task "
              "expression_fraction_by_age_region_cellcycle. CSV: uses column "
              "'gene' if present, otherwise the first column; text: one gene "
              "per line."),
    )
    parser.add_argument(
        "--subfolder-name", default=None,
        help="Optional subfolder under --output-dir.",
    )
    parser.add_argument(
        "--cell-class-col", default="CellClass",
        help="Column in adata.obs with the proliferating/differentiating class labels.",
    )
    parser.add_argument(
        "--cell-cycle-score-col", default="cell_cycle_score",
        help="Column in adata.obs with the cell-cycle score.",
    )
    parser.add_argument(
        "--cell-cycle-threshold", type=float, default=0.004,
        help="cell_cycle_score > this = 'cycling' cell.",
    )
    parser.add_argument(
        "--sym-col", default="Gene",
        help="Column in adata.var containing gene symbols.",
    )
    parser.add_argument(
        "--plot-metric", default="Proliferating_Cycling",
        choices=CELL_BUCKETS,
        help="Which cell-bucket fraction to plot on the y-axis.",
    )
    parser.add_argument(
        "--exclude-regions", nargs="*", default=list(DEFAULT_EXCLUDE_REGIONS),
        help="Region values to skip when plotting (default: Brain Head).",
    )
    parser.add_argument(
        "--min-gene-fraction", type=float, default=None,
        help=("Cell counts as 'expressing' if it expresses ≥ ceil(fraction × n_genes) "
              "of the target genes. e.g. 0.5 = at least half. Mutually exclusive "
              "with --min-gene-count. If neither is set, ≥1 gene (any) is used."),
    )
    parser.add_argument(
        "--min-gene-count", type=int, default=None,
        help=("Absolute minimum number of target genes a cell must express to count. "
              "Overrides --min-gene-fraction when both are provided."),
    )
    return parser


def _load_gene_list(path):
    """CSV ('gene' col or first col) or plain-text (one per line)."""
    if path is None:
        return []
    p = Path(path)
    if p.suffix.lower() == ".csv":
        df = pd.read_csv(p)
        col = "gene" if "gene" in df.columns else df.columns[0]
        genes = df[col].dropna().astype(str).str.strip().tolist()
    else:
        with open(p) as f:
            genes = [line.strip() for line in f if line.strip()]
    seen, uniq = set(), []
    for g in genes:
        if g and g not in seen:
            seen.add(g)
            uniq.append(g)
    return uniq


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
    elif args.task == "gene_region_boxplot":
        result = create_gene_expression_boxplot_by_region(
            h5ad_path=args.h5ad_path,
            genes=args.genes or [],
            output_dir=args.output_dir,
            chemistry=chemistry,
            chemistry_col=args.chemistry_col,
            region_col=args.region_col,
            region_group_col=args.region_group_col,
        )
    else:
        genes = args.genes or _load_gene_list(args.gene_list)
        if not genes:
            raise SystemExit(
                "expression_fraction_by_age_region_cellcycle needs a gene list. "
                "Pass --genes A B C  or  --gene-list path/to/genes.csv"
            )
        result = create_expression_fraction_by_age_region_cellcycle(
            h5ad_path=args.h5ad_path,
            genes=genes,
            output_dir=args.output_dir,
            subfolder_name=args.subfolder_name,
            chemistry=chemistry,
            chemistry_col=args.chemistry_col,
            age_col=args.age_col,
            region_col=args.region_col,
            cell_class_col=args.cell_class_col,
            cell_cycle_score_col=args.cell_cycle_score_col,
            cell_cycle_threshold=args.cell_cycle_threshold,
            sym_col=args.sym_col,
            plot_metric=args.plot_metric,
            exclude_regions=args.exclude_regions,
            min_gene_fraction=args.min_gene_fraction,
            min_gene_count=args.min_gene_count,
        )

    print("Saved dataset analysis results:")
    for key, value in result.items():
        print(f"  {key}: {value}")
