"""
gene_boxplots_umaps.py

Simple gene-visualization module. Given a gene list and one or more comparison
columns from adata.obs, produce:

  1. One boxplot per gene per comparison column
     (expression on y-axis, column values on x-axis; expressing cells only by default)
  2. One UMAP per gene (colored by expression)
  3. One UMAP per comparison column (colored by category)

Usage:
  python modules/gene_boxplots_umaps.py \
    --h5ad-path data/human_dev.h5ad \
    --gene-list data/genes/my_genes.csv \
    --comparison-columns Region CellClass \
    --chemistry v3 \
    --subfolder-name my_gene_summary

Outputs (under results/additional_analyses/gene_boxplots_umaps/<subfolder>/):
  boxplots/
    <ColumnA>/
      <GENE>.png
    <ColumnB>/
      ...
  umaps/
    gene_expression/<GENE>.png
    columns/<ColumnA>.png
"""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_OUTPUT_DIR   = Path("results/additional_analyses/gene_boxplots_umaps")
DEFAULT_CHEMISTRY    = "v3"
DEFAULT_CHEMISTRY_COL = "Chemistry"


def sanitize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_") or "unnamed"


def load_gene_list(path: str) -> list[str]:
    """CSV with a 'gene' column (or first column), or plain text (one per line)."""
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


def filter_adata_by_chemistry(adata, chemistry, chemistry_col):
    if chemistry is None:
        return adata
    if chemistry_col not in adata.obs.columns:
        print(f"  [warn] chemistry_col '{chemistry_col}' not in adata.obs — skipping filter.")
        return adata
    n_before = adata.n_obs
    adata = adata[adata.obs[chemistry_col].astype(str) == str(chemistry)].copy()
    print(f"  Chemistry filter '{chemistry}': {n_before:,} → {adata.n_obs:,} cells")
    return adata


def normalize_if_needed(adata):
    """log1p-normalize if the data still looks raw (heuristic: max > 50)."""
    import scipy.sparse as sp
    import scanpy as sc

    X_check = adata.X[:100].toarray() if sp.issparse(adata.X) else adata.X[:100]
    if float(X_check.max()) > 50:
        print("  Normalizing (log1p)...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)


def build_symbol_to_var_map(adata, sym_col):
    """Map gene symbol → adata.var_names entry."""
    if sym_col in adata.var.columns:
        return (
            pd.Series(adata.var_names.values, index=adata.var[sym_col].astype(str))
            .dropna()
            .to_dict()
        )
    return {g: g for g in adata.var_names.astype(str)}


def build_expression_dataframe(adata, genes, comparison_columns, sym_col):
    """One row per cell, one column per found gene, plus each comparison column."""
    import scipy.sparse as sp

    sym2var = build_symbol_to_var_map(adata, sym_col)
    found_genes = [g for g in genes if g in sym2var]
    missing = [g for g in genes if g not in sym2var]
    if missing:
        preview = ", ".join(missing[:5]) + (" ..." if len(missing) > 5 else "")
        print(f"  [warn] {len(missing):,} gene(s) not found: {preview}")
    if not found_genes:
        raise ValueError("None of the requested genes were found in adata.")

    varnames = [sym2var[g] for g in found_genes]
    X = adata[:, varnames].X
    if sp.issparse(X):
        X = X.toarray()

    expr_df = pd.DataFrame(X, columns=found_genes)
    for col in comparison_columns:
        expr_df[col] = adata.obs[col].astype(str).to_numpy()
    return expr_df, found_genes


def plot_gene_boxplot_by_column(
    expr_df, gene, column, out_path, expressing_cells_only=True,
):
    """One boxplot for `gene`: x = values of `column`, y = expression."""
    sub = expr_df[[gene, column]].copy()
    if expressing_cells_only:
        sub = sub[sub[gene] > 0]

    categories = sorted(sub[column].dropna().unique())
    values_per_cat = [
        sub.loc[sub[column] == c, gene].astype(float).to_numpy()
        for c in categories
    ]

    if all(len(v) == 0 for v in values_per_cat):
        return False

    fig, ax = plt.subplots(figsize=(max(4.5, 0.6 * len(categories) + 3), 4.5))
    box = ax.boxplot(
        values_per_cat, patch_artist=True, labels=categories, widths=0.55,
    )
    for patch in box["boxes"]:
        patch.set_facecolor("#4c78a8")
        patch.set_alpha(0.75)
    for median in box["medians"]:
        median.set_color("#222222")
        median.set_linewidth(1.5)

    ax.set_title(gene, fontsize=12, pad=10)
    ax.set_ylabel(
        "log-normalized expression"
        + (" (expressing cells only)" if expressing_cells_only else "")
    )
    ax.set_xlabel(column)
    ax.grid(axis="y", alpha=0.3)
    ax.set_axisbelow(True)

    for tick in ax.get_xticklabels():
        tick.set_rotation(30)
        tick.set_ha("right")

    ns_line = " | ".join(f"{c} n={len(v)}" for c, v in zip(categories, values_per_cat))
    ax.text(
        0.5, 1.02, ns_line, transform=ax.transAxes,
        ha="center", va="bottom", fontsize=7, color="#555555",
    )

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def export_boxplots(
    expr_df, found_genes, comparison_columns, out_dir, expressing_cells_only,
):
    n_saved = 0
    for column in comparison_columns:
        col_dir = out_dir / sanitize_name(column)
        col_dir.mkdir(parents=True, exist_ok=True)
        print(f"  Boxplots × {column}  → {col_dir}")
        for gene in found_genes:
            path = col_dir / f"{sanitize_name(gene)}.png"
            if plot_gene_boxplot_by_column(
                expr_df, gene, column, path,
                expressing_cells_only=expressing_cells_only,
            ):
                n_saved += 1
    return n_saved


def export_umaps(adata, genes, sym_col, comparison_columns, out_dir):
    """Per-gene expression UMAP + one UMAP per comparison column."""
    import scanpy as sc

    if "X_umap" not in adata.obsm:
        print("  [skip] UMAPs — adata.obsm['X_umap'] not present in this h5ad.")
        return 0

    sym2var = build_symbol_to_var_map(adata, sym_col)

    gene_dir = out_dir / "gene_expression"
    gene_dir.mkdir(parents=True, exist_ok=True)
    n_saved = 0
    for gene in genes:
        var = sym2var.get(gene)
        if var is None:
            continue
        sc.pl.umap(
            adata, color=var, title=gene, s=5, frameon=False,
            vmax="p99", show=False,
        )
        fig = plt.gcf()
        fig.savefig(gene_dir / f"{sanitize_name(gene)}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        n_saved += 1

    if comparison_columns:
        col_dir = out_dir / "columns"
        col_dir.mkdir(parents=True, exist_ok=True)
        for column in comparison_columns:
            if column not in adata.obs.columns:
                continue
            sc.pl.umap(adata, color=column, s=5, frameon=False, show=False)
            fig = plt.gcf()
            fig.savefig(col_dir / f"{sanitize_name(column)}.png", dpi=300, bbox_inches="tight")
            plt.close(fig)

    return n_saved


def run(
    h5ad_path,
    gene_list_path,
    comparison_columns,
    output_dir=DEFAULT_OUTPUT_DIR,
    subfolder_name=None,
    sym_col="Gene",
    chemistry=DEFAULT_CHEMISTRY,
    chemistry_col=DEFAULT_CHEMISTRY_COL,
    expressing_cells_only=True,
    skip_umaps=False,
):
    output_dir = Path(output_dir)
    if subfolder_name:
        output_dir = output_dir / sanitize_name(subfolder_name)
    else:
        output_dir = output_dir / sanitize_name(Path(gene_list_path).stem)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading gene list: {gene_list_path}")
    genes = load_gene_list(gene_list_path)
    print(f"  {len(genes)} genes")

    print(f"Loading AnnData: {h5ad_path}")
    import scanpy as sc
    adata = sc.read_h5ad(h5ad_path)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    adata = filter_adata_by_chemistry(adata, chemistry=chemistry, chemistry_col=chemistry_col)

    missing_cols = [c for c in comparison_columns if c not in adata.obs.columns]
    if missing_cols:
        raise KeyError(
            f"comparison columns not found in adata.obs: {missing_cols}\n"
            f"Available: {list(adata.obs.columns)}"
        )

    normalize_if_needed(adata)

    expr_df, found_genes = build_expression_dataframe(
        adata, genes, comparison_columns, sym_col=sym_col,
    )
    print(f"  Genes found: {len(found_genes)}/{len(genes)}")

    boxplots_dir = output_dir / "boxplots"
    boxplots_dir.mkdir(parents=True, exist_ok=True)
    n_boxplots = export_boxplots(
        expr_df, found_genes, comparison_columns, boxplots_dir,
        expressing_cells_only=expressing_cells_only,
    )

    n_umaps = 0
    if not skip_umaps:
        umaps_dir = output_dir / "umaps"
        umaps_dir.mkdir(parents=True, exist_ok=True)
        n_umaps = export_umaps(
            adata, found_genes, sym_col=sym_col,
            comparison_columns=comparison_columns, out_dir=umaps_dir,
        )

    print("\nDone:")
    print(f"  output_dir:  {output_dir}")
    print(f"  boxplots:    {n_boxplots}")
    print(f"  umaps:       {n_umaps}")
    print(f"  found genes: {len(found_genes)}/{len(genes)}")


def build_arg_parser():
    p = argparse.ArgumentParser(
        description=(
            "Export per-gene boxplots and UMAPs for a gene list, grouped by "
            "one or more adata.obs comparison columns."
        )
    )
    p.add_argument("--h5ad-path", required=True, help="Input AnnData .h5ad file.")
    p.add_argument(
        "--gene-list", required=True,
        help="CSV with a 'gene' column (or first column), or plain text (one gene per line).",
    )
    p.add_argument(
        "--comparison-columns", nargs="+", required=True,
        help="One or more adata.obs columns to group boxplots by, e.g. Region CellClass.",
    )
    p.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    p.add_argument(
        "--subfolder-name", default=None,
        help="Subfolder under --output-dir. Defaults to the gene-list filename stem.",
    )
    p.add_argument(
        "--sym-col", default="Gene",
        help="Column in adata.var containing gene symbols. Falls back to adata.var_names.",
    )
    p.add_argument(
        "--chemistry", default=DEFAULT_CHEMISTRY,
        help="Chemistry value to keep from adata.obs. Use 'None' to disable.",
    )
    p.add_argument(
        "--chemistry-col", default=DEFAULT_CHEMISTRY_COL,
        help="Column in adata.obs holding chemistry labels.",
    )
    p.add_argument(
        "--include-all-cells", action="store_true",
        help="Include cells with expression == 0 in boxplots (default: expressing cells only).",
    )
    p.add_argument(
        "--skip-umaps", action="store_true",
        help="Only produce boxplots — skip the UMAP step (fast when you don't have a UMAP embedding).",
    )
    return p


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    run(
        h5ad_path=args.h5ad_path,
        gene_list_path=args.gene_list,
        comparison_columns=args.comparison_columns,
        output_dir=args.output_dir,
        subfolder_name=args.subfolder_name,
        sym_col=args.sym_col,
        chemistry=None if args.chemistry == "None" else args.chemistry,
        chemistry_col=args.chemistry_col,
        expressing_cells_only=not args.include_all_cells,
        skip_umaps=args.skip_umaps,
    )
