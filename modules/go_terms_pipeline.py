import argparse
import re
import sys
from pathlib import Path

import pandas as pd
import gseapy as gp

DEFAULT_GENE_SETS = [
    "GO_Biological_Process_2021",
    "GO_Molecular_Function_2021",
    "GO_Cellular_Component_2021",
    "KEGG_2016",
    "KEGG_2021_Human",
    "WikiPathway_2021_Human",
]


def derive_outdir(summary_csv: str | Path) -> Path:
    # Expected path: .../{run_name}/data/enrichment_results/GSEA/GSEA_final_summary.csv
    run_name = Path(summary_csv).parents[3].name
    return Path("results/GO_terms") / f"{run_name}_GO_enrichment"


def get_gene_lists_from_gsea_summary(
    summary_csv: str | Path,
    column_conditions: list[str] | None = None,
) -> dict:
    df = pd.read_csv(summary_csv)

    if column_conditions:
        df = df[df["column"].isin(column_conditions)]
        if df.empty:
            raise ValueError(
                f"No rows matched column_conditions={column_conditions}. "
                f"Available columns: {pd.read_csv(summary_csv)['column'].unique().tolist()}"
            )

    df["condition"] = df["column"] + " - " + df["condition"]

    gene_lists = {}
    for cond, sub in df.groupby("condition"):
        genes = (
            sub["Lead_genes"]
            .dropna()
            .str.split(";")
            .explode()
            .str.strip()
            .unique()
        )
        gene_lists[cond] = list(genes)

    return gene_lists


def run_go_enrichment_per_condition(
    gene_lists: dict,
    outdir: str | Path,
    organism: str = "Human",
    gene_sets: list | None = None,
    cutoff: float = 0.05,
    min_genes: int = 5,
) -> dict:
    if gene_sets is None:
        gene_sets = DEFAULT_GENE_SETS

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    results_by_condition = {}

    for condition, genes in gene_lists.items():
        genes = [g for g in genes if isinstance(g, str) and g.strip()]
        genes = list(dict.fromkeys(genes))

        safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", condition).strip("_")
        cond_outdir = outdir / safe_name
        cond_outdir.mkdir(parents=True, exist_ok=True)

        if len(genes) < min_genes:
            print(f"Skipping {condition} (too few genes: {len(genes)})")
            results_by_condition[condition] = pd.DataFrame()
            continue

        print(f"Running enrichment: {condition} (n_genes={len(genes)})")

        enr = gp.enrichr(
            gene_list=genes,
            organism=organism,
            gene_sets=gene_sets,
            outdir=None,
            cutoff=cutoff,
            no_plot=True,
        )

        res = enr.results if enr is not None and hasattr(enr, "results") else pd.DataFrame()
        results_by_condition[condition] = res

        res.to_csv(cond_outdir / "enrichr_results.csv", index=False)
        pd.Series(genes, name="gene").to_csv(cond_outdir / "leading_genes_used.csv", index=False)

        if not res.empty and "Adjusted P-value" in res.columns:
            sig = res[res["Adjusted P-value"] <= cutoff].copy()
            sig.to_csv(cond_outdir / f"significant_hits_fdr_{cutoff}.csv", index=False)
            if sig.empty:
                print(f"  No significant terms at cutoff={cutoff}")
        else:
            print("  No results returned (or unexpected format).")

    return results_by_condition


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Run GO / pathway enrichment (Enrichr ORA) per condition from a GSEA summary CSV."
    )
    parser.add_argument(
        "--summary-csv",
        required=True,
        help="Path to GSEA_final_summary.csv containing Lead_genes, column, and condition columns.",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help=(
            "Output directory for enrichment results. "
            "If omitted, derived automatically as results/GO_terms/{run_name}_GO_enrichment "
            "where run_name is taken from the summary CSV path."
        ),
    )
    parser.add_argument(
        "--column-conditions",
        nargs="+",
        default=None,
        metavar="COLUMN",
        help=(
            "One or more values from the 'column' field of the summary CSV to include. "
            "E.g. --column-conditions CellCyclePhase Region. "
            "If omitted, all columns are included."
        ),
    )
    parser.add_argument(
        "--organism",
        default="Human",
        help="Organism passed to gseapy Enrichr (default: Human).",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=0.05,
        help="Adjusted p-value cutoff for significant hits (default: 0.05).",
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=5,
        help="Minimum number of genes required to run enrichment for a condition (default: 5).",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()

    outdir = Path(args.outdir) if args.outdir else derive_outdir(args.summary_csv)
    print(f"Output directory: {outdir}")

    gene_lists = get_gene_lists_from_gsea_summary(
        args.summary_csv,
        column_conditions=args.column_conditions,
    )
    print(f"Found {len(gene_lists)} conditions in {args.summary_csv}")

    run_go_enrichment_per_condition(
        gene_lists,
        outdir=outdir,
        organism=args.organism,
        cutoff=args.cutoff,
        min_genes=args.min_genes,
    )

    print(f"\nDone. Results saved to: {outdir}")
    sys.exit(0)
