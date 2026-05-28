"""
GSEA enrichment with Tau pre-filtering.

Drop-in replacement for run_gsea() in search_enrichment_gsea.py.
For each condition:
  1. Load GES scores (prerank metric).
  2. Load global Tau scores for the matching column (from tau_pipeline.py output).
  3. Retain only genes at or above `tau_percentile` (default: 90th → top 10%).
  4. Apply the GES score threshold within that filtered set.
  5. Run prerank GSEA.

Results are written to GSEA_tau_filtered/ instead of GSEA/ so they coexist
with unfiltered runs.
"""

import argparse
import contextlib
import datetime
import shutil
import sys
from pathlib import Path

import gseapy as gp
import pandas as pd
import statsmodels.stats.multitest as smm
import yaml

from search_enrichment_gsea import plot_enhanced_gsea


@contextlib.contextmanager
def _log_to_file(log_path: Path):
    """Tee all stdout (print output) to log_path while keeping it on the terminal."""
    class _Tee:
        def __init__(self, *streams):
            self._streams = streams
        def write(self, s):
            for st in self._streams:
                st.write(s)
        def flush(self):
            for st in self._streams:
                st.flush()
        @property
        def encoding(self):
            return getattr(self._streams[0], "encoding", "utf-8")

    orig = sys.stdout
    with open(log_path, "w", encoding="utf-8") as fh:
        sys.stdout = _Tee(orig, fh)
        try:
            yield
        finally:
            sys.stdout = orig


# ---------------------------------------------------------------------------
# Tau helpers
# ---------------------------------------------------------------------------

def load_tau_scores(tau_scores_dir: Path, column: str) -> pd.DataFrame | None:
    """Load tau_scores_{column}.csv produced by tau_pipeline.py."""
    path = Path(tau_scores_dir) / f"tau_scores_{column}.csv"
    if not path.exists():
        print(f"⚠️ Tau scores not found for column '{column}': {path}")
        return None
    return pd.read_csv(path)


def top_tau_gene_set(tau_df: pd.DataFrame, percentile: float = 90.0) -> set:
    """
    Return the set of genes at or above `percentile` of the Tau distribution.
    E.g. percentile=90 keeps the top 10% most cell-type-specific genes.
    """
    threshold = tau_df["tau"].quantile(percentile / 100.0)
    return set(tau_df.loc[tau_df["tau"] >= threshold, "gene"])


def tau_gene_set_by_score(tau_df: pd.DataFrame, score_cutoff: float) -> set:
    """Return genes with an absolute tau score >= score_cutoff (e.g. 0.5)."""
    return set(tau_df.loc[tau_df["tau"] >= score_cutoff, "gene"])


# ---------------------------------------------------------------------------
# Main enrichment function
# ---------------------------------------------------------------------------

def run_gsea_tau_filtered(
    ges_score_path,
    tau_scores_dir,
    gmt_file,
    column_conditions,
    ges_score_threshold,
    out_folder,
    figs_folder,
    tau_percentile: float = 90.0,
    tau_score_cutoff: float | None = None,
    metadata_dir: Path | None = None,
) -> Path | None:
    """
    Run GSEA prerank enrichment, restricting the gene universe by Tau specificity.

    Two mutually exclusive filter modes:
      - tau_score_cutoff (preferred when set): keep genes with tau >= this absolute value.
      - tau_percentile: keep genes at or above the Nth percentile (default 90 → top 10%).

    Parameters
    ----------
    ges_score_path : path-like
        Folder containing per-condition GES CSVs
        (data/ges_spec_{column}_{condition}.csv).
    tau_scores_dir : path-like
        Folder containing tau_scores_{column}.csv files
        (data/ subfolder of a tau_pipeline.py run).
    gmt_file : str
        Path to the .gmt gene-set file.
    column_conditions : dict
        {column: [condition, ...]} mapping.
    ges_score_threshold : float
        Minimum GES score; genes below this are excluded after the Tau filter.
    out_folder : path-like
        Root folder for enrichment result CSVs.
    figs_folder : path-like
        Root folder for figures.
    tau_percentile : float
        Tau quantile cutoff (default 90 → keep genes with tau >= 90th percentile).
        Ignored when tau_score_cutoff is set.
    tau_score_cutoff : float or None
        Absolute tau score threshold (e.g. 0.5). When set, overrides tau_percentile.

    Returns
    -------
    Path to the summary CSV, or None if no results were produced.
    """
    tau_scores_dir = Path(tau_scores_dir)
    results_folder = Path(out_folder) / "GSEA_tau_filtered"
    fig_dir = Path(figs_folder) / "GSEA_tau_filtered"
    results_folder.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    gmt_name = Path(gmt_file).stem
    print(f"\nGMT set: {gmt_name}")
    if tau_score_cutoff is not None:
        print(f"Tau filter: score cutoff (tau >= {tau_score_cutoff})\n")
    else:
        top_pct = 100.0 - tau_percentile
        print(f"Tau filter: top {top_pct:.0f}% (tau >= {tau_percentile:.0f}th percentile)\n")

    # Parse the GMT to know which genes are in the disease gene set.
    # GMT format: gene_set_name \t description \t gene1 \t gene2 \t ...
    disease_genes: set = set()
    try:
        with open(gmt_file) as _f:
            for line in _f:
                parts = line.strip().split("\t")
                disease_genes.update(parts[2:])  # skip name + description
    except Exception:
        pass  # non-fatal; overlap print will be skipped below

    if metadata_dir is not None:
        metadata_dir = Path(metadata_dir)
        filtered_dir = metadata_dir / "filtered_genes"
        filtered_dir.mkdir(exist_ok=True)
    else:
        filtered_dir = None

    # Filtered-out genes stored per column (disease + expressing-after-tau)
    # and per condition (expressing-after-GES).  Keys: column or "column__condition".
    # Tau gene set is computed per column, so disease and tau-filtered expressing
    # genes are the same for every condition within a column.
    _safe = lambda s: str(s).replace(" ", "_").replace("/", "-")

    _filtered_disease_genes: dict[str, set] = {}   # keyed by column
    _filtered_expr_tau:      dict[str, set] = {}   # keyed by column
    _filtered_expr_ges:      dict[str, set] = {}   # keyed by "column__condition"

    # Cache disease-gene ratio per column: "kept/total"
    _disease_ratio_cache: dict[str, str] = {}

    # Cache tau gene sets per column to avoid repeated I/O and computation
    tau_cache: dict[str, set | None] = {}

    summary_rows = []

    for column, condition_list in column_conditions.items():
        if column not in tau_cache:
            tau_df = load_tau_scores(tau_scores_dir, column)
            if tau_df is None:
                tau_cache[column] = None
                _disease_ratio_cache[column] = "N/A"
            else:
                if tau_score_cutoff is not None:
                    gene_set = tau_gene_set_by_score(tau_df, score_cutoff=tau_score_cutoff)
                    print(
                        f"Column '{column}': {len(gene_set)} genes retained "
                        f"(tau >= {tau_score_cutoff})"
                    )
                else:
                    gene_set = top_tau_gene_set(tau_df, percentile=tau_percentile)
                    print(
                        f"Column '{column}': {len(gene_set)} genes retained "
                        f"(tau >= {tau_percentile:.0f}th percentile)"
                    )
                tau_cache[column] = gene_set
                if disease_genes:
                    n_disease_in = len(disease_genes & gene_set)
                    _disease_ratio_cache[column] = f"{n_disease_in}/{len(disease_genes)}"
                    _filtered_disease_genes[column] = disease_genes - gene_set
                    # Also track which expressing genes are dropped by the tau filter
                    # (shared across all conditions in this column).
                    _filtered_expr_tau[column] = set()
                    print(f"  Disease genes retained: {n_disease_in} / {len(disease_genes)}")
                else:
                    _disease_ratio_cache[column] = "N/A"

        top_tau_set = tau_cache[column]

        for condition in condition_list:
            print(f"\nRunning Tau-filtered GSEA: {column} → {condition}")

            ges_path = Path(ges_score_path) / "data" / f"ges_spec_{column}_{condition}.csv"
            if not ges_path.exists():
                print(f"  ⚠️ Missing GES file: {ges_path} — skipping.")
                continue

            ges_df = pd.read_csv(ges_path)
            n_ges_total   = len(ges_df)
            all_ges_genes = set(ges_df["gene"])

            cond_key = f"{column}__{_safe(condition)}"

            # --- Tau filter ---
            if top_tau_set is not None:
                ges_df = ges_df[ges_df["gene"].isin(top_tau_set)].copy()
                n_ges_after_tau = len(ges_df)
                genes_after_tau = set(ges_df["gene"])
                # Accumulate per-column (tau set is the same for all conditions).
                _filtered_expr_tau.setdefault(column, set()).update(all_ges_genes - genes_after_tau)
                expr_tau_ratio  = f"{n_ges_after_tau}/{n_ges_total}"
                print(f"  Tau filter:  {n_ges_total} → {n_ges_after_tau} genes")
            else:
                print("  ⚠️ No Tau scores available for this column — Tau filter skipped.")
                genes_after_tau = all_ges_genes
                expr_tau_ratio  = f"{n_ges_total}/{n_ges_total}"

            # --- GES threshold ---
            n_before_ges    = len(ges_df)
            ges_df          = ges_df[ges_df["ges_score"] > ges_score_threshold]
            genes_after_ges = set(ges_df["gene"])
            # Per condition: GES scores differ per condition.
            _filtered_expr_ges[cond_key] = genes_after_tau - genes_after_ges
            expr_ges_ratio  = f"{len(genes_after_ges)}/{n_before_ges}"
            print(f"  GES filter:  {n_before_ges} → {len(genes_after_ges)} genes (threshold={ges_score_threshold})")

            if ges_df.empty:
                print(f"  ⚠️ No genes pass both filters — skipping.")
                continue

            ranking = ges_df.set_index("gene")["ges_score"]

            cond_dir = results_folder / column / str(condition)
            raw_dir = cond_dir / "gsea_raw"
            cond_dir.mkdir(parents=True, exist_ok=True)

            try:
                gsea_res = gp.prerank(
                    rnk=ranking,
                    gene_sets=gmt_file,
                    outdir=str(raw_dir),
                    min_size=2,
                    max_size=2500,
                    seed=6,
                )
            except LookupError:
                print("  ⚠️ Not enough genes overlap with gene sets — skipping.")
                gsea_res = None

            if gsea_res is None:
                continue

            out_csv = cond_dir / "gsea_results.csv"
            gsea_res.res2d.to_csv(out_csv)
            print(f"  ✔ Saved results → {out_csv}")

            top = gsea_res.res2d.iloc[0]
            summary_rows.append({
                "column": column,
                "condition": condition,
                "term": top["Term"],
                "NES": top["NES"],
                "NOM p-val": top["NOM p-val"],
                "FDR q-val": top["FDR q-val"],
                "FWER p-val": top["FWER p-val"],
                "Tag %": top["Tag %"],
                "Gene %": top["Gene %"],
                "Lead_genes": top["Lead_genes"],
                "tau_filter_mode":       "score"      if tau_score_cutoff is not None else "percentile",
                "tau_score_cutoff":      tau_score_cutoff,
                "tau_percentile_cutoff": None         if tau_score_cutoff is not None else tau_percentile,
                "disease_genes_after_tau_filtering":                 _disease_ratio_cache.get(column, "N/A"),
                "filtered_out_disease_genes_after_tau_location":     str(filtered_dir / f"disease_after_tau__{_safe(column)}.txt")     if filtered_dir else "",
                "expressing_genes_after_tau_filtering":              expr_tau_ratio,
                "filtered_out_expressing_genes_after_tau":           str(filtered_dir / f"expressing_after_tau__{_safe(column)}.txt")  if filtered_dir else "",
                "expressing_genes_after_ges_threshold_filtering":    expr_ges_ratio,
                "filtered_out_expressing_genes_after_ges_threshold": str(filtered_dir / f"expressing_after_ges__{cond_key}.txt")       if filtered_dir else "",
            })

            term = gsea_res.res2d["Term"].iloc[0]
            if term in gsea_res.results:
                plot_out = plot_enhanced_gsea(gsea_res, term, condition, raw_dir)
                print(f"  ✔ Saved plot → {plot_out}")
            else:
                print(f"  ⚠️ Term '{term}' not in result dict — skipping plot.")

    # Write filtered-out gene lists — one file per column or per condition.
    if filtered_dir is not None:
        def _write_genes(path: Path, genes: set) -> None:
            path.write_text(
                "\n".join(sorted(genes)) + ("\n" if genes else ""),
                encoding="utf-8",
            )
            print(f"  📄 filtered_genes/{path.name}  ({len(genes)} genes)")

        for col, genes in _filtered_disease_genes.items():
            _write_genes(filtered_dir / f"disease_after_tau__{_safe(col)}.txt", genes)

        for col, genes in _filtered_expr_tau.items():
            _write_genes(filtered_dir / f"expressing_after_tau__{_safe(col)}.txt", genes)

        for ck, genes in _filtered_expr_ges.items():
            _write_genes(filtered_dir / f"expressing_after_ges__{ck}.txt", genes)

    if not summary_rows:
        print("\n⚠️ No GSEA results produced — summary not written.")
        return None

    final_summary = pd.DataFrame(summary_rows)

    _, pvals_corr = smm.multipletests(
        final_summary["FDR q-val"].astype(float), method="fdr_bh"
    )[:2]
    final_summary["FDR q-val (BH corrected)"] = pvals_corr

    summary_path = results_folder / "GSEA_tau_filtered_summary.csv"
    final_summary.to_csv(summary_path, index=False)

    print(f"\n📄 Summary written: {summary_path}")
    print("🎉 Tau-filtered GSEA finished.")
    return summary_path


# ---------------------------------------------------------------------------
# Config loader + pipeline orchestrator
# ---------------------------------------------------------------------------

def load_config(config_path: str) -> dict:
    """
    Load and resolve paths from an enrichment-format YAML config.

    Required keys (same as enrichment config):
        run_name, ndd_gene_modules_folder_root, output_folder,
        ges_results_folder, gene_list_path, gmt_folder,
        gsea.min_ges_score_threshold, column_conditions_for_gsea

    Extra keys for tau-filtered runs:
        tau_scores_dir   — path to the data/ folder from a tau_pipeline.py run
        tau_percentile   — Tau quantile cutoff (default 90 → top 10%)
    """
    config_path = Path(config_path).resolve()
    with open(config_path) as f:
        config = yaml.safe_load(f)

    root = Path(config["ndd_gene_modules_folder_root"]).resolve()

    def resolve(p):
        p = Path(p)
        return p if p.is_absolute() else (root / p).resolve()

    for key in ("output_folder", "ges_results_folder", "gene_list_path", "gmt_folder"):
        config[key] = resolve(config[key])

    if "tau_scores_dir" not in config:
        raise ValueError(
            "Config is missing required key 'tau_scores_dir' "
            "(path to the data/ folder of a tau_pipeline.py run)."
        )
    config["tau_scores_dir"] = resolve(config["tau_scores_dir"])
    config["_config_path"] = config_path
    return config


def run_tau_filtered_pipeline(
    config_path: str,
    gene_list_path: str | None = None,
    tau_score_cutoff: float | None = None,
) -> Path | None:
    """
    Full tau-filtered GSEA pipeline driven by a YAML config file.
    Creates a dated output folder, builds the GMT if needed, and
    calls run_gsea_tau_filtered().

    Parameters
    ----------
    config_path : str
        Path to the YAML config file.
    gene_list_path : str, optional
        Override for gene_list_path in the config.
    tau_score_cutoff : float, optional
        Absolute tau score threshold (e.g. 0.5). When set, overrides
        the percentile-based filter from the config.
        Folder name will use 'tauscore{value}' instead of 'tau{pct}'.
    """
    from get_gmt import save_to_gmt

    config = load_config(config_path)

    if gene_list_path is not None:
        config["gene_list_path"] = Path(gene_list_path).resolve()

    run_name = config.get("run_name", Path(config["gene_list_path"]).stem + " tau filtered")
    tau_percentile = float(config.get("tau_percentile", 90.0))
    ges_threshold = config["gsea"]["min_ges_score_threshold"]
    column_conditions = config.get("column_conditions_for_gsea", {})

    date_str = datetime.datetime.now().strftime("%Y%m%d")
    if tau_score_cutoff is not None:
        tau_tag = f"tauscore{tau_score_cutoff}"
    else:
        tau_tag = f"tau{int(tau_percentile)}"
    run_dir = (
        Path(config["output_folder"])
        / f"{run_name}_{tau_tag}_threshold_{ges_threshold}_{date_str}"
    )
    enr_dir = run_dir / "data" / "enrichment_results"
    fig_dir = run_dir / "data" / "enrichment_figures"
    metadata_dir = run_dir / "metadata"
    for d in (enr_dir, fig_dir, metadata_dir):
        d.mkdir(parents=True, exist_ok=True)

    summary_path = None
    with _log_to_file(metadata_dir / "pipeline_output.log"):
        print(f"📋 Log: {metadata_dir / 'pipeline_output.log'}")

        src = Path(config["_config_path"]).resolve()
        dst = (metadata_dir / src.name).resolve()
        if src != dst:
            shutil.copy2(src, dst)

        gmt_folder = Path(config["gmt_folder"])
        gmt_folder.mkdir(parents=True, exist_ok=True)
        gmt_out = gmt_folder / (Path(config["gene_list_path"]).stem + ".gmt")
        if not gmt_out.exists():
            print(f"Creating GMT: {gmt_out}")
            save_to_gmt(str(config["gene_list_path"]), gmt_out)
        else:
            print(f"GMT already exists: {gmt_out}")

        print(f"\n{'='*52}")
        print("  Tau-filtered GSEA Pipeline")
        print(f"{'='*52}")
        print(f"• Run name:        {run_name}")
        print(f"• GES folder:      {config['ges_results_folder']}")
        print(f"• Tau scores dir:  {config['tau_scores_dir']}")
        if tau_score_cutoff is not None:
            print(f"• Tau score cutoff: {tau_score_cutoff}  (absolute)")
        else:
            print(f"• Tau percentile:  {tau_percentile}  (top {100 - tau_percentile:.0f}%)")
        print(f"• GES threshold:   {ges_threshold}")
        print(f"• GMT:             {gmt_out}")
        print(f"• Output:          {run_dir}")
        print(f"{'='*52}\n")

        summary_path = run_gsea_tau_filtered(
            ges_score_path=config["ges_results_folder"],
            tau_scores_dir=config["tau_scores_dir"],
            gmt_file=str(gmt_out),
            column_conditions=column_conditions,
            ges_score_threshold=ges_threshold,
            out_folder=enr_dir,
            figs_folder=fig_dir,
            tau_percentile=tau_percentile,
            tau_score_cutoff=tau_score_cutoff,
            metadata_dir=metadata_dir,
        )

        if summary_path is not None:
            _make_tau_gsea_summary_plots(summary_path, fig_dir, run_name)

    return summary_path


def _make_tau_gsea_summary_plots(summary_path: Path, figs_folder: Path, run_name: str) -> None:
    """Create per-column enrichment bar charts from the tau-filtered GSEA summary CSV."""
    from create_figs_ges_for_presentation import plot_bar_chart as plot_ges

    summary_df = pd.read_csv(summary_path)
    out_dir = figs_folder / "GSEA_tau_filtered"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("\n📊 Generating tau-filtered GSEA summary plots")
    for column_name, df in summary_df.groupby("column"):
        output_path = str(out_dir / f"{column_name}_enrichment.png")
        plot_ges(df, output_path, run_name, column_name)
        print(f"  ✔ Saved enrichment summary plot → {output_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run GSEA enrichment with Tau pre-filtering. "
            "Accepts the same YAML config as the enrichment pipeline, "
            "plus 'tau_scores_dir' and optionally 'tau_percentile'."
        )
    )
    parser.add_argument("config", help="Path to the YAML config file.")
    parser.add_argument(
        "--gene-list",
        default=None,
        metavar="PATH",
        help=(
            "Path to the gene list CSV. Overrides 'gene_list_path' in the config. "
            "Useful when running multiple configs with a single gene list."
        ),
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    run_tau_filtered_pipeline(args.config, gene_list_path=args.gene_list)
    sys.exit(0)
