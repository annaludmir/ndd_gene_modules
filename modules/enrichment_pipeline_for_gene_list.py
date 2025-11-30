import yaml
from pathlib import Path
import os

from get_gmt import save_to_gmt
from search_enrichment_gsea import run_gsea
from deseq_calculations import main as run_deseq
from create_figs_ges import plot_bar_chart as plot_ges
from create_figs_deseq import plot_bar_chart as plot_deseq


# =============================================================
# CONFIG LOADER — resolves relative paths cleanly
# =============================================================
def load_config(config_path: str):
    config_path = Path(config_path).resolve()
    with open(config_path) as f:
        config = yaml.safe_load(f)

    root = Path(config["ndd_gene_modules_folder_root"]).resolve()

    # Resolve path fields relative to the project folder
    def resolve(p):
        p = Path(p)
        return p if p.is_absolute() else (root / p)

    config["output_folder"] = resolve(config["output_folder"])
    config["gmt_folder"] = resolve(config["gmt_folder"])

    if config["deseq"].get("enabled", False):
        config["deseq"]["pseudobulk_folder"] = resolve(config["deseq"]["pseudobulk_folder"])
        config["deseq"]["gene_names"] = resolve(config["deseq"]["gene_names"])

    return config


# =============================================================
# STEP A — GMT builder (uses output folders from config)
# =============================================================
def create_gmt_file(gmt_folder: Path, gene_list_file: str, make: bool = True) -> Path:
    """
    Ensure a GMT file exists for the given gene list.
    Returns the path to the GMT file.
    """
    gmt_folder = Path(gmt_folder)
    gmt_folder.mkdir(parents=True, exist_ok=True)

    gmt_out = gmt_folder / (Path(gene_list_file).stem + ".gmt")
    if gmt_out.exists():
        print(f"📂 GMT already exists: {gmt_out}")
        return gmt_out

    if not make:
        raise FileNotFoundError("GMT creation disabled, but GMT file does not exist.")

    print(f"🧬 Creating GMT: {gmt_out}")
    save_to_gmt(gene_list_file, gmt_out.stem, gmt_out)
    return gmt_out


# =============================================================
# STEP B — GSEA enrichment run (NO PLOTTING)
# =============================================================
def run_gsea_enrichment(config, gmt_file: Path) -> Path:
    """
    Run GSEA enrichment only (no plots).
    Returns the folder where GSEA results were written.
    """
    if not config["gsea"]["enabled"]:
        print("🔕 GSEA disabled in config")
        return None

    out = config["output_folder"] / f"GSEA_minScore_{config['gsea']['min_ges_score_threshold']}"
    out.mkdir(parents=True, exist_ok=True)

    print("\n🔥 Running GSEA enrichment")
    print(f"→ GMT file:     {gmt_file}")
    print(f"→ Output folder {out}")
    print(f"→ Threshold     {config['gsea']['min_ges_score_threshold']}\n")

    # NOTE: signature of run_gsea should match what you implemented there.
    run_gsea(
        gmt_file=str(gmt_file),
        out_folder=str(out) + "/",
        data_type=config["data_type"],
        min_score=config["gsea"]["min_ges_score_threshold"],
        column_conditions=config.get("column_conditions_for_gsea", {})
    )

    print("✅ GSEA enrichment finished.")
    return out


def run_gsea_plots(results_folder: str | Path):
    """
    Run GSEA plotting only, given a results folder.
    This is intended to be called as a separate Nextflow process.
    """
    results_folder = Path(results_folder)
    print(f"\n📊 Generating GSEA plots from: {results_folder}")
    plot_ges(results_folder=str(results_folder) + "/")


# =============================================================
# STEP C — DESEQ enrichment run (NO PLOTTING)
# =============================================================
def run_deseq_enrichment(config, gene_list_file: str) -> Path | None:
    """
    Run DESeq2 enrichment (no plots).
    Returns the folder where DESeq results were written.
    """
    if not config["deseq"]["enabled"]:
        print("🔕 DESeq disabled in config")
        return None

    out = config["output_folder"] / "DESEQ_results"
    out.mkdir(parents=True, exist_ok=True)

    print("\n🔥 Running DESeq2 differential enrichment\n")

    run_deseq(
        hsg_file=gene_list_file,
        psb_data_folder=str(config["deseq"]["pseudobulk_folder"]),
        gene_names=str(config["deseq"]["gene_names"]),
        out_path=str(out) + "/",
        data_type=config["data_type"],
    )

    print("✅ DESeq enrichment finished.")
    return out


def run_deseq_plots(results_folder: str | Path):
    """
    Run DESeq plotting only, given a results folder.
    This is intended to be called as a separate Nextflow process.
    """
    results_folder = Path(results_folder)
    print(f"\n📊 Generating DESeq plots from: {results_folder}")
    plot_deseq(results_folder=str(results_folder) + "/")


# =============================================================
# MAIN CONVENIENCE FUNCTION (still runs everything in Python)
# =============================================================
def run_gene_list_pipeline(config_path: str, gene_list_file: str):
    """
    Convenience function to run the whole pipeline from Python:
    - load config
    - create GMT
    - run GSEA enrichment (if enabled)
    - run DESeq enrichment (if enabled)
    NO plotting here; use run_gsea_plots / run_deseq_plots separately.
    """
    config = load_config(config_path)

    print("\n===================================================")
    print(" 🔥 GENE LIST PIPELINE — CONFIGURATION SUMMARY")
    print("===================================================")
    print(f"Root directory:  {config['ndd_gene_modules_folder_root']}")
    print(f"Data type:       {config['data_type']}")
    print(f"Output folder:   {config['output_folder']}")
    print(f"GMT folder:      {config['gmt_folder']}")
    print("===================================================\n")

    # ---------- Step A: GMT ----------
    gmt_file = create_gmt_file(
        gmt_folder=config["gmt_folder"],
        gene_list_file=gene_list_file,
    )

    # ---------- Step B: GSEA (NO plots) ----------
    gsea_out = run_gsea_enrichment(config, gmt_file)

    # ---------- Step C: DESEQ (NO plots) ----------
    deseq_out = run_deseq_enrichment(config, gene_list_file)

    print("\n🎉 PIPELINE COMPLETE — enrichment steps finished.")
    print("   For plots, call run_gsea_plots(...) and/or run_deseq_plots(...) separately.\n")
