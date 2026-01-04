import yaml
from pathlib import Path
import os
import shutil
import datetime
import pandas as pd

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
    config["gene_list_path"] = resolve(config["gene_list_path"])
    config["ges_results_folder"] = resolve(config["ges_results_folder"])

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
    save_to_gmt(gene_list_file, gmt_out)
    return gmt_out


# =============================================================
# STEP B — GSEA enrichment run (NO PLOTTING)
# =============================================================
def run_gsea_enrichment(config, gmt_file: Path, enr_results_dir, figs_folder) -> Path:
    """
    Run GSEA enrichment only (no plots).
    Returns the folder where GSEA results were written.
    """
    if config["analysis_mode"] not in ["gsea","both"]:
        print("🔕 GSEA disabled in config")
        return None


    print("\n🔥 Running GSEA enrichment")
    print(f"→ GMT file:            {gmt_file}")
    print(f"→ Output folder        {enr_results_dir}")
    print(f"→ GES score threshold  {config['gsea']['min_ges_score_threshold']}\n")

    full_results_csv_path = run_gsea(
        ges_score_path = config["ges_results_folder"],
        gmt_file = str(gmt_file),
        column_conditions = config.get("column_conditions_for_gsea", {}),
        ges_score_threshold = config["gsea"]["min_ges_score_threshold"],
        out_folder = enr_results_dir,
        figs_folder = figs_folder)

    print("✅ GSEA enrichment finished.")
    return full_results_csv_path


def run_gsea_plots(results_path: str, output_folder, run_name):
    """
    Run GSEA plotting only, given a results folder.
    This is intended to be called as a separate Nextflow process.
    """
    print(f"\n📊 Generating GSEA plots from: {results_path}")
    ges_results_df=pd.read_csv(results_path)
    #divide into different columns and run each file saparatly with the barplot
    different_columns_dfs = [d for _, d in ges_results_df.groupby("column")]
    for df in different_columns_dfs:
      column_name = df.get("column").iat[0]
      output_path = str(output_folder) + f"/GSEA/{column_name}_enrichment.png"
      plot_ges(df, output_path, run_name, column_name)

# =============================================================
# STEP C — DESEQ enrichment run (NO PLOTTING)
# =============================================================
def run_deseq_enrichment(config, gene_list_file: str, enr_results_dir) -> Path | None:
    """
    Run DESeq2 enrichment (no plots).
    Returns the folder where DESeq results were written.
    """
    if config["analysis_mode"] not in ["deseq","both"]:
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
def run_gene_list_pipeline(config_path: str):
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
    print(f"Run name:              {config['run_name']}")
    print(f"Root directory:        {config['ndd_gene_modules_folder_root']}")
    print(f"Output folder:         {config['output_folder']}")
    print(f"Gene list:             {config['gene_list_path']}")
    print(f"GMT folder:            {config['gmt_folder']}")
    print(f"GES Score results:     {config['ges_results_folder']}")
    print(f"Analysis mode:         {config['analysis_mode']}")
    if config['analysis_mode'] in ('gsea', 'both'):
      print(f"GES score threshold:   {config['gsea']['min_ges_score_threshold']}")
      print(f"Column conditions:     {config['column_conditions_for_gsea']}")
    print("===================================================\n")

    ges_results_folder = config['ges_results_folder']
    gene_list_name = config['gene_list_path']
    ges_score_threshold = config['gsea']['min_ges_score_threshold']
  
    # --------------------------
    # Build top-level run folder
    # --------------------------
    date_str = datetime.datetime.now().strftime("%Y%m%d")

    run_name = f"{config['run_name']}_threshold_{ges_score_threshold}_{date_str}"
    run_dir = Path(config['output_folder']) / run_name

    metadata_dir = run_dir / "metadata"
    data_dir = run_dir / "data"
    enr_results_dir = data_dir / "enrichment_results"
    fig_dir = data_dir / "enrichment_figures"
    add_fig_dir = data_dir / "additional_figures"

    # create all necessary folders
    run_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    enr_results_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    add_fig_dir.mkdir(parents=True, exist_ok=True)

    # Copy YAML config to metadata
    src = Path(config_path).resolve()
    dst = (metadata_dir / src.name).resolve()
    
    # Only copy if different paths (or if you want: if different inode)
    if src != dst:
        shutil.copy2(src, dst)
    else:
        print(f"Config already in metadata: {dst} (skipping copy)")
      
    # ---------- Step A: GMT ----------
    gmt_file = create_gmt_file(
        gmt_folder=config["gmt_folder"],
        gene_list_file=config['gene_list_path']
    )

    # ---------- Step B: GSEA ----------
    gsea_out = run_gsea_enrichment(config, gmt_file, enr_results_dir, fig_dir)

    # ---------- Step C: DESEQ ----------
    deseq_out = run_deseq_enrichment(config, config['gene_list_path'], enr_results_dir)

    # ---------- Step D: GSEA Plots ----------
    run_gsea_plots(gsea_out, fig_dir, config['run_name'])
  
    print("\n🎉 PIPELINE COMPLETE — enrichment steps finished.")
    
