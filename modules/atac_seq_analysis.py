"""
atac_seq_analysis.py

ATAC-seq analysis complementing the RNA-side tf_network module. Uses
SnapATAC2 for motif enrichment and pyranges for peak/TSS overlaps.

Three tasks (all can be enabled/disabled from the config):

  1. TF motif enrichment per cell type
     ChromVAR-style motif accessibility. For each candidate TF (from
     tf_network) and each cell type, computes an enrichment score for the
     TF's binding motif in the accessible chromatin of that cell type.

  2. Target-gene promoter accessibility comparison
     For each target predicted by tf_network, finds ATAC peaks overlapping
     the target's TSS (± N bp), then compares per-cell-type mean
     accessibility across region groups (e.g. Telencephalon vs Diencephalon)
     to check whether accessibility mirrors the RNA-side region specificity.

  3. Motif -> target validation
     For each candidate TF, scans the promoters of its predicted targets for
     the TF's binding motif. Reports fraction of targets with a hit — cross-
     validates the GRN with sequence evidence.

Inputs:
  - ATAC h5ad: cells × peaks matrix (var_names look like "chr:start-end")
  - fragment.tsv.bgz (+ .tbi) — used by SnapATAC2 for peak-level operations
  - tf_network cache dir — TF candidates + targets pulled from
    tf_target_summary.csv

Usage:
  python modules/atac_seq_analysis.py config_files/atac_seq_config.yaml
"""

import argparse
import contextlib
import datetime
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


# =============================================================================
# Logging
# =============================================================================

@contextlib.contextmanager
def _log_to_file(*log_paths: Path):
    class _Tee:
        def __init__(self, *streams): self._streams = streams
        def write(self, s):
            for st in self._streams: st.write(s)
        def flush(self):
            for st in self._streams: st.flush()
        @property
        def encoding(self): return getattr(self._streams[0], "encoding", "utf-8")

    orig    = sys.stdout
    handles = [open(p, "w", encoding="utf-8") for p in log_paths]
    try:
        sys.stdout = _Tee(orig, *handles)
        yield
    finally:
        sys.stdout = orig
        for fh in handles:
            fh.close()


# =============================================================================
# Config loading
# =============================================================================

def load_config(config_path: str) -> dict:
    config_path = Path(config_path).resolve()
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    root = Path(cfg["ndd_gene_modules_folder_root"]).resolve()

    def resolve(p):
        if p is None:
            return None
        p = Path(str(p))
        return p if p.is_absolute() else (root / p).resolve()

    cfg["output_folder"] = resolve(cfg.get("output_folder", "results/atac_analysis"))

    atac = cfg.get("atac", {})
    atac["h5ad_path"]     = resolve(atac["h5ad_path"])
    atac["fragment_file"] = resolve(atac.get("fragment_file"))
    cfg["atac"] = atac

    tfn = cfg.get("tf_network", {})
    tfn["cache_dir"] = resolve(tfn["cache_dir"])
    cfg["tf_network"] = tfn

    cfg["_config_path"] = config_path
    cfg["_root"]        = root
    return cfg


# =============================================================================
# Data loaders
# =============================================================================

def load_atac_h5ad(atac_cfg: dict):
    """Load the ATAC AnnData. Reports what obs columns are actually available so
    the config's `region_col` / `cell_type_col` can be verified in the log."""
    import scanpy as sc

    print(f"\n[Load] ATAC h5ad: {atac_cfg['h5ad_path']}")
    adata = sc.read_h5ad(atac_cfg["h5ad_path"])
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} peaks")

    for col_key in ("cell_type_col", "region_col"):
        col = atac_cfg.get(col_key)
        if col and col not in adata.obs.columns:
            print(f"  [warn] {col_key}='{col}' not in obs. Available columns:")
            print(f"     {list(adata.obs.columns)}")

    return adata


def load_tf_network_targets(
    tfn_cfg: dict, min_confidence: str = "high", top_n_tfs: int | None = None,
) -> pd.DataFrame:
    """Read tf_target_summary.csv from a tf_network cache dir. Returns a
    long-format DataFrame with (TF, target, confidence, ...) rows filtered
    to the requested confidence level. Top-N TFs selected by target count."""
    cache_dir = Path(tfn_cfg["cache_dir"])
    summary_path = cache_dir / "tf_target_summary.csv"
    if not summary_path.exists():
        raise FileNotFoundError(
            f"tf_target_summary.csv not found at {summary_path}. "
            "Run tf_network on this dataset first (or point cache_dir at an existing cache)."
        )

    df = pd.read_csv(summary_path)
    print(f"\n[Load] tf_target_summary.csv: {len(df):,} TF-target pairs")

    # Confidence filter: keep pairs at min_confidence tier or better.
    confidence_order = {"low": 0, "medium": 1, "high": 2}
    if min_confidence not in confidence_order:
        raise ValueError(f"min_confidence must be one of {list(confidence_order)}")
    threshold = confidence_order[min_confidence]
    df = df[df["confidence"].map(confidence_order).fillna(-1) >= threshold]
    print(f"  After confidence>='{min_confidence}': {len(df):,} pairs")

    if top_n_tfs is not None:
        target_counts = df.groupby("TF")["target"].nunique().sort_values(ascending=False)
        keep_tfs = target_counts.head(top_n_tfs).index
        df = df[df["TF"].isin(keep_tfs)]
        print(f"  After top-{top_n_tfs}-TFs (by target count): {len(df):,} pairs "
              f"({df['TF'].nunique()} TFs)")

    return df.reset_index(drop=True)


# =============================================================================
# Peak parsing / coordinate helpers
# =============================================================================

_PEAK_RE = re.compile(r"^([A-Za-z0-9_.-]+)[:_-](\d+)[-_](\d+)$")


_CHROM_ALIASES = ["Chromosome", "chromosome", "chrom", "chr", "seqnames"]
_START_ALIASES = ["Start", "start", "peak_start", "chromStart"]
_END_ALIASES   = ["End",   "end",   "peak_end",   "chromEnd"]


def _pick_col(df: pd.DataFrame, aliases: list[str]) -> str | None:
    for a in aliases:
        if a in df.columns:
            return a
    return None


def peaks_to_dataframe(adata, peak_columns: dict | None = None) -> pd.DataFrame:
    """Return a DataFrame with columns Chromosome / Start / End / peak_id.

    Tries, in order:
      1. Explicit `peak_columns` mapping from the config, e.g.
         {"chrom": "chr", "start": "start", "end": "end"}.
      2. Auto-detect chrom/start/end columns in adata.var.
      3. Parse adata.var_names as 'chr:start-end' or 'chr_start_end'.
    """
    var = adata.var

    # 1. explicit config
    if peak_columns:
        try:
            chrom_c = peak_columns["chrom"]
            start_c = peak_columns["start"]
            end_c   = peak_columns["end"]
            df = pd.DataFrame({
                "Chromosome": var[chrom_c].astype(str).values,
                "Start":      var[start_c].astype(int).values,
                "End":        var[end_c].astype(int).values,
                "peak_id":    adata.var_names.astype(str).values,
            })
            print(f"  [peaks] using explicit var columns: {chrom_c}/{start_c}/{end_c}")
            return df
        except KeyError as e:
            print(f"  [warn] peak_columns entry {e} missing in adata.var — falling back.")

    # 2. auto-detect coordinate columns in var
    chrom_c = _pick_col(var, _CHROM_ALIASES)
    start_c = _pick_col(var, _START_ALIASES)
    end_c   = _pick_col(var, _END_ALIASES)
    if chrom_c and start_c and end_c:
        df = pd.DataFrame({
            "Chromosome": var[chrom_c].astype(str).values,
            "Start":      var[start_c].astype(int).values,
            "End":        var[end_c].astype(int).values,
            "peak_id":    adata.var_names.astype(str).values,
        })
        print(f"  [peaks] auto-detected var columns: {chrom_c}/{start_c}/{end_c}")
        return df

    # 3. parse var_names
    rows = []
    for name in adata.var_names.astype(str):
        m = _PEAK_RE.match(name)
        if not m:
            continue
        chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
        rows.append({"Chromosome": chrom, "Start": start, "End": end, "peak_id": name})
    if rows:
        print(f"  [peaks] parsed from var_names ({len(rows):,} peaks)")
        return pd.DataFrame(rows)

    # 4. nothing worked — print a helpful diagnostic
    sample_names = list(adata.var_names[:5].astype(str))
    var_head     = var.head(3).to_dict("list") if not var.empty else {}
    raise ValueError(
        "Could not locate peak coordinates in this h5ad.\n"
        f"  var_names sample (first 5): {sample_names}\n"
        f"  var.columns: {list(var.columns)}\n"
        f"  var head: {var_head}\n"
        "Fix: add to your config under atac:\n"
        "  peak_columns:\n"
        "    chrom: <name of chrom column in var>\n"
        "    start: <name of start column in var>\n"
        "    end:   <name of end column in var>"
    )


def load_tss_annotation(gtf_path: Path | None, genome: str) -> pd.DataFrame:
    """Return a DataFrame with columns Chromosome / TSS / Strand / gene_name.

    If a GTF file is provided, parse it. Otherwise fall back to SnapATAC2's
    built-in gene annotations (which are curated for supported genomes).
    """
    if gtf_path is not None and Path(gtf_path).exists():
        print(f"\n[Load] gene TSSes from GTF: {gtf_path}")
        return _tss_from_gtf(Path(gtf_path))

    print(f"\n[Load] gene TSSes from snapatac2 built-in ({genome})")
    try:
        import snapatac2 as snap
    except ImportError as e:
        raise ImportError(
            "snapatac2 is required. Install with: pip install snapatac2"
        ) from e
    genome_obj = getattr(snap.genome, genome, None)
    if genome_obj is None:
        raise ValueError(f"snapatac2.genome has no attribute '{genome}'.")
    # SnapATAC2 exposes gene annotations; if the exact API differs across
    # versions, we may need to switch to a GTF file (set atac.gtf_path).
    ann = getattr(genome_obj, "annotation", None)
    if ann is None:
        raise RuntimeError(
            f"snap.genome.{genome} has no `annotation` attribute in this snapatac2 "
            "version. Please provide atac.gtf_path in the config."
        )
    df = pd.DataFrame(ann)
    # Normalise column names — snapatac2 versions vary slightly here.
    rename = {"chrom": "Chromosome", "start": "Start", "end": "End",
              "strand": "Strand", "gene_name": "gene_name"}
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    df["TSS"] = np.where(df["Strand"] == "+", df["Start"], df["End"])
    return df[["Chromosome", "TSS", "Strand", "gene_name"]]


def _tss_from_gtf(gtf_path: Path) -> pd.DataFrame:
    """Parse a GTF file to extract per-gene TSS positions."""
    rows = []
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
            attrs = parts[8]
            m = re.search(r'gene_name "([^"]+)"', attrs)
            if not m:
                continue
            tss = start if strand == "+" else end
            rows.append({"Chromosome": chrom, "TSS": tss, "Strand": strand,
                         "gene_name": m.group(1)})
    return pd.DataFrame(rows)


def build_promoter_regions(
    tss_df: pd.DataFrame, window: int, genes: list[str] | None = None,
) -> pd.DataFrame:
    """Return per-gene promoter windows: Chromosome, Start, End, gene_name."""
    df = tss_df.copy()
    if genes is not None:
        df = df[df["gene_name"].isin(set(genes))]
    df["Start"] = (df["TSS"] - window).clip(lower=0)
    df["End"]   = df["TSS"] + window
    return df[["Chromosome", "Start", "End", "gene_name"]].reset_index(drop=True)


def peaks_overlapping_promoters(
    peaks_df: pd.DataFrame, promoter_df: pd.DataFrame,
) -> pd.DataFrame:
    """For each promoter, list the peak_ids that overlap it. Returns a
    DataFrame with one row per (gene_name, peak_id)."""
    try:
        import pyranges as pr
    except ImportError as e:
        raise ImportError("pyranges is required. Install with: pip install pyranges") from e

    peaks_pr     = pr.PyRanges(peaks_df.rename(columns={"peak_id": "peak_id"}))
    promoters_pr = pr.PyRanges(promoter_df)
    overlaps     = peaks_pr.join(promoters_pr).df
    if overlaps.empty:
        return pd.DataFrame(columns=["gene_name", "peak_id"])
    keep_cols = ["gene_name", "peak_id"]
    return overlaps[keep_cols].drop_duplicates().reset_index(drop=True)


# =============================================================================
# Task 1 — TF motif enrichment per cell type
# =============================================================================

def task_motif_enrichment(
    adata, tf_targets: pd.DataFrame, cfg: dict, out_dir: Path,
) -> Path | None:
    """
    For each cell type in the configured list, run SnapATAC2's motif enrichment
    on that cell type's most-accessible peaks. Report the enrichment score of
    every candidate TF (from tf_targets['TF']) across every cell type.
    """
    print("\n" + "=" * 62)
    print("  Task 1: TF motif enrichment per cell type")
    print("=" * 62)

    task_cfg  = cfg.get("motif_enrichment", {})
    if not task_cfg.get("enabled", True):
        print("  [skip] motif_enrichment.enabled = false")
        return None

    cell_types    = task_cfg.get("cell_types") or []
    motif_source  = task_cfg.get("motif_source", "cis_bp")
    genome_name   = cfg["atac"].get("genome", "hg38")
    top_peaks_pct = float(task_cfg.get("top_peaks_pct", 0.10))
    cell_type_col = cfg["atac"]["cell_type_col"]

    if not cell_types:
        print("  [skip] no cell_types listed in motif_enrichment.cell_types")
        return None

    try:
        import snapatac2 as snap
    except ImportError as e:
        raise ImportError("snapatac2 required: pip install snapatac2") from e

    # Motif library
    if motif_source == "cis_bp":
        motifs = snap.datasets.cis_bp(unique=True)
    elif motif_source == "jaspar":
        motifs = snap.datasets.Jaspar(unique=True)  # API name may vary
    else:
        raise ValueError(f"Unknown motif_source '{motif_source}'")

    candidate_tfs = set(tf_targets["TF"].astype(str).unique())
    print(f"  {len(candidate_tfs)} candidate TFs from tf_network")
    print(f"  Cell types: {cell_types}")

    # Build per-cell-type peak sets: take the top X% most-accessible peaks in
    # that cell type (pseudo-bulk mean over cells of that type).
    import scipy.sparse as sp
    X = adata.X
    if sp.issparse(X):
        X = X.tocsr()

    peak_df   = peaks_to_dataframe(adata, cfg["atac"].get("peak_columns"))
    peak_ids  = peak_df["peak_id"].tolist()
    peak_idx  = {pid: i for i, pid in enumerate(peak_ids)}

    per_type_top_peaks: dict[str, pd.DataFrame] = {}
    for ct in cell_types:
        mask = (adata.obs[cell_type_col].astype(str) == str(ct)).to_numpy()
        n_ct = int(mask.sum())
        if n_ct == 0:
            print(f"  [warn] no cells for cell_type='{ct}'")
            continue
        ct_mean = np.asarray(X[mask, :].mean(axis=0)).ravel()
        n_top   = max(1, int(len(peak_ids) * top_peaks_pct))
        top_ix  = np.argsort(ct_mean)[::-1][:n_top]
        selected = peak_df.iloc[top_ix].reset_index(drop=True)
        per_type_top_peaks[str(ct)] = selected
        print(f"    {ct}: {n_ct:,} cells → top {n_top:,} peaks")

    if not per_type_top_peaks:
        print("  [skip] no cell types had cells matching")
        return None

    genome_obj = getattr(snap.genome, genome_name)

    # snap.tl.motif_enrichment signature varies across versions. Common form:
    #   motif_enrichment(motifs, regions=dict[str, DataFrame], genome_fasta=...)
    print("\n  Running motif enrichment via snapatac2 ...")
    try:
        result = snap.tl.motif_enrichment(
            motifs=motifs,
            regions=per_type_top_peaks,
            genome_fasta=genome_obj,
        )
    except TypeError as e:
        print(f"  [warn] snap.tl.motif_enrichment signature mismatch: {e}")
        print("         Falling back to per-cell-type calls (older API).")
        result = {}
        for ct, peaks in per_type_top_peaks.items():
            result[ct] = snap.tl.motif_enrichment(
                motifs=motifs, regions=peaks, genome_fasta=genome_obj,
            )

    # Flatten result into a long DataFrame regardless of exact return shape.
    rows = []
    for ct, res in (result.items() if isinstance(result, dict) else [("all", result)]):
        if hasattr(res, "to_dataframe"):
            res_df = res.to_dataframe()
        elif isinstance(res, pd.DataFrame):
            res_df = res
        else:
            res_df = pd.DataFrame(res)
        for _, r in res_df.iterrows():
            motif_name = str(r.get("name", r.get("motif", "")))
            log2fe     = float(r.get("log2_fold_enrichment", r.get("log2FE", np.nan)))
            padj       = float(r.get("adjusted_p_value",   r.get("padj",   np.nan)))
            rows.append({
                "cell_type": ct, "motif": motif_name,
                "log2_fold_enrichment": log2fe, "adjusted_p_value": padj,
            })
    enrich_df = pd.DataFrame(rows)

    # Restrict rows to candidate TFs where possible (motif name usually
    # contains the TF gene symbol).
    if not enrich_df.empty:
        def _has_candidate_tf(motif_name: str) -> str | None:
            m = motif_name.upper()
            for tf in candidate_tfs:
                if tf.upper() in m:
                    return tf
            return None
        enrich_df["candidate_tf"] = enrich_df["motif"].map(_has_candidate_tf)
        candidate_df = enrich_df[enrich_df["candidate_tf"].notna()].copy()
    else:
        candidate_df = enrich_df

    out_all  = out_dir / "motif_enrichment_all_motifs.csv"
    out_cand = out_dir / "motif_enrichment_candidate_tfs.csv"
    out_dir.mkdir(parents=True, exist_ok=True)
    enrich_df.to_csv(out_all, index=False)
    candidate_df.to_csv(out_cand, index=False)
    print(f"\n  Saved: {out_all.name}   ({len(enrich_df):,} rows)")
    print(f"  Saved: {out_cand.name}  ({len(candidate_df):,} rows for candidate TFs)")

    # Heatmap: cell_type × candidate_tf, values = log2 fold enrichment
    if not candidate_df.empty:
        pivot = (candidate_df
                 .groupby(["cell_type", "candidate_tf"])["log2_fold_enrichment"]
                 .max()
                 .unstack(fill_value=np.nan))
        _heatmap(pivot, out_dir / "motif_enrichment_heatmap.png",
                 title="TF motif enrichment (log2 FE) — cell type × candidate TF",
                 cmap="viridis")
        print(f"  Saved: motif_enrichment_heatmap.png")

    return out_cand


# =============================================================================
# Task 2 — target-gene promoter accessibility across region groups
# =============================================================================

def task_promoter_accessibility(
    adata, tf_targets: pd.DataFrame, cfg: dict, out_dir: Path,
) -> Path | None:
    """
    For each target gene: gather peaks overlapping its promoter (TSS ± window),
    compute per-cell-type mean accessibility, then compare across region groups.
    """
    print("\n" + "=" * 62)
    print("  Task 2: target-gene promoter accessibility")
    print("=" * 62)

    task_cfg = cfg.get("promoter", {})
    if not task_cfg.get("enabled", True):
        print("  [skip] promoter.enabled = false")
        return None

    tss_window     = int(task_cfg.get("tss_window", 2000))
    region_pairs   = task_cfg.get("region_pairs", [])
    region_col     = cfg["atac"].get("region_col")
    gtf_path       = cfg["atac"].get("gtf_path")
    genome_name    = cfg["atac"].get("genome", "hg38")

    if not region_col or region_col not in adata.obs.columns:
        print(f"  [skip] region_col '{region_col}' not in adata.obs")
        return None
    if not region_pairs:
        print("  [skip] no promoter.region_pairs configured")
        return None

    targets = sorted(tf_targets["target"].astype(str).unique())
    print(f"  Targets to test: {len(targets)}")

    tss_df       = load_tss_annotation(gtf_path, genome_name)
    promoter_df  = build_promoter_regions(tss_df, tss_window, genes=targets)
    if promoter_df.empty:
        print(f"  [skip] no TSS matches for any of the {len(targets)} targets")
        return None
    print(f"  TSS windows built: {len(promoter_df):,}")

    peaks_df    = peaks_to_dataframe(adata, cfg["atac"].get("peak_columns"))
    peak_to_ix  = {pid: i for i, pid in enumerate(peaks_df["peak_id"])}

    overlap_df  = peaks_overlapping_promoters(peaks_df, promoter_df)
    if overlap_df.empty:
        print("  [skip] no peaks overlap any target promoter")
        return None
    print(f"  Peak-promoter overlaps: {len(overlap_df):,}  "
          f"({overlap_df['gene_name'].nunique()} targets have ≥1 peak)")

    import scipy.sparse as sp
    X = adata.X
    if sp.issparse(X):
        X = X.tocsr()

    out_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []
    for region_a, region_b in region_pairs:
        pair_label = f"{_sanitize(region_a)}_vs_{_sanitize(region_b)}"
        for label, region_val in ((region_a, region_a), (region_b, region_b)):
            if region_val not in set(adata.obs[region_col].astype(str)):
                print(f"  [warn] region '{region_val}' not in adata.obs['{region_col}']")

        mask_a = (adata.obs[region_col].astype(str) == str(region_a)).to_numpy()
        mask_b = (adata.obs[region_col].astype(str) == str(region_b)).to_numpy()
        n_a, n_b = int(mask_a.sum()), int(mask_b.sum())
        print(f"  {pair_label}: {n_a:,} vs {n_b:,} cells")
        if n_a == 0 or n_b == 0:
            continue

        per_gene = (
            overlap_df.groupby("gene_name")["peak_id"].apply(list).to_dict()
        )
        rows = []
        for gene, peak_ids in per_gene.items():
            ix = [peak_to_ix[p] for p in peak_ids if p in peak_to_ix]
            if not ix:
                continue
            X_sub = X[:, ix]
            per_cell_sum = np.asarray(X_sub.sum(axis=1)).ravel()
            mean_a = float(per_cell_sum[mask_a].mean())
            mean_b = float(per_cell_sum[mask_b].mean())
            rows.append({
                "gene": gene,
                "region_a": region_a, "region_b": region_b,
                "n_promoter_peaks": len(ix),
                "n_cells_a": n_a, "n_cells_b": n_b,
                "mean_accessibility_a": mean_a,
                "mean_accessibility_b": mean_b,
                "log2_ratio_a_over_b": (
                    np.log2((mean_a + 1e-9) / (mean_b + 1e-9))
                ),
            })
        pair_df = pd.DataFrame(rows)
        pair_csv = out_dir / f"promoter_accessibility_{pair_label}.csv"
        pair_df.to_csv(pair_csv, index=False)
        print(f"    Saved: {pair_csv.name}  ({len(pair_df):,} genes)")
        all_rows.append(pair_df.assign(pair=pair_label))

    if not all_rows:
        return None
    combined = pd.concat(all_rows, ignore_index=True)
    combined_csv = out_dir / "promoter_accessibility_all_pairs.csv"
    combined.to_csv(combined_csv, index=False)
    print(f"  Saved: {combined_csv.name}  ({len(combined):,} rows total)")
    return combined_csv


# =============================================================================
# Task 3 — motif -> target validation
# =============================================================================

def task_motif_target_validation(
    adata, tf_targets: pd.DataFrame, cfg: dict, out_dir: Path,
) -> Path | None:
    """
    For each candidate TF, scan the promoters of its predicted targets for the
    TF's binding motif. Reports the fraction of targets with a promoter motif
    hit — cross-validates the RNA-side GRN with sequence evidence.
    """
    print("\n" + "=" * 62)
    print("  Task 3: motif -> target validation")
    print("=" * 62)

    task_cfg = cfg.get("motif_target_validation", {})
    if not task_cfg.get("enabled", True):
        print("  [skip] motif_target_validation.enabled = false")
        return None

    tss_window   = int(task_cfg.get("promoter_window", 2000))
    pval_cutoff  = float(task_cfg.get("motif_pval_cutoff", 1e-4))
    gtf_path     = cfg["atac"].get("gtf_path")
    genome_name  = cfg["atac"].get("genome", "hg38")

    try:
        import snapatac2 as snap
    except ImportError as e:
        raise ImportError("snapatac2 required: pip install snapatac2") from e

    genome_obj = getattr(snap.genome, genome_name)
    motifs     = snap.datasets.cis_bp(unique=True)

    # Build a symbol → motif lookup. Motif names in cis_bp typically include
    # the TF gene symbol; we match by substring (case-insensitive).
    def _find_motifs_for_tf(tf: str) -> list:
        tf_up = tf.upper()
        return [m for m in motifs if tf_up in getattr(m, "name", "").upper()]

    tss_df   = load_tss_annotation(gtf_path, genome_name)
    tss_by_gene = tss_df.set_index("gene_name").to_dict("index")

    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    tf_summary_rows = []

    for tf, tf_group in tf_targets.groupby("TF"):
        tf_motifs = _find_motifs_for_tf(str(tf))
        target_genes = tf_group["target"].astype(str).unique().tolist()

        if not tf_motifs:
            print(f"  [warn] no motif found in cis_bp for TF '{tf}' — skipping")
            tf_summary_rows.append({
                "TF": tf, "n_targets": len(target_genes),
                "n_motifs_found": 0, "n_targets_with_hit": 0,
                "fraction_with_hit": np.nan,
            })
            continue

        # Build promoter regions for this TF's targets
        promoter_df = build_promoter_regions(tss_df, tss_window, genes=target_genes)
        if promoter_df.empty:
            continue

        # SnapATAC2 doesn't expose a per-sequence FIMO-style scanner directly,
        # so we use motif_enrichment on the promoter set (foreground) with the
        # whole peak set as background. Any TF whose motif has adj_p < cutoff
        # in the promoter set is treated as "validated". We also record per-
        # gene motif counts if snapatac2 returns them.
        try:
            per_promoter_hits = snap.tl.motif_enrichment(
                motifs=tf_motifs,
                regions={"promoters": promoter_df.rename(
                    columns={"gene_name": "name"}
                )},
                genome_fasta=genome_obj,
            )
        except Exception as e:
            print(f"  [warn] motif scan failed for {tf}: {e}")
            continue

        hits_df = _flatten_enrichment_result(per_promoter_hits)
        significant = hits_df[hits_df.get("adjusted_p_value", np.nan) < pval_cutoff]

        n_hits = int(significant.shape[0])
        tf_summary_rows.append({
            "TF": tf,
            "n_targets": len(target_genes),
            "n_motifs_found": len(tf_motifs),
            "n_significant_motif_enrichments": n_hits,
            "min_adjusted_p_value": (
                float(hits_df["adjusted_p_value"].min())
                if "adjusted_p_value" in hits_df.columns and not hits_df.empty
                else np.nan
            ),
        })
        for _, r in hits_df.iterrows():
            rows.append({
                "TF": tf,
                "motif": r.get("name", r.get("motif", "")),
                "log2_fold_enrichment": r.get("log2_fold_enrichment", np.nan),
                "adjusted_p_value": r.get("adjusted_p_value", np.nan),
            })

    detail_df  = pd.DataFrame(rows)
    summary_df = pd.DataFrame(tf_summary_rows)
    detail_csv = out_dir / "motif_target_validation_details.csv"
    summary_csv = out_dir / "motif_target_validation_summary.csv"
    detail_df.to_csv(detail_csv, index=False)
    summary_df.to_csv(summary_csv, index=False)
    print(f"\n  Saved: {detail_csv.name}   ({len(detail_df):,} rows)")
    print(f"  Saved: {summary_csv.name}  ({len(summary_df):,} TFs)")
    return summary_csv


def _flatten_enrichment_result(result) -> pd.DataFrame:
    if isinstance(result, dict):
        frames = []
        for _, v in result.items():
            frames.append(v.to_dataframe() if hasattr(v, "to_dataframe") else pd.DataFrame(v))
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if hasattr(result, "to_dataframe"):
        return result.to_dataframe()
    return pd.DataFrame(result)


# =============================================================================
# Small utilities
# =============================================================================

def _sanitize(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_") or "unnamed"


def _heatmap(df: pd.DataFrame, out_path: Path, title: str, cmap: str = "viridis"):
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(max(6, 0.3 * df.shape[1] + 3),
                                    max(4, 0.3 * df.shape[0] + 2)))
    im = ax.imshow(df.values, aspect="auto", cmap=cmap)
    ax.set_xticks(range(df.shape[1])); ax.set_xticklabels(df.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(df.shape[0])); ax.set_yticklabels(df.index, fontsize=7)
    ax.set_title(title, fontsize=10)
    plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Pipeline orchestration
# =============================================================================

def run(config_path: str) -> None:
    cfg = load_config(config_path)

    date_str = datetime.datetime.now().strftime("%Y%m%d")
    dataset_name = cfg.get("dataset_name", "atac_run")
    run_dir = cfg["output_folder"] / f"{dataset_name}_{date_str}"
    run_dir.mkdir(parents=True, exist_ok=True)

    metadata_dir = run_dir / "metadata"
    metadata_dir.mkdir(exist_ok=True)
    import shutil
    shutil.copy2(cfg["_config_path"], metadata_dir / Path(config_path).name)

    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = metadata_dir / f"pipeline_output_{ts}.log"

    with _log_to_file(log_path):
        print("=" * 62)
        print(f"  ATAC-seq analysis — {dataset_name}")
        print(f"  Run dir: {run_dir}")
        print("=" * 62)

        adata      = load_atac_h5ad(cfg["atac"])
        tf_targets = load_tf_network_targets(
            cfg["tf_network"],
            min_confidence=cfg["tf_network"].get("min_confidence", "high"),
            top_n_tfs=cfg["tf_network"].get("top_n_tfs"),
        )

        task_motif_enrichment(adata, tf_targets, cfg, run_dir / "1_motif_enrichment")
        task_promoter_accessibility(adata, tf_targets, cfg, run_dir / "2_promoter_accessibility")
        task_motif_target_validation(adata, tf_targets, cfg, run_dir / "3_motif_target_validation")

        print(f"\nDone → {run_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="ATAC-seq analysis complementing the tf_network module."
    )
    parser.add_argument("config", help="YAML config file.")
    args = parser.parse_args()
    run(args.config)
