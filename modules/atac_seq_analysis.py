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

# NumPy 2.0 removed `np.asfarray`; snapatac2 still uses it in its internal
# BH-correction helper. Restore a compat shim before snapatac2 is imported.
if not hasattr(np, "asfarray"):
    np.asfarray = lambda a, dtype=np.float64: np.asarray(a, dtype=dtype)


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


def _normalize_chrom(name: str) -> str:
    """Ensure UCSC-style ('chr1', 'chrX', 'chrM') naming to match GENCODE FASTA.
    Ensembl uses '1', 'X', 'MT' — we prepend 'chr' and rename 'MT' → 'chrM'."""
    s = str(name)
    if s.startswith("chr"):
        return s
    if s in ("MT", "Mt", "mt"):
        return "chrM"
    return f"chr{s}"


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
                "Chromosome": [_normalize_chrom(c) for c in var[chrom_c].astype(str).values],
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
            "Chromosome": [_normalize_chrom(c) for c in var[chrom_c].astype(str).values],
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
        rows.append({
            "Chromosome": _normalize_chrom(chrom),
            "Start": start, "End": end, "peak_id": name,
        })
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

    If a GTF file is provided, parse it. Otherwise try to obtain a GTF from
    SnapATAC2's built-in `Genome` object (this downloads on first use and
    caches locally).
    """
    if gtf_path is not None and Path(gtf_path).exists():
        print(f"\n[Load] gene TSSes from GTF: {gtf_path}")
        return _tss_from_gtf(Path(gtf_path))

    print(f"\n[Load] gene TSSes from snapatac2 built-in ({genome})")
    try:
        import snapatac2 as snap
    except ImportError as e:
        raise ImportError("pip install snapatac2") from e

    genome_obj = getattr(snap.genome, genome, None)
    if genome_obj is None:
        raise ValueError(f"snapatac2.genome has no attribute '{genome}'.")

    # Try, in order:
    #   1. .fetch_annotations()  — returns a path to a cached GTF
    #   2. .annotation as a str/Path pointing to a GTF
    #   3. .annotation as an iterable of records → DataFrame
    fetched_gtf = None
    if hasattr(genome_obj, "fetch_annotations"):
        try:
            fetched_gtf = Path(genome_obj.fetch_annotations())
        except Exception as e:
            print(f"  [warn] genome.fetch_annotations() failed: {e}")

    if fetched_gtf is None:
        ann = getattr(genome_obj, "annotation", None)
        if isinstance(ann, (str, Path)):
            fetched_gtf = Path(ann)

    if fetched_gtf is not None and fetched_gtf.exists():
        print(f"  [snapatac2 GTF] {fetched_gtf}")
        return _tss_from_gtf(fetched_gtf)

    ann = getattr(genome_obj, "annotation", None)
    if ann is None:
        raise RuntimeError(
            f"Could not obtain gene annotations for '{genome}' from snapatac2. "
            "Please set atac.gtf_path in the config to a GENCODE/Ensembl GTF."
        )
    try:
        df = pd.DataFrame(list(ann))
    except Exception as e:
        raise RuntimeError(
            f"Failed to parse snap.genome.{genome}.annotation "
            f"(type={type(ann).__name__}): {e}. "
            "Please set atac.gtf_path in the config to a GENCODE/Ensembl GTF."
        )
    rename = {"chrom": "Chromosome", "start": "Start", "end": "End",
              "strand": "Strand", "gene_name": "gene_name"}
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    df["TSS"] = np.where(df["Strand"] == "+", df["Start"], df["End"])
    df["Chromosome"] = df["Chromosome"].map(_normalize_chrom)
    return df[["Chromosome", "TSS", "Strand", "gene_name"]]


_GENE_FEATURE_TYPES = ("gene", "transcript", "mRNA")
# GTF attribute keys where a gene symbol may live, in order of preference.
_GENE_NAME_KEYS = ("gene_name", "Name", "gene_symbol", "symbol")
# Fallback if no symbol column exists — Ensembl IDs.
_GENE_ID_KEYS   = ("gene_id", "ID")


def _tss_from_gtf(gtf_path: Path) -> pd.DataFrame:
    """Parse a GTF (plain or gzipped) into Chromosome / TSS / Strand / gene_name.

    Tolerates GTFs where the feature column is "transcript" instead of "gene",
    and GTFs whose attribute uses "Name"/"gene_symbol" instead of "gene_name".
    Falls back to Ensembl gene_id if no symbol is present.
    """
    import gzip
    opener = gzip.open if str(gtf_path).endswith(".gz") else open

    rows = []
    feature_type_counts: dict[str, int] = {}
    sample_lines: list[str] = []
    sample_attrs_keys: set[str] = set()

    def _find_key(attrs: str, keys) -> str | None:
        for k in keys:
            # GTF: key "value";
            m = re.search(rf'\b{re.escape(k)} "([^"]+)"', attrs)
            if m:
                return m.group(1)
            # GFF3: key=value; (values are semicolon-terminated or end-of-line)
            m = re.search(rf'(?:^|[;\s]){re.escape(k)}=([^;]+)', attrs)
            if m:
                return m.group(1).strip()
        return None

    with opener(gtf_path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            feature = parts[2]
            feature_type_counts[feature] = feature_type_counts.get(feature, 0) + 1
            if len(sample_lines) < 3:
                sample_lines.append(line.rstrip("\n"))
                # Extract attr keys once for diagnostic
                for m in re.finditer(r'(\w+) "', parts[8]):
                    sample_attrs_keys.add(m.group(1))

            if feature not in _GENE_FEATURE_TYPES:
                continue
            chrom, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
            attrs = parts[8]

            name = _find_key(attrs, _GENE_NAME_KEYS) or _find_key(attrs, _GENE_ID_KEYS)
            if not name:
                continue
            tss = start if strand == "+" else end
            rows.append({
                "Chromosome": _normalize_chrom(chrom),
                "TSS": tss, "Strand": strand, "gene_name": name,
            })

    if not rows:
        top_features = sorted(feature_type_counts.items(), key=lambda x: -x[1])[:8]
        raise ValueError(
            f"Parsed 0 gene rows from GTF: {gtf_path}\n"
            f"  Feature-type counts (top 8): {top_features}\n"
            f"  Attribute keys seen: {sorted(sample_attrs_keys)}\n"
            f"  Sample lines (up to 3):\n    "
            + "\n    ".join(sample_lines[:3])
            + "\n"
            f"  Recognized feature types: {_GENE_FEATURE_TYPES}\n"
            f"  Recognized name keys:     {_GENE_NAME_KEYS + _GENE_ID_KEYS}\n"
            "Fix: point atac.gtf_path at a GENCODE-format GTF, or extend the "
            "recognized keys above."
        )

    df = pd.DataFrame(rows)
    # Deduplicate — one TSS per gene_name (keep the first, i.e. the primary).
    df = df.drop_duplicates(subset=["gene_name"], keep="first").reset_index(drop=True)
    print(f"  [gtf] parsed {len(df):,} gene TSSes  "
          f"(feature types: {sorted(feature_type_counts, key=lambda k: -feature_type_counts[k])[:5]})")
    return df


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

    # snapatac2's motif_enrichment expects regions as "chr:start-end" strings.
    peak_region_strs = [
        f"{c}:{s}-{e}"
        for c, s, e in zip(peak_df["Chromosome"], peak_df["Start"], peak_df["End"])
    ]

    per_type_top_peaks: dict[str, list[str]] = {}
    for ct in cell_types:
        mask = (adata.obs[cell_type_col].astype(str) == str(ct)).to_numpy()
        n_ct = int(mask.sum())
        if n_ct == 0:
            print(f"  [warn] no cells for cell_type='{ct}'")
            continue
        ct_mean = np.asarray(X[mask, :].mean(axis=0)).ravel()
        n_top   = max(1, int(len(peak_ids) * top_peaks_pct))
        top_ix  = np.argsort(ct_mean)[::-1][:n_top]
        per_type_top_peaks[str(ct)] = [peak_region_strs[i] for i in top_ix]
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
        for ct, region_strs in per_type_top_peaks.items():
            result[ct] = snap.tl.motif_enrichment(
                motifs=motifs, regions=region_strs, genome_fasta=genome_obj,
            )

    # Flatten result into a long DataFrame regardless of exact return shape.
    rows = []
    printed_debug = False
    for ct, res in (result.items() if isinstance(result, dict) else [("all", result)]):
        res_df = _to_pandas(res)

        if not printed_debug and not res_df.empty:
            print(f"  [debug] motif_enrichment columns for '{ct}': {list(res_df.columns)}")
            print(f"  [debug] index name: {res_df.index.name}  "
                  f"first index value: {res_df.index[0]!r}")
            print(f"  [debug] first row: {res_df.iloc[0].to_dict()}")
            printed_debug = True

        for idx, r in res_df.iterrows():
            motif_name = _extract_motif_name(idx, r)
            tf_name    = _extract_str(r, ["family", "TF", "tf", "gene_name", "gene_symbol", "name"])
            log2fe     = _extract_numeric(r, [
                "log2_fold_enrichment", "log2FE", "log2 fold enrichment",
                "log2(fold change)", "log2 fold change", "log2FoldChange",
                "log2FC", "log2 fold-enrichment", "fold enrichment",
            ])
            padj       = _extract_numeric(r, [
                "adjusted_p_value", "adjusted p-value", "adjusted p value",
                "padj", "FDR", "q-value", "qvalue", "adj_pvalue",
            ])
            pval       = _extract_numeric(r, [
                "p-value", "p_value", "pvalue", "P-value", "P value",
            ])
            rows.append({
                "cell_type": ct, "motif": motif_name,
                "tf_name": tf_name,
                "log2_fold_enrichment": log2fe,
                "p_value": pval,
                "adjusted_p_value": padj,
            })
    enrich_df = pd.DataFrame(rows)

    # Restrict rows to candidate TFs. Try the dedicated tf_name column first
    # (more reliable), fall back to substring-matching the motif ID.
    if not enrich_df.empty:
        candidate_upper = {tf.upper() for tf in candidate_tfs}
        def _match(row) -> str | None:
            tf = str(row.get("tf_name", "") or "").upper()
            if tf in candidate_upper:
                # return the original-case TF from the candidate set
                for c in candidate_tfs:
                    if c.upper() == tf:
                        return c
            motif = str(row.get("motif", "") or "").upper()
            for c in candidate_tfs:
                if c.upper() in motif:
                    return c
            return None
        enrich_df["candidate_tf"] = enrich_df.apply(_match, axis=1)
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
        pivot = _order_heatmap_by_enrichment(pivot)
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

    # ── All-regions heatmap (mirrors the layout of Task 1) ──────────────────
    if task_cfg.get("all_regions_heatmap", True):
        _accessibility_heatmap_all_regions(
            adata, overlap_df, peak_to_ix, X,
            region_col=region_col,
            regions=task_cfg.get("heatmap_regions"),
            exclude_regions=task_cfg.get(
                "heatmap_exclude_regions", ["Brain", "Head"]
            ),
            top_n_genes=int(task_cfg.get("heatmap_top_n_genes", 50)),
            gene_ranking=task_cfg.get("heatmap_gene_ranking", "variance"),
            out_dir=out_dir,
        )

    if not all_rows:
        return None
    combined = pd.concat(all_rows, ignore_index=True)
    combined_csv = out_dir / "promoter_accessibility_all_pairs.csv"
    combined.to_csv(combined_csv, index=False)
    print(f"  Saved: {combined_csv.name}  ({len(combined):,} rows total)")
    return combined_csv


def _accessibility_heatmap_all_regions(
    adata, overlap_df, peak_to_ix, X,
    region_col: str,
    regions: list | None,
    exclude_regions: list,
    top_n_genes: int,
    gene_ranking: str,
    out_dir: Path,
) -> None:
    """Compute per-region mean promoter accessibility for every target gene and
    save a full-matrix CSV + a top-N-genes heatmap (region × gene). The heatmap
    reuses `_order_heatmap_by_enrichment` so rows/columns are ordered to make
    region-specific patterns visually obvious (same convention as Task 1)."""
    obs_region = adata.obs[region_col].astype(str)
    all_regions = list(pd.unique(obs_region))

    exclude_set = set(map(str, exclude_regions or ()))
    if regions:
        target_regions = [r for r in regions if r in all_regions and r not in exclude_set]
    else:
        target_regions = [r for r in all_regions if r not in exclude_set]

    if not target_regions:
        print(f"  [skip] all-regions heatmap: no regions to plot (available: {all_regions})")
        return

    print(f"\n  All-regions heatmap: {len(target_regions)} region(s) × "
          f"{overlap_df['gene_name'].nunique():,} genes")

    # Per-gene per-cell accessibility score (sum of overlapping-peak values).
    per_gene = overlap_df.groupby("gene_name")["peak_id"].apply(list).to_dict()

    # region × gene DataFrame
    region_masks = {r: (obs_region == r).to_numpy() for r in target_regions}
    for r, mask in region_masks.items():
        if int(mask.sum()) == 0:
            print(f"    [warn] region '{r}' has 0 cells — will be all-NaN row.")

    data = {r: {} for r in target_regions}
    for gene, peak_ids in per_gene.items():
        ix = [peak_to_ix[p] for p in peak_ids if p in peak_to_ix]
        if not ix:
            continue
        per_cell_sum = np.asarray(X[:, ix].sum(axis=1)).ravel()
        for r in target_regions:
            mask = region_masks[r]
            data[r][gene] = (
                float(per_cell_sum[mask].mean()) if mask.any() else np.nan
            )

    mat = pd.DataFrame(data).T  # rows = regions, cols = genes
    mat.index.name = region_col
    mat.columns.name = "gene"

    full_csv = out_dir / "promoter_accessibility_region_x_gene.csv"
    mat.to_csv(full_csv)
    print(f"  Saved: {full_csv.name}  ({mat.shape[0]} regions × {mat.shape[1]:,} genes)")

    # Rank genes so the heatmap shows the most informative ones.
    if gene_ranking == "variance":
        scores = mat.var(axis=0, skipna=True)
    elif gene_ranking == "range":
        scores = mat.max(axis=0, skipna=True) - mat.min(axis=0, skipna=True)
    elif gene_ranking == "max":
        scores = mat.max(axis=0, skipna=True)
    else:
        raise ValueError(f"gene_ranking must be one of variance | range | max; got {gene_ranking}")

    top_genes = scores.dropna().sort_values(ascending=False).head(top_n_genes).index
    if len(top_genes) == 0:
        print("  [skip] heatmap: no genes with non-NaN scores")
        return
    sub = mat[top_genes]

    sub = _order_heatmap_by_enrichment(sub)
    _heatmap(
        sub, out_dir / "promoter_accessibility_region_x_gene_heatmap.png",
        title=(f"Promoter accessibility — {region_col} × top-{len(top_genes)} genes "
               f"(ranked by {gene_ranking})"),
        cmap="viridis",
    )
    print(f"  Saved: promoter_accessibility_region_x_gene_heatmap.png")


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
        promoter_region_strs = [
            f"{c}:{s}-{e}"
            for c, s, e in zip(promoter_df["Chromosome"],
                               promoter_df["Start"],
                               promoter_df["End"])
        ]
        try:
            per_promoter_hits = snap.tl.motif_enrichment(
                motifs=tf_motifs,
                regions={"promoters": promoter_region_strs},
                genome_fasta=genome_obj,
            )
        except Exception as e:
            print(f"  [warn] motif scan failed for {tf}: {e}")
            continue

        raw_hits = _flatten_enrichment_result(per_promoter_hits)

        # Normalize column names (snapatac2 uses "p-value" / "adjusted p-value";
        # motif name is often the DataFrame index).
        hits_rows = []
        for idx, r in raw_hits.iterrows():
            hits_rows.append({
                "motif": _extract_motif_name(idx, r),
                "log2_fold_enrichment": _extract_numeric(r, [
                    "log2_fold_enrichment", "log2FE", "log2 fold enrichment",
                    "log2(fold change)", "log2 fold change", "log2FoldChange", "log2FC",
                ]),
                "p_value": _extract_numeric(r, [
                    "p-value", "p_value", "pvalue", "P-value", "P value",
                ]),
                "adjusted_p_value": _extract_numeric(r, [
                    "adjusted_p_value", "adjusted p-value", "adjusted p value",
                    "padj", "FDR", "q-value", "qvalue", "adj_pvalue",
                ]),
            })
        hits_df = pd.DataFrame(hits_rows)

        if hits_df.empty:
            tf_summary_rows.append({
                "TF": tf, "n_targets": len(target_genes),
                "n_motifs_found": len(tf_motifs),
                "n_significant_motif_enrichments": 0,
                "min_adjusted_p_value": np.nan,
            })
            continue

        significant = hits_df[hits_df["adjusted_p_value"] < pval_cutoff]
        tf_summary_rows.append({
            "TF": tf,
            "n_targets": len(target_genes),
            "n_motifs_found": len(tf_motifs),
            "n_significant_motif_enrichments": int(significant.shape[0]),
            "min_adjusted_p_value": float(hits_df["adjusted_p_value"].min(skipna=True)),
        })
        for _, r in hits_df.iterrows():
            rows.append({
                "TF": tf,
                "motif": r["motif"],
                "log2_fold_enrichment": r["log2_fold_enrichment"],
                "p_value": r["p_value"],
                "adjusted_p_value": r["adjusted_p_value"],
            })

    detail_df  = pd.DataFrame(rows)
    summary_df = pd.DataFrame(tf_summary_rows)
    detail_csv = out_dir / "motif_target_validation_details.csv"
    summary_csv = out_dir / "motif_target_validation_summary.csv"
    detail_df.to_csv(detail_csv, index=False)
    summary_df.to_csv(summary_csv, index=False)
    print(f"\n  Saved: {detail_csv.name}   ({len(detail_df):,} rows)")
    print(f"  Saved: {summary_csv.name}  ({len(summary_df):,} TFs)")

    # ── Per-pair MOODS scan (each (TF, target) gets its own row) ───────────
    if task_cfg.get("per_pair_scoring", True):
        _per_pair_motif_scan(
            adata=adata,
            tf_targets=tf_targets,
            tss_df=tss_df,
            motifs=motifs,
            find_motifs_for_tf=_find_motifs_for_tf,
            genome_obj=genome_obj,
            cfg=cfg,
            out_dir=out_dir,
        )

    return summary_csv


# ---------------------------------------------------------------------------
# Per-pair MOODS scan (Task 3, additional output)
# ---------------------------------------------------------------------------

def _extract_pwm(motif) -> np.ndarray | None:
    """Extract a 4×L PWM (probability or count matrix) from a snapatac2 motif.
    Returns None if the format is unrecognizable."""
    for attr in ("counts", "probability_matrix", "matrix", "pwm"):
        if hasattr(motif, attr):
            mat = np.asarray(getattr(motif, attr), dtype=float)
            if mat.ndim == 2:
                # Some libraries store L×4, some 4×L. Normalize to 4×L.
                if mat.shape[0] != 4 and mat.shape[1] == 4:
                    mat = mat.T
                if mat.shape[0] == 4:
                    return mat
    # dict-like access
    if hasattr(motif, "__getitem__"):
        for k in ("matrix", "pwm", "counts"):
            try:
                mat = np.asarray(motif[k], dtype=float)
                if mat.ndim == 2:
                    if mat.shape[0] != 4 and mat.shape[1] == 4:
                        mat = mat.T
                    if mat.shape[0] == 4:
                        return mat
            except Exception:
                continue
    return None


def _pwm_to_moods_log_odds(pwm: np.ndarray, bg=(0.25, 0.25, 0.25, 0.25),
                           pseudocount: float = 0.01) -> list[list[float]]:
    """Convert a 4×L count/probability matrix to a log-odds matrix (rows ACGT)
    in the list-of-lists format MOODS.scan expects."""
    m = np.asarray(pwm, dtype=float)
    # If it's counts, normalize per-column to probabilities.
    col_sums = m.sum(axis=0, keepdims=True)
    col_sums = np.where(col_sums > 0, col_sums, 1.0)
    ppm = m / col_sums
    ppm = np.clip(ppm + pseudocount, 1e-9, None)
    bg_arr = np.asarray(bg).reshape(-1, 1)
    lo = np.log2(ppm / bg_arr)
    return lo.tolist()


def _resolve_fasta_path(genome_obj, cfg) -> Path | None:
    """Find the reference FASTA path snapatac2 uses. Config override wins."""
    explicit = cfg.get("atac", {}).get("fasta_path")
    if explicit and Path(explicit).exists():
        return Path(explicit)
    for attr in ("fasta_file", "fasta", "genome_fasta"):
        v = getattr(genome_obj, attr, None)
        if callable(v):
            try:
                v = v()
            except Exception:
                v = None
        if v and Path(str(v)).exists():
            return Path(str(v))
    # snapatac2 caches under ~/.cache/snapatac2/. Try that.
    cache = Path.home() / ".cache" / "snapatac2"
    for pattern in ("*.fa.gz.decomp", "*.fa", "*.fasta"):
        for p in cache.glob(pattern):
            return p
    return None


def _peaks_by_chrom(peaks_df: pd.DataFrame) -> dict[str, np.ndarray]:
    """Return {chrom: sorted 2-col array of (start, end)} for O(log n) overlap
    queries via numpy searchsorted."""
    out = {}
    for chrom, sub in peaks_df.groupby("Chromosome"):
        starts_ends = sub[["Start", "End"]].to_numpy()
        starts_ends = starts_ends[starts_ends[:, 0].argsort()]
        out[chrom] = starts_ends
    return out


def _hit_in_any_peak(chrom: str, hit_start: int, hit_end: int,
                     peaks_by_chrom: dict) -> bool:
    """True if [hit_start, hit_end) overlaps any peak on `chrom`."""
    arr = peaks_by_chrom.get(chrom)
    if arr is None or len(arr) == 0:
        return False
    # Peaks whose start < hit_end AND end > hit_start.
    idx = np.searchsorted(arr[:, 0], hit_end, side="right")
    for i in range(idx - 1, -1, -1):
        p_start, p_end = arr[i, 0], arr[i, 1]
        if p_end <= hit_start:
            # Since peaks are sorted by start, we can break once we've gone
            # past overlap range on both sides. Cap the backward scan for
            # very large peak sets.
            if arr[max(0, i - 500), 1] <= hit_start:
                break
            continue
        return True
    return False


def _per_pair_motif_scan(
    adata, tf_targets, tss_df, motifs, find_motifs_for_tf,
    genome_obj, cfg, out_dir,
) -> Path | None:
    """MOODS per-sequence scan: one row per (TF, target).

    Reports the best MOODS score across the TF's cis_bp motif variants at each
    target promoter, plus a flag for whether the best hit falls inside an
    accessible ATAC peak.
    """
    task_cfg = cfg.get("motif_target_validation", {})
    pval     = float(task_cfg.get("moods_pvalue_cutoff", 1e-4))
    tss_win  = int(task_cfg.get("promoter_window", 2000))

    print("\n" + "-" * 62)
    print("  Task 3b: per-pair MOODS scan")
    print("-" * 62)

    try:
        import MOODS.scan
        import MOODS.tools
    except ImportError as e:
        raise ImportError(
            "MOODS-python required for per-pair scoring: pip install MOODS-python\n"
            "(or disable via motif_target_validation.per_pair_scoring: false)"
        ) from e

    try:
        import pyfaidx
    except ImportError as e:
        raise ImportError("pyfaidx required: pip install pyfaidx") from e

    fasta_path = _resolve_fasta_path(genome_obj, cfg)
    if fasta_path is None:
        raise RuntimeError(
            "Could not locate reference FASTA. Set atac.fasta_path in the config "
            "to the same file snapatac2 uses (usually under ~/.cache/snapatac2/)."
        )
    print(f"  FASTA: {fasta_path}")
    fasta = pyfaidx.Fasta(str(fasta_path))

    # Peak coords for accessibility overlay (same auto-detection Task 1/2 use).
    peaks_df = peaks_to_dataframe(adata, cfg["atac"].get("peak_columns"))
    peaks_by_chrom = _peaks_by_chrom(peaks_df)
    print(f"  Peak index built: {sum(len(v) for v in peaks_by_chrom.values()):,} peaks "
          f"across {len(peaks_by_chrom)} chromosomes")

    bg = (0.25, 0.25, 0.25, 0.25)
    rows = []
    n_skipped_pwm = 0
    n_skipped_seq = 0

    for tf, tf_group in tf_targets.groupby("TF"):
        tf_motifs = find_motifs_for_tf(str(tf))
        if not tf_motifs:
            continue

        # Build PWMs + thresholds for this TF's motif variants.
        pwms:   list[list[list[float]]] = []
        motif_names: list[str]          = []
        for m in tf_motifs:
            pwm = _extract_pwm(m)
            if pwm is None:
                n_skipped_pwm += 1
                continue
            lo = _pwm_to_moods_log_odds(pwm, bg=bg)
            pwms.append(lo)
            motif_names.append(str(getattr(m, "name", f"motif_{len(motif_names)}")))
        if not pwms:
            continue

        thresholds = [MOODS.tools.threshold_from_p(m, list(bg), pval) for m in pwms]

        target_genes = tf_group["target"].astype(str).unique().tolist()
        promoter_df  = build_promoter_regions(tss_df, tss_win, genes=target_genes)
        if promoter_df.empty:
            continue

        for _, prow in promoter_df.iterrows():
            chrom = str(prow["Chromosome"])
            start = int(prow["Start"])
            end   = int(prow["End"])
            gene  = str(prow["gene_name"])

            try:
                seq = str(fasta[chrom][start:end].seq).upper()
            except (KeyError, ValueError):
                n_skipped_seq += 1
                continue
            if len(seq) < max(len(pwm[0]) for pwm in pwms):
                n_skipped_seq += 1
                continue

            hit_lists = MOODS.scan.scan_dna(seq, pwms, list(bg), thresholds)

            best_score = -np.inf
            best_motif_name = ""
            best_in_peak = False
            n_hits = 0
            n_hits_in_peak = 0

            for motif_i, hits in enumerate(hit_lists):
                width = len(pwms[motif_i][0])
                for pos, score in hits:
                    n_hits += 1
                    hit_start = start + int(pos)
                    hit_end   = hit_start + width
                    in_peak   = _hit_in_any_peak(chrom, hit_start, hit_end, peaks_by_chrom)
                    if in_peak:
                        n_hits_in_peak += 1
                    if float(score) > best_score:
                        best_score      = float(score)
                        best_motif_name = motif_names[motif_i]
                        best_in_peak    = in_peak

            rows.append({
                "TF": tf,
                "target": gene,
                "n_motif_hits": n_hits,
                "n_motif_hits_in_accessible_peak": n_hits_in_peak,
                "best_motif_score": (best_score if np.isfinite(best_score) else np.nan),
                "best_motif": best_motif_name,
                "best_hit_in_accessible_peak": bool(best_in_peak) if n_hits > 0 else False,
                "validated_by_motif_and_accessibility": bool(n_hits_in_peak > 0),
            })

    if n_skipped_pwm:
        print(f"  [warn] {n_skipped_pwm} motif(s) had unrecognizable PWM format and were skipped.")
    if n_skipped_seq:
        print(f"  [warn] {n_skipped_seq} promoter(s) had missing/short sequence and were skipped.")

    if not rows:
        print("  [skip] no per-pair scores produced.")
        return None

    pair_df  = pd.DataFrame(rows).sort_values(["TF", "best_motif_score"],
                                              ascending=[True, False])
    pair_csv = out_dir / "motif_target_pair_scores.csv"
    pair_df.to_csv(pair_csv, index=False)
    print(f"  Saved: {pair_csv.name}  "
          f"({len(pair_df):,} pairs; validated: "
          f"{int(pair_df['validated_by_motif_and_accessibility'].sum()):,})")
    return pair_csv


def _extract_motif_name(idx, row) -> str:
    """Motif name may live in the row's index (typical for snapatac2) or in
    a name/motif/id column. Try each in turn."""
    if idx is not None and not isinstance(idx, (int, np.integer)):
        return str(idx)
    for col in ("name", "motif", "motif_name", "id", "TF"):
        if col in row.index and pd.notna(row[col]):
            return str(row[col])
    return ""


def _extract_str(row, keys) -> str:
    """First matching key coerced to str. Empty string if none."""
    for k in keys:
        if k in row.index:
            v = row[k]
            if pd.notna(v):
                return str(v)
    return ""


def _extract_numeric(row, keys) -> float:
    """First matching key in `keys` present in `row`, coerced to float. NaN if none."""
    for k in keys:
        if k in row.index:
            v = row[k]
            if pd.notna(v):
                try:
                    return float(v)
                except (TypeError, ValueError):
                    pass
    return float("nan")


def _to_pandas(x) -> pd.DataFrame:
    """Convert whatever snapatac2 gave us into a pandas DataFrame with
    preserved column names. Handles polars, pandas, or generic mappings."""
    if isinstance(x, pd.DataFrame):
        return x
    if hasattr(x, "to_pandas"):        # polars.DataFrame
        return x.to_pandas()
    if hasattr(x, "to_dataframe"):     # some snapatac2 result objects
        return x.to_dataframe()
    return pd.DataFrame(x)


def _flatten_enrichment_result(result) -> pd.DataFrame:
    if isinstance(result, dict):
        frames = [_to_pandas(v) for v in result.values()]
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return _to_pandas(result)


# =============================================================================
# Small utilities
# =============================================================================

def _sanitize(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_") or "unnamed"


def _order_heatmap_by_enrichment(df: pd.DataFrame) -> pd.DataFrame:
    """Reorder rows and columns to make enrichment patterns visually obvious.

    - Rows sorted by max enrichment across columns (descending) so the most-
      active cell type is on top.
    - Columns grouped by which row they peak in (in the sorted row order), then
      sorted by descending value within each block. Produces a diagonal /
      block-diagonal structure: TFs preferred by row 0 come first (strongest
      first), then row 1, etc.
    """
    if df.empty:
        return df

    # Row order: strongest cell type first.
    row_scores = df.max(axis=1, skipna=True)
    row_order  = row_scores.sort_values(ascending=False, na_position="last").index
    df = df.reindex(index=row_order)

    # Column order: which row does each column peak in?
    #   idxmax gives the row label (e.g. 'glioblast').
    #   Convert to position via row_order.get_loc for sorting.
    row_pos    = {r: i for i, r in enumerate(row_order)}
    peak_row   = df.idxmax(axis=0, skipna=True).map(row_pos).fillna(len(row_order))
    peak_value = df.max(axis=0, skipna=True)
    col_order  = (
        pd.DataFrame({"peak_row": peak_row, "peak_value": peak_value})
        .sort_values(by=["peak_row", "peak_value"], ascending=[True, False])
        .index
    )
    return df.reindex(columns=col_order)


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
