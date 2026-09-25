"""
cell_cycle_atac_replication.py

Cell-cycle annotation for scATAC that does not go through gene sets at all.

The gene-set approach asks whether cell-cycle gene promoters are more open in
a cell. They are not: promoter accessibility barely oscillates through the
cycle, which is why that route measures nothing in practice. This module uses
two signals that are native to the assay instead.

  S phase  — while a cell replicates it carries two copies of the genome it
             has already copied and one copy of the rest. Replication is
             spatially organized into megabase-scale domains, so this makes
             large-bin fragment coverage OVERDISPERSED relative to Poisson.
             G1 and G2M cells have uniform copy number and are not. This needs
             no replication-timing reference: it only assumes replication is
             regional, not which regions go first.

  G2M      — a cell that has finished S carries twice the DNA, so at matched
             capture efficiency it yields about twice the unique fragments.

Together:
             dispersion   DNA content
  G1/G0      low          1x
  S          high         between 1x and 2x
  G2M        low          2x

The dispersion statistic is depth-corrected. Writing x_ij for the fragments
cell i puts in bin j, n_i for its total and p_j for the population's share in
that bin,

  chi2_i = sum_j (x_ij - n_i p_j)^2 / (n_i p_j)      (= B under pure Poisson)
  D_i    = chi2_i / B                                 (index of dispersion)
  phi_i  = (D_i - 1) * B / n_i

phi is the variance of relative copy number across bins. It is ~0 for a cell
with uniform copy number and ~0.11 for a cell exactly halfway through S
(half the genome at 2 copies, half at 1), and unlike D it does not grow with
sequencing depth — so it is comparable between shallow and deep cells.

Because postmitotic neurons cannot cycle, a per-cell-type breakdown is a
built-in positive/negative control, and this module always prints one.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

PHASE_ORDER_REPLICATION = ["G1/G0", "S", "G2M"]

# phi expected at mid-S: half the genome at 2 copies, half at 1, so relative
# copy number is {4/3, 2/3} and its variance is 1/9.
PHI_MID_S = 1.0 / 9.0

DEFAULT_BIN_SIZE       = 1_000_000   # replication domains are ~0.4-0.8 Mb
DEFAULT_MIN_BIN_PCT    = 5.0         # drop the lowest-coverage bins (gaps, centromeres)
DEFAULT_S_THRESHOLD    = 0.02        # phi above the cell-type baseline
DEFAULT_G2M_THRESHOLD  = 0.6         # log2 DNA content above the cell-type baseline
DEFAULT_S_Z_THRESHOLD  = 3.0         # sd of phi's Poisson null a cell must clear
MIN_G2M_SEPARATION     = 3.0         # sd between the 1x and 2x DNA modes to call G2M

_AUTOSOME_RE = re.compile(r"^(?:chr)?(\d+)$")


# ---------------------------------------------------------------------------
# Genome bins
# ---------------------------------------------------------------------------

def select_autosomes(chrom_sizes) -> list[str]:
    """Autosomes only, in numeric order.

    X and Y are excluded on purpose: their copy number differs between male
    and female cells, which would show up as a constant per-cell deviation
    from the population average and inflate dispersion for everyone.
    """
    out = []
    for c in chrom_sizes:
        m = _AUTOSOME_RE.match(str(c))
        if m:
            out.append((int(m.group(1)), str(c)))
    return [c for _, c in sorted(out)]


BIN_SIZE_CHOICES = [500_000, 1_000_000, 2_000_000, 5_000_000, 10_000_000]
TARGET_FRAGS_PER_BIN = 8.0


def choose_bin_size(median_depth: float, genome_bp: float,
                    target_per_bin: float = TARGET_FRAGS_PER_BIN) -> int:
    """Pick a bin width giving roughly `target_per_bin` fragments per bin.

    Bin width trades two errors against each other. Narrow bins resolve
    replication domains but leave so few fragments per bin that Poisson noise
    swamps the copy-number signal; wide bins are quiet but average the domain
    structure away. Empirically the sweet spot is 6-12 fragments per bin, and
    since depth varies by dataset the right width does too.
    """
    if not np.isfinite(median_depth) or median_depth <= 0:
        return DEFAULT_BIN_SIZE
    raw = genome_bp * target_per_bin / median_depth
    return min(BIN_SIZE_CHOICES, key=lambda c: abs(np.log(c / raw)))


def build_genome_bins(chrom_sizes: dict, bin_size: int,
                      chroms: list[str] | None = None) -> list[str]:
    """Fixed-width bins as 'chrom:start-end'. A trailing bin shorter than half
    the bin width is dropped rather than left to look under-covered."""
    chroms = chroms if chroms is not None else select_autosomes(chrom_sizes)
    regions = []
    for c in chroms:
        size = int(chrom_sizes[c])
        for start in range(0, size, bin_size):
            end = min(start + bin_size, size)
            if end - start >= bin_size // 2:
                regions.append(f"{c}:{start}-{end}")
    return regions


# ---------------------------------------------------------------------------
# Dispersion
# ---------------------------------------------------------------------------

MIN_GROUP_CELLS = 30


def coverage_dispersion(X, min_bin_pct: float = DEFAULT_MIN_BIN_PCT, groups=None):
    """Per-cell overdispersion of large-bin coverage.

    Returns (phi, D, n, n_bins_kept, phi_null_sd). `phi` is the depth-corrected
    copy-number variance described in the module docstring; cells with no
    fragments get NaN. Bins below `min_bin_pct` of the population coverage
    distribution are dropped — assembly gaps and centromeres are empty for
    every cell and would otherwise contribute a fixed, meaningless chi-square
    term.

    `groups` (cell-type labels) is important, not cosmetic. The expected bin
    profile p_j has to come from cells of the same type: a rare cell type with
    a distinctive accessibility profile deviates systematically from the
    population average, and that deviation lands in phi as if it were copy
    number. Estimated within type, phi measures only the within-type variation
    that replication actually causes. Groups smaller than MIN_GROUP_CELLS fall
    back to the global profile, since their own would be too noisy.
    """
    import scipy.sparse as sp

    X = X.tocsr() if sp.issparse(X) else np.asarray(X, dtype=float)
    bin_total = np.asarray(X.sum(axis=0)).ravel()
    nonzero = bin_total > 0
    if not nonzero.any():
        raise RuntimeError("Every genome bin is empty — check the fragment import.")
    cutoff = np.percentile(bin_total[nonzero], min_bin_pct)
    keep = bin_total > max(cutoff, 0.0)
    if keep.sum() < 100:
        raise RuntimeError(
            f"Only {int(keep.sum())} bins survived coverage filtering; "
            f"dispersion needs many bins. Try a larger --bin-size.")

    Xk = X[:, keep]
    n = np.asarray(Xk.sum(axis=1)).ravel().astype(float)
    B = int(keep.sum())
    safe_n = np.where(n > 0, n, 1.0)

    n_cells = Xk.shape[0]
    if groups is None:
        blocks = [("<all>", np.arange(n_cells))]
    else:
        g = pd.Series(groups).astype(str).to_numpy()
        blocks, small = [], []
        for name in pd.unique(g):
            idx = np.flatnonzero(g == name)
            (blocks if len(idx) >= MIN_GROUP_CELLS else small).append((name, idx))
        if small:
            pooled = np.concatenate([idx for _, idx in small])
            print(f"    [note] {len(small)} group(s) under {MIN_GROUP_CELLS} cells "
                  f"({len(pooled):,} cells) pooled for profile estimation")
            blocks.append(("<small groups>", pooled))

    D = np.full(n_cells, np.nan)
    phi = np.full(n_cells, np.nan)
    phi_null_sd = np.full(n_cells, np.nan)

    for _, idx in blocks:
        Xg = Xk[idx]
        tot = np.asarray(Xg.sum(axis=0)).ravel()
        grand = tot.sum()
        if grand <= 0:
            continue
        p = tot / grand
        # Bins empty across the whole group contribute nothing: x_ij is 0 there
        # for every cell in it, so the chi-square term is 0/0. Zeroing 1/p makes
        # that explicit, and those bins drop out of the degrees of freedom.
        inv_p = np.where(p > 0, 1.0 / np.where(p > 0, p, 1.0), 0.0)
        B_eff = int((p > 0).sum())
        if B_eff < 2:
            continue

        # chi2_i = (1/n_i) * sum_j x_ij^2 / p_j - n_i, using sum_j p_j = 1 and
        # sum_j x_ij = n_i. One sparse matvec against the squared matrix.
        sq = Xg.multiply(Xg) @ inv_p if sp.issparse(Xg) else (Xg ** 2) @ inv_p
        sq = np.asarray(sq).ravel()

        n_g = safe_n[idx]
        chi2_g = sq / n_g - n[idx]
        D[idx] = chi2_g / B_eff
        phi[idx] = (D[idx] - 1.0) * B_eff / n_g
        phi_null_sd[idx] = np.sqrt(2.0 * B_eff) / n_g

    # phi_null_sd above comes from chi2 ~ chi-square with B_eff dof under the
    # Poisson null, so sd(D) = sqrt(2/B) and sd(phi) = sqrt(2B)/n. It falls
    # with depth, which is why a single absolute cutoff on phi over-calls S in
    # shallow cells and under-calls it in deep ones — the classifier forms a
    # per-cell z from it instead.
    D[n <= 0] = np.nan
    phi[n <= 0] = np.nan
    phi_null_sd[n <= 0] = np.nan
    return phi, D, n, B, phi_null_sd


# ---------------------------------------------------------------------------
# Baselines
# ---------------------------------------------------------------------------

def _grouped_median(values: np.ndarray, groups: pd.Series | None) -> np.ndarray:
    """Median of `values` within each group, broadcast back to every cell.
    Falls back to a single global median when no grouping is given."""
    if groups is None:
        return np.full(len(values), float(np.nanmedian(values)))
    df = pd.DataFrame({"v": values, "g": pd.Series(groups).astype(str).values})
    med = df.groupby("g", observed=True)["v"].transform("median")
    filled = med.to_numpy(dtype=float)
    # a group whose median is undefined falls back to the global one
    bad = ~np.isfinite(filled)
    if bad.any():
        filled[bad] = float(np.nanmedian(values))
    return filled


def baseline_corrected(values: np.ndarray, groups) -> np.ndarray:
    """`values` minus their within-group median.

    Both signals need this. Cell types differ in how heterogeneous their
    chromatin is (a phi floor) and in how many fragments they yield (a DNA
    offset), and neither difference is cell cycle. Subtracting the group
    median measures each cell against others of its own kind. It assumes most
    cells in a group are not in S / G2M, which holds in brain tissue.
    """
    return values - _grouped_median(values, groups)


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def classify_phases(s_score: np.ndarray, s_z: np.ndarray, dna_score: np.ndarray,
                    s_threshold: float, s_z_threshold: float,
                    g2m_threshold: float, call_g2m: bool = True) -> np.ndarray:
    """S wins over G2M: a cell still replicating is in S even once it has
    picked up appreciable extra DNA. Everything else is G1/G0 — this method
    cannot separate G0 from G1, since neither has replicated.

    An S call needs both a large enough effect (`s_score` over `s_threshold`,
    so the overdispersion is big enough to be replication) and enough evidence
    for it (`s_z` over `s_z_threshold`, so a shallow cell's counting noise
    cannot produce it). Requiring only the first over-calls S in low-depth
    cells; requiring only the second over-calls it in very deep ones.
    """
    valid = np.isfinite(s_score)
    is_s = valid & (s_score > s_threshold) & np.isfinite(s_z) & (s_z > s_z_threshold)
    is_g2m = (call_g2m & valid & ~is_s
              & np.isfinite(dna_score) & (dna_score > g2m_threshold))

    phases = np.full(len(s_score), None, dtype=object)
    phases[valid] = "G1/G0"
    phases[is_g2m] = "G2M"
    phases[is_s] = "S"
    return phases


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def _pct(a, q):
    return float(np.nanpercentile(a, q))


def print_signal_tables(phi, s_score, dna_score, n_fragment,
                        s_threshold, g2m_threshold):
    print("\n  Per-cell signal distribution:")
    print(f"    {'signal':16s}  {'threshold':>10s}  {'p10':>10s}  {'median':>10s}  "
          f"{'p90':>10s}  {'p99':>10s}")
    for name, arr, th in (
        ("phi (raw)",       phi,       float("nan")),
        ("s_score",         s_score,   s_threshold),
        ("log2 DNA",        dna_score, g2m_threshold),
        ("n_fragment",      n_fragment, float("nan")),
    ):
        th_s = "-" if not np.isfinite(th) else f"{th:.4g}"
        print(f"    {name:16s}  {th_s:>10s}  {_pct(arr,10):>10.4g}  "
              f"{_pct(arr,50):>10.4g}  {_pct(arr,90):>10.4g}  {_pct(arr,99):>10.4g}")
    print(f"\n    mid-S would give phi ~ {PHI_MID_S:.3f}; "
          f"observed p99 of s_score is {_pct(s_score, 99):.4g}")


def print_dna_resolution(dna_score, s_score, s_threshold):
    """Can DNA content separate G2M at all?

    The 2x of G2M is a 1.0 shift in log2 DNA. If technical capture variation
    among non-replicating cells is already that wide, the DNA axis cannot
    resolve G2M no matter where the threshold sits, and the G2M calls will be
    the tail of the depth distribution rather than a phase.
    """
    quiet = np.isfinite(s_score) & (s_score <= s_threshold) & np.isfinite(dna_score)
    if quiet.sum() < 100:
        print("\n  [warn] too few non-replicating cells to judge DNA-content resolution.")
        return False
    spread = float(np.nanstd(dna_score[quiet]))
    sep = 1.0 / spread if spread > 0 else float("inf")
    print(f"\n  DNA-content resolution: sd(log2 DNA) among non-replicating cells "
          f"= {spread:.2f}  ->  G2M's 2x is {sep:.1f} sd away")
    if sep < MIN_G2M_SEPARATION:
        print(f"    [warn] below the {MIN_G2M_SEPARATION:.0f} sd needed for usable G2M calls. "
              f"Capture variation, not")
        print(f"           DNA content, dominates the fragment count here, so a G2M call "
              f"would just be")
        print(f"           the deep tail of the depth distribution. G2M calls are DISABLED; "
              f"those cells")
        print(f"           are reported as G1/G0. Pass --force-g2m to override.")
        return False
    return True


def print_consistency_check(phases, dna_score):
    """S cells should sit between G1 and G2M in DNA content.

    This is the strongest internal check available: it is implied by the
    biology and is not something noise in phi would produce, because phi and
    fragment count are computed from different aspects of the data.
    """
    print("\n  Consistency check — median log2 DNA by called phase:")
    med = {}
    for phase in PHASE_ORDER_REPLICATION:
        sel = (phases == phase) & np.isfinite(dna_score)
        med[phase] = float(np.nanmedian(dna_score[sel])) if sel.sum() else float("nan")
        print(f"    {phase:8s}  n={int(sel.sum()):>9,}  median log2 DNA = {med[phase]:+.3f}")
    g1, s = med.get("G1/G0", np.nan), med.get("S", np.nan)
    if np.isfinite(g1) and np.isfinite(s):
        if s > g1:
            print(f"    S sits {s - g1:+.3f} above G1/G0 in DNA content, as replication implies.")
        else:
            print(f"    [warn] S sits {s - g1:+.3f} relative to G1/G0 — replicating cells "
                  f"should carry MORE DNA, not less.")
            print(f"           The dispersion signal is probably not replication.")


def print_cell_type_table(phases, phi, dna_score, cell_types, max_rows: int = 30):
    """Per-cell-type breakdown. Postmitotic neurons are a built-in negative
    control and progenitors a positive one, so this table is the main evidence
    that the signal is cell cycle rather than something technical."""
    if cell_types is None:
        print("\n  [note] no cell-type column given — skipping the per-cell-type control.")
        print("         Pass --cell-type-col to get the strongest available validation.")
        return
    df = pd.DataFrame({
        "cell_type": pd.Series(cell_types).astype(str).values,
        "phase": phases,
        "phi": phi,
        "dna": dna_score,
    })
    rows = []
    for ct, grp in df.groupby("cell_type", observed=True):
        n = len(grp)
        rows.append({
            "cell_type": ct, "n": n,
            "pct_S":   100.0 * (grp["phase"] == "S").sum() / n,
            "pct_G2M": 100.0 * (grp["phase"] == "G2M").sum() / n,
            "median_phi": float(np.nanmedian(grp["phi"])),
        })
    tab = pd.DataFrame(rows).sort_values("pct_S", ascending=False)
    print("\n  Phase composition by cell type (progenitors should cycle, "
          "neurons should not):")
    print(f"    {'cell_type':45s}  {'n':>9s}  {'%S':>7s}  {'%G2M':>7s}  {'med phi':>9s}")
    for _, r in tab.head(max_rows).iterrows():
        print(f"    {r['cell_type'][:45]:45s}  {int(r['n']):>9,}  "
              f"{r['pct_S']:>7.2f}  {r['pct_G2M']:>7.2f}  {r['median_phi']:>9.4f}")
    if len(tab) > max_rows:
        print(f"    ... {len(tab) - max_rows} more cell type(s)")
    return tab


def print_phase_distribution(phases, n_obs, n_missing=0):
    counts = pd.Series(phases).value_counts()
    print("\n  Phase distribution (ATAC, replication-based):")
    for phase in PHASE_ORDER_REPLICATION:
        n = int(counts.get(phase, 0))
        print(f"    {phase:8s}  {n:>10,}  ({100.0 * n / n_obs if n_obs else 0:5.2f}%)")
    if n_missing:
        print(f"    {'(no data)':8s}  {n_missing:>10,}")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def annotate_cell_cycle_atac_replication(
    adata,
    fragments_path: str | Path,
    genome: str = "hg38",
    bin_size: int | None = None,
    min_bin_pct: float = DEFAULT_MIN_BIN_PCT,
    cell_type_col: str | None = None,
    s_threshold: float = DEFAULT_S_THRESHOLD,
    s_z_threshold: float = DEFAULT_S_Z_THRESHOLD,
    g2m_threshold: float = DEFAULT_G2M_THRESHOLD,
    force_g2m: bool = False,
    tempdir: str | Path | None = None,
    n_jobs: int = 8,
):
    """Annotate `adata` with replication-based cell-cycle calls, in place.

    Adds obs columns:
        cc_phi           raw coverage overdispersion
        cc_s_score       phi minus its within-cell-type median
        cc_s_z           cc_s_score in units of phi's Poisson null spread
        cc_log2_dna      log2 fragment count minus its within-cell-type median
        n_fragment       total fragments in the retained autosomal bins
        CellCyclePhase   G1/G0 | S | G2M  (NaN for cells with no fragments)

    Cells absent from the fragments file get NaN throughout.
    """
    import tempfile
    import snapatac2 as snap
    from atac_seq_analysis import get_snapatac2_genome
    from cell_cycle_annotation import _detect_fragments_chrom_style, _to_ensembl_chrom

    bin_size_auto = bin_size is None
    fragments_path = Path(fragments_path)
    if not fragments_path.exists():
        raise FileNotFoundError(f"fragments file not found: {fragments_path}")

    groups = None
    if cell_type_col is not None:
        if cell_type_col not in adata.obs.columns:
            raise KeyError(
                f"cell_type_col '{cell_type_col}' not in adata.obs. "
                f"Available: {list(adata.obs.columns)[:25]}")
        groups_all = adata.obs[cell_type_col].astype(str)
    else:
        groups_all = None
        print("  [note] no --cell-type-col given; baselines will be global rather "
              "than per cell type, which is less robust.")

    genome_obj = get_snapatac2_genome(genome)
    frag_style = _detect_fragments_chrom_style(fragments_path)
    print(f"  Chromosome naming in fragments file: '{frag_style}'")
    to_frag_chrom = _to_ensembl_chrom if frag_style == "ensembl" else (lambda c: c)

    chrom_sizes = genome_obj.chrom_sizes
    chrom_sizes = dict(chrom_sizes() if callable(chrom_sizes) else chrom_sizes)
    chrom_sizes = {to_frag_chrom(k): v for k, v in chrom_sizes.items()}

    autosomes = select_autosomes(chrom_sizes)
    genome_bp = float(sum(int(chrom_sizes[c]) for c in autosomes))

    work_dir = Path(tempdir) if tempdir else Path(tempfile.mkdtemp(prefix="cc_repl_"))
    work_dir.mkdir(parents=True, exist_ok=True)

    barcodes = adata.obs_names.astype(str).tolist()
    print(f"  Importing fragments for {len(barcodes):,} whitelisted cells "
          f"(n_jobs={n_jobs}, tempdir={work_dir}) ...")
    frag = snap.pp.import_data(
        str(fragments_path),
        chrom_sizes=chrom_sizes,
        whitelist=barcodes,
        min_num_fragments=0,
        sorted_by_barcode=False,
        tempdir=str(work_dir),
        n_jobs=n_jobs,
    )
    print(f"    {frag.n_obs:,}/{len(barcodes):,} cells found in fragments file")
    if frag.n_obs == 0:
        raise RuntimeError(
            "No whitelisted barcodes found in the fragments file. Check that "
            "adata.obs_names match the fragment barcodes (column 4).")

    if bin_size is None:
        depth_hint = (float(np.median(frag.obs["n_fragment"].values))
                      if "n_fragment" in frag.obs.columns else float("nan"))
        bin_size = choose_bin_size(depth_hint, genome_bp)
        print(f"  Auto bin size: {bin_size:,} bp "
              f"(median {depth_hint:,.0f} fragments/cell, targeting "
              f"~{TARGET_FRAGS_PER_BIN:.0f} per bin)")

    bins = build_genome_bins(chrom_sizes, bin_size, autosomes)
    print(f"  {len(bins):,} genome bins of {bin_size:,} bp across "
          f"{len(autosomes)} autosomes (X/Y excluded: copy number varies by sex)")

    print("  Counting fragments per genome bin ...")
    bin_mat = snap.pp.make_peak_matrix(frag, use_rep=bins)
    print(f"    {bin_mat.n_obs:,} cells x {bin_mat.n_vars:,} bins")

    reg_obs_pre = pd.Index(bin_mat.obs_names.astype(str))
    pos_in_adata_pre = pd.Index(adata.obs_names.astype(str)).get_indexer(reg_obs_pre)
    if (pos_in_adata_pre < 0).any() and groups_all is not None:
        raise RuntimeError("bin matrix contains cells absent from adata; cannot "
                           "assign cell-type groups.")

    print("  Computing coverage dispersion ...")
    dispersion_groups = (groups_all.to_numpy()[pos_in_adata_pre]
                         if groups_all is not None else None)
    phi, D, n_frag, n_bins_kept, phi_null_sd = coverage_dispersion(
        bin_mat.X, min_bin_pct, groups=dispersion_groups)
    per_bin = float(np.nanmedian(n_frag)) / n_bins_kept
    print(f"    {n_bins_kept:,} bins retained after coverage filtering "
          f"(median {np.nanmedian(n_frag):,.0f} fragments/cell, "
          f"{per_bin:.1f} per bin)")
    if per_bin < 4.0:
        wider = choose_bin_size(float(np.nanmedian(n_frag)), genome_bp)
        print(f"    [warn] under ~4 fragments/bin, Poisson noise hides all but "
              f"deep-S cells, so S")
        print(f"           recall will be low. {wider:,} bp bins would suit this "
              f"depth better.")

    # baselines are computed on the cells that have fragments, in bin_mat order
    reg_obs = reg_obs_pre
    pos_in_adata = pos_in_adata_pre
    keep = pos_in_adata >= 0
    groups = (groups_all.to_numpy()[pos_in_adata[keep]] if groups_all is not None else None)

    log2_dna_raw = np.full(len(n_frag), np.nan)
    good = n_frag > 0
    log2_dna_raw[good] = np.log2(n_frag[good])

    sub_phi, sub_dna = phi[keep], log2_dna_raw[keep]
    s_score_sub = baseline_corrected(sub_phi, groups)
    dna_score_sub = baseline_corrected(sub_dna, groups)
    with np.errstate(divide="ignore", invalid="ignore"):
        s_z_sub = s_score_sub / phi_null_sd[keep]

    def _align(sub):
        arr = np.full(adata.n_obs, np.nan)
        arr[pos_in_adata[keep]] = sub
        return arr

    aligned = {
        "cc_phi":      _align(sub_phi),
        "cc_s_score":  _align(s_score_sub),
        "cc_log2_dna": _align(dna_score_sub),
        "cc_s_z":      _align(s_z_sub),
        "n_fragment":  _align(n_frag[keep]),
    }
    n_missing = int(np.isnan(aligned["cc_s_score"]).sum())
    if n_missing:
        print(f"  [note] {n_missing:,} cells have no fragments -> NaN score / NaN phase")

    print_signal_tables(aligned["cc_phi"], aligned["cc_s_score"],
                        aligned["cc_log2_dna"], aligned["n_fragment"],
                        s_threshold, g2m_threshold)
    median_null = float(np.nanmedian(phi_null_sd[keep]))
    print(f"    phi's Poisson null spread is {median_null:.4g} at the median depth; "
          f"an S call needs")
    print(f"    both s_score > {s_threshold:g} and s_z > {s_z_threshold:g} "
          f"(i.e. > {s_z_threshold * median_null:.4g} at that depth).")

    # the resolution check needs preliminary S calls, so classify twice
    prelim = classify_phases(aligned["cc_s_score"], aligned["cc_s_z"],
                             aligned["cc_log2_dna"], s_threshold,
                             s_z_threshold, g2m_threshold, call_g2m=False)
    resolvable = print_dna_resolution(
        aligned["cc_log2_dna"], aligned["cc_s_score"], s_threshold)
    call_g2m = bool(resolvable or force_g2m)
    if force_g2m and not resolvable:
        print("    [note] --force-g2m given; making G2M calls anyway.")

    phases = classify_phases(
        aligned["cc_s_score"], aligned["cc_s_z"], aligned["cc_log2_dna"],
        s_threshold, s_z_threshold, g2m_threshold, call_g2m=call_g2m)
    print_consistency_check(phases, aligned["cc_log2_dna"])
    print_cell_type_table(phases, aligned["cc_phi"], aligned["cc_log2_dna"],
                          groups_all.to_numpy() if groups_all is not None else None)
    print_phase_distribution(phases, adata.n_obs, n_missing)

    for col, arr in aligned.items():
        adata.obs[col] = arr.astype(np.float32)
    adata.obs["CellCyclePhase"] = pd.Categorical(
        phases, categories=PHASE_ORDER_REPLICATION, ordered=False)
    adata.uns["cell_cycle_scoring"] = {
        "scoring": "replication",
        "bin_size": int(bin_size),
        "bin_size_auto": bool(bin_size_auto),
        "n_bins_kept": int(n_bins_kept),
        "s_threshold": float(s_threshold),
        "s_z_threshold": float(s_z_threshold),
        "g2m_threshold": float(g2m_threshold),
        "g2m_called": bool(call_g2m),
        "cell_type_col": cell_type_col,
    }
    return adata
