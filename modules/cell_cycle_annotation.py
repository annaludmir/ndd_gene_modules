"""
cell_cycle_annotation.py

Annotate an h5ad with per-cell cell-cycle information. Adds two obs columns:

  - cycling_score  (float in [0, 1]) — fraction of UMIs contributed by
    cc_genes_human across the cell's total UMI count.
  - CellCyclePhase (categorical str)  — one of:
        Non-cycling  if cycling_score < 0.4%
        G1 / S / G2M if cycling_score >= 0.4% and the phase-specific UMI
                     fraction passes its threshold (highest ratio to
                     threshold wins when multiple phases pass)
        Post-M       if cycling_score >= 0.4% but no phase passes

Thresholds (from the spec):
    cycling_score        > 0.4%   (0.004)
    G1 gene UMI fraction > 0.2% (0.0002)
    S  gene UMI fraction > 0.2% (0.002)
    G2M gene UMI fraction > 3% (0.03)

Usage:
  python modules/cell_cycle_annotation.py \
    --h5ad-input  data/input.h5ad \
    --h5ad-output data/input_with_cc.h5ad

ATAC-seq has two modes (--modality atac):
  * peak-based (default): sum accessibility of called peaks overlapping each
    gene set's promoters (TSS ± --tss-window) from adata.X.
  * fragment-based (--fragments fragments.tsv.bgz): count Tn5 fragments
    directly in each gene's promoter window with snapatac2, independent of
    peak calling. Denominator is the cell's total fragment count.

Scoring (--scoring, applies to every modality):
  * fraction (default): the scheme described above — share of a cell's total
    signal falling in the gene set. Simple, but the resulting value scales with
    how many features the set has and with the set's baseline signal level, so
    absolute thresholds do not transfer between gene sets or between assays.
  * background: score each set against control features matched on mean signal
    (the scanpy `score_genes` / AUCell idea), i.e.

        score = mean(set features) - mean(matched control features)

    on depth-normalized values, reported as a z-score across cells. Using a
    mean rather than a sum removes the set-size dependence; subtracting the
    matched background removes the baseline level. Thresholds are then in
    standard deviations and are comparable across sets.

    This mode also reports, per set, a `signal` ratio: the cell-to-cell spread
    of the real score divided by that of a null score built from a random set
    of the same size and signal level. signal <= 1 means the set carries no
    more structure than a random set, so its phase calls are not interpretable
    no matter where the threshold is put.
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Gene sets
# ---------------------------------------------------------------------------

CC_GENES_HUMAN = np.array([
    'ABHD3', 'AC016205.1', 'AC073529.1', 'AC084033.3', 'AC087632.1',
    'AC091057.6', 'AC097534.2', 'AC099850.2', 'AC135586.2', 'ACAA2',
    'ACADM', 'ACP1', 'ACTL6A', 'ACYP1', 'ADCY3', 'ADD3', 'ADK', 'AHCY',
    'AKIRIN2', 'AKR7A2', 'AL359513.1', 'AL449266.1', 'AL513165.2',
    'ANAPC11', 'ANLN', 'ANP32A', 'ANP32B', 'ANP32E', 'AP001347.1',
    'APOLD1', 'ARHGAP11A', 'ARHGEF39', 'ARID1A', 'ARL6IP1', 'ARL6IP6',
    'ARMC1', 'ASF1B', 'ASPM', 'ASRGL1', 'ATAD2', 'ATAD5',
    'ATP1B3', 'AURKA', 'AURKB', 'BANF1', 'BARD1', 'BAZ1A', 'BAZ1B',
    'BIRC5', 'BLM', 'BORA', 'BRCA1', 'BRCA2', 'BRD8', 'BRIP1', 'BTG3',
    'BUB1', 'BUB1B', 'BUB3', 'C11orf58', 'C19orf48', 'C1orf112',
    'C1orf35', 'C21orf58', 'C5orf34', 'CACYBP', 'CAMTA1',
    'CARHSP1', 'CBX1', 'CBX3', 'CBX5', 'CCAR1', 'CCDC14', 'CCDC167',
    'CCDC18', 'CCDC34', 'CCDC77', 'CCNA1', 'CCNA2', 'CCNB1', 'CCNB2',
    'CCNE2', 'CCNF', 'CCT4', 'CCT5', 'CDC20', 'CDC25B', 'CDC25C',
    'CDC27', 'CDC45', 'CDC6', 'CDC7', 'CDCA2', 'CDCA3', 'CDCA4',
    'CDCA7L', 'CDCA8', 'CDK1', 'CDK19', 'CDK2', 'CDK4', 'CDK5RAP2',
    'CDKAL1', 'CDKN1B', 'CDKN2C', 'CDKN3', 'CDT1', 'CENPA', 'CENPC',
    'CENPE', 'CENPF', 'CENPH', 'CENPI', 'CENPJ', 'CENPK', 'CENPL',
    'CENPM', 'CENPN', 'CENPO', 'CENPP', 'CENPQ', 'CENPU', 'CENPW',
    'CENPX', 'CEP112', 'CEP128', 'CEP135', 'CEP192', 'CEP295', 'CEP55',
    'CEP57', 'CEP57L1', 'CEP70', 'CETN3', 'CFAP20', 'CFL2', 'CGGBP1',
    'CHAF1A', 'CHCHD2', 'CHEK1', 'CHEK2', 'CHRAC1', 'CIP2A', 'CIT',
    'CKAP2', 'CKAP2L', 'CKAP5', 'CKLF', 'CKS1B', 'CKS2', 'CLSPN',
    'CMC2', 'CMSS1', 'CNIH4', 'CNN3', 'CNTLN', 'CNTRL', 'COA1',
    'COMMD4', 'COX8A', 'CSE1L', 'CTCF', 'CTDSPL2', 'CWF19L2', 'CYB5B',
    'CYCS', 'DACH1', 'DBF4', 'DBF4B', 'DBI', 'DCAF7', 'DCP2', 'DCXR',
    'DDAH2', 'DDX39A', 'DDX46', 'DEK', 'DEPDC1', 'DEPDC1B', 'DESI2',
    'DHFR', 'DIAPH3', 'DKC1', 'DLEU2', 'DLGAP5', 'DNA2', 'DNAJB1',
    'DNAJC9', 'DNMT1', 'DPM1', 'DR1', 'DSCC1', 'DSN1', 'DTL', 'DTYMK',
    'DUSP16', 'DUT', 'DYNLL1', 'DYRK1A', 'E2F3', 'E2F7', 'E2F8',
    'ECT2', 'EED', 'EEF1D', 'EID1', 'EIF1AX', 'EIF2S2', 'EIF4A3',
    'EIF4E', 'EIF5', 'EMC9', 'ENAH', 'ENO1', 'ENY2', 'ERH', 'ESCO2',
    'EWSR1', 'EXOSC8', 'EZH2', 'FAM111B', 'FAM122B', 'FAM72C',
    'FAM72D', 'FAM83D', 'FANCB', 'FANCD2', 'FANCI', 'FANCL', 'FBL',
    'FBXL5', 'FBXO5', 'FDPS', 'FDX1', 'FEN1', 'FGFR1OP', 'FILIP1L',
    'FOXM1', 'FUS', 'FUZ', 'FXR1', 'FZR1', 'G2E3', 'G3BP1', 'GABPB1',
    'GAS2L3', 'GEMIN2', 'GEN1', 'GGCT', 'GGH', 'GINS2',
    'GLO1', 'GMNN', 'GMPS', 'GNG5', 'GPBP1', 'GPSM2', 'GTSE1', 'H1FX',
    'H2AFV', 'H2AFX', 'H2AFY', 'H2AFZ', 'HACD3', 'HADH', 'HAT1',
    'HAUS1', 'HAUS6', 'HAUS8', 'HDAC2', 'HDGF', 'HELLS', 'HES1',
    'HINT1', 'HIRIP3', 'H1-1', 'H1-2', 'H1-3', 'H2BC9',
    'HIST1H4C', 'H2AC6', 'HJURP', 'HMG20B', 'HMGA1', 'HMGA2',
    'HMGB1', 'HMGB2', 'HMGB3', 'HMGN1', 'HMGN2', 'HMGN3', 'HMGN5',
    'HMGXB4', 'HMMR', 'HNRNPA0', 'HNRNPA1', 'HNRNPA2B1', 'HNRNPA3',
    'HNRNPAB', 'HNRNPC', 'HNRNPD', 'HNRNPDL', 'HNRNPF', 'HNRNPH3',
    'HNRNPK', 'HNRNPLL', 'HNRNPM', 'HNRNPU', 'HNRNPUL1', 'HP1BP3',
    'HPF1', 'HSD17B11', 'HSP90B1', 'HSPA13', 'HSPA1B',
    'HSPB11', 'HSPD1', 'HSPE1', 'HYLS1', 'IDH2', 'IFT122', 'IGF2BP3',
    'IKBIP', 'ILF2', 'ILF3', 'ILVBL', 'IMMP1L', 'INCENP', 'IPO5',
    'IQGAP3', 'ISCA2', 'ISOC1', 'ITGAE', 'ITGB3BP', 'JADE1', 'JPT1',
    'KATNBL1', 'KCTD9', 'KIAA0586', 'KIF11', 'KIF14', 'KIF15',
    'KIF18A', 'KIF18B', 'KIF20A', 'KIF20B', 'KIF22', 'KIF23', 'KIF2C',
    'KIF4A', 'KIF5B', 'KIFC1', 'KMT5A', 'KNL1', 'KNSTRN', 'KPNA2',
    'KPNB1', 'LARP7', 'LBR', 'LCORL', 'LIG1', 'LIN52',
    'LINC01224', 'LINC01572', 'LMNB1', 'LMNB2', 'LRR1', 'LSM14A',
    'LSM2', 'LSM3', 'LSM4', 'LSM5', 'LSM6', 'LSM7', 'LSM8', 'LUC7L2',
    'MAD2L1', 'MAGI1', 'MAGOH', 'MAGOHB', 'MAPK1IP1L', 'MAPRE1',
    'MARCKS', 'MASTL', 'MBNL2', 'MCM10', 'MCM2', 'MCM3', 'MCM4',
    'MCM5', 'MCM7', 'MED30', 'MELK', 'MGME1', 'MIS18A',
    'MIS18BP1', 'MKI67', 'MMS22L', 'MND1', 'MNS1', 'MORF4L2',
    'MPHOSPH9', 'MRE11', 'MRPL18', 'MRPL23', 'MRPL47', 'MRPL51',
    'MRPL57', 'MRPS34', 'MTFR2', 'MYBL2', 'MYEF2', 'MZT1', 'MZT2B',
    'NAA38', 'NAA50', 'NAE1', 'NAP1L1', 'NAP1L4', 'NASP', 'NCAPD2',
    'NCAPD3', 'NCAPG', 'NCAPG2', 'NCAPH', 'NCL', 'NDC1', 'NDC80',
    'NDE1', 'NDUFA6', 'NDUFAF3', 'NDUFS6', 'NEDD1', 'NEIL3', 'NEK2',
    'NELFE', 'NENF', 'NFATC3', 'NFYB', 'NIPBL', 'NMU', 'NONO', 'NOP56',
    'NOP58', 'NRDC', 'NSD2', 'NSMCE2', 'NSMCE4A', 'NUCKS1', 'NUDC',
    'NUDCD2', 'NUDT1', 'NUDT15', 'NUDT21', 'NUDT5', 'NUF2', 'NUP107',
    'NUP35', 'NUP37', 'NUP50', 'NUP54', 'NUSAP1', 'ODC1', 'ODF2',
    'OIP5', 'ORC6', 'PA2G4', 'PAICS', 'PAIP2', 'PAK4', 'PAPOLA',
    'PARP1', 'PARPBP', 'PAXX', 'PBK', 'PCBD2', 'PCBP2', 'PCM1', 'PCNA',
    'PCNP', 'PDS5B', 'PHF19', 'PHF5A', 'PHGDH', 'PHIP', 'PIF1',
    'PIMREG', 'PIN1', 'PLCB1', 'PLGRKT', 'PLIN3', 'PLK1',
    'PLK4', 'PMAIP1', 'PNISR', 'PNN', 'PNRC2', 'POC1A', 'POLD2',
    'POLD3', 'POLE2', 'POLQ', 'POLR2C', 'POLR2D', 'POLR2G', 'POLR2J',
    'POLR2K', 'POLR3K', 'PPIA', 'PPIG', 'PPIH', 'PPP1CC', 'PPP2R3C',
    'PPP2R5C', 'PPP6R3', 'PRC1', 'PRDX3', 'PRIM1', 'PRIM2',
    'PRPF38B', 'PRPSAP1', 'PRR11', 'PSIP1', 'PSMA3', 'PSMA4', 'PSMB2',
    'PSMB3', 'PSMC3', 'PSMC3IP', 'PSMD10', 'PSMD14', 'PSMG2', 'PSRC1',
    'PTBP1', 'PTGES3', 'PTMA', 'PTMS', 'PTTG1', 'PUF60', 'RAB8A',
    'RACGAP1', 'RAD21', 'RAD51AP1', 'RAD51B', 'RAD51C', 'RAN',
    'RANBP1', 'RANGAP1', 'RASSF1', 'RBBP4', 'RBBP8', 'RBL1', 'RBM17',
    'RBM39', 'RBM8A', 'RBMX', 'RCC1', 'RDX', 'REEP4', 'RFC1', 'RFC2',
    'RFC3', 'RFC4', 'RFWD3', 'RHEB', 'RMI2', 'RNASEH2B', 'RNASEH2C',
    'RNF138', 'RNF168', 'RNF26', 'RNPS1', 'RPA1', 'RPA3', 'RPL35',
    'RPL39L', 'RPLP0', 'RPLP1', 'RPLP2', 'RPN2', 'RPP30', 'RPS15',
    'RPS16', 'RPS20', 'RPS21', 'RPSA', 'RRM1', 'RSRC1', 'RSRC2',
    'RTKN2', 'RUVBL2', 'SAC3D1', 'SAE1', 'SAP18', 'SAPCD2', 'SCAF11',
    'SCLT1', 'SDHAF3', 'SELENOK', 'SEM1', 'SEPHS1',
    'SERBP1', 'SET', 'SF1', 'SF3B2', 'SFPQ', 'SGO1', 'SGO2',
    'SHCBP1', 'SINHCAF', 'SIVA1', 'SKA1', 'SKA2', 'SKA3', 'SLBP',
    'SLC20A1', 'SLC25A3', 'SLTM', 'SMC1A', 'SMC2', 'SMC3', 'SMC4',
    'SMC5', 'SMCHD1', 'SNAPC1', 'SNRNP25', 'SNRNP40', 'SNRNP70',
    'SNRPA', 'SNRPA1', 'SNRPB', 'SNRPC', 'SNRPD1', 'SNRPD2', 'SNRPD3',
    'SNRPE', 'SNRPF', 'SNRPG', 'SON', 'SPAG5', 'SPATA5', 'SPC25',
    'SPCS2', 'SPDL1', 'SREK1', 'SRI', 'SRP9', 'SRRM1', 'SRSF1',
    'SRSF10', 'SRSF11', 'SRSF2', 'SRSF3', 'SRSF4', 'SRSF7', 'SSB',
    'SSBP1', 'SSNA1', 'SSRP1', 'ST13', 'STAG1', 'STIL', 'STIP1',
    'STK17B', 'STK3', 'STOML2', 'SUGP2', 'SUMO1', 'SUMO3', 'SUPT16H',
    'SUV39H2', 'SUZ12', 'SYNE2', 'TACC3', 'TBC1D31', 'TBC1D5', 'TDP1',
    'TEAD1', 'TEX30', 'TFDP1', 'THRAP3', 'TICRR', 'TIMELESS', 'TIMM10',
    'TK1', 'TMED5', 'TMEM106C', 'TMEM237', 'TMEM60', 'TMEM97', 'TMPO',
    'TMSB15A', 'TOP1', 'TOP2A', 'TPI1', 'TPR', 'TPRKB', 'TPX2',
    'TRA2B', 'TRAIP', 'TROAP', 'TTC28', 'TTF2', 'TTK', 'TXNDC12', 'TYMS',
    'UBA2', 'UBB', 'UBE2C', 'UBE2D2', 'UBE2D3', 'UBE2I', 'UBE2N',
    'UBE2S', 'UBE2T', 'UHRF1', 'UNG', 'UQCC2', 'UQCC3', 'UQCRC1',
    'UQCRFS1', 'USP1', 'VBP1', 'VDAC3', 'VEZF1', 'VRK1', 'WAPL',
    'WDHD1', 'WDPCP', 'WDR34', 'WDR76', 'XPO1', 'XRCC4', 'XRCC5',
    'XRCC6', 'YAP1', 'YBX1', 'YEATS4', 'Z94721.1', 'ZFP36L1', 'ZGRF1',
    'ZMYM1', 'ZNF22', 'ZNF367', 'ZNF43', 'ZNF704', 'ZNF83', 'ZRANB3',
    'ZSCAN16-AS1', 'ZWINT'
], dtype=object)

G1_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
    "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP",
    "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
    "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
    "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]

S_GENES = [
    "H2AC14", "H2AC17", "H1-3", "H4C3", "HIST1H2AJ", "HIST1H2AM",
    "HIST1H1D", "HIST1H4C",
]

G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3",
    "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1",
    "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "HJURP",
    "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
    "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
    "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2",
    "G2E3", "GAS2L3", "CBX5", "CENPA",
]

# Thresholds for --scoring fraction (fractions, not percentages).
CYCLING_THRESHOLD = 0.004     # 0.4%
G1_THRESHOLD      = 0.002   # 0.2%
S_THRESHOLD       = 0.002   # 0.2%
G2M_THRESHOLD     = 0.03    # 3%

# Thresholds for --scoring background, in standard deviations of the
# background-normalized score. Comparable across gene sets by construction.
BG_CYCLING_THRESHOLD = 1.0
BG_G1_THRESHOLD      = 1.0
BG_S_THRESHOLD       = 1.0
BG_G2M_THRESHOLD     = 1.0

# Control-matching parameters for --scoring background.
BG_N_BINS            = 25     # equal-count signal bins features are matched within
BG_N_CTRL_PER_REGION = 50     # controls drawn per set feature
BG_SEED              = 0
BG_N_BACKGROUND_GENES = 5000  # [--fragments] non-cc genes whose promoters form the pool

PHASE_ORDER = ["Non-cycling", "G1", "S", "G2M", "Post-M"]
SCORING_MODES = ("fraction", "background")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_ENSEMBL_RE = re.compile(r"^ENS[A-Z]*G\d+(\.\d+)?$")


def _looks_like_ensembl(names: list[str], sample_size: int = 20) -> bool:
    sample = list(names)[:sample_size]
    if not sample:
        return False
    return sum(1 for s in sample if _ENSEMBL_RE.match(str(s))) >= max(1, len(sample) // 2)


def _symbol_to_var_map(adata, sym_col: str) -> dict[str, str]:
    """Return {gene_symbol: var_names_entry}. Uses adata.var[sym_col] if the
    column exists AND contains non-Ensembl-looking values; otherwise falls back
    to adata.var_names."""
    if sym_col in adata.var.columns:
        symbols = adata.var[sym_col].astype(str).tolist()
        if not _looks_like_ensembl(symbols):
            return dict(zip(symbols, adata.var_names.astype(str)))
        print(f"  [warn] adata.var['{sym_col}'] looks like Ensembl IDs.")
    if _looks_like_ensembl(adata.var_names.astype(str).tolist()):
        print(
            "  [warn] adata.var_names look like Ensembl IDs. Cell-cycle gene "
            "matching will likely find 0 genes. Convert your h5ad's var_names "
            "or add a symbol column before running."
        )
    return {g: g for g in adata.var_names.astype(str)}


def _gene_indices(gene_list, sym2var: dict[str, str],
                  var_pos: dict[str, int], label: str) -> np.ndarray:
    """Positional indices into adata.var for genes present in the AnnData.
    Dedupes and reports missing genes."""
    unique_requested = list(dict.fromkeys(gene_list))   # preserve order + dedupe
    ix = []
    found_names, missing = [], []
    for g in unique_requested:
        var_name = sym2var.get(g)
        if var_name is not None and var_name in var_pos:
            ix.append(var_pos[var_name])
            found_names.append(g)
        else:
            missing.append(g)
    print(f"  {label:15s}  {len(found_names):>4}/{len(unique_requested):>4} genes found in adata")
    if missing and len(missing) <= 10:
        print(f"                   missing: {missing}")
    return np.asarray(ix, dtype=int)


def _per_cell_fraction(X, gene_ix: np.ndarray, total_umi: np.ndarray) -> np.ndarray:
    """(sum X[:, gene_ix] per cell) / total_umi. Zero total → 0."""
    if len(gene_ix) == 0:
        return np.zeros(X.shape[0])
    import scipy.sparse as sp
    if sp.issparse(X):
        sub_sum = np.asarray(X[:, gene_ix].sum(axis=1)).ravel()
    else:
        sub_sum = np.asarray(X[:, gene_ix]).sum(axis=1)
    denom = np.where(total_umi > 0, total_umi, 1.0)
    return sub_sum / denom


# ---------------------------------------------------------------------------
# Background-normalized scoring
#
# The fraction scheme above answers "what share of this cell's signal sits in
# the set?", which is dominated by how many features the set has and by their
# baseline level. Here each set is instead compared against control features
# matched on mean signal, so the score reflects only how the set deviates from
# a comparable background in that cell.
# ---------------------------------------------------------------------------

def _rate_matrix(X, total):
    """Depth-normalize counts: entry (i, j) is the share of cell i's total
    signal that falls in feature j. Cells with zero total become all-zero."""
    import scipy.sparse as sp
    safe  = np.where(total > 0, total, 1.0)
    inv   = np.where(total > 0, 1.0 / safe, 0.0)
    if sp.issparse(X):
        return sp.diags(inv) @ X.tocsr()
    return np.asarray(X, dtype=float) * inv[:, None]


def _weighted_feature_mean(X, ix, n_features):
    """Per-cell mean over the features in `ix`. Repeated indices act as
    weights, which is what keeps a matched control pool weighted by the set's
    bin composition rather than by bin size."""
    if len(ix) == 0:
        return np.zeros(X.shape[0])
    w = np.bincount(np.asarray(ix, dtype=int), minlength=n_features).astype(float)
    total_w = w.sum()
    if total_w == 0:
        return np.zeros(X.shape[0])
    w /= total_w
    return np.asarray(X @ w).ravel()


def _signal_bins(mean_signal, n_bins: int) -> np.ndarray:
    """Assign each feature to one of `n_bins` equal-count bins ordered by mean
    signal. Rank-based rather than value-based, so sparse ATAC features (where
    most means pile up near zero) still spread across bins."""
    n = len(mean_signal)
    if n == 0:
        return np.asarray([], dtype=int)
    order = np.argsort(mean_signal, kind="mergesort")
    rank  = np.empty(n, dtype=int)
    rank[order] = np.arange(n)
    return np.minimum((rank * n_bins) // n, n_bins - 1)


def _matched_controls(bin_id, pools, set_ix, n_ctrl_per_region, rng, with_null=False):
    """Draw control features matched to `set_ix` on signal bin.

    For each set feature we sample `n_ctrl_per_region` controls from that
    feature's own bin. With `with_null`, one extra feature per set feature is
    drawn in the same call and split off, giving a null set that is the same
    size as the real set, matched the same way, and guaranteed disjoint from
    the controls it will be compared against.

    Returns (ctrl_ix, null_ix) as concatenated index arrays; duplicates across
    set features are intentional and carry the matching weights.
    """
    ctrl, null = [], []
    extra = 1 if with_null else 0
    for j in np.asarray(set_ix, dtype=int):
        pool = pools.get(int(bin_id[j]))
        if pool is None or len(pool) == 0:
            continue
        k = min(n_ctrl_per_region + extra, len(pool))
        pick = rng.choice(pool, size=k, replace=False)
        if with_null and k > 1:
            null.append(pick[:1])
            ctrl.append(pick[1:])
        else:
            ctrl.append(pick)
    ctrl_ix = np.concatenate(ctrl) if ctrl else np.asarray([], dtype=int)
    null_ix = np.concatenate(null) if null else np.asarray([], dtype=int)
    return ctrl_ix, null_ix


def _background_normalized_scores(
    X, total, set_indices: dict[str, np.ndarray], background_mask: np.ndarray,
    n_bins: int = BG_N_BINS, n_ctrl_per_region: int = BG_N_CTRL_PER_REGION,
    seed: int = BG_SEED,
) -> dict[str, dict]:
    """Score each gene set against signal-matched control features.

    `background_mask` marks features eligible as controls — everything that is
    not part of any cell-cycle set. Returns, per set name:
        z             z-scored score across cells (what gets thresholded)
        raw           mean(set) - mean(matched controls), depth-normalized
        signal_ratio  sd(raw) / sd(null), where null replaces the set with a
                      size- and signal-matched random set. <= 1 means the set
                      carries no more cell-to-cell structure than chance.
        n_features    features in the set
        n_ctrl        distinct control features drawn
    """
    n_features  = X.shape[1]
    rate        = _rate_matrix(X, total)
    mean_signal = np.asarray(rate.mean(axis=0)).ravel()
    bin_id      = _signal_bins(mean_signal, n_bins)

    eligible = np.asarray(background_mask, dtype=bool)
    pools    = {b: np.flatnonzero(eligible & (bin_id == b)) for b in range(n_bins)}
    n_pool   = int(eligible.sum())
    if n_pool == 0:
        raise RuntimeError(
            "No features are eligible as controls. Background-normalized "
            "scoring needs features outside the cell-cycle sets to compare "
            "against (for --fragments, raise --n-background-genes)."
        )
    rng = np.random.default_rng(seed)

    out = {}
    for name, set_ix in set_indices.items():
        set_ix = np.asarray(set_ix, dtype=int)
        if len(set_ix) == 0:
            zeros = np.zeros(X.shape[0])
            out[name] = {"z": zeros, "raw": zeros, "signal_ratio": float("nan"),
                         "n_features": 0, "n_ctrl": 0}
            continue

        ctrl_ix, null_ix = _matched_controls(
            bin_id, pools, set_ix, n_ctrl_per_region, rng, with_null=True)
        ctrl_mean = _weighted_feature_mean(rate, ctrl_ix, n_features)
        raw  = _weighted_feature_mean(rate, set_ix,  n_features) - ctrl_mean
        null = _weighted_feature_mean(rate, null_ix, n_features) - ctrl_mean

        sd      = float(np.nanstd(raw))
        null_sd = float(np.nanstd(null))
        z = (raw - float(np.nanmean(raw))) / (sd if sd > 0 else 1.0)
        out[name] = {
            "z": z,
            "raw": raw,
            "signal_ratio": (sd / null_sd) if null_sd > 0 else float("nan"),
            "n_features": int(len(set_ix)),
            "n_ctrl": int(len(np.unique(ctrl_ix))),
        }

    print(f"  Control pool: {n_pool:,} features in {n_bins} signal-matched bins "
          f"({n_ctrl_per_region} controls per set feature)")
    return out


def _print_background_table(scores: dict[str, dict], thresholds: dict[str, float]):
    def _pct(a, q): return float(np.nanpercentile(a, q))
    print("\n  Background-normalized scores (z across cells):")
    print(f"    {'set':13s}  {'threshold':>9s}  {'n_feat':>7s}  "
          f"{'p10':>9s}  {'median':>9s}  {'p90':>9s}  {'max':>9s}  {'signal':>7s}")
    weak = []
    for name, res in scores.items():
        z = res["z"]
        print(f"    {name:13s}  {thresholds[name]:>9.3g}  {res['n_features']:>7,}  "
              f"{_pct(z,10):>9.3g}  {_pct(z,50):>9.3g}  {_pct(z,90):>9.3g}  "
              f"{float(np.nanmax(z)):>9.3g}  {res['signal_ratio']:>7.2f}")
        if not np.isnan(res["signal_ratio"]) and res["signal_ratio"] <= 1.1:
            weak.append((name, res["signal_ratio"]))
    for name, ratio in weak:
        print(f"    [warn] '{name}' signal={ratio:.2f} — no more cell-to-cell structure "
              f"than a size- and signal-matched random set.")
        print(f"           Phase calls driven by this set are not interpretable.")


# ---------------------------------------------------------------------------
# Shared plumbing
# ---------------------------------------------------------------------------

def _background_mask_from_sets(n_features: int, set_indices: dict) -> np.ndarray:
    """Features eligible as controls: everything in no cell-cycle set."""
    used = np.zeros(n_features, dtype=bool)
    for ix in set_indices.values():
        ix = np.asarray(ix, dtype=int)
        if len(ix):
            used[ix] = True
    return ~used


def _resolve_thresholds(scoring: str, cycling, g1, s, g2m) -> dict[str, float]:
    """Map the four threshold arguments onto the score names, filling in the
    defaults for `scoring` where an argument is None."""
    defaults = (
        (BG_CYCLING_THRESHOLD, BG_G1_THRESHOLD, BG_S_THRESHOLD, BG_G2M_THRESHOLD)
        if scoring == "background" else
        (CYCLING_THRESHOLD, G1_THRESHOLD, S_THRESHOLD, G2M_THRESHOLD)
    )
    given = (cycling, g1, s, g2m)
    names = ("cycling_score", "g1_frac", "s_frac", "g2m_frac")
    return {n: (d if v is None else float(v))
            for n, v, d in zip(names, given, defaults)}


def _score_sets(X, total, set_indices: dict, scoring: str, thresholds: dict,
                background_mask=None, n_bins=BG_N_BINS,
                n_ctrl_per_region=BG_N_CTRL_PER_REGION, seed=BG_SEED):
    """Compute the four per-cell scores under `scoring` and print the matching
    diagnostic table. Returns {score_name: per-cell array}."""
    if scoring not in SCORING_MODES:
        raise ValueError(f"scoring must be one of {SCORING_MODES}, got {scoring!r}")

    if scoring == "fraction":
        print("  Computing per-cell fractions ...")
        values = {n: _per_cell_fraction(X, ix, total) for n, ix in set_indices.items()}
        _print_fraction_table([(n, values[n], thresholds[n]) for n in set_indices])
        return values

    print("  Computing background-normalized scores ...")
    if background_mask is None:
        background_mask = _background_mask_from_sets(X.shape[1], set_indices)
    scores = _background_normalized_scores(
        X, total, set_indices, background_mask,
        n_bins=n_bins, n_ctrl_per_region=n_ctrl_per_region, seed=seed)
    _print_background_table(scores, thresholds)
    return {n: res["z"] for n, res in scores.items()}


def _write_score_obs(adata, scoring: str, values: dict, tss_window=None):
    """Write the score columns. `cycling_score` keeps its name in both modes —
    a fraction under `fraction`, a z-score under `background` — and uns records
    which, so downstream code can tell them apart."""
    adata.obs["cycling_score"] = np.asarray(values["cycling_score"], dtype=np.float32)
    for key, col in (("g1_frac", "cc_score_g1"),
                     ("s_frac",  "cc_score_s"),
                     ("g2m_frac", "cc_score_g2m")):
        adata.obs[col] = np.asarray(values[key], dtype=np.float32)
    adata.uns["cell_cycle_scoring"] = {
        "scoring": scoring,
        "score_units": "z-score" if scoring == "background" else "fraction of total signal",
    }
    if tss_window is not None:
        adata.uns["cell_cycle_scoring"]["tss_window"] = int(tss_window)


def _report_phases(phases, n_obs: int, title: str, n_missing: int = 0):
    counts = pd.Series(phases).value_counts()
    print(f"\n  {title}")
    for phase in PHASE_ORDER:
        n = int(counts.get(phase, 0))
        pct = 100.0 * n / n_obs if n_obs else 0.0
        print(f"    {phase:12s}  {n:>10,}  ({pct:5.2f}%)")
    if n_missing:
        print(f"    {'(no data)':12s}  {n_missing:>10,}")


# ---------------------------------------------------------------------------
# Main annotation routine
# ---------------------------------------------------------------------------

def annotate_cell_cycle(
    adata,
    sym_col: str = "Gene",
    scoring: str = "fraction",
    cycling_threshold: float | None = None,
    g1_threshold: float | None  = None,
    s_threshold: float | None   = None,
    g2m_threshold: float | None = None,
    n_bins: int = BG_N_BINS,
    n_ctrl_per_region: int = BG_N_CTRL_PER_REGION,
    seed: int = BG_SEED,
):
    """Add `cycling_score` and `CellCyclePhase` columns to adata.obs (in place).

    Under `scoring="background"` each gene set is scored against expression-
    matched control genes instead of as a share of the cell's UMIs, and the
    thresholds are in standard deviations. See the module docstring.
    """
    import scipy.sparse as sp

    thresholds = _resolve_thresholds(
        scoring, cycling_threshold, g1_threshold, s_threshold, g2m_threshold)

    print(f"  Building symbol -> var lookup (sym_col='{sym_col}')")
    sym2var = _symbol_to_var_map(adata, sym_col)
    var_pos = {v: i for i, v in enumerate(adata.var_names.astype(str))}

    print("  Resolving gene sets:")
    cc_ix  = _gene_indices(CC_GENES_HUMAN, sym2var, var_pos, "cc_genes")
    g1_ix  = _gene_indices(G1_GENES,       sym2var, var_pos, "G1")
    s_ix   = _gene_indices(S_GENES,        sym2var, var_pos, "S")
    g2m_ix = _gene_indices(G2M_GENES,      sym2var, var_pos, "G2M")

    if len(cc_ix) == 0:
        raise RuntimeError(
            "None of the cell-cycle gene set was found in adata.var. "
            "Check --sym-col and the gene naming convention in the h5ad."
        )

    X = adata.X
    if sp.issparse(X):
        X = X.tocsr()

    print("  Computing UMI totals ...")
    total_umi   = np.asarray(X.sum(axis=1)).ravel()
    set_indices = {"cycling_score": cc_ix, "g1_frac": g1_ix,
                   "s_frac": s_ix, "g2m_frac": g2m_ix}
    values = _score_sets(X, total_umi, set_indices, scoring, thresholds,
                         n_bins=n_bins, n_ctrl_per_region=n_ctrl_per_region,
                         seed=seed)

    phases = _classify_phases(
        values["cycling_score"], values["g1_frac"], values["s_frac"],
        values["g2m_frac"], thresholds["cycling_score"], thresholds["g1_frac"],
        thresholds["s_frac"], thresholds["g2m_frac"],
    )

    _write_score_obs(adata, scoring, values)
    adata.obs["CellCyclePhase"] = pd.Categorical(
        phases, categories=PHASE_ORDER, ordered=False,
    )
    _report_phases(phases, adata.n_obs, "Phase distribution:")
    return adata


# ---------------------------------------------------------------------------
# ATAC-seq version — same logic, but "UMI fraction" becomes "peak accessibility
# fraction over cell-cycle gene promoters (TSS ± window)".
# ---------------------------------------------------------------------------

def _peak_indices_for_gene_set(
    gene_list, tss_df, peaks_df, peak_to_ix, tss_window: int, label: str,
) -> np.ndarray:
    """Positional indices into adata.var of peaks overlapping any promoter
    window (TSS ± tss_window) of a gene in `gene_list`. Prints a two-stage
    diagnostic so you can see whether losses happen at TSS match or peak
    overlap."""
    from atac_seq_analysis import build_promoter_regions, peaks_overlapping_promoters

    unique_requested = list(dict.fromkeys(gene_list))
    prom_df = build_promoter_regions(tss_df, tss_window, genes=unique_requested)
    n_by_tss = int(prom_df["gene_name"].nunique()) if not prom_df.empty else 0

    if prom_df.empty:
        missing = unique_requested
        print(f"  {label:15s}  0/{len(unique_requested):>4} matched a TSS; "
              f"missing genes: {missing[:12]}{' ...' if len(missing) > 12 else ''}")
        return np.asarray([], dtype=int)

    overlap_df = peaks_overlapping_promoters(peaks_df, prom_df)
    n_by_peak  = int(overlap_df["gene_name"].nunique()) if not overlap_df.empty else 0
    peak_ids   = overlap_df["peak_id"].astype(str).unique().tolist() if not overlap_df.empty else []
    ix         = [peak_to_ix[p] for p in peak_ids if p in peak_to_ix]

    print(f"  {label:15s}  "
          f"TSS-matched: {n_by_tss:>3}/{len(unique_requested):<3}  "
          f"→ peak-overlap: {n_by_peak:>3}  "
          f"→ pooled peaks: {len(ix):>5,}")

    if n_by_tss > 0 and n_by_peak == 0:
        genes_in_prom = set(prom_df["gene_name"].astype(str))
        print(f"    [note] {n_by_tss} promoter(s) built but 0 ATAC peaks overlap them.")
        print(f"           TSS-matched genes: {sorted(genes_in_prom)}")
    elif len(unique_requested) - n_by_tss > 0:
        matched_set = set(prom_df["gene_name"].astype(str))
        missing = [g for g in unique_requested if g not in matched_set]
        print(f"    [note] {len(missing)} gene(s) missing from GTF: "
              f"{missing[:10]}{' ...' if len(missing) > 10 else ''}")
    return np.asarray(ix, dtype=int)


def annotate_cell_cycle_atac(
    adata,
    gtf_path: str | Path | None = None,
    genome: str = "hg38",
    tss_window: int = 2000,
    peak_columns: dict | None = None,
    scoring: str = "fraction",
    cycling_threshold: float | None = None,
    g1_threshold: float | None  = None,
    s_threshold: float | None   = None,
    g2m_threshold: float | None = None,
    n_bins: int = BG_N_BINS,
    n_ctrl_per_region: int = BG_N_CTRL_PER_REGION,
    seed: int = BG_SEED,
):
    """Add `cycling_score` and `CellCyclePhase` obs columns to an ATAC h5ad.

    For each cell-cycle gene set, finds every ATAC peak overlapping a
    (TSS ± tss_window) promoter of a gene in the set. Under `scoring="fraction"`
    the "cycling score" is the total accessibility of those peaks divided by
    the cell's total ATAC signal; under `scoring="background"` it is the set's
    mean accessibility minus that of accessibility-matched control peaks, as a
    z-score. Same priority ordering (G2M > S > G1 > Post-M) as the RNA path.
    """
    import scipy.sparse as sp
    from atac_seq_analysis import peaks_to_dataframe, load_tss_annotation

    thresholds = _resolve_thresholds(
        scoring, cycling_threshold, g1_threshold, s_threshold, g2m_threshold)

    print("  Parsing peak coordinates from adata.var ...")
    peaks_df   = peaks_to_dataframe(adata, peak_columns)
    peak_to_ix = {pid: i for i, pid in enumerate(peaks_df["peak_id"].astype(str))}
    print(f"    {len(peaks_df):,} peaks parsed")

    print(f"  Loading TSS annotations (genome={genome}, gtf_path={gtf_path}) ...")
    tss_df = load_tss_annotation(gtf_path, genome)
    print(f"    {len(tss_df):,} gene TSSes")

    print(f"  Building per-gene-set peak indices (TSS ± {tss_window} bp):")
    cc_ix  = _peak_indices_for_gene_set(CC_GENES_HUMAN, tss_df, peaks_df, peak_to_ix, tss_window, "cc_genes")
    g1_ix  = _peak_indices_for_gene_set(G1_GENES,       tss_df, peaks_df, peak_to_ix, tss_window, "G1")
    s_ix   = _peak_indices_for_gene_set(S_GENES,        tss_df, peaks_df, peak_to_ix, tss_window, "S")
    g2m_ix = _peak_indices_for_gene_set(G2M_GENES,      tss_df, peaks_df, peak_to_ix, tss_window, "G2M")

    if len(cc_ix) == 0:
        raise RuntimeError(
            "0 peaks overlap any cell-cycle gene promoter. Common causes: wrong "
            "gtf_path/genome, peak var_names in a different coordinate style, or "
            "gene names in the GTF using Ensembl IDs. Check the log above."
        )

    X = adata.X
    if sp.issparse(X):
        X = X.tocsr()

    total = np.asarray(X.sum(axis=1)).ravel()
    set_indices = {"cycling_score": cc_ix, "g1_frac": g1_ix,
                   "s_frac": s_ix, "g2m_frac": g2m_ix}
    values = _score_sets(X, total, set_indices, scoring, thresholds,
                         n_bins=n_bins, n_ctrl_per_region=n_ctrl_per_region,
                         seed=seed)

    phases = _classify_phases(
        values["cycling_score"], values["g1_frac"], values["s_frac"],
        values["g2m_frac"], thresholds["cycling_score"], thresholds["g1_frac"],
        thresholds["s_frac"], thresholds["g2m_frac"],
    )

    _write_score_obs(adata, scoring, values, tss_window=tss_window)
    adata.obs["CellCyclePhase"] = pd.Categorical(
        phases, categories=PHASE_ORDER, ordered=False,
    )
    _report_phases(phases, adata.n_obs, "Phase distribution (ATAC):")
    return adata


# ---------------------------------------------------------------------------
# ATAC-seq, fragment-based version — count Tn5 fragments in promoter windows
# (TSS ± window) directly from a fragments.tsv(.bgz) file via snapatac2, so
# genes without a called peak still contribute. Denominator per cell is the
# total number of fragments (snapatac2's `n_fragment`).
# ---------------------------------------------------------------------------

_MAIN_CHROM_RE = re.compile(r"^(chr)?([0-9]+|X|Y)$")


def _chrom_style(name: str) -> str | None:
    """'ucsc' for chr1/chrX, 'ensembl' for 1/X, None for contigs like GL000009.2."""
    m = _MAIN_CHROM_RE.match(str(name))
    if m is None:
        return None
    return "ucsc" if m.group(1) else "ensembl"


def _to_ensembl_chrom(name: str) -> str:
    s = str(name)
    if s == "chrM":
        return "MT"
    return s[3:] if s.startswith("chr") else s


def _open_text(path):
    import gzip
    return gzip.open(path, "rt") if str(path).endswith((".gz", ".bgz")) else open(path, "rt")


def _detect_fragments_chrom_style(fragments_path, max_lines: int = 5_000_000) -> str:
    """Scan the fragments file until a main chromosome (1..22/X/Y, with or
    without 'chr') is seen. Unplaced contigs (GL*/KI*) are often first in
    coordinate-sorted files, so we keep scanning past them."""
    with _open_text(fragments_path) as f:
        for i, line in enumerate(f):
            if line.startswith("#"):
                continue
            style = _chrom_style(line.split("\t", 1)[0])
            if style is not None:
                return style
            if i >= max_lines:
                break
    print("  [warn] Could not detect chromosome naming from fragments file; assuming 'ucsc'.")
    return "ucsc"


def _promoter_regions_for_genes(tss_df, genes, window: int, to_frag_chrom, allowed_chroms):
    """Build merged TSS ± window intervals per gene, restricted to `genes` and
    to chromosomes present in `allowed_chroms` (after renaming with
    `to_frag_chrom`). Returns (regions, region_gene): region strings
    "chrom:start-end" in snapatac2 format and the gene of each region."""
    df = tss_df[tss_df["gene_name"].isin(set(genes))].copy()
    df["Chromosome"] = df["Chromosome"].astype(str).map(to_frag_chrom)
    df = df[df["Chromosome"].isin(allowed_chroms)]
    df["Start"] = (df["TSS"].astype(int) - window).clip(lower=0)
    df["End"]   = df["TSS"].astype(int) + window

    regions, region_gene = [], []
    for (gene, chrom), grp in df.groupby(["gene_name", "Chromosome"], sort=True):
        cur_s = cur_e = None
        for st, en in sorted(zip(grp["Start"], grp["End"])):
            if cur_s is None:
                cur_s, cur_e = st, en
            elif st <= cur_e:            # overlapping / touching → merge
                cur_e = max(cur_e, en)
            else:
                regions.append(f"{chrom}:{cur_s}-{cur_e}"); region_gene.append(gene)
                cur_s, cur_e = st, en
        if cur_s is not None:
            regions.append(f"{chrom}:{cur_s}-{cur_e}"); region_gene.append(gene)
    return regions, region_gene


def _parse_region(r: str) -> tuple[str, int, int]:
    """'chrom:start-end' -> (chrom, start, end). rsplit so chrom may contain ':'."""
    chrom, span = str(r).rsplit(":", 1)
    start, end = span.split("-")
    return chrom, int(start), int(end)


def _drop_regions_overlapping(regions, region_gene, blocked_regions):
    """Drop windows that overlap any window in `blocked_regions`.

    With a wide --tss-window a background gene can sit close enough to a
    cell-cycle gene that their promoter windows share fragments. Leaving those
    in the control pool would put set signal into the background and flatten
    the very contrast the score is measuring, so they are removed.
    """
    import bisect

    by_chrom: dict[str, list[list[int]]] = {}
    for r in blocked_regions:
        c, s, e = _parse_region(r)
        by_chrom.setdefault(c, []).append([s, e])
    for c, ivs in by_chrom.items():          # merge into disjoint, sorted intervals
        ivs.sort()
        merged = []
        for s, e in ivs:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        by_chrom[c] = merged
    starts = {c: [iv[0] for iv in ivs] for c, ivs in by_chrom.items()}

    keep_r, keep_g = [], []
    for r, g in zip(regions, region_gene):
        c, s, e = _parse_region(r)
        ivs = by_chrom.get(c)
        if ivs:
            # last interval starting before this window ends; intervals are
            # disjoint and sorted, so if that one clears `s`, all earlier do too
            j = bisect.bisect_left(starts[c], e) - 1
            if j >= 0 and ivs[j][1] > s:
                continue
        keep_r.append(r)
        keep_g.append(g)
    return keep_r, keep_g


def _classify_phases(cycling_score, g1_frac, s_frac, g2m_frac,
                     cycling_threshold, g1_threshold, s_threshold, g2m_threshold):
    """Same rule as the RNA / peak paths (priority: G2M > S > G1 > Post-M).
    NaN scores (cells without data) get None."""
    n = len(cycling_score)
    valid = ~np.isnan(cycling_score)
    is_cycling = valid & (cycling_score > cycling_threshold)
    passes_g1  = np.nan_to_num(g1_frac)  > g1_threshold
    passes_s   = np.nan_to_num(s_frac)   > s_threshold
    passes_g2m = np.nan_to_num(g2m_frac) > g2m_threshold

    phases = np.full(n, None, dtype=object)
    phases[valid] = "Non-cycling"
    phases[is_cycling] = "Post-M"
    phases[is_cycling & passes_g1]  = "G1"
    phases[is_cycling & passes_s]   = "S"
    phases[is_cycling & passes_g2m] = "G2M"
    return phases


def _print_fraction_table(rows):
    def _pct(a, q): return float(np.nanpercentile(a, q))
    print("\n  Per-cell fraction distribution (helps calibrate thresholds):")
    print(f"    {'set':13s}  {'threshold':>10s}  "
          f"{'p10':>10s}  {'median':>10s}  {'p90':>10s}  {'max':>10s}")
    for name, arr, th in rows:
        print(f"    {name:13s}  {th:>10.4g}  "
              f"{_pct(arr,10):>10.4g}  {_pct(arr,50):>10.4g}  "
              f"{_pct(arr,90):>10.4g}  {float(np.nanmax(arr)):>10.4g}")


def annotate_cell_cycle_fragments(
    adata,
    fragments_path: str | Path,
    gtf_path: str | Path | None = None,
    genome: str = "hg38",
    tss_window: int = 2000,
    scoring: str = "fraction",
    cycling_threshold: float | None = None,
    g1_threshold: float | None  = None,
    s_threshold: float | None   = None,
    g2m_threshold: float | None = None,
    n_bins: int = BG_N_BINS,
    n_ctrl_per_region: int = BG_N_CTRL_PER_REGION,
    n_background_genes: int = BG_N_BACKGROUND_GENES,
    seed: int = BG_SEED,
    tempdir: str | Path | None = None,
    n_jobs: int = 8,
):
    """Add `cycling_score` and `CellCyclePhase` obs columns to an ATAC h5ad
    using a fragments file instead of the peak matrix.

    Steps: snapatac2 import_data (restricted to adata.obs_names) → count
    fragments in TSS ± tss_window windows of the cell-cycle genes
    (make_peak_matrix on explicit regions, so it works across snapatac2
    versions) → per-set score → same priority as the other paths. Cells absent
    from the fragments file get NaN score and a NaN phase.

    Under `scoring="background"` the same promoter windows are also built for
    `n_background_genes` randomly chosen non-cell-cycle genes and counted in the
    same pass. Those windows are the control pool the cell-cycle sets are
    matched against, so the background is promoters scored exactly like the
    sets — not peaks or bins of a different kind.
    """
    import tempfile
    import scipy.sparse as sp
    import snapatac2 as snap
    from atac_seq_analysis import get_snapatac2_genome, load_tss_annotation

    thresholds = _resolve_thresholds(
        scoring, cycling_threshold, g1_threshold, s_threshold, g2m_threshold)

    fragments_path = Path(fragments_path)
    if not fragments_path.exists():
        raise FileNotFoundError(f"fragments file not found: {fragments_path}")

    genome_obj = get_snapatac2_genome(genome)

    # --- chromosome naming: fragments file decides; chrom_sizes + regions follow
    frag_style = _detect_fragments_chrom_style(fragments_path)
    print(f"  Chromosome naming in fragments file: '{frag_style}'")
    to_frag_chrom = _to_ensembl_chrom if frag_style == "ensembl" else (lambda c: c)

    chrom_sizes = genome_obj.chrom_sizes
    chrom_sizes = dict(chrom_sizes() if callable(chrom_sizes) else chrom_sizes)
    chrom_sizes = {to_frag_chrom(k): v for k, v in chrom_sizes.items()}

    work_dir = Path(tempdir) if tempdir else Path(tempfile.mkdtemp(prefix="cc_frag_"))
    work_dir.mkdir(parents=True, exist_ok=True)

    # --- promoter windows (TSS ± window) for the cell-cycle genes only, merged
    #     per gene so a fragment is not counted twice across transcripts
    tss_df = load_tss_annotation(gtf_path, genome)
    all_genes = list(dict.fromkeys([*CC_GENES_HUMAN, *G1_GENES, *S_GENES, *G2M_GENES]))
    regions, region_gene = _promoter_regions_for_genes(
        tss_df, all_genes, tss_window, to_frag_chrom, set(chrom_sizes))
    print(f"  {len(regions):,} promoter windows for "
          f"{len(set(region_gene)):,}/{len(all_genes):,} cell-cycle genes (TSS ± {tss_window} bp)")
    if not regions:
        raise RuntimeError("No promoter windows built; check the annotation / gene names.")

    # --- background promoter windows: the control pool for --scoring background.
    #     Same construction as above on non-cell-cycle genes, so controls differ
    #     from the sets only in which genes they cover.
    cc_universe = set(all_genes)
    if scoring == "background":
        candidates = sorted(set(tss_df["gene_name"].astype(str)) - cc_universe)
        if not candidates:
            raise RuntimeError("No non-cell-cycle genes in the annotation to use as controls.")
        rng_genes = np.random.default_rng(seed)
        if n_background_genes and len(candidates) > n_background_genes:
            bg_genes = list(rng_genes.choice(candidates, size=n_background_genes, replace=False))
        else:
            bg_genes = candidates
        bg_regions, bg_region_gene = _promoter_regions_for_genes(
            tss_df, bg_genes, tss_window, to_frag_chrom, set(chrom_sizes))
        n_raw = len(bg_regions)
        bg_regions, bg_region_gene = _drop_regions_overlapping(
            bg_regions, bg_region_gene, regions)
        print(f"  {len(bg_regions):,} background promoter windows for "
              f"{len(set(bg_region_gene)):,} non-cell-cycle genes "
              f"({n_raw - len(bg_regions):,} dropped for overlapping a cell-cycle window)")
        if not bg_regions:
            raise RuntimeError(
                "Every background window overlaps a cell-cycle window — nothing "
                "left as a control. Lower --tss-window or raise --n-background-genes.")
        regions     = regions + bg_regions
        region_gene = region_gene + bg_region_gene

    # --- import fragments for the cells in adata only
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
    n_found = frag.n_obs
    print(f"    {n_found:,}/{len(barcodes):,} cells found in fragments file")
    if n_found == 0:
        raise RuntimeError(
            "No whitelisted barcodes found in the fragments file. Check that "
            "adata.obs_names match the fragment barcodes (column 4) and that "
            "chromosome names match chrom_sizes (see log above).")

    # --- count fragments in each promoter window (cells × regions)
    print("  Counting fragments in promoter windows ...")
    region_mat = snap.pp.make_peak_matrix(frag, use_rep=regions)
    print(f"    {region_mat.n_obs:,} cells × {region_mat.n_vars:,} regions")

    if "n_fragment" in frag.obs.columns:
        total = np.asarray(frag.obs["n_fragment"].values, dtype=float)
    else:
        raise RuntimeError("frag.obs has no 'n_fragment' column; cannot compute per-cell totals.")

    # make_peak_matrix normally preserves frag's cell order, but the totals are
    # the scoring denominator, so align them explicitly rather than assume it.
    frag_obs, reg_obs = pd.Index(frag.obs_names.astype(str)), pd.Index(region_mat.obs_names.astype(str))
    if not frag_obs.equals(reg_obs):
        pos_in_frag = frag_obs.get_indexer(reg_obs)
        if (pos_in_frag < 0).any():
            raise RuntimeError("region matrix contains cells absent from the fragment object.")
        total = total[pos_in_frag]

    X = region_mat.X
    X = X.tocsr() if sp.issparse(X) else np.asarray(X)
    # region_mat.var_names follow the order of `regions`; map each column to its gene
    var_gene = np.asarray(region_gene, dtype=object)
    if region_mat.n_vars != len(regions):
        region_to_gene = dict(zip(regions, region_gene))
        var_gene = np.asarray(
            [region_to_gene.get(r) for r in region_mat.var_names.astype(str)], dtype=object)

    def _set_indices(gene_list, label):
        wanted = set(gene_list)
        ix = np.flatnonzero(np.isin(var_gene, list(wanted)))
        n_genes = len(set(var_gene[ix]))
        print(f"  {label:15s}  {n_genes:>4}/{len(wanted):>4} genes with promoter windows "
              f"({len(ix):,} regions)")
        return ix

    print("  Resolving gene sets:")
    set_indices = {
        "cycling_score": _set_indices(CC_GENES_HUMAN, "cc_genes"),
        "g1_frac":       _set_indices(G1_GENES,       "G1"),
        "s_frac":        _set_indices(S_GENES,        "S"),
        "g2m_frac":      _set_indices(G2M_GENES,      "G2M"),
    }
    if len(set_indices["cycling_score"]) == 0:
        raise RuntimeError("None of the cell-cycle genes have promoter windows; "
                           "check the annotation's gene_name attribute.")

    # Controls are the background-gene windows only: any window whose gene is
    # outside the cell-cycle universe (and that mapped to a gene at all).
    background_mask = np.asarray(
        [g is not None and g not in cc_universe for g in var_gene], dtype=bool)

    values = _score_sets(X, total, set_indices, scoring, thresholds,
                         background_mask=background_mask, n_bins=n_bins,
                         n_ctrl_per_region=n_ctrl_per_region, seed=seed)

    # --- align back to adata.obs (cells missing from fragments → NaN)
    pos_in_adata = pd.Index(adata.obs_names.astype(str)).get_indexer(reg_obs)
    keep = pos_in_adata >= 0
    aligned = {}
    for k, v in values.items():
        arr = np.full(adata.n_obs, np.nan)
        arr[pos_in_adata[keep]] = v[keep]
        aligned[k] = arr
    n_missing = int(np.isnan(aligned["cycling_score"]).sum())
    if n_missing:
        print(f"  [note] {n_missing:,} cells in adata have no fragments → NaN score / NaN phase")

    phases = _classify_phases(
        aligned["cycling_score"], aligned["g1_frac"], aligned["s_frac"], aligned["g2m_frac"],
        thresholds["cycling_score"], thresholds["g1_frac"],
        thresholds["s_frac"], thresholds["g2m_frac"],
    )

    _write_score_obs(adata, scoring, aligned, tss_window=tss_window)
    adata.obs["CellCyclePhase"] = pd.Categorical(
        phases, categories=PHASE_ORDER, ordered=False,
    )
    total_frag = np.full(adata.n_obs, np.nan)
    total_frag[pos_in_adata[keep]] = total[keep]
    adata.obs["n_fragment"] = total_frag.astype(np.float32)

    _report_phases(phases, adata.n_obs,
                   "Phase distribution (ATAC, fragment-based):", n_missing)
    return adata


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Annotate an h5ad with cycling_score + CellCyclePhase obs columns."
    )
    parser.add_argument("--h5ad-input",  required=True, help="Input .h5ad path.")
    parser.add_argument("--h5ad-output", required=True,
                        help="Output .h5ad path. Can be same as --h5ad-input to overwrite.")
    parser.add_argument("--modality", choices=("rna", "atac"), default="rna",
                        help="rna: peak-independent UMI-based scoring. "
                             "atac: peak-accessibility over cell-cycle gene promoters.")
    parser.add_argument("--overwrite-ok", action="store_true",
                        help="Silently allow --h5ad-output to equal --h5ad-input.")

    # RNA-specific
    parser.add_argument("--sym-col", default="Gene",
                        help="[RNA] adata.var column with gene symbols. Falls back to var_names.")

    # ATAC-specific
    parser.add_argument("--gtf-path", default=None,
                        help="[ATAC] Path to GTF/GFF3 for TSS lookup. If omitted, "
                             "snapatac2's built-in annotations are used.")
    parser.add_argument("--genome", default="hg38",
                        help="[ATAC] Genome name for snapatac2 built-in TSS lookup.")
    parser.add_argument("--tss-window", type=int, default=2000,
                        help="[ATAC] ± bp around TSS to consider promoter peaks.")
    parser.add_argument("--fragments", default=None,
                        help="[ATAC] fragments.tsv(.bgz) file. If given, score cells by "
                             "fragments in promoter windows (snapatac2) instead of adata.X peaks.")
    parser.add_argument("--tempdir", default=None,
                        help="[ATAC --fragments] scratch dir for snapatac2 import (large!).")
    parser.add_argument("--n-jobs", type=int, default=8,
                        help="[ATAC --fragments] threads for snapatac2 import_data.")

    # ATAC: which signal to score. `promoter` is the gene-set route above;
    # `replication` ignores gene sets and uses copy-number effects instead.
    parser.add_argument("--atac-signal", choices=("promoter", "replication"),
                        default="promoter",
                        help="promoter: accessibility over cell-cycle gene promoters. "
                             "replication: S phase from megabase-scale coverage "
                             "overdispersion (replicating cells carry 2 copies of what "
                             "they have copied and 1 of the rest) and G2M from DNA "
                             "content. Requires --fragments.")
    parser.add_argument("--bin-size", type=int, default=None,
                        help="[--atac-signal replication] genome bin width in bp. "
                             "Default: chosen from the data's depth.")
    parser.add_argument("--min-bin-pct", type=float, default=None,
                        help="[--atac-signal replication] drop bins below this "
                             "percentile of population coverage (gaps, centromeres).")
    parser.add_argument("--cell-type-col", default=None,
                        help="[--atac-signal replication] obs column with cell-type "
                             "labels. Baselines become per-cell-type, and the run "
                             "prints a phase-by-cell-type control table.")
    parser.add_argument("--s-dispersion-threshold", type=float, default=None,
                        help="[--atac-signal replication] overdispersion above the "
                             "cell-type baseline needed for an S call.")
    parser.add_argument("--s-z-threshold", type=float, default=None,
                        help="[--atac-signal replication] how many sd of the Poisson "
                             "null that overdispersion must also clear.")
    parser.add_argument("--g2m-dna-threshold", type=float, default=None,
                        help="[--atac-signal replication] log2 DNA content above the "
                             "cell-type baseline needed for a G2M call.")
    parser.add_argument("--force-g2m", action="store_true",
                        help="[--atac-signal replication] make G2M calls even when "
                             "capture variation is too wide to resolve them.")

    # Scoring
    parser.add_argument("--scoring", choices=SCORING_MODES, default="fraction",
                        help="fraction: share of the cell's total signal in the gene set "
                             "(thresholds are absolute fractions, and depend on set size). "
                             "background: set mean minus signal-matched control mean, as a "
                             "z-score (thresholds are in SDs and comparable across sets).")
    parser.add_argument("--n-bins", type=int, default=BG_N_BINS,
                        help=f"[--scoring background] signal bins features are matched "
                             f"within. Default {BG_N_BINS}.")
    parser.add_argument("--n-ctrl-per-region", type=int, default=BG_N_CTRL_PER_REGION,
                        help=f"[--scoring background] controls drawn per set feature. "
                             f"Default {BG_N_CTRL_PER_REGION}.")
    parser.add_argument("--n-background-genes", type=int, default=BG_N_BACKGROUND_GENES,
                        help=f"[--scoring background, --fragments] non-cell-cycle genes whose "
                             f"promoters form the control pool. Default {BG_N_BACKGROUND_GENES}.")
    parser.add_argument("--seed", type=int, default=BG_SEED,
                        help=f"[--scoring background] seed for control sampling. Default {BG_SEED}.")

    # Thresholds — unset means "the default for --scoring". Under `fraction`
    # those defaults are RNA-calibrated and do not transfer to ATAC; check the
    # distribution table printed at run time and tune. Under `background` they
    # are in SDs, so 1.0 means "1 SD above this dataset's mean".
    parser.add_argument("--cycling-threshold", type=float, default=None,
                        help=f"cycling_score > this = cycling. Default "
                             f"{CYCLING_THRESHOLD} (fraction) / {BG_CYCLING_THRESHOLD} (background).")
    parser.add_argument("--g1-threshold",  type=float, default=None,
                        help=f"G1 cutoff. Default {G1_THRESHOLD} / {BG_G1_THRESHOLD}.")
    parser.add_argument("--s-threshold",   type=float, default=None,
                        help=f"S cutoff. Default {S_THRESHOLD} / {BG_S_THRESHOLD}.")
    parser.add_argument("--g2m-threshold", type=float, default=None,
                        help=f"G2M cutoff. Default {G2M_THRESHOLD} / {BG_G2M_THRESHOLD}.")

    args = parser.parse_args()

    use_replication = (args.modality == "atac" and args.atac_signal == "replication")
    if use_replication and not args.fragments:
        parser.error("--atac-signal replication requires --fragments: copy number "
                     "is read from genome-wide coverage, not from the peak matrix.")

    in_path  = Path(args.h5ad_input).resolve()
    out_path = Path(args.h5ad_output).resolve()

    if in_path == out_path and not args.overwrite_ok:
        print(f"[warn] output == input ({out_path}). Overwriting in place. "
              f"Pass --overwrite-ok to silence this warning.")

    print(f"Loading: {in_path}")
    import scanpy as sc
    adata = sc.read_h5ad(str(in_path))
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} features  (modality={args.modality})")

    common = dict(
        scoring=args.scoring,
        cycling_threshold=args.cycling_threshold,
        g1_threshold=args.g1_threshold,
        s_threshold=args.s_threshold,
        g2m_threshold=args.g2m_threshold,
        n_bins=args.n_bins,
        n_ctrl_per_region=args.n_ctrl_per_region,
        seed=args.seed,
    )

    if use_replication:
        from cell_cycle_atac_replication import (
            annotate_cell_cycle_atac_replication,
            DEFAULT_MIN_BIN_PCT, DEFAULT_S_THRESHOLD,
            DEFAULT_S_Z_THRESHOLD, DEFAULT_G2M_THRESHOLD,
        )
        annotate_cell_cycle_atac_replication(
            adata,
            fragments_path=args.fragments,
            genome=args.genome,
            bin_size=args.bin_size,
            min_bin_pct=(DEFAULT_MIN_BIN_PCT if args.min_bin_pct is None
                         else args.min_bin_pct),
            cell_type_col=args.cell_type_col,
            s_threshold=(DEFAULT_S_THRESHOLD if args.s_dispersion_threshold is None
                         else args.s_dispersion_threshold),
            s_z_threshold=(DEFAULT_S_Z_THRESHOLD if args.s_z_threshold is None
                           else args.s_z_threshold),
            g2m_threshold=(DEFAULT_G2M_THRESHOLD if args.g2m_dna_threshold is None
                           else args.g2m_dna_threshold),
            force_g2m=args.force_g2m,
            tempdir=args.tempdir,
            n_jobs=args.n_jobs,
        )
    elif args.modality == "rna":
        annotate_cell_cycle(adata, sym_col=args.sym_col, **common)
    elif args.fragments:
        annotate_cell_cycle_fragments(
            adata,
            fragments_path=args.fragments,
            gtf_path=args.gtf_path,
            genome=args.genome,
            tss_window=args.tss_window,
            n_background_genes=args.n_background_genes,
            tempdir=args.tempdir,
            n_jobs=args.n_jobs,
            **common,
        )
    else:
        annotate_cell_cycle_atac(
            adata,
            gtf_path=args.gtf_path,
            genome=args.genome,
            tss_window=args.tss_window,
            **common,
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nSaving: {out_path}")
    adata.write_h5ad(str(out_path))
    print("Done.")


if __name__ == "__main__":
    main()
