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

# Thresholds (fractions, not percentages).
CYCLING_THRESHOLD = 0.004     # 0.4%
G1_THRESHOLD      = 0.002   # 0.2%
S_THRESHOLD       = 0.002   # 0.2%
G2M_THRESHOLD     = 0.03    # 3%

PHASE_ORDER = ["Non-cycling", "G1", "S", "G2M", "Post-M"]


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
# Main annotation routine
# ---------------------------------------------------------------------------

def annotate_cell_cycle(adata, sym_col: str = "Gene"):
    """Add `cycling_score` and `CellCyclePhase` columns to adata.obs (in place)."""
    import scipy.sparse as sp

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

    print("  Computing UMI totals and per-cell fractions ...")
    total_umi     = np.asarray(X.sum(axis=1)).ravel()
    cycling_score = _per_cell_fraction(X, cc_ix,  total_umi)
    g1_frac       = _per_cell_fraction(X, g1_ix,  total_umi)
    s_frac        = _per_cell_fraction(X, s_ix,   total_umi)
    g2m_frac      = _per_cell_fraction(X, g2m_ix, total_umi)

    # Classification (priority: G2M > S > G1 > Post-M).
    # A cell that passes multiple phase thresholds is assigned the latest phase
    # in the cycle order, since later markers are more informative.
    is_cycling = cycling_score > CYCLING_THRESHOLD

    passes_g1  = g1_frac  > G1_THRESHOLD
    passes_s   = s_frac   > S_THRESHOLD
    passes_g2m = g2m_frac > G2M_THRESHOLD

    phases = np.full(adata.n_obs, "Non-cycling", dtype=object)
    # Cycling cells default to Post-M; then overwritten by G1 → S → G2M so
    # G2M ends up with the highest priority.
    phases[is_cycling] = "Post-M"
    phases[is_cycling & passes_g1]  = "G1"
    phases[is_cycling & passes_s]   = "S"
    phases[is_cycling & passes_g2m] = "G2M"

    adata.obs["cycling_score"]  = cycling_score.astype(np.float32)
    adata.obs["CellCyclePhase"] = pd.Categorical(
        phases, categories=PHASE_ORDER, ordered=False,
    )

    # Log summary.
    counts = pd.Series(phases).value_counts()
    total  = adata.n_obs
    print("\n  Phase distribution:")
    for phase in PHASE_ORDER:
        n = int(counts.get(phase, 0))
        pct = 100.0 * n / total if total else 0.0
        print(f"    {phase:12s}  {n:>10,}  ({pct:5.2f}%)")
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
    parser.add_argument("--sym-col", default="Gene",
                        help="adata.var column with gene symbols. Falls back to var_names.")
    parser.add_argument("--overwrite-ok", action="store_true",
                        help="Silently allow --h5ad-output to equal --h5ad-input.")
    args = parser.parse_args()

    in_path  = Path(args.h5ad_input).resolve()
    out_path = Path(args.h5ad_output).resolve()

    if in_path == out_path and not args.overwrite_ok:
        print(f"[warn] output == input ({out_path}). Overwriting in place. "
              f"Pass --overwrite-ok to silence this warning.")

    print(f"Loading: {in_path}")
    import scanpy as sc
    adata = sc.read_h5ad(str(in_path))
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    annotate_cell_cycle(adata, sym_col=args.sym_col)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nSaving: {out_path}")
    adata.write_h5ad(str(out_path))
    print("Done.")


if __name__ == "__main__":
    main()
