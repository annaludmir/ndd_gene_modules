"""
neuron_type_annotation.py

Annotate the neurons of an h5ad as excitatory or inhibitory from marker-gene
expression, and save the result as a new h5ad. Adds one obs column,
`NeuronType` (categorical), for cells whose `CellClass` is "Neuron"; every
other cell is left empty (NaN).

Per neuron, the UMIs of each marker set are summed. Each set has the
neurotransmitter machinery, plus lineage transcription factors that switch on
earlier, so immature neurons (common in fetal tissue) are not left unassigned:
  excitatory: SLC17A7 (VGLUT1), SLC17A6 (VGLUT2), SLC17A8 (VGLUT3);
              NEUROD2, NEUROD6, TBR1 (pallial glutamatergic);
              EOMES, NEUROG2, NEUROD1 (glutamatergic lineage commitment)
  inhibitory: GAD1, GAD2, SLC32A1 (VGAT), SLC6A5 (GlyT2, glycinergic);
              DLX1, DLX2, DLX5, DLX6, LHX6 (forebrain GABAergic),
              PAX2 (mid/hindbrain GABAergic);
              ASCL1, GSX2, NKX2-1, PTF1A (GABAergic lineage commitment)
The lineage-commitment genes are progenitor markers that newborn neurons
still carry. Only neurons are classified either way.

and the cell is called
  Excitatory  if only excitatory markers are detected, or excitatory UMIs are
              >= --min-ratio x inhibitory UMIs
  Inhibitory  the same, the other way round
  Mixed       both detected, neither dominating by --min-ratio
  Unassigned  neither detected (common in sparse scRNA-seq)

adata.X is taken as UMI counts, as in cell_cycle_annotation.py. Both sets'
UMIs come from the same cell, so library size cancels out of the comparison.

A second column, `NeuronMaturity`, calls the same neurons Immature (young),
Mature, Intermediate or Unassigned from two more marker sets:
  immature: DCX, SOX11, SOX4, DPYSL3, CD24, NNAT (newborn / migrating neurons)
  mature:   RBFOX3 (NeuN), SYP, SYT1, SNAP25, NEFL, NEFM (synaptic, neurofilament)
These genes differ widely in expression level (SNAP25 dwarfs most immature
markers), so raw sums would not compare. Per neuron, each marker is depth-
normalized (counts per 10k) and divided by its mean over all neurons, so every
gene weighs the same; a set's score is the mean over its genes. Then
  neuron_maturity_score = log2((mature + eps) / (immature + eps))
and Mature if score >= log2(--maturity-min-ratio), Immature if <= -that,
Intermediate in between, Unassigned if no marker of either set is detected.
Development is a continuum, so Intermediate is expected to be common.

Usage:
  python modules/neuron_type_annotation.py \
    --h5ad-input  data/human_dev_without_week_5.h5ad \
    --h5ad-output data/human_dev_without_week_5_neuron_type.h5ad
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from cell_cycle_annotation import _gene_indices, _symbol_to_var_map


EXCITATORY_MARKERS = [
    "SLC17A7", "SLC17A6", "SLC17A8",          # neurotransmitter: VGLUT1/2/3
    "NEUROD2", "NEUROD6", "TBR1",             # lineage TFs: pallial glutamatergic
    "EOMES", "NEUROG2", "NEUROD1",            # progenitor commitment: glutamatergic
]
INHIBITORY_MARKERS = [
    "GAD1", "GAD2", "SLC32A1", "SLC6A5",      # neurotransmitter: GAD67/65, VGAT, GlyT2
    "DLX1", "DLX2", "DLX5", "DLX6", "LHX6",   # lineage TFs: forebrain GABAergic
    "PAX2",                                   # lineage TF: mid/hindbrain GABAergic
    "ASCL1", "GSX2", "NKX2-1", "PTF1A",       # progenitor commitment: GABAergic
]

IMMATURE_MARKERS = ["DCX", "SOX11", "SOX4", "DPYSL3", "CD24", "NNAT"]
MATURE_MARKERS   = ["RBFOX3", "SYP", "SYT1", "SNAP25", "NEFL", "NEFM"]

MATURITY_COL           = "NeuronMaturity"
MATURITY_SCORE_COL     = "neuron_maturity_score"
MATURITY_ORDER         = ["Immature", "Intermediate", "Mature", "Unassigned"]
DEFAULT_MATURITY_RATIO = 1.5

NEURON_TYPE_COL   = "NeuronType"
NEURON_TYPE_ORDER = ["Excitatory", "Inhibitory", "Mixed", "Unassigned"]

DEFAULT_CELL_CLASS_COL = "CellClass"
DEFAULT_NEURON_CLASS   = "Neuron"
DEFAULT_MIN_RATIO      = 2.0


def _set_umis(X, rows: np.ndarray, gene_ix: np.ndarray) -> np.ndarray:
    """Per-cell sum of X[rows, gene_ix] (X sparse or dense)."""
    return np.asarray(X[rows][:, gene_ix].sum(axis=1)).ravel()


def _neuron_rows(adata, cell_class_col: str, neuron_class: str) -> np.ndarray:
    """Row indices of cells with obs[cell_class_col] == neuron_class."""
    if cell_class_col not in adata.obs.columns:
        raise KeyError(f"obs column '{cell_class_col}' not found; set --cell-class-col. "
                       f"Available: {list(adata.obs.columns)}")
    rows = np.flatnonzero((adata.obs[cell_class_col].astype(str) == str(neuron_class)).to_numpy())
    print(f"  {len(rows):,}/{adata.n_obs:,} cells with {cell_class_col} == '{neuron_class}'")
    if len(rows) == 0:
        raise ValueError(f"No cells with {cell_class_col} == '{neuron_class}'. Values: "
                         f"{sorted(adata.obs[cell_class_col].astype(str).unique())[:20]}")
    return rows


def _marker_sets(adata, sym_col: str, sets: dict[str, list[str]]) -> dict[str, np.ndarray]:
    """{name: var indices} for each marker set; an entirely missing set is an error."""
    print(f"  Resolving marker genes (sym_col='{sym_col}'):")
    sym2var = _symbol_to_var_map(adata, sym_col)
    var_pos = {v: i for i, v in enumerate(adata.var_names.astype(str))}
    out = {name: _gene_indices(genes, sym2var, var_pos, name) for name, genes in sets.items()}
    if any(len(ix) == 0 for ix in out.values()):
        raise RuntimeError("An entire marker set is missing from adata.var; check --sym-col "
                           "and the gene naming convention in the h5ad.")
    return out


def _report(labels, order, col: str) -> None:
    counts = pd.Series(labels).value_counts()
    print(f"\n  {col} among {len(labels):,} neurons:")
    for label in order:
        n = int(counts.get(label, 0))
        print(f"    {label:12s}  {n:>10,}  ({100.0 * n / len(labels):5.2f}%)")


def classify_neurons(exc: np.ndarray, inh: np.ndarray, min_ratio: float) -> np.ndarray:
    """Excitatory / Inhibitory / Mixed / Unassigned from per-cell marker UMIs."""
    labels = np.full(len(exc), "Unassigned", dtype=object)
    both = (exc > 0) & (inh > 0)
    labels[(exc > 0) & (inh == 0)] = "Excitatory"
    labels[(inh > 0) & (exc == 0)] = "Inhibitory"
    labels[both] = "Mixed"
    labels[both & (exc >= min_ratio * inh)] = "Excitatory"
    labels[both & (inh >= min_ratio * exc)] = "Inhibitory"
    return labels


def annotate_neuron_type(
    adata,
    sym_col: str = "Gene",
    cell_class_col: str = DEFAULT_CELL_CLASS_COL,
    neuron_class: str = DEFAULT_NEURON_CLASS,
    min_ratio: float = DEFAULT_MIN_RATIO,
):
    """Add `NeuronType` to adata.obs (in place) for cells whose
    obs[cell_class_col] == neuron_class; other cells get NaN."""
    rows = _neuron_rows(adata, cell_class_col, neuron_class)
    ix = _marker_sets(adata, sym_col, {"excitatory": EXCITATORY_MARKERS,
                                       "inhibitory": INHIBITORY_MARKERS})
    exc_ix, inh_ix = ix["excitatory"], ix["inhibitory"]

    import scipy.sparse as sp
    X = adata.X.tocsr() if sp.issparse(adata.X) else adata.X
    exc = _set_umis(X, rows, exc_ix)
    inh = _set_umis(X, rows, inh_ix)
    labels = classify_neurons(exc, inh, min_ratio)

    neuron_type = np.full(adata.n_obs, None, dtype=object)
    neuron_type[rows] = labels
    adata.obs[NEURON_TYPE_COL] = pd.Categorical(neuron_type, categories=NEURON_TYPE_ORDER)
    adata.uns["neuron_type_annotation"] = {
        "excitatory_markers": [str(g) for g in EXCITATORY_MARKERS],
        "inhibitory_markers": [str(g) for g in INHIBITORY_MARKERS],
        "cell_class_col": cell_class_col,
        "neuron_class": str(neuron_class),
        "min_ratio": float(min_ratio),
    }

    _report(labels, NEURON_TYPE_ORDER, NEURON_TYPE_COL)
    return adata


def annotate_neuron_maturity(
    adata,
    sym_col: str = "Gene",
    cell_class_col: str = DEFAULT_CELL_CLASS_COL,
    neuron_class: str = DEFAULT_NEURON_CLASS,
    min_ratio: float = DEFAULT_MATURITY_RATIO,
):
    """Add `NeuronMaturity` and `neuron_maturity_score` to adata.obs (in
    place) for neurons; other cells get NaN. See the module docstring."""
    import scipy.sparse as sp
    rows = _neuron_rows(adata, cell_class_col, neuron_class)
    ix = _marker_sets(adata, sym_col, {"immature": IMMATURE_MARKERS, "mature": MATURE_MARKERS})

    X = adata.X.tocsr() if sp.issparse(adata.X) else adata.X
    total = np.asarray(X[rows].sum(axis=1)).ravel()
    scale = 1e4 / np.where(total > 0, total, 1.0)

    def set_score(gene_ix):
        sub = X[rows][:, gene_ix]
        sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
        cp10k = sub * scale[:, None]
        gene_mean = cp10k.mean(axis=0)
        keep = gene_mean > 0                      # genes never seen in neurons carry no weight
        rel = cp10k[:, keep] / gene_mean[keep]
        detected = (sub[:, keep] > 0).any(axis=1)
        return (rel.mean(axis=1) if keep.any() else np.zeros(len(rows))), detected

    imm, imm_det = set_score(ix["immature"])
    mat, mat_det = set_score(ix["mature"])
    eps = 1e-3
    score = np.log2((mat + eps) / (imm + eps))
    cut = np.log2(min_ratio)

    labels = np.full(len(rows), "Intermediate", dtype=object)
    labels[score >= cut]  = "Mature"
    labels[score <= -cut] = "Immature"
    none = ~imm_det & ~mat_det
    labels[none] = "Unassigned"
    score[none] = np.nan

    col = np.full(adata.n_obs, None, dtype=object); col[rows] = labels
    adata.obs[MATURITY_COL] = pd.Categorical(col, categories=MATURITY_ORDER)
    sc_all = np.full(adata.n_obs, np.nan, dtype=np.float32); sc_all[rows] = score
    adata.obs[MATURITY_SCORE_COL] = sc_all
    adata.uns["neuron_maturity_annotation"] = {
        "immature_markers": [str(g) for g in IMMATURE_MARKERS],
        "mature_markers":   [str(g) for g in MATURE_MARKERS],
        "cell_class_col": cell_class_col,
        "neuron_class": str(neuron_class),
        "min_ratio": float(min_ratio),
    }

    q = np.nanpercentile(score, [10, 25, 50, 75, 90])
    print(f"  {MATURITY_SCORE_COL} (log2 mature/immature) p10/25/50/75/90: "
          + " / ".join(f"{v:.2f}" for v in q) + f"   cutoffs: ±{cut:.2f}")
    _report(labels, MATURITY_ORDER, MATURITY_COL)
    return adata


def main():
    p = argparse.ArgumentParser(
        description="Add excitatory / inhibitory NeuronType and Immature / Mature "
                    "NeuronMaturity obs columns to an h5ad.")
    p.add_argument("--h5ad-input", required=True, help="Input .h5ad path.")
    p.add_argument("--h5ad-output", default=None,
                   help="Output .h5ad path. Default: <input>_neuron_type.h5ad next to the input.")
    p.add_argument("--sym-col", default="Gene",
                   help="adata.var column with gene symbols. Falls back to var_names.")
    p.add_argument("--cell-class-col", default=DEFAULT_CELL_CLASS_COL,
                   help=f"obs column identifying neurons. Default '{DEFAULT_CELL_CLASS_COL}'.")
    p.add_argument("--neuron-class", default=DEFAULT_NEURON_CLASS,
                   help=f"Value of --cell-class-col marking neurons. Default '{DEFAULT_NEURON_CLASS}'.")
    p.add_argument("--min-ratio", type=float, default=DEFAULT_MIN_RATIO,
                   help=f"When both marker sets are detected, one must have at least this "
                        f"many times the other's UMIs to win; otherwise 'Mixed'. "
                        f"Default {DEFAULT_MIN_RATIO}.")
    p.add_argument("--maturity-min-ratio", type=float, default=DEFAULT_MATURITY_RATIO,
                   help=f"Mature / Immature when one set's normalized score is at least this "
                        f"many times the other's; otherwise 'Intermediate'. "
                        f"Default {DEFAULT_MATURITY_RATIO}.")
    args = p.parse_args()

    in_path  = Path(args.h5ad_input).resolve()
    out_path = (Path(args.h5ad_output).resolve() if args.h5ad_output
                else in_path.with_name(f"{in_path.stem}_neuron_type.h5ad"))
    if out_path == in_path:
        p.error("--h5ad-output must differ from --h5ad-input; the result is saved as a new file.")

    import scanpy as sc
    print(f"Loading: {in_path}")
    adata = sc.read_h5ad(str(in_path))
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    annotate_neuron_type(adata, sym_col=args.sym_col, cell_class_col=args.cell_class_col,
                         neuron_class=args.neuron_class, min_ratio=args.min_ratio)
    annotate_neuron_maturity(adata, sym_col=args.sym_col, cell_class_col=args.cell_class_col,
                             neuron_class=args.neuron_class,
                             min_ratio=args.maturity_min_ratio)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nSaving: {out_path}")
    adata.write_h5ad(str(out_path))
    print("Done.")


if __name__ == "__main__":
    main()
