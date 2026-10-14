"""
neuron_type_annotation.py

Annotate the neurons of an h5ad as excitatory or inhibitory from marker-gene
expression, and save the result as a new h5ad. Adds one obs column,
`NeuronType` (categorical), for cells whose `CellClass` is "Neuron"; every
other cell is left empty (NaN).

Per neuron, the UMIs of each marker set are summed:
  excitatory: SLC17A7 (VGLUT1), SLC17A6 (VGLUT2), SLC17A8 (VGLUT3)
  inhibitory: GAD1, GAD2, SLC32A1 (VGAT)

and the cell is called
  Excitatory  if only excitatory markers are detected, or excitatory UMIs are
              >= --min-ratio x inhibitory UMIs
  Inhibitory  the same, the other way round
  Mixed       both detected, neither dominating by --min-ratio
  Unassigned  neither detected (common in sparse scRNA-seq)

adata.X is taken as UMI counts, as in cell_cycle_annotation.py. Both sets'
UMIs come from the same cell, so library size cancels out of the comparison.

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


EXCITATORY_MARKERS = ["SLC17A7", "SLC17A6", "SLC17A8"]
INHIBITORY_MARKERS = ["GAD1", "GAD2", "SLC32A1"]

NEURON_TYPE_COL   = "NeuronType"
NEURON_TYPE_ORDER = ["Excitatory", "Inhibitory", "Mixed", "Unassigned"]

DEFAULT_CELL_CLASS_COL = "CellClass"
DEFAULT_NEURON_CLASS   = "Neuron"
DEFAULT_MIN_RATIO      = 2.0


def _set_umis(X, rows: np.ndarray, gene_ix: np.ndarray) -> np.ndarray:
    """Per-cell sum of X[rows, gene_ix] (X sparse or dense)."""
    return np.asarray(X[rows][:, gene_ix].sum(axis=1)).ravel()


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
    if cell_class_col not in adata.obs.columns:
        raise KeyError(f"obs column '{cell_class_col}' not found; set --cell-class-col. "
                       f"Available: {list(adata.obs.columns)}")
    is_neuron = (adata.obs[cell_class_col].astype(str) == str(neuron_class)).to_numpy()
    rows = np.flatnonzero(is_neuron)
    print(f"  {len(rows):,}/{adata.n_obs:,} cells with {cell_class_col} == '{neuron_class}'")
    if len(rows) == 0:
        raise ValueError(f"No cells with {cell_class_col} == '{neuron_class}'. Values: "
                         f"{sorted(adata.obs[cell_class_col].astype(str).unique())[:20]}")

    print(f"  Resolving marker genes (sym_col='{sym_col}'):")
    sym2var = _symbol_to_var_map(adata, sym_col)
    var_pos = {v: i for i, v in enumerate(adata.var_names.astype(str))}
    exc_ix  = _gene_indices(EXCITATORY_MARKERS, sym2var, var_pos, "excitatory")
    inh_ix  = _gene_indices(INHIBITORY_MARKERS, sym2var, var_pos, "inhibitory")
    if len(exc_ix) == 0 or len(inh_ix) == 0:
        raise RuntimeError("An entire marker set is missing from adata.var; check --sym-col "
                           "and the gene naming convention in the h5ad.")

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

    counts = pd.Series(labels).value_counts()
    print(f"\n  {NEURON_TYPE_COL} among {len(rows):,} neurons:")
    for label in NEURON_TYPE_ORDER:
        n = int(counts.get(label, 0))
        print(f"    {label:11s}  {n:>10,}  ({100.0 * n / len(rows):5.2f}%)")
    return adata


def main():
    p = argparse.ArgumentParser(
        description="Add an excitatory / inhibitory NeuronType obs column to an h5ad.")
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

    out_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nSaving: {out_path}")
    adata.write_h5ad(str(out_path))
    print("Done.")


if __name__ == "__main__":
    main()
