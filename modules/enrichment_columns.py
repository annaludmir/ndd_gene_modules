"""
enrichment_columns.py

Read entries of an enrichment config's `column_conditions_for_gsea`. A column is
either a plain list of conditions, scored against the config's top-level
`ges_results_folder`:

  CellClass:
    - "Radial glia"
    - "Neuron"

or a block with its own GES run, for columns scored in a separate GES run:

  NeuronType:
    ges_results_folder: "results/ges_score_results/ges_score_for_cortex_neuron_type_20261026"
    conditions:
      - "Excitatory"
      - "Inhibitory"
"""

from pathlib import Path


def column_conditions(spec) -> list:
    """The condition list of one column entry (list, or block with `conditions`)."""
    if isinstance(spec, dict):
        if "conditions" not in spec:
            raise ValueError(f"column block needs a 'conditions' list: {spec}")
        return list(spec["conditions"])
    return list(spec or [])


def column_ges_folder(spec, default):
    """The GES run folder for one column entry: its own `ges_results_folder`
    if it is a block that sets one, else `default` (the config's top-level one)."""
    if isinstance(spec, dict) and spec.get("ges_results_folder"):
        return spec["ges_results_folder"]
    return default


def resolve_column_ges_folders(column_spec: dict, root: Path) -> None:
    """Resolve per-column `ges_results_folder` paths against `root` in place and
    fail early, naming the column, if one does not exist."""
    for column, spec in column_spec.items():
        if not (isinstance(spec, dict) and spec.get("ges_results_folder")):
            continue
        p = Path(spec["ges_results_folder"])
        p = p if p.is_absolute() else Path(root) / p
        if not p.is_dir():
            raise FileNotFoundError(
                f"ges_results_folder for column '{column}' not found: {p} — set it to "
                f"an existing (dated) GES run folder under results/ges_score_results/.")
        spec["ges_results_folder"] = p
