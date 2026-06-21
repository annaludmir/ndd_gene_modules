#!/usr/bin/env python3
"""
Patch pyscenic + dask source files for NumPy >= 1.24 and Dask compatibility.

Fix 1 — NumPy aliases (np.object, np.bool, …) removed in 1.24:
  pyscenic/transform.py uses np.object at module level.  Because pyscenic
  uses Dask workers (subprocesses), every worker re-imports transform.py
  fresh — so the fix must live in the file, not in a main-process monkey-patch.

Fix 2 — dask from_delayed rejects generators:
  pyscenic/prune.py builds a Dask graph via an internal alias of
  dask.dataframe.from_delayed, passing a generator expression.  Newer Dask
  calls len() on the first argument and raises TypeError.  We fix this in
  dask's own _delayed.py (the exact line that fails) rather than trying to
  parse pyscenic's dynamic alias.

Run ONCE on the cluster (no arguments needed):
  python running_scripts/fix_pyscenic_numpy_compat.py

Backups are written next to each patched file as <name>.py.bak
"""

import re
import sys
from pathlib import Path

CONDA_ENV  = Path("/miridan-data/annaludmir/conda-envs/jupyter-scanpy_new")
SITE_PKGS  = CONDA_ENV / "lib/python3.10/site-packages"
PYSCENIC_DIR = SITE_PKGS / "pyscenic"

# Dask file where the TypeError actually originates
DASK_DELAYED_FILE = SITE_PKGS / "dask" / "dataframe" / "dask_expr" / "io" / "_delayed.py"


# ---------------------------------------------------------------------------
# Fix 1 — np.<alias> → <alias>  (NumPy removed these in 1.24)
# \b word-boundary keeps np.object_ / np.float64 / np.str_ untouched.
# ---------------------------------------------------------------------------

NP_ALIASES = ["object", "bool", "int", "float", "complex", "str"]

def _patch_numpy_aliases(text: str) -> tuple[str, list[str]]:
    changes = []
    for alias in NP_ALIASES:
        new_text, n = re.subn(rf"\bnp\.{alias}\b", alias, text)
        if n:
            changes.append(f"  np.{alias} → {alias}  ({n} occurrence{'s' if n > 1 else ''})")
            text = new_text
    return text, changes


# ---------------------------------------------------------------------------
# Fix 2 — dask from_delayed: accept generators by converting to list first.
# Target: dask/dataframe/dask_expr/io/_delayed.py  line ~122
#
#   BEFORE:
#       if len(dfs) == 0:
#
#   AFTER:
#       if not hasattr(dfs, "__len__"):
#           dfs = list(dfs)
#       if len(dfs) == 0:
# ---------------------------------------------------------------------------

_DASK_OLD = "    if len(dfs) == 0:"
_DASK_NEW  = (
    '    if not hasattr(dfs, "__len__"):\n'
    "        dfs = list(dfs)\n"
    "    if len(dfs) == 0:"
)

def _patch_dask_from_delayed(text: str) -> tuple[str, list[str]]:
    if _DASK_OLD not in text:
        return text, []
    if _DASK_NEW in text:          # already patched
        return text, []
    return text.replace(_DASK_OLD, _DASK_NEW, 1), [
        "  from_delayed: added generator→list coercion before len() check"
    ]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def patch_file(path: Path, fixers) -> bool:
    if not path.exists():
        print(f"  [skip] not found: {path}")
        return False
    original = path.read_text(encoding="utf-8")
    text = original
    all_changes = []
    for fixer in fixers:
        text, changes = fixer(text)
        all_changes.extend(changes)
    if text == original:
        return False
    backup = path.with_suffix(".py.bak")
    backup.write_text(original, encoding="utf-8")
    path.write_text(text, encoding="utf-8")
    print(f"\nPatched: {path}")
    for c in all_changes:
        print(c)
    print(f"  Backup: {backup.name}")
    return True


def main():
    if not PYSCENIC_DIR.exists():
        sys.exit(f"pyscenic not found at: {PYSCENIC_DIR}")

    patched = 0

    # Fix 1: numpy aliases in all pyscenic .py files
    for py_file in sorted(PYSCENIC_DIR.glob("*.py")):
        if py_file.suffix == ".bak":
            continue
        if patch_file(py_file, fixers=[_patch_numpy_aliases]):
            patched += 1

    # Fix 2: dask from_delayed generator handling
    if patch_file(DASK_DELAYED_FILE, fixers=[_patch_dask_from_delayed]):
        patched += 1

    if patched == 0:
        print("No changes needed — files are already compatible.")
    else:
        print(f"\nDone — {patched} file(s) patched.")
        print("Re-run tf_network.py; both fixes are now permanent in the conda env.")


if __name__ == "__main__":
    main()
