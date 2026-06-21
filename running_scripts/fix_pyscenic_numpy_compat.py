#!/usr/bin/env python3
"""
Patch pyscenic source files for NumPy >= 1.24 and Dask compatibility.

NumPy 1.24 removed the deprecated aliases np.object, np.bool, np.int,
np.float, np.complex, np.str.  pyscenic's transform.py uses np.object at
module level, which crashes every Dask worker process that imports it.

Run ONCE on the cluster (no arguments needed):
  python running_scripts/fix_pyscenic_numpy_compat.py

Backups are written next to each patched file as <name>.py.bak
"""

import re
import sys
from pathlib import Path

CONDA_ENV = Path("/miridan-data/annaludmir/conda-envs/jupyter-scanpy_new")
SITE_PKGS  = CONDA_ENV / "lib/python3.10/site-packages"
PYSCENIC_DIR = SITE_PKGS / "pyscenic"

# ---------------------------------------------------------------------------
# Fix 1 — np.<alias> → <alias>  (NumPy removed these in 1.24)
# Regex uses \b word-boundary so np.object_ / np.float64 are left untouched.
# ---------------------------------------------------------------------------

NP_ALIASES = ["object", "bool", "int", "float", "complex", "str"]

def _patch_numpy_aliases(text: str) -> tuple[str, list[str]]:
    changes = []
    for alias in NP_ALIASES:
        pattern = rf"\bnp\.{alias}\b"
        new_text, n = re.subn(pattern, alias, text)
        if n:
            changes.append(f"  np.{alias} → {alias}  ({n} occurrence{'s' if n > 1 else ''})")
            text = new_text
    return text, changes


# ---------------------------------------------------------------------------
# Fix 2 — from_delayed(generator, ...) → from_delayed(list(generator), ...)
# Newer dask.dataframe.from_delayed requires a sequence, not a generator.
# ---------------------------------------------------------------------------

def _patch_from_delayed_generator(text: str) -> tuple[str, list[str]]:
    """
    Wraps bare generator expressions passed as the first argument of
    from_delayed() in list(…).  Handles the single-line form:

        from_delayed(
            (expr for x in seq),
    """
    changes = []

    # Pattern: from_delayed(\n    (  — a generator on the very next line
    pattern = r"(from_delayed\(\s*\n\s*)(\()"
    def wrap_gen(m):
        return m.group(1) + "list(" + m.group(2)

    new_text, n = re.subn(pattern, wrap_gen, text)
    if n:
        # We opened list(  — now find the matching closing ) of the generator
        # and insert a closing ) after it.  Track depth to find the close paren.
        result = []
        i = 0
        depth_list = 0          # tracks open list( insertions we need to close
        inside_list_wrap = False
        while i < len(new_text):
            ch = new_text[i]
            # Detect our inserted "list(" marker
            if new_text[i:i+5] == "list(" and not inside_list_wrap:
                inside_list_wrap = True
                depth_list = 0
                result.append(new_text[i:i+5])
                i += 5
                continue
            if inside_list_wrap:
                if ch == "(":
                    depth_list += 1
                elif ch == ")":
                    if depth_list == 0:
                        # This ) closes the generator expression — add our )
                        result.append(ch)   # close the generator (
                        result.append(")")  # close our list(
                        inside_list_wrap = False
                        i += 1
                        continue
                    depth_list -= 1
            result.append(ch)
            i += 1
        new_text = "".join(result)
        changes.append(f"  from_delayed(generator) → from_delayed(list(generator))  ({n} site{'s' if n > 1 else ''})")

    return new_text, changes


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def patch_file(path: Path, fixers) -> bool:
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
    for py_file in sorted(PYSCENIC_DIR.glob("*.py")):
        if py_file.stem.endswith(".bak"):
            continue
        changed = patch_file(
            py_file,
            fixers=[_patch_numpy_aliases, _patch_from_delayed_generator],
        )
        if changed:
            patched += 1

    if patched == 0:
        print("No changes needed — pyscenic files are already compatible.")
    else:
        print(f"\nDone — {patched} file(s) patched.")
        print("Re-run tf_network.py normally; the Dask workers will now import the fixed files.")


if __name__ == "__main__":
    main()
