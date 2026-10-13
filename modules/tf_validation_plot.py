"""
tf_validation_plot.py

Network figures for tf_validation_report's annotated CSVs. Same layouts and
TF selection as tf_network's figures, and the target node colour is still the
normalized prediction score, but every TF -> target edge is drawn as a bundle
of parallel lines, one per cell type in which ATAC-seq validates that pair:

  - one solid line per validating cell type of the chosen cell class
  - not validated in any of them -> a single thin grey dotted line

By default only validated nodes are shown: TFs without a single edge
validated in the drawn cell types are dropped before the top-N TFs are picked,
and targets without a validated edge to a shown TF are dropped. Grey edges
between the remaining nodes stay. --show-unvalidated-nodes keeps everything.

A gene list with more than MAX_TFS_PER_FIGURE (5) TFs to show is split, in
TF rank order, into the fewest figures of 2-3 TFs each (sizes as even as
possible), saved as `..._part<i>of<n>.png`; within each figure, targets need
a validated edge to one of that figure's TFs.

Only the cell types of one class from --cell-class-json are drawn (default:
progenitor_and_stem_cells from config_files/atac_cell_classes.json). The
class's cell types are matched to the CSV's `validated_<cell_type>` columns by
the same name sanitizing tf_validation_report uses, and each keeps the colour
of its position in the JSON list, so colours are stable across figures even
when a cell type is missing from a run. Cell-phase columns
(`validated_phase_*`) are ignored.

With --motif-hits-csv (the ATAC run's global
`3_motif_target_validation/motif_target_pair_scores.csv`), grey edges whose
pair has a TF motif hit in an accessible peak of the target's promoter when
all peaks are considered (validated_by_motif_and_accessibility) get a star at
the edge midpoint: accessible, but not among the top peaks of any drawn cell
type. Coloured edges aren't starred - the per-cell-type peaks are a subset of
all peaks, so they are always globally validated too.

With --literature-dir, pairs are also tagged from the curated literature
CSVs in that folder (column "Evidence assessment"), with a diamond 3/4 of the
way to the target:
  - "Evidence found (direct)"            -> filled diamond
  - "Plausible, no pair-specific ..."    -> hollow diamond
Other assessments are left untagged. Evidence belongs to the (TF, target)
pair, not to a gene list, so every CSV in the folder is pooled and a pair is
tagged in every figure it appears in, whatever the CSV is named. If two CSVs
disagree, direct wins over plausible.

Inputs:
  --validation-dir <path>
     A tf_validation subfolder holding `tf_targets_<gene_list>.csv` files, or
     the results/tf_validation root, in which case every subfolder is plotted.

  --cell-class-json <path> / --cell-class <name>
     {class: [ATAC cell type, ...]} and the class to draw.

  --motif-hits-csv <path>   (optional)
     motif_target_pair_scores.csv from the matching atac_seq_analysis run.

  --literature-dir <path>   (optional)
     Folder of literature CSVs (tf, target, "Evidence assessment"), e.g.
     results/tf_validation/literature_evidence.

Output:
  <subfolder>/figures/tf_validation_network_<gene_list>.png, or
  <subfolder>/figures/tf_validation_network_<gene_list>_part<i>of<n>.png

Usage:
  python modules/tf_validation_plot.py \
    --validation-dir results/tf_validation/cortex_v3_20260802_vs_atac_20260816 \
    --plot-style connected --top-n-tfs 8 \
    --cell-class progenitor_and_stem_cells \
    --motif-hits-csv results/atac_analysis/<run>/3_motif_target_validation/motif_target_pair_scores.csv
"""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from tf_network import _TF_CMAP
from tf_validation_report import (
    AXES, TARGET_KEY_CANDIDATES, TF_KEY_CANDIDATES, VALIDATED_COL, _pick_col, _sanitize,
)


DEFAULT_CELL_CLASS_JSON = (Path(__file__).resolve().parent.parent
                           / "config_files" / "atac_cell_classes.json")
DEFAULT_CELL_CLASS = "progenitor_and_stem_cells"

# Categorical palette, assigned by position in the class's JSON list; classes
# longer than 8 fall back to the light tab20 variants.
CELL_TYPE_COLORS = [
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
    "#e87ba4", "#008300", "#4a3aa7", "#e34948",
]
CELL_TYPE_FALLBACK = [matplotlib.colormaps["tab20"](i) for i in range(1, 20, 2)]

UNVALIDATED_STYLE = dict(color="#bbbbbb", lw=0.7, linestyle=(0, (1, 2)))
EDGE_LW = 1.1
MOTIF_HIT_STYLE = dict(marker="*", markersize=9, color="#111111",
                       markeredgecolor="white", markeredgewidth=0.6, linestyle="none")
MOTIF_HIT_LABEL = "not validated, accessible motif hit (all peaks)"

LITERATURE_COL    = "Evidence assessment"
_LIT_MARKER = dict(marker="D", markersize=6, markeredgecolor="#111111",
                   markeredgewidth=1.0, linestyle="none")
LITERATURE_STYLES = {
    "direct":    dict(_LIT_MARKER, color="#111111"),
    "plausible": dict(_LIT_MARKER, color="white"),
}
LITERATURE_LABELS = {
    "direct":    "literature: direct evidence",
    "plausible": "literature: plausible",
}


# ---------------------------------------------------------------------------
# Cell class and the annotated CSV
# ---------------------------------------------------------------------------

def load_cell_class(json_path: str | Path, cell_class: str) -> list[str]:
    """The ATAC cell types of `cell_class` in a {class: [cell type, ...]} JSON."""
    with open(json_path) as f:
        classes = json.load(f)
    if cell_class not in classes:
        raise KeyError(f"cell class '{cell_class}' not in {json_path}; "
                       f"available: {list(classes)}")
    return [str(c) for c in classes[cell_class]]


def cell_type_columns(df: pd.DataFrame) -> list[str]:
    """The CSV's per-cell-type `validated_<slug>` columns (phase and summary
    columns excluded)."""
    reserved = {spec[k] for spec in AXES.values()
                for k in ("any_col", "count_col", "list_col")}
    phase_prefix = AXES["cell_phase"]["col_prefix"]
    ct_prefix    = AXES["cell_type"]["col_prefix"]
    return [c for c in df.columns
            if c.startswith(ct_prefix) and not c.startswith(phase_prefix)
            and c not in reserved]


def detect_groups(df: pd.DataFrame, cell_types: list[str]) -> dict[str, str]:
    """{cell type: CSV column} for the class's cell types present in the CSV,
    in class order."""
    present = set(cell_type_columns(df))
    prefix  = AXES["cell_type"]["col_prefix"]
    cols    = {ct: f"{prefix}{_sanitize(ct)}" for ct in cell_types}
    return {ct: col for ct, col in cols.items() if col in present}


def build_line_specs(cell_types: list[str]) -> dict[str, dict]:
    """{cell type: line style}, coloured by position in the class list."""
    palette = CELL_TYPE_COLORS + CELL_TYPE_FALLBACK
    if len(cell_types) > len(palette):
        raise ValueError(f"{len(cell_types)} cell types but only "
                         f"{len(palette)} colours available.")
    return {ct: dict(color=c, linestyle="-") for ct, c in zip(cell_types, palette)}


def edge_groups(df: pd.DataFrame, groups: dict[str, str]) -> list[list[str]]:
    """Per row, the cell types (in class order) that validate the pair."""
    if not groups:
        return [[] for _ in range(len(df))]
    names  = list(groups)
    matrix = df[list(groups.values())].fillna(False).astype(bool).to_numpy()
    return [[names[i] for i in np.where(row)[0]] for row in matrix]


def load_accessible_motif_pairs(path: str | Path) -> set[tuple[str, str]]:
    """(TF, target) pairs, upper-cased, with validated_by_motif_and_accessibility
    = True in an ATAC run's global motif_target_pair_scores.csv, i.e. >=1
    motif hit inside an accessible peak (all peaks). Pairs that were never
    scanned (TF without a motif, promoter not found) are absent."""
    df = pd.read_csv(path)
    tf_col     = _pick_col(df, TF_KEY_CANDIDATES)
    target_col = _pick_col(df, TARGET_KEY_CANDIDATES)
    hit = df[df[VALIDATED_COL].fillna(False).astype(bool)]
    print(f"Accessible motif hits: {len(hit):,} / {len(df):,} scanned pairs ({path})")
    return set(zip(hit[tf_col].astype(str).str.upper(),
                   hit[target_col].astype(str).str.upper()))


def _literature_tag(assessment) -> str:
    """Map an "Evidence assessment" value to "direct", "plausible" or "".
    "No direct evidence" must not count as direct, hence the exact prefixes."""
    a = str(assessment).strip().lower()
    if a.startswith("evidence found (direct)"):
        return "direct"
    if a.startswith("plausible"):
        return "plausible"
    return ""


def load_literature_tags(path: Path) -> dict[tuple[str, str], str]:
    """{(TF, target) upper-cased: "direct" | "plausible"} from one
    literature CSV. Unrecognized assessments are logged and left untagged."""
    df = pd.read_csv(path)
    tf_col     = _pick_col(df, TF_KEY_CANDIDATES)
    target_col = _pick_col(df, TARGET_KEY_CANDIDATES)
    tags = df[LITERATURE_COL].map(_literature_tag)

    other = sorted(df.loc[tags == "", LITERATURE_COL].dropna().astype(str).unique())
    known_untagged = {"no direct evidence"}
    other = [o for o in other if o.strip().lower() not in known_untagged]
    if other:
        print(f"  [note] {path.name}: untagged assessments {other}")

    keep = tags != ""
    print(f"  Literature ({path.name}): {int((tags == 'direct').sum())} direct, "
          f"{int((tags == 'plausible').sum())} plausible of {len(df)} pairs")
    return dict(zip(zip(df.loc[keep, tf_col].astype(str).str.upper(),
                        df.loc[keep, target_col].astype(str).str.upper()),
                    tags[keep]))


def load_literature_dir(literature_dir: str | Path) -> dict[tuple[str, str], str]:
    """Pool load_literature_tags over every CSV in `literature_dir`. Files
    without tf / target / assessment columns are skipped; when two files tag
    the same pair differently, direct wins over plausible."""
    literature_dir = Path(literature_dir)
    csvs = sorted(literature_dir.glob("*.csv"))
    print(f"\nLiterature evidence: {len(csvs)} CSV(s) in {literature_dir}")
    pooled: dict[tuple[str, str], str] = {}
    for path in csvs:
        try:
            tags = load_literature_tags(path)
        except KeyError as e:
            print(f"  [skip] {path.name}: {e}")
            continue
        for pair, tag in tags.items():
            if pooled.get(pair, tag) != tag:
                print(f"  [note] {pair[0]}->{pair[1]}: {pooled[pair]} vs {tag} "
                      f"({path.name}); keeping direct")
                tag = "direct"
            pooled[pair] = tag
    n_direct = sum(t == "direct" for t in pooled.values())
    print(f"  -> {n_direct} direct, {len(pooled) - n_direct} plausible pairs in total")
    return pooled


def select_top_tfs(df: pd.DataFrame, top_n: int) -> list[str]:
    """Same ranking as tf_network: most targets, then highest mean score."""
    return (
        df.groupby("tf")
        .agg(n_targets=("target", "count"), mean_score=("normalized_score", "mean"))
        .sort_values(["n_targets", "mean_score"], ascending=False)
        .head(top_n)
        .index.tolist()
    )


def select_plot_pairs(
    df: pd.DataFrame, groups: dict[str, str], top_n: int, validated_only: bool,
) -> tuple[pd.DataFrame, list[str]]:
    """The rows to draw (with a `_groups` column) and the top TFs.

    With `validated_only`, TFs without a validated edge are dropped before
    ranking, so the figure still gets up to `top_n` TFs, and then targets
    without a validated edge to one of those TFs are dropped. Unvalidated
    edges between remaining nodes are kept."""
    df = df.assign(_groups=edge_groups(df, groups))
    validated = df["_groups"].map(bool)
    if validated_only:
        df = df[df["tf"].isin(set(df.loc[validated, "tf"]))]
    top_tfs = select_top_tfs(df, top_n)
    df = df[df["tf"].isin(top_tfs)]
    if validated_only:
        df = df[df["target"].isin(set(df.loc[df["_groups"].map(bool), "target"]))]
    return df.copy(), top_tfs


MAX_TFS_PER_FIGURE = 5
TFS_PER_PART       = 3


def split_tfs(tfs: list[str]) -> list[list[str]]:
    """Keep up to MAX_TFS_PER_FIGURE TFs in one figure; otherwise split them,
    in rank order, into the fewest parts of at most TFS_PER_PART with sizes as
    even as possible (12 -> 3+3+3+3, 7 -> 3+2+2)."""
    if len(tfs) <= MAX_TFS_PER_FIGURE:
        return [list(tfs)]
    n_parts = -(-len(tfs) // TFS_PER_PART)
    return [[str(t) for t in part] for part in np.array_split(tfs, n_parts)]


def restrict_to_tfs(plot_df: pd.DataFrame, tfs: list[str], validated_only: bool) -> pd.DataFrame:
    """Rows of `plot_df` for `tfs`; with `validated_only`, also drop targets
    without a validated edge to one of these TFs."""
    sub = plot_df[plot_df["tf"].isin(tfs)]
    if validated_only:
        sub = sub[sub["target"].isin(set(sub.loc[sub["_groups"].map(bool), "target"]))]
    return sub


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------

def _draw_validation_edge(ax, start, end, groups, line_specs, node_r,
                          motif_hit=False, literature=""):
    """Draw one TF -> target edge as parallel lines, one per validating group.
    The bundle is kept narrower than the target node it ends on. An
    unvalidated edge with an accessible motif hit gets a star at its midpoint; a
    literature tag gets a diamond 3/4 of the way to the target."""
    d    = end - start
    dist = np.linalg.norm(d)
    if dist < 1e-6:
        return
    unit = d / dist
    perp = np.array([-unit[1], unit[0]])

    if literature:
        p = start + 0.75 * d
        ax.plot(p[0], p[1], zorder=2.5, **LITERATURE_STYLES[literature])

    if not groups:
        ax.plot([start[0], end[0]], [start[1], end[1]],
                zorder=2, solid_capstyle="butt", **UNVALIDATED_STYLE)
        if motif_hit:
            mid = (start + end) / 2
            ax.plot(mid[0], mid[1], zorder=2.5, **MOTIF_HIT_STYLE)
        return

    n       = len(groups)
    spacing = min(0.25 * node_r, 1.6 * node_r / (n - 1)) if n > 1 else 0.0
    offsets = (np.arange(n) - (n - 1) / 2) * spacing
    for off, g in zip(offsets, groups):
        s, e = start + perp * off, end + perp * off
        ax.plot([s[0], e[0]], [s[1], e[1]], lw=EDGE_LW, zorder=2,
                solid_capstyle="butt", dash_capstyle="butt", **line_specs[g])


def _add_legends(fig, groups, line_specs, has_unvalidated, title,
                 show_motif_hits=False, show_literature=False):
    """One legend below the figure (savefig's tight bbox grows to include
    it): the drawn cell types, then the edge marks."""
    handles = [Line2D([], [], lw=2.0, **line_specs[ct]) for ct in groups]
    labels  = list(groups)
    if has_unvalidated:
        handles.append(Line2D([], [], **dict(UNVALIDATED_STYLE, lw=1.5)))
        labels.append("not validated")
    if show_motif_hits and has_unvalidated:
        handles.append(Line2D([], [], **MOTIF_HIT_STYLE))
        labels.append(MOTIF_HIT_LABEL)
    if show_literature:
        for tag, style in LITERATURE_STYLES.items():
            handles.append(Line2D([], [], **style))
            labels.append(LITERATURE_LABELS[tag])
    fig.legend(handles, labels, title=title, loc="upper center",
               bbox_to_anchor=(0.5, 0.0), fontsize=8, title_fontsize=9,
               frameon=False, handlelength=3.0, ncol=2 if len(handles) > 6 else 1)


def _add_colorbar(fig, ax, **kw):
    sm = plt.cm.ScalarMappable(cmap=_TF_CMAP, norm=plt.Normalize(0, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, **kw)
    cbar.set_label("Normalized\nPrediction score", fontsize=9)
    cbar.set_ticks([0, 0.5, 1])


def plot_validation_network_connected(
    plot_df: pd.DataFrame,
    top_tfs: list[str],
    fig_path: Path,
    title: str,
    groups: dict[str, str],
    line_specs: dict[str, dict],
    legend_title: str,
) -> None:
    """Single connected graph (shared targets appear once) of `top_tfs` and
    the rows of `plot_df` (from select_plot_pairs), mirroring
    tf_network.plot_tf_network_connected - same node insertion order and
    spring-layout seed."""
    import networkx as nx

    G = nx.DiGraph()
    for tf in top_tfs:
        G.add_node(tf, node_type="tf")

    target_max_score: dict[str, float] = {}
    for _, row in plot_df.iterrows():
        tgt   = row["target"]
        score = float(row["normalized_score"])
        if tgt not in target_max_score or score > target_max_score[tgt]:
            target_max_score[tgt] = score
        if tgt not in G:
            G.add_node(tgt, node_type="gene")
        G.add_edge(row["tf"], tgt, weight=score, groups=row["_groups"],
                   motif_hit=bool(row.get("_motif_hit", False)),
                   literature=row.get("_literature", ""))

    pos = nx.spring_layout(G, seed=42, k=3.0, iterations=150)

    fig, ax = plt.subplots(figsize=(14, 11))
    ax.set_aspect("equal")
    ax.axis("off")

    tf_nodes   = [n for n in G.nodes if G.nodes[n]["node_type"] == "tf"]
    gene_nodes = [n for n in G.nodes if G.nodes[n]["node_type"] == "gene"]

    tf_r   = 0.07
    gene_r = 0.035

    # ── Edges ────────────────────────────────────────────────────────────────
    for u, v, data in G.edges(data=True):
        pu, pv = np.array(pos[u]), np.array(pos[v])
        d      = pv - pu
        dist   = np.linalg.norm(d)
        if dist < 1e-6:
            continue
        unit = d / dist
        _draw_validation_edge(ax, pu + unit * tf_r, pv - unit * gene_r,
                              data["groups"], line_specs, gene_r,
                              data["motif_hit"], data["literature"])

    # ── Gene (target) nodes ──────────────────────────────────────────────────
    for gene in gene_nodes:
        p     = np.array(pos[gene])
        color = _TF_CMAP(target_max_score.get(gene, 0.0))
        ax.add_patch(plt.Circle(p, gene_r, color=color, ec="none", zorder=3))

        # Push label away from the TF centroid
        tf_preds = list(G.predecessors(gene))
        if tf_preds:
            centroid  = np.mean([pos[t] for t in tf_preds], axis=0)
            push      = p - centroid
            push_norm = np.linalg.norm(push)
            push_dir  = push / push_norm if push_norm > 1e-6 else np.array([1.0, 0.0])
        else:
            push_dir = np.array([1.0, 0.0])

        lp = p + push_dir * (gene_r + 0.028)
        ha = "left" if push_dir[0] >= 0 else "right"
        ax.text(lp[0], lp[1], gene, ha=ha, va="center",
                fontsize=7, fontweight="bold", color="#222222", zorder=5)

    # ── TF hub nodes ─────────────────────────────────────────────────────────
    for tf in tf_nodes:
        p = np.array(pos[tf])
        ax.add_patch(plt.Circle(p, tf_r, color="white", ec="#333333", lw=2.0, zorder=4))
        ax.text(p[0], p[1], tf, ha="center", va="center",
                fontsize=9, fontweight="bold", color="#111111", zorder=6)

    _add_colorbar(fig, ax, shrink=0.28, pad=0.02, aspect=18, anchor=(1.0, 0.5))

    all_pos = np.array(list(pos.values()))
    margin  = 0.30
    ax.set_xlim(all_pos[:, 0].min() - margin, all_pos[:, 0].max() + margin)
    ax.set_ylim(all_pos[:, 1].min() - margin, all_pos[:, 1].max() + margin)

    has_unvalidated = any(not g for g in plot_df["_groups"])
    fig.suptitle(title, fontsize=12)
    plt.tight_layout()
    _add_legends(fig, groups, line_specs, has_unvalidated, legend_title,
                 "_motif_hit" in plot_df.columns, "_literature" in plot_df.columns)
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fig_path.name}")


def plot_validation_network_hub_spoke(
    plot_df: pd.DataFrame,
    top_tfs: list[str],
    fig_path: Path,
    title: str,
    groups: dict[str, str],
    line_specs: dict[str, dict],
    legend_title: str,
) -> None:
    """One hub-and-spoke subplot per TF in `top_tfs`, with the rows of
    `plot_df` (from select_plot_pairs), mirroring tf_network.plot_tf_network."""
    n_tfs   = len(top_tfs)
    ncols   = min(3, n_tfs)
    nrows   = (n_tfs + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5.5, nrows * 5.5),
                             squeeze=False)
    axes_flat = axes.flatten()

    has_unvalidated = False
    for idx, tf in enumerate(top_tfs):
        ax  = axes_flat[idx]
        sub = plot_df[plot_df["tf"] == tf].sort_values("normalized_score", ascending=False)
        targets  = sub["target"].tolist()
        scores   = sub["normalized_score"].tolist()
        validity = sub["_groups"].tolist()
        motif    = (sub["_motif_hit"].tolist() if "_motif_hit" in sub.columns
                    else [False] * len(sub))
        lit      = (sub["_literature"].tolist() if "_literature" in sub.columns
                    else [""] * len(sub))
        has_unvalidated |= any(not g for g in validity)
        n = len(targets)

        radius      = max(0.9, 0.15 * n)
        hub_radius  = min(0.22, radius * 0.22)
        node_radius = min(0.09, radius * 0.09)
        center      = np.array([0.0, 0.0])

        angles     = np.linspace(0, 2 * np.pi, n, endpoint=False)
        target_pos = center + radius * np.stack([np.cos(angles), np.sin(angles)], axis=1)

        for pos, score, gene, grp, hit, tag in zip(
                target_pos, scores, targets, validity, motif, lit):
            unit = (pos - center) / np.linalg.norm(pos - center)
            _draw_validation_edge(ax, center + unit * hub_radius,
                                  pos - unit * node_radius,
                                  grp, line_specs, node_radius, hit, tag)
            ax.add_patch(plt.Circle(pos, node_radius, color=_TF_CMAP(score),
                                    ec="none", zorder=3))
            label_pos = center + (radius + node_radius + 0.05) * unit
            ax.text(label_pos[0], label_pos[1], gene,
                    ha="left" if pos[0] >= center[0] else "right", va="center",
                    fontsize=7, fontweight="bold", color="#222222", zorder=5)

        ax.add_patch(plt.Circle(center, hub_radius, color="white",
                                ec="#333333", lw=1.8, zorder=4))
        ax.text(center[0], center[1], tf, ha="center", va="center",
                fontsize=8, fontweight="bold", color="#111111", zorder=6)

        margin = radius + node_radius + 0.4
        ax.set_xlim(-margin, margin)
        ax.set_ylim(-margin, margin)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"{tf}  ({n} targets)", fontsize=8, pad=4)

    for idx in range(n_tfs, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    fig.suptitle(title, fontsize=12, y=1.01)
    # Colorbar after tight_layout, which otherwise lays a subplot over it.
    plt.tight_layout()
    _add_colorbar(fig, axes_flat[:n_tfs].tolist(), shrink=0.35, pad=0.04, aspect=20)
    _add_legends(fig, groups, line_specs, has_unvalidated, legend_title,
                 "_motif_hit" in plot_df.columns, "_literature" in plot_df.columns)
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fig_path.name}")


PLOT_STYLES = {
    "connected": plot_validation_network_connected,
    "hub_spoke": plot_validation_network_hub_spoke,
}


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def plot_validation_run(
    run_dir: Path,
    plot_style: str = "connected",
    top_n: int = 8,
    motif_hit_pairs: set[tuple[str, str]] | None = None,
    literature_tags: dict[tuple[str, str], str] | None = None,
    cell_class_json: str | Path = DEFAULT_CELL_CLASS_JSON,
    cell_class: str = DEFAULT_CELL_CLASS,
    validated_only: bool = True,
) -> None:
    """Plot every tf_targets_*.csv in one tf_validation subfolder, drawing the
    cell types of `cell_class`. `motif_hit_pairs` (from
    load_accessible_motif_pairs) adds the stars; `literature_tags` (from
    load_literature_dir) adds the literature diamonds."""
    csvs = sorted(run_dir.glob("tf_targets_*.csv"))
    if not csvs:
        print(f"  [skip] no tf_targets_*.csv in {run_dir}")
        return

    fig_dir = run_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    plot_fn = PLOT_STYLES[plot_style]

    cell_types   = load_cell_class(cell_class_json, cell_class)
    line_specs   = build_line_specs(cell_types)
    legend_title = f"Validated in {cell_class.replace('_', ' ')}"

    print(f"\nPlotting {len(csvs)} CSV(s) from {run_dir}")
    # Every CSV of one run has the same columns, so log the matching once.
    header  = pd.read_csv(csvs[0], nrows=0)
    groups  = detect_groups(header, cell_types)
    missing = [ct for ct in cell_types if ct not in groups]
    print(f"  {cell_class}: drawing {len(groups)}/{len(cell_types)} cell types "
          f"({', '.join(groups) or 'none'})")
    if missing:
        print(f"  [note] no validated_* column for: {', '.join(missing)} "
              f"(not scored in the ATAC run, or < 50 cells)")

    for csv_path in csvs:
        gene_list = csv_path.stem.removeprefix("tf_targets_")
        df = pd.read_csv(csv_path)
        if df.empty:
            print(f"  [skip] {csv_path.name}: empty")
            continue
        if not detect_groups(df, cell_types):
            print(f"  [skip] {csv_path.name}: none of the {cell_class} cell types "
                  f"has a validated_* column")
            continue
        keys = list(zip(df["tf"].astype(str).str.upper(),
                        df["target"].astype(str).str.upper()))
        if motif_hit_pairs is not None:
            df["_motif_hit"] = [k in motif_hit_pairs for k in keys]
        if literature_tags:
            df["_literature"] = [literature_tags.get(k, "") for k in keys]
            n_tagged = int((df["_literature"] != "").sum())
            print(f"  {gene_list}: {n_tagged} pair(s) with literature evidence")
        groups = detect_groups(df, cell_types)
        plot_df, top_tfs = select_plot_pairs(df, groups, top_n, validated_only)

        # Clear this gene list's figures from earlier runs, so a split that
        # changed size doesn't leave stale parts next to the new ones.
        for old in fig_dir.glob(f"tf_validation_network_{gene_list}.png"):
            old.unlink()
        for old in fig_dir.glob(f"tf_validation_network_{gene_list}_part*of*.png"):
            old.unlink()

        if plot_df.empty:
            print(f"  [skip] {gene_list}: no validated pairs to draw")
            continue

        parts = split_tfs(top_tfs)
        title = f"TF Regulatory Networks — {gene_list} (ATAC validation)"
        if len(parts) > 1:
            print(f"  {gene_list}: {len(top_tfs)} TFs -> {len(parts)} figures "
                  f"({'+'.join(str(len(p)) for p in parts)})")
        for i, part in enumerate(parts, 1):
            suffix   = f"_part{i}of{len(parts)}" if len(parts) > 1 else ""
            part_ttl = (f"{title}\npart {i}/{len(parts)}: {', '.join(part)}"
                        if len(parts) > 1 else title)
            plot_fn(
                restrict_to_tfs(plot_df, part, validated_only),
                part,
                fig_dir / f"tf_validation_network_{gene_list}{suffix}.png",
                title=part_ttl,
                groups=groups,
                line_specs=line_specs,
                legend_title=legend_title,
            )


def run(
    validation_dir: str,
    plot_style: str = "connected",
    top_n: int = 8,
    motif_hits_csv: str | None = None,
    literature_dir: str | None = None,
    cell_class_json: str | Path = DEFAULT_CELL_CLASS_JSON,
    cell_class: str = DEFAULT_CELL_CLASS,
    validated_only: bool = True,
) -> None:
    root = Path(validation_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"validation dir not found: {root}")

    # A single run folder, or the tf_validation root holding many of them.
    run_dirs = ([root] if any(root.glob("tf_targets_*.csv"))
                else sorted(d for d in root.iterdir()
                            if d.is_dir() and any(d.glob("tf_targets_*.csv"))))
    if not run_dirs:
        raise FileNotFoundError(f"No tf_targets_*.csv files under {root}")

    motif_hit_pairs = load_accessible_motif_pairs(motif_hits_csv) if motif_hits_csv else None
    literature_tags = load_literature_dir(literature_dir) if literature_dir else None
    for d in run_dirs:
        plot_validation_run(d, plot_style=plot_style, top_n=top_n,
                            motif_hit_pairs=motif_hit_pairs,
                            literature_tags=literature_tags,
                            cell_class_json=cell_class_json, cell_class=cell_class,
                            validated_only=validated_only)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_arg_parser():
    p = argparse.ArgumentParser(
        description=(
            "Plot tf_validation_report's annotated CSVs as TF networks whose "
            "edges show the cell types of one cell class that validate each "
            "TF-target pair."
        )
    )
    p.add_argument("--validation-dir", default="results/tf_validation",
                   help="A tf_validation subfolder, or the root holding several "
                        "(default: results/tf_validation).")
    p.add_argument("--plot-style", choices=sorted(PLOT_STYLES), default="connected",
                   help="Layout, as tf_network's query.plot_style (default: connected).")
    p.add_argument("--top-n-tfs", type=int, default=8,
                   help="Top N TFs to show, as tf_network's query.top_n_tfs (default: 8).")
    p.add_argument("--motif-hits-csv", default=None,
                   help="motif_target_pair_scores.csv from the ATAC run's "
                        "3_motif_target_validation/. Stars unvalidated (grey) "
                        "edges whose pair has a motif hit in an accessible "
                        "peak (validated_by_motif_and_accessibility, all peaks).")
    p.add_argument("--literature-dir", default=None,
                   help="Folder of literature CSVs, pooled across gene lists. "
                        "Tags pairs with direct (filled diamond) or plausible "
                        "(hollow diamond) literature evidence.")
    p.add_argument("--cell-class-json", default=str(DEFAULT_CELL_CLASS_JSON),
                   help="JSON of {cell class: [ATAC cell type, ...]} "
                        "(default: config_files/atac_cell_classes.json).")
    p.add_argument("--cell-class", default=DEFAULT_CELL_CLASS,
                   help=f"Class whose cell types are drawn (default: {DEFAULT_CELL_CLASS}).")
    p.add_argument("--show-unvalidated-nodes", action="store_true",
                   help="Keep TFs and targets without a validated edge (by "
                        "default they are dropped).")
    return p


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    run(args.validation_dir, plot_style=args.plot_style, top_n=args.top_n_tfs,
        motif_hits_csv=args.motif_hits_csv, literature_dir=args.literature_dir,
        cell_class_json=args.cell_class_json, cell_class=args.cell_class,
        validated_only=not args.show_unvalidated_nodes)
