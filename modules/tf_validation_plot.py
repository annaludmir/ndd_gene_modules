"""
tf_validation_plot.py

Network figures for tf_validation_report's annotated CSVs. Same layouts and
TF selection as tf_network's figures, and the target node colour is still the
normalized prediction score, but every TF -> target edge is drawn as a bundle
of parallel lines, one per group in which ATAC-seq validates that pair:

  - cell types   -> solid lines, one colour per major cell class
  - cell phases  -> dashed lines, one colour per phase
  - not validated anywhere -> a single thin grey dotted line

With --motif-hits-csv (the ATAC run's global
`3_motif_target_validation/motif_target_pair_scores.csv`), grey edges whose
pair has a TF motif hit in an accessible peak of the target's promoter when
all peaks are considered (validated_by_motif_and_accessibility) get a star at
the edge midpoint: accessible, but not among the top peaks of any cell type /
phase. Coloured edges aren't starred - the per-group peaks are a subset of
all peaks, so they are always globally validated too.

With --literature-dir, pairs are also tagged from the curated literature
CSV of their gene list, `<gene_list>_TF_Validation_with_evidence.csv`
(column "Evidence assessment"), with a diamond 3/4 of the way to the target:
  - "Evidence found (direct)"            -> filled diamond
  - "Plausible, no pair-specific ..."    -> hollow diamond
Other assessments are left untagged. Gene lists without a literature CSV are
plotted without tags.

The ATAC cell types are too many to tell apart by colour, so they are folded
into major classes (CELL_CLASS_OF below); a pair is validated in a class if
it is validated in any of the class's cell types. Cell types missing from the
mapping are plotted on their own and reported in the log.

Colours are assigned in a fixed order (CELL_CLASS_ORDER, then phase cycle
order), so a class / phase keeps its colour across all gene-list figures.

Inputs:
  --validation-dir <path>
     A tf_validation subfolder holding `tf_targets_<gene_list>.csv` files, or
     the results/tf_validation root, in which case every subfolder is plotted.

  --motif-hits-csv <path>   (optional)
     motif_target_pair_scores.csv from the matching atac_seq_analysis run.

  --literature-dir <path>   (optional)
     Folder of <gene_list>_TF_Validation_with_evidence.csv files, e.g.
     results/tf_validation/literature_evidence.

Output:
  <subfolder>/figures/tf_validation_network_<gene_list>.png

Usage:
  python modules/tf_validation_plot.py \
    --validation-dir results/tf_validation/cortex_v3_20260802_vs_atac_20260816 \
    --plot-style connected --top-n-tfs 8 \
    --motif-hits-csv results/atac_analysis/<run>/3_motif_target_validation/motif_target_pair_scores.csv
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from tf_network import _TF_CMAP
from tf_validation_report import (
    AXES, TARGET_KEY_CANDIDATES, TF_KEY_CANDIDATES, VALIDATED_COL, _pick_col,
)


# Major cell classes the ATAC cell types are folded into, keyed by the
# sanitized slug used in the CSV column names (matched case-insensitively).
CELL_CLASS_OF = {
    "radial_glial_cell":                          "Radial glia / progenitors",
    "neural_progenitor_cell":                     "Radial glia / progenitors",
    "progenitor_cell":                            "Radial glia / progenitors",
    "neuroblast_sensu_nematoda_and_protostomia":  "Neuroblast",
    "glutamatergic_neuron":                       "Neuron",
    "interneuron":                                "Neuron",
    "purkinje_cell":                              "Neuron",
    "glycinergic_neuron":                         "Neuron",
    "dopaminergic_neuron":                        "Neuron",
    "serotonergic_neuron":                        "Neuron",
    "sensory_neuron_of_dorsal_root_ganglion":     "Neuron",
    "glioblast":                                  "Glioblast",
    "oligodendrocyte_precursor_cell":             "Oligodendrocyte lineage",
    "committed_oligodendrocyte_precursor":        "Oligodendrocyte lineage",
    "oligodendrocyte":                            "Oligodendrocyte lineage",
    "schwann_cell":                               "Neural crest",
    "endothelial_cell":                           "Vascular",
    "pericyte":                                   "Vascular",
    "vascular_associated_smooth_muscle_cell":     "Vascular",
    "vascular_leptomeningeal_cell":               "Vascular",
    "microglial_cell":                            "Immune",
}
CELL_CLASS_ORDER = [
    "Radial glia / progenitors", "Neuroblast", "Neuron", "Glioblast",
    "Oligodendrocyte lineage", "Neural crest", "Vascular", "Immune",
]

# Cell classes: categorical palette in CELL_CLASS_ORDER; unmapped cell types
# fall back to the light tab20 variants.
CELL_TYPE_COLORS = [
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
    "#e87ba4", "#008300", "#4a3aa7", "#e34948",
]
CELL_TYPE_FALLBACK = [matplotlib.colormaps["tab20"](i) for i in range(1, 20, 2)]

# Cell phases: a separate set of hues that don't repeat the cell-type ones,
# in cycle order (Non-cycling -> G1 -> S -> G2M -> Post-M).
PHASE_COLORS = ["#222222", "#17becf", "#8c564b", "#bcbd22", "#1f3b73"]
PHASE_FALLBACK = ["#7f7f7f", "#b5651d", "#5f9ea0"]

LINESTYLES = {"cell_type": "-", "cell_phase": (0, (4, 2))}
UNVALIDATED_STYLE = dict(color="#bbbbbb", lw=0.7, linestyle=(0, (1, 2)))
EDGE_LW = 1.1
MOTIF_HIT_STYLE = dict(marker="*", markersize=9, color="#111111",
                       markeredgecolor="white", markeredgewidth=0.6, linestyle="none")
MOTIF_HIT_LABEL = "not validated, accessible motif hit (all peaks)"

LITERATURE_SUFFIX = "_TF_Validation_with_evidence.csv"
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

AXIS_LEGEND_TITLE = {
    "cell_type": "Validated in cell class",
    "cell_phase": "Validated in cell phase",
}


# ---------------------------------------------------------------------------
# Reading the annotated CSV
# ---------------------------------------------------------------------------

def detect_slugs(df: pd.DataFrame) -> dict[str, list[str]]:
    """Return {axis: [group_slug, ...]} from the per-group boolean columns,
    in CSV column order. Axes absent from the CSV are omitted."""
    reserved = {spec[k] for spec in AXES.values()
                for k in ("any_col", "count_col", "list_col")}
    phase_prefix = AXES["cell_phase"]["col_prefix"]
    ct_prefix    = AXES["cell_type"]["col_prefix"]

    groups = {"cell_type": [], "cell_phase": []}
    for c in df.columns:
        if c in reserved:
            continue
        if c.startswith(phase_prefix):
            groups["cell_phase"].append(c[len(phase_prefix):])
        elif c.startswith(ct_prefix):
            groups["cell_type"].append(c[len(ct_prefix):])
    return {a: g for a, g in groups.items() if g}


def unmapped_cell_types(slugs: dict[str, list[str]]) -> list[str]:
    return [s for s in slugs.get("cell_type", []) if s.lower() not in CELL_CLASS_OF]


def detect_groups(df: pd.DataFrame) -> dict[str, dict[str, list[str]]]:
    """Return {axis: {legend_label: [group_slug, ...]}}. Cell types are folded
    into their major class (unmapped ones keep their own entry, after the
    classes); phases map one-to-one, in cycle order."""
    slugs  = detect_slugs(df)
    groups = {}
    if "cell_type" in slugs:
        classes: dict[str, list[str]] = {}
        for c in CELL_CLASS_ORDER:
            members = [s for s in slugs["cell_type"] if CELL_CLASS_OF.get(s.lower()) == c]
            if members:
                classes[c] = members
        for s in unmapped_cell_types(slugs):
            classes[s.replace("_", " ")] = [s]
        groups["cell_type"] = classes
    if "cell_phase" in slugs:
        groups["cell_phase"] = {s: [s] for s in slugs["cell_phase"]}
    return groups


def build_line_specs(groups: dict[str, dict[str, list[str]]]) -> dict[tuple[str, str], dict]:
    """{(axis, label): {color, linestyle}} for every legend entry."""
    palettes = {
        "cell_type": CELL_TYPE_COLORS + CELL_TYPE_FALLBACK,
        "cell_phase": PHASE_COLORS + PHASE_FALLBACK,
    }
    specs = {}
    for axis, labels in groups.items():
        palette = palettes[axis]
        if len(labels) > len(palette):
            raise ValueError(
                f"{len(labels)} {AXES[axis]['label']} groups but only "
                f"{len(palette)} colours available - extend CELL_CLASS_OF."
            )
        for label, color in zip(labels, palette):
            specs[(axis, label)] = dict(color=color, linestyle=LINESTYLES[axis])
    return specs


def edge_groups(df: pd.DataFrame, groups: dict[str, dict[str, list[str]]]) -> list[list[tuple[str, str]]]:
    """Per row, the (axis, label) groups that validate the pair - cell
    classes first, then phases, each in legend order. A class validates the
    pair if any of its cell types does."""
    keys, flags = [], []
    for axis, labels in groups.items():
        prefix = AXES[axis]["col_prefix"]
        for label, members in labels.items():
            keys.append((axis, label))
            cols = [f"{prefix}{m}" for m in members]
            flags.append(df[cols].fillna(False).astype(bool).any(axis=1).to_numpy())
    matrix = np.column_stack(flags)
    return [[keys[i] for i in np.where(row)[0]] for row in matrix]


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


def find_literature_csv(literature_dir: Path, gene_list: str) -> Path | None:
    """`<gene_list>_TF_Validation_with_evidence.csv`, matched case-insensitively."""
    want = f"{gene_list}{LITERATURE_SUFFIX}".lower()
    for p in literature_dir.glob("*.csv"):
        if p.name.lower() == want:
            return p
    return None


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


def select_top_tfs(df: pd.DataFrame, top_n: int) -> list[str]:
    """Same ranking as tf_network: most targets, then highest mean score."""
    return (
        df.groupby("tf")
        .agg(n_targets=("target", "count"), mean_score=("normalized_score", "mean"))
        .sort_values(["n_targets", "mean_score"], ascending=False)
        .head(top_n)
        .index.tolist()
    )


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


def _add_legends(fig, groups, line_specs, has_unvalidated, show_motif_hits=False,
                 show_literature=False):
    """One legend per validation axis, side by side below the figure
    (savefig's tight bbox grows to include them)."""
    entries = []
    for axis, names in groups.items():
        handles = [Line2D([], [], lw=2.0, **line_specs[(axis, n)]) for n in names]
        labels  = [n.replace("_", " ") for n in names]
        # Columns are sized by the groups alone, so the edge-mark rows
        # appended below don't push a short legend into two columns.
        ncol = 2 if len(names) > 7 else 1
        entries.append((AXIS_LEGEND_TITLE[axis], handles, labels, ncol))
    # Edge marks that aren't a group go at the end of the last legend, so
    # they don't need a box of their own.
    if has_unvalidated:
        style = dict(UNVALIDATED_STYLE, lw=1.5)
        entries[-1][1].append(Line2D([], [], **style))
        entries[-1][2].append("not validated")
    if show_motif_hits and has_unvalidated:
        entries[-1][1].append(Line2D([], [], **MOTIF_HIT_STYLE))
        entries[-1][2].append(MOTIF_HIT_LABEL)
    if show_literature:
        for tag, style in LITERATURE_STYLES.items():
            entries[-1][1].append(Line2D([], [], **style))
            entries[-1][2].append(LITERATURE_LABELS[tag])

    placements = (
        [("upper center", 0.5)] if len(entries) == 1
        else [("upper right", 0.49), ("upper left", 0.51)]
    )
    for (title, handles, labels, ncol), (loc, x) in zip(entries, placements):
        fig.legend(handles, labels, title=title, loc=loc,
                   bbox_to_anchor=(x, 0.0), fontsize=8, title_fontsize=9,
                   frameon=False, handlelength=3.0,
                   ncol=ncol)


def _add_colorbar(fig, ax, **kw):
    sm = plt.cm.ScalarMappable(cmap=_TF_CMAP, norm=plt.Normalize(0, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, **kw)
    cbar.set_label("Normalized\nPrediction score", fontsize=9)
    cbar.set_ticks([0, 0.5, 1])


def plot_validation_network_connected(
    df: pd.DataFrame,
    fig_path: Path,
    top_n: int = 8,
    title: str = "TF Regulatory Networks",
) -> None:
    """Single connected graph (shared targets appear once), mirroring
    tf_network.plot_tf_network_connected - same node insertion order and
    spring-layout seed, so node positions match the original figure."""
    import networkx as nx

    groups     = detect_groups(df)
    line_specs = build_line_specs(groups)

    top_tfs = select_top_tfs(df, top_n)
    plot_df = df[df["tf"].isin(top_tfs)].copy()
    plot_df["_groups"] = edge_groups(plot_df, groups)

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
    _add_legends(fig, groups, line_specs, has_unvalidated,
                 "_motif_hit" in df.columns, "_literature" in df.columns)
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fig_path.name}")


def plot_validation_network_hub_spoke(
    df: pd.DataFrame,
    fig_path: Path,
    top_n: int = 8,
    title: str = "TF Regulatory Networks",
) -> None:
    """One hub-and-spoke subplot per TF, mirroring tf_network.plot_tf_network."""
    groups     = detect_groups(df)
    line_specs = build_line_specs(groups)

    top_tfs = select_top_tfs(df, top_n)
    n_tfs   = len(top_tfs)
    ncols   = min(3, n_tfs)
    nrows   = (n_tfs + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5.5, nrows * 5.5),
                             squeeze=False)
    axes_flat = axes.flatten()

    has_unvalidated = False
    for idx, tf in enumerate(top_tfs):
        ax  = axes_flat[idx]
        sub = df[df["tf"] == tf].sort_values("normalized_score", ascending=False)
        targets  = sub["target"].tolist()
        scores   = sub["normalized_score"].tolist()
        validity = edge_groups(sub, groups)
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
    _add_legends(fig, groups, line_specs, has_unvalidated,
                 "_motif_hit" in df.columns, "_literature" in df.columns)
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
    literature_dir: str | Path | None = None,
) -> None:
    """Plot every tf_targets_*.csv in one tf_validation subfolder.
    `motif_hit_pairs` (from load_accessible_motif_pairs) adds the stars;
    `literature_dir` adds literature tags for gene lists that have a CSV."""
    csvs = sorted(run_dir.glob("tf_targets_*.csv"))
    if not csvs:
        print(f"  [skip] no tf_targets_*.csv in {run_dir}")
        return

    fig_dir = run_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    plot_fn = PLOT_STYLES[plot_style]

    print(f"\nPlotting {len(csvs)} CSV(s) from {run_dir}")
    # Every CSV of one run has the same columns, so log the folding once.
    header = pd.read_csv(csvs[0], nrows=0)
    for label, members in detect_groups(header).get("cell_type", {}).items():
        print(f"  {label}: {', '.join(members)}")
    unmapped = unmapped_cell_types(detect_slugs(header))
    if unmapped:
        print(f"  [note] not in CELL_CLASS_OF, plotted on their own: {', '.join(unmapped)}")

    for csv_path in csvs:
        gene_list = csv_path.stem.removeprefix("tf_targets_")
        df = pd.read_csv(csv_path)
        if df.empty:
            print(f"  [skip] {csv_path.name}: empty")
            continue
        if not detect_groups(df):
            print(f"  [skip] {csv_path.name}: no validated_* columns")
            continue
        keys = list(zip(df["tf"].astype(str).str.upper(),
                        df["target"].astype(str).str.upper()))
        if motif_hit_pairs is not None:
            df["_motif_hit"] = [k in motif_hit_pairs for k in keys]
        if literature_dir:
            lit_csv = find_literature_csv(Path(literature_dir), gene_list)
            if lit_csv is None:
                print(f"  [note] no {gene_list}{LITERATURE_SUFFIX} — no literature tags")
            else:
                tags = load_literature_tags(lit_csv)
                df["_literature"] = [tags.get(k, "") for k in keys]
        plot_fn(
            df,
            fig_dir / f"tf_validation_network_{gene_list}.png",
            top_n=top_n,
            title=f"TF Regulatory Networks — {gene_list} (ATAC validation)",
        )


def run(
    validation_dir: str,
    plot_style: str = "connected",
    top_n: int = 8,
    motif_hits_csv: str | None = None,
    literature_dir: str | None = None,
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
    for d in run_dirs:
        plot_validation_run(d, plot_style=plot_style, top_n=top_n,
                            motif_hit_pairs=motif_hit_pairs,
                            literature_dir=literature_dir)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_arg_parser():
    p = argparse.ArgumentParser(
        description=(
            "Plot tf_validation_report's annotated CSVs as TF networks whose "
            "edges show the cell types (solid) and cell phases (dashed) that "
            "validate each TF-target pair."
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
                   help="Folder of <gene_list>_TF_Validation_with_evidence.csv "
                        "files. Tags pairs with direct (filled diamond) or "
                        "plausible (hollow diamond) literature evidence.")
    return p


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    run(args.validation_dir, plot_style=args.plot_style, top_n=args.top_n_tfs,
        motif_hits_csv=args.motif_hits_csv, literature_dir=args.literature_dir)
