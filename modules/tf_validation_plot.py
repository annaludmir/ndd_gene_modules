"""
tf_validation_plot.py

Network figures for tf_validation_report's annotated CSVs. Same layouts and
TF selection as tf_network's figures, and the target node colour is still the
normalized prediction score, but every TF -> target edge is drawn as a bundle
of parallel lines, one per group in which ATAC-seq validates that pair:

  - cell types   -> solid lines, one colour per cell type
  - cell phases  -> dashed lines, one colour per phase
  - not validated anywhere -> a single thin grey dotted line

Group colours are assigned in CSV column order, which is identical across
every CSV of one tf_validation run, so a cell type / phase keeps its colour
across all gene-list figures of that run.

Inputs:
  --validation-dir <path>
     A tf_validation subfolder holding `tf_targets_<gene_list>.csv` files, or
     the results/tf_validation root, in which case every subfolder is plotted.

Output:
  <subfolder>/figures/tf_validation_network_<gene_list>.png

Usage:
  python modules/tf_validation_plot.py \
    --validation-dir results/tf_validation/cortex_v3_20260802_vs_atac_20260816 \
    --plot-style connected --top-n-tfs 8
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
from tf_validation_report import AXES


# Cell types: categorical palette in fixed order; extras (>8 cell types) fall
# back to the light tab20 variants so every cell type still gets its own colour.
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

AXIS_LEGEND_TITLE = {
    "cell_type": "Validated in cell type",
    "cell_phase": "Validated in cell phase",
}


# ---------------------------------------------------------------------------
# Reading the annotated CSV
# ---------------------------------------------------------------------------

def detect_groups(df: pd.DataFrame) -> dict[str, list[str]]:
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


def build_line_specs(groups: dict[str, list[str]]) -> dict[tuple[str, str], dict]:
    """{(axis, slug): {color, linestyle}} for every validation group."""
    palettes = {
        "cell_type": CELL_TYPE_COLORS + CELL_TYPE_FALLBACK,
        "cell_phase": PHASE_COLORS + PHASE_FALLBACK,
    }
    specs = {}
    for axis, slugs in groups.items():
        palette = palettes[axis]
        if len(slugs) > len(palette):
            raise ValueError(
                f"{len(slugs)} {AXES[axis]['label_plural']} but only "
                f"{len(palette)} colours available."
            )
        for slug, color in zip(slugs, palette):
            specs[(axis, slug)] = dict(color=color, linestyle=LINESTYLES[axis])
    return specs


def edge_groups(df: pd.DataFrame, groups: dict[str, list[str]]) -> list[list[tuple[str, str]]]:
    """Per row, the (axis, slug) groups that validate the pair - cell types
    first, then phases, each in column order."""
    keys, cols = [], []
    for axis, slugs in groups.items():
        for slug in slugs:
            keys.append((axis, slug))
            cols.append(f"{AXES[axis]['col_prefix']}{slug}")
    matrix = df[cols].fillna(False).to_numpy(dtype=bool)
    return [[keys[i] for i in np.where(row)[0]] for row in matrix]


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

def _draw_validation_edge(ax, start, end, groups, line_specs, node_r):
    """Draw one TF -> target edge as parallel lines, one per validating group.
    The bundle is kept narrower than the target node it ends on."""
    d    = end - start
    dist = np.linalg.norm(d)
    if dist < 1e-6:
        return
    unit = d / dist
    perp = np.array([-unit[1], unit[0]])

    if not groups:
        ax.plot([start[0], end[0]], [start[1], end[1]],
                zorder=2, solid_capstyle="butt", **UNVALIDATED_STYLE)
        return

    n       = len(groups)
    spacing = min(0.25 * node_r, 1.6 * node_r / (n - 1)) if n > 1 else 0.0
    offsets = (np.arange(n) - (n - 1) / 2) * spacing
    for off, g in zip(offsets, groups):
        s, e = start + perp * off, end + perp * off
        ax.plot([s[0], e[0]], [s[1], e[1]], lw=EDGE_LW, zorder=2,
                solid_capstyle="butt", dash_capstyle="butt", **line_specs[g])


def _add_legends(fig, groups, line_specs, has_unvalidated):
    """One legend per validation axis, side by side below the figure
    (savefig's tight bbox grows to include them)."""
    entries = []
    for axis, slugs in groups.items():
        handles = [Line2D([], [], lw=2.0, **line_specs[(axis, s)]) for s in slugs]
        labels  = [s.replace("_", " ") for s in slugs]
        entries.append((AXIS_LEGEND_TITLE[axis], handles, labels))
    if has_unvalidated:
        # Appended to the last legend so it doesn't need a box of its own.
        style = dict(UNVALIDATED_STYLE, lw=1.5)
        entries[-1][1].append(Line2D([], [], **style))
        entries[-1][2].append("not validated")

    placements = (
        [("upper center", 0.5)] if len(entries) == 1
        else [("upper right", 0.49), ("upper left", 0.51)]
    )
    for (title, handles, labels), (loc, x) in zip(entries, placements):
        fig.legend(handles, labels, title=title, loc=loc,
                   bbox_to_anchor=(x, 0.0), fontsize=8, title_fontsize=9,
                   frameon=False, handlelength=3.0,
                   ncol=2 if len(handles) > 6 else 1)


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
        G.add_edge(row["tf"], tgt, weight=score, groups=row["_groups"])

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
                              data["groups"], line_specs, gene_r)

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
    _add_legends(fig, groups, line_specs, has_unvalidated)
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
        has_unvalidated |= any(not g for g in validity)
        n = len(targets)

        radius      = max(0.9, 0.15 * n)
        hub_radius  = min(0.22, radius * 0.22)
        node_radius = min(0.09, radius * 0.09)
        center      = np.array([0.0, 0.0])

        angles     = np.linspace(0, 2 * np.pi, n, endpoint=False)
        target_pos = center + radius * np.stack([np.cos(angles), np.sin(angles)], axis=1)

        for pos, score, gene, grp in zip(target_pos, scores, targets, validity):
            unit = (pos - center) / np.linalg.norm(pos - center)
            _draw_validation_edge(ax, center + unit * hub_radius,
                                  pos - unit * node_radius,
                                  grp, line_specs, node_radius)
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
    _add_legends(fig, groups, line_specs, has_unvalidated)
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

def plot_validation_run(run_dir: Path, plot_style: str = "connected", top_n: int = 8) -> None:
    """Plot every tf_targets_*.csv in one tf_validation subfolder."""
    csvs = sorted(run_dir.glob("tf_targets_*.csv"))
    if not csvs:
        print(f"  [skip] no tf_targets_*.csv in {run_dir}")
        return

    fig_dir = run_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    plot_fn = PLOT_STYLES[plot_style]

    print(f"\nPlotting {len(csvs)} CSV(s) from {run_dir}")
    for csv_path in csvs:
        gene_list = csv_path.stem.removeprefix("tf_targets_")
        df = pd.read_csv(csv_path)
        if df.empty:
            print(f"  [skip] {csv_path.name}: empty")
            continue
        if not detect_groups(df):
            print(f"  [skip] {csv_path.name}: no validated_* columns")
            continue
        plot_fn(
            df,
            fig_dir / f"tf_validation_network_{gene_list}.png",
            top_n=top_n,
            title=f"TF Regulatory Networks — {gene_list} (ATAC validation)",
        )


def run(validation_dir: str, plot_style: str = "connected", top_n: int = 8) -> None:
    root = Path(validation_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"validation dir not found: {root}")

    # A single run folder, or the tf_validation root holding many of them.
    run_dirs = ([root] if any(root.glob("tf_targets_*.csv"))
                else sorted(d for d in root.iterdir()
                            if d.is_dir() and any(d.glob("tf_targets_*.csv"))))
    if not run_dirs:
        raise FileNotFoundError(f"No tf_targets_*.csv files under {root}")

    for d in run_dirs:
        plot_validation_run(d, plot_style=plot_style, top_n=top_n)


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
    return p


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    run(args.validation_dir, plot_style=args.plot_style, top_n=args.top_n_tfs)
