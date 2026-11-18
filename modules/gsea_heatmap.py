"""
gsea_heatmap.py

Square heatmap of GSEA results from one or more batch summaries: a row per
gene list, a column per condition (column_condition_value), with conditions
grouped under their column (column_condition_title) along the bottom.

Each cell is a square coloured by NES (blue = negative, red = positive); its
size and marks give the significance, from the BH-corrected FDR q
(FDR_qval_BH) of the batch summary:
  large, black border, *   |NES| > --nes-cutoff and q < 0.05
  large, *                 q < 0.05
  large                    q < --trend-q (default 0.25, the usual GSEA
                           exploratory cutoff)
  small                    q >= --trend-q
  empty                    gene list not tested in that condition

Inputs:
  --batch-summaries  batch summary CSVs. If the same gene list x column x
                     condition is in more than one, append =Label to each
                     path (e.g. summary.csv=Cortex); the group names then get
                     the label, e.g. "CellClass (Cortex)".
  --gene-lists       gene lists to show, in this order: the file name as in
                     gene_list_file_name (with or without .csv), optionally
                     =Row label. Default: every gene list in the summaries.
  --columns          only these column_condition_title values, in this order.

Usage:
  python modules/gsea_heatmap.py \
    --batch-summaries results/enrichment_results/batch_summary_A.csv=Cortex \
                      results/enrichment_results/batch_summary_B.csv=All_Layers \
    --gene-lists microcephaly_genes "autism_sfari=Autism (SFARI)" \
    --output results/additional_analyses/gsea_heatmap/heatmap.png
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

CMAP = LinearSegmentedColormap.from_list("nes", ["#3a54a4", "#ffffff", "#c3272b"])
CELL = 0.38          # inches per cell
BIG, SMALL = 0.78, 0.46   # square side, fraction of the cell
GROUP_GAP = 0.5      # gap between column groups, in cells


def split_label(item: str) -> tuple[str, str | None]:
    if "=" in item:
        a, b = item.split("=", 1)
        return a.strip(), b.strip()
    return item.strip(), None


def stem(name: str) -> str:
    return Path(str(name)).stem if str(name).lower().endswith((".csv", ".txt", ".tsv")) else str(name)


def load(summaries: list[str]) -> pd.DataFrame:
    frames = []
    for item in summaries:
        path, label = split_label(item)
        df = pd.read_csv(path)
        df["_summary"] = label or Path(path).stem
        frames.append(df)
        print(f"  {path}: {len(df):,} rows, {df['gene_list_file_name'].nunique()} gene lists")
    df = pd.concat(frames, ignore_index=True)
    df["_list"] = df["gene_list_file_name"].map(stem)
    df["column_condition_title"] = df["column_condition_title"].astype(str)
    df["column_condition_value"] = df["column_condition_value"].astype(str)
    key = ["_list", "column_condition_title", "column_condition_value"]
    if df.duplicated(key).any():
        if df.duplicated(key + ["_summary"]).any():
            raise ValueError("A gene list x column x condition is twice in the same batch "
                             "summary; pass summaries that do not overlap.")
        print("  Same gene list x condition in several summaries: groups get the summary label")
        df["column_condition_title"] = df["column_condition_title"] + " (" + df["_summary"] + ")"
    return df


def plot(df: pd.DataFrame, rows: list[tuple[str, str]], groups: list[tuple[str, list[str]]],
         nes_cutoff: float, trend_q: float, vmax: float, title: str | None,
         out: Path) -> None:
    lookup = df.set_index(["_list", "column_condition_title", "column_condition_value"])

    # x position of every column, with a gap between groups
    xs, x = [], 0.0
    for g, conds in groups:
        xs.append([x + i for i in range(len(conds))])
        x += len(conds) + GROUP_GAP
    width = x - GROUP_GAP
    nrow = len(rows)

    fig, ax = plt.subplots(figsize=((width + 9) * CELL + 2.5, max(nrow, 6) * CELL + 2.4))
    norm = Normalize(-vmax, vmax)
    for r, (lst, _) in enumerate(rows):
        for (g, conds), gx in zip(groups, xs):
            for cond, cx in zip(conds, gx):
                ax.add_patch(Rectangle((cx, r), 1, 1, fill=False, ec="#9a9a9a", lw=0.6))
                if (lst, g, cond) not in lookup.index:
                    continue
                rec = lookup.loc[(lst, g, cond)]
                if isinstance(rec, pd.DataFrame):
                    rec = rec.iloc[0]
                nes, q = rec["NES"], rec["FDR_qval_BH"]
                if pd.isna(nes):
                    continue
                sig_q = pd.notna(q) and q < 0.05
                strong = sig_q and abs(nes) > nes_cutoff
                big = pd.notna(q) and q < trend_q
                side = BIG if big else SMALL
                off = (1 - side) / 2
                ax.add_patch(Rectangle((cx + off, r + off), side, side,
                                       fc=CMAP(norm(np.clip(nes, -vmax, vmax))),
                                       ec="black" if strong else "#9a9a9a",
                                       lw=1.6 if strong else 0.6))
                if sig_q:
                    ax.text(cx + 0.5, r + 0.62, "*", ha="center", va="center", fontsize=11)

    # legend to the right of the heatmap, in data units (1 = one cell)
    lx = width + 1.0
    ax.text(lx, -0.1, "NES", fontsize=9, va="bottom")
    nseg = 40
    for k in range(nseg):
        v = -vmax + (k + 0.5) * 2 * vmax / nseg
        ax.add_patch(Rectangle((lx + k * 3 / nseg, 0.0), 3 / nseg, 0.45,
                               fc=CMAP(norm(v)), ec="none", clip_on=False))
    ax.add_patch(Rectangle((lx, 0.0), 3, 0.45, fill=False, ec="#9a9a9a", lw=0.6, clip_on=False))
    for v, t in ((-vmax, f"{-vmax:g}"), (0, "0"), (vmax, f"{vmax:g}")):
        ax.text(lx + 1.5 + 1.5 * v / vmax, 0.55, t, ha="center", va="top", fontsize=8)
    legend = [(f"|NES| > {nes_cutoff:g} & FDR < 0.05", BIG, "black", 1.6, True),
              ("FDR < 0.05", BIG, "#9a9a9a", 0.6, True),
              (f"FDR < {trend_q:g}", BIG, "#9a9a9a", 0.6, False),
              (f"FDR \u2265 {trend_q:g}", SMALL, "#9a9a9a", 0.6, False)]
    for k, (txt, side, ec, lw, star) in enumerate(legend):
        y = 1.6 + k * 1.1
        off = (1 - side) / 2
        ax.add_patch(Rectangle((lx, y), 1, 1, fill=False, ec="#9a9a9a", lw=0.6, clip_on=False))
        ax.add_patch(Rectangle((lx + off, y + off), side, side, fc="#c8c8c8", ec=ec, lw=lw,
                               clip_on=False))
        if star:
            ax.text(lx + 0.5, y + 0.62, "*", ha="center", va="center", fontsize=11)
        ax.text(lx + 1.3, y + 0.5, txt, va="center", fontsize=8.5)

    ax.set_xlim(-0.1, width + 0.1)
    ax.set_ylim(nrow + 0.05, -0.05)
    ax.set_aspect("equal")
    ax.set_yticks([r + 0.5 for r in range(nrow)], [lab for _, lab in rows], fontsize=9)
    ticks = [c + 0.5 for gx in xs for c in gx]
    ax.set_xticks(ticks, [c for _, conds in groups for c in conds], rotation=90, fontsize=9)
    ax.tick_params(length=0)
    ax.tick_params(axis="x", pad=12)   # room for the group line
    for sp in ax.spines.values():
        sp.set_visible(False)
    # a line under each group, and the group name under its condition labels
    for (g, conds), gx in zip(groups, xs):
        ax.plot([gx[0] + 0.05, gx[-1] + 0.95], [nrow + 0.25] * 2, color="black", lw=0.8,
                clip_on=False)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    low = min(t.get_window_extent(renderer).y0 for t in ax.get_xticklabels())
    y_axes = ax.transAxes.inverted().transform((0, low))[1]
    for (g, conds), gx in zip(groups, xs):
        ax.annotate(g, xy=((gx[0] + gx[-1] + 1) / 2, y_axes),
                    xycoords=("data", "axes fraction"), xytext=(0, -6),
                    textcoords="offset points", ha="center", va="top", fontsize=11,
                    fontweight="bold", annotation_clip=False)
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold", loc="left", pad=12)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out} (+ .pdf)")


def main():
    ap = argparse.ArgumentParser(description="GSEA NES heatmap from batch summaries.")
    ap.add_argument("--batch-summaries", nargs="+", required=True,
                    help="Batch summary CSVs, each optionally =Label.")
    ap.add_argument("--gene-lists", nargs="+", default=None,
                    help="Gene lists (file name, optional =Row label), in row order. "
                         "Default: all.")
    ap.add_argument("--gene-list-file", default=None,
                    help="Or a text file with one gene list per line (same NAME[=Label] form).")
    ap.add_argument("--columns", nargs="+", default=None,
                    help="Only these column_condition_title values, in this order.")
    ap.add_argument("--nes-cutoff", type=float, default=1.5,
                    help="|NES| for the black border (with FDR q < 0.05). Default 1.5.")
    ap.add_argument("--trend-q", type=float, default=0.25,
                    help="FDR q below which the square is large (no star). Default 0.25.")
    ap.add_argument("--vmax", type=float, default=None,
                    help="Colour scale limit (symmetric). Default: max |NES|, rounded up.")
    ap.add_argument("--no-n", action="store_true", help="Do not add (n=genes) to row labels.")
    ap.add_argument("--title", default=None)
    ap.add_argument("--output", required=True, help="Output PNG (a PDF is saved next to it).")
    args = ap.parse_args()

    print("Loading batch summaries")
    df = load(args.batch_summaries)

    items = list(args.gene_lists or [])
    if args.gene_list_file:
        items += [l.strip() for l in Path(args.gene_list_file).read_text().splitlines()
                  if l.strip() and not l.startswith("#")]
    if not items:
        items = list(dict.fromkeys(df["_list"]))
    rows, missing = [], []
    for it in items:
        name, label = split_label(it)
        name = stem(name)
        sub = df[df["_list"] == name]
        if sub.empty:
            missing.append(name)
            continue
        if label is None:
            desc = sub.get("gene_list_description")
            label = (str(desc.dropna().iloc[0]) if desc is not None and desc.notna().any()
                     else name)
        n = pd.to_numeric(sub["num_genes_in_gene_list"], errors="coerce").dropna()
        if not args.no_n and len(n):
            label = f"{label}\n(n={int(n.iloc[0]):,})"
        rows.append((name, label))
    if missing:
        print(f"  [warning] not in any batch summary, left out: {missing}")
    if not rows:
        raise ValueError("None of the gene lists are in the batch summaries. Names there: "
                         f"{sorted(df['_list'].unique())[:30]}")

    df = df[df["_list"].isin([r[0] for r in rows])]
    titles = list(dict.fromkeys(df["column_condition_title"]))
    if args.columns:
        wanted = [c for c in args.columns if c in titles]
        # with summary labels the titles are "Col (Label)": match on the prefix too
        wanted += [t for c in args.columns for t in titles
                   if t.startswith(c + " (") and t not in wanted]
        if not wanted:
            raise ValueError(f"--columns not found. Columns: {titles}")
        titles = wanted
    groups = [(t, list(dict.fromkeys(df.loc[df["column_condition_title"] == t,
                                            "column_condition_value"])))
              for t in titles]
    vmax = args.vmax or float(np.ceil(np.nanmax(np.abs(df["NES"].astype(float)))))
    print(f"  {len(rows)} gene lists x {sum(len(c) for _, c in groups)} conditions "
          f"in {len(groups)} columns")
    plot(df, rows, groups, args.nes_cutoff, args.trend_q, vmax, args.title, Path(args.output))


if __name__ == "__main__":
    main()
