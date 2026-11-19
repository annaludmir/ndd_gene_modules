"""
leading_edge_overlap.py

Compare the GSEA leading-edge genes of one gene list in two conditions of the
same column (default CellClass: Neuron vs Neuroblast), per batch summary
(e.g. All / Cortex x V3 / V2).

For each case (gene list x batch summary, from --cases):
  - Venn diagram of the two leading edges: shared genes and genes in only one.
    In the Venn-only figures, below it the top --n-top genes of each part
    (A only / both / B only) by GES score in that condition (both: mean of
    the two conditions), score(s) in brackets.
  - With --es-plots for the case's category (default: the "only" categories):
    running enrichment score (GSEA, weighted as in gseapy prerank) of the
    shared genes and of each condition's exclusive genes, along that
    condition's GES ranking, with the hit positions below. Panel A: ranked
    in condition A (shared + A-only), panel B: ranked in condition B
    (shared + B-only), panel C: the Venn diagram.

The ranking of a condition is the GES file GSEA used:
<GES run>/data/ges_spec_<column>_<condition>.csv, the GES run being the
folder of the batch summary's ges_config_path, filtered to ges_score >
min_ges_score_threshold of the run's enrichment config (as in
search_enrichment_gsea.py).

Cases file (CSV): gene_list, summary, category
  gene_list  gene_list_file_name without .csv (e.g. ASD_HC)
  summary    a batch summary label from --batch-summaries (e.g. All V3)
  category   free text, used for grouping, e.g. both / neuron_only /
             neuroblast_only

Outputs (results/additional_analyses/leading_edge_overlap/<label>/):
  <category>/<gene_list>__<summary>.png   Venn (+ ES panels)
  leading_edge_overlap.csv                per case: NES / FDR of both
                                          conditions, counts and the genes
                                          in shared / A-only / B-only

Usage:
  python modules/leading_edge_overlap.py \
    --batch-summaries "results/.../batch_summary_all_layers.csv=All V3" ... \
    --cases config_files/leading_edge_cases_neuron_neuroblast.csv \
    --label autism_neuron_vs_neuroblast
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.patches import Circle

DEFAULT_OUTPUT_DIR = Path("results/additional_analyses/leading_edge_overlap")
COLORS = {"shared": "#d07a2a", "A": "#2f62e0", "B": "#2e9461"}


def stem(name) -> str:
    s = str(name)
    return s[:-4] if s.lower().endswith(".csv") else s


def split_genes(s) -> list[str]:
    if pd.isna(s) or not str(s).strip():
        return []
    return [g.strip() for g in str(s).replace(",", ";").split(";") if g.strip()]


# ---------------------------------------------------------------- inputs
def load_summaries(items: list[str]) -> pd.DataFrame:
    frames = []
    for it in items:
        path, label = (it.split("=", 1) + [None])[:2] if "=" in it else (it, None)
        path, label = path.strip(), (label or Path(path).stem).strip()
        df = pd.read_csv(path)
        df["_summary"] = label
        frames.append(df)
        print(f"  {label}: {path} ({len(df):,} rows)")
    df = pd.concat(frames, ignore_index=True)
    df["_list"] = df["gene_list_file_name"].map(stem)
    return df


def ges_run_folder(ges_config_path) -> Path | None:
    """The GES run folder holding data/ges_spec_*.csv, from its config path."""
    if pd.isna(ges_config_path):
        return None
    p = Path(str(ges_config_path))
    for d in [p] + list(p.parents):
        if (d / "data").is_dir() and any((d / "data").glob("ges_spec_*.csv")):
            return d
    return None


def gsea_threshold(enrichment_config_path):
    try:
        cfg = yaml.safe_load(open(enrichment_config_path))
        return cfg.get("gsea", {}).get("min_ges_score_threshold")
    except Exception:
        return None


_rank_cache: dict = {}


def ranking(row: pd.Series, column: str, condition: str, default_threshold):
    """GES ranking (gene -> score, descending) GSEA used for this condition."""
    folder = ges_run_folder(row.get("ges_config_path"))
    if folder is None:
        raise FileNotFoundError(f"GES run folder not found from ges_config_path "
                                f"{row.get('ges_config_path')}")
    thr = gsea_threshold(row.get("enrichment_config_path"))
    thr = default_threshold if thr is None else thr
    key = (str(folder), column, condition, thr)
    if key not in _rank_cache:
        f = folder / "data" / f"ges_spec_{column}_{condition}.csv"
        g = pd.read_csv(f)
        if thr is not None:
            g = g[g["ges_score"] > thr]
        _rank_cache[key] = (g.set_index("gene")["ges_score"].astype(float)
                             .sort_values(ascending=False))
    return _rank_cache[key]


# ---------------------------------------------------------------- GSEA curve
def running_es(rank: pd.Series, genes: list[str], weight: float = 1.0):
    """GSEA running enrichment score (Subramanian 2005, as gseapy prerank)."""
    hit = rank.index.isin(genes)
    n, nh = len(rank), int(hit.sum())
    if nh == 0 or nh == n:
        return None, None, None
    w = np.abs(rank.to_numpy()) ** weight
    p_hit = np.where(hit, w, 0) / w[hit].sum()
    p_miss = np.where(~hit, 1.0, 0) / (n - nh)
    es = np.cumsum(p_hit - p_miss)
    return es, np.flatnonzero(hit), nh


# ---------------------------------------------------------------- plots
def venn(ax, a: set, b: set, name_a: str, name_b: str, title: str) -> None:
    r, d = 1.0, 0.75
    for x, col in ((-d / 2, COLORS["A"]), (d / 2, COLORS["B"])):
        ax.add_patch(Circle((x, 0), r, fc=col, alpha=0.35, ec=col, lw=1.5))
    ax.text(-d / 2 - 0.55, 0, f"{len(a - b)}", ha="center", va="center",
            fontsize=16, fontweight="bold")
    ax.text(0, 0, f"{len(a & b)}", ha="center", va="center", fontsize=20,
            fontweight="bold")
    ax.text(d / 2 + 0.55, 0, f"{len(b - a)}", ha="center", va="center",
            fontsize=16, fontweight="bold")
    ax.text(-d / 2 - 0.35, r + 0.12, f"{name_a}\n(total: {len(a)})", ha="center",
            va="bottom", color=COLORS["A"], fontsize=10, fontweight="bold")
    ax.text(d / 2 + 0.35, r + 0.12, f"{name_b}\n(total: {len(b)})", ha="center",
            va="bottom", color=COLORS["B"], fontsize=10, fontweight="bold")
    ax.set_xlim(-2.0, 2.0); ax.set_ylim(-1.3, 1.9)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=11, fontweight="bold")


def es_panel(ax_es, ax_hit, rank, sets: list[tuple[str, list[str], str]], title: str):
    xmax = 0
    for name, genes, col in sets:
        es, pos, nh = running_es(rank, genes)
        if es is None:
            continue
        peak = es[np.argmax(np.abs(es))]
        ax_es.plot(np.arange(len(es)), es, color=col, lw=2.2,
                   label=f"{name} (n={nh}) (peak={peak:.2f})")
        ax_hit.vlines(pos, 0, 1, color=col, lw=1.2)
        xmax = max(xmax, pos.max())
    xmax = int(xmax * 1.03) + 10 if xmax else len(rank)
    for ax in (ax_es, ax_hit):
        ax.set_xlim(0, xmax)
        ax.spines[["top", "right"]].set_visible(False)
    ax_es.axhline(0, color="#888", lw=0.8, ls="--")
    ax_es.grid(True, ls="--", alpha=0.4)
    ax_es.set_ylabel("Enrichment Score (ES)", fontweight="bold")
    ax_es.set_title(title, fontsize=11, fontweight="bold", loc="left")
    ax_es.legend(loc="upper left", fontsize=9)
    ax_es.tick_params(labelbottom=False)
    ax_hit.set_yticks([]); ax_hit.set_ylim(0, 1)
    ax_hit.set_ylabel("Hits", fontweight="bold")
    ax_hit.set_xlabel(f"Rank in ordered gene list (of {len(rank):,} genes)",
                      fontweight="bold")


def significant(rec) -> bool:
    """The batch summary's is_significant (|NES| > 1.5 and FDR <= 0.05)."""
    v = rec.get("is_significant")
    if pd.notna(v) and str(v).strip() != "":
        return str(v).strip().lower() in ("true", "1", "yes")
    nes, q = rec.get("NES"), rec.get("FDR_qval_BH")
    return bool(pd.notna(nes) and abs(nes) > 1.5 and pd.notna(q) and q <= 0.05)


def top_genes(genes, ranks: list, lead_orders: list, n: int) -> list[str]:
    """Top n of `genes`: by GES score (mean over the given rankings, missing = lowest),
    or, without rankings, by mean position in the leading-edge lists."""
    genes = list(genes)
    if not genes:
        return []
    if all(r is not None for r in ranks):
        score = {g: np.mean([r.get(g, r.min()) for r in ranks]) for g in genes}
        order = sorted(genes, key=lambda g: -score[g])
        fmt = lambda g: f"{g} ({' / '.join(f'{r.get(g, np.nan):.2f}' for r in ranks)})"
    else:
        pos = [{g: i for i, g in enumerate(o)} for o in lead_orders]
        order = sorted(genes, key=lambda g: np.mean([q.get(g, len(q)) for q in pos]))
        fmt = str
    return [fmt(g) for g in order[:n]]


def gene_table(ax, cols: list[tuple[str, list[str], str]], note: str) -> None:
    ax.axis("off")
    for k, (head, genes, color) in enumerate(cols):
        x = (k + 0.5) / len(cols)
        ax.text(x, 1.0, head, ha="center", va="top", fontsize=9.5, fontweight="bold",
                color=color, transform=ax.transAxes)
        ax.text(x, 0.86, "\n".join(genes) if genes else "-", ha="center", va="top",
                fontsize=8.5, family="monospace", transform=ax.transAxes, linespacing=1.35)
    ax.text(0.5, -0.04, note, ha="center", va="top", fontsize=7.5, color="#666",
            transform=ax.transAxes)


def case_figure(gl, summary, a, b, cond_a, cond_b, rank_a, rank_b, with_es, path,
                lead_a=(), lead_b=(), n_top: int = 10, stats=None):
    """stats: {condition: (NES, FDR, significant)} for the ES figure titles."""
    shared, a_only, b_only = sorted(a & b), sorted(a - b), sorted(b - a)
    stats = stats or {}
    if with_es:
        fig = plt.figure(figsize=(17, 6))
        gs = fig.add_gridspec(2, 3, height_ratios=[3, 1], width_ratios=[1, 1, 0.8],
                              hspace=0.08, wspace=0.25)
        for k, (cond, rank, own, key) in enumerate(((cond_a, rank_a, a_only, "A"),
                                                    (cond_b, rank_b, b_only, "B"))):
            title = f"{'AB'[k]}. Ranked within {cond} ({gl}, {summary})"
            if cond in stats:
                nes, fdr, sig = stats[cond]
                title += (f"\n{'SIGNIFICANT' if sig else 'not significant'} in {cond}: "
                          f"NES {nes:.2f}, FDR {fdr:.3g}")
            ax_es = fig.add_subplot(gs[0, k])
            es_panel(ax_es, fig.add_subplot(gs[1, k]), rank,
                     [("Shared", shared, COLORS["shared"]),
                      (f"{cond} only", own, COLORS[key])], title)
            if cond in stats and stats[cond][2]:
                ax_es.title.set_color("#b2182b")
        venn(fig.add_subplot(gs[:, 2]), a, b, cond_a, cond_b,
             f"Leading-edge overlap: {cond_a} vs {cond_b}\n({gl}, {summary})")
        sig_in = [c for c in (cond_a, cond_b) if c in stats and stats[c][2]]
        fig.suptitle(f"{gl}, {summary}: enrichment significant in "
                     f"{' and '.join(sig_in) if sig_in else 'neither'}",
                     fontsize=14, fontweight="bold", y=1.04)
    else:
        # Venn, and below it the top leading genes of each part
        fig = plt.figure(figsize=(8, 7.6))
        gs = fig.add_gridspec(2, 1, height_ratios=[4.5, 3.0], hspace=0.05)
        venn(fig.add_subplot(gs[0]), a, b, cond_a, cond_b,
             f"Leading-edge overlap: {cond_a} vs {cond_b}\n({gl}, {summary})")
        if n_top:
            have = rank_a is not None and rank_b is not None
            cols = [(f"{cond_a} only", top_genes(a - b, [rank_a], [lead_a], n_top), COLORS["A"]),
                    ("Both", top_genes(a & b, [rank_a, rank_b], [lead_a, lead_b], n_top),
                     COLORS["shared"]),
                    (f"{cond_b} only", top_genes(b - a, [rank_b], [lead_b], n_top), COLORS["B"])]
            note = (f"Top {n_top} by GES score (in brackets; Both: {cond_a} / {cond_b}, "
                    f"ordered by their mean)" if have else
                    f"Top {n_top} by position in the leading edge (GES rankings not found)")
            gene_table(fig.add_subplot(gs[1]), cols, note)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description="Leading-edge overlap of two conditions.")
    p.add_argument("--batch-summaries", nargs="+", required=True,
                   help="Batch summary CSVs, each =Label (e.g. path=All V3).")
    p.add_argument("--cases", required=True,
                   help="CSV with columns gene_list, summary, category.")
    p.add_argument("--column", default="CellClass")
    p.add_argument("--condition-a", default="Neuron")
    p.add_argument("--condition-b", default="Neuroblast")
    p.add_argument("--es-plots", nargs="*", default=None,
                   help="Categories that get the ES panels. Default: every category "
                        "with 'only' in its name. Pass 'all' for every case.")
    p.add_argument("--n-top", type=int, default=10,
                   help="Venn-only figures: top genes listed per part (0 = none).")
    p.add_argument("--default-threshold", type=float, default=0.0,
                   help="GES threshold when the run's enrichment config can't be read.")
    p.add_argument("--label", required=True)
    p.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    args = p.parse_args()

    out_dir = Path(args.output_dir) / args.label
    out_dir.mkdir(parents=True, exist_ok=True)
    print("Batch summaries")
    df = load_summaries(args.batch_summaries)
    df = df[df["column_condition_title"].astype(str) == args.column]
    cases = pd.read_csv(args.cases, skipinitialspace=True)
    cases = cases.apply(lambda c: c.astype(str).str.strip())
    print(f"  {len(cases)} cases from {args.cases}")

    rows = []
    for c in cases.itertuples():
        gl, summ, cat = stem(c.gene_list), c.summary, c.category
        sub = df[(df["_list"] == gl) & (df["_summary"] == summ)]
        recs = {}
        for cond in (args.condition_a, args.condition_b):
            r = sub[sub["column_condition_value"].astype(str) == cond]
            recs[cond] = r.iloc[0] if len(r) else None
        if any(v is None for v in recs.values()):
            print(f"  [skip] {gl} / {summ}: missing "
                  f"{[k for k, v in recs.items() if v is None]} in the batch summary")
            continue
        ra, rb = recs[args.condition_a], recs[args.condition_b]
        a, b = set(split_genes(ra["lead_genes"])), set(split_genes(rb["lead_genes"]))
        with_es = (args.es_plots is not None and ("all" in args.es_plots or cat in args.es_plots)) \
            or (args.es_plots is None and "only" in cat.lower())
        rank_a = rank_b = None
        try:
            rank_a = ranking(ra, args.column, args.condition_a, args.default_threshold)
            rank_b = ranking(rb, args.column, args.condition_b, args.default_threshold)
        except (FileNotFoundError, KeyError) as e:
            rank_a = rank_b = None
            print(f"  [note] {gl} / {summ}: GES rankings not found ({e})"
                  + ("; no ES panels" if with_es else ""))
            with_es = False
        safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in f"{gl}__{summ}")
        path = out_dir / cat / f"{safe}.png"
        case_figure(gl, summ, a, b, args.condition_a, args.condition_b, rank_a, rank_b,
                    with_es, path, split_genes(ra["lead_genes"]),
                    split_genes(rb["lead_genes"]), args.n_top,
                    {args.condition_a: (ra["NES"], ra["FDR_qval_BH"], significant(ra)),
                     args.condition_b: (rb["NES"], rb["FDR_qval_BH"], significant(rb))})
        print(f"  {cat:16s} {gl:22s} {summ:10s} shared {len(a & b):3d} | "
              f"{args.condition_a}-only {len(a - b):3d} | {args.condition_b}-only "
              f"{len(b - a):3d}{'  (+ES)' if with_es else ''}")
        A, B = args.condition_a, args.condition_b
        rows.append({
            "category": cat, "gene_list": gl, "summary": summ,
            f"NES_{A}": ra["NES"], f"FDR_{A}": ra["FDR_qval_BH"],
            f"NES_{B}": rb["NES"], f"FDR_{B}": rb["FDR_qval_BH"],
            f"n_lead_{A}": len(a), f"n_lead_{B}": len(b),
            "n_shared": len(a & b), f"n_{A}_only": len(a - b), f"n_{B}_only": len(b - a),
            "jaccard": len(a & b) / len(a | b) if a | b else np.nan,
            "shared": ";".join(sorted(a & b)), f"{A}_only": ";".join(sorted(a - b)),
            f"{B}_only": ";".join(sorted(b - a)), "figure": str(path),
        })
    pd.DataFrame(rows).to_csv(out_dir / "leading_edge_overlap.csv", index=False)
    print(f"\nOutputs: {out_dir}")


if __name__ == "__main__":
    main()
