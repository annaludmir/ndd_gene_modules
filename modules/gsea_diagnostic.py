"""
gsea_diagnostic.py

Why does a gene list get a given NES in one condition and not another? For a
gene list x batch summary, and each condition of one column (default
CellClass: Neuron, Neuroblast), recompute the GSEA statistics on the same
GES ranking GSEA used and show what drives the NES:

  NES = ES / mean ES of random gene sets of the same size (same sign).

So a high ES can still give a low NES when random gene sets also reach a
high ES on that ranking, e.g. when a few genes at the top have very large
GES scores (weighted GSEA, p = 1: every hit moves the curve in proportion
to its score).

Per gene list x summary, one figure with a row per condition:
  1. Ranked GES scores (log scale) with the gene list's genes marked; the
     share of the total score held by the top 1% / 10% of genes.
  2. Running ES, weighted (p = 1, as GSEA ran) and classic (p = 0).
  3. Null ES of random gene sets of the same size, weighted, with the
     observed ES: mean null ES (the NES denominator), NES, nominal p.
  4. The same for p = 0 (does the result depend on the weighting?).
The batch summary's NES / FDR (gseapy) is in the row title for comparison;
the recomputed NES should be close to it (different random sets).

Inputs: the batch summaries (as for gsea_heatmap / leading_edge_overlap),
with either --cases (CSV with gene_list, summary columns, e.g.
config_files/leading_edge_cases_neuron_neuroblast.csv) or --gene-lists and
--summaries (every combination).
Ranking: <GES run>/data/ges_spec_<column>_<condition>.csv of the summary's
GES run, filtered to ges_score > min_ges_score_threshold of the enrichment
config. Gene list: gene_list_path of that enrichment config ('gene' column),
else --gene-list-folder/<gene list>.csv.

Outputs (results/additional_analyses/gsea_diagnostic/<label>/):
  <gene_list>__<summary>.png
  gsea_diagnostic.csv   per gene list x summary x condition x weight: genes
                        ranked / in the set, score concentration, ES, mean
                        null ES, NES, nominal p, batch summary NES / FDR

Usage:
  python modules/gsea_diagnostic.py \
    --batch-summaries "results/.../batch_summary_all_layers.csv=All V3" \
    --gene-lists ASD_SFARI_syndromic --summaries "All V3" \
    --conditions Neuron Neuroblast --label asd_syndromic
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from leading_edge_overlap import load_summaries, ranking, stem  # noqa: E402

DEFAULT_OUTPUT_DIR = Path("results/additional_analyses/gsea_diagnostic")
COL_W1, COL_W0, COL_HIT = "#1b7f3a", "#7b3294", "#d62728"


# ---------------------------------------------------------------- inputs
def gene_list_genes(row: pd.Series, gl: str, folder: str | None) -> list[str]:
    """Genes of the gene list: the enrichment config's gene_list_path, else
    <folder>/<gl>.csv."""
    cands = []
    try:
        cfg = yaml.safe_load(open(row["enrichment_config_path"]))
        p = Path(str(cfg.get("gene_list_path", "")))
        root = cfg.get("ndd_gene_modules_folder_root")
        cands += [p] + ([Path(root) / p] if root and not p.is_absolute() else [])
    except Exception:
        pass
    if folder:
        cands.append(Path(folder) / f"{gl}.csv")
    for p in cands:
        if p.is_file() and p.stem == gl:
            g = pd.read_csv(p)
            col = "gene" if "gene" in g.columns else g.columns[0]
            return g[col].dropna().astype(str).str.strip().unique().tolist()
    raise FileNotFoundError(f"gene list {gl} not found (tried {[str(c) for c in cands]})")


# ---------------------------------------------------------------- GSEA statistics
def es_curves(scores: np.ndarray, hits: np.ndarray, weight: float) -> np.ndarray:
    """Running ES for each row of the boolean hit matrix `hits` (sets x genes)."""
    n = scores.size
    w = np.abs(scores) ** weight
    hw = hits * w
    nh = hits.sum(1, keepdims=True)
    p_hit = hw / hw.sum(1, keepdims=True)
    p_miss = (~hits) / (n - nh)
    return np.cumsum(p_hit - p_miss, axis=1)


def peak(curves: np.ndarray) -> np.ndarray:
    """ES: the running ES furthest from 0 (signed), per row."""
    i = np.abs(curves).argmax(1)
    return curves[np.arange(len(curves)), i]


def null_es(scores: np.ndarray, k: int, weight: float, n_perm: int, seed: int,
            chunk: int = 200) -> np.ndarray:
    """ES of n_perm random gene sets of size k on the same ranking (gene-set
    permutation, as gseapy prerank)."""
    rng = np.random.default_rng(seed)
    n, out = scores.size, []
    for start in range(0, n_perm, chunk):
        m = min(chunk, n_perm - start)
        hits = np.zeros((m, n), dtype=bool)
        idx = np.argsort(rng.random((m, n)), axis=1)[:, :k]
        np.put_along_axis(hits, idx, True, axis=1)
        out.append(peak(es_curves(scores, hits, weight)))
    return np.concatenate(out)


def nes_p(es: float, null: np.ndarray) -> tuple[float, float, float]:
    """(mean null ES of the same sign, NES, nominal p) as in GSEA."""
    same = null[null >= 0] if es >= 0 else null[null < 0]
    if same.size == 0:
        return np.nan, np.nan, np.nan
    mean = float(np.abs(same).mean())
    p = float(((np.abs(same) >= abs(es)).sum() + 0) / same.size)
    return mean, float(np.sign(es) * abs(es) / mean), p


def concentration(scores: np.ndarray) -> dict:
    s = np.sort(np.abs(scores))[::-1]
    tot = s.sum()
    top = lambda f: float(s[:max(1, int(round(f * s.size)))].sum() / tot)
    return {"top1pct_score_share": top(0.01), "top10pct_score_share": top(0.10),
            "max_score": float(s[0]), "median_score": float(np.median(s))}


# ---------------------------------------------------------------- plot
def fmt_p(p, n_perm) -> str:
    """Nominal p; 0 means no random set was as extreme: < 1/n_perm."""
    return f"< {1 / n_perm:.0e}" if p == 0 else f"= {p:.3g}"


def null_panel(ax, null, es, mean, nes, p, color, weight):
    same = null[null >= 0] if es >= 0 else null[null < 0]
    ax.hist(same, bins=40, color=color, alpha=0.45)
    ax.axvline(es, color="black", lw=2, label=f"observed ES {es:.3f}")
    ax.axvline(np.sign(es) * mean, color=color, lw=1.6, ls="--",
               label=f"mean null ES {mean:.3f}")
    ax.set_title(f"Null, p = {weight:g}: NES = {es:.3f} / {mean:.3f} = {nes:.2f}, "
                 f"nominal p {fmt_p(p, null.size)}", fontsize=9.5)
    ax.set_xlabel("ES of random gene sets (same size)")
    ax.set_ylabel("Random sets")
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)


def figure(gl, summary, results, n_perm, path):
    nrow = len(results)
    fig, axes = plt.subplots(nrow, 4, figsize=(24, 4.6 * nrow), squeeze=False,
                             gridspec_kw={"width_ratios": [1.1, 1.2, 1, 1]})
    for r, res in enumerate(results):
        cond, rank, pos = res["condition"], res["rank"], res["pos"]
        ax = axes[r, 0]
        x = np.arange(1, len(rank) + 1)
        ax.plot(x, rank.to_numpy(), color="#444", lw=1.4)
        ax.scatter(pos + 1, rank.to_numpy()[pos], s=8, color=COL_HIT, zorder=3,
                   label=f"{gl} genes (n={len(pos)})")
        ax.set_yscale("log")
        c = res["conc"]
        ax.set_title(f"{cond}: ranked GES scores ({len(rank):,} genes)\n"
                     f"top 1% of genes hold {c['top1pct_score_share']:.0%} of the total score, "
                     f"top 10% {c['top10pct_score_share']:.0%}", fontsize=9.5)
        ax.set_xlabel("Rank"); ax.set_ylabel("GES score (log)")
        ax.legend(fontsize=8, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)

        ax = axes[r, 1]
        for wt, col in ((1, COL_W1), (0, COL_W0)):
            st = res[wt]
            ax.plot(x, st["curve"], color=col, lw=2,
                    label=f"p = {wt} ({'weighted, as GSEA ran' if wt else 'classic'}): "
                          f"ES {st['es']:.3f}, NES {st['nes']:.2f}")
        ax.vlines(pos + 1, -0.02, 0.0, transform=ax.get_xaxis_transform(),
                  color=COL_HIT, lw=0.6, clip_on=False)
        ax.axhline(0, color="#999", ls="--", lw=0.8)
        b = res["batch"]
        ax.set_title(f"{cond}: running ES ({gl}, {summary})\nbatch summary (gseapy): "
                     f"NES {b[0]:.3f}, FDR {b[1]:.3g}", fontsize=9.5, fontweight="bold")
        ax.set_xlabel("Rank"); ax.set_ylabel("Enrichment score")
        ax.legend(fontsize=8, frameon=False, loc="best")
        ax.grid(alpha=0.3)
        ax.spines[["top", "right"]].set_visible(False)

        for k, (wt, col) in enumerate(((1, COL_W1), (0, COL_W0))):
            st = res[wt]
            null_panel(axes[r, 2 + k], st["null"], st["es"], st["mean"], st["nes"],
                       st["p"], col, wt)
    fig.suptitle(f"{gl}, {summary}: what drives the NES ({n_perm} random gene sets "
                 f"per condition)", fontsize=13, fontweight="bold", y=1.0)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved: {path}")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="GSEA NES diagnostic per condition.")
    ap.add_argument("--batch-summaries", nargs="+", required=True,
                    help="Batch summary CSVs, each =Label.")
    ap.add_argument("--cases", default=None,
                    help="CSV with gene_list, summary columns (other columns ignored).")
    ap.add_argument("--gene-lists", nargs="+", default=None)
    ap.add_argument("--summaries", nargs="+", default=None,
                    help="Summary labels for --gene-lists. Default: all.")
    ap.add_argument("--column", default="CellClass")
    ap.add_argument("--conditions", nargs="+", default=["Neuron", "Neuroblast"])
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=6)
    ap.add_argument("--default-threshold", type=float, default=0.0,
                    help="GES threshold when the enrichment config can't be read.")
    ap.add_argument("--gene-list-folder", default="data/genes/final_genes_to_run_on",
                    help="Where to find <gene list>.csv if the config's path is missing.")
    ap.add_argument("--label", required=True)
    ap.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    args = ap.parse_args()

    print("Batch summaries")
    df = load_summaries(args.batch_summaries)
    df = df[df["column_condition_title"].astype(str) == args.column]
    if args.cases:
        c = pd.read_csv(args.cases, skipinitialspace=True).astype(str)
        pairs = [(stem(g.strip()), s.strip()) for g, s in zip(c["gene_list"], c["summary"])]
    elif args.gene_lists:
        summ = args.summaries or list(dict.fromkeys(df["_summary"]))
        pairs = [(stem(g), s) for g in args.gene_lists for s in summ]
    else:
        ap.error("pass --cases or --gene-lists")
    pairs = list(dict.fromkeys(pairs))
    out_dir = Path(args.output_dir) / args.label
    rows = []
    for gl, summ in pairs:
        print(f"\n{gl} / {summ}")
        sub = df[(df["_list"] == gl) & (df["_summary"] == summ)]
        results = []
        for cond in args.conditions:
            r = sub[sub["column_condition_value"].astype(str) == cond]
            if r.empty:
                print(f"  [skip] {cond}: not in the batch summary")
                continue
            rec = r.iloc[0]
            try:
                rank = ranking(rec, args.column, cond, args.default_threshold)
                genes = gene_list_genes(rec, gl, args.gene_list_folder)
            except (FileNotFoundError, KeyError) as e:
                print(f"  [skip] {cond}: {e}")
                continue
            hit = rank.index.isin(genes)
            pos, k = np.flatnonzero(hit), int(hit.sum())
            if k < 2:
                print(f"  [skip] {cond}: {k} genes of the list in the ranking")
                continue
            scores = rank.to_numpy()
            res = {"condition": cond, "rank": rank, "pos": pos,
                   "conc": concentration(scores),
                   "batch": (rec.get("NES", np.nan), rec.get("FDR_qval_BH", np.nan))}
            for wt in (1, 0):
                curve = es_curves(scores, hit[None, :], wt)[0]
                es = float(peak(curve[None, :])[0])
                null = null_es(scores, k, wt, args.n_perm, args.seed)
                mean, nes, p = nes_p(es, null)
                res[wt] = {"curve": curve, "es": es, "null": null, "mean": mean,
                           "nes": nes, "p": p}
                rows.append({"gene_list": gl, "summary": summ, "condition": cond,
                             "weight": wt, "n_ranked": len(rank), "n_in_set": k,
                             **res["conc"], "ES": es, "mean_null_ES": mean, "NES": nes,
                             "nominal_p": p, "batch_NES": res["batch"][0],
                             "batch_FDR_qval_BH": res["batch"][1]})
                print(f"  {cond:12s} p={wt}: ES {es:.3f}, mean null {mean:.3f}, "
                      f"NES {nes:.2f}, p {fmt_p(p, args.n_perm)}"
                      + (f"   (batch summary NES {res['batch'][0]:.3f})" if wt else ""))
            results.append(res)
        if results:
            safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in f"{gl}__{summ}")
            figure(gl, summ, results, args.n_perm, out_dir / f"{safe}.png")
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_dir / "gsea_diagnostic.csv", index=False)
    print(f"\nOutputs: {out_dir}")


if __name__ == "__main__":
    main()
