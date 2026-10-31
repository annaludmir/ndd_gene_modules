import marimo

__generated_with = "0.17.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import re
    import h5py
    import numpy as np
    import pandas as pd
    try:
        from anndata.io import read_elem            # anndata >= 0.11
    except ImportError:
        from anndata.experimental import read_elem  # older anndata

    CORTEX_PATH = "/miridan-data/annaludmir/ndd_gene_modules/data/Cortex_EMX1_louvain3_passedQC_PostM_rev1.h5ad"
    HUMAN_PATH = "/miridan-data/annaludmir/ndd_gene_modules/data/human_dev_without_week_5_cc_annotated_fraction.h5ad"

    def read_obs(path):
        with h5py.File(path, "r") as f:
            return pd.DataFrame(read_elem(f["obs"]))

    obs_cx = read_obs(CORTEX_PATH)
    obs_hd = read_obs(HUMAN_PATH)
    obs_cx.shape, obs_hd.shape
    return obs_cx, obs_hd, pd, re


@app.cell
def _(obs_cx, obs_hd, pd, re):
    def _describe(o, name):
        return pd.DataFrame({
            "file": name,
            "example cell IDs": [", ".join(map(str, o.index[:3]))],
            "phase-like columns": [", ".join(c for c in o.columns if re.search("phase|cycle", c, re.I))],
            "donor-like columns": [", ".join(c for c in o.columns if re.search("donor", c, re.I))],
        })

    pd.concat([_describe(obs_cx, "cortex"), _describe(obs_hd, "human_dev (new)")], ignore_index=True)
    return


@app.cell
def _(obs_cx, obs_hd, pd):
    PHASE_CX = "CellCyclePhase"   # phase column in the cortex file
    PHASE_HD = "CellCyclePhase"   # <- set to the NEW phase column in the human_dev file
    DONOR_CX = "Donor"
    DONOR_HD = "donor_id"
    CLASS_CX = "CellClass"

    def _norm_donor(s):
        return s.astype(str).str.replace(r"^(XHU|XDD):\d+:(\d+)$", r"\1:\2", regex=True)

    def _cell_key(o, donor_col):
        ids = pd.Series(o.index.astype(str), index=o.index)
        barcode = ids.str.extract(r"([ACGT]{16})", expand=False).fillna(ids)
        return _norm_donor(o[donor_col]) + "|" + barcode

    cx = pd.DataFrame({"key": _cell_key(obs_cx, DONOR_CX), "phase_cortex": obs_cx[PHASE_CX].astype(str),
                       "cell_class": obs_cx[CLASS_CX].astype(str), "age_cortex": obs_cx.get("Age")})
    hd = pd.DataFrame({"key": _cell_key(obs_hd, DONOR_HD), "phase_human_dev": obs_hd[PHASE_HD].astype(str),
                       "age_human_dev": obs_hd.get("Age")})

    dup_cx, dup_hd = cx.key.duplicated(keep=False), hd.key.duplicated(keep=False)
    matched = cx[~dup_cx].merge(hd[~dup_hd], on="key", how="inner")

    pd.Series({
        "cortex cells": len(cx),
        "cortex cells with an ambiguous key (dropped)": int(dup_cx.sum()),
        "human_dev cells with an ambiguous key (dropped)": int(dup_hd.sum()),
        "cells matched by donor + barcode": len(matched),
        "% of cortex cells matched": round(100 * len(matched) / len(cx), 1),
        "cell IDs identical as written (cross-check)": int(obs_cx.index.isin(obs_hd.index).sum()),
        "matched cells whose age disagrees": int((matched.age_cortex.astype(str) != matched.age_human_dev.astype(str)).sum()),
    })
    return (matched,)


@app.cell
def _(matched, pd):
    def _norm_phase(s):
        s = s.str.strip().replace({"PostM": "Post-M", "post-m": "Post-M", "nan": "NA"})
        return s

    matched["pc"] = _norm_phase(matched.phase_cortex)
    matched["ph"] = _norm_phase(matched.phase_human_dev)

    counts = pd.crosstab(matched.pc, matched.ph, margins=True)
    row_pct = (pd.crosstab(matched.pc, matched.ph, normalize="index") * 100).round(1)

    agreement = pd.Series({
        "overall agreement (%)": round(100 * (matched.pc == matched.ph).mean(), 1),
        **{f"agreement in {c} (%)": round(100 * (g.pc == g.ph).mean(), 1)
           for c, g in matched.groupby("cell_class")},
    })

    counts, row_pct, agreement
    return


@app.cell
def _():
    return


@app.cell
def _():
    import anndata as ad
    a = ad.read_h5ad("/miridan-data/annaludmir/ndd_gene_modules/data/Cortex_EMX1_louvain3_passedQC_PostM_rev1_neuron_type.h5ad", backed="r")
    return (a,)


@app.cell
def _(a, pd):
    n = a.obs[a.obs["NeuronType"].notna()].copy()
    n["unassigned"] = n["NeuronType"] == "Unassigned"
    print(n.groupby("Age")["unassigned"].mean().round(2))                   # falls with age → immaturity
    print(n.groupby("Chemistry")["unassigned"].mean().round(2))             # higher in v2 → dropout
    umi = "TotalUMI"              # whichever depth column exists
    print(n.groupby(pd.qcut(n[umi], 5))["unassigned"].mean().round(2))      # falls with depth → dropout
    return


@app.cell
def _():
    # 1. Select specific columns (using .copy() to avoid SettingWithCopyWarning)
    selected_cols = ["gene_list_file_name", "num_genes_in_gene_list", "column_condition_title",
        "column_condition_value", "NES", "FDR_qval_BH_reported", "lead_genes", "num_of_lead_genes"]

    return (selected_cols,)


@app.cell
def _(pd, selected_cols):
    summary_v3_all_layer_neuron__maturity = '/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/batch_summary_20261031_final_genes_to_run_on_all_layers_neuron_maturity.csv'
    df_v3_all_layer_neuron_maturity  = pd.read_csv(summary_v3_all_layer_neuron__maturity)

    filtered_df_v3_all_layer_neuron_maturity = df_v3_all_layer_neuron_maturity.query("is_significant == True and NES > 0")[selected_cols].copy()

    # 2. Add the constant column
    filtered_df_v3_all_layer_neuron_maturity['chemistry'] = 'v3'
    filtered_df_v3_all_layer_neuron_maturity['layer'] = 'all_layers'
    return (filtered_df_v3_all_layer_neuron_maturity,)


@app.cell
def _(pd, selected_cols):
    summary_v2_all_layer_neuron_maturity = '/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/batch_summary_20261031_final_genes_to_run_on_all_layers_neuron_maturity_v2.csv'
    df_v2_all_layer_neuron_maturity  = pd.read_csv(summary_v2_all_layer_neuron_maturity)

    filtered_df_v2_all_layer_neuron_maturity = df_v2_all_layer_neuron_maturity.query("is_significant == True and NES > 0")[selected_cols].copy()

    # 2. Add the constant column
    filtered_df_v2_all_layer_neuron_maturity['chemistry'] = 'v2'
    filtered_df_v2_all_layer_neuron_maturity['layer'] = 'all_layers'
    return (filtered_df_v2_all_layer_neuron_maturity,)


@app.cell
def _(filtered_df_v2_all_layer_neuron_maturity):
    filtered_df_v2_all_layer_neuron_maturity
    return


@app.cell
def _(pd, selected_cols):
    summary_v3_cortex_neuron_maturity = '/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/batch_summary_20261031_final_genes_to_run_on_cortex_neuron_maturity.csv'
    df_v3_cortex_neuron_maturity  = pd.read_csv(summary_v3_cortex_neuron_maturity)

    filtered_df_v3_cortex_neuron_maturity = df_v3_cortex_neuron_maturity.query("is_significant == True and NES > 0")[selected_cols].copy()

    # 2. Add the constant column
    filtered_df_v3_cortex_neuron_maturity['chemistry'] = 'v3'
    filtered_df_v3_cortex_neuron_maturity['layer'] = 'cortex'
    return


@app.cell
def _(pd, selected_cols):
    summary_v2_cortex_neuron_maturity = '/miridan-data/annaludmir/ndd_gene_modules/results/enrichment_results/batch_summary_20261031_final_genes_to_run_on_cortex_neuron_maturity_v2.csv'
    df_v2_cortex_neuron_maturity  = pd.read_csv(summary_v2_cortex_neuron_maturity)

    filtered_df_v2_cortex_neuron_maturity = df_v2_cortex_neuron_maturity.query("is_significant == True and NES > 0")[selected_cols].copy()

    # 2. Add the constant column
    filtered_df_v2_cortex_neuron_maturity['chemistry'] = 'v2'
    filtered_df_v2_cortex_neuron_maturity['layer'] = 'cortex'
    return


@app.cell
def _(
    filtered_df_v2_all_layer_neuron_maturity,
    filtered_df_v3_all_layer_neuron_maturity,
    pd,
):
    # 1. Concatenate the DataFrames along rows
    neuron_maturity = pd.concat([filtered_df_v3_all_layer_neuron_maturity, filtered_df_v2_all_layer_neuron_maturity], ignore_index=True)

    # 2. Sort by specific columns
    sorted_neuron_maturity = neuron_maturity.sort_values(
        by=["gene_list_file_name"], ascending=[True]
    )
    return (sorted_neuron_maturity,)


@app.cell
def _(sorted_neuron_maturity):
    sorted_neuron_maturity[["gene_list_file_name", "layer", "chemistry", "column_condition_title",
        "column_condition_value", "NES", "FDR_qval_BH_reported", "lead_genes", "num_genes_in_gene_list", "num_of_lead_genes"]]
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
