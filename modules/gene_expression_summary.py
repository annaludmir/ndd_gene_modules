import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc


def collapse_region(region: str):
    """Map detailed regions to meta-regions."""
    if pd.isna(region):
        return np.nan

    forebrain = {"Forebrain", "Telencephalon", "Diencephalon"}
    hindbrain = {"Hindbrain", "Cerebellum", "Pons", "Medulla"}
    midbrain = {"Midbrain"}

    if region in forebrain:
        return "Forebrain"
    elif region in hindbrain:
        return "Hindbrain"
    elif region in midbrain:
        return "Midbrain"
    else:
        return np.nan  # ignore Brain / Head / others


def summarize_gene_expression_top_region(
    adata,
    genes,
    region_col="Region",
    cell_filter_col="CellClass",
    cell_filter_val="Radial glia",
    sym_col="Gene",
    expr_threshold=0,
):
    """
    Returns:
      1) long_summary: one row per gene x meta-region
      2) final_summary: one row per gene with top region and region-specific stats
    """

    print("Filtering cells...")
    adata_f = adata[adata.obs[cell_filter_col] == cell_filter_val].copy()

    print("Normalizing...")
    sc.pp.normalize_total(adata_f)
    sc.pp.log1p(adata_f)

    # Map regions to meta-regions
    print("Collapsing regions...")
    adata_f.obs["meta_region"] = adata_f.obs[region_col].map(collapse_region)

    # Keep only the 3 meta-regions
    adata_f = adata_f[adata_f.obs["meta_region"].notna()].copy()

    # gene symbol -> var name
    sym2var = (
        pd.Series(adata_f.var_names.values, index=adata_f.var[sym_col].astype(str))
        .dropna()
        .to_dict()
    )

    found_genes = [g for g in genes if g in sym2var]
    missing = [g for g in genes if g not in sym2var]

    print(f"Found {len(found_genes)} genes")
    print(f"Missing {len(missing)} genes")
    if missing:
        print("Missing symbols:", missing)

    varnames = [sym2var[g] for g in found_genes]

    # Extract expression matrix only for selected genes
    X = adata_f[:, varnames].X
    if sp.issparse(X):
        X = X.toarray()

    meta_regions = adata_f.obs["meta_region"].values
    long_results = []

    print("Calculating summaries...")

    for i, gene in enumerate(found_genes):
        vals = X[:, i]

        tmp = pd.DataFrame({
            "meta_region": meta_regions,
            "expr": vals
        })

        # all cells per region
        grouped_all = tmp.groupby("meta_region")["expr"]

        n_total = grouped_all.size().reset_index(name="n_total_cells")
        frac_expr = grouped_all.apply(lambda x: (x > expr_threshold).mean()).reset_index(name="fraction_expressing")

        # expressing cells only
        tmp_expr = tmp[tmp["expr"] > expr_threshold].copy()
        grouped_expr = tmp_expr.groupby("meta_region")["expr"]

        summary_expr = grouped_expr.agg(
            mean_expr="mean",
            median_expr="median",
            std_expr="std",
            min_expr="min",
            max_expr="max",
            n_expressing_cells="count"
        ).reset_index()

        # merge
        summary = n_total.merge(frac_expr, on="meta_region", how="left")
        summary = summary.merge(summary_expr, on="meta_region", how="left")

        # region score: intensity * prevalence
        summary["region_score"] = summary["mean_expr"] * summary["fraction_expressing"]

        summary["gene"] = gene
        long_results.append(summary)

    long_summary = pd.concat(long_results, ignore_index=True)

    # ---------- Build one-row-per-gene final summary ----------
    final_rows = []

    wanted_regions = ["Forebrain", "Midbrain", "Hindbrain"]

    for gene, sub in long_summary.groupby("gene"):
        sub = sub.set_index("meta_region").reindex(wanted_regions)

        # scores for ranking
        scores = sub["region_score"].fillna(0)

        top_region = scores.idxmax()
        top_score = scores.max()

        sorted_scores = scores.sort_values(ascending=False)
        second_region = sorted_scores.index[1] if len(sorted_scores) > 1 else np.nan
        second_score = sorted_scores.iloc[1] if len(sorted_scores) > 1 else np.nan

        enrichment_score = (
            top_score / second_score
            if pd.notna(second_score) and second_score > 0
            else np.nan
        )

        row = {
            "gene": gene,
            "top_region": top_region,
            "top_region_score": top_score,
            "second_region": second_region,
            "second_region_score": second_score,
            "enrichment_score_top_vs_second": enrichment_score,
        }

        for region in wanted_regions:
            row[f"{region}_mean_expr"] = sub.loc[region, "mean_expr"]
            row[f"{region}_median_expr"] = sub.loc[region, "median_expr"]
            row[f"{region}_std_expr"] = sub.loc[region, "std_expr"]
            row[f"{region}_fraction_expressing"] = sub.loc[region, "fraction_expressing"]
            row[f"{region}_n_total_cells"] = sub.loc[region, "n_total_cells"]
            row[f"{region}_n_expressing_cells"] = sub.loc[region, "n_expressing_cells"]
            row[f"{region}_region_score"] = sub.loc[region, "region_score"]

        final_rows.append(row)

    final_summary = pd.DataFrame(final_rows)

    # Sort by how strongly region-specific the gene is
    final_summary = final_summary.sort_values(
        ["top_region", "enrichment_score_top_vs_second"],
        ascending=[True, False]
    )

    return long_summary, final_summary


if __name__ == "__main__":
    genes = [
        'AARS1','ACBD6','ADARB1','AFG2B','AGMO','AKT3','ANKLE2','AP4B1','AP4E1','AP4M1','AP4S1','ARCN1',
        'ARF3','ARPC4','ASPM','ATP11A','ATP1A2','ATP6V0A1','ATP6V0C','ATP9A','ATR','ATRIP','ATRX','BLM',
        'BPTF','BRCA2','BRD4','BRIP1','BUB1','BUB1B','CAMK2B','CAMSAP1','CASK','CCDC88A','CCND2','CDC6',
        'CDK5RAP2','CDK6','CDT1','CENPE','CENPF','CEP135','CEP152','CEP55','CEP57','CEP63','CHAMP1','CHKA',
        'CIT','CKAP2L','COASY','COG3','COPB1','COPB2','CPAP','CPSF3','CREBBP','CRIPT','CSNK2A1','CTCF',
        'CTNNB1','CTSF','CTU2','DDX11','DHCR7','DIAPH1','DNA2','DNMT3A','DOHH','DONSON','DPM1','DPP6',
        'DROSHA','DYNC1I2','DYRK1A','EFTUD2','EIF2S3','EIF5A','EOMES','ERCC4','ERCC5','ERCC6','ERCC8',
        'EXOC7','FANCA','FANCB','FANCC','FANCD2','FANCE','FANCF','FANCG','FANCI','FANCL','FANCM','FBRSL1',
        'FILIP1','FOXG1','FRA10AC1','GINS2','GINS3','GMNN','GPT2','GRM7','GTF2E2','HDAC8','HHAT',
        'HIKESHI','HIST1H4C','HMGB1','HPDL','IARS1','IER3IP1','IGF1','IGF1R','INTS11','KIF11','KIF14',
        'KIF1BP','KMT2B','KNL1','LAGE3','LARP7','LHX2','LIG4','LMNB1','LMNB2','MCM7','MCPH1','MECP2',
        'MED11','MED17','METTL5','MFSD2A','MINPP1','MORC2','MPLKIP','MRE11','MSMO1','MYCN','NAA20','NAPB',
        'NARS1','NBN','NCAPD2','NCAPD3','NCAPH','NDE1','NHEJ1','NIN','NIPBL','NSD2','NSMCE2','NSRP1','NUF2',
        'NUP107','NUP188','NUP214','ORC1','ORC4','ORC6','OSGEP','PALB2','PCDH12','PCDHGC4','PCLO','PCNT',
        'PDCD6IP','PDHA1','PHC1','PLAA','PLK4','PNKP','POC1A','POGZ','PPFIBP1','PPIL1','PPP1R15B','PPP1R35',
        'PQBP1','PRIM1','PRUNE1','PSMC3','PTPN23','PUF60','PUS7','QARS','RAD21','RAD50','RAD51','RAD51C',
        'RBBP8','RING1','RMI1','RNU4-2','RNU4ATAC','RPL10','RRP7A','RTTN','RUSC2','SARS1','SASS6','SLC1A4',
        'SLC25A19','SLC38A3','SLC9A6','SLF2','SLX4','SMARCA5','SMC1A','SMC3','SMC5','SMG8','SMO','STAMBP',
        'STIL','SVBP','TAF13','TCF4','TMX2','TNPO2','TOP3A','TP53RK','TPR','TPRKB','TRA2B','TRAIP',
        'TRAPPC10','TRAPPC12','TRAPPC6B','TRAPPC9','TRIO','TRIP13','TRMT1','TRMT10A','TRRAPC14','TSEN15',
        'TSEN54','TTC5','TTI1','TUBG1','TUBGCP2','TUBGCP3','TUBGCP4','TUBGCP6','UBA5','UBE3A','UFC1','UFM1',
        'UGP2','UNC80','VPS50','VRK1','WDFY3','WDR11','WDR37','WDR4','WDR62','WDR73','WLS','XRCC4','YIF1B',
        'YIPF5','ZEB2','ZNF335','ZNF526','ZNF668','ZNHIT3','ZPR1'
    ]

    print("Uploading data...")
    adata = sc.read_h5ad("/miridan-data/annaludmir/ndd_gene_modules/data/human_dev.h5ad")

    long_summary, final_summary = summarize_gene_expression_top_region(
        adata,
        genes,
        region_col="Region",
        cell_filter_col="CellClass",
        cell_filter_val="Radial glia",
        sym_col="Gene",
        expr_threshold=0
    )

    # optional: save the detailed per-region table too
    long_summary.to_csv(
        "/miridan-data/annaludmir/ndd_gene_modules/results/additional_analyses/"
        "microcephaly_genes_region_long_summary_radial_glia_above_0.csv",
        index=False
    )

    # save the final one-row-per-gene table
    final_summary.to_csv(
        "/miridan-data/annaludmir/ndd_gene_modules/results/additional_analyses/"
        "microcephaly_genes_top_region_summary_radial_glia.csv",
        index=False
    )