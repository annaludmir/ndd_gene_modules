import marimo

__generated_with = "0.17.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    import scanpy as sc
    import pandas as pd
    return pd, sc


@app.cell
def _(sc):
    def add_chemistry_column(
        in_h5ad: str,
        out_h5ad: str,
        v3_sample_ids,
        v2_sample_ids,
        sample_col: str = "sample_id",
        unknown_label: str = "unknown",
        overwrite: bool = True,
    ):
        """
        Add adata.obs['Chemistry'] based on sample_id membership in v3/v2 lists and save.
        """

        adata = sc.read_h5ad(in_h5ad)

        if sample_col not in adata.obs.columns:
            raise KeyError(f"'{sample_col}' not found in adata.obs columns: {list(adata.obs.columns)}")

        # Make sure sample IDs are comparable (avoid int/str mismatch)
        sample_ids = adata.obs[sample_col].astype(str)
        v3_set = set(map(str, v3_sample_ids))
        v2_set = set(map(str, v2_sample_ids))

        overlap = v3_set & v2_set
        if overlap:
            raise ValueError(f"These sample_ids appear in BOTH v3 and v2 lists: {sorted(list(overlap))[:20]}")

        if (not overwrite) and ("Chemistry" in adata.obs.columns):
            raise ValueError("Chemistry column already exists and overwrite=False.")

        # Default
        adata.obs["Chemistry"] = unknown_label

        # Assign
        adata.obs.loc[sample_ids.isin(v3_set), "Chemistry"] = "v3"
        adata.obs.loc[sample_ids.isin(v2_set), "Chemistry"] = "v2"

        # Quick sanity report
        print("Chemistry value counts:")
        print(adata.obs["Chemistry"].value_counts(dropna=False))

        # Save
        adata.write_h5ad(out_h5ad)
        print(f"Saved: {out_h5ad}")

        return adata  # optional, in case you want it in-memory too
    return (add_chemistry_column,)


@app.cell
def _(add_chemistry_column):
    # ---- Example usage ----
    in_h5ad = "/miridan-data/annaludmir/ndd_gene_modules/data/human_dev-GRCh38-3.0.0_all_layers_from_github.h5ad"
    out_h5ad = "/miridan-data/annaludmir/ndd_gene_modules/data/human_dev-GRCh38-3.0.0_all_layers_from_github_with_chemistry.h5ad"

    v3_sample_ids = [b'10X181_3',b'10X177_4',b'10X177_5',b'10X177_6',b'10X178_1',b'10X178_5',b'10X178_2',b'10X178_3',b'10X178_4',b'10X179_4',b'10X298_1',b'10X298_2',b'10X298_6',b'10X298_7',b'10X298_8',b'10X298_3',b'10X298_4',b'10X298_5',b'10X287_4',b'10X287_5',b'10X287_6',b'10X287_7',b'10X287_8',b'10X288_1',b'10X288_2',b'10X288_3',b'10X288_4',b'10X288_5',b'10X288_6',b'10X288_7',b'10X288_8',b'10X167_1',b'10X167_2',b'10X167_7',b'10X167_8',b'10X168_4',b'10X168_5',b'10X168_6',b'10X169_7',b'10X169_8',b'10X170_3',b'10X170_4',b'10X167_5',b'10X167_6',b'10X168_1',b'10X168_2',b'10X168_3',b'10X169_1',b'10X169_2',b'10X170_1',b'10X170_2',b'10X169_5',b'10X169_6',b'10X169_3',b'10X169_4',b'10X200_1',b'10X200_2',b'10X199_5',b'10X199_6',b'10X199_7',b'10X199_8',b'10X196_3',b'10X196_4',b'10X197_3',b'10X197_4',b'10X198_3',b'10X198_4',b'10X196_5',b'10X196_6',b'10X197_5',b'10X197_6',b'10X198_5',b'10X198_6',b'10X196_1',b'10X196_2',b'10X197_1',b'10X197_2',b'10X198_1',b'10X198_2',b'10X199_1',b'10X199_2',b'10X199_3',b'10X199_4',b'10X200_3',b'10X302_5',b'10X302_6',b'10X302_7',b'10X302_8',b'10X302_1',b'10X302_2',b'10X302_3',b'10X302_4',b'10X163_3',b'10X163_4',b'10X164_3',b'10X164_4',b'10X164_5',b'10X163_5',b'10X163_6',b'10X163_7',b'10X163_8',b'10X164_6',b'10X164_7',b'10X163_1',b'10X163_2',b'10X164_1',b'10X164_2',b'10X165_1',b'10X165_2',b'10X165_3',b'10X165_4',b'10X208_1',b'10X208_2',b'10X209_1',b'10X208_7',b'10X208_8',b'10X209_3',b'10X208_4',b'10X208_3',b'10X209_2',b'10X208_5',b'10X208_6',b'10X210_1',b'10X210_2',b'10X211_1',b'10X210_7',b'10X210_8',b'10X211_4',b'10X211_5',b'10X210_3',b'10X210_4',b'10X211_2',b'10X211_3',b'10X207_1',b'10X207_2',b'10X207_3',b'10X210_5',b'10X210_6',b'10X187_1',b'10X187_2',b'10X188_1',b'10X187_3',b'10X187_4',b'10X188_2',b'10X187_5',b'10X187_6',b'10X188_3',b'10X185_7',b'10X185_8',b'10X186_4',b'10X187_7',b'10X187_8',b'10X188_4',b'10X185_4',b'10X185_3',b'10X186_2',b'10X185_1',b'10X185_2',b'10X186_1',b'10X185_5',b'10X185_6',b'10X186_3',b'10X212_3',b'10X212_4',b'10X212_1',b'10X212_2',b'10X212_5',b'10X212_6',b'10X213_3',b'10X213_4',b'10X213_2',b'10X213_1',b'10X231_1',b'10X231_2',b'10X231_5',b'10X231_6',b'10X257_1',b'10X257_2',b'10X257_3',b'10X257_4',b'10X257_5',b'10X257_6',b'10X252_2',b'10X252_3',b'10X252_4',b'10X254_7',b'10X254_8',b'10X255_7',b'10X254_5',b'10X254_6',b'10X255_6',b'10X254_2',b'10X255_1',b'10X255_2',b'10X255_3',b'10X251_7',b'10X251_8',b'10X252_1',b'10X252_5',b'10X252_6',b'10X252_7',b'10X254_3',b'10X254_4',b'10X255_4',b'10X255_5',b'10X258_3',b'10X258_4',b'10X258_1',b'10X258_2',b'10X262_5']
    v2_sample_ids = [b'10X147_1',b'10X147_2',b'10X148_1',b'10X148_2',b'10X110_4',b'10X110_5',b'10X110_6',b'10X109_4',b'10X109_5',b'10X109_1',b'10X109_2',b'10X109_3',b'10X109_1',b'10X109_7',b'10X109_7',b'10X118_1',b'10X118_2',b'10X118_3',b'10X118_4',b'10X119_1',b'10X119_2',b'10X119_3',b'10X119_4',b'10X119_5',b'10X119_6',b'10X119_7',b'10X119_8',b'10X118_5',b'10X118_6',b'10X118_7',b'10X118_8',b'10X154_4',b'10X152_3',b'10X152_4',b'10X154_1',b'10X154_2',b'10X154_3',b'10X152_1',b'10X152_2',b'10X154_5',b'10X154_6',b'10X154_7',b'10X152_5',b'10X152_6',b'10X152_7',b'10X152_8',b'10X122_1',b'10X122_2',b'10X122_3',b'10X123_7',b'10X123_8',b'10X122_7',b'10X122_8',b'10X123_1',b'10X123_2',b'10X123_3',b'10X123_4',b'10X123_5',b'10X123_6',b'10X122_4',b'10X122_5',b'10X122_6',b'10X124_3',b'10X124_4',b'10X125_3',b'10X125_4',b'10X124_5',b'10X124_6',b'10X125_5',b'10X125_6',b'10X126_7',b'10X126_8',b'10X126_1',b'10X126_2',b'10X126_3',b'10X126_4',b'10X124_1',b'10X124_2',b'10X125_1',b'10X125_2',b'10X124_7',b'10X124_8',b'10X111_1',b'10X111_2',b'10X112_1',b'10X112_2',b'10X89_1',b'10X89_2',b'10X89_5',b'10X89_6',b'10X89_3',b'10X89_4',b'10X89_7',b'10X89_8',b'10X96_5',b'10X96_6',b'10X96_1',b'10X96_2',b'10X97_1',b'10X97_2',b'10X99_5',b'10X99_6',b'10X98_1',b'10X98_2',b'10X98_3',b'10X98_4',b'10X99_1',b'10X99_2',b'10X96_7',b'10X96_8',b'10X97_3',b'10X97_4',b'10X96_3',b'10X96_4',b'10X99_3',b'10X99_4',b'10X98_5',b'10X98_6',b'10X156_7',b'10X156_8',b'10X155_1',b'10X155_3',b'10X156_3',b'10X156_4',b'10X156_5',b'10X156_6',b'10X156_1',b'10X156_2',b'10X157_2',b'10X157_1',b'10X105_1',b'10X105_2',b'10X105_3',b'10X103_1',b'10X103_2',b'10X104_6',b'10X104_7',b'10X104_8',b'10X104_3',b'10X104_4',b'10X104_5',b'10X103_3',b'10X104_1',b'10X104_2',b'10X132_1',b'10X132_2',b'10X132_3',b'10X132_4',b'10X132_5',b'10X132_6',b'10X132_7',b'10X132_8',b'10X143_1',b'10X143_2',b'10X114_3',b'10X114_4',b'10X115_4',b'10X115_5',b'10X115_6',b'10X114_1',b'10X115_1',b'10X115_2',b'10X115_3',b'10X114_2',b'10X114_5',b'10X114_6',b'10X114_7',b'10X114_8',b'10X115_7',b'10X115_8',b'10X116_5',b'10X116_6',b'10X116_7',b'10X116_8',b'10X116_1',b'10X116_2',b'10X116_3',b'10X116_4',b'10X92_1',b'10X92_2',b'10X92_3',b'10X92_4',b'10X101_5',b'10X101_6',b'10X101_7',b'10X101_8',b'10X102_5',b'10X102_6',b'10X102_7',b'10X102_8',b'10X101_1',b'10X101_2',b'10X101_3',b'10X101_4',b'10X102_1',b'10X102_2',b'10X102_3',b'10X102_4']

    add_chemistry_column(in_h5ad, out_h5ad, v3_sample_ids, v2_sample_ids)
    return


@app.cell
def _(mo):
    mo.md(r"""
    #Remove bytes
    """)
    return


@app.cell
def _(pd):

    def decode_obs_bytes(adata):
        """
        Convert any byte or mixed byte/string columns in adata.obs into proper UTF-8 strings.
        Works for categorical and object dtype columns.
        """
        for col in adata.obs.columns:
            s = adata.obs[col]

            # Skip numeric columns
            if pd.api.types.is_numeric_dtype(s):
                continue

            # Decode object / categorical columns
            if s.dtype == "object" or pd.api.types.is_categorical_dtype(s):
                adata.obs[col] = (
                    s.astype(str)
                     .str.replace("^b'", "", regex=True)
                     .str.replace("'$", "", regex=True)
                )

        return adata
    return (decode_obs_bytes,)


@app.cell
def _(decode_obs_bytes, sc):
    #upload data
    adata=sc.read_h5ad("/miridan-data/annaludmir/ndd_gene_modules/data/human_dev-GRCh38-3.0.0_all_layers_from_github_with_chemistry.h5ad")
    adata_decoded = decode_obs_bytes(adata)
    adata_decoded.write_h5ad("/miridan-data/annaludmir/ndd_gene_modules/data/human_dev-GRCh38-3.0.0_all_layers_from_github_with_chemistry_no_bytes.h5ad")
    return (adata_decoded,)


@app.cell
def _(adata_decoded):
    adata_decoded.obs
    return


@app.cell
def _(sc):
    # load file
    adata_full = sc.read_h5ad("/miridan-data/annaludmir/ndd_gene_modules/data/human_dev.h5ad")

    # inspect unique values first (recommended)
    print(adata_full.obs["Age"].unique())

    return (adata_full,)


@app.cell
def _(adata_full):
    # remove week 5
    adata_filtered = adata_full[adata_full.obs["Age"] != "5"].copy()

    # save new file
    adata_filtered.write("/miridan-data/annaludmir/ndd_gene_modules/data/human_dev_without_week_5.h5ad")

    print("Saved successfully")
    print(adata_filtered.shape)
    return


@app.cell
def _(adata_full):
    for i in sorted(adata_full.obs["Age"].unique()):
        print(i)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
