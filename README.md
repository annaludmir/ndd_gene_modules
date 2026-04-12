# ndd_gene_modules

Analysis code for gene expression specificity (GES) scoring and downstream enrichment on developmental human brain single-cell data.

The main workflow in this repo is:

1. compute GES tables for one or more cell groups from an `.h5ad` dataset
2. run GSEA enrichment on those ranked GES tables using a gene list from a config file
3. save summary tables and figures in date-stamped result folders

## Repository layout

### `modules/`

Core analysis scripts.

- `specificity_score_calculations.py`: main GES pipeline. Reads a YAML config, loads the AnnData object, optionally filters by chemistry, creates derived group columns, computes GES scores, and writes one CSV per target condition.
- `enrichment_pipeline_for_gene_list.py`: main enrichment pipeline. Reads a YAML config, creates a GMT file from a gene list if needed, runs GSEA on the GES outputs, and creates summary plots.
- `search_enrichment_gsea.py`: lower-level GSEA execution with `gseapy.prerank`, plus per-term enrichment plots.
- `get_gmt.py`: converts a gene list CSV into GMT format.
- `deseq_calculations.py`: DESeq2 / Fisher-based enrichment path for pseudobulk data.
- `create_figs_ges.py` and `create_figs_ges_for_presentation.py`: bar plots for GSEA summary results.
- `create_figs_deseq.py`: bar plots for DESeq enrichment results.
- `gene_expression_summary.py`: helper analysis for summarizing where selected genes are most enriched across broad brain regions.
- `dataset_analysis_helper.py`: helper analysis utilities for basic dataset summaries such as sample counts by age.
- `plot_umaps.py`, `dot_plots.py`, `tsne_plots.py`: visualization helpers for genes and leading-edge sets.
- `enrichment_cal_lists_loop.py`: helper utilities for running enrichment across multiple gene lists.

### `config_files/`

YAML configs for both steps of the pipeline:

- `ges_score_*.yaml`: configs for the GES scoring step
- `enrichment_*.yaml`: configs for the enrichment step

### `notebooks/`

Interactive analysis and figure-generation notebooks/scripts.

### `results/`

Local output folder for generated tables and plots.

### `project_description.md`

Long-form project background, biological motivation, and notes on older analysis flows.

## Expected inputs

The config files assume the project root is set with:

```yaml
ndd_gene_modules_folder_root: "/path/to/ndd_gene_modules"
```

Most paths inside the YAML files are then resolved relative to that root.

Typical inputs are:

- an `.h5ad` single-cell dataset
- a gene list CSV with a `gene` column
- optionally, pseudobulk files for the DESeq branch

## Main configs

### GES config

The GES pipeline uses configs like `config_files/ges_score_cortex_config.yaml`.

Important fields:

- `name_of_run`: base name for the output folder
- `data_path`: `.h5ad` input
- `output_folder`: parent output directory
- `column_conditions`: groups to score
- `expression_threshold`: minimum fraction of target cells expressing a gene
- `chemistry`: optional chemistry filter, for example `v3`
- `normalize_data`: whether to run normalization and log transform

`column_conditions` can be either:

- a direct list of values already present in `adata.obs`
- a derived column created from a boolean expression or a set of rules

### Enrichment config

The enrichment pipeline uses configs like `config_files/enrichment_cortex_config.yaml`.

Important fields:

- `run_name`: base name for the enrichment run
- `output_folder`: parent output directory
- `ges_results_folder`: path to the GES run generated in step 1
- `gene_list_path`: gene list CSV
- `gmt_folder`: where the GMT file should be created or reused
- `analysis_mode`: `gsea`, `deseq`, or `both`
- `gsea.min_ges_score_threshold`: filter applied before GSEA
- `column_conditions_for_gsea`: which GES result files to use

## Full GES -> GSEA workflow

### 1. Run the GES scoring step

Pick or edit one of the `ges_score_*.yaml` files, then run:

```bash
python modules/specificity_score_calculations.py config_files/ges_score_cortex_config.yaml
```

This creates a dated run folder like:

```text
results/ges_score_results/ges_score_for_cortex_YYYYMMDD/
```

Inside that folder:

- `metadata/`: copied config, data hash, and an `edited_adata.h5ad` with derived columns
- `data/`: one file per target, for example `ges_spec_CellClass_Radial glia.csv`

### 2. Point the enrichment config to the new GES output

Open the matching enrichment config and set:

```yaml
ges_results_folder: "results/ges_score_results/ges_score_for_cortex_YYYYMMDD"
gene_list_path: "data/genes/microcephaly_genes.csv"
gmt_folder: "data/genes/"
analysis_mode: "gsea"
```

Make sure `column_conditions_for_gsea` matches the columns and targets you scored in step 1. The enrichment step looks for files in:

```text
<ges_results_folder>/data/ges_spec_<column>_<condition>.csv
```

### 3. Run the enrichment step

```bash
python modules/enrichment_pipeline_for_gene_list.py config_files/enrichment_cortex_config.yaml
```

What this does:

- loads the enrichment YAML
- creates `<gene_list_name>.gmt` in `gmt_folder` if it does not already exist
- runs GSEA on each requested column/condition
- saves a combined summary CSV
- creates per-column summary bar plots

### 4. Inspect the outputs

The enrichment run creates a dated folder like:

```text
results/enrichment_results/Microcephaly Cortex_threshold_1_YYYYMMDD/
```

Important subfolders:

- `metadata/`: copy of the enrichment config
- `data/enrichment_results/GSEA/`: per-condition GSEA tables and the combined summary
- `data/enrichment_figures/GSEA/`: summary bar plots

The main combined output is usually:

```text
data/enrichment_results/GSEA/GSEA_final_summary.csv
```

## Minimal example

If you want to run the cortex microcephaly workflow end to end:

```bash
python modules/specificity_score_calculations.py config_files/ges_score_cortex_config.yaml
python modules/enrichment_pipeline_for_gene_list.py config_files/enrichment_cortex_config.yaml
```

Before the second command, update `ges_results_folder` in `config_files/enrichment_cortex_config.yaml` so it points to the dated GES output created by the first command.

## Notes

- The gene list file used for GMT creation should contain a `gene` column.
- GES output folders are date-stamped, so the enrichment config usually needs to be updated after each new GES run.
- Several configs in `config_files/` are specialized for cortex, all-layers, proliferating cells, and cell-cycle analyses.
- `project_description.md` contains broader biological and historical context, including older scripts in `old_stuff/`.

## Gene expression summary example

The module `modules/gene_expression_summary.py` reads the leading genes for one selected GSEA summary condition, summarizes their expression across Forebrain/Midbrain/Hindbrain, and writes:

- a long per-gene-per-region CSV
- a final one-row-per-gene CSV
- a per-gene Wilcoxon comparison CSV for the user-selected region versus the other regions

The final summary also includes a one-sided paired Wilcoxon signed-rank test comparing `top_region_score` versus `second_region_score` across genes, along with an FDR/BH-corrected p-value. The per-gene Wilcoxon CSV compares the user-selected region against each other region using log-normalized expression and only cells with expression greater than `0` for that gene. If you provide an extra summary file and condition list, the output adds boolean columns such as `S_leading_gene` and `G2M_leading_gene`, and also writes one group-level Wilcoxon summary CSV with rows for `all`, `G1 & S`, and `S & G2M` when those flags are available. That summary reports the top region for each gene set and Wilcoxon results for the user-selected region versus each other available region, plus FDR/BH-corrected p-values across those group comparisons. By default the module filters the AnnData input to `Chemistry == v3`; use `--chemistry` to change that or `--chemistry None` to disable chemistry filtering.

Example call:

```bash
python modules/gene_expression_summary.py \
  --h5ad-path data/human_dev.h5ad \
  --gsea-summary-file path/to/GSEA_final_summary.csv \
  --condition Forebrain \
  --leading-gene-summary-file path/to/GSEA_final_summary.csv \
  --leading-gene-conditions S G2M G1 \
  --wilcoxon-region Forebrain \
  --chemistry v3 \
  --subfolder-name forebrain_gene_expression_summary
```

By default, outputs are written under:

```text
results/additional_analyses/gene_expression_summary/
```

## Early/Mid/Late GO heatmap example

The module `modules/early_late_go_heatmap.py` creates Early/Mid/Late heatmaps for the leading genes of a selected GSEA condition. It writes one GO-grouped heatmap with GO labels under the gene blocks and one expression-pattern-ordered heatmap, using these developmental stages:

- `Early`: weeks 5.5 to 7.0
- `Mid`: weeks 7.1 to 8.9
- `Late`: weeks 9.0 to 14.0

Example call:

```bash
python modules/early_late_go_heatmap.py \
  --gsea-summary-file path/to/GSEA_final_summary.csv \
  --condition Radial_glia \
  --go-term-file path/to/go_enrichment.csv \
  --h5ad-path data/human_dev.h5ad \
  --subfolder-name radial_glia_cell_cycle
```

The script reads the `Lead_genes` value from the row where the summary file `condition` column matches the provided `--condition`. By default it filters the AnnData input to `Chemistry == v3`; use `--chemistry` to choose a different chemistry or `--chemistry None` to disable that filter.

By default, outputs are written to:

```text
results/time_analysis/early_late/
```

If `--subfolder-name` is provided, outputs are written under:

```text
results/time_analysis/early_late/<your_subfolder_name>/
```

## Pseudotime leading-gene heatmap example

The module `modules/pseudotime_leading_genes_heatmap.py` creates pseudotime-style heatmaps for the leading genes of a selected GSEA condition. It reads the matching row from the GSEA summary file, uses the row `column` and `condition` values to filter cells, and plots gene-expression changes across ordered ages. If you provide `--go-term-file`, it writes both a GO-grouped heatmap and an expression-pattern-ordered heatmap. By default it filters the AnnData input to `Chemistry == v3`; use `--chemistry` to choose a different chemistry or `--chemistry None` to disable that filter.

Example call:

```bash
python modules/pseudotime_leading_genes_heatmap.py \
  --gsea-summary-file path/to/GSEA_final_summary.csv \
  --condition Forebrain \
  --h5ad-path data/human_dev_without_week_5.h5ad \
  --go-term-file path/to/enrichr_results.csv \
  --subfolder-name microcephaly_forebrain_pseudotime
```

By default, outputs are written to:

```text
results/time_analysis/pseudotime/
```

If `--subfolder-name` is provided, outputs are written under:

```text
results/time_analysis/pseudotime/<your_subfolder_name>/
```

## Leading-gene correlation example

The module `modules/leading_gene_condition_correlations.py` scans a folder tree for GSEA summary files, compares all leading-gene condition pairs, and writes a CSV of non-self condition pairs with correlation values above `0.5` and their overlap p-values, plus a heatmap of the strongest overlaps. Heatmap cells with Jaccard scores above `0.3` are annotated with overlap p-values.

Example call:

```bash
python modules/leading_gene_condition_correlations.py \
  results/enrichment_results/ \
  --top-n 25
```

By default, outputs are written under:

```text
results/correlations/<input_folder_name>/
```
