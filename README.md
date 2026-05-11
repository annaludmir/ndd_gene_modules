# ndd_gene_modules

Analysis code for gene expression specificity (GES) scoring and downstream enrichment on developmental human brain single-cell data.

The main workflow in this repo is:

1. Compute GES tables for one or more cell groups from an `.h5ad` dataset
2. Run GSEA enrichment on those ranked GES tables using a gene list from a config file
3. Save summary tables and figures in date-stamped result folders

An alternative Tau-filtered GSEA workflow is also supported:

1. Compute global Tau specificity scores for all genes (from the same `.h5ad`)
2. Run GSEA enrichment, restricting the gene universe to the top 10% most cell-type-specific genes by Tau before preranking by GES

Downstream modules can also extract GO term or OMIM disease associations from the leading genes produced by either GSEA workflow.

## Repository layout

### `modules/`

Core analysis scripts.

**GES pipeline**

- `specificity_score_calculations.py`: main GES pipeline. Reads a YAML config, loads the AnnData object, optionally filters by chemistry, creates derived group columns, computes GES scores, and writes one CSV per target condition. Also contains `calculate_tau()` for computing per-target Tau specificity.
- `enrichment_pipeline_for_gene_list.py`: main enrichment pipeline. Reads a YAML config, creates a GMT file from a gene list if needed, runs GSEA on the GES outputs, and creates summary plots.
- `search_enrichment_gsea.py`: lower-level GSEA execution with `gseapy.prerank`, plus per-term enrichment plots.
- `search_enrichment_deseq.py`: lower-level DESeq enrichment execution.
- `get_gmt.py`: converts a gene list CSV into GMT format.
- `deseq_calculations.py`: DESeq2 / Fisher-based enrichment path for pseudobulk data.

**Tau-filtered GSEA pipeline**

- `tau_pipeline.py`: computes global Tau specificity scores for every gene across all groups in each condition column. Accepts the same YAML config as the GES pipeline. Saves one `tau_scores_{column}.csv` per column. Run this before `search_enrichment_gsea_tau_filtered.py`.
- `search_enrichment_gsea_tau_filtered.py`: drop-in replacement for `run_gsea()` that first filters genes to those in the top 10% of Tau (configurable), then applies the GES threshold, and runs prerank GSEA. Results are written to `GSEA_tau_filtered/` so they coexist with standard GSEA outputs.

**Downstream term extraction**

- `go_terms_pipeline.py`: runs Enrichr ORA (via `gseapy`) on leading genes extracted from a `GSEA_final_summary.csv`. Queries GO Biological Process, GO Molecular Function, GO Cellular Component, KEGG, and WikiPathways by default. Saves significant-hit tables per condition.
- `omim_pipeline.py`: same structure as `go_terms_pipeline.py` but queries OMIM Disease and OMIM Expanded gene sets instead of GO/pathway libraries. Outputs go to `results/OMIM/`.

**Visualization**

- `create_figs_ges.py` and `create_figs_ges_for_presentation.py`: bar plots for GSEA summary results.
- `create_figs_deseq.py`: bar plots for DESeq enrichment results.
- `plot_umaps.py`, `dot_plots.py`, `tsne_plots.py`: visualization helpers for genes and leading-edge sets.

**Analysis helpers**

- `gene_expression_summary.py`: summarizes where selected genes are most enriched across broad brain regions, with optional Wilcoxon comparisons and per-gene boxplots.
- `dataset_analysis_helper.py`: basic dataset summaries such as CellID counts by age, with optional chemistry and region filtering.
- `enrichment_cal_lists_loop.py`: batch runner that iterates every `*.csv` in a gene-lists folder, runs the standard GSEA enrichment pipeline for each file using a shared base config, optionally prunes non-significant run folders, and writes a single `batch_summary_{date}.csv` collecting NES and FDR results for all gene lists and condition columns.
- `enrichment_cal_lists_loop_tau_comparison.py`: extended batch runner that, for each gene-list CSV in a folder, runs all 12 combinations of scope (cortex / cell_phase / all_layers) × chemistry (v2 / v3) × tau filtering (with / without). Produces one `{gene_list_stem}_batch_summary_tau_vs_v2_v3_{date}.csv` per gene list with columns `scope`, `chemistry`, and `tau_filtered` for direct cross-variant comparison.
- `early_late_go_heatmap.py`: Early/Mid/Late developmental heatmaps for leading genes from a selected GSEA condition.
- `pseudotime_leading_genes_heatmap.py`: pseudotime-style heatmaps for leading genes across ordered ages.
- `leading_gene_condition_correlations.py`: scans enrichment result trees, computes pairwise Jaccard overlaps for all leading-gene condition pairs, and writes a correlation heatmap.

### `config_files/`

YAML configs for both steps of the pipeline:

- `ges_score_*.yaml`: configs for the GES scoring step (and Tau pipeline, which accepts the same format)
- `enrichment_*.yaml`: configs for the enrichment step. Variants exist for:
  - scope: `cortex`, `cortex_cell_phase`, `all_layers`
  - chemistry: plain (v3) and `_v2` suffix
  - tau filtering: plain and `_tau_filtered` suffix (six tau-filtered configs exist, covering both chemistries and all three scopes)

### `notebooks/`

Interactive analysis and figure-generation notebooks/scripts.

### `running_scripts/`

Cluster-oriented SLURM launchers for common runs.

- `python_gsea_enrichment.sh`: submit one enrichment run from one fixed config.
- `python_gsea_enrichment_all_configs.sh`: run one gene list across the main enrichment configs.
- `python_gsea_enrichment_folder_all_configs.sh`: run every `*.csv` gene list in a chosen folder across the main enrichment configs.
- `python_gsea_tau_filtered.sh`: run tau-filtered GSEA for one gene list across all three scope configs (cortex, cell_phase, all_layers).
- `python_gsea_tau_comparison.sh`: run the full 12-variant batch comparison (v2/v3 × tau/no-tau × all scopes) for every gene list in a folder, producing one summary CSV per gene list.
- `python_go_terms.sh`: run GO/pathway Enrichr ORA on leading genes from a GSEA summary CSV.
- `python_omim.sh`: run OMIM disease Enrichr ORA on leading genes from a GSEA summary CSV.
- `python_gene_expression_summary.sh`, `python_dataset_analysis_helper.sh`, `python_early_late.sh`, `python_pseudotime.sh`, `python_correlation_matrix.sh`: cluster launchers for the matching analysis helper modules.
- `python_ges_score_all_layers_v2*.sh`: cluster launchers for the GES scoring step across dataset subsets.

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

The GES pipeline (and the Tau pipeline) use configs like `config_files/ges_score_cortex_config.yaml`.

Important fields:

- `name_of_run`: base name for the output folder
- `data_path`: `.h5ad` input
- `output_folder`: parent output directory
- `column_conditions`: groups to score
- `expression_threshold`: minimum fraction of target cells expressing a gene
- `chemistry`: optional chemistry filter, for example `v3`
- `normalize_data`: whether to run normalization and log transform
- `tau_agg`: aggregation method used by the Tau pipeline (`mean` or `median`, default `mean`)

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

## Full GES → GSEA workflow

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

Make sure `column_conditions_for_gsea` matches the columns and targets scored in step 1. The enrichment step looks for files in:

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

### 3b. Run all main enrichment configs for one gene list on the cluster

If you want the same gene list to be evaluated with the standard `Cortex`, `Cell Phase`, and `All Layers` enrichment configs, use:

```bash
sbatch running_scripts/python_gsea_enrichment_all_configs.sh data/genes/my_gene_list.csv
```

The script derives the base run name from the CSV filename by default. You can also pass a custom base name as the second argument.

### 3c. Run all main enrichment configs for a whole folder of gene lists on the cluster

If you want to process a folder of gene lists in one batch, use:

```bash
sbatch running_scripts/python_gsea_enrichment_folder_all_configs.sh data/genes/my_gene_lists_folder
```

This script:

- finds every `*.csv` file in the folder
- runs `enrichment_cortex_config.yaml`, `enrichment_cortex_cell_phase_config.yaml`, and `enrichment_all_layers_config.yaml` for each file
- uses each CSV stem as the default run-name base

You can optionally add a shared prefix for all generated run names:

```bash
sbatch running_scripts/python_gsea_enrichment_folder_all_configs.sh data/genes/my_gene_lists_folder "AH MPRA"
```

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

## Tau-filtered GSEA workflow

This alternative workflow restricts the GSEA gene universe to the top 10% most cell-type-specific genes by Tau score before preranking by GES. This combines two complementary specificity signals:

- **Tau** filters out broadly expressed genes (keeps tissue-specific ones)
- **GES** then ranks the surviving genes for enrichment

### 1. Run the GES scoring step (same as standard workflow)

```bash
python modules/specificity_score_calculations.py config_files/ges_score_cortex_config.yaml
```

### 2. Compute global Tau scores

Use the same GES config file:

```bash
python modules/tau_pipeline.py config_files/ges_score_cortex_config.yaml
```

This creates a dated folder like:

```text
results/ges_score_results/ges_score_for_cortex_YYYYMMDD_tau/
```

Inside `data/`:

- `tau_scores_{column}.csv` — one file per condition column with columns: `gene`, `tau`, `max_group`, `mean_expr_{group}` for each group

`tau` ranges from 0 (uniform across all groups) to 1 (expressed in exactly one group).

### 3. Run the Tau-filtered GSEA enrichment

The simplest way to run all three scopes for one gene list on the cluster:

```bash
sbatch running_scripts/python_gsea_tau_filtered.sh data/genes/my_gene_list.csv
```

This runs `cortex`, `cell_phase`, and `all_layers` tau-filtered configs in sequence. Each scope creates its own dated output folder:

```text
results/enrichment_results/microcephaly_cortex_tau_filtered_tau90_threshold_1_YYYYMMDD/
results/enrichment_results/microcephaly_cell_phase_tau_filtered_tau90_threshold_1_YYYYMMDD/
results/enrichment_results/Microcephaly All Data (Without Week 5) tau filtered_tau90_threshold_1_YYYYMMDD/
```

Inside each folder the tau-filtered results land in:

```text
data/enrichment_results/GSEA_tau_filtered/   ← summary CSV and per-condition result CSVs
data/enrichment_figures/GSEA_tau_filtered/   ← per-column enrichment bar plots
```

You can also call the pipeline directly from Python:

Call `run_gsea_tau_filtered` from `search_enrichment_gsea_tau_filtered.py` in place of the standard `run_gsea`. It can be called directly from Python or integrated into a custom enrichment config runner:

```python
from search_enrichment_gsea_tau_filtered import run_gsea_tau_filtered

run_gsea_tau_filtered(
    ges_score_path="results/ges_score_results/ges_score_for_cortex_YYYYMMDD",
    tau_scores_dir="results/ges_score_results/ges_score_for_cortex_YYYYMMDD_tau/data",
    gmt_file="data/genes/microcephaly_genes.gmt",
    column_conditions={"CellClass": ["Radial glia", "Neuron"], "Region": ["Forebrain"]},
    ges_score_threshold=1,
    out_folder=enr_results_dir,
    figs_folder=fig_dir,
    tau_percentile=90.0,   # keep genes with tau >= 90th percentile (top 10%)
)
```

Results are written to `GSEA_tau_filtered/` under the enrichment output folder, so they coexist with standard `GSEA/` outputs. The summary CSV includes a `tau_percentile_cutoff` column for traceability.

## Cross-variant batch comparison (v2 / v3 × tau / no-tau)

To compare enrichment results across chemistries and tau-filtering for a folder of gene lists, use the `enrichment_cal_lists_loop_tau_comparison.py` pipeline.

### What it runs

For each gene-list CSV in the input folder, it executes 12 enrichment variants:

| Scope       | Chemistry | Tau filtered |
|-------------|-----------|--------------|
| cortex      | v3        | no           |
| cortex      | v2        | no           |
| cortex      | v3        | yes          |
| cortex      | v2        | yes          |
| cell_phase  | v3        | no           |
| cell_phase  | v2        | no           |
| cell_phase  | v3        | yes          |
| cell_phase  | v2        | yes          |
| all_layers  | v3        | no           |
| all_layers  | v2        | no           |
| all_layers  | v3        | yes          |
| all_layers  | v2        | yes          |

Each variant creates its own dated result folder under `results/enrichment_results/`.

### Batch summary output

One CSV is written per gene list:

```text
results/enrichment_results/{gene_list_stem}_batch_summary_tau_vs_v2_v3_{YYYYMMDD}.csv
```

The CSV has the same columns as the standard `batch_summary_{date}.csv` produced by `enrichment_cal_lists_loop.py`, plus three extra columns that identify the variant:

| Column        | Values                               |
|---------------|--------------------------------------|
| `scope`       | `cortex`, `cell_phase`, `all_layers` |
| `chemistry`   | `v2`, `v3`                           |
| `tau_filtered`| `True`, `False`                      |

### Running on the cluster

```bash
sbatch running_scripts/python_gsea_tau_comparison.sh data/genes/my_gene_lists_folder
```

### Note on v2 tau scores

The v2 tau-filtered configs currently point to the same `tau_scores_dir` as their v3 counterparts. Update the `tau_scores_dir` field in the three `*_tau_filtered_v2.yaml` config files once v2-specific tau scores have been computed with `tau_pipeline.py` on the v2 dataset.

## GO term and OMIM extraction

Both modules read a `GSEA_final_summary.csv` produced by either GSEA workflow, extract leading genes per condition, and run Enrichr ORA via `gseapy`.

### GO term enrichment

```bash
python modules/go_terms_pipeline.py \
  --summary-csv results/enrichment_results/my_run/data/enrichment_results/GSEA/GSEA_final_summary.csv
```

Default gene set libraries queried: `GO_Biological_Process_2021`, `GO_Molecular_Function_2021`, `GO_Cellular_Component_2021`, `KEGG_2021_Human`, `WikiPathway_2021_Human`.

Output folder is derived automatically as `results/GO_terms/{run_name}_GO_enrichment/`. Override with `--outdir`.

Optional arguments:

```
--column-conditions CellCyclePhase Region   # restrict to specific columns
--cutoff 0.05                               # adjusted p-value cutoff
--min-genes 5                               # skip conditions with fewer genes
```

### OMIM disease enrichment

```bash
python modules/omim_pipeline.py \
  --summary-csv results/enrichment_results/my_run/data/enrichment_results/GSEA/GSEA_final_summary.csv
```

Default gene set libraries queried: `OMIM_Disease`, `OMIM_Expanded`.

Output folder is derived automatically as `results/OMIM/{run_name}_OMIM_enrichment/`. Override with `--outdir`.

Accepts the same optional arguments as `go_terms_pipeline.py`, plus:

```
--gene-sets OMIM_Disease OMIM_Expanded   # override the queried libraries
```

Both modules write per-condition subfolders containing:

- `enrichr_results.csv`: full Enrichr results table
- `leading_genes_used.csv`: the gene list submitted
- `significant_hits_fdr_{cutoff}.csv`: terms passing the adjusted p-value cutoff

## Minimal example (standard workflow)

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
- a region-pair significant-gene count PNG with 3 divergent horizontal bar plots
- an optional folder of significant per-gene expression boxplots

The final summary also includes a one-sided paired Wilcoxon signed-rank test comparing `top_region_score` versus `second_region_score` across genes, along with an FDR/BH-corrected p-value. The per-gene Wilcoxon CSV compares the user-selected region against each other region using log-normalized expression and only cells with expression greater than `0` for that gene. That file includes both `region_of_interest` and `region_of_comparison` columns. If you pass `--wilcoxon-compare-all-region-pairs`, the per-gene CSV instead contains every ordered region-vs-region comparison, so with 3 regions you get 6 rows per gene. The PNG output always summarizes all 3 unordered region pairs, with one panel per pair and counts split by `combined_leading_gene_group`; each side counts genes with `wilcoxon_fdr_bh_pvalue < 0.05` for that region as the region of interest. If you also pass `--export-significant-gene-boxplots`, the module writes a `significant_gene_boxplots/` folder containing one 3-region boxplot per gene for each requested region of interest that is significant against both other regions. Those boxplots include only cells with expression greater than `0` for the plotted gene. If you provide an extra summary file and condition list, the output adds boolean columns such as `S_leading_gene` and `G2M_leading_gene`, and also writes one group-level Wilcoxon summary CSV with rows for `all`, `G1 & S`, and `S & G2M` when those flags are available. That summary reports the top region for each gene set and Wilcoxon results for the user-selected region versus each other available region, plus FDR/BH-corrected p-values across those group comparisons. By default the module filters the AnnData input to `Chemistry == v3`; use `--chemistry` to change that or `--chemistry None` to disable chemistry filtering.

Example call:

```bash
python modules/gene_expression_summary.py \
  --h5ad-path data/human_dev.h5ad \
  --gsea-summary-file path/to/GSEA_final_summary.csv \
  --condition Forebrain \
  --leading-gene-summary-file path/to/GSEA_final_summary.csv \
  --leading-gene-conditions S G2M G1 \
  --wilcoxon-region Forebrain \
  --wilcoxon-compare-all-region-pairs \
  --export-significant-gene-boxplots \
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
