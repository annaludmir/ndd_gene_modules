# description of the project pipeline and work

this project explores human specific genes and human accelerated regions genes in human brain develoment cortex data
contect wise, the project icludes:

- **data**- h5ad files, gene lists. disease gene lists, speceficity tables, deseq results tables, gtf files etc..
- **codes**- python scripts that analyse and create the results data this include: ges score creation, and the enrichment pipeline
- **notebooks**- more interactive python scripts that includes data analysis, plots and more with description of the gene lists creation and results overview.
- **results**- results tables for enrichments and plots, some results are displayed or explored in the notebooks

## navigate
[data description and file locations](#data-description-and-file-locations)<br>
[gene lists used in the project](#gene-lists)<br>
[explanation about the project](#so-lets-begin-from-the-start)<br>
[speceficity](#speceficity)<br>
[important steps in ges calculation](#important-steps-in-ges-caculation)<br>
[code to create new ges tables](#code-to-create-new-ges-tables)<br>
[differential expression analysis](#differential-expression)<br>
[important steps in deg calculation](#important-steps-in-deg-calculation)<br>
[enrichmnet strategy](#enrichment_strategy)<br>
[the pipeline](#the-pipeline)<br>
[scripts that run the pipeline](#scripts-that-run-the-pipeline)<br>
[ges enrichment pipeline](#ges-enrichmnet-pipeline)<br>
[dge enrichment pipeline](#dge-enrichment-pipeline)<br>
[more notebooks to know about](#more-notebooks-to-know-about)<br>
[more scripts to know about](#more-scripts-to-know-about)<br>
[more data to know about](#more-data-to-know-about)

## data description and file locations

### data- 
the most important data files that have been the source for most of the work:
#### h5ad
<span style="color:red">**cloud location-h5ad_files directory**<span>

* the big h5ad file with all scRNA data- <span style="color:blue">/scratch200/reutj/data/human_dev_GRCh38-3.0.0_all_layers.h5ad</span><br>
updated- <span style="color:blue">/scratch200/reutj/data/updated_data_all.h5ad</span>

* cortex data- <span style="color:blue">/scratch200/reutj/data/cortex_adata.h5ad</span><br>
updated (with new criteria and slight modifications)- <span style="color:blue">/scratch200/reutj/data/updated_cortex_data_hg38.h5ad</span>

#### speceficity tables- 
<span style="color:red">**cloud location-spec_tables_new directory**</span><br>
tables of genes specific to each cell type that specified for the cortex or big data, those were caculated using the ges score for the normalized data counts. each tablke incluse the genes with their ges_score, mean exp in the target cell and in the other cell types compared based on the condition (cell type, cell cycle..).

those tables are found in **/scratch200/data/spec_tables_new**

#### deseq2 pseudobulks 
<span style="color:red">**cloud location-deseq2_pseudobulks
directory**</span><br>
differential analysis was done in deseq2 using pydeseq which stores the dds object for each big comparison group (cell cycle, cell class, etc..) in an h5ad file format. those are found in - **/scratch200/data/deseq/** under cortex or data_all folders

#### notebooks
<span style="color:red">**cloud location-notebooks and scripts directories**</span><br>
all the jupyter notebboks and python codes and bash codes that were used in this project are stored in the **notebooks** folder

#### out files
<span style="color:red">**cloud location-out_files directory**</span><br>
log files of some of the sbatch scripts that I run , mostly ges and fisher out files. and also some full out files with how a full pipeline run should look.

#### conda environment setting
most of the scripts and work was done in the jupyter-scanpy environment later change to jupyter-scanpy new. you cannsee the list of packages in the environment with their versions in /scratch200/reutj/jupyter_scanpy_conda_list.txt or simply just **jupyter_scanpy_conda_list.txt** in the cloud.

#### results
all the results from the ges and deseq enrichments and plots created by the main pipline. saved in the **enricment_results** directory

### gene lists
<span style="color:red">**cloud location-hsg_gene_lists directory**<span><br>
<span style="color:blue">**cluster location-/scratch200/reutj/data/hsg_gene_lists directory**<span><br>

we have gene lists from 3 types :
- **hsg**- those are human specific genes were sequence changes and duplications cause creation of unique human genes.
- **har**- genes that are assocuated with human accelareate regions identifies by proximity r functional assyas such as Hi-C to connect hAR to theier proposed gene targets.
- disease genes- these include lists of genes implicated in microcephaly, intellectal disability panel and autism.

**each of this lists also have gmt file format to support gsea analysis on that list**

* the lists on the folder that are not listed below are the basis to create the used lists here I am going to present just the list that I have chceked in the project

##### hsg lists
- <span style="color:blue">updated_hsg_list </span>- list of all genes with human specific chnages such as duplications, de-novo and sequence changes not including har. 
- <span style="color:blue">dup_genes_list</span>-list of just the duplicated genes 
- <span style="color:blue">seq_change_hsg_list</span> just the genes that went through sequence changes or de-novo- 

those lists are based of of smaller lists taken from different papers - de-novo, hsg, duplications and the smaller lists and creation is described in the relevant notebook:**creation_of_hsg_gene_list**

##### har lists
here we have different kinds of lists taken from diffenr paper

1. <span style="color:blue">updated_har_list</span>- first list we did analysis on with more than 1600 genes taken from [Won. et a](https://www.nature.com/articles/s41467-019-10248-3) , those genes were captured based on HI-C data from fetal brain in the cortical plate and germinal zones from mid-gestation embryos. note that some of the old measdurments was doen on the old list, since then i have added some modifiations to unidentifies gene names to increase the number of identifies gene sin the list. th old list is located in har_associated_genes.xlsx file in the hsg_gene_lists folder. code to update to new list in the notebook- **update_har_genes_list**.
2. <span style="color:blue">har700_df</span>- 700+ genes that were associted to har based on HI-C from different tissue types (mostly fibroblasts) and not from brain taken from [Walsh. et al](https://pmc.ncbi.nlm.nih.gov/articles/PMC5063026/) 
3. <span style="color:blue">har2000_df</span>-  2000+ genes associated by genomic proximity to HARs. taken also from [Walsh. et al](https://pmc.ncbi.nlm.nih.gov/articles/PMC5063026/). 
4. <span style="color:blue">har_genes_nsc_df</span>- genes associated to HAR in neural stem cells based on Chi-C new analysis. I redid the analysis to associate  the har to genes and 5kb upstream based on their results with some filtering. sourch [Noonan new cell article](https://www.cell.com/cell/fulltext/S0092-8674(25)00036-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867425000364%3Fshowall%3Dtrue). path-
creation of this list is described in the notebook - **new_har_gene_lists**
5. <span style="color:blue">har_genes_neurons_df</span>- genes associated to HAR in neuronal organoids based on the same new article . I redid the analysis to associate the har to genes and 5kb upstream based on their results with some filtering. 
creation of this list is also described in the notebook -<span style="color:blue">new_har_gene_lists <span>

* how i did the filtering to get to the last 2 gene lists?
the paper contains about 4000 genes associated to HAR in nsc or neurons cultures based on comprehensive Hi-c analysis the researchers chose areas capturing the HAR that include also the ares from the sides - secondary and teriatry. and look for genomic interactions with threshold of 4.6 (common threshold is 5 for the strength of interactions) to areas that unclude cofind genes and also 5kb upstream to them. to filter to more specific results I took only interactions that were in the prinmary region that contain the HAR with strength of intertaction >5 and used the GTF gencode annotation to assocaute all coding genes including 5kb upstream form the TSS as they  

##### disease lists

- <span style="color:blue">microcephaly_genes</span>- genes implicated in primary microcephaly taken fron panelapp genes. source [here](https://panelapp.genomicsengland.co.uk/panels/162/). 
- <span style="color:blue">id_list_updated</span>- genes implicated in intellectual disability and NDDs taken also from panelapp. [sourch here](https://panelapp.genomicsengland.co.uk/panels/285/). 
- <span style="color:blue">/scratch200/reutj/data/hsg_gene_lists/autism_genes.csv</span>- autism (divide also genes ranked as top autism genes (1) in a saparate **autism strong** list) genes takem form SFARI site. the list contains genws assoited to autism bu different levels of cofidence 1 being the most related and then 2 and 3 , the list also contains genes that cause syndromic disease when mutated that contain an autism symptom along with others (marked as "s" in the table) . source [here](https://gene.sfari.org/). 

* creation of those lists is described in the notebook -**create_disease_gene_lists**

## so lets begin from the start

basically the purpose is to look for enrichments of the human specific genes or HAR genes in the single cell data. to start we wanted to focuse on individual cell types each time and look for enrichments there hoping to find interesting patterns and genes. we started of by using two gene lists:

* hsg- the big hsg list that we strated which is comprised of de-novo, duplications and gene with sequence chnages taken from different sources and papers described in the `creation of the hsg gene lists` notebook the original list is in **data/hsg_gene_lists/updated_human_specific_genes3.csv** since i have gone over the list updated annotaions and so on and the updated list is in the same filder under the name - updated_hsg_list. and this list is very similar to the previous list.

* har-
1600 genes taken from [Won. et a](https://www.nature.com/articles/s41467-019-10248-3) described above.
general overview of those lists in the big data (mean exp, looking at different regions) is found in the notebook- `general_overview_on_data`

general overview of those in the cortex data is found in the `general_expression_hsg_in_cortex` notebook.

## what we wanted to do?

As i wrote here. our idea was to look for genes that identify specific cell type by their expression and then see if hsg or har genes are over represented in those lists. for that we used 2 approaches 

### speceficity
specific genes are genes that are uniquly or mostly expressed in one type of cell and do not express or express very little in other cell types. the requirment from a true specific gene is to be expressed in the target cell population, even expressed in a little fraction of the target cell population (here we set threshold of 5% expression fraction) and to have mean expression greater in the target cell then all the other cell types. 

* note that this is different from diffrentially expressed genes where other cell tyopes can express the gene and the level of expression is what changes, also diffrentially expressed genes usually expressed in large fraction of the target cell while here we get a lot of genes hich express only in a fraction of the target cell type

#### so how we calculate speceficity?

this is done by using ges score. shown here: ![ges_equation](../data/ges_equation.png)

- if image is not shown edit the path to the ges_equation.png file path. in the cloud- more_data/ges_equation.png

**Pt** is the fraction of cells in the data that belong to the target cell type. what the equation does is to divided the **mean exp** of the said gene , multiplied by (1-Pt) to **weighthed mean exp** of the other cell types meaning a sum of each cell type mean exp of the gene multiplied by its fraction in the data. that way all other cell types are considered in the calculation with each given the proportionate weight. the (1-Pt) correct the answer so that it is centered in 1:

* ges score>1 - genes that are specific to the cell target (the biggere the score the more the gene is specific)
* ges score <1 - genes that are more specific to other cell types and not specific to the target cell type.

#### important steps in ges caculation

1. data is uploaded - cortex/ entire brain. **in cortex data the adata is further filtered to V3 chemistry** since it contains more gene expression and cells and considered better than V2 and to not create a batch effect.
2. basic filtering to the data to cells with more then 200 expressing genes and genes that express in at least 3 cells
``` python
   sc.pp.filter_cells(adata, min_genes=200)
   sc.pp.filter_genes(adata, min_cells=3)
``` 
3. data is **normalized** using
``` python
   sc.pp.normalize_total(adata)
   sc.pp.log1p(adata)
```
4. for each target cell only genes that express in more than 5% of the target cell population will be considered in the ges table
5. ges score is caculated to each gene, with more detsils such as mean exp in each cell type..
6. *optional* - permutations can be done on random cell type lables to create ges scores on all the genes and the p value is caculated by the number of real ges score that are  smaller than the permutated. I tried with 500 permutations but most p values were completely 0 or 1 so it dosent provide much value (almost 8000 genes with p_val<0.05). even with padj.

**the resulted file**- table with each gene- gene name, the ges score , expression in the cell target and more details..

#### code to create new ges tables

* the creation of the ges scores was done using the **ges_score_corrected_no_permutations.py** script in the notebooks/scripts directory.
* bash scripts to run the ges_creation script-<br>
  for the cortex data - **ges_sbatch.sh**<br>
  for the big data- **ges_sbatch_big.sh**

##### usage:
```bash
python -u ges_score_corrected_no_permutations.py [column list] [conditions list] 
```
*column list*- list of columns names in the adata_obs that contain criteria of cell types to create ges scores on 
*conditions list*- list of conditions inside those columns that we want the ges score to be created

### differential expression

diffrentially expressed genes are genes that are expressed more in the target cell than the other types of cells. to caulcate deg in single cell there are many ways. I chose to use pseudobulks which takes groups of cells and sum them to behave like bulk RNA then use DESEQ2 to do the analysis. 

we have 4 pseudobulks,2 for the big data and 2 for cortex data:
- big data- cell class
- big data- region
- cortex- cell class
- cortex- cell cycle

*those pseudobulks are in /scratch200/data/deseq2/ in the cortex or data_all files or deseq_pseudobulks directory in the cloud

then from the pseudobulks pydeseq2 analysis can create deg table for each comparison between the condition in the pseudobulks (each condition against al the others, inter-comparison between one cell type and another is also posssible but i havent done it here)

* the creation of those deg pseudobulks is described in the notebook **deseq2_differential_expression**
* the creation of tables is with the enrichment analysis (which we will talk about later) in the **deseq_calculations.py** script


#### important steps in deg calculation

1.data is uploaded, cortex data was filtered to V3 chemistry and then wwnt through basic filtering using 
``` python
   sc.pp.filter_cells(adata, min_genes=200)
   sc.pp.filter_genes(adata, min_cells=3)
``` 
2. data is grouped to Donor and comparison criteria .for example:
``` python
   adata_v3.obs.groupby(['Donor','CellClass']).size() 
```
then **if the size of the group > 30** each of this subsections is sliced into individual h5ad. also if specific donor sample have most of its groups with very little number of cells or extreme outliers (like erythrocytes in big data - cell class). those are dropped

3. all those h5ad are concated to one h5ad. 
4. edit the obs of the new object to inclucde individual comparisons (each cell against others) and sub categories that invlove more than one cell type (for example proliferating) and this is the pseudobulks
5. the obs from the pseudobulks are the basis for the design to the deg analysis and the peseudobulk also behaves as the DDS object. to start the analysis each time different criteria is chosen within the pseudobulk and PYDESQ2 analysis is done.
6. results from the dds table arfe saved in the summery object. and file is saved to csv in the deseq2 directory.
7.  to filter significant genes -the threshold for significance results was set to **padj <0.05 , log2FC >1.5 , basemean >5**

* tables and examples is in the **deseq2_differential_expression** notebook
* overview of top genes in some groups in the **deseq2_differential_expression** and **top_sig_genmes_inspection** notebooks

### enrichment strategy 

the enrichment analysis idea is to take the hsg lists and to see if they are ennriched among the top deg or specific genes to get the idea about the functions of those hsg and which cells and stages in development they take part in. and also to get specific genes that are in the top of those enricment.
strategy is divןded to 2:

1. in speceficity list- with gsea to view the distribution of the genes in the list in the specificity score. high enrichmnet meaning a lot of genes appearing in the top places of the list.
2. in dge list- using fisher exact test to see if the proportion of hsg genes inside the significant dge list is higher than expected.

the enrichment analysis is done using the pipeline. description..

## the pipeline

the overall pipeline description is this:

      Gene list
         ↓
    ENRICHMENT TYPE
     ↙      ↘
    GES      DGE
     ↓        ↓
    gsea    fisher
     ↓        ↓
    figure   figure

<span style="color:blue">**main pipeline script- gene_list_pipeline**</span>

### run parameters:
##### usage:
```bash
python -u gene_list_piepline.py gene_list out_folder [ges_enrichment | dge_enrichment][data_all | cortex]  extra_args
```

1. **gene_list**- complete path to the gene list csv file. the gene list should be a dataframe with 2 neccesary columns- gene- gene symbols. and ens- matching gene ens annnotation. **important!** column names should be gene and ens with small letters
2. **out_folder**- out folder path where the results should be located. in my data the path is <span style="color:blue"> /scratch200/reutj/enrichment_results/ </span>
3. **choice** - whether to do ges enrichmnet analysis or sge. written like that - *ges_enrichment | dge_enrichment*
4. **data_type**- wheter to perform the analysis on the cortex data o the bog data. written like this- *data_all | cortex*
5. **extra_args** - additional arguments based on the *choice*. <br><br>
**in ges_enrichmnet:**
- folder with gmt_paths - should be the folder where the gene lists is located. in my data <span style="color:blue"> /scratch200/reutj/data/hsg_gene_lists/ </span><br>
- conditions to check- a list of conditions to check ges scores on, *notice* that the conditions should match the naming in the spec table folder directory spec_tables_new. so for example- use "radialglia" as the file name is ges_spec_V3_radialglia.csv and not "Radial glia" as written in the adta.obs.<br><br>
**in dge_enrichment:**
- pseudobulks data file folder - *note* that the pseudobulks should inside directory name that matches the data type - meaning cortex or data/all
- file with gene names of the genes in data that express in more than 100 cells. <br>
in the big data-
<span style="color:blue">/scratch200/reutj/data/all_data_ens_above100.txt</span> <br>
in the cortex data- <span style="color:blue">/scratch200/reutj/data/cortex_over_100_genes.txt</span>

### scripts that run the pipeline
for ges enrichment:
1. run_pipeline_ges.sh - for the cortex data
2. run_pipeline_ges_big.sh - for the big brain data
for dge enrichmnet:
1. run_pipeline_deseq.sh - for the cortex data
2. run_pipeline_deseq_big.sh - for the big brain data

those scripts contain all the parmeters needed and conditions and to run them on different gene lists we just need to change the gene_list name 

### ges enrichment pipeline
the pipeline contain 3 steps:
#### gmt file creation
creation of .gmt file of the gene symbols included in the provided gene list , the file is then saved in the directory of the gene list. if gmt file already exusts the program continue to the next level.<br>
<span style="color:blue">python file</span> - get_gmt.py <br>
<span style="color:blue">parameters</span>- the gene list path
#### gsea on the ges table calculation
here the program takes the hsg gene names from the gmt file and all the ges scores files with the condition specificed in the extra_args. Then, to make the enrichment more accurate **ges tables are filtered only to genes with ges score >1**  meaning genes that tend to or are more specific to the target cell. and uses *gseapy* to look for the positions of the given gene lists genes in the ges score ordered table (high ges scores on top) each calculation provides **NES score** which represent the strength off the enrichment based on the distribution of the genes at the top section of the ges list. the program also provides lists of hsg lead genes in the ges table with their percentage among all the hsg genes that were identified within each ges table. and also provides p-value for the calculation. gseaplots showing the distribution of the hsg within the ges tables with representation of the enrichment plots is drawn to each calculation an saved as well withh all the gseapy analysis files and results. and the results from all ges tables are organizes in a dataframe and adj p-value to correct for miltiple testing is performed on all of the other tests.

* it is important to note that beacuse of the relative small number of tests performed and the similary between some ges tables (for example- proliferating cells contain NPC , radial glia , and ipc which all have ges tables) the padj is not so good at reflecting the significane of the enrichment therefor NES score is a better measure. with those measures- <br>  **good enrichment- NES score >1.5 <br> strong enrichment- NES score >2**

* more information on this step is found in [important steps in ges calculation](#important-steps-in-ges-caculation) section.

<span style="color:blue">python file</span> - search_enrichmnet_new.py <br>
<span style="color:blue">parameters</span>- 
1. the gene list path- taken from the gene_list parameter in the pipeline parameters 
2. the gmt folder path, should be the one where the hsg gene list is located- extra_args[1] in the pipeline parameters.
3. condition to check- list of condition/cell_types to calculate enrichment on. the condition names should be exactly as written in the file names of the ges score tables- extra_args[2] in the pipeline parameters.
4. out_folder- the folder where results would be stored - same out_folder that is specified in the pipeline parameters.
5. data_type- whether the ges_scores tables that the program runs or are form the bigger data - "data_all" or cortex data- "cortex". this parameter is taken fron the data_type parameter of the pipeline. 

* notice the program dont take as a aparmeter the location of the ges files folder (spec_tables_new),and its written in the code itself as the path written in my user. so if you have another ges tables and folders it should be changed inside the script.

<span style="color:blue">outputs</span>- 
1. raw ges results- folder with the raw results and caculation steps, including the gme gene set used, preranking and intial results. saved into this specific directory:<br>

    [out_folder]/ges_results_above_1/[data_type]/[gene_list_name]/[condition]/gsea_raw/

2. gsea results table for each cell type. located in the same folder as gsea_raw folder as *gsea_results.csv*

3. gsea plot- plot for the gene list enrichmtnrt with nes score , gene list distrbution within the ges table list and gsea enrichment plot all in one plot. saved in the gsea_raw folder as:<br>

    custom_gseaplot_[gene_list_name].png

4. summerized table with the results for all measurments done for the specific gene list includes nes socre, top genes, p val amd padj (fdr-bh) for each condtion. save in this path:<br>

    [out_folder]/ges_results_above_1/[data_type]/[gene_list_name]/full_results.csv

#### gsea enrichment reaults plotting

the program takes the full_results.csv from the previous step to present barplot produced for each gene list representing the nes scores. negative enrichment scores (NES<0-  meaning gene list os enriched in the button of the ges list) are presented in red. dashed line at 1.5 gives the thresdhold for significant results and above 2 to very strong results. for each data_type 2 plots are created. <br>
**big data**- regions plot and cell class.<br>
**cortex**- cell class and cell cycle plots

* note that the plots are desinged specificaly to fit all the cell types and conditions in an order i set so if performeng ges caculations on other condtion orless condition the plot won't work , to change that you can go into the sciprt and change the line numbers and their order based on the line numbers from the input gsea results table to modify that to the cell_types and results youwant to display.

<span style="color:blue">python file</span> - create_figs_ges.py <br>
<span style="color:blue">parameters</span>- 
1. the gene list path- taken from the gene_list parameter in the pipeline parameters 
2. out_folder- taken from the out_folder parameter in the pipeline parameters
3. data_type- data_all/cortex, taken from the data_type parameter in the pipeline parameters <br>

<span style="color:blue">output</span><br>the figure as png file into this path:<br>

[out_folder]/figs/[gene_list_name]/[data_type]/[ges_results]/[comparison_type]_enrichment.png

*comparison type*- cell class/ cell cycle / regions.

### dge enrichmnet pipeline
this pipeline contain 2 steps:

#### fisher exact calculation on deseq results-

1. the script takes a file with list of genes expressed in the relevant data_type filtered to genes that express in more than 100 genes.
2. the genes in the hsg gene list intersected with this list to create **'genes_in_data'** list.
3. this script takes the pseudobulks relevant to the data type chosen (data_all/cortex) and performes **differential analysis** using *pydeseq2* on all the columns in the psueobulk.obs that contain comparisons of 1 cell type vs others (all column except the first one that contains inter comparisons radialglia vs IPC , radialglia vs neuron etc..).
4. the **dds results** are then saved into csv file. they contain gene names (in the big data gene ens) with deg information - meanbase, log2FC, pvalue , padj..
5. significant upregulated results for deseq results  list are set with those thresholds-
```python
de[(de.padj<0.05)&(de.log2FoldChange>1.5)&(de.baseMean>5)]
```
6. fisher test is done on the 'genes_in_data' with the each deseq2 results table. this is done using **fisher_exact** function from scipy.stats. background for the total amount of genes is the number of fenes in the above 100 file. this function provide for each calculation the names of the genes in the intersection between the lists, and the p_value.
7. results from all comparison are summed up to one sumersing table and p-adj (fdr-bh) to all the comparisons is caculated and saved to csv. 

* more information on this step is found in [important steps in deg calculation](#important-steps-in-deg-caculation) section.

<span style="color:blue">python file</span> - deseq_calculations.py <br>
<span style="color:blue">parameters</span>- 
1. the gene list path- taken from the gene_list parameter in the pipeline parameters
2. psb_data_folder- the folder with the pseudoblk h5ad files, note that the pseudobulk reside in one parent directory split to cortex and data_all direcories that contain the relvant psb to each data_type. here the parent directory should be specified. this is passed from the pipeline extra_args[1] argument.
3. gene_names- path to a txt file that contain the gene names of genes that express in more than 100 cells in the relavent data. the txt file format is of each gene name in a line.  
4. out_folder- taken from the out_folder parameter in the pipeline parameters
5. data_type- data_all/cortex, taken from the data_type parameter in the pipeline parameters<br>

<span style="color:blue">output</span>- 
1. differntially expressed gene table from the pydeseq2 analysis.sorted by padj ascending and saved in <br>

[out_folder]/deseq_fisher_results/[data_type]/[gene_list_name]/[pseusdobulk_name]/[comparison- condition_vs_other].csv

2. full results table with intersection gene names (note that in the data_all ens and not gene symbols are used, I havent fixed that yet but its easy to get gene symbols from ens when you have the list), num of intersecting genes, p_val and p_adj to each comparison. saved in the same directory as the dds results as - <br>

    full_results.csv

* **note** that I forgot to add .csv ending to the deseq results file names so i modified the script now but all the results from before are saved without the csv ending. in the cloud. in the cluster itself I have changed the ending. also the new script was uploaded under deseq_calculation_new.py with this small corection (not critical). 

#### plotting the results

very similar to the ges plotting.the program takes the full_results.csv from the previous step to present barplot produced for each gene list representing the -log(10) for the pvalues. dashed line at -log10(0.05) sets the thresdhold for significant results . for each data_type 2 plots are created based on the pseudobulk datas. <br>
**big data**- regions plot and cell class.<br>
**cortex**- cell class and cell cycle plots 

<span style="color:blue">python file</span> - create_figs_deseq.py <br>
<span style="color:blue">parameters</span>- 
1. the gene list path- taken from the gene_list parameter in the pipeline parameters 
2. out_folder- taken from the out_folder parameter in the pipeline parameters
3. data_type- data_all/cortex, taken from the data_type parameter in the pipeline parameters <br>

<span style="color:blue">output</span><br>the figure as png file into this path:<br>

[out_folder]/figs/[gene_list_name]/[data_type]/[deseq_results]/[comparison_type]_enrichment.png

* note that the plots are desinged specificaly to fit all the cell types and conditions in an order i set so if performeng ges caculations on other condtion orless condition the plot won't work , to change that you can go into the sciprt and change the line numbers and their order based on the line numbers from the input gsea results table to modify that to the cell_types and results youwant to display.

#### more notebooks to know about
<span style="color:blue">check_hsg_exp</span>- check hsg genes (old notebook so its the old hsg gene list and i thing also the old har list) expression in the cortex data- if they are expressed and in how many cells. also divided to celltypes.. <br>
<span style="color:blue">human specific genes search</span>- blast analysis on the 400 last bases of the 3utr of the hsg genes (the area that was sequenced) to check where it would be mapped to the genome using gtf data and check the level of similary and changes between very similar areas. this was done beacuse if the large number of duplicated hsg genes that can be mapped to more than one location and confuse the reads quantificantion.with that analysis i got the **hsg_with_paralogs.csv** file in the hsg_lists folder that shos for each hsg query its matches with percent identity, only genes with paralogs are shown in the table.<br>  
<span style="color:blue">ges_results_intersection</span>- venn plots and intersection of the lead genes provided for the different cell types in one gsea calculation. this was doen to see if the ges_genes are similar with similar gene_types and if the lead genes from the hsg match to them as well or shows some specificty to each contion tested<br>
<span style="color:blue">final_results_plots</span>- barplots of the old results from gsea enrichment and deseq enrichment analysis<br>
<span style="color:blue">bunch of staff</span>- as the name suggests include many things, for example the code to update and add criteria such as proliferating/differentiating to the adata files. 


#### more scripts to know about
* <span style="color:blue">speceficity_script</span>, ges_manwhitney - old ways to calculate speceficity, not used forward..
* <span style="color:blue">ges_scores_forreal</span> - the former ges table cauclation script which almost accurate just missing the (1-Pt) part in the equation. the results provided are basicaly the sameas when looking at the corrected results just a little shofted in the numbers and cannot be distinguished by ges score >1.
* <span style="color:blue">search_enrichment</span>- very similar to search_enrichmnt new .was used also to calculte gsea on ges tables.
* <span style="color:blue">search enrichment deseq</span>- script to calculate gsea analysis on full deseq tables (ordered by padj ascending) , another way to find locations of hsg genes inside the dge lists, the thing is that in deseq encrimnet we can prvide resonable thereshold for singnificance so tho smethod which looks at all the genes is probably less efficient.
* <span style="color:blue">bed_to_df</span>- script to convert bd files to datframe used to create the new har gene lists. comes with bed_to_df.sh bash file.

**note**- we can also see bash files that are accompanied to the scripts in the pipeline and outside of it and runs them. (we can understand what scripts the .sh runs by its name of by lookimng inside the file). 

#### more data to know about
<span style="color:red">**cloud location- more_data directory**</span><br>
there are more data files that are important in the project:
1. all_data_ens_above_100 and all_data_genes_above100- list of genes in the brain data that expressed in more than 100 cells. needed for the deseq_fisher analysis.
2. cortex_ens_above_100 and cortex_over_100_genes- list of genes in the cortex data that expressed in more than 100 cells. needed for the deseq_fisher analysis.
3. gencode.v43.annotation.gtf - gtf annotaion (used for extracting the new har gene lists)
4. genes.bed, genes_with5kb.bed- bed files with gene annotations and gens+5kb upstream the tss used for extracting the new har gene lists.
5. ges_equation-figure of the ges_equation
6. all interacrtion files- the base to create the new har lists
7. exp_df- tables of hsg and har (the old lists) expression in the data divided to cell types including mean exp and number of cells in eachcell type
   



   






















