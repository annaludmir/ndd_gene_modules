import scanpy as sc
import numpy as np
import pandas as pd
import sys
import time
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
print('imports done')

col = sys.argv[1]
condition = sys.argv[2]

print(f'starting running for {col} : {condition}')
start_time=time.time()

#upload the data
adata_cortex=sc.read_h5ad("/scratch200/reutj/data/updated_cortex_data_hg38.h5ad")

print(f"data loaded in {time.time()-start_time:.2f} seconds.")

#data_preprocessing
prepro_time=time.time()
adata_cortex.layers["counts"] = adata_cortex.X.copy()
#filter cells with less than 200 genes expressd anf=d genes that express in less than 3 cells
sc.pp.filter_cells(adata_cortex, min_genes=200)
sc.pp.filter_genes(adata_cortex, min_cells=3)
sc.pp.normalize_total(adata_cortex)
sc.pp.log1p(adata_cortex)
# Store the full data set in 'raw' as log-normalised data for statistical testing
adata_cortex.raw = adata_cortex
print('preprocessing done')

#turn the normalized count matrix to array for better visualization and divide based on v2 and v3 chemistry
v2_adata=adata_cortex[adata_cortex.obs['Chemistry']=='v2',:]
print('data copied to v2 array')
v3_adata=adata_cortex[adata_cortex.obs['Chemistry']=='v3',:]
print('data copied to v3 array')
prepro_time=time.time()
print(f"preprocessig and data separation done in {time.time()-prepro_time:.2f} seconds")


def manwhitney_test_per_gene(adata, target_cell_type, condition_col):
    """
    Perform Wilcoxon rank-sum test for each gene between the target cell type
    and all other cell types individually.

    Parameters:
        adata: AnnData object containing single-cell data.
        target_cell_type: The cell type to compare against others.
        condition_col: Column in `adata.obs` specifying cell types.

    Returns:
        pd.DataFrame with columns: 'gene', 'comparison_cell_type', 'p_value'
    """
    genes=adata.var_names.tolist()
    cell_types = adata.obs[condition_col].unique()
    cell_types = [ct for ct in cell_types if ct != target_cell_type]

    results_df=[]
    gene_time=time.time()
    # Iterate through genes
    for i,gene in enumerate(genes):
        if i%10==0:
            print(f"finished {i} genes caculations in {time.time()-gene_time:.2f} seconds")
            gene_time=time.time()
        target_expr = adata[:, gene].X[adata.obs[condition_col] == target_cell_type].toarray().flatten()
        results_df.append({"gene":gene})
        for other_cell_type in cell_types:
            other_expr = adata[:, gene].X[adata.obs[condition_col] == other_cell_type].toarray().flatten()

            # Perform Wilcoxon rank-sum test
            _, p_value = mannwhitneyu(target_expr, other_expr, alternative='two-sided')
            
            # Store the result
            results_df.append({
                "gene": gene,
                "comparison_cell_type": other_cell_type,
                "p_value": p_value,
            })

    # Convert results to DataFrame
    print("finished manwhitney calculations")
    results_df = pd.DataFrame(results)
    results_df = df.pivot_table(index=['gene'], columns='comparison_cell_type', values='p_value')
    for cell_type in cell_types:
        print("performing multiple correction ", cell_type) 
        results_df[f"adj_p_val_{cell_type}"]= multipletests(results_df[cell_type], method='fdr_bh')[1]

    return results_df

def calculate_ges(adata,condition_col, target_cell_type, chemistry, num_permutations=1000):
    """
    Calculate the Generalized Expression Specificity (GES) score and perform permutation testing.
    
    Parameters:
    - adata: AnnData object (single-cell data)
    - target_cell_type: str, cell type to calculate specificity for
    - num_permutations: int, number of permutations for p-value calculation
    
    Returns:
    - result_df: DataFrame with gene names, GES, p-values, and adjusted p-values
    """
    ges_time=time.time()
    cell_types = adata.obs[condition_col].unique()

    # Subset the data for the target cell type
    target_mask = adata.obs[condition_col] == target_cell_type
    target_data = adata[target_mask]  # Subset for cells in the target cell type

    # Calculate the percentage of cells expressing each gene in the target cell type
    gene_expression_percentage = (target_data.X > 0).mean(axis=0).A1

    # Filter genes with at least `expression_threshold` proportion of expression
    expressed_genes_mask = gene_expression_percentage >= 0.05
    filtered_adata = adata[:, expressed_genes_mask]  # Subset genes
    print("filtered data shape ", filtered_adata.shape)

    # Update the variable names to include only expressed genes
    filtered_gene_names = filtered_adata.var_names

    # Recompute masks for filtered data
    filtered_target_data = filtered_adata[target_mask]
    filtered_other_data = filtered_adata[~target_mask]

    # Calculate mean expression
    mean_expression_cond = filtered_target_data.X.mean(axis=0).A1
    mean_expression_other = filtered_other_data.X.mean(axis=0).A1

    # Calculate percentages of expressing cells
    percent_expressed_target = (filtered_target_data.X > 0).mean(axis=0).A1
    percent_expressed_other = (filtered_other_data.X > 0).mean(axis=0).A1

    # Calculate GES
    epsilon=1e-10
    ges = mean_expression_cond / (mean_expression_other + epsilon)

    # Combine scores into a DataFrame
    ges_df = pd.DataFrame({
        "gene": filtered_gene_names,
        "GES": ges,
        "mean_expression_target": mean_expression_cond,
        "mean_expression_other": mean_expression_other,
        "per_expressed_target": percent_expressed_target,
        "per_expressed_other": percent_expressed_other
    })

    ges_df = ges_df.sort_values(by="GES", ascending=False)
    ges_df.to_csv(f"ges_only_{chemistry}_{condition_col}_{target_cell_type}.csv")
    print(f"calculated ges scores in {time.time()-ges_time:.2f} seconds")

    #ges_df['p_value']=p_values
    print('starting manwhiteney test')
    wilc_time=time.time()
    p_values= manwhitney_test_per_gene(filtered_adata,target_cell_type,condition_col)
    p_values['pass_p_value']=[all(p_values.iloc[i,5:]<0.05) for i in range(p_values.shape[0])]
    print(f"finished p_values testing in {time.time()-wilc_time:.2f} seconds")
    ges_df=pd.merge(ges_df,p_values,on="gene")
    
    return ges_df


# Usage example:
ges_results_v3 = calculate_ges(v3_adata,col,condition,'v3')
print(f"finished calcultaing ges for v3, total time: {time.time()-start_time()}")
ges_results_v3.to_csv(f"ges_spec_v3_{col}_{condition}.csv")
print('saved v3 results to csv')

v2_time=time.time()
ges_results_v2 = calculate_ges(v2_adata,col,condition,'v2')
print(f"finished calculating ges for v2, total time: {time.time()-v2_time()}")
ges_results_v2.to_csv(f"ges_spec_v2_{col}_{condition}.csv")
print ('saved v2 results to csv')
print('DONE')

