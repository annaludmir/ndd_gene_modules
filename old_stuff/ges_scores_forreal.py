import scanpy as sc
import numpy as np
import pandas as pd
import sys
import time
from statsmodels.stats.multitest import multipletests
print('imports done')

col = sys.argv[1].strip('[]').split(',')
condition = sys.argv[2].strip('[]').split(',')

print(f'starting running')
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
#v2_adata=adata_cortex[adata_cortex.obs['Chemistry']=='v2',:]
#print('data copied to v2 array')
v3_adata=adata_cortex[adata_cortex.obs['Chemistry']=='v3',:]
print('data copied to v3 array')
prepro_time=time.time()
print(f"preprocessig and data separation done in {time.time()-prepro_time:.2f} seconds")

def calculate_ges(adata,condition_col, target_cell_type, chemistry):
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
    other_cell_types = [ct for ct in cell_types if ct != target_cell_type]

    target_idx = adata.obs[condition_col] == target_cell_type
    filtered_data=adata[target_idx]
    # Mean expression in the target cell type
    mean_target_expr = filtered_data.X.mean(axis=0).A1
    
    # Weighted mean expression in other cell types
    total_cells = adata.shape[0]
    weighted_sum = np.zeros(adata.shape[1])
    per_exp_other={}
    for ct in other_cell_types:
        ct_idx = adata.obs[condition_col] == ct
        ct_fraction = ct_idx.sum() / total_cells
        ct_mean_expr = adata.X[ct_idx].mean(axis=0).A1
        weighted_sum += ct_fraction * ct_mean_expr
        if len(other_cell_types)>1:
            per_exp_other[f"per_exp_{ct}"]= (adata.X[ct_idx]>0).mean(axis=0).A1

    # Calculate GES scores
    epsilon=1e-10
    ges_scores = mean_target_expr / (weighted_sum+epsilon)

    # Recompute masks for filtered data -other
    filtered_other_data = adata[~target_idx]
    # Calculate mean expression
    mean_expression_other = filtered_other_data.X.mean(axis=0).A1
    # Calculate percentages of expressing cells
    percent_expressed_target = (filtered_data.X >0).mean(axis=0).A1
    percent_expressed_other = (filtered_other_data.X > 0).mean(axis=0).A1

    # Create DataFrame with results
    results_df = pd.DataFrame({
        "gene": adata.var_names,
        "ges_score": ges_scores,
        "mean_expression_target": mean_target_expr,
        "mean_expression_other": mean_expression_other,
        "per_expressed_target":percent_expressed_target ,
        "per_expressed_other":percent_expressed_other
    })
    
    if len(per_exp_other)>0:
        per_expressed_df=pd.DataFrame.from_dict(per_exp_other)
        per_expressed_df["gene"]=adata.var_names
        results_df=results_df.merge(per_expressed_df, on= "gene")

    print(f"finished_ges_caculations in {time.time()-ges_time:.2f} seconds")

    return results_df

def calculate_ges_with_permutations(adata,condition_col,target_cell_type, chemistry, expression_threshold=0.05, n_permutations=500):
       
    # Subset the data for the target cell type
    target_mask = adata.obs[condition_col] == target_cell_type
    target_data = adata[target_mask]  # Subset for cells in the target cell type

    # Calculate the number of cells expressing each gene in the target cell type
    gene_expression_num = (target_data.X > 0).sum(axis=0).A1 

    # Filter genes with at least `expression_threshold` proportion of expression
    expressed_genes_mask = gene_expression_num >= (target_data.shape[0]*expression_threshold)
    filtered_adata = adata[:, expressed_genes_mask] # Subset genes that expressed in more than 5% of the target population

    actual_results = calculate_ges(filtered_adata, condition_col, target_cell_type, chemistry)
    actual_results.to_csv(f"ges_only_{chemistry}_{condition_col}_{target_cell_type}.csv")
    actual_scores = actual_results["ges_score"].values
    # Initialize array to store permuted scores
    
    permuted_scores = np.zeros((n_permutations, len(actual_scores)))
    p_values=np.zeros(len(actual_scores))
    
    # Perform permutations
    for i in range(n_permutations):
        print("permutation ", i+1)
       # Shuffle the condition column
        shuffled_conditions = np.random.permutation(filtered_adata.obs[condition_col])
        filtered_adata.obs[condition_col] = shuffled_conditions
        
        # Recalculate GES scores for the shuffled data
        permuted_results = calculate_ges(filtered_adata, condition_col, target_cell_type,chemistry)
        permuted_scores[i] = permuted_results["ges_score"].values
        p_values+=np.array((actual_scores<permuted_scores[i]).astype(int).tolist())
            
    # Add p-values to results
    actual_results["p_values"] = p_values/n_permutations
    #subset the data for gens with more than 5% expression in the cell population
    actual_results=actual_results.loc[actual_results["per_expressed_target"]>=expression_threshold,:]
    actual_results["adj_p_val"]= multipletests(actual_results["p_values"], method='fdr_bh')[1]
    
    return actual_results

for i in col:
    cell_types=adata_cortex.obs[i].unique().tolist()
    for j in condition:
        if j in cell_types:
            print(f"starting ges calculations on {i} : {j}")
            v3_time=time.time()
            ges_results_v3 = calculate_ges_with_permutations(v3_adata,i,j,'v3')
            print(f"finished calcultaing ges for v3, total time: {time.time()-v3_time:.2f} seconds")
            ges_results_v3.to_csv(f"ges_spec_v3_{j}.csv")
            print('saved v3 results to csv')
            #v2_time=time.time()
            #ges_results_v2 = calculate_ges_with_permutations(v2_adata,i,j,'v2')
            #print(f"finished calculating ges for v2, total time: {time.time()-v2_time:.2f} seconds")
            #ges_results_v2.to_csv(f"ges_spec_v2_{j}.csv")
            #print ('saved v2 results to csv')

print('DONE')
