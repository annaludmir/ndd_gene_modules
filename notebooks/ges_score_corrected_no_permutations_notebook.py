import marimo

__generated_with = "0.17.8"
app = marimo.App()


@app.cell
def _():
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import sys
    import time
    from statsmodels.stats.multitest import multipletests
    print('imports done')

    col = ['region_general']
    condition = ['radialglia','Cerebellum','Diencephalon','Fibroblast','Glioblast','Forebrain_general','hindbrain_general','Immune','Midbrain_general','Neural_crest','IPC','Neuroblast','Neuron','Oligo','Forebrain','Hindbrain','Placodes','Pons','Medulla','Telencephalon','Vascular']
    data_type = "data_all"
    output_folder = "/miridan-data/annaludmir/data/ges_results/"


    print(f'starting running')
    start_time=time.time()


    if data_type=="cortex":
        adata=sc.read_h5ad("/miridan-data/annaludmir/data/updated_cortex_data_hg38.h5ad")
        print('data uploaded')
        adata=adata[adata.obs['Chemistry']=='v3',:]#extract only v3
        print('data copied to v3 array')
    else:#all data
        adata=sc.read_h5ad("/miridan-data/annaludmir/data/updated_data_all.h5ad")
        # adata=adata[~adata.obs.Region.isin(["b'Head'","b'Brain'"])]#substract head and brain
        print('data uploaded')

    print(f"data loaded in {time.time()-start_time:.2f} seconds.")

    #data_preprocessing
    prepro_time=time.time()
    #adata.layers["counts"] = adata.X.copy()
    #filter cells with less than 200 genes expressd anf=d genes that express in less than 3 cells
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    # Store the full data set in 'raw' as log-normalised data for statistical testing
    #adata.raw = adata
    print(f"preprocessing done in {time.time()-prepro_time:.2f} seconds")


    def calculate_ges(adata,condition_col, target_cell_type, chemistry,data_type):
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
        print(type(adata.X))
        total_cells = adata.shape[0]
        target_idx = adata.obs[condition_col] == target_cell_type
        target_fraction = target_idx.sum() / total_cells
        filtered_data=adata[target_idx]
        # Mean expression in the target cell type
        mean_target_expr = filtered_data.X.mean(axis=0).A1
        weighted_target_expr=(1-target_fraction)*mean_target_expr

        # Weighted mean expression in other cell types
        cell_types = adata.obs[condition_col].unique()
        other_cell_types = [ct for ct in cell_types if ct != target_cell_type]
        weighted_sum = np.zeros(adata.shape[1])
        per_exp_other={}
        for ct in other_cell_types:
            print(ct)
            ct_idx=adata.obs[condition_col]==ct
            ct_idx=np.array(ct_idx)
            ct_fraction = ct_idx.sum() / total_cells
            ct_mean_expr = adata[ct_idx].X.mean(axis=0).A1
            weighted_sum += ct_fraction * ct_mean_expr
            if len(other_cell_types)>1:
                per_exp_other[f"per_exp_{ct}"]= (adata.X[ct_idx]>0).mean(axis=0).A1

        # Calculate GES scores
        epsilon=1e-10
        ges_scores = weighted_target_expr / (weighted_sum+epsilon)

        # Recompute masks for filtered data -other
        filtered_other_data = adata[~target_idx]
        # Calculate mean expression
        mean_expression_other = filtered_other_data.X.mean(axis=0).A1
        # Calculate percentages of expressing cells
        percent_expressed_target = (filtered_data.X >0).mean(axis=0).A1
        percent_expressed_other = (filtered_other_data.X > 0).mean(axis=0).A1

        if data_type=="cortex":
            genes=adata.var_names
        else: #in the big data gene symbols are in var_names["accession"]
            genes=adata.var["Gene"]

        # Create DataFrame with results
        results_df = pd.DataFrame({
            "gene": genes,
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

    def calculate_ges_with_permutations(adata,condition_col,target_cell_type, chemistry, data_type, expression_threshold=0.05, n_permutations=500):

        # Subset the data for the target cell type
        target_mask = adata.obs[condition_col] == target_cell_type
        target_data = adata[target_mask]  # Subset for cells in the target cell type

        # Calculate the number of cells expressing each gene in the target cell type
        gene_expression_num = (target_data.X > 0).sum(axis=0).A1 

        # Filter genes with at least `expression_threshold` proportion of expression
        expressed_genes_mask = gene_expression_num >= (target_data.shape[0]*expression_threshold)
        filtered_adata = adata[:, expressed_genes_mask] # Subset genes that expressed in more than 5% of the target population

        actual_results = calculate_ges(filtered_adata, condition_col, target_cell_type, chemistry, data_type)
        actual_results.to_csv(f"ges_only_{chemistry}_{target_cell_type}.csv")
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
            permuted_results = calculate_ges(filtered_adata, condition_col, target_cell_type,chemistry,data_type)
            permuted_scores[i] = permuted_results["ges_score"].values
            p_values+=np.array((actual_scores<permuted_scores[i]).astype(int).tolist())

        # Add p-values to results
        actual_results["p_values"] = p_values/n_permutations
        #subset the data for gens with more than 5% expression in the cell population
        actual_results=actual_results.loc[actual_results["per_expressed_target"]>=expression_threshold,:]
        actual_results["adj_p_val"]= multipletests(actual_results["p_values"], method='fdr_bh')

        return actual_results

    for i in col:
        cell_types=adata.obs[i].unique().tolist()
        cell_types = [c.decode('utf-8') if isinstance(c, bytes) else c for c in cell_types]
        print(cell_types)
        for j in condition:
            print(j)
            if j in cell_types:
                print(f"starting ges calculations on {i} : {j}")

                target_mask = adata.obs[i] == j
                target_data = adata[target_mask]  # Subset for cells in the target cell type

                # Calculate the number of cells expressing each gene in the target cell type
                gene_expression_num = (target_data.X > 0).sum(axis=0).A1 

                # Filter genes with at least `expression_threshold` proportion of expression
                expressed_genes_mask = gene_expression_num >= (target_data.shape[0]*0.05)
                filtered_adata = adata[:, expressed_genes_mask] # Subset genes that expressed in more than 5% of the target population

                actual_results = calculate_ges(filtered_adata,i,j,"v3",data_type)
                actual_results = actual_results.sort_values("ges_score", ascending=False)
                if data_type=="cortex":
                    actual_results.to_csv(f"{output_folder}ges_spec_v3_{j}.csv")
                else:
                    actual_results.to_csv(f"{output_folder}ges_spec_big_{j[2:-1]}.csv")
                print('saved results to csv')

                #permutations
                #v3_time=time.time()
                #ges_results_v3 = calculate_ges_with_permutations(v3_adata,i,j,'v3',data_type)
                #print(f"finished calcultaing ges for v3, total time: {time.time()-v3_time:.2f} seconds")
                #ges_results_v3.to_csv(f"ges_spec_v3_{j}.csv")


    print('DONE')
    return


if __name__ == "__main__":
    app.run()
