import scanpy as sc
import pandas as pd
import numpy as np
import sys

col = sys.argv[1]
condition = sys.argv[2]
inp_k=sys.argv[3]
#from scipy import stats

print(f'starting running for {col} : {condition} with k= {inp_k}') 

print('imports done')

#upload the data
adata_cortex=sc.read_h5ad("/scratch200/reutj/data/updated_cortex_data_hg38.h5ad")

print('data loaded')

#data_preprocessing
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
v2_array=v2_adata.X.toarray()
print('data copied to v2 array')
v3_adata=adata_cortex[adata_cortex.obs['Chemistry']=='v3',:]
v3_array=v3_adata.X.toarray()
print('data copied to v3 array')

gene_sym_adata=adata_cortex.var_names.to_list()

#specficity matrix
def get_specificity_adjusted(adata,counts_array,obs_col,condition,k,c,alpha):
    ##this function takes the data with count matrix and gets speceficty scores + other important measures
    ## depending on certain condition. spec score is sombined from exp_ration_mean changes betwen the
    ## condition and other cells , ratio chnages of the percent of cells expressed from the condition or other cells.
    ## we used here sigmoid function tocontrol the effect of percent_other in low valies and **2 to percent condition
    ## to give it more weight, also gene expression in the condition cells is weighthed here using the log function
    specificity={"gene":[],"mean_exp":[],"exp_ratio":[],"per_exp_con":[],"per_exp_other":[],"spec_score":[]}
    n_condition=sum(adata.obs[obs_col]==condition)
    n_other=len(adata.obs)-n_condition
    count=0
    for gene in gene_sym_adata:
        gene_exp_con=counts_array[adata.obs[obs_col]==condition,adata.var_names==gene]
        count+=1
        if count%1000==0:
            print(count) # to keep track of progress
        if (np.count_nonzero(gene_exp_con))>=(n_condition/10): #maybe avoid that and consider all?
            #limit to cellswith some expression - more than 1/10 cells based on the population size expressing
            #to avoid unnecessary genes
            gene_exp_o=counts_array[adata.obs[obs_col]!=condition,adata.var_names==gene]
            mean_exp=np.mean(counts_array[:,adata.var_names==gene])
            mean_exp_con=(np.sum(gene_exp_con))/n_condition
            mean_exp_o=(np.sum(gene_exp_o))/n_other
            per_exp_con=(np.count_nonzero(gene_exp_con))/n_condition
            per_exp_o=(np.count_nonzero(gene_exp_o))/n_other
            adjusted_percent_other = 1 / (1 + np.exp(-k * (per_exp_o - c)))
            specificity_score= (mean_exp_con/mean_exp_o)*((per_exp_con**alpha)/(adjusted_percent_other))*np.log1p(mean_exp_con)
            specificity["gene"].append(gene)
            specificity["mean_exp"].append(mean_exp)
            specificity["exp_ratio"].append(mean_exp_con/mean_exp_o)
            specificity["per_exp_con"].append(per_exp_con)
            specificity["per_exp_other"].append(per_exp_o)
            specificity["spec_score"].append(specificity_score)
    return specificity

def get_per_exp_o(adata,counts_array,obs_col,condition):
    count=0
    per_exp_o=[]
    n_condition=sum(adata.obs[obs_col]==condition)
    print("n_condition: ",n_condition)
    n_other=len(adata.obs)-n_condition
    for gene in gene_sym_adata:
        num_cells_exp=np.count_nonzero(counts_array[:,adata.var_names==gene])
        count+=1
        if count%1000==0:
            print(count) # to keep track of progress
        if num_cells_exp >(n_condition/10):
            gene_exp_o=counts_array[adata.obs[obs_col]!=condition,adata.var_names==gene]
            per_exp_o.append((np.count_nonzero(gene_exp_o))/n_other)
    return np.mean(per_exp_o)

print('extracting per_exp_other v2')
v2_per_o=get_per_exp_o(v2_adata,v2_array,col,condition)
print(v2_per_o)
print('extracting per_exp_other v3')
v3_per_o=get_per_exp_o(v3_adata,v3_array,col,condition)
print(v3_per_o)

print('getting speceficity scores') 
spec_v2=get_specificity_adjusted(v2_adata,v2_array,col,condition,int(inp_k),v2_per_o,2)
print('speceficty obtained for v2')
spec_v3=get_specificity_adjusted(v3_adata,v3_array,col,condition,int(inp_k),v3_per_o,2)
print('speceficty obtained for v3')

#create dataframe
spec_v2_df=pd.DataFrame(spec_v2)
spec_v2_df["spec_score"]=round(spec_v2_df["spec_score"],2)
spec_v2_df=spec_v2_df.sort_values(by='spec_score', ascending=False).reset_index()
spec_v2_df.to_csv(f"adj_spec_v2_{col}_{condition}_{inp_k}.csv")
print('saved v2 spec to csv')

spec_v3_df=pd.DataFrame(spec_v3)
spec_v3_df["spec_score"]=round(spec_v3_df["spec_score"],2)
spec_v3_df=spec_v3_df.sort_values(by='spec_score', ascending=False).reset_index()
spec_v3_df.to_csv(f"adj_spec_v3_{col}_{condition}_{inp_k}.csv")
print('saved v3 spec to csv')

print('DONE')



