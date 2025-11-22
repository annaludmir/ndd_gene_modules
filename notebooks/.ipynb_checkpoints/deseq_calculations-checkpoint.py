import scanpy as sc
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import fisher_exact
import sys
import time
print('imports done')

psb_data_path=sys.argv[1]
out_path=sys.argv[2]
start_time=time.time()

def fisher_test_hsg(hsg_genes,sig_table,n_hsg,n_genes=30251):
    hsg_intersection=list(set(sig_table.index) & set(hsg_genes))
    table=np.array([[len(hsg_intersection),n_hsg-len(hsg_intersection)],
                    [len(sig_table)-len(hsg_intersection),n_genes-len(sig_table)-len(hsg_intersection)]])
    res = fisher_exact(table, alternative='two-sided')
    return (hsg_intersection , res.pvalue)
    
#upload data
adata_time=time.time()
adata=sc.read_h5ad("/scratch200/reutj/data/updated_cortex_data_hg38_2.h5ad")
adata_v3=adata[adata.obs.Chemistry=='v3',:].copy()
print('data loaded')
#basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print(f'data loaded and preprocessed in {time.time()-adata_time:.2f} seconds')

#get gene symbols that appear in the data
gene_sym_adata=adata.var_names.to_list()
hsg=pd.read_csv("/scratch200/reutj/data/hsg_gene_lists/updated_human_specific_genes3.csv")
hsg_gene_in_data=[i for i in hsg["gene"] if i in gene_sym_adata]
har_genes=pd.read_excel("/scratch200/reutj/data/hsg_gene_lists/har_associated_genes.xlsx",header=1)
har_genes_list=[i for i in har_genes["HAR-associated genes.1"] if i in gene_sym_adata]

#start the process
print('reading psb data')
psb_data=sc.read_h5ad(psb_data_path)
psb_data.obs.columns = psb_data.obs.columns.str.replace("_", "-", regex=True)
counts=pd.DataFrame(psb_data.X,columns=psb_data.var_names)
for col in psb_data.obs.columns[1:]:
    col_time=time.time()
    print(f'starting deseq analysis on {col}')
    dds=DeseqDataSet(counts=counts,metadata=psb_data.obs,design_factors=['sample',col],refit_cooks=True)
    dds.deseq2()
    print(f'finished_creating_dds_obsect_for_{col} in {time.time()-col_time:.2f} seconds ')
    #now ro create all the posibble combinations we would look at all the posibble pairs
    criteria_to_check=psb_data.obs[col].unique().tolist()
    pairs=[(a, b) for idx, a in enumerate(criteria_to_check) for b in criteria_to_check[idx + 1:]]
    for i in pairs:
        pair_time=time.time()
        print(f'starting analysis for comparison: {i}')
        just_up=False#check only upregulated or also downregulated genes
        reverse_order=False
        if i[0]=="other":#if there is condition vs other we would just want to test that direction
            ds = DeseqStats(dds,contrast=(col,i[1],"other"),alpha=0.05,cooks_filter=True,independent_filter=True)
            reverse_order=True
            just_up=True
        elif i[1]=="other":
            ds = DeseqStats(dds,contrast=(col,i[0],"other"),alpha=0.05,cooks_filter=True,independent_filter=True)
            just_up=True
        else:
            ds = DeseqStats(dds,contrast=(col,i[0],i[1]),alpha=0.05,cooks_filter=True,independent_filter=True)
        #summerize scores
        print('summarizing scores')
        ds.summary()
        de=ds.results_df
        de=de.sort_values(by='padj')
        #save to file
        out_path2=out_path+f"deseq_results_{i[0]}_{i[1]}"
        de.to_csv(out_path2)
        print(f'saved results to {out_path2}')
        de_sig_up=de[(de.padj<0.05)&(de.log2FoldChange>1.5)&(de.baseMean>5)]
        print("length sig results: ", len(de_sig_up))
        fisher_up_hsg=fisher_test_hsg(hsg_gene_in_data,de_sig_up,n_hsg=443)
        if reverse_order:
            print(f'printing fisher exact results for {i[1]} vs {i[0]}')
        else:
            print(f'printing fisher exact results for {i[0]} vs {i[1]}')
        print("hsg genes results: ",fisher_up_hsg)
        fisher_up_har=fisher_test_hsg(har_genes_list,de_sig_up,n_hsg=710)
        print("har genes results: ", fisher_up_har)
        if not just_up:
            de_sig_down=de[(de.padj<1e-3)&(de.log2FoldChange<-1.5)&(de.baseMean>5)]
            fisher_down_hsg=fisher_test_hsg(hsg_gene_in_data,de_sig_down,n_hsg=443)
            print("length sig results: ", len(de_sig_down))
            print(f'printing fisher exact results for {i[1]} vs {i[0]}')
            print("hsg genes results: ",fisher_down_hsg)
            fisher_down_har=fisher_test_hsg(har_genes_list,de_sig_down,n_hsg=710)
            print("har genes results: ",fisher_down_har)
        print(f'finished analysis for pair {i} in {time.time()-pair_time:.2f} seconds')

print(f'finished all in {time.time()-start_time:.2f} seconds')
print('DONE')
        
        
