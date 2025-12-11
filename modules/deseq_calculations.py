import scanpy as sc
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import fisher_exact
import sys
import time
import os
import statsmodels.stats.multitest as smm

#fisher exact function for computing fisher of hsg list interction with sig deseq2 gene table
def fisher_test_hsg(hsg_genes,sig_table,n_hsg,n_genes):
    hsg_intersection=list(set(sig_table.index) & set(hsg_genes))
    table=np.array([[len(hsg_intersection),n_hsg-len(hsg_intersection)],
                    [len(sig_table)-len(hsg_intersection),n_genes-len(sig_table)-len(hsg_intersection)]])
    res = fisher_exact(table, alternative='two-sided')
    return (hsg_intersection , res.pvalue)

#start the process
def run_fisher(psb_data_path,out_path,hsg_gene_in_data,n_genes):
    psb_time=time.time()
    print('reading psb data from: ', psb_data_path)
    psb_data=sc.read_h5ad(psb_data_path)
    psb_name=psb_data_path.split("/")[-1][:-5]
    out_path2=out_path+psb_name+"/"
    os.makedirs(out_path2,exist_ok=True)#make sure path exists..
    counts=pd.DataFrame(psb_data.X,columns=psb_data.var_names)
    comparison,len_intersection,intersecting_genes,p_value=[],[],[],[]
    for col in psb_data.obs.columns[2:]:
        col_time=time.time()
        print(f'starting deseq analysis on {col}')
        dds=DeseqDataSet(counts=counts,metadata=psb_data.obs,design_factors=['sample',col],refit_cooks=True)
        dds.deseq2()
        print(f'finished_creating_dds_obsect_for_{col} in {time.time()-col_time:.2f} seconds ')
        #now create all the posibble combinations we would look at all the posibble pairs
        criteria_to_check=psb_data.obs[col].unique().tolist()
        pairs=[(a, b) for idx, a in enumerate(criteria_to_check) for b in criteria_to_check[idx + 1:]]
        for i in pairs:
            pair_time=time.time()
            print(f'starting analysis for comparison: {i}')
            just_up=False#check only upregulated or also downregulated genes
            if i[0]=="other":#if there is condition vs other we would just want to test that direction
                ds = DeseqStats(dds,contrast=(col,i[1],"other"),alpha=0.05,cooks_filter=True,independent_filter=True)
                comparison.append(f"{i[1]}_vs_{i[0]}")
                out_path3=os.path.join(out_path2,f"{i[1]}_vs_{i[0]}.csv")
                just_up=True
            elif i[1]=="other":
                ds = DeseqStats(dds,contrast=(col,i[0],"other"),alpha=0.05,cooks_filter=True,independent_filter=True)
                comparison.append(f"{i[0]}_vs_{i[1]}")
                out_path3=os.path.join(out_path2,f"{i[0]}_vs_{i[1]}.csv")
                just_up=True
            else:
                ds = DeseqStats(dds,contrast=(col,i[0],i[1]),alpha=0.05,cooks_filter=True,independent_filter=True)
                comparison.append(f"{i[0]}_vs_{i[1]}")
                out_path3=os.path.join(out_path2,f"{i[0]}_vs{i[1]}.csv")
            #summerize scores
            print('summarizing scores')
            ds.summary()
            de=ds.results_df
            de=de.sort_values(by='padj')
            #save to file
            de.to_csv(out_path3)
            print(f'saved results to {out_path3}')
            de_sig_up=de[(de.padj<0.05)&(de.log2FoldChange>1.5)&(de.baseMean>5)]
            fisher_up_hsg=fisher_test_hsg(hsg_gene_in_data,de_sig_up,n_hsg=len(hsg_gene_in_data),n_genes=n_genes)
            print(f'printing fisher exact results for {i[0]} vs {i[1]} with hsg genes')
            intersecting_genes.append(fisher_up_hsg[0])
            len_intersection.append(len(fisher_up_hsg[0]))
            p_value.append(float(fisher_up_hsg[1]))
            print('num of sig genes', len(de_sig_up))
            print(fisher_up_hsg)
            if not just_up:
                de_sig_down=de[(de.padj<0.05)&(de.log2FoldChange<-1.5)&(de.baseMean>5)]
                fisher_down_hsg=fisher_test_hsg(hsg_gene_in_data,de_sig_down,n_hsg=len(hsg_gene_in_data),n_genes=n_genes)
                print(f'printing fisher exact results for {i[1]} vs {i[0]} with hsg genes')
                comparison.append(f"{i[1]}_vs_{i[0]}")
                len_intersection.append(len(fisher_down_hsg[0]))
                intersecting_genes.append(fisher_down_hsg[0])
                p_value.append(float(fisher_down_hsg[1]))
                print('num of sig genes', len(de_sig_down))
                print(fisher_down_hsg)
            print(f'finished analysis for pair {i} in {time.time()-pair_time:.2f} seconds')
    print(f'finished all fisher_run for psb in {time.time()-psb_time:.2f} seconds')
    results=pd.DataFrame(data={'comparison':comparison,'intersect_genes':intersecting_genes,'num_intersect_genes':len_intersection,'pval':p_value})
    reject, pvals_corr = smm.multipletests(p_value, method='fdr_bh')[:2]
    results["p-adj"]= pvals_corr
    final_outpath=os.path.join(out_path2, "full_results.csv")
    print("final outpath: ", final_outpath)
    results.to_csv(final_outpath)
         

def main():
    start_time=time.time()
    #gene_names_list=
    print('uploading hsg file: ', hsg_file)
    hsg=pd.read_csv(hsg_file)
    hsg_list_name=hsg_file.split("/")[-1][:-4]
    # Open the file in read mode
    gene_list_names=[]
    n_genes=0
    with open(gene_names, 'r') as file:
        for line in file:
            gene_list_names.append(line.strip())
    if data_type=="data_all":
        hsg_gene_in_data=[i for i in hsg["ens"] if i in gene_list_names]#use ens here cause deseq results are in ens
        n_genes=28587
    else:#cortex
        hsg_gene_in_data=[i for i in hsg["gene"] if i in gene_list_names]
        n_genes=23883
    #add deseq_results so the results will go to the specifies folder
    out_folder=out_path+"deseq_fisher_results/"+data_type+"/"+hsg_list_name+"/"
    #get clean gmt names out of the gme folders.. change!
    psb_data_folder_dt=psb_data_folder+data_type+"/"
    psb_datas=[file for file in os.listdir(psb_data_folder_dt) if file.endswith(".h5ad")]
    for psb_data in psb_datas:
        print("psb_data: ",psb_data)
        psb_data_path=psb_data_folder_dt+psb_data
        run_fisher(psb_data_path,out_folder,hsg_gene_in_data,n_genes)
    print(f'finished all fisher_run in {time.time()-start_time:.2f} seconds')
    print('DONE')

if __name__=="__main__":
    main()

        
        
