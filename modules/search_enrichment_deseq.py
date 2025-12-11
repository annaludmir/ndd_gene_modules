import scanpy as sc
import numpy as np
import pandas as pd
import sys
import os
import gseapy as gp
import statsmodels.stats.multitest as smm

#gsea function..
def run_gsea(gmt_file,conditions,out_folder):
    condition_cell,gsea_result,num_lead,per_lead,NES=[],[],[],[],[]
    for i in conditions:
        print(f"testing condition {i} on {gmt_file} list")
        condition_cell.append(i)
        #sort ges_scores_table
        diff_scores= pd.read_csv(f"/scratch200/reutj/data/deseq2/deseq_results_{i}.csv")
        diff_scores_up=diff_scores[diff_scores.log2FoldChange>0.58]#fc>1.5
        diff_scores_up=diff_scores_up.sort_values(by='padj')
        #prepare table for gsea
        diff_scores_en=diff_scores_up.stat
        diff_scores_en.index=diff_scores_up.Gene
        #prepare paths
        outpath=out_folder+gmt_name+"/"+i+"/"
        print("outpath: ",outpath)
        outraw=outpath+"gsea_raw.csv"
        outgsea=outpath+"gsea_results.csv"
        #gsea    
        try:
            gsea = gp.prerank(rnk=diff_scores_en,
               gene_sets=gmt,  
               outdir=outraw,  # Output directory
               min_size=1,
               max_size=2000,
               permutation_type='geneset',
               permutation_num=10000,
               format='pdf')
            gsea_result.append(gsea.res2d.loc[0,"FDR q-val"])
            num_lead.append(gsea.res2d.loc[0,"Tag %"])
            per_lead.append(gsea.res2d.loc[0,"Gene %"])
            NES.append(gsea.res2d.loc[0,"NES"])
            print(gsea.res2d.loc[0,"FDR q-val"])
            gsea.res2d.to_csv(outgsea)
            if gsea.res2d.loc[0,"FDR q-val"]<0.05:
                print(f"hsg result for: {i} is significant!")
            else:
                print(f"hsg result for: {i} is non significant")
       except:
           print(f"could not find enrichment in this gene set for condition {i}")

    results=pd.DataFrame(data={"condition":condition_cell,"num_genes_in_lead":num_lead,"%genes_in_lead":per_lead,"nes_score":NES,"gsea_pval":gsea_result})
    reject, pvals_corr = smm.multipletests(gsea_result, method='fdr_bh')[:2]
    results["p-adj"]= pvals_corr
    final_outpath=out_folder+gmt_name+"/"+"full_results.csv"
    print("final outpath: ", final outpath)
    results.to_csv(final_outpath)

if __name__=="__main__":
    #add ges_results so the results will go to the specifies folder
    out_folder=out_folder+"/"+"ges_results"+"/"
    #get clean gmt names out of the gme folders
    gmt_names=[file.split("/")[-1][:-4] for file in os.listdir(gmt_folder) if file.endswith(".gmt")]
    for gmt_file in gmt_names:
        run_gsea(gmt_file,conditions,out_folder)
    print('DONE')
