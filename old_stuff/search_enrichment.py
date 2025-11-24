import scanpy as sc
import numpy as np
import pandas as pd
import sys
import os
import gseapy as gp
import statsmodels.stats.multitest as smm

hsg_list=sys.argv[1]
gmt_folder=sys.argv[2]
conditions=sys.argv[3].strip('[]').split(',')
out_folder=sys.argv[4]

#gsea function..
def run_gsea(gmt_file,conditions,out_folder):
    condition_cell,gsea_result,num_lead,per_lead,NES=[],[],[],[],[]
    gmt_name=gmt_file.split("/")[-1][:-4]
    for i in conditions:
        print(f"testing condition {i} on {gmt_name} list")
        if os.path.exists(f"/scratch200/reutj/data/spec_score_tables/ges_spec_v3_updated_{i}.csv"):
            ges_df=pd.read_csv(f"/scratch200/reutj/data/spec_score_tables/ges_spec_v3_updated_{i}.csv")
        elif os.path.exists(f"/scratch200/reutj/data/spec_score_tables/ges_spec_v3_{i}.csv"):
            ges_df=pd.read_csv(f"/scratch200/reutj/data/spec_score_tables/ges_spec_v3_{i}.csv")
        else:
            print(f"condition {i} is not a valid condition name")
            sys.exit()
        #take only with ges_score above 1
        ges_df=ges_df[ges_df.ges_score>=1]
        print("number of genes with ges>1", len(ges_df))
        #prepare table for gsea
        ges_scores_en=ges_df.ges_score
        ges_scores_en.index=ges_df.gene
        #prepare paths
        outpath=out_folder+gmt_name+"/"+i+"/"
        print("outpath: ",outpath)
        outraw=outpath+"gsea_raw"
        outgsea=outpath+"gsea_results.csv"
        #gsea    
        try:
            gsea = gp.prerank(rnk=ges_scores_en,
               gene_sets=gmt_file,  
               outdir=outraw,  # Output directory
               min_size=5,
               seed=6,
               max_size=2500,
               permutation_type='geneset',
               permutation_num=10000,
               format='pdf')
            gsea_result.append(gsea.res2d.loc[0,"FDR q-val"])
            num_lead.append(gsea.res2d.loc[0,"Tag %"])
            per_lead.append(gsea.res2d.loc[0,"Gene %"])
            NES.append(gsea.res2d.loc[0,"NES"])
            condition_cell.append(i)
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
    print("final outpath: ", final_outpath)
    results.to_csv(final_outpath)

if __name__=="__main__":
    #add ges_results so the results will go to the specifies folder
    out_folder=out_folder+"ges_results_above1"+"/"
    gene_list_name=hsg_list.split("/")[-1][:-4]
    gmt_file=gmt_folder+f"{gene_list_name}.gmt"
    print(f"starting ges enrichment analysis on {gmt_file}")
    run_gsea(gmt_file,conditions,out_folder)
    print('DONE')
