import scanpy as sc
import numpy as np
import pandas as pd
import sys
import os
import gseapy as gp
import statsmodels.stats.multitest as smm 
import matplotlib.pyplot as plt

gene_list=sys.argv[1]
gmt_folder=sys.argv[2]
conditions=sys.argv[3].strip('[]').split(',')
out_folder=sys.argv[4]
data_type=sys.argv[5]

def plot_enhanced_gsea(gsea_res, term, cell_type, output_dir):
    """
    Create an enhanced GSEA plot and save it to a file.
    
    Parameters:
    -----------
    gsea_res : GSEApy result object
        The result object from GSEApy prerank analysis
    term : str
        The gene set term to plot
    output_dir : str
        Directory to save the output figure
        
    Returns:
    --------
    fig_path : str
        Path to the saved figure
    """
    # Get the enrichment results for the term
    enrichment_results = gsea_res.results[term]
    
    # Create an enhanced custom GSEA plot
    plt.figure(figsize=(10, 8))
    
    # Plot the enrichment score profile with better styling
    plt.subplot(211)
    plt.plot(enrichment_results.get('RES'), color='forestgreen', linewidth=2)
    plt.axhline(y=0, color='gray', linestyle='--', alpha=0.7)
    plt.ylabel('Enrichment Score', fontsize=12)
    plt.title(f"{term} enrichment in {cell_type} ges scores\nNES={enrichment_results.get('nes'):.3f} FDR={enrichment_results.get('fdr'):.3e}", 
              fontsize=14, fontweight='bold')
    plt.grid(alpha=0.3)
    
    # Add vertical range to match GSEA standard plots
    max_es = max(abs(min(enrichment_results.get('RES'))), abs(max(enrichment_results.get('RES'))))
    plt.ylim(-max_es*1.1, max_es*1.1)
    
    # Plot the hits with a more visible representation
    plt.subplot(212)
    hit_indices = enrichment_results.get('hits')
    # Create a rug plot for the hits
    y = [1] * len(hit_indices)
    plt.plot([0, len(enrichment_results.get('RES'))], [1, 1], color='black', linewidth=0.5, alpha=0.5)
    plt.scatter(hit_indices, y, color='red', s=15, marker='|', alpha=0.8)
    
    # Add a heatmap-style gradient for the rank metric
    ax2 = plt.twinx()
    xs = range(len(gsea_res.ranking))
    ys = [0] * len(xs)
    ax2.scatter(xs, ys, c=gsea_res.ranking, cmap='coolwarm', s=10, marker='_', alpha=0.7)
    ax2.set_ylim(-0.5, 0.5)
    ax2.set_yticks([])
    
    plt.xlabel('Rank in Ordered Dataset', fontsize=12)
    plt.yticks([])
    
    # Add some labels - with adjusted position for "Hits" to avoid the frame
    plt.annotate('Hits', xy=(0.01, 0.95), xycoords='axes fraction', fontsize=10, 
                 fontweight='bold', color='red')
    plt.annotate('Ranked list', xy=(0.01, 0.02), xycoords='axes fraction', fontsize=10, 
                 fontweight='bold', color='blue')
    
    # Improve overall appearance
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.3)
    
    # Save figure with higher quality
    fig = plt.gcf()
    fig_path = os.path.join(output_dir, f"custom_gseaplot_{term}.png")
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return fig_path

#gsea function..
def run_gsea(gmt_file,conditions,out_folder,data_type):
    condition_cell,gsea_result,num_lead,per_lead,NES,lead_genes=[],[],[],[],[],[]
    gmt_name=gmt_file.split("/")[-1][:-4]
    for i in conditions:
        print(f"testing condition {i} on {gmt_name} list")
        ges_df=pd.read_csv(check_ges_path(i,data_type))
        print("uploaded ges file")
        #prepare table for gsea
        ges_df=ges_df[ges_df.ges_score>1]
        print(f"number of genes in {i} ges table with ges>1 is {ges_df.shape[0]}")
        #prepare ges list
        ges_scores_en=ges_df.ges_score
        ges_scores_en.index=ges_df.gene
        #prepare paths
        outpath=out_folder+gmt_name+"/"+i+"/"
        print("outpath: ",outpath)
        outraw=outpath+"gsea_raw"
        outgsea=outpath+"gsea_results.csv"
        #gsea                
        gsea_res = gp.prerank(rnk=ges_scores_en,
               gene_sets=gmt_file,  
               outdir=outraw,  # Output directory
               min_size=5,
               max_size=2500,
               seed=6,
               permutation_type='geneset',
               permutation_num=10000,
               format='pdf')
        #make a plot for the enrichment and save to file
        # Extract term (ensure it's a valid term)
        term = gsea_res.res2d.Term.iloc[0]  # Get the first term from res2d
        # Check if the term exists in results
        if term in gsea_res.results:
            print(term)
            fig_path = plot_enhanced_gsea(gsea_res, term, i,outraw)
            print(f"Plot saved to: {fig_path}")
    
            #other plotting options- gseaplot - does not work
            # Call the plotting function directly gseaplot... (dont forget import gseaplot for that)
            #gseaplot(rank_metric=gsea_res.ranking,
             #term=term, 
             #hits=enrichment_results.get('hits'),
             #nes=enrichment_results.get('nes'),
             #pval=enrichment_results.get('pval'),
             #fdr=enrichment_results.get('fdr'),
             #RES=enrichment_results.get('RES'),
             #figsize=(8, 6))

            #simple gsea plot - does not work
            #gsea_plot= gsea_res.plot(term)
            #print(f"Plot object: {gsea_plot}")  #
            #fig=plt.gcf()
            #plt.tight_layout()
            #plt.draw()
            #plt.pause(1)
            #fig=plt.gcf()
            #print(plt.get_fignums())
            #print(f"Figure dimensions: {fig.get_size_inches()}") 
            #fig_path = f"{outraw}/gseaplot_{term}.png"
            #fig.savefig(fig_path, dpi=300, bbox_inches='tight')
            #plt.close()
        else:
            print(f"Term {term} not found in the results.")
        #use the other results to build the final dataframe..
        condition_cell.append(i)
        gsea_result.append(gsea_res.res2d.loc[0,"FDR q-val"])
        num_lead.append(gsea_res.res2d.loc[0,"Tag %"])
        lead_genes.append(gsea_res.res2d.loc[0,"Lead_genes"])
        per_lead.append(gsea_res.res2d.loc[0,"Gene %"])
        NES.append(gsea_res.res2d.loc[0,"NES"])
        print(gsea_res.res2d.loc[0,"FDR q-val"])
        gsea_res.res2d.to_csv(outgsea)
        if gsea_res.res2d.loc[0,"FDR q-val"]<0.05:
            print(f"gene list result for: {i} is significant!")
        else:
            print(f"gene list result for: {i} is non significant")
        #except: #was written before with try statement
            #print(f"could not find enrichment in this gene set for condition {i}")
    results=pd.DataFrame(data={"condition":condition_cell,"lead_genes":lead_genes,"num_genes_in_lead":num_lead,"%genes_in_lead":per_lead,"nes_score":NES,"gsea_pval":gsea_result})
    reject, pvals_corr = smm.multipletests(gsea_result, method='fdr_bh')[:2]
    results["p-adj"]= pvals_corr
    final_outpath=out_folder+gmt_name+"/"+"full_results.csv"
    print("final outpath: ", final_outpath)
    results.to_csv(final_outpath)

def check_ges_path(condition,data_type):
    if data_type=="cortex":
        if os.path.exists(f"/miridan-data/annaludmir/data/ges_results/ges_spec_v3_{condition}.csv"):
            path=f"/miridan-data/annaludmir/data/ges_results/ges_spec_v3_{condition}.csv"
        else:
            print(f"condition {condition} is not a valid condition name")
            sys.exit()
    elif data_type=="data_all":
        if os.path.exists(f"/miridan-data/annaludmir/data/ges_results/ges_spec_big_{condition}.csv"):
            path=f"/miridan-data/annaludmir/data/ges_results/ges_spec_big_{condition}.csv"
        else:
            print(f"condition {condition} is not a valid condition name")
            sys.exit()
    else:
        print("data_type is not correct give one of this values: data_all or cortex")
        sys.exit()
    return path


if __name__=="__main__":
    #add ges_results so the results will go to the specifies folder
    out_folder=out_folder+"ges_results_above1"+"/"+data_type+"/"
    gene_list_name=gene_list.split("/")[-1][:-4]
    gmt_file=gmt_folder+f"{gene_list_name}.gmt"
    print(f"starting ges enrichment analysis on {gmt_file}")
    run_gsea(gmt_file,conditions,out_folder,data_type)
    print('DONE')
