import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

def plot_bar_chart(deseq_results_df, output_file, gene_list_name,comparison_criteria):
    x_pos = np.arange(len(deseq_results_df)) * 0.9

    fig, ax = plt.subplots(figsize=(4,4))
    
    ax.bar(x_pos,-np.log10(deseq_results_df["p-adj"]),color='#1f77b4',edgecolor='black',linewidth=0.7,width=0.7, alpha=0.8)
    ax.set_xticks(x_pos)
    xlabels=["_".join(i.split("_")[:-2]) for i in deseq_results_df.comparison]
    ax.set_xticklabels(xlabels, fontsize=10,rotation=80)
    ax.set_ylim(0, max(-np.log10(deseq_results_df["p-adj"])) + 0.2)
    ax.axhline(y=-np.log10(0.05), color="#555555", linestyle="--",linewidth=1.5)
    ax.set_title(f"fisher enrichment for {gene_list_name} - {comparison_criteria} deseq results ",pad=14,fontsize=11)
    ax.set_ylabel('-log10(padj)',fontsize=10)
    plt.tight_layout()

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

if __name__=="__main__":
    gene_list_name=hsg_list.split("/")[-1][:-4]
    print("starting plotting deseq results for hsg list: ", gene_list_name)
    #upload ges_results for the gene list
    deseq_results_path=output_path+"deseq_fisher_results/"+data_type+"/"+gene_list_name+"/"
    os.makedirs(deseq_results_path,exist_ok=True)#make sure path exists
    #create figure output directory
    figs_outpath=output_path+"figs/"+gene_list_name+"/"+data_type+"/deseq_results/"
    os.makedirs(figs_outpath,exist_ok=True)#make sure path exists..
    for directory in os.listdir(deseq_results_path):
        joined_path=os.path.join(deseq_results_path,directory,"full_results.csv")
        deseq_results_df=pd.read_csv(joined_path)
        print("uploaded deseq results from: ",  joined_path)
        if "cellclass" in directory:  
            #make figure output path
            fig_full_path=figs_outpath+"cell_class_enrichment.png"
            if data_type=="data_all":
                deseq_results_df=deseq_results_df.iloc[[1,6,7,9,10,0,8,5,4,2,3],:]#reorder rows for better graph representation
            else:#cortex
                deseq_results_df=deseq_results_df.iloc[[6,5,7,3,0,1,2,4],:]#reorder rows for better representation in the graph
            plot_bar_chart(deseq_results_df, fig_full_path, gene_list_name,"cell class")
            print("plotting cell class deseq results, output: ", fig_full_path)
        elif "cellcycle" in directory:
            deseq_results_df=deseq_results_df.iloc[[2,3,4,5,1,0],:]#reorder rows for better representation in the graph
            fig_full_path=figs_outpath+"cell_cycle_enrichment.png"
            plot_bar_chart(deseq_results_df, fig_full_path, gene_list_name,"cell cycle")
            print("plotting cell cycle deseq results, output: ", fig_full_path)
        elif "region" in directory:
            fig_full_path=figs_outpath+"brain_region_enrichment.png"
            deseq_results_df=deseq_results_df.iloc[[8,2,9,3,0,1,7,4,5,6],:]#reorder rows for better representation in the graph
            plot_bar_chart(deseq_results_df, fig_full_path, gene_list_name,"region")
            print("plotting region deseq results, output: ", fig_full_path)


    print('DONE')

    
