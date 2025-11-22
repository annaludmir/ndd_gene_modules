import marimo

__generated_with = "0.17.8"
app = marimo.App(width="full")


@app.cell
def _():
    import sys
    import subprocess

    gene_list="/miridan-data/annaludmir/data/genes/sysndd_id_abnormal_facial_shape_abnormal_heart_morphology_for_pipeline.csv"
    out_folder="/miridan-data/annaludmir/data/enrichment_results/sysndd_id_abnormal_facial_shape_abnormal_heart_morphology_for_pipeline"
    choice = "ges_enrichment"
    data_type="cortex"
    extra_args = ["/miridan-data/annaludmir/data/genes/","['radialglia','proliferating','differentiating','IPC','NPCs','G1','S','G2M','PostM','Neuron','Glioblast','Neuroblast','differentiating_non_cycling','proliferating_non_cycling']"] 
              
    if choice == "ges_enrichment":#ges_enrichment
        result = subprocess.run(["python","/miridan-data/annaludmir/notebooks/get_gmt.py",gene_list]) #creates gmt file in the same folder if it does not exist
        result = subprocess.run(["python", "/miridan-data/annaludmir/notebooks/search_enrichment_new.py",gene_list ,extra_args[0],extra_args[1] , out_folder, data_type])  # Runs ges_script
        print(result.stderr)
        result = subprocess.run(["python", "/miridan-data/annaludmir/notebooks/create_figs_ges.py",gene_list,out_folder,data_type]) 
        print(result.stderr)
    elif choice == "dge_enrichment": #run_deseq
        subprocess.run(["python", "/miridan-data/annaludmir/notebooks/deseq_calculations.py",gene_list, extra_args[0],extra_args[1], out_folder,data_type])  # Runs deseq fisher script
        subprocess.run(["python", "create_figs_deseq.py",gene_list,out_folder,data_type])  
    else:
        print(f"Invalid option: {choice}")
        sys.exit(1)
    return


if __name__ == "__main__":
    app.run()
