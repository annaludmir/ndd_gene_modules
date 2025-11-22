import sys
import subprocess

if len(sys.argv) < 6:
    print("Usage: python run_script.py gene_list.csv out_folder [ges_enrichment | dge_enrichment][all_data | cortex] [additional argumants]", 
    "if ges_enrichment is written add arguments: 1- folder with gmt_paths 2-conditions to check" ,
    " if dge_enrichment is written additional args are: 1- pseudobulks data file folder, 2- gene names of the genes indata that express in more than 100 cells" )
    sys.exit(1)

gene_list=sys.argv[1]
out_folder=sys.argv[2]
choice = sys.argv[3]
data_type=sys.argv[4]
extra_args = sys.argv[5:]

if choice == "ges_enrichment":#ges_enrichment
    subprocess.run(["python","/miridan-data/annaludmir/notebooks/get_gmt.py",gene_list])#creates gmt file in the same folder if it does not exist
    subprocess.run(["python", "/miridan-data/annaludmir/notebooks/search_enrichment_new.py",gene_list ,extra_args[0],extra_args[1] , out_folder, data_type])  # Runs ges_script
    subprocess.run(["python", "create_figs_ges.py",gene_list,out_folder,data_type]) 
elif choice == "dge_enrichment": #run_deseq
    subprocess.run(["python", "/miridan-data/annaludmir/notebooks/deseq_calculations.py",gene_list, extra_args[0],extra_args[1], out_folder,data_type])  # Runs deseq fisher script
    subprocess.run(["python", "create_figs_deseq.py",gene_list,out_folder,data_type])  
else:
    print(f"Invalid option: {choice}")
    sys.exit(1)
