#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=python
#SBATCH --mem=500G
# #SBATCH --mem=250G
#SBATCH --account=miridan-users_v2  
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err
#SBATCH --time=1-00:00:00
# # public
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public
# # miris partition
# #SBATCH --partition=gpu-miridan-pool 
# #SBATCH --qos=owner
# #SBATCH --gres=gpu:0
# # deprecated
# # SBATCH --account=public-users_v2

set -euo pipefail

#module load mamba/mamba1.4.2-environmentally
module load mamba/mamba-1.5.8
#mamba activate /scratch200/reutj/conda-envs/jupyter-scanpy
mamba activate /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new
cd /miridan-data/annaludmir/ndd_gene_modules/

mamba run -p /miridan-data/annaludmir/conda-envs/jupyter-scanpy_new python modules/early_late_go_heatmap.py \
  --leading-genes FOXG1 CDK6 HPDL BRIP1 CEP152 CIT BLM RAD51 KIF11 KNL1 FANCI NCAPH CDT1 ORC1 POC1A CDC6 BUB1B DNA2 CKAP2L NCAPD2 ZEB2 BUB1 STIL CENPF PLK4 ASPM FANCB BRCA2 WDR62 CENPE LMNB2 CEP55 FANCD2 HHAT FANCA CCDC88A GINS2 MRE11 NDE1 CCND2 DDX11 TRAIP RBBP8 FANCL TRIP13 GINS3 CEP135 NUF2 FANCE GMNN TUBGCP3 CDK5RAP2 ORC6 KIF14 VRK1 MCM7 LMNB1 SMC1A DIAPH1 PCNT NUP188 NUP107 FANCG GPT2 NCAPD3 TCF4 NUP214 SASS6 PUS7 SMO RAD50 TTI1 SMC5 ZNF526 RAD21 ATP1A2 TOP3A NSD2 FBRSL1 SLF2 ERCC5 OSGEP SMC3 EFTUD2 ANKLE2 \
  --go-term-file GO_terms/microcephaly_GO_enrichment/Region_-_Forebrain/enrichr_results.csv \
  --h5ad-path data/human_dev_without_week_5.h5ad \
  --subfolder-name microcephaly_forebrain_leading_genes
echo "Python exit code: $rc"
exit $rc
