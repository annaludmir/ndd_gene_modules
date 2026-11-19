#!/bin/bash
#SBATCH --mail-user=annaludmir@mail.tau.ac.il
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=cyclevi_setup
#SBATCH --account=miridan-users_v2
#SBATCH --partition=gpu-miridan-pool
#SBATCH --qos=owner
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0-02:00:00
#SBATCH --output=/miridan-data/annaludmir/jobs_output/%j.out
#SBATCH --error=/miridan-data/annaludmir/jobs_output/%j.err

# One-time setup of the CycleVI environment, as a job (mamba is not
# available on the login node). Runs on a GPU node so the last step can
# check that PyTorch sees the GPU.
#   sbatch running_scripts/setup_cyclevi_env.sh
# Rerunning is safe: an existing env is reused and pip skips what is installed.

set -euo pipefail

ENV=/miridan-data/annaludmir/conda-envs/cyclevi
export PYTHONNOUSERSITE=1   # ignore packages in ~/.local

module load mamba/mamba-1.5.8

echo "== GPU on this node"
nvidia-smi || echo "[warning] nvidia-smi failed"

echo "== 1. create the env"
if [[ -d "$ENV" ]]; then
  echo "  $ENV exists, reusing it"
else
  mamba create -y -p "$ENV" python=3.11
fi

echo "== 2. install CycleVI (pulls PyTorch, scvi-tools, scanpy)"
mamba run -p "$ENV" python -m pip install --upgrade pip
mamba run -p "$ENV" python -m pip install cyclevi

echo "== 3. check the install and the GPU"
mamba run -p "$ENV" python - <<'PY'
import torch, scvi, scanpy, cyclevi
from importlib.metadata import version
print("cyclevi", version("cyclevi"), "| scvi-tools", scvi.__version__,
      "| torch", torch.__version__, "| torch built for CUDA", torch.version.cuda)
ok = torch.cuda.is_available()
print("GPU visible to torch:", ok, "|", torch.cuda.get_device_name(0) if ok else "-")
if not ok:
    print("[problem] torch does not see the GPU: compare 'CUDA Version' in nvidia-smi "
          "above with the CUDA version torch was built for.")
PY
echo "Done."
