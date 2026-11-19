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
# Rerunning is safe: an existing env is reused, torch is reinstalled for
# TORCH_CUDA, and pip skips the rest of what is installed.

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

# PyTorch from PyPI is built for the newest CUDA (13.x), which needs NVIDIA
# driver >= 580. The miridan GPU nodes have driver 570 (CUDA 12.8), so torch
# comes from PyTorch's CUDA 12.8 index instead. Change TORCH_CUDA if the
# driver changes ('CUDA Version' in nvidia-smi above): cu126, cu128, cu130 ...
TORCH_CUDA="${TORCH_CUDA:-cu128}"

echo "== 2. PyTorch for ${TORCH_CUDA}"
mamba run -p "$ENV" python -m pip install --upgrade pip
# a torch built for another CUDA would count as installed: remove it first
mamba run -p "$ENV" python -m pip uninstall -y torch triton || true
mamba run -p "$ENV" python -m pip install torch --index-url "https://download.pytorch.org/whl/${TORCH_CUDA}"

echo "== 3. CycleVI (keeps the torch above; pulls scvi-tools, scanpy)"
mamba run -p "$ENV" python -m pip install cyclevi

echo "== 4. check the install and the GPU"
# (written to a file: mamba run does not pass a script on stdin to python)
CHECK=$(mktemp --suffix=.py)
cat > "$CHECK" <<'PY'
import torch, scvi
from importlib.metadata import version
print("cyclevi", version("cyclevi"), "| scvi-tools", scvi.__version__,
      "| torch", torch.__version__, "| built for CUDA", torch.version.cuda)
ok = torch.cuda.is_available()
print("GPU visible to torch:", ok, "|", torch.cuda.get_device_name(0) if ok else "-")
if ok:
    x = torch.randn(1000, 1000, device="cuda")
    print("GPU computation ok:", (x @ x).shape)
else:
    print("[problem] torch does not see the GPU: compare 'CUDA Version' in nvidia-smi "
          "with the CUDA version torch was built for.")
PY
mamba run -p "$ENV" python "$CHECK"
rm -f "$CHECK"
echo "Done."
