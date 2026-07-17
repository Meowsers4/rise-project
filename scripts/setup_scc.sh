#!/bin/bash -l
#
# setup_scc.sh -- one-time environment build + verification for the BU SCC.
# Run on an SCC login node (env build) then on a GPU node (the GPU check).
#
#   bash scripts/setup_scc.sh            # build env + import checks (login node)
#   qrsh -l gpus=1 -l gpu_c=7.0          # grab a GPU, then:
#   bash scripts/setup_scc.sh --gpu      # OpenMM GPU/platform check
#
# Override module/conda locations for your SCC if the defaults are wrong:
#   CONDA_MODULE=miniconda CUDA_MODULE=cuda CONDA_BASE=$HOME/miniconda3 SOD1_ENV=sod1-fep
set -euo pipefail

CONDA_MODULE="${CONDA_MODULE:-miniconda}"
CUDA_MODULE="${CUDA_MODULE:-cuda}"
SOD1_ENV="${SOD1_ENV:-sod1-fep}"

echo "== module load =="
module load "${CONDA_MODULE}" || echo "WARN: 'module load ${CONDA_MODULE}' failed -- set CONDA_MODULE"
module load "${CUDA_MODULE}"  || echo "WARN: 'module load ${CUDA_MODULE}' failed -- set CUDA_MODULE"

source "$(conda info --base)/etc/profile.d/conda.sh"

if [[ "${1:-}" != "--gpu" ]]; then
  echo "== create env from env/environment.yml (skips if it exists) =="
  conda env list | grep -qE "\b${SOD1_ENV}\b" || conda env create -f env/environment.yml -n "${SOD1_ENV}"
  conda activate "${SOD1_ENV}"

  echo "== import checks =="
  python - <<'PY'
import importlib
for m in ["openmm", "openmmtools", "perses", "pymbar", "mdtraj", "pdbfixer", "yaml", "numpy"]:
    try:
        mod = importlib.import_module(m)
        print(f"  OK   {m:12} {getattr(mod, '__version__', '')}")
    except Exception as e:
        print(f"  FAIL {m:12} {e}")
PY
  echo "== pipeline unit tests (no GPU) =="
  python -m pytest tests/ -q
  echo
  echo "Next: grab a GPU (qrsh -l gpus=1 -l gpu_c=7.0) and run: bash scripts/setup_scc.sh --gpu"
else
  conda activate "${SOD1_ENV}"
  echo "== OpenMM platform / GPU check =="
  python - <<'PY'
import openmm
from openmm import Platform
print("OpenMM", openmm.__version__)
for i in range(Platform.getNumPlatforms()):
    print("  platform:", Platform.getPlatform(i).getName())
try:
    print("  CUDA platform available:", Platform.getPlatformByName("CUDA") is not None)
except Exception as e:
    print("  CUDA platform NOT available:", e)
PY
  nvidia-smi --query-gpu=name,compute_cap,memory.total --format=csv || echo "nvidia-smi unavailable"
fi
