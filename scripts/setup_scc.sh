#!/bin/bash -l
#
# setup_scc.sh -- one-time environment build + verification for the BU SCC.
# Run on an SCC login node (env build) then on a GPU node (the GPU check).
#
#   bash scripts/setup_scc.sh            # build env + import checks (login node)
#   qrsh -l gpus=1 -l gpu_c=7.0          # grab a GPU, then:
#   bash scripts/setup_scc.sh --gpu      # GROMACS GPU check
#
# Module names default from config (cluster.conda_module / cluster.gromacs_module).
# Override for a different SCC if needed:
#   CONDA_MODULE=miniconda/25.3.1 GROMACS_MODULE=gromacs/2025.3 SOD1_ENV=sod1-fep
set -euo pipefail

CONFIG="config/pipeline.yaml"
# Read a cluster.<key> from config without needing python/yaml (pre-env bootstrap).
cfg_get() { grep -E "^[[:space:]]*$1:" "${CONFIG}" | head -1 | sed -E 's/^[^:]+:[[:space:]]*//; s/[[:space:]]*#.*$//; s/^"//; s/"$//'; }

# Precedence: env-var override > config value.
CONDA_MODULE="${CONDA_MODULE:-$(cfg_get conda_module)}"
GROMACS_MODULE="${GROMACS_MODULE:-$(cfg_get gromacs_module)}"
GROMACS_PREREQ="${GROMACS_PREREQ:-$(cfg_get gromacs_prereq_modules)}"
CUDA_MODULE="${CUDA_MODULE:-$(cfg_get cuda_module)}"   # optional; gromacs usually pulls CUDA
GMX="${GMX:-$(cfg_get gmx_binary)}"
SOD1_ENV="${SOD1_ENV:-sod1-fep}"

: "${CONDA_MODULE:?set cluster.conda_module in ${CONFIG} or CONDA_MODULE}"
: "${GROMACS_MODULE:?set cluster.gromacs_module in ${CONFIG} or GROMACS_MODULE}"

echo "== module load =="
module load "${CONDA_MODULE}"   || echo "WARN: 'module load ${CONDA_MODULE}' failed -- set CONDA_MODULE"
# Lmod prerequisites first, IN ORDER -- the SCC's gromacs/2025.3 is an OpenMPI build and
# refuses to load without openmpi already present (cluster.gromacs_prereq_modules).
for _m in ${GROMACS_PREREQ}; do
  module load "${_m}" || echo "WARN: prerequisite 'module load ${_m}' failed"
done
module load "${GROMACS_MODULE}" || echo "WARN: 'module load ${GROMACS_MODULE}' failed -- set GROMACS_MODULE"
[[ -n "${CUDA_MODULE}" ]] && { module load "${CUDA_MODULE}" || echo "WARN: cuda module ${CUDA_MODULE} failed"; }

# MPI builds ship gmx_mpi, non-MPI builds ship gmx. Auto-detect unless pinned in config.
if [[ -z "${GMX}" ]]; then
  for _cand in gmx gmx_mpi; do
    command -v "${_cand}" >/dev/null 2>&1 && { GMX="${_cand}"; break; }
  done
fi
if [[ -n "${GMX}" ]]; then
  echo "  gmx binary: ${GMX} -> $(command -v "${GMX}")"
  echo "  >> record this as cluster.gmx_binary in ${CONFIG} if auto-detection is wrong"
else
  echo "  WARN: neither 'gmx' nor 'gmx_mpi' on PATH after loading ${GROMACS_MODULE}"
fi

source "$(conda info --base)/etc/profile.d/conda.sh"

if [[ "${1:-}" != "--gpu" ]]; then
  echo "== create env from env/environment.yml (skips if it exists) =="
  conda env list | grep -qE "\b${SOD1_ENV}\b" || conda env create -f env/environment.yml -n "${SOD1_ENV}"
  conda activate "${SOD1_ENV}"

  echo "== import checks =="
  python - <<'PY'
import importlib
for m in ["openmm", "pdbfixer", "pmx", "pymbar", "alchemlyb", "pandas", "yaml", "numpy"]:
    try:
        mod = importlib.import_module(m)
        print(f"  OK   {m:12} {getattr(mod, '__version__', '')}")
    except Exception as e:
        print(f"  FAIL {m:12} {e}")
PY
  echo "== pmx mutation force fields available =="
  python - <<'PY'
import os
try:
    import pmx
    d = os.path.join(os.path.dirname(pmx.__file__), "data", "mutff")
    print("  ", sorted(os.listdir(d)) if os.path.isdir(d) else f"no mutff dir at {d}")
except Exception as e:
    print("   pmx unavailable:", e)
PY
  echo "== pipeline unit tests (no GPU) =="
  python -m pytest tests/ -q
  echo
  echo "Next: grab a GPU (qrsh -l gpus=1 -l gpu_c=7.0) and run: bash scripts/setup_scc.sh --gpu"
else
  conda activate "${SOD1_ENV}"
  echo "== GROMACS / GPU check =="
  if [[ -n "${GMX}" ]]; then
    "${GMX}" --version 2>/dev/null \
      | grep -Ei "GROMACS version|Precision|GPU support|SIMD|CUDA (driver|runtime)" \
      || echo "  WARN: '${GMX} --version' produced nothing"
  else
    echo "  WARN: no gmx binary found -- cannot check GROMACS/GPU"
  fi
  nvidia-smi --query-gpu=name,compute_cap,memory.total --format=csv || echo "nvidia-smi unavailable"
fi
