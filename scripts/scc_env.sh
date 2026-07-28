# scc_env.sh -- SOURCE this to get a working SOD1-FEP shell on the BU SCC.
#
#   source scripts/scc_env.sh
#
# Loads every module the pipeline needs (in the order Lmod requires) and activates the
# conda env, reading all names from config/pipeline.yaml so the interactive shell, the
# setup script and the SGE array tasks cannot drift apart -- a mismatch between them is
# exactly how a 864-task array dies 864 times.
#
# Override any of these from the environment before sourcing:
#   CONDA_MODULE  GROMACS_MODULE  GROMACS_PREREQ  CUDA_MODULE  SOD1_ENV  GMX
#
# Sourced, not executed: it must mutate YOUR shell. It deliberately does not `set -e`,
# so a failed module load warns instead of killing an interactive session.

_SCC_CONFIG="${SCC_CONFIG:-config/pipeline.yaml}"
if [[ ! -f "${_SCC_CONFIG}" ]]; then
  echo "scc_env.sh: ${_SCC_CONFIG} not found -- source this from the repo root." >&2
  return 1 2>/dev/null || exit 1
fi

# Read a cluster.<key> without needing python/yaml (this runs before conda activate).
_scc_cfg() {
  grep -E "^[[:space:]]*$1:" "${_SCC_CONFIG}" | head -1 \
    | sed -E 's/^[^:]+:[[:space:]]*//; s/[[:space:]]*#.*$//; s/^"//; s/"$//'
}

CONDA_MODULE="${CONDA_MODULE:-$(_scc_cfg conda_module)}"
GROMACS_MODULE="${GROMACS_MODULE:-$(_scc_cfg gromacs_module)}"
GROMACS_PREREQ="${GROMACS_PREREQ:-$(_scc_cfg gromacs_prereq_modules)}"
CUDA_MODULE="${CUDA_MODULE:-$(_scc_cfg cuda_module)}"
GMX="${GMX:-$(_scc_cfg gmx_binary)}"
SOD1_ENV="${SOD1_ENV:-sod1-fep}"

# Lmod defines `module` as a shell FUNCTION via a profile script. A non-interactive
# shell -- notably `qrsh ... bash -lc '...'` and some batch contexts -- can start without
# it, and then every `module load` silently no-ops and GROMACS later dies with
# "error while loading shared libraries: libmkl_core.so.2". Bootstrap it explicitly.
if ! type module >/dev/null 2>&1; then
  # NOTE the globs: the SCC's Lmod lives in a VERSIONED directory
  # (LMOD_PKG=/usr/local/lmod/8.7.12), and LMOD_PKG itself is not exported into the
  # qrsh shell -- so both a bare .../lmod/init/bash and "${LMOD_PKG}/init/bash" miss it.
  # cluster.lmod_init pins the path if this discovery ever fails.
  for _init in "$(_scc_cfg lmod_init)" \
               "${LMOD_PKG:-/nonexistent}/init/bash" \
               /etc/profile.d/lmod.sh /etc/profile.d/modules.sh \
               /usr/local/lmod/*/init/bash /usr/share/lmod/*/init/bash \
               /share/module.8/lmod/*/init/bash; do
    if [[ -n "${_init}" && -f "${_init}" ]]; then
      # shellcheck disable=SC1090
      if source "${_init}" 2>/dev/null && type module >/dev/null 2>&1; then
        echo "scc_env.sh: bootstrapped Lmod from ${_init}" >&2
        break
      fi
    fi
  done
fi

_have_module=1
if ! type module >/dev/null 2>&1; then
  _have_module=0
  echo "ERROR: Lmod's 'module' command is unavailable and could not be bootstrapped." >&2
  echo "       Nothing will load and GROMACS will fail on missing MKL/CUDA libraries." >&2
  echo "       Set cluster.lmod_init in config/pipeline.yaml to your Lmod init script" >&2
  echo "       (look for \$LMOD_PKG/init/bash in a login shell), or run interactively." >&2
fi

# if-statements rather than `[[ ... ]] && cmd`: this file is sourced by scripts running
# under `set -e`, where a false test at the head of an && list can abort the caller.
if [[ "${_have_module}" == "1" ]]; then
  if [[ -n "${CONDA_MODULE}" ]]; then
    module load "${CONDA_MODULE}" || echo "WARN: module load ${CONDA_MODULE} failed" >&2
  fi
  # Runtime prerequisites for the CUDA GROMACS build (MKL via intel, CUDA, gcc), in order.
  for _m in ${GROMACS_PREREQ}; do
    module load "${_m}" || echo "WARN: module load ${_m} failed" >&2
  done
  if [[ -n "${GROMACS_MODULE}" ]]; then
    module load "${GROMACS_MODULE}" || echo "WARN: module load ${GROMACS_MODULE} failed" >&2
  fi
  if [[ -n "${CUDA_MODULE}" ]]; then
    module load "${CUDA_MODULE}" || echo "WARN: module load ${CUDA_MODULE} failed" >&2
  fi
fi

if [[ -n "${CONDA_MODULE}" ]]; then
  if source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null; then
    conda activate "${SOD1_ENV}" || echo "WARN: could not activate ${SOD1_ENV}" >&2
  else
    echo "WARN: conda profile script not found" >&2
  fi
fi

# Repo root on PYTHONPATH so `python -m src.fep.window` works from anywhere.
export PYTHONPATH="${PWD}:${PYTHONPATH:-}"

# GMXLIB: pdb2gmx only finds pmx's mutation force fields if this points at the directory
# holding them. Derived from the installed pmx so it follows the env, not a hardcoded
# path (the engine sets this too; exporting it here makes interactive gmx calls work).
if python -c "import pmx" >/dev/null 2>&1; then
  GMXLIB="$(python -c "import pmx, os; print(os.path.join(os.path.dirname(pmx.__file__), 'data', 'mutff45'))" 2>/dev/null)"
  if [[ -d "${GMXLIB}" ]]; then export GMXLIB; else unset GMXLIB; fi
fi

# A specific GROMACS build (e.g. the CUDA tree, which the module does NOT put on PATH)
# is selected by sourcing its GMXRC -- that sets PATH/LD_LIBRARY_PATH for that install.
GMX_GMXRC="${GMX_GMXRC:-$(_scc_cfg gmx_gmxrc)}"
if [[ -n "${GMX_GMXRC}" ]]; then
  if [[ -f "${GMX_GMXRC}" ]]; then
    source "${GMX_GMXRC}" || echo "WARN: sourcing ${GMX_GMXRC} failed" >&2
  else
    echo "WARN: cluster.gmx_gmxrc=${GMX_GMXRC} not found" >&2
  fi
fi

if [[ -z "${GMX}" ]]; then
  for _c in gmx gmx_mpi; do
    if command -v "${_c}" >/dev/null 2>&1; then GMX="${_c}"; break; fi
  done
fi
export GMX
if ! command -v "${GMX:-gmx}" >/dev/null 2>&1; then
  echo "WARN: '${GMX:-gmx}' not on PATH -- check cluster.gromacs_module" >&2
fi

unset _m _c
