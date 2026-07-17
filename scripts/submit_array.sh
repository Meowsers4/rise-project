#!/bin/bash -l
#
# submit_array.sh -- SGE array submission helper for Stage 3 (alchemical FEP) on the
# BU SCC (Open Grid Scheduler / SGE). One array of N tasks PER VARIANT, where
# N = legs x lambda_windows x replicates (2 x 18 x 5 = 180 with the current config).
# SGE_TASK_ID (1..N) decodes to (leg, window, replicate); the variant is passed in.
#
# Usage:
#   mkdir -p logs/fep                       # -o dir must exist at submit time
#   qsub -v VARIANT=A4V scripts/submit_array.sh
#
# One window == one small MD sim on ONE GPU. Parallelism is at the JOB level: the
# scheduler runs as many array tasks concurrently as fair-share allows. GPUs are
# requested by compute capability (gpu_c), NOT by model. The scheduler sets
# CUDA_VISIBLE_DEVICES -- src.fep.window must use it as-is and never hardcode it.
#
# NOTE: the -t range below must equal legs*windows*replicates. It is validated at
# runtime against config/pipeline.yaml; if you change the config, resubmit with
# `qsub -t 1-<N> ...` or update the directive.
#
#$ -P rise-batteries
#$ -l gpus=1
#$ -l gpu_c=7.0
#$ -l h_rt=12:00:00
#$ -pe omp 8
#$ -N sod1_fep
#$ -j y
#$ -o logs/fep/
#$ -t 1-180

set -euo pipefail

: "${VARIANT:?set VARIANT, e.g. qsub -v VARIANT=A4V scripts/submit_array.sh}"

CONFIG="config/pipeline.yaml"

# ---- read fan-out geometry from config (no hardcoded params; CLAUDE.md rule 4) --
read_fep() { python -c "import yaml,sys; print(yaml.safe_load(open('${CONFIG}'))['fep'][sys.argv[1]])" "$1"; }
N_WINDOWS=$(read_fep lambda_windows)
N_REPS=$(read_fep replicates)
mapfile -t LEGS < <(python -c "import yaml; print('\n'.join(yaml.safe_load(open('${CONFIG}'))['fep']['legs']))")
N_TASKS=$(( ${#LEGS[@]} * N_WINDOWS * N_REPS ))

# Guard: the SGE array size must match the config-derived fan-out.
if [[ "${SGE_TASK_LAST:-$N_TASKS}" -ne "$N_TASKS" ]]; then
  echo "ERROR: array size ${SGE_TASK_LAST:-?} != legs*windows*reps=${N_TASKS}." >&2
  echo "       Resubmit with: qsub -t 1-${N_TASKS} -v VARIANT=${VARIANT} $0" >&2
  exit 2
fi

# ---- decode SGE_TASK_ID (1..N) -> (leg, window, replicate) ----------------------
idx=$(( SGE_TASK_ID - 1 ))                     # 0-based
rep=$(( idx % N_REPS ))
win=$(( (idx / N_REPS) % N_WINDOWS ))
leg_i=$(( idx / (N_REPS * N_WINDOWS) ))
leg="${LEGS[$leg_i]}"

echo "task ${SGE_TASK_ID}: variant=${VARIANT} leg=${leg} window=${win} rep=${rep}"

# ---- activate the conda env (override via env vars for your SCC; see env/modules.md) --
# Pass FEP_MOCK=1 (e.g. `qsub -v VARIANT=A4V,FEP_MOCK=1 ...`) for a scheduler shake-out
# that runs the synthetic window instead of the (GPU-only) Perses window.
[[ -n "${CUDA_MODULE:-}" ]] && module load "${CUDA_MODULE}"
source "${CONDA_BASE:-$HOME/miniconda3}/etc/profile.d/conda.sh"
conda activate "${SOD1_ENV:-sod1-fep}"

echo "host=$(hostname) CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>} mock=${FEP_MOCK:-0}"

python -m src.fep.window \
    --variant "${VARIANT}" \
    --leg "${leg}" \
    --window "${win}" \
    --rep "${rep}" \
    --config "${CONFIG}" \
    --out "results/fep/${VARIANT}/${leg}/w${win}_r${rep}.npz"
