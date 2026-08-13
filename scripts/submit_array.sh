#!/bin/bash -l
#
# submit_array.sh -- SGE array submission helper for Stage 3 (alchemical FEP) on the
# BU SCC (Open Grid Scheduler / SGE). One array of N tasks PER VARIANT, where
# N = fep.legs x fep.lambda_windows x fep.replicates, read from config/pipeline.yaml
# below -- never restate the product in prose (it drifts). The #$ -t default must equal
# N; the runtime guard recomputes N from config and refuses to launch if the submitted
# array size disagrees. SGE_TASK_ID (1..N) decodes to (leg, window, replicate); the
# variant is passed in.
#
# Usage (submit from the REPO ROOT -- the job runs with -cwd, see below):
#   cd <repo> && mkdir -p logs/fep          # -o dir must exist at submit time
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
# Pin the GPU MODEL, not just a capability floor. gpu_c has SGE relation `<=`, i.e. it is
# a MINIMUM -- so a card NEWER than this GROMACS build also satisfies it. scc-702 carries
# an RTX PRO 6000 Blackwell (compute 12.0); GROMACS 2025.3 was built against CUDA 12.8 and
# ships no sm_120 kernels, so the driver JIT-compiles the embedded PTX and fails with
# `CUDA error #218 (cudaErrorInvalidPtx)` -- after mdrun has already started, so it looks
# like a run failure rather than a placement problem. Killed A4V unfolded/w17_r0 on
# 2026-08-11 while the other 107 windows ran fine on L40S. 32 L40S nodes are available, so
# pinning costs no meaningful queue time and makes window timings comparable.
#$ -l gpu_type=L40S
#$ -l h_rt=12:00:00
#$ -pe omp 8
#$ -N sod1_fep
#$ -j y
# Rerunnable: src.fep.window exits 99 when mdrun cannot get a usable GPU on this host,
# and Grid Engine only honours 99-means-reschedule for a rerunnable task. Safe here
# because a rescheduled task resumes from prod.cpt and assert_resumable() refuses to
# resume across a protocol change.
#$ -r y
#$ -cwd
#$ -o logs/fep/
# MUST equal fep.legs * fep.lambda_windows * fep.replicates (2*18*3). SGE cannot template
# this from config, so the guard below re-derives it and aborts on a mismatch. If you raise
# fep.replicates to 5, change this to 1-180 IN THE SAME COMMIT or every task exits 2.
#$ -t 1-108

# NOT -u: GROMACS's GMXRC (sourced via scc_env.sh) reads unbound $shell/$GMXLDLIB and
# would abort every array task before python runs. Batch must match interactive.
set -eo pipefail

: "${VARIANT:?set VARIANT, e.g. qsub -v VARIANT=A4V scripts/submit_array.sh}"
# -u is off (GMXRC reads unbound vars), so an unset SGE_TASK_ID would arithmetic to
# idx=-1 -> rep=-1 and cheerfully write w0_r-1.npz. This script is only ever a task body.
: "${SGE_TASK_ID:?not running as an SGE array task -- submit with qsub, do not run directly}"

# -cwd above runs the task in the SUBMIT directory. Without it SGE starts in $HOME and
# every relative path here (CONFIG, results/, -o logs/fep/) silently points at the wrong
# place -- the task then dies before it ever reaches python. Submit from the repo root.
CONFIG="config/pipeline.yaml"
if [[ ! -f "${CONFIG}" ]]; then
  echo "ERROR: ${CONFIG} not found in $(pwd). Submit from the repo root:" >&2
  echo "       cd <repo> && mkdir -p logs/fep && qsub -v VARIANT=${VARIANT} $0" >&2
  exit 2
fi

# ---- modules + conda env, from the ONE shared helper -----------------------------
# scripts/scc_env.sh loads cluster.gromacs_prereq_modules -> cluster.gromacs_module ->
# cluster.conda_module and activates the env, all read from config. Sourcing it here (and
# from an interactive shell) means the batch job and your terminal cannot drift apart --
# this script previously loaded CUDA and never GROMACS, which would have failed all 864
# tasks identically.
# Pass FEP_MOCK=1 (`qsub -v VARIANT=A4V,FEP_MOCK=1 ...`) for a scheduler shake-out that
# runs the synthetic window instead of the real GPU one.
source scripts/scc_env.sh
echo "host=$(hostname) CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>} mock=${FEP_MOCK:-0}"

# Direct qsub submission must enforce the same control gate as the Snakemake DAG. Gate
# controls are allowed to produce the data needed by the barrier; every other variant
# requires an existing, passing validation report before spending GPU time.
python - "${VARIANT}" "${CONFIG}" <<'PY'
import json
import sys
from pathlib import Path

import yaml

variant, config_path = sys.argv[1:]
config_file = Path(config_path).resolve()
cfg = yaml.safe_load(config_file.read_text())
subset = cfg["validation"]["gate_subset"]
if variant in subset:
    raise SystemExit(0)

report = Path(cfg["validation"]["outputs"]["gate_report"])
if not report.is_absolute():
    report = config_file.parent.parent / report
try:
    passed = json.loads(report.read_text()).get("passed") is True
except (FileNotFoundError, json.JSONDecodeError, OSError) as exc:
    print(
        f"ERROR: cannot submit non-gate variant {variant}: validation gate report "
        f"{report} is missing or invalid ({exc}). Run the gate first.",
        file=sys.stderr,
    )
    raise SystemExit(2)
if not passed:
    print(
        f"ERROR: validation gate has not passed; refusing FEP for {variant}. "
        f"See {report}.",
        file=sys.stderr,
    )
    raise SystemExit(2)
PY

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

python -m src.fep.window \
    --variant "${VARIANT}" \
    --leg "${leg}" \
    --window "${win}" \
    --rep "${rep}" \
    --config "${CONFIG}" \
    --out "results/fep/${VARIANT}/${leg}/w${win}_r${rep}.npz"
