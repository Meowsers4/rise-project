#!/bin/bash -l
# F64A folded only: 20 lambda states x 3 retained paths = 60 tasks.
#$ -P rise-batteries
#$ -l gpus=1
#$ -l gpu_c=7.0
#$ -l gpu_type=L40S
#$ -l h_rt=12:00:00
#$ -pe omp 8
#$ -N f64a_ext_6ns
#$ -j y
#$ -r y
#$ -cwd
#$ -o logs/f64a_extension/
#$ -t 1-60
#$ -tc 8

set -eo pipefail
: "${SGE_TASK_ID:?submit with qsub; do not execute the task body directly}"
[[ -f pilots/f64a_sampling_extension_v1/config.yaml ]] || {
  echo "ERROR: submit from the repository root" >&2
  exit 2
}
source scripts/scc_env.sh

# Production is scientifically conditional on a real continuation first light, not merely
# on the files having been staged. Every task checks the same tiny NPZ before requesting MD.
python -m src.fep.f64a_extension \
  --config pilots/f64a_sampling_extension_v1/config.yaml check-first-light

N_WINDOWS=$(python -c "import yaml; print(yaml.safe_load(open('pilots/f64a_sampling_extension_v1/config.yaml'))['lambda_windows'])")
N_REPS=$(python -c "import yaml; print(yaml.safe_load(open('pilots/f64a_sampling_extension_v1/config.yaml'))['replicates'])")
N_TASKS=$((N_WINDOWS * N_REPS))
if [[ "${SGE_TASK_LAST:-$N_TASKS}" -ne "${N_TASKS}" ]]; then
  echo "ERROR: array size ${SGE_TASK_LAST:-?} != windows*repeats=${N_TASKS}" >&2
  exit 2
fi

idx=$((SGE_TASK_ID - 1))
rep=$((idx % N_REPS))
window=$((idx / N_REPS))
echo "host=$(hostname) gpu=${CUDA_VISIBLE_DEVICES:-unset} window=${window} rep=${rep}"
python -m src.fep.f64a_extension \
  --config pilots/f64a_sampling_extension_v1/config.yaml \
  run --window "${window}" --rep "${rep}"
