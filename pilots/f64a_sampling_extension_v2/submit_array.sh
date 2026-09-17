#!/bin/bash -l
# Guarded submit overrides h_rt using the recorded aggregate budget.
#$ -P rise-batteries
#$ -l gpus=1
#$ -l gpu_c=7.0
#$ -l gpu_type=L40S
#$ -l h_rt=12:00:00
#$ -pe omp 8
#$ -N f64a_v2_6ns
#$ -j y
#$ -r y
#$ -cwd
#$ -o logs/f64a_continuation/
#$ -t 1-60
#$ -tc 8
set -eo pipefail
[[ -f pilots/f64a_sampling_extension_v2/config.yaml ]] || exit 2
[[ "${SGE_TASK_ID:-undefined}" =~ ^[0-9]+$ ]] || exit 2
(( SGE_TASK_ID >= 1 && SGE_TASK_ID <= 60 )) || exit 2
source scripts/scc_env.sh
read -r n_windows n_reps < <(python -c "import yaml; c=yaml.safe_load(open('pilots/f64a_sampling_extension_v2/config.yaml')); print(c['lambda_windows'], c['replicates'])")
(( n_windows * n_reps == 60 )) || exit 2
idx=$((SGE_TASK_ID - 1))
python -m src.fep.f64a_continuation \
  --config pilots/f64a_sampling_extension_v2/config.yaml \
  run --window "$((idx / n_reps))" --rep "$((idx % n_reps))"
