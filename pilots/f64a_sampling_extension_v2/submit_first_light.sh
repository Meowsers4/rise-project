#!/bin/bash -l
# Submit only through: python -m src.fep.f64a_continuation submit --first-light
#$ -P rise-batteries
#$ -l gpus=1
#$ -l gpu_c=7.0
#$ -l gpu_type=L40S
#$ -l h_rt=01:00:00
#$ -pe omp 8
#$ -N f64a_v2_light
#$ -j y
#$ -r y
#$ -cwd
#$ -o logs/f64a_continuation/
set -eo pipefail
[[ -f pilots/f64a_sampling_extension_v2/config.yaml ]] || exit 2
source scripts/scc_env.sh
python -m src.fep.f64a_continuation \
  --config pilots/f64a_sampling_extension_v2/config.yaml \
  run --window 0 --rep 0 --first-light
