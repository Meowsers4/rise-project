#!/bin/bash -l
#$ -P rise-batteries
#$ -l gpus=1
#$ -l gpu_c=7.0
#$ -l gpu_type=L40S
#$ -l h_rt=01:00:00
#$ -pe omp 8
#$ -N f64a_ext_light
#$ -j y
#$ -r y
#$ -cwd
#$ -o logs/f64a_extension/

set -eo pipefail
[[ -f pilots/f64a_sampling_extension_v1/config.yaml ]] || {
  echo "ERROR: submit from the repository root" >&2
  exit 2
}
source scripts/scc_env.sh
python -m src.fep.f64a_extension \
  --config pilots/f64a_sampling_extension_v1/config.yaml \
  run --window 0 --rep 0 --first-light
