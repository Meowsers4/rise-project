# Cluster module notes -- BU SCC (SGE/OGS)

Environment setup that lives outside `environment.yml`: SCC modules and licensed
tools. Keep this in sync with `env/environment.yml` and `config/pipeline.yaml`.

## Conda + CUDA

```bash
module load miniconda        # or the site's anaconda module
module load cuda             # match the CUDA the conda openmm build expects
conda activate sod1-fep      # created from env/environment.yml
```

## GPU requests (Stage 3, `scripts/submit_array.sh`)

- Request GPUs by **compute capability**, not model:
  `-l gpus=1 -l gpu_c=7.0` (draws from A40 / A6000 / L40S / V100 / A100 / H200 ...).
- One window = one GPU. The scheduler sets `CUDA_VISIBLE_DEVICES`; use it as-is.
- Keep `-l h_rt` <= 12h to stay eligible for buy-in nodes (`cluster.max_walltime_hr`).
- `-P rise-batteries` charges the correct project.

## Licensed tools (Stage 2 prescreen) -- NOT conda packages

- **FoldX** -- academic license; install/module-load separately, expose on `PATH`.
- **Rosetta** (`cartesian_ddg`) -- license required; module-load or point the
  prescreen wrapper at the binary. Record the path in config, do not hardcode.

> TODO: fill in the exact SCC module names/versions once confirmed on the cluster.
