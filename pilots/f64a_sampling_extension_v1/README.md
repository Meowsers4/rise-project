# Isolated F64A 3→9 ns continuation package

This package continues the exact archived F64A folded trajectories for 6 ns. Read
[`PREREGISTRATION.md`](PREREGISTRATION.md) before submitting. The old general-purpose
`scripts/submit_array.sh` is not used.

## Deployment and provenance

Deploy this arm through Git only. Commit and push the package before touching the SCC;
do not rsync an uncommitted copy. Staging refuses untracked or modified pilot paths and
records both `git rev-parse HEAD` and the SHA-256 of `src/fep/f64a_extension.py`. Every
output NPZ carries the same identities, and analysis rejects a mixed set.

From the repository root on the SCC:

```bash
git pull --ff-only
git status --short -- src/fep/f64a_extension.py pilots/f64a_sampling_extension_v1 \
  workflow/Snakefile config/pipeline.yaml
mkdir -p logs/f64a_extension
source scripts/scc_env.sh

python -m src.fep.f64a_extension \
  --config pilots/f64a_sampling_extension_v1/config.yaml verify-tools

# The SCC retained source is the historical results tree.  The separately recorded
# archive_2026-09-11 path is not present on the SCC; validation below proves this copy
# against the frozen off-cluster archive manifest before staging it.
SOURCE=/projectnb/rise-batteries/bode/rise-project/results/fep/F64A/folded
python -m src.fep.f64a_extension \
  --config pilots/f64a_sampling_extension_v1/config.yaml \
  validate-source --source-folded "$SOURCE"

python -m src.fep.f64a_extension \
  --config pilots/f64a_sampling_extension_v1/config.yaml \
  stage --source-folded "$SOURCE" --first-light

qsub pilots/f64a_sampling_extension_v1/submit_first_light.sh
```

After that job finishes, inspect its log and confirm the isolated output exists:

```bash
qstat -u "$USER"
tail -n 80 logs/f64a_extension/f64a_ext_light.*
python - <<'PY'
import numpy as np
p = 'results/pilots/f64a_sampling_extension_v1/first_light/fep/F64A/folded/w0_r0.npz'
with np.load(p) as d:
    print(d['u_kn_window'].shape, d['protocol'], np.isfinite(d['u_kn_window']).all())
PY
```

The `git status` command above must print nothing. Source validation checks the exact
0–3500 ps time grid and re-parses 500–3500 ps from every archived `dhdl.xvg`; those
energies must exactly equal the checksum-frozen NPZ before any GPU work is staged.

The expected first-light shape is `(20, 3011)` and the finite flag must be `True`. Have
the package enforce the same check; only then stage production and submit the bounded array:

```bash
python -m src.fep.f64a_extension \
  --config pilots/f64a_sampling_extension_v1/config.yaml check-first-light

python -m src.fep.f64a_extension \
  --config pilots/f64a_sampling_extension_v1/config.yaml \
  stage --source-folded "$SOURCE"

qsub pilots/f64a_sampling_extension_v1/submit_array.sh
```

Do not pull or change config/code while the array is active. At most eight L40S jobs run
simultaneously because the script uses `-tc 8`; each task requests one GPU and eight CPU
slots. Failed GPU placement exits 99 and is reschedulable.

After all 60 NPZs exist, run the CPU analysis once:

```bash
find results/pilots/f64a_sampling_extension_v1/production/fep/F64A/folded \
  -maxdepth 1 -name 'w*_r*.npz' | wc -l
python -m src.fep.f64a_extension \
  --config pilots/f64a_sampling_extension_v1/config.yaml validate-outputs
python -m src.fep.f64a_extension \
  --config pilots/f64a_sampling_extension_v1/config.yaml analyze
```

The expected count is 60. `analyze` repeats the complete output validation itself and
refuses to overwrite an existing result.

## Snakemake boundary

`workflow/Snakefile` contains `f64a_extension_first_light`, `f64a_extension_window`, and
`f64a_extension_analysis`, so the dependency graph is explicit. The pilot is intentionally
not added to the default target: `rule all` remains the frozen CPU forensic release.
Archive staging remains an explicit operation because it creates a hashed copy of an external
immutable archive. GPU dispatch uses the dedicated SGE scripts, matching the existing Stage 3
one-window-per-task scheduler contract and its `-tc 8`/exit-99 rescheduling behavior.
