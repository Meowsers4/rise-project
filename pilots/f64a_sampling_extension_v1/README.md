# Isolated F64A 3→9 ns continuation package

**Completed 2026-09-16:** production array 7589742 produced 60 validated outputs.
The mean folded-leg estimate shifted by −0.933 kcal/mol and repeat range fell 71%,
without establishing convergence at 9 ns. The supplied SCC report and interpretation
are preserved in [the result record](../../docs/f64a_extension_results.md).
GPU instructions below are execution history, **not** a request to resubmit.

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
energies must exactly equal the checksum-frozen NPZ before any GPU work is staged. After
restart, 500–3499 ps must remain exact. GROMACS may regenerate the 3500 ps checkpoint
boundary; that value is audited but excluded in favor of the frozen source column, and
only exact-grid continuation samples beginning at 3501 ps are admitted.

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

## Authorized CPU follow-up (not yet run on real extended data)

The [block-analysis design](../../docs/f64a_block_analysis_design.md) is fixed after the
primary report, before follow-up execution. It compares disjoint 3 ns blocks and the two
halves of the final block, retaining each window's t0/g/selected counts and each signed
adjacent discrepancy. The original GPU config and primary JSON are not modified.

After job 7589742 has completely left the queue, deploy the committed CPU package and run:

```bash
source scripts/scc_env.sh
python -m src.analysis.f64a_blocks \
  --config pilots/f64a_sampling_extension_v1/cpu_blocks.yaml
```

The new immutable output is `production/analysis/f64a_blocks_v1.json` below the pilot root.
The analyzer requires the original cluster report to agree with the preserved supplied
report, fingerprints all 60 production matrices, and checks reproduction of its first
and final 3 ns blocks before writing any result.

For explicit CPU-only Snakemake execution, exclude the historical GPU producers:

```bash
snakemake f64a_extension_blocks --snakefile workflow/Snakefile --cores 1 \
  --allowed-rules f64a_extension_blocks
```

Missing completed inputs are fatal; this command cannot schedule GPU work. Do not rerun
the old primary analysis or first-light gate under the later CPU code identity.

## CPU selection sensitivity (authorized 2026-09-16)

The completed disjoint-block report and its interpretation are preserved in
[`docs/f64a_block_analysis_results.md`](../../docs/f64a_block_analysis_results.md).
The next CPU-only package compares adaptive trimming with fixed 0% and 25% cutoffs
on the final two halves. See
[`docs/f64a_selection_sensitivity_design.md`](../../docs/f64a_selection_sensitivity_design.md).
Run after transferring the committed update, with no array in flight:

```bash
source scripts/scc_env.sh
python -m src.analysis.f64a_sensitivity \
  --config pilots/f64a_sampling_extension_v1/cpu_sensitivity.yaml
```

This writes `production/analysis/f64a_selection_sensitivity_v1.json` under the pilot
root. It never launches GPU work or overwrites the earlier reports. Keep the original
`cpu_blocks.yaml`, GPU pilot config and base pipeline config unchanged; do not rerun
the completed block command after updating analysis code.
