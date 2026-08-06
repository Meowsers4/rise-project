# AGENTS.md — SOD1 FEP variant pipeline

Operating instructions for Codex on this repository. Read `README.md` for the
full design; this file is the short list of rules that must shape every action.

<!-- Maintainer note: block-level HTML comments are stripped before injection, so
     they cost no context. Use them for notes to humans, not to Codex. -->

## What this project is
A GPU-cluster pipeline computing ΔΔG of folding for SOD1 (UniProt P00441) variants
via alchemical FEP, validated against experimental controls, to triage
uncharacterized ALS variants. It is a hybrid methods/target computational-biophysics
project. `config/pipeline.yaml` is the single source of parameters; `data/variants.csv`
is the single source of truth for the panel.

## Non-negotiable rules
1. **Apo-first.** v1 simulates the apo, disulfide-reduced form. Do NOT switch to
   holo/metal-bound. If holo seems needed, STOP and ask — it is a separate
   parameterized arm, not a config flag or a refactor.
2. **The validation gate is go/no-go.** Do not run FEP on uncharacterized variants
   until FEP reproduces the experimental ΔΔG of the positive controls above
   `validation.min_pearson`. Never bypass or lower this threshold to "make progress".
3. **Mature numbering.** Variant names use mature (153-residue) numbering; the
   precursor has an extra N-terminal Met. The offset lives ONLY in
   `pipeline.yaml:project.mature_offset`. Never hardcode residue indices anywhere.
   `tests/test_numbering.py` guards this — keep it passing.
4. **No parameters in code.** Read everything from `config/pipeline.yaml`. If a value
   is missing or `TODO`, ask the user rather than inventing a default.
5. **Replicates and error bars are mandatory** for any free-energy result. A ΔΔG
   without an uncertainty and a cycle-closure check is not a result.

## Before writing implementation code
Resolve the open decisions in `README.md §9` WITH THE USER first — especially:
- FEP framework (openmm_perses | gromacs_pmx | amber_ti) — shapes all of `src/fep/`.
- GPU count, partition name, walltime — sizes the panel and SGE array chunking.
Ask these interactively; do not pick defaults and proceed.

## How to work here
- Prefer editing existing modules under `src/<stage>/` over adding new top-level code.
- Every stage is a Snakemake rule; wire new work into the DAG in `workflow/Snakefile`,
  don't create standalone run scripts.
- Stage 3 is submitted as SGE job arrays (`qsub -t`) on the BU SCC — one task per
  `(variant, leg, window, replicate)`. Keep that granularity — it's the whole point.
- GPUs are requested via `-l gpus=1 -l gpu_c=7.0` and assigned by the scheduler
  through `CUDA_VISIBLE_DEVICES`; never hardcode device IDs (one window = one GPU).
- SCC modules are versioned (no bare names): `miniconda/25.3.1`, `cuda/12.8` (stay on
  the 12.x line for OpenMM). These live in `config/pipeline.yaml:cluster.{conda,cuda}_module`;
  scripts read them from there (env-var override allowed).
- Checkpoint long GPU jobs; assume windows will die and must resume.
- `scripts/submit_array.sh` must run `set -eo pipefail`, NOT `set -euo pipefail`:
  GROMACS's `GMXRC` (sourced via `scc_env.sh`) reads unbound `$shell`/`$GMXLDLIB`,
  and `set -u` aborts every array task before python runs. Batch must match
  interactive.
- Never hand-edit the SCC copy of the repo. Fix files locally, commit and push,
  then `git restore` the touched paths and `git pull` on the SCC — ad-hoc SCC
  edits diverge and make pulls conflict.
- Honor `cluster.trajectory_retention` — do not persist raw frames when the config
  says estimates only (disk fills fast).

## Verification
- Run `pytest tests/` after touching numbering, panel parsing, or config loading.
- Dry-run the DAG with `snakemake -n` after changing any rule's inputs/outputs.
- When array tasks fail, grep `logs/fep/` for python `Traceback|ToolError` AND bash
  `unbound variable` errors — the latter abort before python ever runs.
- Report any command you could not run rather than assuming it passed.

## Style
- Python 3.11+, type hints, docstrings on public functions.
- Write rules as actions, not vibes: name the file/function to touch, not "be careful".
- If a rule here conflicts with something you infer from the code, surface the conflict
  to the user instead of silently choosing.
