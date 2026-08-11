# CLAUDE.md — SOD1 FEP variant pipeline

Operating instructions for Claude Code on this repository. Read `README.md` for the
full design; this file is the short list of rules that must shape every action.

<!-- Maintainer note: block-level HTML comments are stripped before injection, so
     they cost no context. Use them for notes to humans, not to Claude. -->

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
   `validation.min_pearson` AND within `validation.max_rmse_kcal` AND with a median
   |cycle closure| under `validation.max_median_cycle_closure_kcal`. All three are
   **pre-registered** (fixed 2026-08-07, before any gate evaluation) — never bypass,
   lower, or "temporarily relax" one to make progress; that destroys the only thing
   pre-registration buys. `validation.pivot_pearson` is a *different* number: falling
   below it means reframe per README §10, not retune. The gate subset is
   **charge-neutral only**: a net-charge change carries a PME finite-size artifact that
   does not cancel between the folded and unfolded legs (different box sizes), and the
   pmx engine has no counterion co-alchemy. Charge-changing variants get their own
   **sub-gate** once that is solved — see rule 6.
3. **Mature numbering.** Variant names use mature (153-residue) numbering; the
   precursor has an extra N-terminal Met. The offset lives ONLY in
   `pipeline.yaml:project.mature_offset`. Never hardcode residue indices anywhere.
   `tests/test_numbering.py` guards this — keep it passing.
4. **No parameters in code.** Read everything from `config/pipeline.yaml`. If a value
   is missing or `TODO`, ask the user rather than inventing a default.
5. **Replicates and error bars are mandatory** for any free-energy result. A ΔΔG
   without an uncertainty and a cycle-closure check is not a result. `fep.replicates`
   is pending a 3 → 5 raise (README §9). Do NOT raise it while an array is in flight:
   every task re-reads this config and `submit_array.sh` aborts when
   `legs*windows*replicates` disagrees with `SGE_TASK_LAST`. Raise it between arrays;
   r3/r4 are then additive (72 new tasks per variant), not a rerun.
6. **Four claims are load-bearing** (README §2.3): C1 charge-changing coverage,
   C2 apo-2SH convergence, C3 DMS-independent VUS triage, C4 concordance/discordance.
   Weakening one is a scope change needing the user's sign-off, not a refactor. Two
   consequences bind day-to-day work: convergence diagnostics must be **written at run
   time** to `results/convergence/<variant>.json` (C2 cannot be claimed retroactively
   from uninstrumented runs), and charge-changing variants are no longer merely deferred
   — covering them is the project's strongest methods claim. Never write "first",
   "novel FEP protocol", or "we show SOD1 variants are destabilizing": all three are
   occupied by Wells 2021, which used **this same toolchain** (README §11.2).

## Stage 3 stack — RESOLVED, do not relitigate
**GROMACS + pmx.** Not OpenMM+Perses. A research spike (2026-07-21) proved perses
hard-requires OpenEye for protein-mutation atom mapping in every version — including the
"OpenEye-free" lower-level path — confirmed by a GPU smoke test dying at
`_construct_atom_map -> import openeye.oechem`. The SCC has no OpenEye and no licence.
`src/fep/perses_engine.py` is kept only as dead reference; the live engine is
`src/fep/pmx_engine.py`, selected by `fep.framework: gromacs_pmx` through the dispatch
map in `src/fep/window.py`.

Per (variant, leg) the engine runs: `pdb2gmx` → `pmx mutate` → `pdb2gmx` → `pmx gentop`
→ box/solvate/genion → EM, cached behind a `SYSTEM_READY` marker; then per window a
generated `.mdp`, `grompp`, `mdrun`, and `alchemlyb` on `dhdl.xvg` → the same
`u_kn_window` NPZ schema `src/fep/analyze.py` has always consumed. Force field is
`amber99sb-star-ildn-mut` (pdb2gmx finds it only via `GMXLIB`).

**pmx must be new enough to read a PDB.** Require a revision at or after upstream
`866f34cf0` (2022-03-30, "do not relabel chain IDs if the chain already has an ID").
In `2.0+38.ga2311b9` the `self.atoms.append(a)` in `Model.__readPDBTER` is mis-indented
into an inner `if (a.chain_id==' ') or ...` branch, so a PDB whose atoms carry a non-blank
chain id — i.e. anything `pdb2gmx` writes — parses to ZERO atoms. `make_chains()` then
appends its still-`None` accumulator and dies with `'NoneType' object has no attribute
'model'`. Blanking the chain id does NOT work around it: the same block sets
`bNewChain = False`, so only the first atom would ever land. Verified on the SCC 2026-08-01.

**The force-field directory name is version dependent** — `data/mutff45` on 2.0-era
snapshots (where a flat legacy `data/mutff` also exists and is NOT the one to use),
`data/mutff` on current develop, which dropped mutff45 and is now the maintained set.
Never hardcode either: `pmx_engine.pmx_ff_dir()` probes both and picks whichever holds
`fep.pmx_forcefield`, and `scripts/scc_env.sh` calls that same function to set `GMXLIB`.

**pdb2gmx order matters:** pmx's own help requires its input to have come from pdb2gmx.
Feeding it a raw PDBFixer file makes pmx parse zero residues and die in `make_chains()` —
the SAME traceback as the stale-pmx bug above, so check the pmx revision before assuming
it is an input-ordering problem.

## Invariants GROMACS will silently undo
Stage 1 enforces these on an OpenMM topology; the GROMACS path needs its OWN guard:
- **Reduced disulfide.** `pdb2gmx` re-detects Cys57–Cys146 by SG–SG distance and rebuilds
  it. The engine answers `-ss` prompts "n" per cysteine (`fep.keep_disulfide_reduced`).
- **Which residue gets mutated.** PDBFixer/pdb2gmx RENUMBER. In the capped tripeptide the
  site is the *middle* residue (3), not `mature_pos` — sending `mature_pos` mutates a
  neighbour and yields a plausible, wrong ΔΔG. `pmx_engine.mutation_resid()` reads the id
  from the actual file and verifies the residue name against the panel's wild type.

## First-light traps — the first real GROMACS+pmx window (2026-08-06)
Until 2026-08-06 no Stage 3 window had ever run; the smoke test found these in
order. All are encoded in `pmx_engine.py` / `submit_array.sh` — do not
reintroduce them:
1. `pmx mutate` takes `--script`, not `-script`. A single dash is silently
   dropped by argparse, pmx falls into its interactive residue picker and blocks
   on stdin (`Enter residue number:`). The engine captures stdout, so it looks
   like a hang.
2. No `-ignh` on the hybrid `pdb2gmx` (pmx's own tutorial omits it). With it,
   pdb2gmx discards pmx's hydrogens and rebuilds them from `mutres.hdb`, which
   has no entry for hybrid residues (A2V, …) → "N missing atoms" fatal. Keep
   `-ignh` only on the WT pass.
3. `OMP_NUM_THREADS` must equal `-ntomp` for mdrun (GROMACS ≥ 2025). SCC shells
   export `OMP_NUM_THREADS=1`. Set it on the mdrun subprocess env from
   `cluster.mdrun_ntomp`; exporting it in `scc_env.sh` does not stick (later
   module/GMXRC sourcing resets it).
4. Resume flags are conditional: `-cpi X.cpt -append` fails on a fresh run. Pass
   them only when the checkpoint exists.
5. System build is serialized with an flock — all array tasks of a (variant, leg)
   start at once and would `rmtree` each other's half-built system.
6. `set -u` kills array tasks: GROMACS's `GMXRC` references unbound variables
   (`$shell`, `$GMXLDLIB`), so `submit_array.sh` must use `set -eo pipefail`
   (batch must match interactive).

Cost/monitoring: a real folded window is ~15 min on one GPU (A4V w0_r0 = 892.8 s);
unfolded ~2–3 min. An A4V array (108 windows) is ~16 GPU-hours. `h_rt=12:00:00`
is per-task, never hit. Progress = count `results/fep/<variant>/**/w*_r*.npz`.
SGE log files appear when a task STARTS; grep them for python Tracebacks AND bash
`unbound variable` errors.

## Minimisation must be lambda-aware
`grompp` without `free-energy = yes` builds the system from **A-state parameters alone**.
In a pmx hybrid topology the appearing atoms are dummies in state A — no charge, no LJ —
so minimisation relaxes real atoms *into* them. Switching them on at a high lambda then
explodes the system on step 1: LINCS reports bond lengths in the millions and the GPU
faults with `cudaErrorIllegalAddress` (a symptom, not the cause — do not go looking for a
CUDA bug). A4V folded w17 (λ=1.0) died this way in all three replicates while every
unfolded window survived, because position 4 is buried and the tripeptide is not.

`_minim_mdp(cfg, window)` therefore emits the free-energy + soft-core block, and
`run_pmx_window` minimises **per window at that window's λ** before production, starting
from the shared `system/em.gro`. The shared build still calls `_minim_mdp(cfg)` with no
window — solvation and ions only care about the physical topology.
`tests/test_fep.py:test_window_minimisation_is_lambda_aware` guards this.

Corollary: a window that merely *survives* is not evidence the protocol is right. The
high-λ windows that did not crash were started from the same A-state-minimised structure.

## Never `git pull` the SCC while an array is in flight
Array tasks read the repo code and `config/pipeline.yaml` when EACH TASK STARTS, not at
submit time. Pulling mid-array gives tasks 1..n one protocol and n+1..108 another, and
the result is a complete-looking dataset that MBAR will happily combine. Pull only
between arrays. `run_pmx_window` now hashes the production `.mdp` (minus seed and
`init-lambda-state`) into each window's NPZ as `protocol`, and
`analyze._check_single_protocol` refuses to analyse a variant whose windows disagree —
the settings-level counterpart to the engine-level provenance rule below.

## Results must carry provenance
Every window records the engine that produced it; `analyze.py` refuses to mix engines or
accept an unlabelled window; `validate.py` gates only on `provenance == fep.framework`.
This exists because fabricated ΔΔG files (experimental values + noise) once produced a
`validation_gate.json` reading "passed, pearson 0.994". Never hand-write a result file.

## How to work here
- Prefer editing existing modules under `src/<stage>/` over adding new top-level code.
- Every stage is a Snakemake rule; wire new work into the DAG in `workflow/Snakefile`,
  don't create standalone run scripts.
- Stage 3 is submitted as SGE job arrays (`qsub -t`) on the BU SCC — one task per
  `(variant, leg, window, replicate)`. Keep that granularity — it's the whole point.
- GPUs are requested via `-l gpus=1 -l gpu_c=7.0` and assigned by the scheduler
  through `CUDA_VISIBLE_DEVICES`; never hardcode device IDs (one window = one GPU).
- `source scripts/scc_env.sh` is the ONE way to get a working shell (interactive or
  batch); `submit_array.sh` sources the same file so the two cannot drift. It bootstraps
  Lmod, loads the module chain, activates conda, and exports `PYTHONPATH`/`GMXLIB`.
- SCC modules are versioned (no bare names) and live in `config/pipeline.yaml:cluster.*`.
  Three traps, all already encoded in config — don't rediscover them:
  1. `module load gromacs/2025.3` gives the **CPU-only** build. The CUDA build is a
     separate tree (`install_gpu`), selected via `cluster.gmx_gmxrc`. It needs
     `cuda/12.8 intel/2024.0 gcc/12.2.0` at runtime (it links MKL), not openmpi.
  2. `/etc/profile.d` skips Lmod in **non-interactive** shells (`qrsh ... bash -lc`, SGE
     tasks), so `module` is undefined. `cluster.lmod_init` pins
     `/usr/local/lmod/8.7.12/init/bash`; a stale broken 7.8 tree sorts ahead of it in globs.
  3. SCC home has a 10 GB quota — conda envs must live under `/projectnb`.
  4. `-l gpu_c=7.0` is a **minimum** (SGE relation `<=`), so it also admits cards NEWER
     than the GROMACS build. `scc-702`'s RTX PRO 6000 Blackwell (compute 12.0) has no
     sm_120 kernels in GROMACS 2025.3/CUDA 12.8, so mdrun starts, JIT-compiles the
     embedded PTX, and dies with `CUDA error #218 (cudaErrorInvalidPtx)`. `submit_array.sh`
     therefore also pins `-l gpu_type=L40S` (32 nodes; all timings are measured there).
     Do NOT add this error to the mdrun retry signatures — it is deterministic per node,
     so a retry on the same host fails identically and hides a placement problem.
- Checkpoint long GPU jobs; assume windows will die and must resume.
- Honor `cluster.trajectory_retention` — do not persist raw frames when the config
  says estimates only (disk fills fast).

## Verification
- Run `pytest tests/` after touching numbering, panel parsing, or config loading.
- Dry-run the DAG with `snakemake -n` after changing any rule's inputs/outputs.
- Report any command you could not run rather than assuming it passed.
- When diagnosing array failures, grep `logs/fep/` for both `Traceback|ToolError`
  and `unbound variable` — bash env errors abort before python ever runs.

## Style
- Python 3.11+, type hints, docstrings on public functions.
- Write rules as actions, not vibes: name the file/function to touch, not "be careful".
- If a rule here conflicts with something you infer from the code, surface the conflict
  to the user instead of silently choosing.
