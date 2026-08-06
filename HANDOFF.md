# HANDOFF — SOD1 FEP pipeline, Stage 3 first light (2026-08-06)

Paste the block below to the next agent taking over this work.

---

You are taking over the SOD1 FEP pipeline (BU SCC, GROMACS+pmx). Repo: local
`/Users/bodebosell/sod1fep`, origin `github.com/Meowsers4/rise-project`, branch
`main`. SCC copy: `/projectnb/rise-batteries/bode/rise-project` (login `scc1.bu.edu`).
Read `README.md`, `AGENTS.md`, `CLAUDE.md` first.

**STATE:** Stage 3 works end-to-end for the first time (2026-08-06). Verified on
the SCC: `pmx mutate`, both `pdb2gmx` passes, `pmx gentop`, solvate/genion, EM,
`grompp` and GPU `mdrun` all run; a real A4V folded window completed in 892.8 s
(~15 min), `u_kn(18,151)`, provenance `gromacs_pmx`. Committed + pushed: `pmx
mutate --script` (not `-script`), no `-ignh` on the hybrid `pdb2gmx`,
`OMP_NUM_THREADS` forced to `cluster.mdrun_ntomp` on mdrun subprocesses, resume
flags only when `prod.cpt` exists, and `fcntl.flock` serializing the
per-(variant, leg) system build.

**THE BLOCKER — fixed locally, NOT yet committed/pushed:**
`scripts/submit_array.sh` ran `set -euo pipefail`; under `set -u`, GROMACS's
`GMXRC` (sourced via `scc_env.sh`) aborts every array task with
`shell: unbound variable` / `GMXLDLIB: unbound variable`, so tasks died before
python ran and no `.npz` was produced. It is now `set -eo pipefail` with a
comment saying why. `bash -n` clean, `PYTHONPATH=. pytest -q tests/` = 82 passed.
Array job **7048729** (A4V, 108 tasks) is still running-but-failing on the old
code — `qdel` it before resubmitting.

**TASKS:**

1. Commit and push the `scripts/submit_array.sh` change (already made locally).
2. On SCC: `git restore src/fep/pmx_engine.py scripts/submit_array.sh`
   (discard local hand-edits) then `git pull`.
3. Remove stale test npz so progress counts cleanly:
   `rm -f results/fep/A4V/folded/w0_r0.npz results/fep/A4V/unfolded/w0_r0.npz`
4. Submit A4V array: `cd /projectnb/rise-batteries/bode/rise-project &&
   mkdir -p logs/fep && qsub -v VARIANT=A4V scripts/submit_array.sh`
5. Monitor: progress = `find results/fep/A4V -name "w*_r*.npz" | wc -l`
   (108 = done). `qstat -u bodeb`. Check `logs/fep` for python Tracebacks AND
   bash `unbound variable` errors.
6. When 108 npz exist: `python -m src.fep.analyze --variant A4V --config
   config/pipeline.yaml --out results/fep/A4V/ddg.json`. Require converged,
   `cycle_closure_kcal < 1.0`, provenance `gromacs_pmx`.
7. Only then scale: submit the other 7 `gate_subset` variants
   [F64A, I18V, I113T, G93A, G93S, I149A, G93V], one `qsub -v VARIANT=<name>`
   each. 8 x 108 = 864 tasks total.
8. Run the validation gate (`src/analysis/validate.py`). GO/NO-GO per AGENTS.md:
   `min_pearson` 0.6 AND `max_rmse_kcal` 1.5 on the charge-neutral subset.
   Never lower thresholds.

**NON-NEGOTIABLE (AGENTS.md / CLAUDE.md):** apo-first; the gate is go/no-go;
mature numbering via `project.mature_offset`; no parameters in code;
replicates + error bars mandatory; only `provenance == gromacs_pmx` counts as
a result.

**FACTS:** folded window ~15 min/GPU, unfolded ~2-3 min; A4V array ~16 GPU-hours
(~4-8 h wall at 2-5 concurrent). `h_rt=12:00:00` is per-task, never hit. Array
tasks request `-l gpus=1 -l gpu_c=7.0 -pe omp 8`; use `CUDA_VISIBLE_DEVICES` as
given. Killed tasks resume via `-cpi`. Never hand-edit SCC files; fix locally,
commit, push, then `git restore` + `git pull` on the SCC.

**DONE WHEN:** A4V `ddg.json` is a real `gromacs_pmx` result; the gate report
(`results/validation_gate.json`) exists; `CLAUDE.md` documents all first-light
traps.
