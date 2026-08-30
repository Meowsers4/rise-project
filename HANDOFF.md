# HANDOFF — SOD1 FEP pipeline

Written for an agent with **no prior context**. Updated 2026-08-30.

Read `README.md` for the scientific design and `CLAUDE.md` for the operating rules. This
file is only: where things stand, what to do next, and what not to break.

---

## 1. Orientation in sixty seconds

The pipeline computes ΔΔG of folding for SOD1 variants by alchemical FEP (GROMACS + pmx)
on the BU SCC, and validates against experimental controls before touching uncharacterized
variants. Stage 3 (the GPU work) is submitted as SGE job arrays, one task per
(variant, leg, window, replicate).

- Local repo: `/Users/bodebosell/sod1fep` · SCC: `/projectnb/rise-batteries/bode/rise-project`
- Origin: `github.com/Meowsers4/rise-project`, branch `main`
- **Never edit files on the SCC.** Fix locally → commit → push → `git pull` on the SCC.
- Environment on the SCC is one command: `source scripts/scc_env.sh` (prompt shows
  `(sod1-fep)`). Do not `conda activate` by hand — you get the env without `GMXLIB` and
  `pdb2gmx` cannot find the pmx force field.

**Current protocol hash: `822108e9db71124d`.** 20 λ-windows, 3 replicates, 2 legs = **120
tasks per variant**, ~23 GPU-hours, ~10 h wall at current throughput. Every window records
this hash; `analyze` refuses to mix protocols within a variant and `evaluate_gate` refuses
to mix them across variants.

---

## 2. Where things stand

### Completed under the current protocol

| variant | ΔΔG | exp | closure (cap 1.0) | usable? |
|---|---|---|---|---|
| F64A | 6.94 ± 0.34 | −0.20 | **1.10** | ❌ not converged |

Full record: [`docs/f64a_20window_result.md`](docs/f64a_20window_result.md).

### Everything else must be rerun

| variant | exp ΔΔG | status |
|---|---|---|
| A4V | 1.62 | 18-window only (2.77 ± 0.30) — **invalidated** |
| G93A | 2.43 | 18-window only (1.38 ± 0.11) — **invalidated** |
| I18V | 0.37 | never run |
| I113T | 1.25 | pre-fix only |
| G93S | 3.70 | never run |
| I149A | 4.05 | never run |
| G93V | 7.00 | never run |

The gate needs **`min_gate_points: 6`** usable variants out of the 8 in
`validation.gate_subset`. F64A is already out, so **6 of these 7 must converge**. One more
failure ends the gate as specified.

Superseded results are preserved because `results/` is gitignored:
- [`docs/prefix_diagnostics.md`](docs/prefix_diagnostics.md) — the first three runs
- [`docs/postfix_18window_results.md`](docs/postfix_18window_results.md) — A4V 2.77, G93A 1.38
- [`docs/f64a_folded_leg_failure.md`](docs/f64a_folded_leg_failure.md) — why the protocol changed

---

## 3. What to do next

```bash
cd /projectnb/rise-batteries/bode/rise-project
qstat -u bodeb            # must be empty before pulling
git status --short        # investigate anything modified; never discard blind
git pull
qsub -v VARIANT=A4V scripts/submit_array.sh
```

**One array at a time.** With ~2 concurrent GPUs a second array adds no throughput and only
creates contention.

A4V first as a control: it gave 2.77 under the 18-window protocol, so it measures what the
new ladder and independent boxes change on a variant that already converged. Then G93A,
then the five untested ones.

### Monitoring

```bash
V=A4V
echo "$(date +%H:%M) | done: $(find results/fep/$V -name 'w*_r*.npz' | wc -l)/120 | live: $(find results/fep/$V -name prod.log -newermt '-2 minutes' | wc -l) | queue: $(qstat -u bodeb | grep -c sod1_fep)"
```

`live` is the honest check — a count that isn't moving with an empty queue means it stopped.

### Analysing a finished variant

```bash
source scripts/scc_env.sh
python -c '
import numpy as np, glob, sys
fs = glob.glob(f"results/fep/{sys.argv[1]}/*/w*_r*.npz")
print(len(fs), {np.load(f)["u_kn_window"].shape for f in fs}, {str(np.load(f)["protocol"]) for f in fs})' A4V
python -m src.fep.analyze --variant A4V --config config/pipeline.yaml --out results/fep/A4V/ddg.json
```

**Always run that first line.** Want `120`, `{(20, 3001)}`, one protocol hash. A run that
silently produced empty windows once got as far as a confident-looking analysis crash.

### When at least 6 have converged

```bash
python -m src.analysis.validate --gate-only --config config/pipeline.yaml \
  --out results/validation_gate.json
```

Pre-registered 2026-08-07, before any gate evaluation. **Never lower one to make progress.**

| criterion | value |
|---|---|
| `min_pearson` | 0.70 |
| `max_rmse_kcal` | 1.5 |
| `max_median_cycle_closure_kcal` | 0.75 |
| `pivot_pearson` | 0.60 — below this, reframe per README §10, do not retune |

---

## 4. Known traps — all previously paid for

| trap | what happens | guard |
|---|---|---|
| `git pull` mid-array | tasks split across two protocols | `analyze._check_single_protocol` |
| stale run dirs after a protocol change | mdrun resumes a complete old checkpoint, writes nothing, windows come out `(n,0)` | `assert_resumable`, `discard_equilibration` |
| raising `fep.replicates` | `#$ -t` must change in the same commit or every task exits 2 | `test_submit_array_task_count_matches_the_config` |
| a host that will not yield its GPU | mdrun dies "no GPU detected" | retry, then exit 99 so SGE reschedules |
| Blackwell nodes | `cudaErrorInvalidPtx` (GROMACS has no sm_120 kernels) | `#$ -l gpu_type=L40S` |
| deleting a variant directory | superseded results vanish (`results/` is gitignored) | record in `docs/` **first** |

---

## 5. Open questions — decisions, not tasks

1. **Three experimental values are suspect.** F64A (−0.20; two Phe→Ala at equal burial
   differ by 2.3 kcal/mol, and it breaks the atom-count trend within its own source),
   G93V (7.00) and G93S (3.70) (the only controls where monomer > dimer, all from
   Stathopulos 2006). F64A and G93V are the two anchors of the gate's dynamic range.
   **Do not drop a control because it produced an inconvenient FEP number** — that is the
   same category of post-hoc adjustment as lowering `min_pearson`. Check the primary
   sources first, and write the argument down before any gate evaluation.
2. **No overlap floor among the gate criteria.** G93A passed every convergence check with
   a minimum adjacent overlap of 0.017 — indistinguishable from F64A's 0.018 — while being
   1.05 kcal/mol wrong. Adding a floor only *tightens*, so it is defensible against the
   pre-registration rule, but it must be a deliberate decision, not a mid-analysis tweak.
3. **`replicates` 3 → 5** (README §9 resolved it at 5). Verified additive: the protocol
   hash is identical at 3 and 5, so r3/r4 can be added later without a rerun. `#$ -t` must
   go to 200 in the same commit.
4. **`-maxwarn 2`** is hardcoded on all four `grompp` calls and `mdout.mdp` is never
   archived, so nobody knows what grompp chose for `nstpcouple`/`nsttcouple`.
5. **Trajectory retention.** `nstxout-compressed = 0` means structural hypotheses cannot be
   tested after the fact. `prod.gro` (final frame) may be enough.

---

## 6. Not started

- **C1 (charge-changing variants)** — the strongest surviving methods claim. 17
  charge-changing positive controls carry experimental ΔΔG and would form the sub-gate's
  validation set. Needs co-alchemical counterions or Rocklin corrections. Unimplemented.
- **C3/C4** — `src/prescreen/run.py:get_backend()` raises `NotImplementedError`; the string
  "concordance" appears nowhere in the code; `data/axakova_dms.csv` does not exist. All
  CPU/data work that costs no GPU time and can proceed in parallel with the gate.
- **Negative controls** — the panel has 54 positive_control and 38 vus, and zero benign
  controls, so `classify_uncharacterized`'s "stable" label has no calibration set.
- **`rule validate` and `rule prep`** — `rule prep`'s output is read by nothing yet gates
  the GPU jobs in the DAG. Neither is on the execution path today (Stage 3 goes through
  `qsub`, not Snakemake), so this is latent rather than blocking.
- **`AGENTS.md` is a stale fork of CLAUDE.md** — an agent routed there gets the pre-audit
  rulebook.

---

## 7. The scientific state, stated plainly

Every variant run so far is wrong by ~1 kcal/mol or more, and the two most complete
diagnostics point at the same thing: **the machinery is sound and the folded leg is not
converged.** Unfolded legs reproduce to 0.05–0.24 kcal/mol with low hysteresis in every
variant; every problem lives in the folded leg.

The sharpest result is not a ΔΔG. It is that **cycle closure cannot distinguish a converged
calculation from a wrong one** — F64A's best-behaved replicate (hysteresis 0.02, the lowest
recorded anywhere in the project) is the one that disagrees with its siblings by 1.2
kcal/mol. That is a direct, checkable statement about the limitation Wells 2021 named and
could not resolve, and it is publishable whether or not the gate passes.

If the gate fails, README §10 pre-commits the project to a methods/sampling-limits result
rather than a retune. That is a real outcome, not a fallback.
