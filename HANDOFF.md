# HANDOFF — SOD1 FEP pipeline

Written for an agent with **no prior context**. Updated 2026-09-14.

Read `README.md` for the scientific design and `CLAUDE.md` for the operating rules. This
file is only: where things stand, what to do next, and what not to break.

> **G93A SS DIAGNOSTIC COMPLETE — 2026-09-14:** All 120 windows passed inventory checks at
> protocol `822108e9db71124d`, and the direct topology gate found C57-C146 in all three
> independently built folded systems. The pre-registered primary endpoint is folded ΔG:
> **9.7493 ± 0.0547 kcal/mol**, a shift of **−0.0169 kcal/mol** from the frozen 2SH mean
> (9.7662). This is squarely inside the pre-declared `|Δ| < 0.3` negligible band: the
> experimental/reference-state mismatch is real but does not explain G93A's ~1.17 kcal/mol
> underprediction. Secondary ΔΔG is **1.2698 ± 0.0618 kcal/mol**; closure is 0.1247,
> minimum adjacent overlap 0.030062, and the run is converged. Preserve the SS tree as
> `G93A_SS_diagnostic` with its manifest, restore `G93A_2SH_baseline` to `G93A`, and leave
> `diag/g93a-ss` without merging it. Full record: `docs/prereg_g93a_disulfide_diagnostic.md`.

> **Read [`docs/stop_rule_reanalysis.md`](docs/stop_rule_reanalysis.md) first (2026-09-12).**
> Every documented estimate below has now been re-derived from the raw archived windows and
> reproduces to floating-point rounding, so the review's "unverified against raw output"
> caveat is discharged — with one exception (A4V's gate point) and one correction (the F64A
> "lowest hysteresis anywhere" claim is false; see §7 below).
>
> **Then [`docs/post_pivot_review.md`](docs/post_pivot_review.md) (2026-09-11).** An
> external review corrects several claims in this file — the protocol hash is now per-leg,
> the ~23 GPU-h/variant figure does not reconcile with window timings, and the F64A
> uncertainty argument in §7 does not hold as stated. It also recommends against gate
> attempt 2. Nothing below has been edited pending the user's decision.

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
V=A4V
python - "$V" <<'EOF'
import glob, sys, numpy as np
for leg in ("folded", "unfolded"):
    fs = glob.glob(f"results/fep/{sys.argv[1]}/{leg}/w*_r*.npz")
    print(leg, len(fs),
          {np.load(f)["u_kn_window"].shape for f in fs},
          {str(np.load(f)["protocol"]) for f in fs})
EOF
python -m src.fep.analyze --variant "$V" --config config/pipeline.yaml --out results/fep/"$V"/ddg.json
```

**Always run that first block, and read it PER LEG.** Want 60 windows, one shape and one
protocol hash *within* each leg. A run that silently produced empty windows once got as far
as a confident-looking analysis crash.

The two legs may legitimately differ, and under the current config they do --
`ns_per_window` and `equilibration_ns` are per-leg, so folded and unfolded hash differently
by design and `analyze` records the pair as `folded=...|unfolded=...`.
[`_check_single_protocol`](src/fep/analyze.py#L279) enforces one protocol per leg, not per
variant. Expected shapes follow from `fep.frames_per_window` and the leg's sampling:

| leg | ns / equil | nstdhdl | records | discarded | kept |
|---|---|---|---|---|---|
| folded | 9.0 / 2.0 | 1500 | 3667 | 667 | **(20, 3000)** |
| unfolded | 3.0 / 0.5 | 500 | 3501 | 500 | **(20, 3001)** |

Everything run before `7db3095` is `(20, 3001)` on both legs at hash `822108e9db71124d`.
Inventory of what is on disk:
[`docs/raw_result_reconciliation.md`](docs/raw_result_reconciliation.md).

**Reconciling an archived result: do not write to `results/fep/<V>/ddg.json`.** That
overwrites the record you are checking against. Send `--out` somewhere else and diff.

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

1. **Reference-state mismatch — see [`docs/reference_state_audit.md`](docs/reference_state_audit.md)
   (2026-09-12), which supersedes the "suspect values" framing that stood here.** Every gate
   control was measured with the **Cys57–Cys146 disulfide intact**; we simulate the reduced
   (2SH) form. Separately, G93S and G93V are **DSC measurements on the apo dimer**, not
   apo-monomer as `variants.csv` labels them — which explains the "monomer > dimer" oddity
   without impugning the measurements. A4V/G93A/I113T re-derive exactly from Lindberg 2005, so
   the panel arithmetic is sound. **All six values with a reachable primary table now reproduce
   exactly** (Lindberg 2005 and Nordlund & Oliveberg 2006, the latter closed 2026-09-12 via
   PMC1502438). F64A's −0.20 is a real measured ΔG *above* pWT (3.07 vs 2.87), so "the
   experiment is wrong" is no longer available; and I149A's apo monomer has negative absolute
   stability (ΔG −1.18), making its 4.05 an extrapolation. Only G93S/G93V remain unchecked —
   Stathopulos is not in PMC and needs the BU library.
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
calculation from a wrong one** — F64A's best-behaved replicate (hysteresis 0.023, lowest
~20% of the archive at rank 10/48) is the one that disagrees with its siblings by 1.20
kcal/mol, the largest such disagreement in the dataset. Over all 48 (variant, leg, replicate)
records, folded-leg hysteresis barely tracks independent-box disagreement at all (Pearson
+0.072, n=24). Re-derived from raw windows 2026-09-12:
[`docs/stop_rule_reanalysis.md`](docs/stop_rule_reanalysis.md).

It was previously written here that r1's hysteresis was "the lowest recorded anywhere in the
project". **That is false** — it is 10th of 48 — and the correction is the reason §6 of that
document exists. Note also the review's caution: two agreeing replicates are not ground truth,
so this shows the diagnostic is uninformative, not which replicate is right.

If the gate fails, README §10 pre-commits the project to a methods/sampling-limits result
rather than a retune. That is a real outcome, not a fallback.
