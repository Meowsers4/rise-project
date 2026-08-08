# HANDOFF — SOD1 FEP pipeline

Current state and the ordered work queue. Updated 2026-08-08.
Supersedes the 2026-08-06 first-light handoff (that blocker is long fixed).

## Where things stand

Stage 3 works end to end. One variant is analysed, two are on the cluster.

| | value | note |
|---|---|---|
| A4V ΔΔG | **+4.02 ± 0.23** kcal/mol | experiment (apo monomer) is **1.62** — off by +2.4 |
| A4V cycle closure | 0.63 (cap 1.0) | reports `converged: true` |
| A4V min adjacent overlap | 0.040 | marginal; low pairs move between replicates → noise, not a coarse ladder |
| per-replicate | 3.75 / 3.84 / 4.48 | spread 0.73 vs a 0.23 error bar |

The variance and hysteresis live entirely in the **folded** leg (−12.56/−12.49/−11.91,
hyst 0.63/0.15/0.32); the unfolded leg is pinned (−16.31/−16.33/−16.39, hyst 0.10 each).

**Running now:** G93A and I113T, 108 tasks each, 3 replicates, **pre-fix protocol**.
These are diagnostics, not gate data. They discriminate three hypotheses:

| | constant offset | dimer-like physics | actually correct |
|---|---|---|---|
| G93A (exp 2.43) | ~4.8 | ~3.0 | 2.43 |
| I113T (exp 1.25) | ~3.7 | ~2.5 | 1.25 |

A4V's apo-**dimer** experimental value is 4.31, and we produced 4.02. A4 sits in the
dimer interface, prep carves out chain A, and 3.5 ns is not enough for the exposed face
to relax. I113T is the discriminating variant (monomer/dimer gap 1.23).

**Do not `git pull` on the SCC until both arrays finish.** Tasks read the repo and config
when each task *starts*, so pulling mid-array splits a variant across two protocols.
`analyze.py` now refuses to combine them, but only for post-fix runs.

## Queue

### 1. When the arrays finish — pull and verify

Already committed and pushed (`da09142`), waiting on a pull:

- `sc-coul = yes` — was defaulting to `no`, a silent deviation from pmx's protocol
- `frames_per_window` 150 → 3000 (`nstdhdl` 20 ps → 1 ps): 3001 usable records/window, was 151
- HIS aliasing (HID/HIE/HIP→HIS) — unblocks H43R, H46R, H71Q, H110Y, H120L
- `-cpt` from `fep.checkpoint_interval_min` — windows previously had no checkpoint
- `.mdp` protocol fingerprint in every NPZ; `analyze` refuses to mix protocols
- conservative decorrelation (max of g from total potential and from neighbour Δu)

After pulling, confirm on one window before committing a full array:
`grep sc-coul results/fep/<v>/<leg>/w0_r0/prod.mdp` and check `u_kn` is `(18, 3001)`.

### 2. Decisions needed before the next production run

| Item | Cost / risk | Recommendation |
|---|---|---|
| **Box padding** `structure.solvent_padding_nm: 1.0` against a 1.0 nm cutoff = **zero** minimum-image margin. Comment says it was written for Stage 1's `Modeller.addSolvent`, now reused for `editconf -d`. | 1.2 nm adds ~20% atoms → ~20% per window | raise to 1.2 |
| **`-maxwarn 2`** hardcoded on all four grompp calls; `mdout.mdp` never archived, so we cannot see what grompp chose for `nstpcouple`/`nsttcouple` | may surface warnings needing interpretation | archive `mdout.mdp`, tighten to `-maxwarn 0` and fix what appears |
| **`cycle_closure_kcal` is not a cycle closure** — it is within-ladder Zwanzig hysteresis, and reads "converged" on a result 2.4 off. It is a *pre-registered* gate criterion. | changing a pre-registered threshold is exactly what pre-registration forbids | ADD an overlap floor as an extra criterion; do not redefine the existing one |
| **`replicates` 3 → 5** (README §9) | +72 tasks/variant, additive not a rerun | do it between arrays — never mid-flight, the `SGE_TASK_LAST` guard aborts every unstarted task |
| **Independent replicate systems** — `rep` is not in the system-build seed, so all replicates share one box and one EM structure; only velocities differ | the expensive fix: real NPT equilibration per replicate | decide from the G93A/I113T outcome |
| **Trajectory retention** — `nstxout-compressed = 0` blocks testing the interface hypothesis. Conflicts with `cluster.trajectory_retention: estimates_only` | disk | check whether `prod.gro` (final frame) is already enough |

### 3. Then re-run the gate

A4V, G93A and I113T were all produced under the pre-fix protocol and cannot be mixed with
post-fix windows. All gate variants must share one protocol.

- gate_subset is 8: F64A, I18V, I113T, A4V, G93A, G93S, I149A, G93V
- consider adding **L38V, G41S, V148G** — already in the panel, charge-neutral, with
  experimental ΔΔG. Takes Wells-10 overlap from 4 to 7. Caveat: the Wells-7 span only
  1.25–3.70 kcal/mol, so report paired agreement (MUE), not correlation alone
- at 2 concurrent GPUs, 11 variants × 5 replicates ≈ 275 GPU-h ≈ 5.7 days.
  Check `qquota -u bodeb` first — concurrency dominates every scope decision

### 4. In parallel — costs no GPU time

- **Axakova 2025 DMS** → `data/axakova_dms.csv`, populate `axakova_class`. Unblocks C3 and
  answers the README §10 pivot question (how many of the 38 VUS are already resolved)
- **AlphaMissense** — precomputed table, 92 lookups
- `src/analysis/concordance.py` + Snakemake rule + `validation.outputs.concordance`
  (the string "concordance" currently appears nowhere in the code)
- `src/prescreen/run.py:get_backend()` still raises `NotImplementedError` — parsers exist,
  nothing can generate input for them

## Standing rules

- Never hand-edit the SCC copy. Fix locally → commit → push → `git restore` + `git pull`.
- Never lower a pre-registered threshold (`min_pearson` 0.70, `max_rmse_kcal` 1.5,
  `max_median_cycle_closure_kcal` 0.75). `pivot_pearson` 0.60 means reframe, not retune.
- Only `provenance == gromacs_pmx` with a single `protocol` counts as a result.
