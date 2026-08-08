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

**CORRECTION (2026-08-08).** An earlier version of this file said A4 sits in the dimer
interface and named I113T the discriminating variant. Both are wrong — measured from
`data/structures/3ECU.pdb`, chain A against its true partner chain B:

| residue | min any-atom → B | min side-chain → B | partner atoms <6 Å of side chain |
|---|---|---|---|
| A4 | 6.38 Å | **7.31 Å** | **0** |
| I113 | 3.75 Å | **3.75 Å** | **15** |
| G93 | 17.17 Å | 17.17 Å | 0 |
| I149 | 4.16 Å | 6.69 Å | 0 |
| F64 / I18 | 7.32 / 6.66 Å | 7.32 / 8.98 Å | 0 |

A4V is a dimer-**stability** variant in the literature, but the CB that grows into a Val
has no partner contact here, so a freshly exposed interface cannot directly explain a
+2.4 kcal/mol error at that side chain. The 4.02-vs-4.31 near-match is likely coincidence.

**I113T is the confounded variant, not the discriminating one** — its side chain is
directly desolvated by the partner subunit, so its monomer/dimer gap is entangled with
the very hypothesis it was meant to test. **G93A is the clean discriminator** (17 Å from
the partner: "dimer-like physics" predicts nothing there, so if G93A also comes back
~2.4 high, the cause is not the interface).

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

### 3b. Found by the 2026-08-08 review, still open

Ranked. Fixed items are in the changelog below, not here.

1. **`rule validate` rewrites the file it declares as an input.** `validate.py:main()`
   rewrites `results/validation_gate.json` even without `--gate-only`, while the Snakefile
   declares it as an input of `rule validate`. That makes `ddg_map.csv` permanently stale
   and, via `gate_dependency`, invalidates all 864 VUS window jobs. Currently masked only
   because `prescreen` raises `NotImplementedError` so the rule never runs. Fix needs a
   decision on gate semantics — do not rush it.
2. **G93V's 7.00 kcal/mol is the gate's leverage point and may be the wrong reference
   state.** G93A comes from Lindberg 2005 (monomer < dimer, like 24 of 33 controls);
   G93S/R/D/V come from Stathopulos 2006 and are the *only* controls where
   monomer > dimer. So the "G93A/G93S/G93V same-site ladder … internal consistency check"
   in `config/pipeline.yaml` crosses two datasets and cannot serve that purpose. With n=8,
   G93V alone sets the top of the −0.2 … 7.0 range and will dominate the gate Pearson.
   Also worth sanity-checking against apo-monomer total unfolding ΔG (~4–5 kcal/mol).
3. **`data/structures/A4V.pdb` on disk is OXIDIZED** — CYS57/CYS146 carry no HG. It
   predates the disulfide fix and `data/structures/` is gitignored, so
   `assert_disulfide_reduced` never ran on it. Inert for FEP (the pmx engine builds its
   own input) but Stage 2 and Stage 4 would consume it. `rm data/structures/*.pdb` except
   3ECU forces a correct rebuild.
4. **The GROMACS leg never verifies the disulfide stayed reduced.** Stage 1 asserts;
   Stage 3 pipes `"n"` answers to `pdb2gmx -ss` open-loop and never checks the resulting
   topology. `_RESNAME_ALIASES` maps CYX→CYS, so a re-formed bond would pass the
   wild-type check too.
5. **The protocol fingerprint does not survive a resume.** It is computed from the
   *current* config at run time, then `mdrun -cpi ... -append` appends to the existing
   `dhdl.xvg`. A window killed under one protocol and resumed under another yields one
   file with two protocols and a single hash claiming the new one. Persist the hash to
   `run_dir/protocol.sha` on first write and refuse a mismatched resume.
6. **`build_capped_tripeptide` builds a collinear ACE cap** (`tripeptide.py:108-113`):
   CH3–C–N and C–N–CA both 180° against ideals of ~117°/~121°, with O placed along an
   arbitrary lab-frame axis, so the unfolded leg is not invariant under rotation of the
   input PDB. Largely cancels within a leg, but it is a real reproducibility defect.
7. **`classify_uncharacterized` ignores `ddg_err`** — a VUS at 1.1 ± 0.9 is labelled
   destabilizing with the same confidence as 5.0 ± 0.2. This is the one place the
   pipeline turns numbers into clinical-adjacent claims; it needs an abstain band.
8. **No negative controls.** The panel is 54 positive_control + 38 vus. README §4 lists
   benign controls as required stratification, and `classify_uncharacterized`'s "stable"
   label has no calibration set without them.
9. **Charge-neutrality of `gate_subset` is not enforced in code** — only a curated list
   and a comment. Two lines in `_validate_gate_subset` reading `charge_change` close it.
10. **Interface proximity is not screened for the gate subset.** I113T is a genuine
    interface residue simulated as a carved-out monomer; nothing flags that.
11. **`AGENTS.md` is a stale fork of CLAUDE.md** (71 lines vs ~200, last touched
    `64f6e5d`). It is missing the pmx stack section, the GROMACS invariants, lambda-aware
    minimisation, the no-pull-mid-array rule, and the charge-neutral gate constraint. Any
    agent routed there gets the pre-audit rulebook.
12. **`snakemake` is absent from `env/environment.yml`** though CLAUDE.md's Verification
    section mandates `snakemake -n` after any rule change.
13. **`rule prep`'s output is read by nothing.** `data/structures/{variant}.pdb` is
    declared an input to `fep_window` and `prescreen`; neither reads it, and it is a
    solvated *mutant* while the FEP path builds its own WT. Dead wiring gating 864 jobs.

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
