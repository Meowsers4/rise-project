# Raw-result inventory and archive — 2026-09-11

Week-1 step 1 of the plan in [`post_pivot_review.md`](post_pivot_review.md): establish that
the raw windows behind every documented result still exist, are internally consistent, and
are preserved outside the live run directories.

**This records structural integrity only.** It does not verify a single ΔΔG.

> **Superseded in part, 2026-09-12.** That re-derivation has since been done — see
> [`stop_rule_reanalysis.md`](stop_rule_reanalysis.md). Every documented estimate except
> A4V's gate point now reproduces from these windows to floating-point rounding, so the
> review's "documentation arithmetic, not a reproduced MBAR estimate" finding no longer
> stands. The structural manifest below is unchanged and still current.

## Archive

| | |
|---|---|
| Cluster | `/projectnb/rise-batteries/bode/archive_2026-09-11/` — 7 variants, ~6.6 GB |
| Off-cluster | `~/sod1fep_archive_2026-09-11/fep/` on the local Mac |
| Verified by | `rsync -an --itemize-changes --exclude A4V results/fep/ <archive>/` → empty output (byte-for-byte complete) |
| A4V, 2026-09-12 | Landed and pulled to the **local Mac** — `fep/A4V` now 60/60 on both legs, plus `convergence/` (8 files) and `archive/` (the 3 ns A4V records + `F64A_f9bded6f07b4abe5`). |
| Still to do on the cluster | Copy A4V into `/projectnb/rise-batteries/bode/archive_2026-09-11/`. The off-cluster copy is currently the only archived one, so A4V is back to single-copy on the SCC side. |

`results/` is gitignored, so before today every raw window in the project existed in exactly
one place. A4V demonstrated the cost of that: see "What was lost" below.

## Manifest

Windows per leg, `u_kn_window` shape, mdp fingerprint, engine. Generated on the SCC by
walking `results/fep/<variant>/<leg>/w*_r*.npz`.

| variant | leg | n | shape | protocol | provenance |
|---|---|---|---|---|---|
| A4V | folded | ~~52~~ **60/60** (2026-09-12) | (20, 3000) | `cf1e632168579261` | gromacs_pmx |
| A4V | unfolded | 60 | (20, 3001) | `822108e9db71124d` | gromacs_pmx |
| F64A | folded | 60 | (20, 3001) | `822108e9db71124d` | gromacs_pmx |
| F64A | unfolded | 60 | (20, 3001) | `822108e9db71124d` | gromacs_pmx |
| G93A | folded | 60 | (20, 3001) | `822108e9db71124d` | gromacs_pmx |
| G93A | unfolded | 60 | (20, 3001) | `822108e9db71124d` | gromacs_pmx |
| G93S | folded | 60 | (20, 3001) | `822108e9db71124d` | gromacs_pmx |
| G93S | unfolded | 60 | (20, 3001) | `822108e9db71124d` | gromacs_pmx |
| G93V | folded | 60 | (20, 3001) | `822108e9db71124d` | gromacs_pmx |
| G93V | unfolded | 60 | (20, 3001) | `822108e9db71124d` | gromacs_pmx |
| I113T | folded | 60 | (20, 3001) | `822108e9db71124d` | gromacs_pmx |
| I113T | unfolded | 60 | (20, 3001) | `822108e9db71124d` | gromacs_pmx |
| I149A | folded | 60 | (20, 3001) | `822108e9db71124d` | gromacs_pmx |
| I149A | unfolded | 60 | (20, 3001) | `822108e9db71124d` | gromacs_pmx |
| I18V | folded | 60 | (20, 3001) | `822108e9db71124d` | gromacs_pmx |
| I18V | unfolded | 60 | (20, 3001) | `822108e9db71124d` | gromacs_pmx |

All seven gate-attempt-1 variants are complete at 60 windows per leg (20 λ × 3 replicates),
carry one fingerprint across both legs, and name one engine. This matches the protocol
stated in [`gate_attempt_1.md`](gate_attempt_1.md) and predates the per-leg sampling split
(`7db3095`, `5beab16`), which is why both legs share a hash.

**F64A is complete.** It was excluded from the gate by the pre-registered convergence
criterion (closure 1.10 > cap 1.0), not for missing data, so the replicate-disagreement
example the LiveCoMS note rests on can be re-derived from raw output.

### Two shapes and two hashes are correct now

A4V is running under the attempt-2 protocol and legitimately differs from the archive:

| leg | ns / equil | nstdhdl | nsteps | records | discarded | kept |
|---|---|---|---|---|---|---|
| folded | 9.0 / 2.0 | 1500 | 5,500,000 | 3667 | 667 | **3000** |
| unfolded | 3.0 / 0.5 | 500 | 1,750,000 | 3501 | 500 | **3001** |

Derived from `fep.frames_per_window: 3000` and the scheduling arithmetic in
`pmx_engine.py`; confirmed against the produced windows.
[`analyze._check_single_protocol`](../src/fep/analyze.py#L279) is **per leg** by design, so
a variant reporting `folded=cf1e63…|unfolded=822108…` is consistent, not contaminated. The
check in `HANDOFF.md` §3 predates this and still says to expect one shape and one hash.

## Run chronology

Preserved by `cp -a` from the source directory mtimes. Recoverable from nothing else.

| variant | last written |
|---|---|
| F64A | 2026-08-30 |
| G93A | 2026-09-01 |
| I149A, I18V | 2026-09-02 |
| G93S, I113T | 2026-09-03 |
| G93V | 2026-09-04 |

G93V last, immediately before the gate was evaluated and recorded on 2026-09-04.

## What was lost

A4V's folded leg reports only `cf1e632168579261`, so the 3 ns folded windows behind its
**3.54 ± 0.39** entry in `gate_attempt_1.md` were overwritten in place by the attempt-2
resubmission. That point can no longer be reconciled against raw output; the rounded value
in the gate table is all that survives. The unfolded leg is untouched.

Consequence for the plan: **archive before resubmitting, not after.** Weeks 2–3 rewrite
F64A and G93V, and both are already preserved above.

## What this does and does not establish

Established:

- Every documented gate-attempt-1 result has its complete raw windows on disk, in two places.
- No variant mixes protocols or engines; none carries an unlabelled window.
- Window counts match `legs × windows × replicates` exactly, so no silently-missing tasks.

Not established — and not addressable by this inventory:

- That any reported ΔΔG, cycle closure, or overlap value re-derives from these windows.
- That the windows encode the intended physical system (topology, redox state, numbering, caps, units).
- That the experimental reference values are comparable to what was simulated — see the separate reference-state audit.

## Next

1. ~~Re-run `src.fep.analyze` per variant against these archived windows~~ — **done 2026-09-12**, everything reproduces to floating-point rounding: [`stop_rule_reanalysis.md`](stop_rule_reanalysis.md) §2, §4.
2. ~~Recalculate F64A's per-leg estimates, per-replicate hysteresis and replicate SEM~~ — **done; the stop rule passes**: [`stop_rule_reanalysis.md`](stop_rule_reanalysis.md) §3.
3. ~~Add A4V to both archives~~ — **done 2026-09-12.** Its gate-attempt-1 point survives only as a run-time JSON record, as predicted above.
4. `docs/reference_state_audit.md` — the remaining week-1 deliverable. Reanalysis cannot reach it: re-deriving an estimate from saved reduced potentials cannot detect a wrong topology, redox state, mapping or reference construct.
