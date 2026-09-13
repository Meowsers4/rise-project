# Stop-rule reanalysis: every documented estimate re-derived from raw windows — 2026-09-12

Week-1 step 2 of [`post_pivot_review.md`](post_pivot_review.md) §5, and the **stop rule** the
four-week recommendation is conditioned on:

> if provenance-consistent reanalysis does not reproduce the documented F64A low-hysteresis /
> replicate-disagreement example, withdraw this recommendation.

**Verdict: the stop rule passes. The four-week plan stands.** Every documented ΔΔG, cycle
closure, per-leg ΔG, hysteresis and overlap **except A4V's gate point** re-derives from the
archived NPZ windows to
floating-point rounding, in **three independent local environments spanning both pymbar
backends** (§1): Δ ≤ 3.1e-12 against the SCC records across all eight variants (≤ 3.2e-12 for
any pairwise comparison), and ≤ 5.7e-15 on every F64A field (3.1e-15 on the per-leg estimates
tabulated in §3). **A4V's gate point is the one exception and cannot be re-derived** — its 3 ns
folded windows were overwritten by the attempt-2 resubmission, so it survives only as a
run-time JSON record (§2, §7). It contributes to the gate r/RMSE, so that statistic is
seven-eighths reproduced and one-eighth attested. The review's finding that "every reported
SOD1 result in this repo is documentation arithmetic, not a reproduced MBAR estimate" is
discharged for every result but that one.

**One documented claim is false and is corrected below** (§5): F64A folded r1's hysteresis is
*not* the lowest recorded anywhere in the project. The finding it was cited for survives, and
is now quantified across 48 records instead of one (§6).

## 1. What was run

| | |
|---|---|
| Input | `~/sod1fep_archive_2026-09-11/fep/` — 8 variants × 2 legs × 20 λ × 3 replicates = 960 windows |
| Code | `src.fep.analyze.analyze_variant` at `024114b`, unmodified |
| Config | `config/pipeline.yaml` at `024114b`, unmodified |
| Machine (a) | local Mac, throwaway venv on Python **3.12.13**, numpy 2.5.3, scipy 1.18.1, pymbar 4.0.3 |
| Machine (b) | same Mac, conda env `rise`: Python **3.13.9**, numpy 2.5.0, scipy 1.18.0, pymbar 4.0.3, no JAX. Agrees with (a) to ≤ 1.3e-12 and with the SCC records to ≤ 1.1e-12. |
| Machine (c) | same Mac, venv identical to (a) but **with JAX 0.11.1** (pymbar's 64-bit JAX backend, as the SCC used). Isolates the backend: (a) and (c) differ in nothing else. Agrees with the SCC records to ≤ 3.1e-12. |
| pymbar backend | Both backends exercised. (a) and (b) are numpy/scipy; (c) is JAX, matching the SCC. See §2.1 — matching the backend did **not** reduce the residual. |
| Compare against | the run-time `results/convergence/<V>.json` written by the original SCC runs (Python 3.11, JAX backend), pulled 2026-09-12 |
| Output | `~/sod1fep_archive_2026-09-11/reanalysis_2026-09-12/` — (a) top level, (b) `rise_py313/`, (c) `jax_py312/` |
| Test suite | 105 passed, 12 skipped (skips need OpenMM) — the documented baseline |

### Why the current config is provenance-consistent with runs from 2026-08-30

`analyze_variant` reads six config values — `lambda_windows`, `replicates`, `legs`,
`temperature_K`, `convergence.max_cycle_closure_kcal` and `decorrelate`. (`fep.framework` is
checked too, but by [`main()`](../src/fep/analyze.py#L433) and `validate.py`, not by
`analyze_variant`.) It does **not** read `ns_per_window`, `equilibration_ns` or
`lambda_vector` — those shape the windows, and the windows already exist.

All six are unchanged since before the F64A run; `legs:` has not changed since the initial
scaffold `a043270`, and `max_cycle_closure_kcal` has been 1.0 since the same commit. Note that
`git log -S legs` flags `927b302` and `7db3095`, but **both are false positives** — each merely
added a comment containing the word "legs" (`legs*windows*replicates`, "both legs"); the `legs:`
key itself is untouched, as `git log -L '/^  legs:/,+5'` confirms.

The only change to `analyze.py` since F64A ran (`7db3095`) converts the protocol check from
per-variant to per-leg. It touches provenance bookkeeping and no numerics. So the code path
that produced today's numbers is the same one that produced the recorded ones.

The archived F64A windows carry `shape (20, 3001)`, protocol `822108e9db71124d`, provenance
`gromacs_pmx`, and `lambda_index == w` on all 120 — matching
[`raw_result_reconciliation.md`](raw_result_reconciliation.md) exactly.

## 2. Reproduction, all eight variants

| variant | recorded ΔΔG | re-derived ΔΔG | recorded closure | re-derived closure | source |
|---|---|---|---|---|---|
| I18V | 0.83 ± 0.21 | **0.8292 ± 0.2074** | 0.18 | **0.1810** | raw NPZ |
| I113T | 2.59 ± 0.23 | **2.5908 ± 0.2337** | 0.22 | **0.2249** | raw NPZ |
| G93A | 1.26 ± 0.06 | **1.2587 ± 0.0612** | 0.04 | **0.0450** | raw NPZ |
| G93S | 1.20 ± 0.07 | **1.2036 ± 0.0652** | 0.12 | **0.1162** | raw NPZ |
| I149A | 4.98 ± 0.16 | **4.9819 ± 0.1567** | 0.77 | **0.7678** | raw NPZ |
| G93V | 2.83 ± 0.08 | **2.8309 ± 0.0833** | 0.67 | **0.6696** | raw NPZ |
| F64A | 6.94 ± 0.34 | **6.9443 ± 0.3388** | 1.10 | **1.1038** | raw NPZ |
| A4V (gate) | 3.54 ± 0.39 | 3.5441 ± 0.3948 | 0.66 | 0.6567 | **run-time record only** |
| A4V (attempt 2) | 3.36 ± 0.08 | **3.3591 ± 0.0827** | — | **0.5054** | raw NPZ |

Seven of the eight re-derive from raw output. **A4V's gate point does not, and cannot**: the
attempt-2 resubmission overwrote its 3 ns folded windows in place. What survives is
`archive/A4V_3ns_ddg.json` + `A4V_3ns_convergence.json`, the run-time JSON the original analysis
wrote — full per-leg detail (folded −13.412/−12.964/−11.952, unfolded −16.365/−16.350/−16.245),
but not re-derivable. That is the cost recorded in `raw_result_reconciliation.md` §"What was
lost", and it is now the single weakest link in the chain of custody.

### 2.1 The residual is solver noise, not a backend difference

Running (c) with JAX — the same backend the SCC used, in an environment otherwise identical to
(a) — was expected to tighten the agreement. It did the opposite:

| comparison | worst Δ | where |
|---|---|---|
| (c) JAX vs SCC — **same backend** | **3.034e-12** | `I18V per_replicate_ddg[2]` |
| (a) non-JAX vs SCC | 1.158e-12 | `I18V per_replicate_ddg[2]` |
| (b) non-JAX vs SCC | 1.024e-12 | `I113T` overlap element |
| (c) vs (a) — **backend isolated** | 1.876e-12 | `I18V per_replicate_ddg[2]` |

So the ~1e-12 residual is **not** attributable to the backend: matching it exactly does not
remove the disagreement, and the same field dominates every comparison. It is accumulated
floating-point noise in the MBAR solve and the subsequent summation. `I18V per_replicate_ddg[2]`
is the worst field in five of six comparisons, which is what one would expect of the
smallest-ΔΔG variant in the set (0.83 kcal/mol) rather than of a systematic backend effect.

This is a decomposition, not a strengthening: the headline bound is set by the noisiest
comparison and therefore *rose* from 1.2e-12 to 3.1e-12 when (c) was added. Recorded that way
deliberately — the looser number is the honest one, and 3e-12 kcal/mol remains ~12 orders of
magnitude below anything the project reports.

One incidental check: pymbar's JAX 64-bit banner is emitted on stderr and
[`_solver_notes`](../src/fep/analyze.py#L157) is supposed to strip it. Under (c) the F64A
`solver_notes` came back `[]`, matching the SCC record exactly — the first time that stripping
has been exercised against genuine JAX output rather than against the SCC's saved result.

A4V under 9 ns folded also carries two protocol hashes by design —
`folded=cf1e632168579261|unfolded=822108e9db71124d` — and `_check_single_protocol` accepted it
per leg, as intended.

## 3. F64A per-leg, recorded vs re-derived

The stop rule's actual subject. Recorded 2026-08-30 on the SCC; re-derived 2026-09-12 locally.

| record | ΔG recorded | ΔG re-derived | hyst recorded | hyst re-derived | min ovl | samples |
|---|---|---|---|---|---|---|
| folded r0 | +1.2496 | **+1.2496** | 0.7866 | **0.7866** | 0.0146 | 9405 |
| folded **r1** | **+2.5059** | **+2.5059** | **0.0227** | **0.0227** | 0.0554 | 7735 |
| folded r2 | +1.3572 | **+1.3572** | 1.1038 | **1.1038** | 0.0289 | 8320 |
| unfolded r0 | −5.3567 | **−5.3567** | 0.1794 | **0.1794** | 0.0745 | 11070 |
| unfolded r1 | −5.1159 | **−5.1159** | 0.1258 | **0.1258** | 0.0881 | 11655 |
| unfolded r2 | −5.2476 | **−5.2476** | 0.1918 | **0.1918** | 0.0291 | 12432 |

Largest disagreement on any field in this table: **3.1e-15**, identical under both local environments; over F64A's top-level scalars it is 1.8e-15 (`replicate_spread_kcal`), and 5.4e-15 across the full 20×20 overlap matrices. Independent-sample counts are identical
(60,617 of 360,120), so `pymbar.timeseries` made the same decorrelation decisions on both
machines.

**The stop-rule example reproduces:**

- folded r1 hysteresis **0.0227** — the lowest of F64A's three folded replicates and of all six F64A records
- folded r1 ΔG **+2.5059** vs siblings +1.2496 / +1.3572 (mean +1.3034)
- disagreement **+1.2026 kcal/mol**

All three figures are bit-identical under all three local environments (§1 (a), (b), (c)), so
the stop-rule example depends on neither the numpy/scipy build nor the pymbar backend.

## 4. Gate arithmetic, recomputed from re-derived values

| | recorded | recomputed |
|---|---|---|
| Pearson r | 0.326 | **0.32601** |
| Spearman ρ | 0.500 | **0.50000** |
| RMSE | 2.123 | **2.12349** |
| SSE | 31.5663 (review) | **31.5646** |
| G93V share of SSE | 17.3889 (review) | **17.3814** |
| other six | 14.1774 (review) | **14.1832** |

The gate verdict is unchanged: FAIL on `pearson` and `rmse`, pivot line triggered. The small
SSE difference is expected — the review computed from the rounded table in
[`gate_attempt_1.md`](gate_attempt_1.md); these come from full-precision re-derived values. It
does not move the corrected G93V requirement (≤ 1.254) by a meaningful amount.

## 5. Correction: "the lowest hysteresis recorded anywhere" is false

Three files state that F64A folded r1's 0.02 is the lowest hysteresis in the project:
`CLAUDE.md:164`, `HANDOFF.md:214`, `docs/f64a_20window_result.md:55`. Ranking all 48
(variant, leg, replicate) records:

| rank | hysteresis | record |
|---|---|---|
| 1 | 0.0017 | G93A / unfolded / r2 |
| 2 | 0.0034 | A4V / unfolded / r1 |
| 3 | 0.0067 | G93S / **folded** / r0 |
| 4 | 0.0098 | G93S / **folded** / r1 |
| 5 | 0.0147 | G93S / unfolded / r0 |
| … | | |
| **10** | **0.0227** | **F64A / folded / r1** |

F64A folded r1 ranks **10th of 48**, and 3rd of 24 among folded legs. G93A unfolded r2 is 13×
lower. The claim was never checked against the full set; it was true within F64A and got
generalised. All three files are corrected in the same commit as this document.

This is a correction of the same kind as items 3, 6 and 7 in `post_pivot_review.md` §7 — a
factual error in committed documentation, found by measurement — not a scope decision.

## 6. What replaces it, and it is stronger

The superlative was decoration. The finding it decorated is intact, and the archive supports a
quantitative version the single anecdote could not.

**F64A folded r1 has the largest replicate disagreement of all 48 records (rank 1/48,
1.2026 kcal/mol) while sitting in the lowest ~20% of the hysteresis distribution (rank
10/48).** That is the same point, stated as a position in a measured distribution rather than
as a superlative.

Across all 48 records, taking each replicate's |ΔG − mean(its two independently-solvated
siblings)| as the disagreement:

| set | n | Pearson | Spearman |
|---|---|---|---|
| all records | 48 | +0.218 | +0.250 (p = 0.087) |
| **folded leg only** | 24 | **+0.072** | **+0.080** |
| unfolded leg only | 24 | +0.383 | +0.357 |

**On the folded leg — the one that limits every result in this project — within-ladder
hysteresis carries essentially no information about whether a replicate agrees with
independently solvated repeats.** The overall correlation is weak and not significant at
n = 48; what little there is comes from the unfolded leg, where both quantities are small and
the tripeptide converges anyway.

Both tails are informative. G93S folded r0 and r1 have the 3rd and 4th lowest hysteresis in
the dataset *and* agree with their siblings (0.11, 0.17) — low hysteresis is not a *predictor*
of disagreement either. The diagnostic is close to uninformative on this leg in both
directions, which is the honest claim.

Folded-leg replicate spread, for context — F64A is the outlier variant, not merely the outlier
replicate:

| variant | folded spread | unfolded spread |
|---|---|---|
| **F64A** | **1.256** | 0.241 |
| I113T | 0.583 | 0.623 |
| I18V | 0.556 | 0.231 |
| I149A | 0.401 | 0.204 |
| A4V (9 ns) | 0.266 | 0.121 |
| G93A | 0.250 | 0.048 |
| G93S | 0.189 | 0.166 |
| G93V | 0.088 | 0.324 |

### The same is true of minimum adjacent overlap

`min_adjacent_overlap` is the other diagnostic `CLAUDE.md` tells us to read. Re-derived per
variant, it does not order the failures either:

| variant | closure | min adj overlap | converged | \|error\| vs exp |
|---|---|---|---|---|
| I18V | 0.181 | **0.0054** | yes | 0.46 |
| I149A | 0.768 | **0.0096** | yes | 0.93 |
| F64A | 1.104 | 0.0146 | **no** | — (excluded) |
| G93A | 0.045 | 0.0249 | yes | 1.17 |
| I113T | 0.225 | 0.0288 | yes | 1.34 |
| G93V | 0.670 | 0.0342 | yes | 4.17 |
| A4V (9 ns) | 0.505 | 0.0384 | yes | 1.74 |
| G93S | 0.116 | 0.0425 | yes | 2.50 |

The two thinnest ladders in the project (I18V 0.0054, I149A 0.0096) belong to the two most
*accurate* variants. Ranking the seven usable variants:

| pair | n | Spearman |
|---|---|---|
| min adjacent overlap vs \|error\| | 7 | **+0.893** (p = 0.007) |
| cycle closure vs \|error\| | 7 | +0.107 (p = 0.819) |

Closure carries no signal, as expected from §6. Overlap correlates **strongly and in the
wrong direction** — wider overlap goes with *larger* error.

**Do not read that as a result.** At n = 7 it is one rearrangement away from nothing, and it
is confounded with mutation site: the position-93 variants carry the project's dominant
systematic error (the 3× compression in `gate_attempt_1.md`) and sit at overlap ranks 1, 3
and 5 of 7. Excluding position 93 leaves n = 4 with a perfect Spearman of +1.000 — which is
what four points do, not evidence. The defensible statement is the negative one: **neither
diagnostic orders these seven variants by accuracy in the intended direction**, and the
dataset is too small and too confounded to say why. Recorded so the pattern is not
rediscovered later and mistaken for a finding.

This is consistent with the endpoint-refinement result in
[`f64a_20window_result.md`](f64a_20window_result.md) — halved hysteresis, unchanged overlap —
which already indicated the dominant error is not ladder density.

This retires the framing in the old `CLAUDE.md` bullet, which contrasted G93A's overlap with
F64A's as though those two were the extremes. On the 20-window data they are mid-pack.

## 7. What this does and does not establish

Established:

- Every documented estimate except A4V's gate point re-derives from raw windows on an independent machine, Python version, numpy version and pymbar backend.
- The decorrelation and MBAR path is deterministic across those changes to floating-point rounding.
- The gate failure (r = 0.326, RMSE 2.123) is a reproduced result, not documentation arithmetic.
- The F64A stop-rule example is real as stated, minus the superlative.
- On the folded leg, hysteresis and independent-box disagreement are near-uncorrelated over 48 records.

Not established, and not addressable by reanalysis:

- That the windows encode the intended physical system. Reproducing an estimate from saved reduced potentials cannot detect a wrong topology, redox state, residue mapping, cap geometry or unit. That is the reference-state audit, still to be written.
- That the experimental references are comparable to what was simulated (`docs/reference_state_audit.md`, not yet written).
- Anything about A4V's gate point beyond its run-time record.
- Generality. 48 records is 8 variants at **6 sites** (4, 18, 64, 93, 113, 149; the seven *usable* gate points span 5, F64A being the excluded one), one target, one force field, one estimator, one water model. The correlations in §6 are a statement about this dataset, not about FEP.
- Which replicate sampled the correct basin. Two agreeing replicates are still not ground truth — the review's caution stands.
- The §6 overlap/closure-vs-error rankings mix protocols: A4V's row is the 9 ns folded rerun, the other six are 3 ns. Six of seven are internally consistent; the seventh is the only A4V point that re-derives from raw output. Stated as a negative result, so the mixing weakens nothing it claims — but it rules the table out as a positive finding on its own.

## 8. Next

1. `docs/reference_state_audit.md` — the remaining week-1 deliverable, and now the binding one.
2. LiveCoMS presubmission inquiry (user's task).
3. The §6 analysis is the LiveCoMS note's central figure. It is also the smallest working instance of candidate 1 in `post_pivot_review.md` §6 — the same computation over public OpenFE repeats is that project's week-1 killer test.
