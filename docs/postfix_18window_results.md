# Post-fix, 18-window results — protocol `f9bded6f07b4abe5`

A4V and G93A, run after the 2026-08-09 protocol fixes (`sc-coul = yes`, `nstdhdl` 20 ps →
1 ps, corrected ACE/NME caps, box padding 1.2 nm) but BEFORE the 2026-08-15 endpoint and
independent-box changes. Recorded here because `results/` is gitignored and both are
invalidated by protocol `822108e9db71124d` — they are 18-window runs and must be rerun.

They are kept for one reason: **they are the only before/after pair that isolates what the
2026-08-09 fixes actually bought.** Pre-fix numbers are in
[`prefix_diagnostics.md`](prefix_diagnostics.md).

## Result

| variant | pre-fix | post-fix | exp | pre-error | post-error |
|---|---|---|---|---|---|
| A4V | 4.02 ± 0.23 | **2.77 ± 0.30** | 1.62 | +2.40 | **+1.15** |
| G93A | 1.21 ± 0.06 | **1.38 ± 0.11** | 2.43 | −1.22 | **−1.05** |

A4V's error halved. G93A's moved 14%, in the right direction (it was under-predicting and
went up). Both `converged: true`. RMSE over the pair 1.10, inside the 1.5 bound; Pearson is
meaningless on two points 0.81 kcal/mol apart against ~1.1 errors.

| | A4V | G93A |
|---|---|---|
| cycle closure (cap 1.0) | 0.49 | 0.08 |
| min adjacent overlap | 0.031 | **0.017** |
| replicate spread | 1.02 | 0.35 |
| quoted `ddg_err` | 0.30 | 0.11 |
| per-replicate ΔΔG | 3.24, 2.22, 2.86 | 1.16, 1.51, 1.45 |
| independent samples | 52,854 / 324,108 | 66,247 / 324,108 |
| `solver_notes` | none | none (JAX banner only) |

## Per-leg detail

| variant | leg | r0 | r1 | r2 | spread | hysteresis | overlap |
|---|---|---|---|---|---|---|---|
| A4V | folded | −13.15 | −13.93 | −13.46 | 0.78 | 0.49 / 0.29 / 0.20 | 0.031–0.064 |
| A4V | unfolded | −16.39 | −16.15 | −16.32 | 0.24 | 0.06 / 0.01 / 0.12 | 0.092–0.105 |
| G93A | folded | +9.75 | +9.90 | +9.88 | 0.15 | 0.08 / 0.02 / 0.04 | **0.017**–0.057 |
| G93A | unfolded | +8.59 | +8.39 | +8.42 | 0.20 | 0.08 / 0.02 / 0.00 | 0.044–0.076 |

## What they establish

1. **The 2026-08-09 fixes were worth roughly 1.25 kcal/mol on A4V and 0.17 on G93A.** Not
   a uniform shift — A4V came down, G93A went up, each toward experiment. That rules out a
   constant calibration offset as the explanation for either.
2. **`converged: true` still does not mean accurate.** G93A reports the cleanest closure in
   the whole project (0.08) and is 1.05 kcal/mol wrong. F64A was caught only because it
   failed on a different axis (closure 2.52).
3. **Cycle closure cannot see a thin ladder.** G93A's minimum adjacent overlap is 0.017 —
   effectively the same as F64A's 0.018 — yet it passes every convergence criterion. There
   is no overlap floor among the gate criteria; adding one remains open in HANDOFF §2.
4. **`ddg_err` understates the error 3–10×** in both. Replicate spread exceeds the quoted
   uncertainty because all three replicates shared one box — the defect
   `independent_replicate_systems` addresses.
5. The folded/unfolded asymmetry holds: every unfolded leg is tight and low-hysteresis;
   every problem is in the folded leg.

## Status

Superseded by `822108e9db71124d` (20 windows, refined λ endpoint, independent replicate
boxes). Not gate data — `evaluate_gate` now checks `protocol` across variants and would
refuse to mix these with post-2026-08-15 runs.
