# Validation gate, attempt 1 — FAILED (2026-09-04)

Protocol `822108e9db71124d`: 20 λ-windows with a refined endpoint, 3 replicates with
independent solvated boxes, 3 ns production + 0.5 ns discarded equilibration per window,
`sc-coul = yes`, 1.2 nm box padding. 120 tasks per variant.

Thresholds pre-registered 2026-08-07, **before any gate evaluation**, and untouched.

```
VALIDATION GATE FAILED (failed on: pearson, rmse; pearson=0.326, rmse=2.123, n=7)
```

| criterion | value | threshold | |
|---|---|---|---|
| Pearson r | **0.326** | ≥ 0.70 | FAIL |
| Spearman ρ | 0.500 | reported | |
| RMSE (kcal/mol) | **2.123** | ≤ 1.5 | FAIL |
| median \|cycle closure\| | 0.220 | ≤ 0.75 | PASS |
| usable points | 7 | ≥ 6 | PASS |
| **pivot line** (README §10) | 0.326 < 0.60 | | **TRIGGERED** |

## The seven points

| variant | exp | FEP | error | closure |
|---|---|---|---|---|
| I18V | 0.37 | 0.83 ± 0.21 | **+0.46** | 0.18 |
| I113T | 1.25 | 2.59 ± 0.23 | **+1.34** | 0.22 |
| A4V | 1.62 | 3.54 ± 0.39 | **+1.92** | 0.66 |
| G93A | 2.43 | 1.26 ± 0.06 | **−1.17** | 0.04 |
| G93S | 3.70 | 1.20 ± 0.07 | **−2.50** | 0.12 |
| I149A | 4.05 | 4.98 ± 0.16 | **+0.93** | 0.77 |
| G93V | 7.00 | 2.83 ± 0.08 | **−4.17** | 0.67 |

F64A (exp −0.20) is excluded: `converged: false`, closure 1.10 against a cap of 1.0. That
exclusion is by a pre-registered criterion, not by its value —
[`f64a_20window_result.md`](f64a_20window_result.md).

## The failure has one address

Every non-glycine error is **positive**; all three position-93 errors are **negative**.

| | exp | FEP |
|---|---|---|
| G93A | 2.43 | 1.26 |
| G93S | 3.70 | 1.20 |
| G93V | 7.00 | 2.83 |

**Experiment spans 4.57 kcal/mol at that site. FEP spans 1.63** — a ~3× compression. The
calculation is not blind (it does rank G93V highest) but it separates G93A from G93S by
0.06 where experiment says 1.27.

Excluding position 93 — **diagnostic only, NOT a legitimate gate move** — the remaining
four give r = 0.935, ρ = **1.000** (perfect rank ordering), RMSE 1.280. Dropping controls
because they are inconvenient is the same category of post-hoc adjustment as lowering
`min_pearson`; it is computed here solely to localise the failure.

## Mechanism

Two opposite biases, both sampling:

- **Non-glycine (+):** the folded leg cannot reorganise around the new side chain within
  3 ns, so insertion looks too costly → over-destabilises. Magnitude tracks repacking
  demand.
- **Glycine (−):** Gly→X is dominated by loss of backbone φ/ψ freedom. Three nanoseconds
  does not sample the conformational space glycine explores, so the entropy penalty is
  underestimated → under-destabilises. Gly→X in constrained backbone geometry is the
  textbook hard case for alchemical FEP.

## Consequences

1. **No VUS may be classified.** `src.fep.window._assert_validation_gate` enforces this in
   code: a non-gate variant will not run while the report says `passed: false`.
2. **README §10 pivot is triggered** (r < 0.60): reframe as a methods/sampling-limits
   result rather than retune. No threshold is revisited.
3. **Gate attempt 2** raises folded-leg sampling 3 → 9 ns and equilibration 0.5 → 2.0 ns
   (`5beab16`), leaving the already-converged unfolded leg untouched. Legitimate because
   pre-registration constrains thresholds, not methods, and the change follows an
   independently diagnosed defect. **This is the last attempt justifiable on those
   grounds** — iterating protocols until one passes is fishing even with fixed thresholds.

### Arithmetic for attempt 2, stated in advance

RMSE ≤ 1.5 over seven points requires total squared error ≤ 15.75. It is currently **31.6,
of which G93V alone contributes 17.4**. G93V must fall from −4.17 to about −1.5 for the
gate to be arithmetically reachable — and it is the variant whose failure mode is least
likely to respond to folded-leg sampling. The expected outcome is improved non-glycine
errors, largely unchanged glycine errors, and a gate that still fails but isolates the
residual to backbone-entropy sampling at glycine sites.

Raw: `results/validation_gate.json`, `results/fep/<variant>/ddg.json`,
`results/convergence/<variant>.json`.
