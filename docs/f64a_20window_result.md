# F64A under the redistributed ladder — protocol `822108e9db71124d` (2026-08-30)

The rerun that [`f64a_folded_leg_failure.md`](f64a_folded_leg_failure.md) called for.
120 windows (20 λ × 3 replicates × 2 legs), refined λ endpoint, and
`independent_replicate_systems: true`.

**Outcome: the convergence defects were real and the fixes measurably helped, but F64A
still fails the gate's convergence criterion.** Cycle closure 2.52 → 1.10 against a cap
of 1.0.

## Result

| | 18-window `f9bded6f07b4abe5` | 20-window `822108e9db71124d` |
|---|---|---|
| ΔΔG | 6.63 ± 0.19 | **6.94 ± 0.34** |
| cycle closure (cap 1.0) | 2.52 | **1.10** |
| `converged` | false | **false** |
| replicate spread | 0.64 | **1.02** |
| min adjacent overlap | 0.018 | 0.015 |
| independent samples | 54,299 / 324,108 | 60,617 / 360,120 |
| `solver_notes` | JAX banner only | none |

Experimental reference is −0.20, but see the caveat below before treating +7.1 as the error.

### Per leg

| leg | r0 | r1 | r2 | spread |
|---|---|---|---|---|
| folded | +1.25 (hyst 0.79, ovl 0.015) | **+2.51 (hyst 0.02, ovl 0.055)** | +1.36 (hyst 1.10, ovl 0.029) |
| unfolded | −5.36 (hyst 0.18) | −5.12 (hyst 0.13) | −5.25 (hyst 0.19) | 0.24 |

Per-replicate ΔΔG: 6.61, 7.62, 6.60.

## What the changes bought

**`independent_replicate_systems` did exactly what it was added for.** The replicate
spread *rose* 0.64 → 1.02 and `ddg_err` rose 0.19 → 0.34. That is the correct direction:
the old agreement was an artefact of three replicates sharing one box and being stuck
behind the same barrier. The new spread is honest. Note it is still larger than the quoted
`ddg_err`, so the error bar continues to understate the true uncertainty.

**The endpoint refinement roughly halved the hysteresis** (2.52 → 1.10) without fixing the
overlap (0.018 → 0.015). Those are different problems: the endpoint hole was costing
forward/reverse consistency; the thin mid-ladder overlap is undersampling and needs
sampling time, not spacing.

## The finding that matters most: low hysteresis does not mean the right basin

| folded replicate | ΔG | hysteresis |
|---|---|---|
| r0 | +1.25 | 0.79 |
| **r1** | **+2.51** | **0.02** |
| r2 | +1.36 | 1.10 |

r1 has **the lowest hysteresis recorded anywhere in this project** — and it is the
replicate that disagrees with the other two by ~1.2 kcal/mol. With independent boxes, each
replicate settled into its own basin and sampled that basin consistently. Cycle closure
measures *within-basin* self-consistency; it cannot see that the basins differ.

This is a single-variant demonstration of what G93A only suggested across variants (closure
0.08, still 1.05 kcal/mol wrong). It is the sharpest available statement of the C2 claim:
**apparent convergence and accuracy are separable, and the standard diagnostic cannot tell
them apart.** It also means `replicate_spread_kcal` is now the more informative field of the
two — but only since independent boxes made it meaningful.

## Additional evidence on the experimental value

Complements the F45A/F64A comparison in the failure doc. Within Nordlund & Oliveberg 2006
alone — one source, one method, one lab:

| variant | side-chain heavy atoms deleted | burial (neighbours <6 Å) | exp ΔΔG |
|---|---|---|---|
| I18V | 1 | 66 | 0.37 |
| I149A | 3 | 72 | 4.05 |
| **F64A** | **6** | **70** | **−0.20** |

F64A deletes twice what I149A does from an equally buried site and is reported at roughly
zero. Extrapolating the other two gives ~1.35 kcal/mol per deleted atom, i.e. ~8 for F64A —
so the FEP value sits on that trend and the experimental value is the outlier.

**This is not grounds to drop F64A**, for the reason already argued in the failure doc:
excluding a control because it produced an inconvenient FEP number is the same category of
post-hoc adjustment as lowering `min_pearson`. It is recorded so that the +7 "error" is not
uncritically attributed to the calculation. The primary values must be checked against
Kumar 2017 Table S1 and Nordlund & Oliveberg 2006 before F64A anchors any gate.

## Status

**Not gate data.** `converged: false` (closure 1.10 > 1.0), so `evaluate_gate` excludes it
on a pre-registered criterion, independent of any question about the reference value.

Consequence: F64A is out, so **6 of the remaining 7 gate variants must converge** to reach
`min_gate_points: 6`. One further failure ends the gate as currently specified.

Raw data: `results/fep/F64A` (120 windows). The superseded 18-window run is archived at
`results/archive/F64A_f9bded6f07b4abe5`.
