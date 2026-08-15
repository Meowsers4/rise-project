# F64A under the post-fix protocol — folded-leg sampling failure (2026-08-15)

First variant run under the **post-fix** protocol (`sc-coul = yes`, 3001 dH/dλ records per
window, corrected ACE/NME caps, `solvent_padding_nm` 1.2), protocol fingerprint
`f9bded6f07b4abe5`, 18 uniform λ-windows, 3 replicates, one shared solvated box per
(variant, leg). Recorded here because `results/` is gitignored and this run is being
archived and rerun — and because it is the cleanest evidence in the project so far that
the *folded leg* of a large deletion is where apo-2SH sampling breaks.

## Result

| | value | |
|---|---|---|
| FEP ΔΔG | **+6.63 ± 0.19** kcal/mol | experiment (apo monomer) **−0.20** → error **+6.83** |
| cycle closure | **2.52** | cap 1.0 → `converged: false` |
| min adjacent overlap | **0.018** | pre-fix runs were 0.040–0.056, already called marginal |
| per-replicate ΔΔG | 6.74 / 6.89 / 6.26 | spread 0.64 against a quoted error of 0.19 |
| independent samples | 54,299 of 324,108 | 108 × 3001 exactly — no window double-appended |

This is by far the largest error recorded: 2.8× A4V's pre-fix +2.40, and on a variant
whose experimental value is ≈0. F64A is one of the near-zero controls whose job is to
catch a method that calls everything destabilizing.

## The failure is entirely in the folded leg

| leg | rep | dG (kcal/mol) | hysteresis | min overlap | fewest samples |
|---|---|---|---|---|---|
| folded | r0 | +1.62 | 1.23 | 0.018 | 27 @ w13 |
| folded | r1 | +1.76 | 1.00 | 0.033 | 14 @ w17 |
| folded | r2 | +1.09 | **2.52** | 0.030 | 16 @ w17 |
| unfolded | r0 | −5.12 | 0.14 | 0.115 | 323 @ w17 |
| unfolded | r1 | −5.14 | 0.26 | 0.115 | 183 @ w17 |
| unfolded | r2 | −5.17 | 0.23 | 0.072 | 227 @ w2 |

The unfolded leg is essentially perfect — three replicates agreeing to **0.05 kcal/mol**,
hysteresis well inside the cap, overlap 4–6× the folded leg's. The capped-tripeptide
reference state is **not** implicated. Every folded replicate is at or over the 1.0
hysteresis cap, and all three collapse at the top of the ladder.

Folded r0 independent samples per window (3001 raw records each, 1 ps apart):

```
w0..w8    735  269  626  473  744  689  694  684  408
w9..w17   554  495  143  398   27  319  207  271   45
```

Correlation time runs **~4 ps at low λ and ~111 ps at w13** — a 25× slowdown that switches
on as the phenylalanine disappears. Windows 13 and 17 carry 27 and 45 independent samples.

## Mechanism

F→A deletes an entire benzyl group — seven heavy atoms, the largest alchemical change
attempted in this project (A4V adds two carbons, G93A one, I113T swaps). Removing a buried
aromatic leaves a cavity that water must fill, and solvent reorganising into a
newly-vacated pocket is slow against a 3 ns window. That predicts damage concentrated in
the folded leg at high λ, which is exactly what the table shows, and it predicts a positive
bias, which is what the +6.83 error is.

## What this establishes

1. **`ddg_err` is not an uncertainty on ΔΔG.** 0.19 against a cycle closure of 2.52, a 13×
   gap. All three replicates shared one solvated box and one minimised structure
   (`stable_seed(variant, leg, "system")` had no `rep` in it), so they were stuck behind
   the same solvent barrier and their spread measured nothing. This is the same conclusion
   as `prefix_diagnostics.md` §5, now with the mechanism localised.
2. **The convergence machinery does fire.** Unlike the pre-fix three — all of which
   reported `converged: true` while being wrong by 0.8–2.4 — F64A was correctly refused.
   The diagnostics are not useless; they were being read at the wrong threshold of damage.
3. **MBAR itself was clean.** `solver_notes` contained only pymbar's JAX import banner and
   no fallback warnings, so the solver converged. The failure is in the input samples, not
   the estimator. (The banner is filtered out as of this change; it made the field
   unreadable without ever hiding anything.)
4. **A large deletion is a different sampling problem from a substitution.** The protocol
   that produced acceptable diagnostics for A4V, G93A and I113T is not adequate for F64A,
   and F64A is in `gate_subset`.

## Consequence for the gate

F64A is the **only** sub-zero point in the 8-variant `gate_subset`:

```
F64A -0.20 | I18V +0.37 | I113T +1.25 | A4V +1.62
G93A +2.43 | G93S +3.70 | I149A +4.05 | G93V +7.00
```

It is the single anchor for the "does not destabilize" end of the gate. Without it the
subset spans +0.37 … +7.00 and the Pearson is dominated by G93V — already flagged in
`HANDOFF.md` as the gate's leverage point and a possibly wrong reference state (G93V comes
from Stathopulos 2006, the only control source where monomer > dimer). A gate that loses
F64A and leans on G93V is measuring one point twice.

## A separate concern: the experimental reference value

Raised here because it was found while diagnosing this run, and it must be on record
*before* the gate is evaluated rather than discovered afterwards.

F64 is genuinely buried — burial rank 45/153 by heavy-atom neighbour count within 8 Å of
the side-chain centroid in `3ECU` chain A, comparable to I149 (99 neighbours) and L106
(102). It is not a surface phenylalanine whose truncation would be expected to cost
nothing. Its closest structural analogue in the same control set disagrees sharply:

| variant | burial | mutation | heavy atoms deleted | exp ΔΔG |
|---|---|---|---|---|
| F45A | 93 | Phe→Ala | 6 | **+2.07** |
| F64A | 94 | Phe→Ala | 6 | **−0.20** |

Two near-identical mutations, 2.3 kcal/mol apart. The compilation has more scatter of this
kind (V81A +0.01 against V29A +2.81, both Val→Ala at comparable burial), so the ±0.3 error
bars carried in `variants.csv` understate what it can support. F64A's −0.20 is the least
corroborated point in `gate_subset`, and it is that subset's **only** sub-zero anchor.

**This is not grounds to remove F64A from the gate.** Dropping a control because it
produced a bad FEP number is the same category of post-hoc adjustment as lowering
`min_pearson`, and it would destroy what pre-registration buys. If F64A is ever excluded it
must be on independently established grounds, argued and written down before a gate
evaluation, with the Kumar 2017 / Nordlund & Oliveberg 2006 primary values checked against
the source. Recorded here as an open question, not a decision.

## Status

**Not a result.** `converged: false`, so `evaluate_gate` excludes it. Superseded by the
rerun under the redistributed-ladder protocol; kept as the justification for that change.
