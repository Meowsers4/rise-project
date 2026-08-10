# Pre-fix diagnostic runs (2026-08-07 → 2026-08-09)

Three variants run under the **pre-fix** Stage 3 protocol, before `sc-coul`, the 20×
sampling increase, the corrected ACE/NME caps, and `solvent_padding_nm` 1.2. Recorded
here because `results/` is gitignored and these numbers are the evidence behind the
claim that **apparent convergence does not imply accuracy** — which is the sharper form
of the limitation Wells 2021 reported for apo-2SH (README §11.1).

Protocol: GROMACS+pmx, 18 λ-windows, 3 replicates, 3 ns/window + 0.5 ns discarded,
`nstdhdl` 10000 (one dH/dλ record per 20 ps), `sc-coul` **off**, box padding 1.0 nm,
collinear ACE/NME caps. Provenance `gromacs_pmx` on every window.

## Result

| variant | FEP ΔΔG | exp (apo monomer) | error | cycle closure | min overlap |
|---|---|---|---|---|---|
| A4V | +4.02 ± 0.23 | 1.62 | **+2.40** | 0.63 | 0.040 |
| G93A | +1.21 ± 0.06 | 2.43 | **−1.22** | 0.09 | 0.043 |
| I113T | +2.08 ± 0.18 | 1.25 | **+0.83** | 0.34 | 0.056 |

RMSE 1.63 kcal/mol (pre-registered bound 1.5). Pearson −0.50, Spearman −0.50 — but the
three experimental values span only 1.18 kcal/mol while the errors are 0.8–2.4, so the
correlation is not measurable on this subset and should not be quoted. The gate subset
spans −0.20 … 7.00 precisely to avoid that. `min_gate_points` is 6, so no gate was
evaluated.

## Per-leg detail

| variant | leg | r0 | r1 | r2 | spread | hysteresis |
|---|---|---|---|---|---|---|
| A4V | folded | −12.56 | −12.49 | −11.91 | **0.65** | 0.63 / 0.15 / 0.32 |
| A4V | unfolded | −16.31 | −16.33 | −16.39 | 0.08 | 0.10 / 0.10 / 0.10 |
| G93A | folded | +9.71 | +9.73 | +9.67 | **0.06** | 0.05 / 0.09 / 0.07 |
| G93A | unfolded | +8.55 | +8.50 | +8.44 | 0.11 | 0.03 / 0.03 / 0.03 |
| I113T | folded | −27.41 | −28.08 | −27.61 | **0.67** | 0.15 / 0.20 / 0.34 |
| I113T | unfolded | −29.78 | −29.83 | −29.74 | 0.09 | 0.00 / 0.07 / 0.28 |

Per-replicate ΔΔG — A4V [3.75, 3.84, 4.48]; G93A [1.15, 1.23, 1.23]; I113T [2.37, 1.75, 2.13].
Independent samples after decorrelation: 13844 / 13531 / 13493 of 16308 raw.

## What these three establish

1. **All three report `converged: true`** (closures 0.63, 0.09, 0.34, all under the 1.0
   cap) **and all three are wrong**, by 0.8–2.4 kcal/mol, across three sites and three
   mutation types. The convergence diagnostics carry no information about accuracy.
2. **Errors are opposite in sign**, so there is no constant offset to calibrate away.
3. Errors are 10–20% of a single leg (ΔΔG is a small difference of two ~9–30 kcal/mol
   leg free energies), which is the scale of force-field and reference-state
   approximations, not of sampling noise.
4. **Replicate spread tracks repacking demand, not correctness.** A4V (buried Ala→Val)
   and I113T (interface Ile→Thr) both need neighbours to move: noisy folded leg, tight
   unfolded leg, and both *overestimate*. G93A (Gly→Ala, minimal steric change) has a
   folded leg tighter than its unfolded one — and *underestimates* by 1.22.
5. Consequently a small spread is **not** evidence of convergence. All replicates share
   one solvated box and one minimised structure (`stable_seed(variant, leg, "system")`
   has no `rep` in it); only velocities differ. If the slow mode is inaccessible in 3 ns,
   every replicate is stuck identically and the spread collapses while the mean stays
   biased. `ddg_err` is therefore not an uncertainty on ΔΔG for any of these three.

## Status

Diagnostics only — **not gate data**, and not mixable with post-fix windows
(`analyze._check_single_protocol` enforces this; these windows predate the fingerprint
and count as `unlabelled`). Raw `ddg.json` and `results/convergence/<v>.json` are kept
on the SCC under `results/prefix_diagnostics/`.
