# Idea 09 — Kinetic Stability from Unfolding Replicates

## Research question

**Which protein variants unfold fastest?** Thermodynamic stability (ΔΔG, the parent repo's
metric) tells you the *depth* of the free-energy well; **kinetic stability** tells you how
*long* a protein stays folded. Two proteins can have the same ΔΔG but very different
unfolding rates — and for aggregation-prone disease proteins (like SOD1 in ALS), the
**time-to-first-unfold** is arguably the more disease-relevant quantity.

This project runs **replicate unfolding simulations** (high temperature or chemical-denaturant
mimics, e.g., pulling the protein apart or heating it) on a panel of variants and measures
the distribution of unfolding times: **survival curves** `P(folded | t)`, mean first-passage
times, and whether kinetic stability ranks variants the same as thermodynamic ΔΔG.

## Why this needs an SGE cluster

```
tasks = variants × replicates  (10 × 10 = 100 GPU runs)
```

- Unfolding is a **rare, stochastic event** — you cannot watch it happen once; you must
  average over many replicate trajectories, each several hundred ns of GPU mdrun. That is
  10–20 GPU-days per variant → the cluster is the only way to get statistically meaningful
  survival curves.
- **1 array task = 1 (variant, seed)** — identical to the parent repo's replicate discipline
  (README rule 5: replicates and error bars are mandatory).
- The fan-out is perfectly load-balanced (every run is the same length and cost).

## Data & tools

- **Engine:** GROMACS + pmx (parent repo stack). For unfolding you typically:
  - run at **high temperature** (e.g., 400–500 K) to accelerate unfolding, and/or
  - add a **bias/pull** (steered MD) for measurable unfolding rates, or
  - use **replica exchange / parallel tempering** for better sampling.
- **Panel:** a set of variants with known experimental ΔΔG (reuse the parent repo's gate
  subset + a few benign controls).
- **Analysis:** define an order parameter (RMSD from native, fraction of native contacts Q,
  radius of gyration); threshold → "unfolded"; compute survival curves and mean first-passage
  time; correlate with ΔΔG.

## Skill prerequisites

- Solid MD literacy (can set up, run, and debug GROMACS).
- Survival-analysis / statistics basics (Kaplan–Meier-style curves, exponential fits).
- Intermediate-advanced; do Idea 02 or 08 first.

## Cluster budget

| Parameter | Value |
|---|---|
| Variants | 8–10 (incl. 2–3 stable controls) |
| Replicates | 10 per variant |
| Per-run | ~200–500 ns GPU (hours to ~12 h) |
| Wall-clock on 8 GPUs | **~1–2 weeks** |

## Milestones

1. Reproduce the parent repo's folding setup for WT; verify the protein stays folded at 300 K
   (the negative control must NOT unfold — otherwise the whole protocol is broken).
2. Calibrate temperature: run 2 replicates of WT at 350/400/450 K; find the temperature where
   WT unfolds on the 100–500 ns timescale (that's your working point).
3. Pick the unfolding metric (Q-fraction is robust) and a "folded" threshold.
4. Array over (variant × seed × temperature if sweeping); save only per-frame order
   parameters (tiny files — `trajectory_retention: estimates_only` discipline).
5. Build survival curves per variant; fit mean first-passage time (or median) ± error.
6. **The science:** correlate kinetic unfolding times vs experimental ΔΔG. Do the kinetic and
   thermodynamic ranks agree? Where do they diverge (the interesting cases)?

## Deliverables

- **Survival curves** `P(folded) vs t` for all variants (with shaded confidence bands).
- Kinetic unfolding time vs experimental ΔΔG scatter — the money figure.
- A "kinetic vs thermodynamic stability" comparison table with the divergent variants flagged.
- A mechanistic note on how the fastest-unfolding variant loses structure (which contacts
  break first).

## Pitfalls

- **Calibrate the driving force first.** If no variant unfolds in your window, or WT unfolds
  instantly, the data is worthless — the temperature/force must be tuned on WT before the
  array.
- **High temperature changes the mechanism.** Unfolding at 450 K may not be the same pathway
  as at body temperature. Report it as "accelerated unfolding", never as "the" unfolding
  pathway, unless you cross-check at a lower T.
- **Survival curves need censoring.** A run that never unfolds within the wall time is
  "right-censored", not a failure — handle it properly (Kaplan–Meier) or your mean times are
  biased short.
- **Order-parameter threshold arbitrariness.** Report results for ≥2 thresholds; the
  conclusions should be robust to a ±10% change.
- **Replicates, replicates, replicates.** 3 is a floor, not a target; for survival curves you
  want 10+.
