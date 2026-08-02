# Idea 19 — Markov State Models of the SOD1 Folding Landscape

## Research question

MD gives you trajectories, but a trajectory is just a time series. **Markov State Models
(MSMs)** turn many short trajectories into a *kinetic network*: a discrete set of metastable
conformational states plus transition probabilities between them, from which you can extract
the **free-energy landscape, mean first-passage times between states, and the dominant
transition paths**. Applied to the parent repo's Stage 4 mechanism MD, the question is sharp:

**Do disease variants change the *kinetic network* of SOD1 — not just its stability — in a way
that single-trajectory analysis misses?** E.g., do mutants populate a near-unfolded metastable
state more often, or make the folding transition slower, or open new pathways toward
aggregation-prone conformations?

## Why this needs an SGE cluster

```
tasks = many short independent MD trajectories (the MSM training data)
```

- MSMs **require many independent trajectories** (dozens to hundreds), because they estimate
  transition probabilities from *counts* of transitions between states — and you need enough
  transitions to estimate the rates. One long trajectory is worse than 50 short ones for MSM
  building (this is exactly the "replicate discipline" of README rule 5, mechanized).
- **1 array task = 1 short trajectory** (say 100–200 ns). 60 trajectories per variant × 4
  variants = 240 GPU jobs — the cluster's embarrassingly-parallel sweet spot, exactly the
  `qsub -t` geometry the parent repo already uses.
- The MSM construction/validation (clustering, lag-time testing, Chapman-Kolmogorov checks) is
  cheap CPU work on the aggregated data — the cluster's job is feeding it.

## Data & tools

- **MD engine:** GROMACS (parent stack). Plain unbiased MD of the apo-reduced monomer (the
  parent repo's species) at 300 K; optionally a second temperature.
- **MSM software (free):** **PyEMMA** (standard, maintained) or **MSMBuilder**. Both do
  featurization → dimension reduction → clustering → MSM → validation → implied timescales.
- **Order parameters:** native contacts (Q), RMSD, radius of gyration; feed to time-lagged
  independent component analysis (tICA) in PyEMMA to get the slow collective coordinates.
- **Interpretation:** stationary distribution, free-energy projection onto tICs, committor
  analysis to find transition states, and per-variant comparison of these.

## Skill prerequisites

- Comfortable running MD (Idea 08/09 level).
- Willing to learn MSM concepts: metastability, implied timescale, Chapman–Kolmogorov test.
- Intermediate-advanced. The analysis is genuinely data-science flavored.

## Cluster budget

| Parameter | Value |
|---|---|
| Trajectories | 60/variant × 4 variants (WT + 3 mutants) = 240 |
| Per-trajectory | ~100–200 ns GPU (~hours) |
| MSM construction | CPU, minutes–hours after collection |
| Wall-clock on 8 GPUs | **~1–2 weeks** |

## Milestones

1. Build one WT system (parent prep); verify a 20 ns test run is stable (no unfolding, sane
   RMSD).
2. Generate the trajectory array (60 per variant, distinct seeds). Use the parent repo's
   `trajectory_retention` policy: keep order parameters per frame, downsample frames.
3. Featurize (Q / contacts / tICA input) and cluster → MSM per variant. **Validate each MSM**
   (implied timescale plot; Chapman–Kolmogorov).
4. Extract per-variant: stationary distribution, free-energy landscape (tICA projection),
   dominant transition paths, mean first-passage time folded→near-unfolded.
5. **Compare variants:** which landscape features change? Do mutants populate a high-energy
   "pre-misfolding" basin more? Is the barrier (rate) different even when ΔΔG (depth) is the
   same?
6. Cross-check a headline claim with a plain-trajectory statistic (e.g., RMSF) so the MSM
   result is anchored.

## Deliverables

- **Per-variant free-energy landscapes** (tICA projections with basin labels) — the money
  figures.
- A kinetic-network diagram (nodes = states, arrows = rates) per variant.
- A comparison table: state populations, barrier heights, mean first-passage times, variant
  vs WT.
- A verdict on whether kinetic (not just thermodynamic) signatures separate the mutants.

## Pitfalls

- **Garbage-in: MSMs need stationary, well-sampled trajectories.** If trajectories are too
  short to revisit states, the transition counts are meaningless. Check implied-timescale
  *convergence with lag time* before trusting anything.
- **Do NOT report a single trajectory's story.** The whole point is ensemble statistics; a
  "path" you see in one trajectory is anecdote until the MSM says it has significant flux.
- **State-count sensitivity.** Results depend on how many clusters you keep; vary it and show
  conclusions are robust.
- **Choice of features changes the landscape.** tICA projections are only as good as the input
  features; validate on a physical order parameter (Q) as a cross-check.
- **Trajectory retention policy.** MSMs need per-frame features but not raw frames; extract
  features at the cluster and delete most frames (disk fills fast — CLAUDE.md rule).
