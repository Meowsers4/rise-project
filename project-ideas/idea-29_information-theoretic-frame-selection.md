# Idea 29 — Information-Theoretic Frame Selection: Which Frames Actually Matter for MBAR?

## Research question

A 3 ns FEP window saved every 20 ps produces 150 frames per window, but **most are redundant**.
Consecutive frames are correlated (decorrelation typically reduces usable samples to ~20–50%),
and some frames carry far more information about the free-energy difference than others (frames
at rare/fast-mixing conformations are worth more). The parent repo already does
`decorrelate: true` (subsample to independent samples). This project asks the sharper
information-theoretic question:

**Which frames, selected optimally, minimize the MBAR estimator variance for a given number of
frames kept?** You compare: (a) standard decorrelation, (b) uniform thinning, (c) greedy
variance-minimizing selection, (d) frames chosen to maximize u_kn diversity. The payoff is a
**retention policy**: with `trajectory_retention: estimates_only`, how few frames can you keep
without inflating the ΔΔG error?

## Why this needs an SGE cluster

```
tasks = windows (real FEP) + per-window selection experiments + MBAR refits
```

- You need *real* windows to test selection (the parent repo's array). The frame-selection
  experiment is a resampling study: for each window, try many selection strategies and many
  budget fractions, refitting MBAR each time — thousands of parallel refits.
- To measure "which frames matter," you need *more than 3 ns* of one window to have enough
  frames to subsample meaningfully (say 6 ns × one variant), which is cluster compute.
- The mock harness (known answers) lets you develop the selection criterion cheaply first.
- The result directly feeds a cluster-ops decision (how much disk/retention you need), which
  is exactly the kind of finding that impresses judges because it's *actionable*.

## Data & tools

- **Data:** u_kn matrices from real (or mock) windows; MBAR via `pymbar`.
- **Selection criteria to test:**
  - Baseline: `pymbar.timeseries` decorrelation (what the repo does).
  - Uniform thinning (every kth frame).
  - **Variance-minimizing selection**: drop the frame whose removal least increases the MBAR
    variance estimate (greedy, or a heuristic like importance weighting by u_kn dispersion).
  - **Representative sampling**: cluster frames in u_kn space, keep one per cluster.
- **Metrics:** final ΔΔG error vs the full-data answer; MBAR variance; per-window sample
  efficiency (effective samples per byte stored).

## Skill prerequisites

- Comfortable with the parent repo's window NPZ schema and `pymbar`.
- Some information-theory / estimator intuition (variance of estimators, importance).
- Intermediate-advanced.

## Cluster budget

| Parameter | Value |
|---|---|
| Windows (long, real) | 1–2 variants × 2 legs × 18 windows (6 ns each) |
| Selection experiments | ~10 strategies × 6 budgets × 18 windows × 2 legs |
| MBAR refits | ~5,000–20,000 (cheap CPU) |
| Wall-clock | **~1 week** (real MD) + **~1 day** (refits) |

## Milestones

1. Generate long real windows (6 ns) for one variant pair (or develop on mock data).
2. Implement the selection strategies; for each, subsample frames at several budgets
   (10%–90%) and refit MBAR → per-strategy error vs budget curve.
3. **The money comparison:** error-vs-frames-kept for all strategies. How far above baseline
   is variance-minimizing selection? How much disk does it save at equal error?
4. Sanity: confirm the full-data answer is reproduced by all strategies at high budget (no
   selection bug).
5. Generalize to a second variant/leg to check the result isn't a fluke.
6. Verdict + a concrete retention-policy recommendation for the parent repo (how much to keep,
   and how to choose which frames).

## Deliverables

- **Error-vs-frames-kept curves** per selection strategy — the money figure.
- The best-strategy retention recipe (keep X% of frames, selected by Y) with the disk/error
  trade-off quantified.
- A short note tying it to `trajectory_retention: estimates_only` and to how MBAR's decorrelation
  actually compares to the fancier strategies.

## Pitfalls

- **Variance estimates are noisy.** Comparing strategies by "MBAR variance" is circular if you
  use the same variance estimator to select frames; use a *held-out* measure (e.g., the final
  ΔΔG against the full-data answer, or bootstrap).
- **Don't over-fit the selection to one window.** Selection that wins on window 3 must win on
  the held-out windows too; report across windows.
- **Consecutive-frame correlation must be respected.** A greedy "keep the informative frames"
  scheme that keeps adjacent frames is just data at full correlation — the whole point is
  that selection ≠ thinning, so show independence (autocorrelation of kept frames).
- **The mock harness's smoothness is not real MD.** Real u_kn has rough, non-Gaussian
  structure; validate the winner on real data before believing it.
- **Disk savings vs correctness.** A policy that saves 80% disk but adds 0.3 kcal/mol error is
  a trade, not a win; report both axes.
