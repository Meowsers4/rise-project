# Idea 25 — Are FEP Error Bars Honest? A Bootstrap Audit of MBAR

## Research question

The parent repo's rule 5 is non-negotiable: "A ΔΔG without an uncertainty and a cycle-closure
check is not a result." But what is the *quality* of that uncertainty? MBAR's error bars are
**asymptotic** (they assume you sampled long enough and that the u_kn overlap is good), and
replicate-based errors assume replicates are a fair sample of protocol noise. When MBAR says
±0.3 kcal/mol but two independent replicas disagree by 1.2 kcal/mol — which is the real error?

The question: **Do the error bars that FEP pipelines report actually cover the truth?** You
audit this with **bootstrap resampling** of the window data (re-blocking trajectories,
resampling windows, resampling replicas) and compare three error estimates on the same data:
(1) MBAR's analytic error, (2) block-bootstrap error, (3) replicate scatter. The result is a
calibration report: **"your pipeline's error bars are [honest / optimistic / pessimistic] by
a factor of X, and here's where they break."**

## Why this needs an SGE cluster

```
tasks = bootstrap resamples × windows (thousands of MBAR refits)
```

- Each bootstrap estimate requires **re-fitting MBAR on resampled data thousands of times** —
  the classic embarrassingly-parallel resampling loop (1 task per resample). On a laptop this
  is slow; on an array it's minutes.
- To audit honestly you need *multiple replicas of every window* (the parent repo's `replicates:
  3`), which means the underlying FEP array is real compute — and you should run **additional
  replicas** (say 5) for a subset so the audit has enough replicate scatter to measure.
- The block-bootstrap variant (resampling trajectory blocks at several block sizes to find the
  convergence plateau) multiplies the refits again — all parallel.

## Data & tools

- **Data:** the parent repo's FEP window NPZs (real runs, 3–5 replicates) plus the mock
  harness for developing the audit.
- **Tools:** `pymbar` (MBAR + its analytic errors), `sklearn`/`scipy` for block bootstrap,
  custom re-blocking (overlapping vs non-overlapping blocks).
- **Error estimators to compare:**
  1. MBAR asymptotic error (what the pipeline reports).
  2. **Block bootstrap** over the time series in each window (accounts for correlation).
  3. **Replicate bootstrap** (resample the 3–5 replicas per window).
  4. **Chain bootstrap** (resample whole replica trajectories) — the most conservative.

## Skill prerequisites

- Comfortable with the parent repo's FEP output and with `pymbar`.
- Solid stats: bootstrap, bias, coverage, the difference between variance and bias.
- Intermediate-advanced (this is a statistics-of-simulation project).

## Cluster budget

| Parameter | Value |
|---|---|
| Windows | 8 variants × 2 legs × 18 windows × 5 reps = 1,440 (real GPU) |
| Bootstrap refits | ~5,000 per audit × a few audit configs |
| Wall-clock | **~1–2 weeks** (FEP) + **hours** (refits) |

## Milestones

1. Run (or reuse) a real FEP panel with ≥5 replicates for a small variant set (a few controls).
2. Implement the four error estimators; develop on mock data (where you know the true error).
3. Compute per-variant ΔΔG ± error under each estimator. **The audit table**: MBAR error vs
   bootstrap error vs replicate scatter, per variant.
4. The calibration check: if replicates are a fair sample of protocol noise, the replicate
   scatter should (approximately) bracket MBAR's error. Where do they diverge?
5. Sensitivity: does the block-bootstrap plateau at a larger error than MBAR reports (the
   classic "MBAR is optimistic" result)?
6. Verdict + a concrete recommendation for the parent repo: **is `replicates: 3` enough?** Does
   the gate need a minimum-bootstrap-error guard?

## Deliverables

- **The audit table** (error estimate per estimator per variant) and a coverage-style figure
  (does the ±1σ interval include the replicate mean).
- A "bias-vs-variance" decomposition of the pipeline's error.
- A written recommendation: how the parent repo should report error bars and whether to raise
  the replicate floor.

## Pitfalls

- **Don't audit on mock data only.** The mock harness has known answers; the point is that
  *real* FEP data has unmodeled correlation, and that's exactly what the bootstrap exposes.
- **Block size choice.** Non-overlapping block bootstrap is unstable at small block counts;
  use several block sizes and report the plateau, not one number.
- **Replicate count is itself a source of noise.** With 3 replicates the "replicate scatter"
  estimate has huge variance; that's a finding (the floor is too low), not an excuse.
- **MBAR's analytic error is not "wrong"** — it's conditional on the model (sampling is
  converged, u_kn correct). Frame the comparison as "under which conditions do they agree,"
  not "MBAR is a lie."
- **Cost control:** the audit's real compute is extra replicates. Run 5 reps on a *subset* (the
  gate variants), not the whole panel.
