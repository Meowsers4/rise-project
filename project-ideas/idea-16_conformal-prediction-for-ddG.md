# Idea 16 — Conformal Prediction for ΔΔG: Error Bars You Can Trust

## Research question

Every ΔΔG in this pipeline ships with an error bar — but how trustworthy is it? MBAR's
asymptotic error assumes the sampling was converged, and replicate-based error bars assume the
replicates are a fair sample of the protocol's noise. **Conformal prediction** (a distribution-
free ML framework) wraps any predictor and produces prediction intervals with a *guaranteed*
coverage probability — no distributional assumptions about the model's errors.

The research question: **Can conformal prediction turn FoldX or FEP ΔΔG predictions into
calibrated "85% of true values fall in this interval" statements that hold in practice?** And
the follow-up that makes it a biophysics project rather than a stats exercise: **do
conformality and MBAR/replicate error bars agree or disagree, and which one is honest when
they differ?**

## Why this needs an SGE cluster

```
tasks = mutations × predictors (FoldX, Rosetta, FEP-window subsets)
```

- Calibrating a conformal interval needs a **calibration set** of size ~1,000+ — exactly the
  saturation-scan arrays from Idea 01. One array task per mutation gives you the calibration
  data in hours instead of weeks.
- **Conformalization with out-of-distribution coverage** requires *many* re-scored data points
  so you can measure empirical coverage in held-out slices (buried vs surface, charged vs
  neutral). The cluster's array of thousands of mutations is what makes the coverage curve
  statistically meaningful.
- For the FEP leg, you need *per-mutation replicate subsets* (1-rep, 2-rep, 3-rep pools) so
  you can measure how interval width and coverage change with compute — a resampling loop that
  is embarrassingly parallel.

## Data & tools

- **Predictor:** FoldX (`BuildModel`, the Idea-01/02 machinery), Rosetta `cartesian_ddg`, or
  the parent repo's FEP `ddg.json`.
- **Labels:** experimental ΔΔG from `data/variants.csv` / ProTherm / a DMS dataset (Idea 05).
- **Method:** `mapie` (Python conformal prediction library — free, MIT) or a ~50-line
  split-conformal implementation you write yourself.
- **Conformal variants to compare:** split-conformal (simple), jackknife+, and
  conformalized quantile regression (CQR) for heteroscedastic intervals.

## Skill prerequisites

- Python + pandas.
- Comfort with the concept of "coverage": `P(true ∈ interval) ≥ 1 − α`.
- Intermediate. The statistics are the hard part, not the compute.

## Cluster budget

| Parameter | Value |
|---|---|
| Calibration mutations | ~1,000–2,000 (FoldX) |
| FEP replicate subsets | 8 gate variants × 3 rep-pools |
| Wall-clock | **~1–2 days** (FoldX arrays) + small GPU |

## Milestones

1. Generate (or reuse) a FoldX saturation scan over a protein with experimental ΔΔG labels.
2. Split data into train / calibration / test (by **position**, not mutation — avoid leakage).
3. Train a base predictor (e.g., gradient boosting on the Idea-01 features), then build
   conformal intervals on the calibration set.
4. Measure **empirical coverage** on the test set overall AND in slices (buried vs surface,
   large vs small substitutions). Do the intervals hit the nominal rate?
5. Repeat on the FEP leg: build conformal intervals from the replicate pools; compare interval
   widths to MBAR errors on the same mutations.
6. **The verdict:** where do MBAR/replicate error bars and conformal intervals disagree? Is
   the disagreement concentrated in particular mutation classes (the real finding)?

## Deliverables

- Coverage-vs-nominal calibration plot (the money figure).
- Interval-width comparison: MBAR vs replicate vs conformal, per mutation class.
- A "which error bar do you trust, and when" recommendation — publishable as a small methods note.

## Pitfalls

- **Leakage again.** If the same mutation sits in calibration and test, coverage is fake. Split
  by position/protein.
- **Conformal ≠ correct.** It guarantees coverage of the *predictor's* errors, not of physics.
  If FoldX is biased +1 kcal/mol everywhere, conformal intervals are honest-but-offset. Report
  bias separately from interval coverage.
- **Heteroscedasticity.** Mutation error is not constant (buried mutations are harder); plain
  split-conformal gives over-wide intervals for easy cases. Use CQR to get width that adapts.
- **Small FEP sets.** 8 gate variants cannot calibrate anything by themselves — conformalize
  the *cheap* predictor at scale, and only *test* conformality on the expensive FEP points.
- **Pre-register the α.** Decide the nominal coverage before running, or you'll be tempted to
  pick α after seeing the plot.
