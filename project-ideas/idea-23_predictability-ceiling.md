# Idea 23 — The Predictability Ceiling: How Much ΔΔG Is Knowable at All?

## Research question

Every predictor (FoldX, Rosetta, FEP, ML) is judged by correlation with experiment. But
**experiment itself is noisy** — the same mutation measured by two labs often differs by
1–2 kcal/mol. This project asks the deeper question:

**What is the information-theoretic ceiling on ΔΔG prediction imposed by experimental noise
itself?** Concretely: use *repeated measurements of the same mutations* (multiple labs /
methods in ProTherm, or the parent repo's `exp_source`/`exp_ddg_err` columns) to estimate the
experimental noise floor σ_exp; then compute the maximum achievable Pearson/MAE of any
predictor against that noisy target. The headline: **"even a perfect predictor can't exceed
R ≈ 0.9 against experiment — so a predictor scoring R=0.85 is closer to perfect than it
looks."**

This is a **meta-science / information-theory** project: the deliverable is a corrected
"intrinsic performance" of predictors, which changes how the parent repo's validation gate
(and the whole field) should be interpreted.

## Why this needs an SGE cluster

```
tasks = resampling / bootstrap jobs over repeated-measurement data
```

- Estimating σ_exp and the noise-corrected ceiling requires **resampling statistics**:
  bootstrapping lab/method replicates and fitting measurement-error models — thousands of
  independent resamples, each an array task.
- You also generate *synthetic* "true ΔΔG + known noise" datasets to calibrate your ceiling
  estimator (a Monte Carlo grid over (noise level, effect size) — embarrassingly parallel).
- Combining datasets (ProTherm, DMS, the parent repo's controls) to estimate per-mutation
  σ_exp is data-engineering that runs naturally as a pipeline.

## Data & tools

- **Repeated-measurement data:** ProTherm (many mutations measured by multiple labs), DMS
  datasets with replicate lanes, and — conveniently — the parent repo's `data/variants.csv`
  (which has both `exp_ddg` and `exp_ddg_err`, plus multiple `exp_source` per variant).
- **The parent repo's gate numbers are the perfect worked example:** A4V, G93A etc. have
  literature scatter; you can measure their σ_exp directly.
- **Methods:** measurement-error modeling (e.g., hierarchical Gaussian model on repeated
  measurements), bootstrap CI on the noise-corrected ceiling, maybe a small Stan/PyMC fit.

## Skill prerequisites

- Solid stats (bootstrap, variance decomposition).
- Python/pandas for the dataset work.
- Intermediate. No MD required — a data-science-first project with biophysics meaning.

## Cluster budget

| Parameter | Value |
|---|---|
| Resampling jobs | ~5,000 bootstrap/MC tasks |
| Per-task | seconds–minutes CPU |
| Wall-clock | **~1 day** |

## Milestones

1. Assemble repeated-measurement ΔΔG data (ProTherm + variants.csv + a DMS).
2. Estimate per-mutation / per-method σ_exp; quantify the distribution of inter-lab scatter.
3. Build a noise model: observed = truth + N(0, σ_exp²); estimate the noise-corrected ceiling
   on Pearson/MAE via simulation.
4. Re-score the parent repo's predictors (FoldX, FEP) against the *noise-corrected* target:
   the "attainable fraction" metric (achieved R / ceiling R).
5. Sensitivity: how does the ceiling change if experiment were 2× cleaner? (A compelling
   "what if labs measured better" plot.)
6. Write it up as a mini meta-analysis: "how to read ΔΔG predictor performance honestly."

## Deliverables

- **The ceiling figure**: maximum achievable Pearson vs σ_exp, with real predictors plotted
  against it — the money figure.
- A noise-floor report: how reproducible are ΔΔG measurements, really (per protein/method)?
- A corrected performance table: predictors' achieved fraction of the ceiling.
- A concrete recommendation for the parent repo's `max_rmse_kcal: 1.5` gate in light of the
  noise floor.

## Pitfalls

- **Don't confuse ceiling with excuse.** A high noise floor explains but does not excuse a
  predictor's errors; report the *fraction achieved*, which punishes genuinely bad predictors.
- **σ_exp must come from repeated measurements of the SAME mutation**, not from correlated
  variants — mis-estimating it poisons everything. Audit the data first.
- **Survivorship bias in ProTherm.** Negative (non-destabilizing) results are under-reported;
  note it and consider a DMS with more complete coverage.
- **Beware the parent repo's own "Kumar 2017" sourcing issue** (citations.md): use only
  confirmed repeated measurements.
- **Keep it honest:** the ceiling is a target, and "perfect predictor R=0.9" is a *statement
  about noise*, not about physics. Say so in one line at the top of the write-up.
