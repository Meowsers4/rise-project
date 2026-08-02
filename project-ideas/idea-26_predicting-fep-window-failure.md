# Idea 26 — Predicting FEP Window Failure Before You Pay for It

## Research question

A large FEP array will have windows that fail: `mdrun` blows up (end-state singularities),
the u_kn overlap is so poor that MBAR gives garbage, or a window never equilibrates. On a
3,000-window array, even a 5% failure rate wastes 150 GPU-jobs' worth of compute — and, worse,
failing windows are often silently excluded, biasing the result. The question:

**Can a classifier predict which (variant, leg, window, replicate) will fail or fail to
converge — from cheap features available BEFORE the expensive mdrun runs?**

Features to test (all cheap): mutation type (buried/surface, charge change), λ-window position,
soft-core parameters, force field, topology complexity (number of perturbed atoms), and a short
10 ps "probe" mdrun's early diagnostics (energy jumps, max force, end-state energy).

## Why this needs an SGE cluster

```
tasks = all windows you will run (the failure labels) + probe runs + training
```

- The *labels* come from running a big real FEP array and recording which windows
  fail/converge poorly — the cluster's normal job, but with a failure-logging twist.
- To build the probe-based features you run a **tiny 10 ps mdrun per window** (cheap, CPU/GPU
  but ~100× cheaper than the real one) — another array, but 100× smaller.
- The classifier training/validation is a small grid (`qsub -t`) on top.
- The payoff is a **scheduling policy**: if a probe predicts failure, skip the window, retry
  with different soft-core/λ, or run it on a different replicate — a "prediction + mitigation"
  loop that only matters at the cluster's scale.

## Data & tools

- **Labels:** failure = mdrun crash / NaN energies, OR "MBAR garbage" = poor overlap or
  implied-uncertainty blow-up; non-convergence = replicate scatter >> MBAR error (Idea 25).
  Collect these from the parent repo's window output with provenance.
- **Features:** mutation/structural descriptors (burial via SASA, charge change, residue pair),
  λ value, soft-core config, perturbed-atom count, probe-run diagnostics (max force, energy
  range, RMSD drift in 10 ps).
- **Models:** logistic regression / gradient boosting / random forest — interpretability
  matters (you want to know *why* a window is at risk).
- **Mitigation evaluation:** if predicted-failure windows are re-run with a fix (e.g.,
  different soft-core, more equilibration), what fraction are rescued?

## Skill prerequisites

- Comfortable running the parent repo's FEP windows and reading their output.
- Python/ML basics (Idea-14 level).
- Intermediate-advanced; this is "applied reliability ML on simulation data."

## Cluster budget

| Parameter | Value |
|---|---|
| Label-generating array | ~500 windows (real) |
| Probe array | ~500 × 10 ps (cheap) |
| Training grid | ~100 configs |
| Wall-clock | **~1–2 weeks** (real windows) + **~1 day** (probes/training) |

## Milestones

1. Run a real FEP panel; log failures and convergence quality per window (this is your label
   set — smallish is fine, ~200–500 windows).
2. Extract the cheap features (including probe runs) for every window.
3. Train the failure classifier; report AUROC and, crucially, **precision at the failure
   fraction** (how many flagged windows are real failures).
4. Feature importance: which cheap feature predicts failure best (my bet: end-state/soft-core
   interplay + charge/burial)?
5. **Mitigation test:** re-run predicted-failures with a documented fix; measure rescue rate.
6. Verdict: at what cluster scale does this pre-screening pay for itself in saved GPU-hours?

## Deliverables

- **Failure-prediction ROC/precision curve** — the money figure.
- Feature-importance table ("what predicts a dead window").
- A rescue-rate result (mitigations that actually save windows).
- A written recommendation for the parent repo's `submit_array.sh`: a two-stage (probe →
  real) submission policy.

## Pitfalls

- **Define "failure" precisely before labeling.** Crash ≠ bad-overlap ≠ non-converged; they
  need different mitigations and different predictors. Do not merge them into one label.
- **Survivorship bias in labels.** If your pipeline already auto-resumes or re-tries, your
  label set is clean by construction — document exactly what "failed" means in the logging.
- **Probe features must be cheap.** The whole point is pre-screening; a probe that costs 1/3
  of the real run defeats the purpose. Keep it ≤ a few % of real cost.
- **Small label sets.** With ~200 failures, don't fit a deep net — a regularized linear/GBM
  model with honest CV is the right tool.
- **Don't over-fit the mitigations.** A fix that rescues 90% of flagged windows in your data
  must be re-tested on windows it wasn't tuned on.
