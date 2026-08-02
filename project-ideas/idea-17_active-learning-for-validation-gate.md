# Idea 17 — Active Learning for the Validation Gate

## Research question

The parent repo's gate rule is: **run FEP on controls first; only after the gate passes do you
spend GPU hours on uncharacterized variants** (README §2.2). But *which* controls should you
run FEP on first? Running all 54 controls is expensive; running the 8 chosen gate variants is
a heuristic. **Active learning** (a.k.a. Bayesian experimental design) answers: *which next
variant maximizes the information gained about the gate decision per GPU-hour spent?*

Concretely: build a cheap surrogate (FoldX/Idea-14 ML) over all controls, then use a query
strategy — uncertainty sampling, expected model change, or maximizing the expected change in
the gate's Pearson/RMSE — to rank which control to run FEP on next. Simulate the process
**in silico**: given a fixed GPU budget, does active FEP selection beat running the gate
variants blindly, and beat random selection?

## Why this needs an SGE cluster

```
tasks = candidate FEP windows (variant × leg × window × rep) × active-learning rounds
```

- Active learning is *iterative*: choose → run expensive FEP → update surrogate → choose
  again. Each round is a batch of GPU windows (the parent repo's array). A 10-round experiment
  with 2 variants/round = 20 variants of FEP — the cluster is the only way to afford the *loop*.
- To *evaluate* active learning you must compare it to baselines (random, greedy, fixed gate
  set) — each baseline is a whole FEP experiment. That's 3–4× the compute of one run, all
  embarrassingly parallel.
- The surrogate retraining is cheap (CPU, one task per round); the FEP is the expensive
  step and the cluster makes it schedulable per-round.

## Data & tools

- **Surrogate:** gradient boosting or GP on the Idea-01/14 features (FoldX ΔΔG + structural
  descriptors), fit on all controls.
- **Expensive oracle:** the parent repo's GROMACS + pmx FEP engine (`src/fep/window.py`).
- **Query strategies to compare:** (a) uncertainty sampling (largest surrogate variance),
  (b) expected improvement of the gate metric, (c) a *fixed* pre-chosen gate set (baseline),
  (d) random selection (baseline).
- **Simulation layer:** replay active learning with *already-computed* FEP labels if they
  exist, or with mock labels (FEP_MOCK) for the harness, before spending real GPU time.

## Skill prerequisites

- Python; the Idea-14 ML stack (boosting, features).
- Understands "gate = Pearson/RMSE vs experimental controls" (config `validation.min_pearson`).
- Intermediate-advanced: this is ML × experimental-design thinking.

## Cluster budget

| Parameter | Value |
|---|---|
| Rounds | 6–10 |
| Variants per round | 2–4 |
| FEP windows per variant | 2 legs × 18 windows × 3 reps = 108 |
| Wall-clock on 8 GPUs | **~2–4 weeks** (or cut variants/windows) |

## Milestones

1. Define the gate metric precisely (Pearson + RMSE on the chosen controls, from config).
2. Build the surrogate on all controls (FoldX/structural features → predicted FEP ΔΔG).
3. **Simulate** active learning with mock FEP labels: do the query strategies differ in how
   quickly the surrogate's ranking of variants stabilizes? (Cheap, no GPU.)
4. Run round 1: FEP the top queried variants (real). Update surrogate. Record gate metric.
5. Repeat rounds; at the end, compare final surrogate ranks vs a *hold-out* set of FEP-verified
   variants (compute a few extra, never-queried variants as the test set).
6. **The verdict:** did active selection reach the gate decision with fewer GPU-hours than the
   baselines? How many controls does the gate *really* need (power analysis)?

## Deliverables

- Learning curves: gate metric (or surrogate error) vs GPU-hours for each query strategy.
- A "query map": which variants did active learning prioritize, and why (structural/chemical
  rationale for the queried mutations).
- A power-style conclusion: the minimum controls needed to make the gate decision confidently.

## Pitfalls

- **The surrogate must not leak the oracle.** Never train the surrogate on variants it's about
  to query. Strict per-round discipline is the whole experiment.
- **Active learning optimizes surrogate uncertainty, not biological truth.** A variant that's
  "uncertain in FoldX space" may be trivially predictable in FEP space. Add a *random* control
  arm so you can measure whether querying was worth it.
- **Gate metric noise.** Pearson on 8 points is noisy; use RMSE as the primary target (it
  averages better), and report both.
- **GPU budget realism.** 10 rounds of real FEP is a lot; simulate first, then run a
  *reduced* version (3 rounds, 2 variants) as the demonstrator.
- **Pre-register the query strategies.** Otherwise the "best strategy" is just whichever you
  tuned after seeing results.
