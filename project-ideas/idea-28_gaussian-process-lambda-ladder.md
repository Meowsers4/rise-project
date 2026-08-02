# Idea 28 — Gaussian-Process Surrogates for the λ-Ladder: Smarter Windows

## Research question

The parent repo runs a **fixed λ-ladder** (18 equally spaced windows). But the free energy as
a function of λ is not uniformly smooth — some λ regions (e.g., where soft-core kicks in) need
denser sampling, others are nearly linear and waste windows. The question:

**Can a Gaussian process (GP) model of the running windows tell you adaptively where to add
(or remove) λ-windows, saving GPU time while keeping the MBAR estimate accurate?**

Concretely: fit a GP to `dG(λ)` (or the u_kn overlap structure) from the windows run so far,
compute an acquisition function (predicted variance per λ, or expected reduction in MBAR
error), and decide the next λ to add. Compare: fixed 18-window ladder vs GP-guided placement
at equal total windows, and at reduced windows.

## Why this needs an SGE cluster

```
tasks = windows (each is a GPU job) + GP fitting per round
```

- The whole point is *sequential*: run some windows → fit GP → choose next λ → run it. That's
  an active-learning loop where the GPU work is windows (arrays) and the decision step is cheap
  CPU.
- Evaluating the approach honestly requires **repeating the experiment** over many mutations
  (each is a mini active-learning episode) and comparing against the fixed ladder — each
  repeat is more GPU windows, all embarrassingly parallel across variants.
- You also need to *re-fit MBAR* many times (per episode, per candidate λ) — a cheap
  resampling-style load, parallelizable.
- The parent repo's mock harness (known-answer windows) is the perfect place to develop the
  GP before any real GPU spend.

## Data & tools

- **Data:** u_kn matrices from real (or mock) windows; MBAR gives dG and its variance at each λ.
- **GP tools:** `gpytorch` or `scikit-learn` GaussianProcessRegressor on dG(λ) (with λ and
  maybe soft-core-region indicator as inputs).
- **Acquisition:** predicted variance (uncertainty sampling), or the expected change in the
  final ΔΔG estimate (the truly Bayesian criterion, more work).
- **Baseline:** the fixed 18-window ladder with identical physics — the control you must match.

## Skill prerequisites

- Comfortable with the parent repo's FEP output and MBAR.
- Some Bayesian-ML comfort (what a GP posterior and acquisition function are).
- Advanced: this is methodology research.

## Cluster budget

| Parameter | Value |
|---|---|
| Variants for the repeat study | 6–8 |
| Windows per variant | ~18 (fixed) vs ~12–18 (adaptive) |
| GP/MBAR refits | ~50–200 per variant (CPU) |
| Wall-clock on 8 GPUs | **~1–2 weeks** |

## Milestones

1. Develop the GP+acquisition on the **mock** harness (where true ΔΔG is known): does adaptive
   placement recover the true value with fewer windows than the fixed ladder?
2. Lock the protocol: acquisition function, prior, stopping rule, and the *evaluation metric*
   (error of final ΔΔG vs fixed ladder at equal window count).
3. Run the real repeat study: for each variant, an adaptive episode and a fixed-ladder control
   (same physics, same total windows).
4. Compare final ΔΔG errors and per-λ sampling patterns.
5. **Interpretation:** where did the GP choose to put extra windows? Does it match the
   soft-core / end-state trouble regions predicted by the physics?
6. Verdict: how much GPU time does adaptive λ-placement save at equal accuracy (or accuracy
   gained at equal cost)?

## Deliverables

- **Fixed vs adaptive comparison** (ΔΔG error vs windows used, per variant) — the money figure.
- The learned sampling pattern: where the GP added windows (plotted on the dG(λ) curve).
- A recommendation for the parent repo: should `lambda_windows: 18` become adaptive?

## Pitfalls

- **The GP must model the right thing.** dG(λ) smoothness is the model; the actual MBAR error
  comes from u_kn *overlap*, which is a different object. Model the overlap/variance, not just
  the mean curve, or you'll add windows in the wrong place.
- **Sequential comparison is noisy.** One adaptive episode vs one fixed ladder is anecdote;
  repeat over many variants and report the distribution of savings.
- **Equal-cost comparison, not equal-accuracy.** "Adaptive used 12 windows and got the same
  error" is the honest claim; "adaptive is more accurate" requires equal window counts.
- **Pre-register the stopping rule.** If the adaptive loop stops when "the answer looks
  converged," you're peeking. Define it before running.
- **Cost of the decision step.** If GP+MBAR refitting per round is more expensive than a
  window, the scheme is pointless; keep the decision step CPU-cheap.
