# Idea 15 — Soft-Core Protocol Benchmark

## Research question

**Which soft-core scheme gives the most accurate and stable alchemical FEP for protein
mutations?** The parent repo's config already anticipates this question
(pipeline.yaml:114-127): the classic **Beutler** soft-core (`sc_alpha 0.3`, `sc_power 1`,
`sc_sigma 0.25`) vs the newer **Gapsys** scheme (PMID 26588970), with a note "switch only
after a window has run cleanly." This project is that systematic comparison:

- Run a panel of representative mutations under **both** soft-core functions (identical
  everything else),
- Measure **accuracy** (ΔΔG vs experimental controls) and **robustness** (convergence speed,
  hysteresis forward/reverse, end-state sampling stability — no end-point singularity blowups),
- Recommend a protocol with numbers, not vibes.

This is a **methods-development** project — the least "flashy" and most publishable: protocol
choice is exactly where free-energy pipelines go wrong silently.

## Why this needs an SGE cluster

```
tasks = mutations × soft-core × legs × windows × replicates
      ≈ 12 × 2 × 2 × 18 × 3 = 2,592 GPU windows
```

- The comparison is only meaningful if every arm is run at the **same cost and same
  replicates** — that's an array of thousands of identical GPU windows. The cluster makes the
  "controlled experiment" affordable.
- Convergence/hysteresis statistics need many windows AND replicates (rule 5 of this repo:
  error bars are mandatory) — the cluster provides the scale.
- The two soft-core arms are independent → perfect `qsub -t` split. You can even run both arms
  of one mutation concurrently and compare as they land.

## Data & tools

- **Engine:** GROMACS + pmx (parent stack). Soft-core is a `.mdp` option set:
  - Beutler: `sc-function = beutler`, `sc-alpha = 0.3`, `sc-power = 1`, `sc-sigma = 0.25`.
  - Gapsys: `sc-function = gapsys`, `sc-gapsys-scale-linpoint-lj = 0.85`,
    `sc-gapsys-scale-linpoint-q = 0.3`, `sc-gapsys-sigma-lj = 0.3` (values already in
    pipeline.yaml).
- **Panel:** 12 charge-neutral mutations spanning substitution types (hydrophobic↔polar,
    small↔large, buried vs surface). Reuse the gate subset (F64A, I18V, I113T, A4V, G93A,
    G93S, I149A, G93V) + 4 added for diversity.
- **Experimental anchors:** `data/variants.csv` ΔΔG (apo-monomer) — the accuracy yardstick.

## Skill prerequisites

- **Advanced, and you should already run the parent repo's FEP successfully.** You are
  modifying a protocol, so you must be able to debug a failing `mdrun` and read `dhdl.xvg`
  critically.

## Cluster budget

| Parameter | Value |
|---|---|
| Mutations × soft-core × legs × windows × reps | 12 × 2 × 2 × 18 × 3 = 2,592 |
| Per-window | ~3–6 h GPU |
| Wall-clock on 8 GPUs | **~2–3 weeks** (or cut to 8 mutations ≈ 1.5–2 weeks) |

## Milestones

1. Reproduce one window under **Beutler** (the config default) on a GPU node — the parent
   repo's smoke mode does this.
2. Flip the `.mdp` to **Gapsys**; run the same window; confirm `grompp` accepts the
   gapsys-specific parameters (they are *different* parameter names — mixing them is a grompp
   error, per config note).
3. Define the measurement protocol **before** the array (pre-register):
   - Accuracy: per-mutation ΔΔG vs experiment (MBAR; parent `analyze.py`).
   - Robustness: forward/reverse hysteresis (`max_cycle_closure_kcal: 1.0` config), overlap
     quality of the MBAR u_kn matrix, end-state energy sanity.
4. Launch the 2-arm array; collect NPZs with the soft-core scheme stamped in provenance.
5. Aggregate: per-scheme MAE/RMSE vs experiment, per-scheme hysteresis distribution,
   per-scheme failure rate (windows that didn't converge).
6. **The verdict:** does Gapsys improve accuracy, robustness, or both? Is the switch worth
   it at the parent repo's current config? (Answer with error bars.)

## Deliverables

- Accuracy table: Beutler vs Gapsys ΔΔG vs experiment (MAE/RMSE, per mutation).
- Robustness figure: hysteresis distributions and overlap-quality per scheme.
- Convergence plot: ΔΔG as a function of sampled time per scheme (which converges faster?).
- A written protocol recommendation for the parent repo (`pipeline.yaml` change proposal).

## Pitfalls

- **Don't mix soft-core parameter sets.** Beutler uses `sc_alpha/sc_power/sc_sigma`; Gapsys
  uses `sc_gapsys_*`. A grompp failure at task 1 of 2,592 is embarrassing — smoke-test both
  `.mdp`s first (this is literally what the parent config comment warns about).
- **End-state singularities.** The classic failure mode soft-core fixes is the appearing/disappearing
  atom at λ=0/1. If a scheme's end states still blow up, that IS a finding — report the
  energies, don't hide them.
- **Provenance.** Stamp the soft-core scheme in every NPZ (parent repo rule), or analyze.py
  cannot tell the arms apart and you will average them by accident.
- **Same budget both arms.** Any wall-time/length difference between Beutler and Gapsys arms
  contaminates the comparison. Identical windows, differing only in `.mdp` soft-core lines.
- **Charge-neutral only** (config rule) — a charge-changing mutation adds a finite-size
  artifact that swamps the soft-core effect.
