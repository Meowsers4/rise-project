# Idea 02 — Consensus Force-Field FEP Benchmark

## Research question

**Do you get a more accurate ΔΔG by averaging alchemical free-energy calculations over several
molecular-mechanics force fields?** Both Gapsys et al. 2016 (PMID 27122231, 762-mutation scan)
and Wells et al. 2021 (PMID 34076319, SOD1) report that averaging ΔΔG over 2–3 force fields
reduces error to near-experimental reproducibility (~0.8–1.0 kcal/mol AUE). This project
**reproduces and quantifies that claim on your own system**: run the same mutations with
GROMACS + pmx under every force field in `pmx/data/mutff45`, and measure whether the
consensus (mean or median) beats any single force field.

The natural test bed is the **SOD1 gate subset already in this repo** (A4V, G93A, G93S, G93V,
F64A, I18V, I113T, I149A) — their experimental ΔΔG are already in `data/variants.csv`, so you
get a validation correlation for free.

## Why this needs an SGE cluster

```
tasks = variants × force_fields × legs × windows × replicates
      = 8 × 4 × 2 × 18 × 3 = 3,456 GPU windows
```

- Each window is a **3 ns GROMACS mdrun on its own GPU** (~3–6 h). That is weeks of continuous
  single-GPU compute → the cluster turns it into ~1–2 weeks of wall-clock on 8 GPUs.
- This is a **controlled, pre-registered experiment**: the only thing that changes between arms
  is the force field. Every arm is embarrassingly parallel.
- The fan-out geometry is *identical* to this repo's Stage 3 (`(variant, leg, window,
  replicate)`), so `scripts/submit_array.sh` is directly reusable — just add a `FF` dimension
  to the decode.

## Data & tools

- **Engine:** GROMACS + pmx — **this repo's existing Stage 3 stack** (`src/fep/pmx_engine.py`).
- **Force fields available** in `pmx/data/mutff45` (config/pipeline.yaml:107-108):
  - `amber99sbmut`
  - `amber99sb-star-ildn-mut` (the current config default — pmx's own benchmark FF)
  - `charmm22starmut`
  - `oplsaamut`
- **Water model:** pmx's `mutff45` ships TIP3P (keep `tip3p`; config already locks this).
- **Experimental controls:** `data/variants.csv` (apo-monomer ΔΔG), the same numbers the
  parent project's validation gate uses.

## Skill prerequisites

- Comfortable with the parent repo: config, `src/fep/window.py`, `submit_array.sh`.
- Understands "ΔΔG = ΔG_folded − ΔG_unfolded" and can read a `dhdl.xvg`.
- Some stats: Pearson correlation, RMSE, mean absolute error.

## Cluster budget

| Parameter | Value |
|---|---|
| Variants × FF × legs × windows × reps | 8 × 4 × 2 × 18 × 3 = 3,456 |
| Per-window | ~3–6 h GPU (3 ns) |
| Wall-clock on 8 GPUs | **~2–3 weeks** (or cut to 8 gate variants × 3 FFs ≈ 2,592 → ~2 weeks) |
| Output per window | ~200 KB NPZ (estimates only — no trajectories) |

**Cost-control lever (legit):** this is a *comparison of methods*, so every arm must be equal
— but you can legitimately cut `lambda_windows` to 18 (keep) and use 3 replicates (the floor).

## Milestones

1. Verify one real (non-mock) window runs under `amber99sb-star-ildn-mut` on a GPU node
   (this repo's `--smoke` mode is designed for exactly this).
2. Extend the engine/config so `fep.pmx_forcefield` is an array task dimension rather than a
   single config value (a small, contained change to `pmx_engine.py` + `submit_array.sh`).
3. Launch the 4-FF × 8-variant array; collect NPZs; run `src/fep/analyze.py` per
   (variant, FF) to get per-FF ΔΔG.
4. **Consensus:** per variant, average ΔΔG across FFs (and compute the median). Propagate
   errors (add variances / N², or block bootstrap across FFs).
5. Validation: correlate per-FF ΔΔG and consensus ΔΔG vs experimental controls.
   Compute Pearson, Spearman, RMSE, MAE per arm.
6. **The money figure:** a plot of RMSE/MAE vs experiment for each FF and for the consensus —
   showing the consensus lands below every single FF (Gapsys' claim).

## Deliverables

- Per-FF ΔΔG table + consensus ΔΔG table (with uncertainties) for the 8 variants.
- Bar chart: MAE/RMSE of each FF vs consensus vs experiment.
- Scatter: consensus ΔΔG vs experimental ΔΔG (with the gate thresholds from config drawn on).
- A verdict: does consensus help here? How big is the gain? Does it agree with Gapsys' 762-mutation
  conclusion at the scale of one protein?

## Pitfalls

- **Charge-changing mutations.** The parent repo deliberately excluded them from the gate
  (PME finite-size artifact, config:161-166). Keep the benchmark charge-neutral; you can
  discuss charge-changing as a separate, documented arm.
- **Provenance.** Every window NPZ is provenance-stamped; `analyze.py` refuses to mix engines.
  When you add an `FF` dimension, stamp it too — otherwise you will accidentally average
  across force fields in one MBAR call.
- **Force-field-specific pdb2gmx quirks.** Each FF in `mutff45` has its own `.ff`; the engine
  resolves `GMXLIB` once. Verify `pdb2gmx` succeeds for every FF before launching the full
  array (one smoke task per FF).
- **Don't tune on the test set.** Decide the consensus rule (mean vs median, which FFs) BEFORE
  looking at the validation numbers, or the result is overfit.
