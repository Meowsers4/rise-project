# Idea 05 — Deep-Mutational-Scan Validation

## Research question

**How well do empirical ΔΔG predictors reproduce a published experimental deep mutational
scan (DMS)?** A DMS measures the functional/stability effect of ~1,000+ single mutations of
one protein in one experiment. Your job: predict those mutations with FoldX/Rosetta and
quantify agreement (Pearson, Spearman, ROC on "destabilizing yes/no"), then dissect WHERE
the predictor fails (buried vs surface, hydrophobic vs charged substitutions).

This is the ideal **"method validation"** project: the experimental data is free, huge, and
public; the computation is trivially parallel; and the scientific content is a rigorous
benchmark plus an error analysis — exactly the shape of a publishable methods paper.

## Why this needs an SGE cluster

```
tasks = ~1,000 mutations (or 19 × N positions)
```

- **1 array task = 1 mutation** → ~1,000 independent FoldX runs (~5 min each) = **~80 CPU
  hours**. Weeks on a laptop, hours on the cluster.
- Same fan-out as Idea 01, but the *validation* is the science rather than a side check. The
  cluster buys you the whole landscape so the benchmark is not cherry-picked.
- You can afford **ensemble scans** (MutateX runs each mutation on multiple relaxed
  structures) and **two predictors** (FoldX + Rosetta cartesian_ddg) for a consensus-vs-single
  comparison — the compute makes the study better, not just faster.

## Data & tools

- **Experimental DMS datasets (free, published):**
  - **Fowler & Fields** comprehensive collection of DMS datasets (Nat Methods 2014).
  - Protein-specific DMS: TEM-1 β-lactamase, GBP1, GFP, SUMO, yeast Hsp90, influenza HA.
  - Choose one with **well-annotated single-mutant effects** (not combinatorial libraries).
- **Predictors:** FoldX `BuildModel` (primary), Rosetta `cartesian_ddg` (cross-check),
  **MutateX** (PMID 35323860) to automate.
- **Structure:** the protein's PDB (or ColabFold prediction if none exists).

## Skill prerequisites

- Python/pandas + basic statistics (correlation, ROC).
- Understanding that "functional score" ≠ "stability ΔΔG" — a key analysis step.

## Cluster budget

| Parameter | Value |
|---|---|
| Mutations | ~1,000 |
| Per-mutation | ~5 min CPU (×2 predictors) |
| Array size | 1 per mutation |
| Wall-clock (30 concurrent) | **~12–24 h** |

## Milestones

1. Pick the DMS dataset; download; parse to a clean `(variant → experimental score)` table.
2. Obtain the protein structure; map DMS residue numbering → PDB numbering (the classic trap).
3. Run FoldX on ~50 mutations locally; sanity-check sign convention and scale.
4. Array over all mutations (FoldX); parallel array over all mutations (Rosetta) if time.
5. Merge: for each mutation, experimental score vs (FoldX ΔΔG, Rosetta ΔΔG, consensus).
6. Metrics: Pearson/Spearman on raw values; **rank correlation** (Spearman often better for
   noisy experimental scores); AUC for binary "destabilizing vs not".
7. Error anatomy: split by burial (SASA), substitution class (hydrophobic↔polar, charge),
   and secondary structure. What correlates with predictor failure?
8. Write it up as a benchmark figure + error table.

## Deliverables

- **Correlation/ROC figure** (predicted vs experimental, both predictors + consensus).
- **Error-anatomy table**: MAE broken down by burial / substitution class.
- A "where you can and cannot trust the predictor" section — the most competition-worthy part.

## Pitfalls

- **DMS scores are often functional, not stability.** A stop-codon or non-native mutation
  scores "low" for expression reasons, not folding. Separate the analysis by mutation class and
  say so.
- **Numbering mismatch** (assay numbering vs PDB) is the #1 silent error — map once, verify on
  known mutations (e.g., mutations with both functional and structural data).
- **Single structure vs ensemble.** A crystal structure represents one conformation; MutateX
  solves this by relaxing. If you only use one structure, your buried-residue predictions will
  be systematically noisy.
- **Leakage.** If you "tune" the predictor on the test DMS (e.g., calibrate the cutoff) the
  benchmark is meaningless. Pre-register: report the predictor's default output; calibration is
  a separate, clearly-labeled arm.
