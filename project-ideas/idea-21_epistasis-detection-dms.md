# Idea 21 — Epistasis Detection: When Is a Double Mutant Not the Sum of Its Parts?

## Research question

A foundational assumption of single-mutant ΔΔG scans is **additivity**: that the effect of a
double mutant equals the sum of its parts (ΔΔG(A+B) = ΔΔG(A) + ΔΔG(B)). This is the assumption
behind nearly every "combinatorial variant" inference, including clinical variant triage. But
**epistasis** — non-additive coupling between mutations — is real, and it is precisely where
single-mutant predictors break.

The research question: **Can we systematically detect and model epistasis from deep mutational
scanning (DMS) data, and connect the detected interactions to physical contacts?** Concretely:
fit an additive model to a DMS with double mutants, measure the residual (epistasis) per pair,
and test whether the *structural distance* between the two mutated positions predicts the
magnitude of epistasis. Do close contacts show stronger epistasis? Do "charge-swap" pairs
show systematic non-additivity?

## Why this needs an SGE cluster

```
tasks = epistasis screening over position pairs (thousands of FoldX double-mutant runs)
```

- Measuring epistasis *in silico* means computing double mutants — **~N² pairs**, far more
  than single-mutant scans. A 200-residue protein has ~20,000 position pairs; each is a FoldX
  double-mutant BuildModel (~5 min) → **~1,500 CPU-hours**, only feasible as an array.
- DMS datasets with double mutants exist (e.g., GB1, TEM-1, protein G, influenza HA) but are
  limited; the cluster lets you *generate* the double-mutant matrix yourself at scale and
  cross-validate against whatever experimental doubles exist.
- The modeling side (hierarchical/epistatic regression, structure-aware kernels) is CPU grid
  work, also parallelizable per fold.

## Data & tools

- **Primary:** FoldX `BuildModel` on single AND double mutants (MutateX supports this; Ideas
  01/05 tooling). A curated protein with known allosteric/interface coupling (e.g., an enzyme
  active site) is ideal.
- **Validation:** experimental DMS with double mutants — **GB1** (Wu et al.) and **protein G
  / TEM-1** are standard free datasets.
- **Analysis:** fit ΔΔG ≈ Σ singles + Σ ε_ij pair-terms; estimate ε_ij; correlate ε_ij with
  structural features (contact distance, same-domain, interface, charge of the pair).
- **Tools:** python/scipy (nested fits), `sklearn` (ridge/ElasticNet for sparse ε), PyRosetta
  or FoldX for generation.

## Skill prerequisites

- Python + stats (model comparison, residual analysis).
- Structural biology basics (contact distance, interface).
- Intermediate; comfortably extends Idea-01/05.

## Cluster budget

| Parameter | Value |
|---|---|
| Double-mutant pairs | ~5,000–20,000 |
| Per-pair | ~5–10 min CPU |
| Array size | 1 per pair (or per position row) |
| Wall-clock (50 concurrent) | **~1–2 days** |

## Milestones

1. Choose the protein + DMS validation set; build the single-mutant baseline (Idea-01 style).
2. Fit the additive model on singles; measure residual structure in doubles.
3. Generate the double-mutant array (only pairs that exist / matter, or all if feasible).
4. Estimate epistasis ε_ij per pair; make the **epistasis heatmap** (position × position).
5. Correlate ε_ij with structural features: distance, interface, charge, secondary structure.
6. Validate: do the *experimental* DMS double-mutant residuals agree with your computed
   epistasis? (This is the publishable check.)
7. Write up: which residue pairs are "hot" for epistasis, and what structural rule explains
   them?

## Deliverables

- **Epistasis heatmap** (pairwise ε_ij) with structural annotation — the money figure.
- A table of the strongest epistatic pairs with their structural context.
- Correlation: computed ε vs experimental ε (if doubles exist) or vs structural distance.
- A one-line "rule" for when additivity fails (e.g., "pairs in contact are non-additive").

## Pitfalls

- **Data volume blows up fast.** N² pairs is the trap; restrict to *interesting* pairs (near
  in structure, or all pairs within a distance cutoff) and state the restriction.
- **Sign conventions.** Define ΔΔG(A+B) − (ΔΔG(A)+ΔΔG(B)) consistently, once, in the README.
- **Experimental doubles are rare.** Most DMS cover singles only; treat the double-mutant DMS
  (GB1 etc.) as precious validation, not a primary label source.
- **Additivity is the null hypothesis, not an assumption.** The interesting result is *where
  it fails*, so design the analysis to detect failure honestly (residual plots, not just a
  correlation).
- **Epistasis ≠ interaction energy.** Computed ΔΔG double-mutant residuals mix true coupling
  with force-field noise; report error bars on ε and treat sign carefully.
