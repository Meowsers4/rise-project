# Idea 14 — ML Surrogate Trained on Cluster-Generated Stability Data

## Research question

**Can a machine-learning model learn the output of expensive physics and predict protein
stability changes instantly?** The parent repo (and Idea 01/02/05) generates ΔΔG values that
cost CPU/GPU hours each. An ML surrogate — trained on a **cluster-generated dataset** — could
predict ΔΔG in milliseconds. This project:

1. **Builds a training dataset with the cluster** (thousands of FoldX/mutation-scan ΔΔG
   values, or hundreds of FEP ΔΔG values with provenance),
2. **Trains ML models** (gradient boosting on sequence/structure features; or a graph neural
   network on the structure graph),
3. **Benchmarks the surrogate vs the physics** on held-out mutations.

The twist that makes this *research* rather than a tutorial: the parent repo (and this repo's
rule 4) already treat **provenance and anti-overfitting as sacred** — apply the same rigor to
the ML split (leakage-free by protein, not by mutation).

## Why this needs an SGE cluster

```
tasks = dataset-generation jobs + hyperparameter sweep (hundreds of training runs)
```

- **Dataset generation is the parallel part:** 3,000 mutation-scans (Idea 01) or an FEP panel
  (Idea 02) — each is an independent array task. This is what makes the training set
  *cheap enough to be large*.
- **The hyperparameter sweep is the second parallel part:** training 200–500 model configs
  (model class × features × depth × regularization × seed) is a classic `qsub -t` grid. Each
  training run is independent.
- The cluster gives you the *whole loop*: generate data (arrays), tune models (arrays),
  score held-out predictions (arrays). A laptop can do this loop once at small scale; the
  cluster makes it a proper study.

## Data & tools

- **Training labels:** your own cluster-generated ΔΔG (FoldX scan from Idea 01 = cheap, huge;
  FEP from Idea 02 = expensive, small, gold-standard). If you need more, public ΔΔG datasets
  exist (S2648, Q3488, Mega-scale DMS), but *your own* provenance-stamped data is the
  differentiator.
- **Features:**
  - Sequence: one-hot AA, AA-index physicochemical properties.
  - Structure: residue burial (SASA), secondary structure, contact count, mutant-vs-WT
    substitution (BLOSUM), distance to the mutation.
  - Graph: a residue-contact graph (node = residue, edges = contacts) → **GNN** input.
- **Models:** `scikit-learn` GradientBoosting / XGBoost (fast, robust) and/or
  `torch_geometric` GNN (state of the art, more setup).
- **Evaluation:** strict nested CV; **split by protein** (never mutate-level leakage) for
  any multi-protein dataset.

## Skill prerequisites

- Solid Python; basics of supervised ML (train/val/test, features, regularization).
- Intermediate-advanced; the stats discipline matters more than the models.

## Cluster budget

| Parameter | Value |
|---|---|
| Dataset generation | Idea 01 array (~3,000 tasks, CPU) or FEP panel (~hundreds of GPU tasks) |
| Hyperparameter sweep | ~300 training runs |
| Per-training-run | ~5–30 min CPU (boosting) / GPU (GNN) |
| Wall-clock | **~1–2 days** (generation) + **~1 day** (tuning) |

## Milestones

1. Generate (or borrow) a labeled ΔΔG dataset with provenance; split into train/val/test
   with a **protein-level** split if multi-protein.
2. Engineer a baseline feature set; train a simple model locally; get a working CV loop.
3. Baseline the target: how well does the *physics* (FoldX/FEP) predict the same held-out
   labels? (The surrogate only matters if it is competitive with or complements physics.)
4. Array the hyperparameter sweep; collect per-config CV metrics; pick the best config on
   validation only.
5. **Blend physics + ML:** does a hybrid (physics prediction as an input feature) beat either
   alone? This is often the strongest result.
6. Error analysis: where does the surrogate fail (rare substitutions, buried sites)?

## Deliverables

- Surrogate vs physics vs hybrid benchmark table (Pearson/MAE on held-out).
- Learning-curve figure: surrogate accuracy vs training-set size (shows the cluster's value —
   more data helps how much?).
- Feature-importance ranking (what matters most: burial? substitution class?).
- A fast prediction tool: the trained model scores new variants in milliseconds.

## Pitfalls

- **Leakage kills this project.** If a mutation appears in both training and test (because
  they're correlated variants of the same position), the metrics lie. Split by protein, and
  for single-protein data, split by position groups carefully.
- **"ML beat physics" is easy to fake.** Make sure the physics baseline is honestly computed
  on the *same* split with the *same* labels; report both, and prefer the hybrid framing.
- **Rare classes.** Some substitutions (e.g., Trp→Gly) are rare in training; report
  per-class errors, not just overall R.
- **Reproducibility:** fix seeds, record feature versions, and (parent-repo rule) stamp every
  prediction with the model + data version that made it.
- **Don't over-claim generalizability.** A surrogate trained on one protein predicts that
  protein best; extrapolation to new folds is a claim you must test, not assume.
