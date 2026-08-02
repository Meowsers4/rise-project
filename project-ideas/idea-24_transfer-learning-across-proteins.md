# Idea 24 — Transfer Learning: Can One Protein's Physics Predict Another's?

## Research question

The parent repo spends GPU hours computing FEP ΔΔG for SOD1. Most stability predictors are
trained on *many* proteins. **Transfer learning** asks: does training an ML model on the
cluster-generated stability data of one protein (or a handful) improve predictions on another
protein with sparse data?

Concretely: for a *new* protein with only ~20 experimental ΔΔG points, is it better to
(a) train a model only on those 20 points, (b) train on 3,000 FoldX points from SOD1, or
(c) fine-tune a generic model with a few of the new protein's points (transfer)? The twist:
**physics tells us transfer should be limited** (stability effects are structure-specific), but
ML transfer often wins on small data. This project measures the actual trade-off curve.

## Why this needs an SGE cluster

```
tasks = generate source-protein labels (array) + train/test across many splits
```

- The "source" data — thousands of labeled mutations for one or more proteins — is generated
  by the saturation-scan arrays (Idea 01). Without the cluster you'd be restricted to small
  public datasets.
- Testing transfer *properly* means evaluating across **many train/test splits and several
  target proteins** (a nested CV grid) — hundreds of independent training runs, each a
  `qsub -t` task.
- You want multiple source proteins (2–4) to see if transfer depends on which source you
  pick — more arrays.

## Data & tools

- **Source labels:** FoldX saturation scans (Idea 01) on 2–4 proteins of differing fold
  similarity to the target.
- **Target proteins:** 2–3 proteins with *experimental* ΔΔG (variants.csv controls, DMS,
  ProTherm) — the "test" side.
- **Models:** gradient boosting / GNN (Idea-14 machinery); feature sets that transfer
  (BLOSUM, burial, SASA — structural features must be recomputed per protein).
- **Framing:** three regimes — no transfer (target-only), direct transfer (train on source,
  score target), fine-tune (pretrain on source, fit on a few target points).

## Skill prerequisites

- Python + the Idea-14 ML stack.
- Comfortable thinking about data regimes (how many target points before transfer stops
  helping).
- Intermediate-advanced.

## Cluster budget

| Parameter | Value |
|---|---|
| Source scans | 2–4 × ~3,000 tasks |
| Train/test grid | ~300–600 runs |
| Wall-clock | **~2–4 days** |

## Milestones

1. Generate FoldX scans for source proteins; collect experimental targets.
2. Define the transfer regimes and the evaluation protocol (fixed splits, pre-registered).
3. Run the grid: for each (source, target, regime, #target-points), train + evaluate.
4. **The transfer curve**: model error on the target vs number of target points, for each
   regime. Where does fine-tuning beat both no-transfer and direct-transfer?
5. Does source–target *structural similarity* (fold, sequence identity) predict transfer
   benefit? (The physical hypothesis.)
6. Verdict: when is it worth spending a cluster to label a source protein for transfer?

## Deliverables

- **Transfer curves** (error vs target-points, per regime and source) — the money figure.
- Transfer-benefit vs structural-similarity scatter.
- A practical recommendation: "to predict protein X cheaply, label protein Y like this."

## Pitfalls

- **Leakage by homolog.** If source and target are homologous, "transfer" is just extra
  training data — still useful, but report sequence identity so the reader knows.
- **Structural features don't transfer as-is.** SASA/burial must be recomputed per target
  structure; a shared-feature mistake is a silent bug.
- **Fine-tune without catastrophe.** Standard ML fine-tuning can destroy pretrained weights
  (catastrophic forgetting); the small-target-data regime is exactly where this bites. Use a
  small fine-tune LR.
- **Don't pick sources that trivially win.** Pre-register source/target pairs; report the
  "boring" cases (similar-fold transfer) as the baseline, not as the finding.
- **Report distributions, not single splits.** One train/test split proves nothing; the grid
  is the experiment.
