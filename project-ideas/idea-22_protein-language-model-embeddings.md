# Idea 22 — Protein Language Model Embeddings as Stability Features

## Research question

Protein language models (pLMs) — neural networks pre-trained on millions of protein sequences
(e.g., ESM-2 by Meta) — produce per-residue embeddings that capture evolutionary and
functional context. A hot (and heavily contested) claim is that **pLM zero-shot scores can
rank mutation effects without any MD or structure**. This project tests that claim on the
parent repo's exact question:

**Do ESM-2 embeddings/log-likelihoods predict SOD1 (and other) stability changes as well as
structure-based physics (FoldX, Rosetta)? And does combining pLM features with physical
features beat either alone?**

The novelty is the *hybrid*: most people either (a) throw physics at it or (b) throw a pLM at
it. You test whether the two sources of signal are redundant or complementary — a question the
literature genuinely argues about.

## Why this needs an SGE cluster

```
tasks = sequence scoring (per mutation) + feature extraction + model sweeps
```

- Scoring every single mutation with a pLM is a GPU pass over the protein sequence; doing it
  for a *family* of proteins (all SOD1 homologs, or a pathogen's proteome) is thousands of
  GPU tasks → the array.
- The pLM itself (ESM-2 ~650M params) is too big for a laptop GPU comfortably; the cluster's
  GPUs are the point.
- The comparison "pLM vs physics vs hybrid" is a model-selection sweep — `qsub -t` grid
  (Idea-14 machinery).

## Data & tools

- **pLM:** **ESM-2** (Meta, open weights; `fair-esm` on PyPI / HuggingFace). Zero-shot score =
  log-likelihood of the mutant sequence under the model (masked-marginal or pseudolikelihood
  approximation — the standard "mutation effect" estimate).
- **Physics baseline:** FoldX (Idea 01) or the parent repo's FEP for a small gold set.
- **Labels:** experimental ΔΔG (variants.csv / DMS / ProTherm).
- **Models:** logistic/GBM on features: [ESM score, ESM embedding distance, FoldX ΔΔG,
  structural descriptors]. Compare single-source vs hybrid.

## Skill prerequisites

- Python + PyTorch basics (to run a pretrained ESM, not train one — training is out of scope).
- Understands embeddings (what a vector per residue means).
- Intermediate-advanced: the ML framing matters as much as the physics.

## Cluster budget

| Parameter | Value |
|---|---|
| Sequences scored | 1 protein's saturation scan (~3,000) to a whole family (~10^4) |
| Per-sequence GPU | seconds–minutes |
| Model sweeps | ~100 configs |
| Wall-clock | **~1–2 days** (scoring) + **~1 day** (models) |

## Milestones

1. Choose a validation set with experimental ΔΔG (variants.csv gate + a DMS).
2. Run ESM-2 zero-shot scoring over the full saturation scan; compute per-mutation score.
3. Baseline: FoldX over the same scan (Idea-01 array, shared if possible).
4. Fit predictors: (a) physics only, (b) ESM only, (c) hybrid. Report Pearson/MAE on
   held-out mutations.
5. **The decisive test:** does the hybrid beat the max of the two singles by a real margin?
   Where does each source of signal win (e.g., ESM good on buried hydrophobic→charged, physics
   good elsewhere)?
6. (Stretch) check on a *second* protein to see if the conclusion generalizes.
7. Verdict: are pLMs and physics redundant or complementary for stability prediction?

## Deliverables

- **Single-source vs hybrid benchmark** (Pearson/MAE) — the money figure.
- Error anatomy: where ESM wins and where physics wins (a per-class table).
- A short, honest discussion of the "pLM vs physics" debate grounded in your numbers.

## Pitfalls

- **Zero-shot pLM scores correlate with function AND stability AND expression.** They are not
  pure stability; a "hit" may be a functional/expression effect. Say this explicitly.
- **Leakage in the pLM.** ESM was pre-trained on all of UniProt — your test mutations are
  *in its training data* (as sequences). That's not leakage the way you can fix; it means
  "zero-shot" is really "memorized-ish." Report it as a caveat, and lean on the *hybrid
  delta* rather than absolute claims.
- **Model size ≠ better.** ESM-2 has many sizes; the small ones run on one GPU and often score
  comparably. Don't pay for the 15B unless you test it.
- **Beware sign conventions** between log-likelihood (higher = more favorable) and ΔΔG
  (higher = destabilizing). Map them once, in the README.
- **Pre-register the split.** Same leakage discipline as Idea 16: split by position/protein.
