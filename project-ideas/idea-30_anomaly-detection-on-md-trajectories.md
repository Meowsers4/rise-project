# Idea 30 — Anomaly Detection on MD Trajectories: Finding the Rare Events

## Research question

Long MD trajectories are mostly boring (the protein jiggles near its native state); the
interesting events — partial unfolding, a transient exposed loop, a near-misfolding excursion,
an intermittent salt-bridge loss — are **rare**. Standard statistics (mean RMSD, RMSF) wash
them out. **Anomaly detection** (unsupervised ML for "this frame is unusual") turns the
question around: instead of defining "interesting" in advance, learn what "normal" looks like
and flag deviations.

The research question: **Can unsupervised anomaly detection on MD trajectories identify
functionally/disease-relevant rare events that conventional order-parameter analysis misses —
and do disease variants change the anomaly budget (how often, how big, where)?**

A clean framing: run plain MD on WT and 2–3 variants of a disease-relevant protein (SOD1!),
train an autoencoder (or one-class SVM / isolation forest) on the WT frames, then ask: do the
variant trajectories trigger anomalies more often, and are the anomalous frames localized to
known disease-relevant regions (loops, interface)?

## Why this needs an SGE cluster

```
tasks = trajectories (the data) + anomaly-scoring passes (per trajectory) + training sweeps
```

- You need **many long trajectories** per system to accumulate rare events (the array pattern;
  e.g., 30–60 × 200 ns per variant). One trajectory almost never contains a rare event.
- Anomaly *scoring* every frame of every trajectory is a separate parallel pass (1 task per
  trajectory).
- Training/validation of the anomaly model (which architecture, which features) is a small
  `qsub -t` sweep.
- The variant comparison (anomaly rates across WT vs mutants) is the science, and it needs the
  trajectory ensemble only the cluster can produce.

## Data & tools

- **Data:** plain unbiased MD (GROMACS, parent stack) of the apo-reduced monomer, WT + 2–3
  variants; 30–60 × 200 ns per system.
- **Features per frame:** contact-map slice, per-residue RMSF window, local SASA, dihedral
  PCA (dPCA) features, or raw coordinates after alignment — choose features that capture
  *conformational* anomaly, not global drift.
- **Anomaly models (free, sklearn/PyTorch):** autoencoder (reconstruction error = anomaly
  score), one-class SVM, isolation forest. Autoencoder is the most interpretable for
  "which residues are anomalous" (per-residue reconstruction error).
- **Validation:** compare flagged frames against physical labels where they exist (e.g., a
  known partial-unfolding event visible in Q-fraction), and check anomaly rates vs variant
  severity.

## Skill prerequisites

- Comfortable running MD (Idea 08/09/19 level) — this project needs real trajectories.
- Python + basic unsupervised ML (what an autoencoder reconstruction error means).
- Intermediate-advanced.

## Cluster budget

| Parameter | Value |
|---|---|
| Trajectories | 30–60 per system × 4 systems = 120–240 GPU runs (~200 ns) |
| Anomaly scoring | 1 task per trajectory (fast, feature-based) |
| Training sweep | ~50 configs |
| Wall-clock on 8 GPUs | **~2–4 weeks** |

## Milestones

1. Build systems (parent prep); run a small pilot (3 trajectories/variant) to validate that
   MD is stable and the feature set is discriminating.
2. Train the anomaly model on WT frames (development set); tune on held-out WT frames.
3. Score all trajectories (parallel pass); define the anomaly-rate metric per trajectory
   (fraction of frames above threshold × their duration).
4. **Compare variants:** anomaly rate, peak anomaly location (per-residue), and duration
   distributions vs WT. Do the disease variants show more/longer anomalies, and are they
   localized to the loops/interface?
5. Cross-validate against a physical order parameter (Q-fraction, salt-bridge monitor): do the
   flagged frames correspond to real structural excursions?
6. Verdict: what did anomaly detection find that RMSF/Q missed?

## Deliverables

- **Anomaly-rate comparison** (WT vs variants, with distributions) — the money figure.
- Per-residue anomaly map (which residues are anomalously flexible/mobile per variant).
- A case-study of the top anomalous event per variant with structural snapshots.
- A written "can anomaly detection triage disease variants faster than hand-picked order
  parameters" conclusion.

## Pitfalls

- **Define "normal" on WT only, then score variants.** If you train on a mix of WT+mutants,
  the model "learns" the disease signal and the comparison is circular.
- **Global drift fakes anomalies.** Align frames and detrend slow global motions (use
  contact/RMSF features, not raw coordinates) or a whole trajectory looks anomalous.
- **Threshold arbitrariness.** Report anomaly rates over a range of thresholds (ROC-style) and
  show the variant ordering is robust.
- **Rare events are rare even in 200 ns.** 30 trajectories per variant is a floor; report the
  variance across trajectories (per-trajectory anomaly rate spread), not just the mean.
- **Don't over-interpret one anomalous frame.** An anomaly is a flag for investigation, not a
  mechanism; follow up the top events with structural analysis and an MSM (Idea 19) if they
  recur.
