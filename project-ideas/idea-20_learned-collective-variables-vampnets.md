# Idea 20 — Learning Collective Variables with Deep Autoencoders (VAMPnets)

## Research question

The success of MD analysis hinges on finding the right **collective variables (CVs)** — the
slow coordinates along which conformational change actually happens. Hand-crafted CVs (RMSD,
radius of gyration, native contacts) are guesses. **Time-lagged autoencoders / VAMPnets**
(Mardt et al., PNAS 2018) *learn* the slow coordinates directly from trajectory data by
maximizing the autocorrelation of the learned features.

The question: **Can a neural network discover the slow collective variables of SOD1's
dynamics that RMSD and Q miss — and do those learned variables reveal variant-specific
differences in the landscape?** This is the data-science heart of Idea 19: the MSM says "here
are states"; the VAMPnet says "here is the *right coordinate* in which to see them."

## Why this needs an SGE cluster

```
tasks = many short trajectories (the training data) + hyperparameter sweep
```

- Same trajectory requirement as Idea 19: many independent short runs (50–100 per system) —
  the array pattern. VAMPnets are *data-hungry*: the more trajectories, the better the learned
  coordinates.
- Training the VAMPnet (a small recurrent/convolutional net over featurized frames) has a
  hyperparameter sweep (lag time, hidden size, layers) — the `qsub -t` grid from Idea 14.
- You need enough data to split into train/validation and to test whether the learned CV
  generalizes to *held-out* trajectories — again, the cluster's scale.

## Data & tools

- **Data:** the same trajectory ensemble as Idea 19 (generate once; share the array!).
- **Software:** **PyEMMA** (has VAMPnet tutorials and `vampnet` wrappers) or a custom
  PyTorch implementation of the time-lagged variational principle (the core is ~150 lines).
- **Baselines to beat:** hand-crafted CVs (Q, RMSD, contacts) and linear tICA.
- **Metrics for "good CV":** VAMP-2 score (how much kinetic variance the CV captures) on
  validation trajectories; how well the learned CV separates the MSM states (Idea 19).

## Skill prerequisites

- Python + PyTorch basics (you can build the net, or use PyEMMA's).
- MSM concepts from Idea 19 are a strong prerequisite.
- Advanced-ish for the ML; the MD data collection is the same as Idea 19.

## Cluster budget

| Parameter | Value |
|---|---|
| Trajectory array | shared with Idea 19 (~240–400 GPU runs) |
| VAMPnet training sweep | ~50–100 configs |
| Wall-clock | **~1–2 weeks** (data) + **~1 day** (training) |

## Milestones

1. Reuse/build the Idea-19 trajectory ensemble (featurized frames only — cheap storage).
2. Implement or run a baseline: linear **tICA** CVs; compute VAMP-2 score on held-out
   trajectories.
3. Train a VAMPnet (time-lagged AE) on train trajectories; score on held-out.
4. Project trajectories onto learned vs tICA vs hand-crafted CVs; compare free-energy
   landscapes (does the learned CV show basins the others smear out?).
5. **Variant comparison:** do the learned CVs for a destabilizing mutant vs WT differ in a
   way that identifies a *pathway* (which residues move along the slow coordinate)?
6. Verdict: did learning beat hand-crafted? Where and why?

## Deliverables

- VAMP-2 scores: hand-crafted vs tICA vs learned CV (a clean, quantitative comparison).
- Landscape projections in the learned coordinates (vs baseline) — show the difference.
- An interpretability pass: which physical descriptors load onto the learned slow coordinate
  (attention/saliency over the input features).
- A written story: "the slow coordinate of SOD1 dynamics is X, and variants change it by Y."

## Pitfalls

- **Overfitting the variational objective.** VAMP-2 can be gamed by degenerate solutions;
  always evaluate on *held-out* trajectories and report generalization, not train score.
- **Learned CVs can be uninterpretable.** Balance with the saliency/loading analysis — a
  black-box CV that beats Q but you can't explain is half a result.
- **Lag-time tyranny.** The "slow" coordinate is defined relative to your lag time; show the
  landscape is stable across lag choices.
- **Same-data reuse.** If you share the Idea-19 trajectories, keep the splits identical or you
  can't compare the two analyses.
- **Don't over-engineer the net.** A small net beats a big one here; the data quality matters
  more than the architecture.
