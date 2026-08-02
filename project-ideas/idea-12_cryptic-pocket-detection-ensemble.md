# Idea 12 — Cryptic-Pocket Detection from an MD Ensemble

## Research question

**Are there transient "cryptic" binding pockets on a protein that the crystal structure
hides?** Many drug targets are considered "undruggable" because their crystal structure shows
no pocket — yet long-timescale MD reveals **transient pockets** that open and close with
conformational dynamics (e.g., the well-known p53–MDM2 and KRAS stories). This project:

1. Runs **many independent MD replicas** of a target protein,
2. Detects pockets in every frame (**fpocket** or similar),
3. Computes **pocket occupancy/volatility statistics** across the ensemble,
4. Ranks "cryptic" pockets by how often and how long they open,
5. (Stretch) docks a fragment library into the top-ranked cryptic pocket to show it is
   druggable.

A great target framing: *"Predicting cryptic pockets on a protein that currently has no known
drug binding site."*

## Why this needs an SGE cluster

```
tasks = replicas (e.g., 32 × 1 µs) + pocket analysis per frame
```

- Cryptic pockets appear on **slow (µs) timescales** and are **rare events** — you need many
  long independent replicas to sample them, not one heroic run. 32 replicas × 1 µs of atomistic
  GPU MD is ~30 GPU-days: the cluster's whole reason to exist.
- **1 array task = 1 replica.** Each replica is independent, resumable, and identical in cost —
  the perfect `qsub -t` shape. The aggregate pocket statistics then combine across all replicas.
- The pocket-detection pass (fpocket over thousands of frames) is also trivially parallel
  (1 task per replica's frame-slice).
- You can iterate: if no pockets open, extend replicas or raise temperature, re-submit the
  array — cheap on the cluster.

## Data & tools

- **Target:** pick a protein with an interesting allostery/undruggable story and a good crystal
  structure. Examples: KRas, TEM-1 β-lactamase, a GPCR (harder), a viral protease, or SOD1
  (does the apo/reduced form open cryptic pockets near the aggregation-prone loops? — ties to
  the parent repo!).
- **Engine:** GROMACS (atomistic; the parent stack). Multiple seeds, fixed T (e.g., 310 K) or
  a couple of temperatures.
- **Pocket detection:** **fpocket** (free, standard) — per-frame pockets with volumes and
  "druggability" scores.
- **Analysis:** cluster pockets spatially across frames (match by residue lining), compute per-pocket
  (open frequency, mean volume, residence time), and make a pocket-volatility heatmap per residue.

## Skill prerequisites

- Solid MD literacy (setup, run, restart, analyze trajectories).
- Python + some structural bio (residue contacts, pocket vocabulary).
- Intermediate; the MD itself is not as hard as the analysis design.

## Cluster budget

| Parameter | Value |
|---|---|
| Replicas | 32 × ~1 µs (2 seeds × 16 starting conformations, or 32 seeds) |
| Per-replica | ~1–2 days GPU |
| Pocket pass | 1 task per replica (frame slice) |
| Wall-clock on 8 GPUs | **~1–2 weeks** |

## Milestones

1. Choose target; fetch structure; build system (parent repo prep is reusable); verify a short
   run is stable.
2. Design the ensemble (seeds, temperatures, optional: different initial conformations from a
   short high-T pre-equilibration).
3. Validate on a **positive control**: a protein with a known cryptic pocket (e.g., a published
   case) — your pipeline must detect it before you trust it on your target.
4. Array over replicas; collect trajectories (retention policy: keep processed pocket lists,
   drop raw frames or downsample).
5. fpocket every Nth frame per replica (parallel); map pockets to residues; cluster into
   persistent pocket identities.
6. Rank pockets: high open-frequency × long residence × reasonable volume = cryptic druggable
   candidate.
7. (Stretch) Dock a fragment library (Idea 04 tooling) into the top cryptic pocket and show it
   engages transiently (short MD).

## Deliverables

- **Pocket-volatility map** (residue × time of pocket opening) — the money figure.
- Ranked cryptic-pockets table (open frequency, mean volume, lining residues, druggability).
- Comparison vs crystal-structure pockets (which are new?).
- If you did the stretch: fragment-docking evidence the cryptic pocket can bind a ligand.

## Pitfalls

- **Positive control first.** Without a known cryptic-pocket protein validating the pipeline,
  "no pocket found" is meaningless (could be protocol failure).
- **Rare-event statistics.** One replica that opens a pocket is an anecdote; report per-pocket
  open frequencies with binomial confidence intervals across replicas.
- **Detection artifacts.** fpocket flags shallow grooves; apply a minimum volume + minimum
  residence threshold and *say what you chose*.
- **"Estimate-only" retention.** 32 µs of trajectories is terabytes; keep pocket/frame features,
  not frames. Your analysis must not need the raw trajectory after feature extraction.
- **Temperature bias.** Raising T opens pockets more; report at physiological T as the headline
  and treat high-T as a discovery screen.
- **Cryptic ≠ druggable.** A pocket that opens rarely may still bind poorly; the docking
  stretch is what turns "a pocket exists" into "a pocket can bind."
